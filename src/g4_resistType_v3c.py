'''
Resistance prediction for gram negative, to work with the combined resistance database curated by Nicole Stoesser

Author
Hang Phan - hang.phan@ndm.ox.ac.uk
13 July 2015
Last update 23/12/2015 - adding copy number variant checking
Last update 15/1/2015 - change concensus sequence module to use GATK, not own script to make it more stable, output sequences with padded 'N'
                      - change refineGeneSearch module to make it more compact
Want to do: aim to perform this for metagenomics data where assemblies are not easily available
   1. first map to gene database, get all read pairs that have gene                
   2. use reads that mapped to make assembly
   3. continue the pipeline as normal 
The pipeline works as follow:
1. velvet assembly from short reads, kmer size = 2/3 read length
2. blastn contigs against resistance allele database
3. look for exact match + protein level exact match between contigs and database, record these
4. exclude exact matches from the reference set for use in the next step of bwa alignment
5. bwa alignment + concensus sequence making using samtools. The reference of the alignment is the reduced set of alleles (by cdhit 90% identity)
6. Matching of concensus sequences against the database to get exactMatch/closest hit
7. Linking resulting matches with database and give predictions

# Software versions for use with this pipeline
        1. samtools 0.1.19-44428cd
        2. bcftools 0.1.19-44428cd
        3. bedtools v2.19.1
        4. bwa 0.7.10-r789
        5. velvet Version 1.0.18
        6. ncbi Blast 2.2.30+
'''

from __future__ import division
import sys, os, pysam, gzip, logging, subprocess, uuid, shutil
import operator, random
import logging.handlers
import time, datetime
import os.path
from subprocess import Popen, PIPE
from optparse import  OptionParser
from math import fabs
import numpy as np
from utility import  * 
from pysam import Fastafile
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from resistDBParser_v3c import CResistDB
_baseDir = '/well/gerton/hangphan/MMM/'

# Set up logging
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger('Log')
ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.setLevel(logging.DEBUG)
KPNE="R00000049"
ECOL="R00000042"
KOXY="R00000129"
speciesDicts={
    KPNE: "Klebsiella pneumoniae",
    KOXY: "Klebsiella oxytoxa",
    ECOL: "Escherichia coli"
    
}

class ResistType(object):
    def __init__(self, args):
        self.args = args
        self.refid = args[0]
        if args[1].endswith('.bam'):
            self.sampleid = args[1].split('/')[-1].split('_R000')[0]
            self.samplePath = args[1]
        else:
            self.sampleid = args[1]
            self.samplePath = 'data/_extensions/bam/{0}_{1}_v3.bam'.format(self.sampleid, self.refid)
        self.suffix=self.sampleid            
        if args[2] != None:
            self.plateName= args[2]
            self.suffix=self.plateName + "/" + self.sampleid 
        else:
            self.plateName = None
        self.isMetaGenomics = args[3]

        self.nSamples, self.sumLen = None, None  #estimate of how mixed the sample is (1 = pure culture, 2 = mixture of two subtypes, >2: metagenomic samples)
        self.meanCov, self.stdCov, self.meanDP = None, None, None  #Cov is coverage for contigs, DP is coverage depth for alignment bam
        self.outputDir = 'resistType/{0}/'.format(self.suffix)
        self.outputStats = 'resistType/{0}/sampleStats.txt'.format(self.suffix)
        self.velvetOutputDir= 'velvetOutput/{0}/'.format(self.suffix)
        self.spadesOutputDir= 'spadesOutput/{0}/'.format(self.suffix)
        self.contigFile = self.velvetOutputDir + 'contigs.fa'
        self.spadesContigFile=self.spadesOutputDir + "contigs.fasta"
        self.refFile = self.outputDir + '/reference.fa'
        self.bamFile = self.outputDir + '/output.bam'
        self.vcfFile = self.outputDir + '/output.vcf'
        self.blastResultFile = self.outputDir + '/blastResults.txt'
        self.outFastaFile = self.outputDir + '/outMatchedSequences.fa'
        self.outFastaFileRefined = self.outputDir + '/outRefinedMatchedSeqs.fa'
        self.tempFastqFile="/well/bag/grn/tmpDir/{0}/reads.fq.gz".format(self.suffix)
        if not os.path.exists(self.outputDir):
            os.mkdir(self.outputDir)
        if not os.path.exists(self.velvetOutputDir):
            os.mkdir(self.velvetOutputDir)
 

        #set resources 
        self.resistdb = _baseDir + "ResistDB/"
        self.resistGeneFastaFull= self.resistdb + 'ResistanceGeneSeqs_fw.fasta'
        self.resistGeneFastaFullPadded= self.resistdb + 'ResistanceGeneSeqs_fw_padded.fa'
        self.resistProteinFastaFull= self.resistdb + 'ResistanceGeneSeqs_aa.fasta'
        self.resistGeneFasta= self.resistdb + 'temp90.fa'
        self.resistGeneInfos = {}
        self.getResistGeneInfos()
        self.resistDBFile = "{0}/ResistDB/ResistanceProfiles.txt".format(_baseDir)

        #set link to tools used in the pipeline
        self.velvethPath = '/well/gerton/hangphan/MMM/bin/velveth'
        self.velvetgPath = '/well/gerton/hangphan/MMM/bin/velvetg'
        self.samtoolsPath = '/apps/well/samtools/0.1.19/bin/samtools'
        self.bedtoolsPath = '/apps/well/bedtools/2.19.1/bedtools'
        self.bwaPath = '/apps/well/bwa/0.7.10/bwa'
        self.makeblastdbProg='/well/gerton/hangphan/MMM/bin/makeblastdb'
        self.blastprog = '/well/gerton/hangphan/MMM/bin/blastn'
        
        #internal variables
        self.filteredGeneList= {}
        self.geneClusters = {}
        self.geneClusterMap = {}

        #resistGenes
        self.resistGenesMatch=[]
        self.resistGenesMismatch={}
        self.promoters=['ampC_promoter', 'P4', 'P5', 'Pa/Pb', 'Pc/Pd', 'P3', "ampCpromoter", 'Pa_Pb', 'Pc_Pd']
        self.chromGenes=['gyrA_ecol', 'gyrB_ecol', 'gyrA_kpne', 'gyrB_kpne', 'parC_ecol', 'parC_kpne', 'parE_ecol', 'parE_kpne']

        self.protDB=SeqIO.index(self.resistProteinFastaFull, "fasta")

        self.fastqFile1 =  "fastq/{0}/reads1.fq.gz".format(self.suffix)
        self.fastqFile2 =  "fastq/{0}/reads2.fq.gz".format(self.suffix)
        self.outputStatsFile = open(self.outputStats, "w")

    def runTest1(self):

        if self.isMetaGenomics:
            self.preFilterStep()
        else:
            if os.path.exists(self.spadesContigFile):
                self.contigFile = self.spadesContigFile

            if not os.path.exists(self.contigFile):
                readLen = self.getReadLength()
                self.runVelvet(readLen)
        self.runBlast()

        self.readResistGeneCluster()
        [listForExemptionFromAlignment, backToReference] = self.getPresentGenes()
        self.refFile  = self.makeTempRefFile(listForExemptionFromAlignment, backToReference) 
        if self.runBWA():
            with open(self.outFastaFile, "w") as f:
                pass
            self.getFastaSeqs()
            self.refineGeneSearch()

        self.makePredictions()
        self.outputStatsFile.close()
        logger.info("Finished resistType")

    def preFilterStep(self):
        '''
        filter reads by mapping to the big resistance gene file, then assemble using spades
        '''
        cmdLine = 'mkdir -p filteredFastq/{0} tmpDir/{0}'.format(self.suffix)
        os.system(cmdLine)
        
        cmdLine = '{0} mem {1} {2} {3} -O 9 -B 3  | python {4}/scripts/filterUnmappedReadpairs.py filteredFastq/{5}/reads1.fq.gz filteredFastq/{5}/reads2.fq.gz '.format(self.bwaPath, self.resistGeneFastaFullPadded, self.fastqFile1, self.fastqFile2, _baseDir, self.suffix)
        logger.info(cmdLine)
        os.system(cmdLine)
        self.spadesPath = 'spades.py'
        cmdLine = '{0} -1 filteredFastq/{1}/reads1.fq.gz -2 filteredFastq/{1}/reads2.fq.gz -o tmpDir/{1}   --careful -t 6  --phred-offset 33'.format(self.spadesPath, self.suffix)
        os.system(cmdLine)
        cmdLine = 'cp tmpDir/{0}/contigs.fasta resistType/{0}/contigs.fasta'.format(self.suffix)
        os.system(cmdLine)
        cmdLine = 'rm -rf tmpDir/{0}/'.format(self.suffix)
        os.system(cmdLine)

        self.contigFile = 'resistType/{0}/contigs.fasta'.format(self.suffix)
        self.fastqFile1 =  "filteredFastq/{0}/reads1.fq.gz".format(self.suffix)
        self.fastqFile2 =  "filteredFastq/{0}/reads2.fq.gz".format(self.suffix)
        


    def convert2Fastq(self):
        if not os.path.exists(self.fastqFile1):
            logger.info("Converting bam to fastq")
            unzipFile1 = self.fastqFile1.split(".gz")[0]
            unzipFile2 = self.fastqFile2.split(".gz")[0]
            commandLine = '{0} sort -n  {1} {2}/in'.format(self.samtoolsPath, self.samplePath, self.outputDir)
            logger.info(commandLine)
            os.system(commandLine)
            commandLine = '{0} bamtofastq -i {1}/in.bam -fq {2} -fq2 {3}'.format(self.bedtoolsPath, self.outputDir, unzipFile1, unzipFile2)
            logger.info(commandLine)
            os.system(commandLine)
            os.system("gzip  -f {0}".format(unzipFile1))
            os.system("gzip  -f {0}".format(unzipFile2))
            os.system("rm {0}/in.bam".format(self.outputDir))
        
    def inProtSeqDB(self, prot):
        for protid in self.protDB:
            if str(self.protDB[protid].seq) == str(prot).strip("*"):
                print "Matched protein sequence", protid
                return protid
        return None


    def makeTempRefFile(self, candidateList, addToRef):
        '''
        make temporary reference file which does not contain the ones in candidateList
        '''
        fullRecords =SeqIO.index(self.resistGeneFasta, "fasta")
        tempRefFile = self.outputDir + "/reference.fa"
        nSeq=""
        for i in range(100):
            nSeq+= "N"

        
        sequenceInRef = set()
        with open(tempRefFile, "w") as fo:
            nRegions = 7
            if self.nSamples ==1:
                contigs = SeqIO.index(self.contigFile, "fasta")
                #pick 7 regions randomly from contigs that are within the coverage range meanCov +- stdCov
                stdContigs = getSTDContigs(contigs, self.meanCov, self.stdCov)
                count = 0
                for i in range (nRegions):
                    randidx = random.randint(0, len(stdContigs)-1)
                    randContigID, randContigLen = stdContigs[randidx]
                    randContigStart = random.randint(0, randContigLen- 2001)
                    randSeq = str(contigs[randContigID].seq)[randContigStart: randContigStart + 2000]
                    randContigID = "REF_" + randContigID + "_" + str(randContigStart)
                    fo.write(">{0}\n{1}{2}{1}\n".format(randContigID, nSeq, randSeq, nSeq ))

            for record in fullRecords:
                if self.refid == KPNE and 'ecol' in record:
                    newRecord = record.split("_")[0] + "_kpne" 
                    sequenceInRef.add(newRecord)
                    continue
                if self.refid == ECOL and 'kpne' in record:
                    newRecord = record.split("_kpne")[0] + "_ecol" 
                    sequenceInRef.add(newRecord)
                    continue

                if (record not in candidateList and record not in sequenceInRef) and record not in sequenceInRef :
                    sequenceInRef.add(record)
                    fo.write(">{0}\n{2}{1}{2}\n".format(record, fullRecords[record].seq,nSeq))
            for gene in addToRef:
                if gene not in sequenceInRef:
                    sequenceInRef.add(gene)
                    fo.write(">{0}\n{2}{1}{2}\n".format(gene, fullRecords[gene].seq,nSeq))
                
        fullRecords.close()
        cmdLine = '{0} index {1}'.format(self.bwaPath, tempRefFile);      os.system(cmdLine)
        cmdLine = '{0} faidx {1}'.format(self.samtoolsPath, tempRefFile); os.system(cmdLine)

        return tempRefFile
        
        

    def runBlast(self):
        blastdb = self.contigFile.split(".fa")[0] 
        cmdLine = '{0} -in {1} -dbtype nucl -out {2}'.format(self.makeblastdbProg, self.contigFile, blastdb)
        #pipe result of this commandline as input to the next bit to extract the dna sequence, set up reference, run bwa-alignment, etc ...)
        p = Popen(cmdLine.split())
        p.wait()
        
        cmdLine = '{0} -gapopen 4 -gapextend 2 -evalue 0.001 -outfmt \"6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq\" -word_size 14 -db {1} -query {2} > {3} '.format(self.blastprog, blastdb, self.resistGeneFastaFull, self.blastResultFile)
        os.system(cmdLine)

    def getPresentGenes(self):
        resistGenes = SeqIO.index(self.resistGeneFastaFull, "fasta")
        genes_contigs = {}
        candidate_genes_presences = {}
        genes = {}
        if not os.path.exists(self.contigFile):
            logger.info("Failed to get an assembly of the samples. Move on to map based approach. ")
            return set()

        contigs = SeqIO.index(self.contigFile , "fasta")

        candidateContigs = set()
        exactContigs = set()
        fin = open(self.blastResultFile ,'r')
        logger.info("Candidates of gene presence with >92% matching in identity and >90% sequence length")
        EXACTMATCH, PARTIALMATCH, PROTMATCH= 2,0, 1#dna matching, partial matching,  protein matching
        contigs_genes = {}
        overlapGeneSet=set()
        for line in fin:
            qseqid, sseqid, pident, qlen, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qseq, sseq = line.split()
            length = int(length)
            qlen = int(qlen)
            qstart = int(qstart)
            qend = int(qend)
            sstart = int(sstart)
            send = int(send)
            pident = float(pident)
            gapopen = int(gapopen)
            geneLen = len(resistGenes[qseqid])
            sprot = ""
            qprot = ""

            if  pident < 90 or  qend - qstart + 1 < qlen * 0.9: #filter by length of alignment
                continue
            isTruncated = 0
            dnaMismatches =None
            newVector=[qseqid, pident, sseq, qseq, dnaMismatches, isTruncated, sprot, qprot, sstart, send]
            sVector=[sseqid, EXACTMATCH, pident, dnaMismatches, sstart, send, sprot, qseq, sseq]
            if qend - qstart +1 == qlen and pident == 100:
                if qseqid not in self.promoters:
                    sseq = Seq(str(sseq).replace("-",""), generic_dna)
                    sprot = sseq.translate(table="Bacterial")
                logger.info("{0}".format("\t".join(map(str, [qseqid, sseqid, length, qlen, geneLen, pident]))))
            else:
                #extend the hit to either end until either reaching the end of genes or end of contigs.
                offset1 = qstart -1
                offset2 = qlen - qend
                newSseq, newQseq = sseq, qseq
                newSstart, newSend = sstart, send
                if offset1 >0 or offset2 >0:
                    lensseq = len(str(contigs[sseqid].seq))
                    if sstart < send: #if hit match forwardly
                        newSstart = max(sstart - offset1 - 1, 0) #-1 because blastn return 1-index coordinate, in here we use 0-index coordinate
                        newSend = min(send + offset2-1, lensseq)
                        
                        isTruncated = newSstart - sstart + send - newSend + offset1 + offset2
                        leftSseq = str(contigs[sseqid].seq)[newSstart:sstart-1]
                        rightSseq = str(contigs[sseqid].seq)[send:newSend]
                        newSseq = leftSseq + sseq + rightSseq

                    if sstart > send: #if hit match reverse complement
                        newSstart = min(sstart + offset1-1, lensseq) 
                        newSend = max(send - offset2-1, 0)
                        isTruncated= sstart - newSstart  + newSend - send + offset1 + offset2
                        newSseq = contigs[sseqid].reverse_complement().seq[lensseq-newSstart: lensseq-newSend+1]
                        if fabs(sstart - newSstart)>=1:
                            leftSseq = newSseq[:int(fabs(sstart - newSstart))-1]
                        else: leftSseq = ""
                        if int(send- newSend) !=0:
                            rightSseq = newSseq[- int(fabs(send - newSend)) :]
                        else:
                            rightSseq = ""
                        newSseq = leftSseq + sseq + rightSseq

                    leftQseq = str(resistGenes[qseqid].seq)[max(qstart - int(fabs(sstart - newSstart)), 0):offset1]
                    rightQseq = str(resistGenes[qseqid].seq)[qend:min(qend + int(fabs(newSend - send)), qlen)]
                    newQseq = leftQseq + qseq + rightQseq
                if qseqid not in self.promoters and isTruncated ==0:
                    if len(str(newQseq).replace("-", "")) %3 ==0 and len(str(newSseq).replace("-", "")) %3 ==0:
                        sseq = Seq(str(newSseq).replace("-",""), generic_dna)
                        qseq = Seq(str(newQseq).replace("-",""), generic_dna)
                        sprot = sseq.translate(table="Bacterial")
                        qprot = qseq.translate(table="Bacterial")

                lenRatio=float(len(str(newSseq).replace("-", "")))/qlen
                pident = pident * lenRatio
                dnaMismatches = self.getMismatches(qseqid, newQseq, newSseq, mode=0, isDNA=1)
                newVector=[qseqid, pident, dnaMismatches, newSseq, newQseq, isTruncated, sprot, qprot, newSstart, newSend]

                if len(dnaMismatches.keys()) - int(mismatch) - offset1 - offset2  > 100:
                    logger.error("Something wrong, too much diffs: {0}".format(line))
                    print len(dnaMismatches.keys()), mismatch
                    print qseqid, len(newQseq), len(newSseq)
                    print newQseq[:100]
                    print newSseq[:100]
                    print newQseq[-100:]
                    print newSseq[-100:]
                if len(sprot) != 0 :
                    if str(sprot) != str(qprot):
                        sVector=[sseqid, PARTIALMATCH, pident, dnaMismatches, newSstart, newSend, sprot, newSseq, newQseq]
                    else:
                        newVector[1] = 100.0
                        sVector= [sseqid, PROTMATCH, pident, dnaMismatches,  newSstart, newSend, sprot, newSseq, newQseq]
                else:
                    sVector = [sseqid, PARTIALMATCH, pident, dnaMismatches, newSstart, newSend, sprot, newSseq, newQseq]

        
            if qseqid not in genes_contigs:
                genes_contigs[qseqid] = [sVector]
            else:
                genes_contigs[qseqid].append(sVector)
                
            if sseqid not in contigs_genes:
                contigs_genes[sseqid] = [newVector]
            else:
                #update partial matching of contigs to genes based on the overlap of gene matching to contigs.
                isOverlap = 0
                for idx, item in enumerate(contigs_genes[sseqid]):
                    overlapSize = min(max(send, sstart), max(item[-2], item[-1]))-max(min(sstart, send), min(item[-1],item[-2]))
                    if  overlapSize > qlen * 0.5 : #if is overlap
                        isOverlap =1
                        if newVector[1] > contigs_genes[sseqid][idx][1]:
                            overlapGeneSet.add(contigs_genes[sseqid][idx][0])
                            contigs_genes[sseqid][idx] = newVector
                        elif newVector[1] == contigs_genes[sseqid][idx][1] and qseqid not in contigs_genes[sseqid][idx]:
                            contigs_genes[sseqid][idx].extend(newVector)
                        else:
                            overlapGeneSet.add(qseqid)

                if isOverlap ==0:
                    contigs_genes[sseqid].append(newVector)
            continue

        fin.close()
        candidateList = set()
        fo =  open(self.outFastaFileRefined, "w")
        #write matched sequences to file
        for gene in sorted(genes_contigs.keys()):
            
            items = genes_contigs[gene] #contain vector of [contigName, matchStatus, pident, dnaMismatches, sstart, send, sprot, qseq, sseq]
            for item in items:
                contigName, matchStatus, pident, dnaMismatches, sstart, send, sprot, qseq, sseq = item
                if matchStatus == EXACTMATCH:
                    fo.write(">{0}\texactMatch:mappedContig:{1}\n{2}\n".format(gene, ":".join(map(str, item[:-2]) ), sseq))
                    self.resistGenesMatch.append([gene, "dnaMatch",pident, dnaMismatches,None, [contigName, sstart, send]])
                    candidateList.add(self.geneClusterMap[gene])
                if matchStatus == PROTMATCH:
                    duplicate = 0
                    for thisGene in genes_contigs.keys():
                        if thisGene != gene and genes_contigs[thisGene][0][5] == send and genes_contigs[thisGene][0][1] == EXACTMATCH:
                            duplicate = 1
                    if duplicate == 0:
                        fo.write(">{0}\tproteinLevelMatch:mappedContig:{1}\n{2}\n".format(gene, ":".join(map(str, item[:-2]) ),sseq))
                        self.resistGenesMatch.append([gene, "protMatch", pident, dnaMismatches, None, [contigName, sstart, send]])
                        candidateList.add(self.geneClusterMap[gene])

        #write close match sequences to file, update snps information along the way
        backToReference=[]
        for contig in contigs_genes:

            for hit in contigs_genes[contig]:
                [gene, pident, dnaMismatches, sseq, qseq] = hit[:5]
                isTruncated = hit[5]
                if isTruncated > 0:
                    backToReference.append(self.geneClusterMap[gene]) #for mapping purpose due to the ambiguity in assembly approach
                    continue
                if pident >=100:#this is already covered by perfect match genes
                    continue
                if pident < 90:# don't consider if sequence identity is low 
                    continue
                logger.info("{0}".format("\t".join(map(str, [gene, contig, len(sseq), pident]))))
                #newVector=[qseqid, pident, dnaMismatches, newSseq, newQseq, sprot, qprot, newSstart, newSend]
                items = genes_contigs[gene]
                isTruncatedProtein=0
                contigName, sstart, send = item[0], item[4], item[5]
                proteinMismatches=None
                for item in items:
                    if gene in self.promoters or "-" in sseq or "-" in qseq: #if is promoter of there is gap in the alignment, indication of frameshift mutations, no point to translate and compare
                        fo.write(">{0}\tcloseHit:mappedContig:{1}\t{2}\n{3}\n".format(gene, ":".join(map(str, item[:-2]) ), item[-1], item[-2]))   
                        if gene not in self.promoters:
                            sprot = str(Seq(str(sseq).replace("-", "")).translate(table="Bacterial"))
                            #print sprot
                            if "*" in sprot:
                                if sprot.index('*') < len(sprot)-1:
                                    isTruncatedProtein = 1
                            else:
                                pass
                        self.resistGenesMismatch[gene]= [pident, dnaMismatches, proteinMismatches, isTruncatedProtein , None, [contigName, sstart, send]]
                        candidateList.add(self.geneClusterMap[gene])
                    else:
                        sseq = Seq(str(sseq).replace("-", ""), generic_dna)
                        sprot = sseq.translate(table="Bacterial")
                        isInProtSeq = self.inProtSeqDB(sprot)
                        sprot = str(sprot)
                        if "*" in sprot:
                            if sprot.index('*')< len(sprot)-1:
                                isTruncatedProtein = 1
                        if isInProtSeq != None:
                            self.resistGenesMatch.append([isInProtSeq, "protMatch", None, None, None, [contigName, sstart, send]]) 
                            fo.write(">{0}\tproteinLevelMatch1:mappedContig:{1}\t{2}\n{3}\n".format(isInProtSeq, ":".join(map(str, item[:-2]) ), item[-1], item[-2]))  
                            continue
                        proteinMismatches= self.getMismatches(gene, qseq, sseq, mode=0, isDNA=0)
                        if dnaMismatches == None:
                            self.resistGenesMatch.append([gene, "dnaMatch", 100.0, None, None, [contigName, sstart, send]])
                            continue
                        if proteinMismatches==None or len(proteinMismatches.keys()) == 0:
                            self.resistGenesMatch.append([gene, "protMatch", pident, dnaMismatches, None, [contigName, sstart, send]])
                            fo.write(">{0}\tproteinLevelMatch:mappedContig:{1}\t{2}\n{3}\n".format(gene, ":".join(map(str, item[:-2]) ), item[-1], item[-2]))  
                        else:
                            self.resistGenesMismatch[gene] = [pident, dnaMismatches, proteinMismatches, isTruncatedProtein, None, [contigName, sstart, send]]
                            fo.write(">{0}\tcloseHit:mappedContig:{1}\t{2}\n{3}\n".format(gene, ":".join(map(str, item[:-2]) ), item[-1], item[-2]))  
                candidateList.add(self.geneClusterMap[gene])
        
        self.estimateCopyNumber()
        return [set(list(candidateList) + list(overlapGeneSet)), backToReference]

    def estimateCopyNumber(self):
        contigs = SeqIO.index(self.contigFile, "fasta")
        self.nSamples, self.sumLen = checkAssemblySize(self.contigFile)
        self.outputStatsFile.write("Estimated number of samples: {0}\n".format(self.nSamples))
        self.outputStatsFile.write("Assembly size: {0}\n".format(self.sumLen))

        if self.nSamples != 1:
            logger.info("Not estimating copy number of resistance genes because the sample is not a pure culture sample")
            return 

        self.meanCov, self.stdCov, rateLen  = getAssemblyMeanCoverage(self.contigFile)  #rateLen is the ratio of sumLen of contigs calculated in meanCov over sumLen of all contigs
        self.outputStatsFile.write("Mean contig coverage: {0}\n".format(self.meanCov))
        self.outputStatsFile.write("Std contig coverage: {0}\n".format(self.stdCov))
        for idx, item in enumerate(self.resistGenesMatch):
            contigName=item[-1][0]
            self.resistGenesMatch[idx][4] = estimateContigCopyNumber(contigName, contigs, self.meanCov,self.stdCov)
            self.resistGenesMatch[idx][-1].extend([self.meanCov, self.stdCov])
        for gene in self.resistGenesMismatch:
            contigName= self.resistGenesMismatch[gene][-1][0]
            self.resistGenesMismatch[gene][-2] = estimateContigCopyNumber(contigName, contigs, self.meanCov, self.stdCov)
            self.resistGenesMismatch[gene][-1].extend([self.meanCov, self.stdCov])
        return


    def readResistGeneCluster(self):
        clusterFile = self.resistGeneFasta + ".clstr"

        genes, repGene = [], None
        with open(clusterFile, "r") as f:
            for line in f:
                if line.startswith(">"):
                    if len(genes) > 0:
                        self.geneClusters[repGene] = genes
                        genes = []
                
                else:
                    geneName = line.strip().split()[2].rstrip(".")[1:]
                    genes.append(geneName)
                    if line.strip().endswith("*"):
                        repGene = geneName

        self.geneClusters[repGene] = genes
        
        for cluster in self.geneClusters:
            for gene in self.geneClusters[cluster]:
                self.geneClusterMap[gene] = cluster


    def runBWA(self):
        if not os.path.exists(self.fastqFile1):
            return 0
        logger.info("Run bwa to align to candidate genes")
        logger.info("Start BWA mapping")
        cmdLine = '{0} mem {1} {2} {3} -O 9 -B 3  | {4} view -F 4 -Shu - |{4} sort - {5}'.format(self.bwaPath, self.refFile, self.fastqFile1, self.fastqFile2, self.samtoolsPath, self.bamFile.split(".bam")[0])
        logger.info(cmdLine)
        os.system(cmdLine)
        cmdLine = '{0} index {1}'.format(self.samtoolsPath, self.bamFile)
        os.system(cmdLine)
        logger.info("Start samtools mpileup ")
        cmdLine = 'samtools mpileup -f {0} -E -M0 -q10   -F0.002 -D -g -S  {1}|bcftools view -c -g -b -A -L -t0.01 -i-1 -p0.5 -Pfull -|bcftools view - >{2}'.format(self.refFile, self.bamFile, self.vcfFile)
        cmdLine = 'samtools mpileup -f {0} -E -A  -D -g -S  {1}|bcftools view -c -g -b -A -L -t0.01 -i-1 -p0.5 -Pfull -|bcftools view - |sed \'s/,X//\'>{2}'.format(self.refFile, self.bamFile, self.vcfFile)
        os.system(cmdLine)
        return 1


    #####################################################
    def getFastaSeqs(self):
        '''
        From the vcf file, make this into sequences of genes, with info on number of low coverage bases, number of Heterozygote calls (to make inference about multiple gene copies) 
        and the average depth 
        '''
        gatk="/apps/well/gatk/3.5-0/GenomeAnalysisTK.jar"
        cmdLine = "rm {0}/reference.dict".format(self.outputDir)
        os.system(cmdLine)
        cmdLine = "java -jar /gpfs0/apps/well/picard-tools/1.111/CreateSequenceDictionary.jar R={0}/reference.fa O={0}/reference.dict".format(self.outputDir)
        os.system(cmdLine)
        cmdLine= "bedtools genomecov -ibam {0}/output.bam -dz > {0}/outputBam.cov".format(self.outputDir)
        os.system(cmdLine)
        
        records = SeqIO.index("{0}/reference.fa".format(self.outputDir), "fasta")
        covInfo = {}
        lenInfo = {}
        nLowDP = {}
        for record in records:
            covInfo[record] = []
            lenInfo[record] = len(str(records[record].seq))
            nLowDP[record] = 0

        with open("{0}/outputBam.cov".format(self.outputDir), "r") as f:
            for line in f:
                cols = line.strip().split()
                coord = int(cols[1])
                cov = int(cols[2])
                if coord < 100 or coord > lenInfo[cols[0]] - 100:
                    continue
                covInfo[cols[0]].append(cov)
                if cov <2:
                    nLowDP[record] +=1
        RefDP = []
        for record in records:
            if record.startswith("REF"):
                avgCov = float(np.sum(np.array(covInfo[record]))) / (lenInfo[record] - 200)
                RefDP .append(avgCov)
            else:
                nLowDP[record] += lenInfo[record]  -200 - len(covInfo[record])
        self.meanDP = np.average(np.array(RefDP)) 
        self.outputStatsFile.write("Mean depth from BAM: {0}\n".format(self.meanDP))
        with  open("{0}/gatkIntervals.list".format(self.outputDir), "w")       as f:
            for record in records:
                if record.startswith("REF"):
                    continue
                if float(nLowDP[record]) > (lenInfo[record] -200) * 0.1:
                    continue
                f.write("{0}\n".format(record))
                self.filteredGeneList[record] = "all"
                

        cmdLine = "java -jar {0} -T FastaAlternateReferenceMaker     -R {1}/reference.fa     -o   {1}/output.fasta   -V {1}/output.vcf     -IUPAC  {1}/output.bam -L {1}/gatkIntervals.list ".format(gatk, self.outputDir)
        os.system(cmdLine)
        fo = open("{0}/output1.fasta".format(self.outputDir), "w") 
        records = SeqIO.parse("{0}/output.fasta".format(self.outputDir), "fasta")  
        for record in records:
            geneName= record.description.split()[1].split(":")[0]
            lowDP = max(nLowDP[geneName], 0)
            seq = str(record.seq)
            nHet = 0
            for thisChar in seq:
                if thisChar not in ['A', 'T', 'G', 'C', 'N']:
                    nHet +=1
            
            meanDP = float(np.sum(np.array(covInfo[geneName])))/(lenInfo[geneName] - 200)
            fo.write(">{0}\tnLowDP_{1}:nHet_{2}:meanDP_{3}\n{4}\n".format(geneName, lowDP, nHet, meanDP, seq))

        return

    #####################################################  
    def makePredictions(self):
        '''
        Given a list of genes, match these with the database and make predictions of antibiotic resistance
        
        # self.resistGenesMatch #list of genes with exact match to the alleles
        # self.resistGenesMismatch #list of genes with partial matchs (have mutation) - only care about those that are chromosomal ,e.g parE, parC, gyrA, gyrB etc. 
        '''

        ResistDB= CResistDB(self.resistDBFile)
        species = {'R00000042':'Ecol', 'R00000049': 'Kpne'}
        thisSpecies='none'
        if self.refid in species:
            thisSpecies = species[self.refid]
        self.outPredictionFile = '{0}/resistancePred.txt'.format(self.outputDir)
        self.outGenotypeFile='{0}/genotypePred.txt'.format(self.outputDir)
        ResistDB.resistancePrediction(self.resistGenesMatch, self.resistGenesMismatch, thisSpecies, self.outPredictionFile, self.sampleid, self.outGenotypeFile)
        


    #####################################################  
    def refineGeneSearch(self):

        hitGenes = []
        closestHits = []
        self.outFastaFile= self.outputDir + "/output1.fasta"
        records1 = SeqIO.index(self.outFastaFile, "fasta")
        if len(records1) ==0:
            return
        fo = open(self.outFastaFileRefined, "a")
        records=SeqIO.index(self.resistGeneFastaFull, "fasta")

        for gene in self.filteredGeneList.keys():
            if gene.startswith("REF"):
                continue
            fn =  "{0}/ResistDB/clusterFas/{1}.fa".format(_baseDir, gene)
            fn = fn.replace("(", "").replace(")", "").replace("\'", "")
            clusterRecords = SeqIO.index(fn, "fasta")
            thisDP = float(records1[gene].description.split()[-1].split(":")[-1].split("_")[-1]) #coverage depth
            copyNumber = None

            if self.meanDP:
                copyNumber = round(thisDP/self.meanDP)
            if os.path.exists(fn):
                outBlastFN = self.outputDir + "/tempBlast.txt"
                cmdLine = '{0} -gapopen 5 -gapextend 2 -evalue 0.001 -outfmt \"6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq\" -word_size 17 -subject {1} -query {2} >{3} '.format(self.blastprog,"{0}/output1.fasta".format(self.outputDir), fn , outBlastFN)
                os.system(cmdLine)
                with open(outBlastFN, "r") as f:
                    hasHit = 0
                    closestHit = []
                    maxIdent= 0
                    ssequence = ""                    
                    bestSeq = []
                    for line in f:

                        cols = line.strip().split()
                        pident = float(cols[2])
                        qseqid, sseqid = cols[0], cols[1]
                        qlen, length = int(cols[3]), int(cols[5])
                        if length < qlen*.9:
                            continue

                        qstart, qend, sstart, send = int(cols[8]), int(cols[9]), int(cols[10]), int(cols[11])
                        qseq = cols[-2]
                        sseq = cols[-1]
                        slen = int(cols[4])
                        if qlen != length:
                            continue
                        
                        dnaMismatches = self.getMismatches(qseqid, qseq, sseq, mode=0, isDNA=1)
                        proteinMismatches= self.getMismatches(qseqid, qseq, sseq, mode=0, isDNA=0)
                        isTruncatedProtein = 0
                        sprot = str(Seq(sseq.replace("-","")).translate(table="Bacterial"))
                        if "*" in sprot:
                            if sprot.index("*") < len(sprot) -1:
                                isTruncatedProtein =1


                        if len(dnaMismatches) ==0 or not dnaMismatches:
                            hasHit =1
                            hitGenes.append(cols[0])
                            if cols[0] not in self.resistGenesMatch:
                                self.resistGenesMatch.append([qseqid, "dnaMatch", 100.0, dnaMismatches, copyNumber, ['mapBased', str(thisDP), str(self.meanDP)]])
                            fo.write(">{0}\texactHit:{1}\n{2}\n".format(cols[0], records1[gene].description, ssequence))
                            break
                        else:
                            if len(proteinMismatches) ==0:
                                hasHit =1
                                hitGenes.append(cols[0])
                                if cols[0] not in self.resistGenesMatch:
                                    self.resistGenesMatch.append([qseqid, "protMatch", pident, dnaMismatches, copyNumber, ['mapBased', str(thisDP), str(self.meanDP)]])
                                    fo.write(">{0}\tproteinMatched1:{1}\n{2}\n".format(cols[0], records1[gene].description, ssequence))
                                break
                        
                            tempVal = float(qlen - len(dnaMismatches)) / qlen 
                            if tempVal > maxIdent and (maxIdent ==0 or length > len(bestSeq[0][0])):
                                maxIdent = tempVal
                                closestHit = [qseqid]
                                bestSeq = [[qseq, sseq, dnaMismatches, proteinMismatches,  isTruncatedProtein]]
                                
                            elif tempVal == maxIdent:
                                closestHit.append(qseqid)
                                bestSeq.append([qseq, sseq, dnaMismatches, proteinMismatches, isTruncatedProtein])
                    if hasHit==0 and len(bestSeq) >0:
                        nN = bestSeq[0][1].count('N')
                        if float(nN)/len(bestSeq[0][1]) >0.05:
                            continue
                        fo.write(">{0}\tcloseHit:{1}\t{2}\n{3}\n".format("|".join(closestHit), records1[gene].description, bestSeq[0][0], bestSeq[0][1]))
                        #sys.stdout.write(">{0}\tcloseHit:{1}\t{2}\n{3}\n".format("|".join(closestHit), records1[gene].description, len(bestSeq[0][0]), len(bestSeq[0][1])))
                        dnaMismatches=bestSeq[0][2]
                        proteinMismatches = bestSeq[0][3]
                        isTruncatedProtein=bestSeq[0][4]
                        if dnaMismatches!=None:
                            self.resistGenesMismatch[closestHit[0]] = [maxIdent*100, dnaMismatches, proteinMismatches, isTruncatedProtein, copyNumber, ['mapBased', str(thisDP), str(self.meanDP)]]
                            
                os.system("rm {0}".format(outBlastFN))
            else:
                fo.write(">{0}\tcloseHit:{1}\n{2}\n".format(gene, records1[gene].description, records1[gene].seq))
        records.close()
        records1.close()
        return

    def getMismatches(self, seqName, qseq, sseq, mode =0, isDNA=0): 
        """
        Given two sequences, try to get the mismatches between them. Translate the sequence if needed.
        """
        qseq = str(qseq)
        sseq = str(sseq)

        mismatches = {}
        for i in range(min(len(sseq), len(qseq))):
            if sseq[i] != qseq[i]:
                mismatches[i] = [qseq[i], sseq[i]]
        if seqName in self.promoters or "-" in sseq + qseq or isDNA:
            return mismatches
        needReverse = 0
        if len(mismatches.keys()) > len(sseq) * 0.4 and len(sseq) < 1000 and mode ==1:
            alignments = pairwise2.align.globalxx(qseq, sseq)[0]
            newqseq, newsseq, nmatch = alignments[0:3]
            if nmatch > len(sseq)*0.8:
                qseq, sseq = newqseq, newsseq
                mismatches = {}
                for i in range(min(len(sseq), len(qseq))):
                    if sseq[i] != qseq[i]:
                        mismatches[i] = [qseq[i], sseq[i]]
            else:
                needReverse = 1
        if  len(mismatches.keys()) > len(sseq) * 0.4 and (mode ==0 or len(sseq) >1000):
            needReverse =1
            mismatches = {}
            sseq = Seq(str(sseq), generic_dna).reverse_complement()
            for i in range(min(len(sseq), len(qseq))):
                if sseq[i] != qseq[i]:
                    mismatches[i] = [qseq[i], sseq[i]]
        if len(mismatches.keys()) > len(sseq) * 0.4:
            logger.debug("{0} sequences not matching".format(seqName))
            logger.debug("{0}\n{1}\n".format(qseq, sseq))
            return None
        if seqName not in self.promoters:
            sseq = Seq(str(sseq).replace("-", ""), generic_dna)
            qseq = Seq(str(qseq).replace("-", ""), generic_dna)

            sprot = sseq.translate(table="Bacterial")
            qprot = qseq.translate(table="Bacterial")
            countStar = str(qprot).count("*")

            if countStar >1:#check reverse complement
                qprot = sseq.reverse_complement().translate(table="Bacterial")
                sprot = qseq.reverse_complement().translate(table="Bacterial")
            mismatches={}
            for i in range(min(len(qprot), len(sprot))):
                if qprot[i] != sprot[i]:
                    mismatches[i] = [qprot[i], sprot[i]]

            return mismatches
        else:
            return mismatches

    #####################################################  
    def getReadLength(self, readnum=20):
        """
        Calculate read length based on first N reads that passed QC.
        """

        if os.path.exists(self.samplePath):
            # extract first 'readnum' reads that have not failed QC and are not PCR duplicates
            cmd1 = '{samtools} view -F 1536 {inbampath}'.format(samtools=self.samtoolsPath, inbampath=self.samplePath)    
            cmd2 = 'head -n {0}'.format(readnum)
            cmd3 = 'cut -f10'
            logger.info('Estimating read length based on first {0} high-quality reads in {1}:'.format(readnum, self.samplePath))
            logger.info(' | '.join((cmd1, cmd2, cmd3)))
            p1 = Popen(cmd1.split(), stdout=PIPE, stderr=PIPE)
            p2 = Popen(cmd2.split(), stdout=PIPE, stderr=PIPE)
            p3 = Popen(cmd3.split(), stdin=p1.stdout, stdout=PIPE, stderr=PIPE)
            p1.stdout.close() # mimic pipe
            p2.stdout.close() 
            stdout, stderr = p3.communicate()    
            if not stdout:
                logger.critical('Failed to get reads from {b}: {e}'.format(b = self.samplePath, e = stderr if stderr else 'None')); sys.exit(1)
            try:
                readLength = np.mean([len(x) for x in stdout.split()])
            except:
                logger.critical('Failed to estimate read length : {e}'.format(e = stderr if stderr else 'None')); sys.exit(1)
            logger.info('Read length is {0}'.format(readLength))
            return readLength
        elif os.path.exists(self.tempFastqFile):
            f=gzip.open(self.tempFastqFile)
            counter = 0
            readLen = []
            for line in f:
                counter +=1
                if counter%4 ==2:
                    readLen.append(len(line.strip()))

                if counter > readnum *4:
                    break
            f.close()
            return float(sum(readLen)) / len(readLen)
        return 151
    def runVelvet(self, readLength):

        tmpdir = 'tmpDir/{0}.velvet'.format(self.sampleid)
        if not os.path.isdir(tmpdir):
            os.mkdir(tmpdir)
        self.tempFastqFile = "fastq/{0}/reads.fq.gz".format(self.suffix)
        if not(os.path.exists(self.tempFastqFile)):
            cmdLine = "python {0}/scripts/interleaveFastq.py {1} {2} {3}".format(_baseDir, self.fastqFile1, self.fastqFile2, self.tempFastqFile)
            os.system(cmdLine)

        k = int(0.67 * readLength)
        k = k - (1-k%2) # make k odd
        logger.info('Starting Velvet on {0} with k={1}'.format(self.samplePath, k))
        cmd1 = '{vh} {tmpdir} {k} -shortPaired  -fastq.gz {infq1} '.format(vh=self.velvethPath, k=k, tmpdir=tmpdir, infq1=self.tempFastqFile)
        logger.info(cmd1)
        p1 = Popen(cmd1.split(), stdout=PIPE, stderr=PIPE)
        stdout, stderr = p1.communicate()    
        if stderr:
            logger.critical('Error while trying to run Velvet: {e}'.format(e = stderr if stderr else 'None')); sys.exit(1)
        cmd2 = '{vg} {tmpdir} -exp_cov auto -cov_cutoff auto -very_clean yes'.format(vg=self.velvetgPath, tmpdir=tmpdir)
        logger.info(cmd2)
        p2 = Popen(cmd2.split(), stdout=PIPE, stderr=PIPE)
        stdout, stderr = p2.communicate()    
        
        tmp_outpath = os.path.join(tmpdir, 'contigs.fa')
        if not os.path.isfile(tmp_outpath):
            logger.critical('Error while trying to run Velvet: {e}'.format(e = stderr if stderr else 'None')); sys.exit(1)
    
        logger.info('Saving contigs: {0}'.format(self.contigFile))
        shutil.move(tmp_outpath, self.contigFile)
        shutil.move(os.path.join(tmpdir, 'Log'), self.velvetOutputDir+ '/Log')
        shutil.rmtree(tmpdir)
        return

    def getResistGeneInfos(self):
        '''
        read information about location (chromosomal or not) and species (to filter genes from presence/absence list and decide to make annotation or not (aa changes)
        '''
        chromInfoFile = self.resistdb + '/chromosomalGenes.txt'
        speciesInfoFile = self.resistdb + '/speciesSpecificGenes.txt'
        
        with open(chromInfoFile, 'r') as f:
            for line in f:
                gene, location = line.strip().split("\t")[:2]
                self.resistGeneInfos[gene] = [location, []]
                
        with open(speciesInfoFile, 'r') as f:
            for line in f:
                gene, species = line.strip().split("\t")[:2]
                species = species.split(",")
                if gene in self.resistGeneInfos:
                    self.resistGeneInfos[gene][1].extend(species)
                else:
                    self.resistGeneInfos[gene] = ["N/A", species]
        return 
######################################################    

    
if __name__ == "__main__":
    usage = "usage: python %prog [options] \n" 
    version = "%prog 0.1"
    
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--refid", dest="refid", type="string", default=None, help="Reference identification")
    parser.add_option("-s", "--sampleName", dest="sampleName", type="string", default=None, help="Sample name")
    parser.add_option("-p", "--plateName", dest="plateName", type="string", default=None, help="Plate name")
    parser.add_option("-m", "--metaGenomics", dest = "metaGenomics", type="int", default = 0, help="Sample is metagenomics data. Will filter reads before assembly.")
    (opts, args) = parser.parse_args()
    
    resistTypeModule= ResistType([opts.refid, opts.sampleName, opts.plateName, opts.metaGenomics])
    resistTypeModule.runTest1()
