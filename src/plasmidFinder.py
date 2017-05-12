from __future__ import division
import sys, os, pysam, gzip, logging, subprocess, uuid, shutil
import operator
import logging.handlers
import time, datetime
import os.path
from subprocess import Popen, PIPE
from optparse import  OptionParser
import numpy as np
from pysam import Fastafile
from Bio import SeqIO
_baseDir="/".join(os.path.dirname(os.path.realpath(sys.argv[0])).split("/")[:-1])
# Set up logging
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger('Log')
ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.setLevel(logging.DEBUG)

class PlasmidFinder(object):
    def __init__(self, args):
        self.args = args
        self.sampleid = args[0]
        self.plasmidfinderDir= _baseDir + "/resources/plasmidfinder/"
        self.alleleDB= self.plasmidfinderDir + "plasmid_database.fsa"
        self.contigFile = args[1]

        self.outDir = "plasmidfinderOutput/{0}/".format(self.sampleid)
        cmdLine = "rm -rf "+ self.outDir
        os.system(cmdLine)
        cmdLine = "mkdir -p " + self.outDir
        os.system(cmdLine)
        self.blastFile = self.outDir + "outputBlast.txt"
        self.outFile = self.outDir + "outPlasmidFinder.txt"
        self.closestMatchedAlleles=[]
        self.matchedAlleles=[]
        
    def runPlasmidFinder(self):
        self.runBlast()
        self.getMatchedAlleles()
        return

    def runBlast(self):
        '''
        run blast of contig file against the database of all alleles 
        '''
        cmdLine = "blastn -query {0} -db {1} -word_size 17  -evalue 0.001  -gapopen 5 -gapextend 2 -culling_limit 1 -outfmt '6 qseqid sseqid pident length slen qstart qend sstart send qseq sseq'|awk '$3>80 && $4/$5>0.90' > {2}".format(self.contigFile, self.alleleDB, self.blastFile)
        print cmdLine
        os.system(cmdLine)
        return


    def getMatchedAlleles(self):
        '''
        read blast result from blastFile and get matched allele or it's closest hit
        '''
        self.matchedAlleles = []
        self.closestMatchedAlleles = []
        contig_allele = {}
        for line in open(self.blastFile):
            cols = line.strip().split()
            qseqid, sseqid = cols[:2]
            pident, length, slen, qstart, qend, sstart, send = map(float, cols[2:9])
            qseq, sseq = cols[9:11]

            if send < sstart:
                temp = send
                send = sstart
                sstart = temp
            lenRatio = length / (send - sstart +1)
            if lenRatio <0:  lenRatio = -lenRatio
            if lenRatio >1: lenRatio = 1/lenRatio
            pident = pident  * lenRatio
            newVector = [sseqid, pident, qstart, qend, qseq, sseq]
            if qseqid not in contig_allele:
                contig_allele[qseqid] = [newVector]
            else:
                #check overlap
                isOverlap = 0
                for idx, item in enumerate(contig_allele[qseqid]):
                    overlapSize = min(qend, item[3]) - max(qstart, item[2])
                    if overlapSize > slen *0.8 :
                        isOverlap = 1
                        if newVector[1] > item[1]:
                            contig_allele[qseqid][idx] = newVector
                if isOverlap ==0:
                    contig_allele[qseqid].append(newVector)
            
        f = open(self.outFile, 'w')
        for item in contig_allele.keys():
            matches = contig_allele[item]
            for match in matches:
                f.write("{0},{1},{2},{3},{4},{5}\n".format(self.sampleid, item, match[0], match[1], int(match[2]),int(match[3])))
                self.matchedAlleles.append(match)
        f.close()
        return


if __name__ == "__main__":
    usage = "usage: python %prog [options] \n" 
    version = "%prog 0.1"
    
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-s", "--sampleName", dest="sampleName", type="string", default=None, help="Sample name")
    parser.add_option("-c", "--contigFile", dest="contigFile", type="string", default=None, help="link to bam file")
    (opts, args) = parser.parse_args()
    
    PlasmidFinderModule= PlasmidFinder([opts.sampleName, opts.contigFile])
    PlasmidFinderModule.runPlasmidFinder()
