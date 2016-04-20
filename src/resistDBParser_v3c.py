import sys
import os
# 0 gene_name
# 1 notes
# 2 alternative_name_1
# 3 alternative_name_2
# 4 alternative_name_3
# 5 alternative_name_4
# 6 location
# 7 accession_number_1
# 8 start_pos
# 9 end_pos
# 10 length
# 11 sequence
# 12 rev_comp
# 13 is_rev_complement
# 14 aa_sequence
# 15 PubmedID
# 16 reference_author
# 17 reference_title
# 18 reference_publication
# 19 pub_host_organism
# 20 summary_resistance_class
# 21 sub_resistance_class
# 22 primers_forward
# 23 primers_reverse
# 24 primers_seq_forward
# 25 primers_seq_reverse
# 26 primers_universal_forward
# 27 primers_universal_reverse
# 28 primer_references
# 29 ambler_class
# 30 bush_jacoby_class
# 31 OXA_family
# 32 phenicol_resistance_group
# 33 penicillin
# 34 amoxicillin_ampicillin
# 35 cephalexin
# 36 cefuroxime
# 37 cefoxitin
# 38 ceftriaxone
# 39 ceftazidime
# 40 cefepime
# 41 co-amoxiclav
# 42 piperacillin-tazobactam
# 43 carbapenems
# 44 aztreonam
# 45 carbenicillin
# 46 oxacillin
# 47 tetracycline
# 48 minocycline
# 49 doxycycline
# 50 tigecycline
# 51 arbekacin
# 52 amikacin
# 53 astromycin_fortimicin
# 54 butirosin
# 55 gentamicin
# 56 hygromycin_B
# 57 isepamicin
# 58 sisomicin
# 59 5-epi-sisomicin
# 60 dibekacin
# 61 netilmicin
# 62 2'-N-ethylnetilmicin
# 63 6'-N-ethylnetilmicin
# 64 tobramycin
# 65 kanamycin
# 66 neomycin
# 67 lividomycin
# 68 paromomycin
# 69 ribostamycin
# 70 spectinomycin
# 71 streptomycin
# 72 fosfomycin
# 73 chloramphenicol
# 74 nalidixic acid
# 75 norfloxacin
# 76 ofloxacin
# 77 ciprofloxacin
# 78 levofloxacin
# 79 moxifloxacin
# 80 trimethoprim
# 81 sulfamethoxazole
# 82 co-trimoxazole
# 83 nitrofurantoin
# 84 notes_on_phenotype
# 85 prediction_confidence_level
# 86 additional_ref_authors
# 87 additional_ref_1
# 88 additional_ref_2
# 89 audit_notes

class geneFeatures(object):
    def __init__(self, line, headerLine):
        cols = line.strip().split("\t")
# 15 PubmedID
# 16 reference_author
# 17 reference_title
# 18 reference_publication
# 19 pub_host_organism
# 20 summary_resistance_class
# 21 sub_resistance_class
# 22 primers_forward
# 23 primers_reverse
# 24 primers_seq_forward
# 25 primers_seq_reverse
# 26 primers_universal_forward
# 27 primers_universal_reverse
# 28 primer_references
# 29 ambler_class
# 30 bush_jacoby_class
        self.sequence = cols[12].strip()
        self.name = cols[0].strip().replace("(", "").replace(")", "").replace("\'", "_prime").replace("/", "_").replace(",", "_").replace(" ", "")
        self.refAuthor=None
        self.refTitle=None
        self.refPub=None
        self.pubHostOrg=None
        self.summaryResClass=None
        self.subResClass=None
        self.primersForw=None
        self.primersRev=None
        self.primerSeqForw=None
        self.primerSeqRev=None
        self.primersUniForw=None
        self.primersUniRev=None
        self.primerRef=None
        self.amblerClass=None
        self.bushJacobyClass=None
        self.OXAFam=None
        self.resistanceProfile={}
        self.specialTreatment=False
        if len(cols) >=21:
            self.summaryResClass= cols[20]
            self.specialTreatment=True 
        else:
            pass  #chromosomal genes - gyrA, gyrB, parC, parE
            #print self.name, len(cols)
        if len(cols) >=32:
            self.pubmedID=cols[15]
            self.refAuthor=cols[16]
            self.refTitle=cols[17]
            self.refPub=cols[18]
            self.pubHostOrg=cols[19]
            self.summaryResClass=cols[20]

            self.subResClass=cols[21]
            self.primersForw=cols[22]
            self.primersRev=cols[23]
            self.primerSeqForw=cols[24]
            self.primerSeqRev=cols[25]
            self.primersUniForw=cols[26]
            self.primersUniRev=cols[27]
            self.primerRef=cols[28]
            self.amblerClass=cols[29]
            self.bushJacobyClass=cols[30]
            self.OXAFam=cols[31]
            headerKeys = headerLine.strip().split("\t")
            self.resistanceProfile={}
            for i in range(33,min(len(cols), 84)):
                if cols[i].strip() != "":
                    self.resistanceProfile[headerKeys[i]]= cols[i].strip()
class CResistDB(object):
    def __init__(self, fileName):
        self.resistanceGenes = {}
        self.headerLine = None
        with open(fileName, "r") as f:
            lines=f.readlines()
            self.headerLine = lines[0]
            headerKeys = self.headerLine.strip().split("\t")
            for line in lines[1:]:
                
                gene = geneFeatures(line, self.headerLine)
                self.resistanceGenes[gene.name]=gene
                

        self.chromGenes={}
        self.chromGenes['gyrA_ecol']=[67,106]    
        self.chromGenes['gyrB_ecol']=[426,464]
        self.chromGenes['parC_ecol']=[47,133]
        self.chromGenes['parE_ecol']=[420,458]
        self.chromGenes['gyrA_kpne']=[67,106]
        self.chromGenes['gyrB_kpne']=[426,464]
        self.chromGenes['parC_kpne']=[47,133]
        self.chromGenes['parE_kpne']=[420,458]
        self.promoters=['P4', 'P5', 'Pa_Pb', 'Pc_Pd', 'ampC_promoter', 'P3', 'ampCpromoter']
    
    def getResistanceGeneFeature(self, geneName):
        if geneName in self.resistanceGenes:
            return self.resistanceGenes[geneName]
        return None
    
    def resistancePrediction(self, exactMatchList, mismatchDic, species, foN, sampleName, genotypeFile):


        fog = open(genotypeFile, 'w')
        fo= open(foN, "w")
        exactMatchList = sorted(exactMatchList)
        exactMatchList1 = [thisGene[0] for thisGene in exactMatchList]
        #TODO: need to check if have duplicated match 

        fo.write("#Antibiotic resistance for {0}\n".format(sampleName))
        fo.write("********** List of presence genes: ***********\n")
        #exactMatchGenes="\t".join(thisGene[0] for thisGene in exactMatchList )
        
        quinolones=['nalidixic acid', 'norfloxacin', 'ofloxacin', 'ciprofloxacin', 'levofloxacin', 'moxifloxacin']
        resistanceProfile = {}
        exceptionGeneList =['murA', 'uhpT', 'glpT', 'uhpA', 'uhpB', 'uhpC', 'ptsI', 'cyaA']
        fo.write("Exact gene matches: \n")
        for gene in exactMatchList:
            fo.write("\t{0}".format(gene[0]))
            #outGline: resistancePredictionProgram, sampleName, geneName, sequence coverage, percentIdentity, dnaMismatches(for exact match, only dnaMatch or protein match), coverage, protein truncated, comments- mapped contigs
            if gene[4]:
                copyNumber = gene[4]
            else:
                copyNumber = "None"
            outGLine= ["resistType", sampleName, gene[0], 100.0, 100.0,"None",copyNumber,"0", "|".join(map(str, gene[-1]+ [gene[1]]))]
            if gene[1] == 'protMatch':
                outGLine[4] = gene[2] #percent identity
                mismatches = gene[3]  #dnaMismatches
                if mismatches != None:
                    mutationList = ["{0}{1}{2}".format(mismatches[k][0], k+1, mismatches[k][1]) for k in mismatches.keys()]
                    mutationLine = "|".join(mutationList)
                    outGLine[5] = mutationLine

            fog.write("{0}\n".format(",".join(map(str, outGLine))))
            fo.write(":{0}|{1}|{2}".format(copyNumber, gene[1], "|".join(map(str, gene[-1]))))
            featureProfile = self.getResistanceGeneFeature(gene[0])
            if featureProfile == None:
                fo.write("\n")
                continue
            fo.write(",{0}\n".format(",".join('{0}:{1}'.format(key,val) for key, val in featureProfile.resistanceProfile.items())))
        
        fo.write("Inexact matches: {0} (sequences with many mismatches likely having indels causing frameshift, or sequence where there are many Ns because it's too  diverged from the reference sequence for the aligner to map to)\n".format(len(mismatchDic.keys()))) #print number of inexact matches
        
        for gene in sorted(mismatchDic.keys()):
            if mismatchDic[gene][4]:
                copyNumber = str(mismatchDic[gene][4])
            else:
                copyNumber = "None"
            fo.write("\t{0}:{1}|{2},{3:.2f}".format(gene,  copyNumber, "|".join(map(str, mismatchDic[gene][-1])), mismatchDic[gene][0])) 
                                               #print gene name, copyNumber variant, percent identity          
            outGLine=["resistType", sampleName, gene, 100.0, mismatchDic[gene][0], "", copyNumber, mismatchDic[gene][3], "|".join(map(str, mismatchDic[gene][-1]))]
            mismatches = mismatchDic[gene][1]
            dnaMutationList=["{0}{1}{2}".format(mismatches[k][0], k+1, mismatches[k][1]) for k in mismatches.keys()]
            outGLine[5]= "|".join(dnaMutationList)
            proteinMutationList = []
            if gene in self.promoters:
                mismatchDic[gene].append(dnaMutationList)
            else:
                mismatches = mismatchDic[gene][2]
                if mismatches:
                    proteinMutationList=["{0}{1}{2}".format(mismatches[k][0], k+1, mismatches[k][1]) for k in mismatches.keys()]
                    mismatchDic[gene].append(proteinMutationList)
            
            fo.write("|{0}|{1}".format(":".join(dnaMutationList), ":".join(proteinMutationList)))
            featureProfile = self.getResistanceGeneFeature(gene)
            if featureProfile == None:
                fo.write("\n")
                continue
            fo.write(",{0}\n".format(",".join('{0}:{1}'.format(key,val) for key, val in featureProfile.resistanceProfile.items())))
            fog.write("{0}\n".format(",".join(map(str, outGLine))))
        hasAmpC=False
        if 'ampC' in exactMatchList1 or ('ampC' in mismatchDic and mismatchDic['ampC'][0]>90):
            hasAmpC = True


        fo.write("********* Resistance prediction ********\n")
        for gene in exactMatchList:
            if gene[0] in self.chromGenes:
                continue
            if gene[0] in self.promoters:
                continue
            featureProfile = self.getResistanceGeneFeature(gene[0])
            if featureProfile == None:
                continue
            if 'TEM' in gene[0] and featureProfile.bushJacobyClass == '2b':
                isIRT = False
                if featureProfile.resistanceProfile['co-amoxiclav'] == 'R' and featureProfile.resistanceProfile['piperacillin-tazobactam'] == 'R':
                    isIRT = True

                thisGeneProfile={}
                for drug in featureProfile.resistanceProfile:
                    resistanceProfile[drug] = featureProfile.resistanceProfile[drug]
                thisGeneProfile= featureProfile.resistanceProfile
                TEM_promoters=['P3', 'P4', 'P5', 'Pa_Pb', 'Pc_Pd']
                promoters =  list(set(TEM_promoters).intersection(set(exactMatchList1)))
                if len(promoters) ==0:
                    continue
                thisPromoter = ''
                #Pairing promoter with TEM genes, only pair when two contig names are different (assemblies are too fragmented to stitch the contigs together) or when TEM and promoter are on the same contigs and are less than 100bp apart
                for promoter in promoters:
                    for x in exactMatchList:
                        if x[0] != promoter:
                            continue
                        if x[-1][0] == 'mapBased':
                            thisPromoter = promoter
                            break
                        contigName, start, end = x[-1][0:3]
                        if contigName != gene[-1][0]:
                            thisPromoter = promoter
                            break
                        if contigName == gene[-1][0]:
                            geneStart, geneEnd = gene[-1][1:3]
                            diff = None
                            if min(start, end) < min(geneStart, geneEnd):
                                diff = min(geneStart, geneEnd) - max(start, end)
                            else:
                                diff = max(geneStart, geneEnd) - min(start, end)
                            if diff < 0: diff = -diff
                            if diff < 100:
                                thisPromoter = promoter
                            #print diff
                            break
                if thisPromoter == '':
                    #print "This TEM gene does not have an exact match promoter associated"
                    continue
                if thisPromoter == 'P3':
                    if not isIRT: 
                        resistanceProfile['amoxicillin_ampicillin'] = 'R'
                        thisGeneProfile['amoxicillin_ampicillin'] = 'R'
                    for drug in ['cephalexin', 'co-amoxiclav', 'piperacillin-tazobactam']:
                        if drug in resistanceProfile and resistanceProfile[drug]!= 'S':
                            pass
                        else:
                            resistanceProfile['cephalexin'] = 'S'             ; thisGeneProfile['cephalexin'] = 'S' 
                            resistanceProfile['co-amoxiclav'] = 'S'           ; thisGeneProfile['co-amoxiclav'] = 'S'
                            resistanceProfile['piperacillin-tazobactam'] = 'S'; thisGeneProfile['piperacillin-tazobactam'] = 'S'
                elif thisPromoter in ['Pa_Pb', 'Pc_Pd']:
                    if not isIRT:
                        resistanceProfile['amoxicillin_ampicillin'] = 'R'; thisGeneProfile['amoxicillin_ampicillin'] = 'R'
                        resistanceProfile['co-amoxiclav'] = 'R'          ; thisGeneProfile['co-amoxiclav'] = 'R'
                    for drug in ['cephalexin', 'piperacillin-tazobactam']:
                        if drug in resistanceProfile and resistanceProfile[drug]!= 'S':
                            pass
                        else:
                            resistanceProfile['cephalexin'] = 'S'             ; thisGeneProfile['cephalexin'] = 'S' 
                            resistanceProfile['piperacillin-tazobactam'] = 'S'; thisGeneProfile['piperacillin-tazobactam'] = 'S'
                elif thisPromoter in ['P4', 'P5'] and not isIRT:
                    for drug in ['amoxicillin_ampicillin', 'cephalexin', 'co-amoxiclav', 'piperacillin-tazobactam']:
                        resistanceProfile[drug] = 'R'
                        thisGeneProfile[drug] = 'R'

                fo.write('{0}, promoter {1}, {2}\n'.format(gene[0], thisPromoter, ",".join('{0}:{1}'.format(key, val) for key,val in sorted(thisGeneProfile.items()))))
                continue

            if gene[0].startswith('ompK'):
                continue
            if gene[0].startswith('oqx') and species == 'Ecol':
                fo.write("{0}, none".format(gene[0]))
                for drug in quinolones:
                    fo.write(", {0}:R".format(drug))
                    resistanceProfile[drug] = 'R'
                fo.write("\n")
                continue
            if gene[0].startswith('dfr') :
                resistanceProfile['trimethoprim']='R'
                fo.write('{0},none, trimethoprim:R\n'.format(gene[0]))
                otherGene = []
                for thisGene in exactMatchList1:
                    if 'sul' in thisGene or 'folP' in thisGene:
                        otherGene.append(thisGene)
                for thisGene in mismatchDic:
                    if 'sul' in thisGene or 'folP' in thisGene and mismatchDic[thisGene][0]>=99: #percent identity
                        otherGene.append(thisGene)
                if len(otherGene) >0:
                    resistanceProfile['co-trimoxazole'] = 'R'
                    fo.write('{0}, {1}, co-trimoxazole:R\n'.format(gene[0], ":".join(otherGene)))
                continue
            if gene[0].startswith('sul') :
                resistanceProfile['sulfamethoxazole']='R'
                fo.write('{0},none, sulfamethoxazole:R\n'.format(gene[0]))
                continue
            if gene[0] == 'folP' :
                resistanceProfile['sulfamethoxazole']='R'
                fo.write('{0},none, sulfamethoxazole:low level R\n'.format(gene[0]))
                continue
            if  gene[0].startswith('qnrB') and 'aac6_prime-Ib-cr' in exactMatchList1:
                fo.write("{0}, aac6_prime-Ib-cr".format(gene[0]))
                for drug in quinolones:
                    fo.write(", {0}:R".format(drug))
                    resistanceProfile[drug] = 'R'
                fo.write("\n")
                continue
            if gene[0] in exceptionGeneList:
                continue

            fo.write("{0}, none, ".format(gene[0])) #none stands for no other associated gene/promoter combination
            if len(featureProfile.resistanceProfile.keys()) ==0 and featureProfile.summaryResClass != None:
                fo.write("{0}:R\n".format(featureProfile.summaryResClass))
            else:
                fo.write("{0}\n".format(",".join('{0}:{1}'.format(key,val) for key, val in featureProfile.resistanceProfile.items())))                

#1.            At least two gyrA amino acid QRDR mutations [Jacoby CID paper]
#2. 1 amino acid mutation in gyrA QRDR + 1 amino acid mutation in parC [Simone Bagel Impact of gyrA and parC Mutations on Quinolone Resistance, Doubling Time, and Supercoiling Degree of Escherichia coli] 
#3. for high level resistance (MIC > 8) 2 amino acid mutations in gyrA + 1 in parC
#4. presence of qnrB variant + aac(6')-Ib-cr (based on data supporting qnrB1 + aac(6')-Ib-cr leading to resistance)
#5. if gene_name = oqx* is present and query species is "ecol", then call resistant (significance of oqx* presence in klebsiella unclear)
#
#fosfomycin: murA (R with Cys115Asp; Asp369Asn; Leu370Ile), uhpT,  glpT, uhpA, uhpB, uhpC, ptsI, cyaA: R with mutations inactivating gene/downregulating expression, important mutations not known

        geneListWithMismatch = "\t".join(mismatchDic.keys())
        for gene in mismatchDic:
            featureProfile=  self.getResistanceGeneFeature(gene)
            if gene in ['ampC_promoter', 'ampCpromoter'] and mismatchDic[gene][0] > 0.9 and hasAmpC:
                mutationList=set( mismatchDic[gene][1])
                ampCMutShortList=set(['C110T', 'T120A'])
                ampCMutLongList=set(['C110T', 'T120A', 'G124A', 'G134A', 'G137A', 'T138G', 'C151T', 'C168T', 'C173T', 'T177G', 'A178T', 'G183A', 'G185A', 'C209T', 'T214C', 'C221T'])
                if len(mutationList & ampCMutShortList)>0:
                    fo.write("ampC|ampCpromoter,{0:.2f}:{1}, ".format(mismatchDic[gene][0], ":".join('{0}'.format(i) for i in mutationList & ampCMutShortList )))
                    fo.write("{0}\n".format(",".join('{0}:{1}'.format(key,val) for key, val in featureProfile.resistanceProfile.items())))
                elif len(mutationList & ampCMutLongList)>0:
                    fo.write("ampC|ampCpromoter,{0:.2f}:{1}\n".format(mismatchDic[gene][0], ":".join('{0}'.format(i) for i in mutationList & ampCMutLongList )))
            if gene == 'folP':
                pass
            if gene == 'murA':
                mutationList = set(mismatchDic[gene][1])
                murAMutList = set(['C115D', 'D369N', 'L370I'] )
                if len(mutationList & murAMutList) > 0:
                    fo.write('murA, {0:.2f}:{1}, fosfomycin:R\n'.format(mismatchDic[gene][0], ":".join('{0}'.format(i) for i in mutationList & murAMutList )))
            if gene in self.chromGenes:
                nQRDR = 0
                QRDRMuts = []
                #print gene, mismatchDic[gene], self.chromGenes[gene]
                for mutation in mismatchDic[gene][1]:
                    if mutation >= self.chromGenes[gene][0] and mutation <= self.chromGenes[gene][1]:
                        nQRDR +=1
                        QRDRMuts.append("{0}{1}{2}".format(mismatchDic[gene][1][mutation][0], mutation+1, mismatchDic[gene][1][mutation][1]))
                nQRDR_parC = 0
                QRDRMuts_parC = []
                parC  = None
                if 'parC_ecol' in mismatchDic:
                    parC = 'parC_ecol'
                if 'parC_kpne' in mismatchDic:
                    parC = 'parC_kpne'
                if parC   != None:
                    for mutation in mismatchDic[parC][1]:
                        if mutation >= self.chromGenes[parC][0] and mutation <= self.chromGenes[parC][1]:
                            nQRDR_parC +=1
                            QRDRMuts_parC.append("{0}{1}{2}".format(mismatchDic[parC][1][mutation][0], mutation+1, mismatchDic[parC][1][mutation][1]))
                
                if 'gyrA' not in gene:
                    if len(QRDRMuts) >0:
                        fo.write("{0}, {1}\n".format(gene, ":".join(QRDRMuts)))    
                    continue
                else:
                    if nQRDR >1 and nQRDR_parC ==0:
                        fo.write("{0}, {1}".format(gene, ":".join(QRDRMuts)))
                        for drug in quinolones:
                            fo.write(", {0}:R".format(drug))
                            resistanceProfile[drug] = 'R'
                        fo.write("\n")
                        continue
                    if nQRDR >1 and nQRDR_parC >0:
                        fo.write("{0}:{1}, {2}:{3}".format(gene, ":".join(QRDRMuts), parC, ":".join(QRDRMuts_parC)))
                        for drug in quinolones:
                            fo.write(", {0}:R(MIC>8)".format(drug))
                            resistanceProfile[drug] = 'R(MIC>8)'
                        fo.write("\n")
                        continue
                    
                    if nQRDR==1 and nQRDR_parC>0:
                        fo.write("{0}:{1}, {2}:{3}".format(gene, ":".join(QRDRMuts), parC, ":".join(QRDRMuts_parC)))
                        for drug in quinolones:
                            fo.write(", {0}:R".format(drug))
                            resistanceProfile[drug] = 'R'
                        fo.write("\n")
                        continue
                    if nQRDR>0:
                        fo.write("{0}, {1}\n".format(gene, ":".join(QRDRMuts)))
            elif gene  in exceptionGeneList:
                #TODO: what kind of profile in here
                continue
            elif gene.startswith("ompK"):
                #TODO
                continue
            elif gene.startswith("dfr") and mismatchDic[gene][0] > 95:
                resistanceProfile['trimethoprim']='R'
                fo.write('{0}, {1:.2f}, trimethoprim:R\n'.format(gene, mismatchDic[gene][0] ))
                otherGene = []
                for thisGene in exactMatchList1:
                    if 'sul' in thisGene or 'folP' in thisGene:
                        otherGene.append(thisGene)
                for thisGene in mismatchDic:
                    if 'sul' in thisGene or 'folP' in thisGene and mismatchDic[thisGene][0]> 99:
                        otherGene.append(thisGene)
                if len(otherGene) >0:
                    resistanceProfile['co-trimoxazole'] = 'R'
                    fo.write('{0}, {1}, co-trimoxazole:R\n'.format(gene, ":".join(otherGene)))
                continue
        if species == 'Kpne':
            resistanceProfile['ampicillin'] = 'R' #clinically do not use ampicillin for Kpne
            fo.write('Kpne,  none, ampicillin: R\n')
        fo.close()
        return

if __name__ == "__main__":
    usage = "usage: python %prog [options] \n" 
    version = "%prog 0.1"
    
    resistDB= CResistDB(sys.argv[1])



   
#genes with special treatment:
#need special treatment ampC C-42T; T-32A; C-11T is this nucleotide change or amino acid change? 
#need special treatment ampC promoter
#
#fusC 15 
#nimC 20
#nimD 20
#need special treatment ompK35: loss of ompK35 leads to resistance to cephalosporins and carbapenems (Role of Klebsiella pneumoniae OmpK35 Porin in Antimicrobial Resistance- Antimicrob. Agents Chemother. October 2003 vol. 47 no. 10) . Nicole, can you curate this? together with ompK36. I have read the paper but to put it into the table I need your expertise.
#Need special treatment ompK36
#SHV  A187T A22V D101H E31K G156D L35Q M129V M211L P226S S14F T18A Y7F

