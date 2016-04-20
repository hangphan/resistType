'''
Author
Hang Phan - hang.phan@ndm.ox.ac.uk
11 Jan  2016
Given a list of samples, look in to the folders containing resistPred.txt and get the resistance gene profile, this version contain the copy number report 
'''
from __future__ import division

import sys, os, pysam, gzip, logging, subprocess, uuid, shutil
import operator
import logging.handlers
import time, datetime
import os.path
from optparse import  OptionParser

_baseDir = '/well/gerton/hangphan/MMM/'

# Set up logging
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger('Log')
ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.setLevel(logging.DEBUG)
MODE1, MODE2, MODE3 = 0,1,2
class makeResistGeneProfiles(object):
    def __init__(self, args):
        self.sampleFile=args[0]
        self.modePred=args[1]
        self.modeGene=args[2]
        self.header=None
        self.drugDict={}
        self.drugList=[]
        self.samples=[]
        self.outputFN1 = self.sampleFile.split(".txt")[0] + "_phenotypePred.csv"
        self.outputFN2 = self.sampleFile.split(".txt")[0] + "_genotypePred.csv"
        self.geneList=['CTX-M', 'TEM', 'TEM-promoter', 'OXA', 'SHV', 'aac', 'aad', 'aph', 'ant', 'KPC', 'CMY', 'NDM', 'VIM', 'tet', 'gyrA']
        self.temPromoter=['P3', 'Pa_Pb', 'Pc_Pd', 'P4', 'P5', 'Pa/Pb', 'Pc/Pd']
    def makeHeaderPred(self): #make the header line for the phenotypic prediction
        header = ""
        
        if self.modePred == MODE1:
            with open("{0}/resources/drugnamesInResistDB.txt".format(_baseDir), "r") as f:
                for idx, line in enumerate(f):
                    self.drugDict[line.strip()] = 1 + idx
                    self.drugList.append( line.strip())
            self.header = ",".join(["sampleID"] + self.drugList)

        if self.modePred == MODE3:
            with open("{0}/HICF/jacResults_header.txt".format(_baseDir), "r") as f:
                for idx, line in enumerate(f):
                    self.drugDict[line.strip()] = 1 + idx
                    self.drugList.append( line.strip())
            self.header = ",".join(["sampleID"] + self.drugList)
    
        if self.modePred == MODE2:
            with open("{0}/resources/headerForHICF.csv".format(_baseDir), "r") as f:
                header = f.readlines()[0]
                cols = header.strip().split(",")
                for idx, item in enumerate(cols):
                    if ( "Routine lab" in item or "Phoenix" in item ) and "SIR" in item:
                        drugName = item.replace("Routine lab ", "").replace("Phoenix ", "").split()[0]
                        if drugName in drugList + ["amoxicillin", "ampicillin", "ertapenem", "meropenem"]:
                            self.drugDict[drugName] = idx
            self.header = header
        with open(self.outputFN1, 'w') as f:
            f.write("{0}\n".format( self.header))

        return None

    def makeHeaderGene(self):
        if self.modeGene==0:
            header= ['sampleID'] + self.geneList + ['Others']
            
        else:
            header = ['sampleID', 'matchState', 'pident', 'geneName', 'resistanceProfile']

        with open(self.outputFN2, 'w') as f:
            f.write("{0}\n".format(",".join(header)))

    def run(self):
    
        with open(self.sampleFile, "r") as f:
            for line in f:
                self.samples.append(line.strip().split()[0])
        self.makeHeaderPred()
        self.makeHeaderGene()
        for sample in self.samples:
            exactList, inexactList, resistanceProfile = self.readResult(sample)
            if exactList ==None and inexactList==None and resistanceProfile ==None:
                continue
            self.writeGenePred(sample, exactList, inexactList)
            self.writePhenoPred(sample, resistanceProfile)
        
    def readResult(self, sample):
        resultFile = "resistType/{0}/resistancePred.txt".format(sample)
        if not os.path.exists(resultFile):
            logger.info("{0}'s resistance prediction not found \n".format(sample))
            return None, None, None
        exactList={}
        inexactList={}
        with open(resultFile, "r") as f:
            lines = f.readlines()
            n=0
            exactMatch = 0
            inexactMatch = 0
            while n < len(lines):
                line = lines[n]
                n+=1
                if line.startswith("#"):
                    continue
                if "List of " in line:
                    continue
                if "Exact gene match" in line:
                    exactMatch = 1
                    continue
                if exactMatch:
                    if "Inexact" in line:
                        exactMatch = 0
                        inexactMatch = 1
                        continue
                    else:
                        cols = line.strip().rstrip().split(",")
                        geneName = cols[0].split(":")[0]
                        if geneName not in exactList:
                            exactList[geneName]= [line.strip()]
                        else:
                            exactList[geneName].append(line.strip())
                        
                if inexactMatch:
                    if "* Resistance prediction" in line:
                        break
                    cols = line.strip().rstrip().split(",")
                    if float(cols[1].split(":")[0]) < 80:
                        continue
                                                        
                    inexactList[cols[0].split(":")[0]]=[line.strip()]
                if "* Resistance prediction" in line:
                    break

        resistanceProfile = {}
        for line in lines[n:]:
            cols = line.strip().replace(", possible", "-possible").split(",")
            if len(cols) <=2:
                continue
            for x in cols[2:]:
                x = x.rstrip().strip()
                y = x.split(":")
                if y[0] not in resistanceProfile:
                    resistanceProfile[y[0]] = y[1]
                else:
                    if resistanceProfile[y[0] ] == 'S':
                        resistanceProfile[y[0]] = y[1]
                    if (resistanceProfile[y[0]] == '( R )' or "?" in resistanceProfile[y[0]]) and (y[1] == "R"  or y[1] == "R(MIC>8)"):
                        resistanceProfile[y[0]] == y[1]

        return exactList, inexactList, resistanceProfile

            
    def writeGenePred(self, sample, exactList, inexactList):
        foG=  open(self.outputFN2, "a")
        nCols=100        
        if self.modeGene ==0:  #report a selected set of genes in separate columns, then aggregate the remaining into one single column, separated by |
            outLine = ['']*(len(self.geneList) + 2)
            outCols=[]*(len(self.geneList) + 2)
            for idx in range(len(self.geneList)+ 2):
                outCols.append([])
            outLine[0] =  sample
            outCols[0] = [sample]

            for item in exactList:
                inList=0
                
                for idx, gene in enumerate(self.geneList):
                    if item.startswith(gene):
                        outCols[idx + 1].append(item)
                        inList=1
                        break
                if item in self.temPromoter:
                    idx = self.geneList.index("TEM-promoter")
                    outCols[idx +1].append(item )
                    inList=1
                if inList==0:
                    outCols[-1].append(item)
            for item in inexactList:
                inList =0
                pident = inexactList[item][0].split(",")[1].split(":")[0]
                #if float(pident) == 100:
                #    pident = "98.0"
                for idx, gene in enumerate(self.geneList):
                    if item.startswith(gene) and float(pident) >80:
                        outCols[idx + 1].append(item + ":" + pident)
                        inList=1
                        break
                if item in self.temPromoter:
                    idx = self.geneList.index("TEM-promoter")
                    outCols[idx +1].append(item + ":" + pident)
                    inList=1
                if inList==0:
                    outCols[-1].append(item + ":"+ pident)

            for idx in range(1, len(outCols)):
                outLine[idx] = "|".join(outCols[idx])
            foG.write("{0}\n".format(",".join(outLine)))

        else: #write one gene per line
            for item in sorted(exactList.keys()):
                lines= exactList[item]
                for line in lines:
                    cols = line.strip().rstrip().split(",")
                    outCols = [""]*nCols
                    outCols[0] = sample
                    outCols[1] = "exactMatch"
                    outCols[2] = "100.0"
                    outCols[3] = cols[0].split(":")[0]
                    if len(cols[0].split(":")) ==1:
                        print sample, line
                    outCols[4] = cols[0].split(":")[1]

                    for idx in range(1, len(cols)):
                        outCols[idx +4] = cols[idx].rstrip().strip()
                    
                    foG.write("{0}\n".format(",".join(outCols)))
            for item in sorted(inexactList.keys()):
                lines= inexactList[item]
                for line in lines:
                    cols = line.strip().rstrip().split(",")
                    outCols = [""]*nCols
                    outCols[0] = sample
                    outCols[1] = "inexactMatch"
                    outCols[2] = cols[1][:20].strip()
                    outCols[3]= cols[0].split(":")[0]
                    outCols[4] = cols[0].split(":")[1]
                    for idx in range(2, len(cols)):
                        outCols[idx +3] = cols[idx].rstrip().strip()
                    foG.write("{0}\n".format(",".join(outCols)))
            foG.close()
        return
    def writePhenoPred(self, sample, resistanceProfile):
        '''
        write the resistance profile prediction to output file, using several options
        '''
        nCols=100
        drugDict=self.drugDict
        line = [""] * nCols
        line[0] = sample
        if self.modePred==MODE1:
            for x in resistanceProfile:
                if x == "ampicillin":
                    line[self.drugDict["amoxicillin_ampicillin"]] = resistanceProfile[x]
                else:
                    line[self.drugDict[x]] = resistanceProfile[x]
        if self.modePred == MODE2:
            for x in resistanceProfile:
                if x in drugDict:
                    line[drugDict[x]] = resistanceProfile[x]
                if x == "amoxicillin_ampicillin":
                    line[drugDict["amoxicillin"]] = resistanceProfile[x]
                    line[drugDict["ampicillin"]] = resistanceProfile[x]
                if x == "carbapenems":
                    line[drugDict["ertapenem"]] = resistanceProfile[x]
                    line[drugDict["meropenem"]] = resistanceProfile[x]
                if x == "aminoglycosides and quinolones":
                    line[drugDict["aminoglycosides"]] = resistanceProfile[x]
                    line[drugDict["quinolones"]] = resistanceProfile[x]
        if self.modePred == MODE3:
            for x in resistanceProfile:
                if x in drugDict:
                    line[drugDict[x]] = resistanceProfile[x]
                if x == "amoxicillin_ampicillin":
                    line[drugDict["ampicillin"]] = resistanceProfile[x]
                if x == "carbapenems":
                    line[drugDict["ertapenem"]] = resistanceProfile[x]
                    line[drugDict["meropenem"]] = resistanceProfile[x]

        with open(self.outputFN1, 'a') as fo:
            fo.write("{0}\n".format(",".join(line)))
        return


if __name__ == "__main__":
    usage = "usage: python %prog [options] \n" 
    version = "%prog 0.1"
    
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-s", "--sampleFile", dest="sampleFile", type="string", default=None, help="File containing sample names")
    parser.add_option("-m", "--predMode", dest="outputPredictionMode", type="int", default=0, help="Mode of output for the resistance prediction \n\t\t\t0: all drug mode \n\t\t\t1:HICF comparison mode \n\t\t\t2: JAC set comparison mode")
    parser.add_option("-g", "--geneMode", dest="outputGeneMode", type="int", default=0, help="Mode of output for the resistance gene\n\t\t\t0: report selected genes mode, one row per sample\n\t\t\t1: report all genes mode, a sample cover several rows, with associated resistance mechanisms")
    (opts, args) = parser.parse_args()
    
    makeResistGeneProfileModule= makeResistGeneProfiles([opts.sampleFile, opts.outputPredictionMode, opts.outputGeneMode])
    makeResistGeneProfileModule.run()
