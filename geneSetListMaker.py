import csv
import scipy
import numpy
import math
from _collections import defaultdict
import sys
import pickle

#can use the FPKM values to look at genes with different expression levels
RNAseqData1 = list(csv.reader(open('Rep1.txt'), delimiter='\t'))
RNAseqData2 = list(csv.reader(open('Rep2.txt'), delimiter='\t'))
RNAseqDict1, RNAseqDict2 = defaultdict(str), defaultdict(str)
RNAseqGenes = [RNAseqData1[i][3] for i in range(1, len(RNAseqData1))]
for i in range(1, len(RNAseqData1)):
    gene1, gene2, FPKM1, FPKM2 = RNAseqData1[i][3], RNAseqData2[i][3], RNAseqData1[i][9], RNAseqData2[i][9]
    RNAseqDict1[gene1] = float(FPKM1)
    RNAseqDict2[gene2] = float(FPKM2)
RNAseqDict = {gene: numpy.mean([RNAseqDict1[gene], RNAseqDict2[gene]]) for gene in RNAseqGenes}
RNAseqData = [(gene, RNAseqDict[gene]) for gene in RNAseqGenes]
RNAseqData.sort(key=lambda x: x[1], reverse=True) #sorted largest FPKM to lowest

#get the top most expressed genes according to RNA-seq data which are also protein coding genes and save to text file
geneSet = [RNAseqData[index][0] for index in range(400) if any([IDtype in RNAseqData[index][0] for IDtype in ['SPAC', 'SPBC', 'SPCC', 'SPAP', 'SPBP', 'SPCP']])]
outputFile = open('geneList.txt', 'w')
for gene in geneSet:
    outputFile.write(gene + '\n')
outputFile.close()

