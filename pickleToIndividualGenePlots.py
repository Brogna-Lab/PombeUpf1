import csv
import scipy
import numpy
import math
from _collections import defaultdict
import matplotlib
import matplotlib.pyplot
import sys
import pickle
import matplotlib.gridspec

pickleFilePathIP = sys.argv[1] #pickleFilePathIP = 'geneCoverages/SB1.pickle'
pickleFilePathInput = sys.argv[2] #pickleFilePathInput = 'geneCoverages/SB10.pickle'
flankingSize = int(sys.argv[3]) #flankingSize = 1000
plotDirectoryPath = sys.argv[4] #plotDirectoryPath = 'plots/'
geneSetFilePath = sys.argv[5] #geneSetFilePath = 'geneList2.txt'
geneSet =  list(csv.reader(open(geneSetFilePath), delimiter='\t'))
geneSet = [gene[0] for gene in geneSet]


geneCoverages = pickle.load(open(pickleFilePathIP, 'rb'))
chromosomes = list(geneCoverages.keys())

IPgeneValues = {geneID: [] for geneID in geneSet}
for chromosome in chromosomes:
    for gene in geneSet:
        if gene in geneCoverages[chromosome].keys():
            IPgeneValues[gene] = numpy.array(geneCoverages[chromosome][gene])


geneCoverages = pickle.load(open(pickleFilePathInput, 'rb'))
inputgeneValues = {geneID: [] for geneID in geneSet}
for chromosome in chromosomes:
    for gene in geneSet:
        if gene in geneCoverages[chromosome].keys():
            inputgeneValues[gene] = numpy.array(geneCoverages[chromosome][gene])


normalisedGeneValues = {gene: IPgeneValues[gene] / inputgeneValues[gene] for gene in geneSet}
xlab = 'distance from TSS (bp)'
ylab = 'normalised ChIP-seq signal'
for gene in geneSet:
    xVals = numpy.array(list(range(-flankingSize, len(normalisedGeneValues[gene]) - flankingSize)))
    yVals = normalisedGeneValues[gene]
    matplotlib.pyplot.figure(figsize=(12, 6))
    matplotlib.pyplot.axvline(0, color='black')
    matplotlib.pyplot.axvline(len(normalisedGeneValues[gene]) - flankingSize * 2 - 1, color='black')
    matplotlib.pyplot.plot(xVals, yVals)
    matplotlib.pyplot.xlabel(xlab)
    matplotlib.pyplot.ylabel(ylab)
    matplotlib.pyplot.savefig('plots/' + gene + '.pdf')
    #matplotlib.pyplot.show()

