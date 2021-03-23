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

gtfFilePath = sys.argv[1] #gtfFilePath = 'pombeGenomeIndex/Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Genes/genes.gtf'
pickleFilePathIP = sys.argv[2] #pickleFilePathIP = 'geneCoverages/SB1.pickle'
pickleFilePathInput = sys.argv[3] #pickleFilePathInput = 'geneCoverages/SB10.pickle'
flankingSize = int(sys.argv[4]) #flankingSize = 1000
plotFilePath = sys.argv[5] #plotFilePath = 'plots/test.pdf'
inputGeneSet = bool(sys.argv[6] == 'True') #inputGeneSet = True / False as to whether a .txt file with desired gene set

if inputGeneSet:
    geneSetFilePath = sys.argv[7] #geneSetFilePath = 'geneList.txt'
    geneSet =  list(csv.reader(open(geneSetFilePath), delimiter='\t'))
    geneSet = [gene[0] for gene in geneSet]

gtfContents = list(csv.reader(open(gtfFilePath), delimiter='\t'))
chromosomes = list(set([line[0] for line in gtfContents]))
genesByChromosome = {chromosome: defaultdict(tuple) for chromosome in chromosomes}
for line in gtfContents:
    if line[2] == 'transcript':
        chromosome = line[0]
        geneID = [entry for entry in line[8].split(';') if 'gene_id' in entry][0].split('"')[1]
        strand = line[6]
        start = int(line[3])
        end = int(line[4])
        genesByChromosome[chromosome][geneID] = (strand, start, end)

maxChromosomeBasePositions = {chromosome: numpy.max([genesByChromosome[chromosome][geneID][2] for geneID in genesByChromosome[chromosome].keys()]) + 20000 for chromosome in chromosomes}

metageneValues = {'IP': {}, 'Input': {}}

#store IP metagene values
geneCoverages = pickle.load(open(pickleFilePathIP, 'rb'))
geneOrder = []
geneValues, upstreamValues, downstreamValues = [], [], []
for chromosome in chromosomes:
        for geneID in genesByChromosome[chromosome]:
            geneLength = len(geneCoverages[chromosome][geneID]) - (2 * flankingSize) #take off the 1000 bp either side of the gene
            interval = geneLength / 100
            intervalBreaks = [interval * index for index in range(100 + 1)]
            intervalStartEndVals = [(math.ceil(intervalBreaks[index]) + (flankingSize + 1), math.floor(intervalBreaks[index + 1]) + 1 + (flankingSize + 1)) for index in range(100)]
            intervalCoverages = [numpy.mean([geneCoverages[chromosome][geneID][position] for position in range(intervalStartEndVals[index][0], intervalStartEndVals[index][1])]) for index in range(len(intervalStartEndVals))]
            upstreamCoverage = geneCoverages[chromosome][geneID][:(flankingSize + 1)] #include TSS in single-base resolution
            downstreamCoverage = geneCoverages[chromosome][geneID][(geneLength + (flankingSize - 1)):] #include TES in single-base resolution
            if not any([math.isnan(value) for value in intervalCoverages]):  # cut out genes with length < 200, which are too small to properly divide into 100 % averaged values
                geneValues.append(intervalCoverages)
                upstreamValues.append(upstreamCoverage)
                downstreamValues.append(downstreamCoverage)
                geneOrder.append(geneID)
metageneValues['IP'] = {'upVals': upstreamValues, 'geneVals': geneValues, 'downVals': downstreamValues}

#store input metagene values
geneCoverages = pickle.load(open(pickleFilePathInput, 'rb'))
geneOrder = []
geneValues, upstreamValues, downstreamValues = [], [], []
for chromosome in chromosomes:
        for geneID in genesByChromosome[chromosome]:
            geneLength = len(geneCoverages[chromosome][geneID]) - (2 * flankingSize) #take off the 1000 bp either side of the gene
            interval = geneLength / 100
            intervalBreaks = [interval * index for index in range(100 + 1)]
            intervalStartEndVals = [(math.ceil(intervalBreaks[index]) + (flankingSize + 1), math.floor(intervalBreaks[index + 1]) + 1 + (flankingSize + 1)) for index in range(100)]
            intervalCoverages = [numpy.mean([geneCoverages[chromosome][geneID][position] for position in range(intervalStartEndVals[index][0], intervalStartEndVals[index][1])]) for index in range(len(intervalStartEndVals))]
            upstreamCoverage = geneCoverages[chromosome][geneID][:(flankingSize + 1)] #include TSS in single-base resolution
            downstreamCoverage = geneCoverages[chromosome][geneID][(geneLength + (flankingSize - 1)):] #include TES in single-base resolution
            if not any([math.isnan(value) for value in intervalCoverages]):  # cut out genes with length < 200, which are too small to properly divide into 100 % averaged values
                geneValues.append(intervalCoverages)
                upstreamValues.append(upstreamCoverage)
                downstreamValues.append(downstreamCoverage)
                geneOrder.append(geneID)
metageneValues['Input'] = {'upVals': upstreamValues, 'geneVals': geneValues, 'downVals': downstreamValues}

#if you want to plot the metagene for a specific set of genes

if inputGeneSet:
    geneIndices = [index for index in range(len(geneOrder)) if geneOrder[index] in geneSet]
else:
    geneIndices = list(range(len(geneOrder)))


selectedInputgeneVals, selectedInputupVals, selectedInputdownVals = [metageneValues['Input']['geneVals'][index] for index in geneIndices], [metageneValues['Input']['upVals'][index] for index in geneIndices], [metageneValues['Input']['downVals'][index] for index in geneIndices]
inputGeneValsArray, inputUpValsArray, inputDownValsArray = numpy.array(selectedInputgeneVals), numpy.array(selectedInputupVals), numpy.array(selectedInputdownVals)
inputPercentageMeans, inputUpMeans, inputDownMeans = numpy.mean(inputGeneValsArray, axis=0), numpy.mean(inputUpValsArray, axis=0), numpy.mean(inputDownValsArray, axis=0)

selectedIPgeneVals, selectedIPupVals, selectedIPdownVals = [metageneValues['IP']['geneVals'][index] for index in geneIndices], [metageneValues['IP']['upVals'][index] for index in geneIndices], [metageneValues['IP']['downVals'][index] for index in geneIndices]
ipGeneValsArray, ipUpValsArray, ipDownValsArray = numpy.array(selectedIPgeneVals), numpy.array(selectedIPupVals), numpy.array(selectedIPdownVals)
ipPercentageMeans, ipUpMeans, ipDownMeans = numpy.mean(ipGeneValsArray, axis=0), numpy.mean(ipUpValsArray, axis=0), numpy.mean(ipDownValsArray, axis=0)

upstreamVals, geneBodyVals, downStreamVals = ipUpMeans / inputUpMeans, ipPercentageMeans / inputPercentageMeans, ipDownMeans / inputDownMeans

fig, axes = matplotlib.pyplot.subplots(1, 3, sharex=False, sharey=True, figsize=(12, 6), gridspec_kw={'width_ratios': [1, 2, 1]})
fig.subplots_adjust(wspace=0, hspace=0)
percentages = [index + 0.5 for index in range(100)]
percentages[0], percentages[-1] = 0, 100
axes[0].plot(list(range(-flankingSize, 1)), upstreamVals)
axes[1].plot(percentages, geneBodyVals)
axes[2].plot(list(range(0, flankingSize + 1)), downStreamVals)
axes[0].set_xlim(-flankingSize, 0)
axes[1].set_xlim(0, 100)
axes[1].xaxis.tick_top()
axes[2].set_xlim(0, flankingSize)
#fig.suptitle(plotTitle)
fig.text(0.25, 0.04, 'distance from TSS (bp)', ha='center')
fig.text(0.5, 0.08, 'position through gene (%)', ha='center')
fig.text(0.75, 0.04, 'distance from TES (bp)', ha='center')
fig.text(0.06, 0.5, 'normalised coverage depth', va='center', rotation='vertical')
matplotlib.pyplot.savefig(plotFilePath)
#fig.show()

