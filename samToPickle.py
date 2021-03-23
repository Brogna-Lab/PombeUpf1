import csv
import scipy
import numpy
import math
from _collections import defaultdict
import matplotlib
import matplotlib.pyplot
import sys
import pickle
samFilePath = sys.argv[1]
gtfFilePath = sys.argv[2]
pickleFilePath = sys.argv[3]
flankingSize = int(sys.argv[4])
gtfContents = list(csv.reader(open(gtfFilePath), delimiter='\t'))

samContents = list(csv.reader(open(samFilePath), delimiter='\t'))
samHeaders = samContents[:8]
samContents = samContents[8:]
normalisingValue = len(samContents)

#now create a dictionary with the information of the number of reads covering every base in the genome for every chromosome from the sam file, then can for each gene according to the transcript
#entry from the gtf file and extract the basewise coverage scores and then make metagenes from that information

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

#now make a dictionary for every base in the genome
coverageDict = {chromosome: {base: 0 for base in range(-20000, maxChromosomeBasePositions[chromosome])} for chromosome in chromosomes}

#now go through and for every read and do +1 coverage for every base between and including the start and end positions
for alignment in samContents:
    chromosome = alignment[2]
    start = int(alignment[3]) - 1
    size = len(alignment[9])
    end = start + size
    for position in range(start, end):
        coverageDict[chromosome][position] += 1 / normalisingValue

geneCoverages = {chromosome: {geneID: None for geneID in genesByChromosome[chromosome].keys()} for chromosome in chromosomes}
for chromosome in chromosomes:
    for geneID in list(geneCoverages[chromosome].keys()):
        strand, start, end = genesByChromosome[chromosome][geneID]
        coordExtremes = [start - flankingSize, end + 1 + flankingSize]
        if strand == '-':
            coordExtremes.reverse() #coordRange = list(reversed(coordRange))
            coordRange = list(range(coordExtremes[0], coordExtremes[1], -1))
        else:
            coordRange = list(range(coordExtremes[0], coordExtremes[1]))
        geneCoverages[chromosome][geneID] = [coverageDict[chromosome][base] for base in coordRange]

pickle.dump(geneCoverages, open(pickleFilePath, 'wb'))


