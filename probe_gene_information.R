library(seqinr)

probe_data <- read.csv('ProbeDataTable.csv', stringsAsFactors = F)

cdf <- read.table('Sp20b_M_v04.CDF', sep = ' ')
genes <- unlist(sapply(strsplit(cdf[, 1], '='), function(gene) {
  if (length(gene) == 2) {
    if (gene[1] == 'Name' & gene[2] != 'NONE') {
      return(gene[2])
    }
  }
}))[-1]
probesByGene <- sapply(genes, function(gene) {
  return(matrix(ncol = 3, nrow = 0))
})
gene <- ''
lastProbeSequence <- ''
go <- FALSE
for (i in 1:nrow(cdf)) {
  print(i)
  if (substr(cdf[i, 1], 1, 4) == 'Name') {
    gene <- strsplit(cdf[i, 1], '=')[[1]][2]
  } else if (substr(cdf[i, 1], 1, 10) == 'CellHeader') {
    go <- TRUE
  } else if (substr(cdf[i, 1], 1, 2) == '[U') {
    go <- FALSE
  } else if (go) {
    probeInfo <- strsplit(strsplit(cdf[i, 1], '=')[[1]][2], '\t')[[1]]
    X <- probeInfo[1]
    Y <- probeInfo[2]
    probeSequence <- probeInfo[3]
    if (probeSequence == lastProbeSequence) {
      #mismatch
    } else {
      #perfect match
      probesByGene[[gene]] <- rbind(probesByGene[[gene]], c(probeSequence, X, Y))
    }
    lastProbeSequence <- probeSequence
  }
}
probeInfo <- matrix(nrow = 0, ncol = 4)
for (gene in names(probesByGene)) {
  print(gene)
  geneProbes <- probesByGene[[gene]]
  geneProbes <- cbind(rep(gene, nrow(geneProbes)), geneProbes)
  probeInfo <- rbind(probeInfo, geneProbes)
}


compliment <- c('A', 'T', 'C', 'G', 'N')
names(compliment) <- c('T', 'A', 'G', 'C', 'N')


reverse_compliment <- function(sequence, compliment) {
  sequence_vector <- strsplit(sequence, c())[[1]]
  complimentary_sequence_vector <- unname(compliment[sequence_vector])
  reverse_complimentary_vector <- rev(complimentary_sequence_vector)
  reverse_compliment_sequence <- paste(reverse_complimentary_vector, collapse = '')
  return(reverse_compliment_sequence)
}

#including the UTRs is important because the TSS is the first base of the 5' UTR and the TES is the last base of the 3' UTR
gene_data <- read.fasta('S.pombeCDS+introns+UTRs/cds+introns+utrs.fa', as.string = T, forceDNAtolower = F) #downloaded from https://www.pombase.org/downloads/genome-datasets
locusIDs <- names(gene_data)

gene_information <- probe_data
gene_information$ProbeSequence <- probeInfo[, 2]
gene_information$GeneLength <- 0
gene_information$FivePrimePos <- 0
gene_information$ThreePrimePos <- 0
probe_setIDs <- unique(gene_information$ProbeSetName)

MatchIndex <- which(probe_setIDs %in% locusIDs)
NoOfProbeGeneMatches <- length(MatchIndex)

NoMatchIndex <- (1:length(probe_setIDs))[-MatchIndex]

for(i in NoMatchIndex) { #to remove probes sets for which we dont have corresponding gene sequence data
  #cat(i)
  NonMatchingProbeID <- probe_setIDs[i]
  gene_information <- gene_information[which(gene_information$ProbeSetName != NonMatchingProbeID), ]
}

for(index in MatchIndex) { #get gene length and index (in the pombe gene sequence data) for each probe set with corresponding info in pombe data
  #cat(index)
  MatchingProbeID <- probe_setIDs[index]
  gene <- gene_data[[MatchingProbeID]][1]
  gene_length <- nchar(gene)
  gene_information$GeneLength[which(gene_information$ProbeSetName == MatchingProbeID)] <- gene_length
}


SenseAntiSense <- c(1, 2)
ProbesOutsideGenesIndices <- list() #vector of probes that are not being found within the gene they are supposed to be inside, this could be because they are slightly outside the gene or because they are not perfect matches to the genomic sequence if its a new assembly to the one that the probes were originally
counter <- 1
ProbeSequences <- gene_information$ProbeSequence
ProbeSetNames <- gene_information$ProbeSetName
GeneLengths <- gene_information$GeneLength
for(i in 1:nrow(gene_information)){
  print(i)
  ProbeSeq <- ProbeSequences[i]
  Gene <- gene_data[[ProbeSetNames[i]]][1]
  GeneRC <- reverse_compliment(Gene, compliment)
  SensePosition <- gregexpr(ProbeSeq, Gene)[[1]][1]
  AntiSensePosition <- gregexpr(ProbeSeq, GeneRC)[[1]][1]
  StrandIndex <- which(c(SensePosition, AntiSensePosition) != -1)
  if(length(StrandIndex) > 0) {
    gene_information$Direction[i] <- StrandIndex
    FirstBasePos <- c(SensePosition, AntiSensePosition)[StrandIndex]
    LastBasePos <- FirstBasePos + 24
    if(StrandIndex == 2) { #to flip the positions the right way around if the probe is on the antisense strand
      GeneLength <- GeneLengths[i]
      FirstBasePos <- GeneLength - FirstBasePos + 1 #so the TSS position is position 1 and the TES is position length of gene
      LastBasePos <- GeneLength - LastBasePos + 1
      gene_information$FivePrimePos[i] <- LastBasePos
      gene_information$ThreePrimePos[i] <- FirstBasePos
    } else {
      gene_information$FivePrimePos[i] <- FirstBasePos
      gene_information$ThreePrimePos[i] <- LastBasePos
    }
  } else {
    print(i)
    ProbesOutsideGenesIndices[[counter]] <- i
    counter <- counter + 1
  }
}
ProbesOutsideGenesIndices <- unlist(ProbesOutsideGenesIndices)

# gene_info_non_matching_probes <- gene_information[ProbesOutsideGenesIndices,]
# write.csv(gene_info_non_matching_probes, 'GeneNonMatchProbeInformation.csv', row.names = F)
gene_information <- gene_information[-ProbesOutsideGenesIndices, ]
write.csv(gene_information, 'ProbeGeneInformation.csv', row.names = F)

