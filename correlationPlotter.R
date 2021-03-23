#the chip chip signal log values, can transform with 2 ^ (signalLogValue) to get the actual signal value for each gene
signalLogValuesList <- readRDS('signalLogValues.Rdata')

fileNames <- read.table('celFilePaths.txt', stringsAsFactors = F)[, 1]
fileNames <- strsplit(fileNames, '/')
fileNames <- unlist(lapply(fileNames, function(entry) {return(entry[2])}))
fileNames <- unique(unlist(strsplit(fileNames, '.CEL')))
fileNames <- gsub('-', '_', fileNames)


#get the actual signal values for each gene
signalValuesList <- sapply(fileNames, function(fileName, signalLogValuesList) {
  return(2 ^ signalLogValuesList[[fileName]])
}, signalLogValuesList, simplify = F)

geneList <- names(signalLogValuesList[["Assynchronous_Upf1_1_04082008"]])


valuesCombiner <- function(values1, values2) { #combines two datasets, scaling the 2nd to the 1st based on their medians
  med1 <- median(values1)
  med2 <- median(values2)
  values2 <- values2 * med1 / med2
  values <- c(values1, values2)
  return(values)
}


controlValues <- valuesCombiner(signalValuesList[["CONTROL"]], signalValuesList[["SB_Pombe_Control_04082008"]])
controlValues <- sapply(unique(names(controlValues)), function(gene, controlValues) {
  return(mean(controlValues[which(names(controlValues) == gene)]))
}, controlValues)

polControlValues <- valuesCombiner(signalValuesList[["IN1"]], signalValuesList[["IN3"]])
polControlValues <- sapply(unique(names(polControlValues)), function(gene, polControlValues) {
  return(mean(polControlValues[which(names(polControlValues) == gene)]))
}, polControlValues)



FPKMValues1 <- read.table('SPombeRNASeqData/Rep1.txt', header = T, stringsAsFactors = F)
FPKMValues2 <- read.table('SPombeRNASeqData/Rep2.txt', header = T, stringsAsFactors = F)
FPKMgenes <- FPKMValues1$gene_id[which(FPKMValues1$gene_id %in% geneList)] #get the genes that we have data for
FPKMValues1 <- sapply(FPKMgenes, function(gene, FPKMValues1) {
  return(FPKMValues1$FPKM[which(FPKMValues1$gene_id == gene)])
}, FPKMValues1) 
FPKMValues2 <- sapply(FPKMgenes, function(gene, FPKMValue2) {
  return(FPKMValues2$FPKM[which(FPKMValues2$gene_id == gene)])
}, FPKMValues2) 
FPKMValues <- sapply(FPKMgenes, function(gene_id, FPKMValues1, FPKMValues2) {
  return(mean(FPKMValues1[[gene_id]], FPKMValues2[[gene_id]]))
}, FPKMValues1, FPKMValues2)
logFPKMValues <- log(FPKMValues + 1)


library(ggplot2)


allPol2Values <- valuesCombiner(signalValuesList[["POL2"]], signalValuesList[["POL3"]])
meanPol2Values <- sapply(unique(names(logFPKMValues)), function(gene, allPol2Values) {
  return(mean(c(allPol2Values[which(names(allPol2Values) == gene)])))
}, allPol2Values)
meanPol2Signal <- sapply(names(logFPKMValues), function(gene, meanPol2Values, polControlValues) {
  return(unname(meanPol2Values[which(names(meanPol2Values) == gene)]) / unname(polControlValues[which(names(polControlValues) == gene)]))
}, meanPol2Values, polControlValues)


asynchSignal <- sapply(names(logFPKMValues), function(gene, signalValuesList, controlValues) {
  return(signalValuesList[["Assynchronous_Upf1_1_04082008"]][[gene]] / controlValues[[gene]])
}, signalValuesList, controlValues)

SPhaseSignal <- sapply(names(logFPKMValues), function(gene, signalValuesList, controlValues) {
  return(signalValuesList[["Upf1_S_phase_2_22032011"]][[gene]] / controlValues[[gene]])
}, signalValuesList, controlValues)


correlationPlotter <- function(xValues, yValues, xLabel, yLabel, yLimits, plotName) {
  sigTest <- cor.test(xValues, yValues, method = 'spearman')
  plotObject <- ggplot(mapping = aes(xValues, yValues)) +
    stat_density_2d(aes(fill = stat(density)), geom = 'raster', contour = FALSE) +
    scale_fill_viridis_c() +
    coord_cartesian(expand = FALSE) +
    geom_point(shape = '.', col = 'white') +
    xlab(xLabel) +
    ylab(yLabel) +
    ggtitle(paste('rho = ', signif(sigTest[["estimate"]], 3), ', p = ', signif(sigTest[["p.value"]], 3), sep = '')) +
    #xlim() +
    ylim(yLimits)
  pdf(plotName)
  print(plotObject)
  dev.off()
}



quantileValue <- quantile(logFPKMValues, probs = 0.1)
selectedIndices <- which(logFPKMValues >= quantileValue)
selectedGenes <- names(selectedIndices)


yLabel <- 'log(FPKM + 1)'
xLabel <- 'S-phase Upf1 ChIP-chip signal'
yValues <- logFPKMValues[selectedIndices]
yLimits <- c(min(yValues), max(yValues))
xValues <- SPhaseSignal[selectedIndices]
plotName <- 'RNAseqVsSphaseUpf1.pdf'
correlationPlotter(xValues, yValues, xLabel, yLabel, yLimits, plotName)
xLabel <- 'Asynchronous Upf1 ChIP-chip signal'
xValues <- asynchSignal[selectedIndices]
plotName <- 'RNAseqVsAsynchUpf1.pdf'
correlationPlotter(xValues, yValues, xLabel, yLabel, yLimits, plotName)
yLabel <- 'Pol2 ChIP-chip signal'
xLabel <- 'S-phase Upf1 ChIP-chip signal'
yValues <- meanPol2Signal[selectedIndices]
yLimits <- c(min(yValues), 3) #for upf1 (s-phase or asynch) vs pol2, when pol2 is on yaxis to cut out the anomolously high value from the plot
xValues <- SPhaseSignal[selectedIndices]
plotName <- 'Pol2VsSphaseUpf1.pdf'
correlationPlotter(xValues, yValues, xLabel, yLabel, yLimits, plotName)
xLabel <- 'Asynchronous Upf1 ChIP-chip signal'
xValues <- asynchSignal[selectedIndices]
plotName <- 'Pol2VsAsynchUpf1.pdf'
correlationPlotter(xValues, yValues, xLabel, yLabel, yLimits, plotName)


#save plotting data as table txt file, representing the values used in the correlation plots
plottingData <- as.data.frame(list('geneID' = selectedGenes, 'log_FPKM' = logFPKMValues[selectedIndices], 'Pol2_ChIP' = meanPol2Signal[selectedIndices], 'Asynchronous_Upf1_ChIP' = asynchSignal[selectedIndices], 'Sphase_Upf1_ChIP' = SPhaseSignal[selectedIndices]))
write.table(x = plottingData, file = 'processedDataTable.txt', quote = F, sep = '\t', row.names = F, col.names = T)


