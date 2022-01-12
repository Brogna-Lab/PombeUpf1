ProbeData <- read.csv('ProbeGeneInformation.csv')
ProbeData$PositionValue <- 0
RowNumber <- nrow(ProbeData)
fivePrimePositions <- ProbeData$FivePrimePos
maxPositions <- ProbeData$GeneLength - 24
positionValues <- fivePrimePositions * 100 / maxPositions
ProbeData$PositionValue <- positionValues


library(PtProcess)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)


PositionMeanFinder <- function(Position, Data, Positions) {
  PositionSpecificIntensities <- Data[which(Positions == Position)]
  AverageIntensity <- mean(PositionSpecificIntensities)
  return(AverageIntensity)
}

PlottingWrapperFunction <- function(Data, PositionMeanFinder, span, Positions, PositionClasses, title) {
  positionRange <- seq(1, 1000, 1) / 10 - 0.05
  mus <- c()
  sigmas <- c()
  stdErrors <- c()
  for (position in positionRange) {
    currentIndices <- which(Positions >= position - 0.05 & Positions < position + 0.05)
    mus <- c(mus, mean(Data[currentIndices]))
    sigmas <- c(sigmas, sd(Data[currentIndices]))
    stdErrors <- c(stdErrors, sd(Data[currentIndices]) / sqrt(length(currentIndices)))
  }
  plot(positionRange, mus, type = 'l', xlab = "Position in gene body (%)", ylab = "Normalised signal intensity (AU)", main = title)
  # plot(positionRange, mus, type = 'l', xlab = "Position in gene body (%)", ylab = "Normalised signal intensity (AU)", main = title, ylim = c(min(mus - stdErrors), max(mus + stdErrors)))
  # lines(positionRange, mus + stdErrors, col = 'red')
  # lines(positionRange, mus - stdErrors, col = 'red')
}

Positions <- ProbeData$PositionValue
controlMeans <- sapply(1:nrow(ProbeData), function(index, ProbeData) {
  return(mean(c(ProbeData$PMupfassynCONTROL[index] / ProbeData$MMupfassynCONTROL[index], ProbeData$PMupfassynSBPombeControl04082008[index] / ProbeData$MMupfassynSBPombeControl04082008[index])))
}, ProbeData)
pol2ControlMeans <- sapply(1:nrow(ProbeData), function(index, ProbeData) {
  return(mean(c(ProbeData$PMRNAPolsChIPsIN1[index] / ProbeData$MMRNAPolsChIPsIN1[index], ProbeData$PMRNAPolsChIPsIN3[index] / ProbeData$MMRNAPolsChIPsIN3[index])))
}, ProbeData)
normalisedAsynchUpf1 <- (ProbeData$PMupfassynAssynchronousUpf1104082008 / ProbeData$MMupfassynAssynchronousUpf1104082008) / controlMeans
normalisedSphaseUpf1 <- (ProbeData$PMupfs1Upf1Sphase222032011 / ProbeData$MMupfs1Upf1Sphase222032011) / controlMeans
normalisedPol2_1 <- (ProbeData$PMRNAPolsChIPsPOL2 / ProbeData$MMRNAPolsChIPsPOL2) / controlMeans
normalisedPol2_2 <- (ProbeData$PMRNAPolsChIPsPOL3 / ProbeData$MMRNAPolsChIPsPOL3) / controlMeans
normalisedPol2 <- sapply(1:length(normalisedPol2_1), function(i) {
  return(mean(c(normalisedPol2_1[i], normalisedPol2_2[i])))
})
span <- 1

Data <- normalisedAsynchUpf1
title <- 'Asynchronous Upf1 signal (PM / MM) normalised by mean of two control signals (PM / MM)'
PlottingWrapperFunction(Data, PositionMeanFinder, span, Positions, PositionClasses, title)

Data <- normalisedSphaseUpf1
title <- 'S-phase Upf1 signal (PM / MM) normalised by mean of two control signals (PM / MM)'
PlottingWrapperFunction(Data, PositionMeanFinder, span, Positions, PositionClasses, title)

Data <- controlMeans
title <- 'Control mean signal (PM / MM)'
PlottingWrapperFunction(Data, PositionMeanFinder, span, Positions, PositionClasses, title)

Data <- normalisedPol2 #
title <- 'Pol2 mean signal (PM / MM) normalised by mean of two Pol2 control signals (PM / MM)'
PlottingWrapperFunction(Data, PositionMeanFinder, span, Positions, PositionClasses, title)

Data <- pol2ControlMeans
title <- 'Pol2 controls mean signal (PM / MM)'
PlottingWrapperFunction(Data, PositionMeanFinder, span, Positions, PositionClasses, title)






