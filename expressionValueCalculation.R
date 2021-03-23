#now that we have the values with the noise corrected for, we can calculate the expression values of each gene
adjustedIntensitiesList <- readRDS('adjustedIntensities.Rdata') 


#convert adjusted intensities back from list to environment structure 
fileNames <- names(adjustedIntensitiesList)
adjustedIntensities <- list2env(sapply(fileNames, function(fileName, adjustedIntensitiesList) {
  return(list2env(adjustedIntensitiesList[[fileName]]))
}, adjustedIntensitiesList, simplify = F))

allCoordsKey <- ls(adjustedIntensities[["Assynchronous_Upf1_1_04082008"]])
allCoords <- list2env(sapply(allCoordsKey, function(coord) {return(as.numeric(strsplit(coord, '-')[[1]]))}, simplify = F))

#so the expression values are calculated as follows, the intensities have already been preprocessed for global background noise, these will be 
#taken and an ideal mismatch value is calculated and subtracted to adjust the PM intensity, these are then log-transformed to stabilise the variance
#then the biweight estimator will provide a robust mean (per gene) of resulting values, with the signal then being the antilog of that mean
probeData <- read.csv('ProbeDataTable.csv', header = T, stringsAsFactors = F)
probeSetNames <- unique(probeData$ProbeSetName)
probeSetIntensities <- sapply(fileNames, function(fileName) {
  return(list2env(sapply(probeSetNames, function(probeSetName) {
    return(new.env())
    }, simplify = F)))
  })
counters <- list2env(sapply(probeSetNames, function(probeSetNames) {
  return(0)
  }, simplify = F))
for (rowIndex in 1:nrow(probeData)) {
  PMcoords <- c(probeData[rowIndex, "PMX"], probeData[rowIndex, "PMY"])
  MMcoords <- c(probeData[rowIndex, "MMX"], probeData[rowIndex, "MMY"])
  PMkey <- paste(PMcoords[1], PMcoords[2], sep = '-')
  MMkey <- paste(MMcoords[1], MMcoords[2], sep = '-')
  probeSetName <- probeData[rowIndex, "ProbeSetName"]
  counters[[probeSetName]] <- counters[[probeSetName]] + 1
  index <- paste(counters[[probeSetName]])
  for (fileName in fileNames) {
    PMsignal <- adjustedIntensities[[fileName]][[PMkey]]
    MMsignal <- adjustedIntensities[[fileName]][[MMkey]]
    probeSetIntensities[[fileName]][[probeSetName]][[index]] <- list("PM" = PMsignal, "MM" = MMsignal) #add the probe pair noise-corrected intensities
  }
}

#get the specific background (SB) for each probe set using the tukey biweight of the log2(PM) - log2(MM) probe values,
#this means the specific background is just a robust estimate of the reliability of the probe set values, since the larger the (SB),
#the larger the PM signals are compared to the MM signals

#first get the probe log ratio values we need for the specific background
probeLogRatios <- list2env(sapply(fileNames, function(fileName, probeSetIntensities, probeSetNames) {
  return(list2env(sapply(probeSetNames, function(probeSetName, fileName, probeSetIntensities) {
    return(list2env(sapply(ls(probeSetIntensities[[fileName]][[probeSetName]]), function(index, probeSetIntensities, fileName, probeSetName) {
      return(log2(probeSetIntensities[[fileName]][[probeSetName]][[index]][["PM"]]) - log2(probeSetIntensities[[fileName]][[probeSetName]][[index]][["MM"]]))
      }, probeSetIntensities, fileName, probeSetName, simplify = F)))
    }, fileName, probeSetIntensities, simplify = F)))
  }, probeSetIntensities, probeSetNames, simplify = F))

#now calculate the specific background using the one-step biweight
C <- 5 #tuning constant
e <- 0.0001 #avoids divisions by zero

oneStepBiweight <- function(xVec, C, e) {
  M <- median(xVec)
  S <- median(abs(xVec - M))
  uVec <- (xVec - M) / (C * S + e)
  wVec <- sapply(uVec, function(u) {
    if (abs(u) <= 1) {
      return((1 - (u ^ 2)) ^ 2)
    } else {
      return(0)
    }
  }) 
  Tbi <- sum((xVec * wVec) / sum(wVec))
  return(Tbi)
}

SB <- list2env(sapply(fileNames, function(fileName, probeLogRatios, probeSetNames, oneStepBiweight) {
  return(list2env(sapply(probeSetNames, function(probeSetName, probeLogRatios, fileName, oneStepBiweight) {
    xVec <- unname(sapply(ls(probeLogRatios[[fileName]][[probeSetName]]), function(index, probeLogRatios, fileName, probeSetName) {
      return(probeLogRatios[[fileName]][[probeSetName]][[index]])
    }, probeLogRatios, fileName, probeSetName))
    return(oneStepBiweight(xVec, C, e))
  }, probeLogRatios, fileName, oneStepBiweight, simplify = F)))
}, probeLogRatios, probeSetNames, oneStepBiweight, simplify = F))


#now get ideal mismatch (IM) values for each probe pair based on whether MM < PM, MM >= PM & SB > contrast, or MM >= PM & SB <= contrast
contrastTau <- 0.03 #for determining how small the SB should be before we start making the IM depend more on the PM
scaleTau <- 10
idealIntensities <- list2env(sapply(fileNames, function(fileName, probeSetNames, probeSetIntensities, SB) {
  return(list2env(sapply(probeSetNames, function(probeSetName, probeSetIntensities, fileName, SB) {
    SBi <- SB[[fileName]][[probeSetName]]
    return(list2env(sapply(ls(probeSetIntensities[[fileName]][[probeSetName]]), function(index, probeSetIntensities, fileName, probeSetName, SBi) {
      PM <- probeSetIntensities[[fileName]][[probeSetName]][[index]][["PM"]]
      MM <- probeSetIntensities[[fileName]][[probeSetName]][[index]][["MM"]]
      if (MM < PM) {
        return(list2env(list('PM' = PM, 'IM' = MM)))
      } else if (SBi > contrastTau) {
        return(list2env(list('PM' = PM, 'IM' = PM / (2 ^ SBi))))
      } else {
        return(list2env(list('PM' = PM, 'IM' = PM / (2 ^ (contrastTau / (1 + (contrastTau - SBi) / scaleTau))))))
      }
    }, probeSetIntensities, fileName, probeSetName, SBi, simplify = F)))
  }, probeSetIntensities, fileName, SB, simplify = F)))
}, probeSetNames, probeSetIntensities, SB, simplify = F))

#now get the probe value, which is the log-transformed difference between PM and IM
d <- 2 ^ (- 20) #value to guarantee numerical stability
probeValues <- list2env(sapply(fileNames, function(fileName, idealIntensities, probeSetNames) {
  return(list2env(sapply(probeSetNames, function(probeSetName, idealIntensities, fileName) {
    return(list2env(sapply(ls(idealIntensities[[fileName]][[probeSetName]]), function(index, idealIntensities, fileName, probeSetName) {
      PM <- idealIntensities[[fileName]][[probeSetName]][[index]][["PM"]]
      IM <- idealIntensities[[fileName]][[probeSetName]][[index]][["IM"]]
      return(log2(max(PM - IM, d)))
    }, idealIntensities, fileName, probeSetName, simplify = F)))
  }, idealIntensities, fileName, simplify = F)))
}, idealIntensities, probeSetNames, simplify = F))

#now get the signal log value, which is the absolute expression value for the probe set, as the one-step biweight estimate of the probe values
signalLogValues <- list2env(sapply(fileNames, function(fileName, probeValues, probeSetNames) {
  return(list2env(sapply(probeSetNames, function(probeSetName, probeValues, fileName) {
    xVec <- unname(sapply(ls(probeValues[[fileName]][[probeSetName]]), function(index, probeValues, fileName, probeSetName) {
      probeValues[[fileName]][[probeSetName]][[index]]
    }, probeValues, fileName, probeSetName))
    return(oneStepBiweight(xVec, C, e))
  }, probeValues, fileName, simplify = F)))
}, probeValues, probeSetNames, simplify = F))

signalLogValuesList <- sapply(fileNames, function(fileName, signalLogValues, probeSetNames) {
  return(sapply(probeSetNames, function(probeSetName, signalLogValues, fileName) {
    return(signalLogValues[[fileName]][[probeSetName]])
  }, signalLogValues, fileName))
}, signalLogValues, probeSetNames, simplify = F)

saveRDS(signalLogValuesList, 'signalLogValues.Rdata')




