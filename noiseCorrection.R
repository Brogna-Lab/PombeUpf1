probeData <- read.csv('ProbeDataTable.csv', stringsAsFactors = F)

fileNames <- read.table('celFilePaths.txt', stringsAsFactors = F)[, 1]
fileNames <- strsplit(fileNames, '/')
fileNames <- unlist(lapply(fileNames, function(entry) {return(entry[2])}))
fileNames <- unique(unlist(strsplit(fileNames, '.CEL')))
fileNames <- gsub('-', '_', fileNames)


#split the arrays into k (default k=16) rectangular zones
k <- 16
nRows <- max(max(probeData$PMX), max(probeData$MMX))
nCols <- max(max(probeData$PMY), max(probeData$MMY))
zoneRows <- nRows/sqrt(k)
zoneCols <- nCols/sqrt(k)
rowIndices <- c((1:sqrt(k) * zoneRows - zoneRows + 1), nRows + 1)
colIndices <- c(sapply(1:sqrt(k), function(index, zoneCols) {if (index %% 2 != 0) {return(floor(zoneCols * index - zoneCols + 1))} else {return(ceiling(zoneCols * index - zoneCols + 1))}}, zoneCols), nCols + 1)
xDims <- sapply(1:sqrt(k), function(index, rowIndices) {return(rowIndices[index + 1] - rowIndices[index])}, rowIndices) #these just show that the rectangles have the correct dimensions
yDims <- sapply(1:sqrt(k), function(index, colIndices) {return(colIndices[index + 1] - colIndices[index])}, colIndices)


#get the zone backgrounds (bZ) and background variabilities (nZ) for every zone by getting the bottom 2% for values in each zone
#sort through the table and put the values into a list named based on the x/y coordinate
orderedProbeData <- probeData[order(probeData$PMY, probeData$PMX), ]


probeDataHash <- lapply(colnames(orderedProbeData), function(colName, orderedProbeData) {
  result <- as.list(orderedProbeData[, colName])
  names(result) <- as.character(1:length(result))
  return(list2env(result))
  }, orderedProbeData)
names(probeDataHash) <- colnames(probeData)
#probeDataHash <- list2env(probeDataHash)


#make hash that has all of the values for every file for each pair of coordinates, first level is all pairs of x and y coords, second is 13 files per pair of coordinates (for each of the to x coords)
fileNamesHash <- sapply(fileNames, function(fileName) {
  return(list('PM' = paste('PM', fileName, sep = ''), 'MM' = paste('MM', fileName, sep = '')))
}, simplify = F)



listNames <- c(paste(orderedProbeData[, "PMX"], orderedProbeData[, "PMY"], sep = "-"), paste(orderedProbeData[, "MMX"], orderedProbeData[, "MMY"], sep = "-"))
# library(unix)
# rlimit_as(24e+10)
valuesHash <- list2env(sapply(listNames, function(coordPair, fileNames) {
  filesList <- vector('list', length(fileNames))
  names(filesList) <- fileNames
  return(NULL)
}, fileNames, simplify = F))


for (index in 1:length(probeDataHash[["PMX"]])) {
  key <- paste(index)
  PMX <- probeDataHash[["PMX"]][[key]]
  MMX <- probeDataHash[["MMX"]][[key]]
  Y <- probeDataHash[["PMY"]][[key]]
  PMcoordPair <- paste(PMX, Y, sep = '-')
  MMcoordPair <- paste(MMX, Y, sep = '-')
  for (fileName in fileNames) {
    valuesHash[[PMcoordPair]][[fileName]] <- probeDataHash[[fileNamesHash[[fileName]][["PM"]]]][[key]]
    valuesHash[[MMcoordPair]][[fileName]] <- probeDataHash[[fileNamesHash[[fileName]][["MM"]]]][[key]]
  }
}

allCoords <- list2env(sapply(ls(valuesHash), function(coordPair) {
  coordValues<- as.numeric(strsplit(coordPair, '-')[[1]])
  names(coordValues) <- c('x', 'y')
  return(coordValues)
}, simplify = F))
allCoordsKey <- ls(allCoords)
zoneCoords <- new.env() #get the coordinates corresponding to each zone
for (rowIndex in 1:sqrt(k)) {
  xMin <- rowIndices[rowIndex]
  xMax <- rowIndices[rowIndex + 1] - 1
  for (colIndex in 1:sqrt(k)) {
    currentCoords <- new.env()
    counter <- 1
    yMin <- colIndices[colIndex]
    yMax <- colIndices[colIndex + 1] - 1
    for (coordPair in allCoordsKey) {
      xCoord <- allCoords[[coordPair]][['x']]
      yCoord <- allCoords[[coordPair]][['y']]
      if (xMin <= xCoord && xCoord <= xMax && yMin <= yCoord && yCoord <= yMax) {
        currentCoords[[paste(counter)]] <- coordPair
        counter <- counter + 1
      }
    }
    zoneCoords[[paste(rowIndex, colIndex, sep = '-')]] <- unname(sapply(ls(currentCoords), function(key, currentCoords) {return(currentCoords[[key]])}, currentCoords))
  }
}

#find the centre point for each zone
zoneCentres <- new.env()
for (rowIndex in 1:sqrt(k)) {
  xCoordinate <- (rowIndices[rowIndex] + rowIndices[rowIndex + 1] - 1) / 2
  for (colIndex in 1:sqrt(k)) {
    yCoordinate <- (colIndices[colIndex] + colIndices[colIndex + 1] - 1) / 2
    zoneCentres[[paste(rowIndex, colIndex, sep = '-')]] <- list(x = xCoordinate, y = yCoordinate)
  }
}

#find the distance between every coordinate and every zone's centre point
twoDimensionalDistance <- function(coordPair1, coordPair2) { #find the 2d euclidean distance between 2 points using their x and y coordinates
  return(sqrt((coordPair1[["x"]] - coordPair2[["x"]])^2 + (coordPair1[["y"]] - coordPair2[["y"]])^2))
}
#also find the weighting that each zone has on each coord based on distance to centre for removing background noise
distances <- new.env()
weightings <- new.env()
smooth <- 100 #default value to ensure distance value is never 0 when we are doing the dividing
zones <- as.vector(sapply(1:sqrt(k), function(i, k) {return(sapply(1:sqrt(k), function(j, i) {return(paste(i, j, sep = '-'))}, i))}, k))
for (zone in zones) {
  distances[[zone]] <- new.env()
  weightings[[zone]] <- new.env()
  for (coord in allCoordsKey) {
    distances[[zone]][[coord]] <- twoDimensionalDistance(zoneCentres[[zone]], allCoords[[coord]])
    weightings[[zone]][[coord]] <- 1/(distances[[zone]][[coord]] + smooth)
  }
}


allZoneValues <- new.env()
bZ <- new.env()
nZ <- new.env()
for (zone in zones) {
  allZoneValues[[zone]] <- new.env()
  bZ[[zone]] <- new.env()
  nZ[[zone]] <- new.env()
  for (fileName in fileNames) {
    zoneValues <- sapply(zoneCoords[[zone]], function(coords, valuesHash, fileName) {return(valuesHash[[coords]][[fileName]])}, valuesHash, fileName)
    allZoneValues[[zone]][[fileName]] <- list2env(as.list(zoneValues))
    orderedValues <- unname(zoneValues[order(zoneValues)])
    bottom2percent <- orderedValues[1:round(length(orderedValues) * 0.02)]
    bZ[[zone]][[fileName]] <- mean(bottom2percent)
    nZ[[zone]][[fileName]] <- sd(bottom2percent)
  }
}


#now get the background and noise values for each cell using the weighted values of the bottom 2% mean and sd of each zone
bZsVecs <- list2env(sapply(fileNames, function(fileName, zones, bZ) {return(unname(sapply(zones, function(zone, bZ, fileName) {return(bZ[[zone]][[fileName]])}, bZ, fileName)))}, zones, bZ, simplify = F))
nZsVecs <- list2env(sapply(fileNames, function(fileName, zones, nZ) {return(unname(sapply(zones, function(zone, nZ, fileName) {return(nZ[[zone]][[fileName]])}, nZ, fileName)))}, zones, nZ, simplify = F))
zoneWeightsVecs <- list2env(sapply(allCoordsKey, function(coord, zones, weightings) {return(unname(sapply(zones, function(zone, coord, weightings) {return(weightings[[zone]][[coord]])}, coord, weightings)))}, zones, weightings, simplify = F))
zoneWeightsSums <- list2env(sapply(allCoordsKey, function(coord, zoneWeightsVecs) {return(sum(zoneWeightsVecs[[coord]]))}, zoneWeightsVecs, simplify = F))

backgrounds <- list2env(sapply(fileNames, function(fileName) {return(new.env())}))
noises <- list2env(sapply(fileNames, function(fileName) {return(new.env())}))
for (fileName in fileNames) {
  for (coord in allCoordsKey) {
    backgrounds[[fileName]][[coord]] <- sum(zoneWeightsVecs[[coord]] * bZsVecs[[fileName]]) / zoneWeightsSums[[coord]]
    noises[[fileName]][[coord]] <- sum(zoneWeightsVecs[[coord]] * nZsVecs[[fileName]]) / zoneWeightsSums[[coord]]
  }
}


#now correct each intensity value for noise getting adjusted values by shifting the intensities down by the local background, howwever must ensure that they do not become
#negative since we want to do calculations on log intensity values, so using floor values based on the variation of the local background

adjustedValues <- list2env(sapply(fileNames, function(fileName) {return(new.env())}))
floorVal <- 0.5 #minimum possible value the intensity could be, to ensure it isnt 0
noiseFrac <- 0.5 #what fraction of local noise we ant to set the threshold to for saying our adjusted intensity has to be at least as big as that fraction of the noise
for (fileName in fileNames) {
  for (coord in allCoordsKey) {
    adjustedValues[[fileName]][[coord]] <- max(max(valuesHash[[coord]][[fileName]], floorVal) - backgrounds[[fileName]][[coord]], noises[[fileName]][[coord]] * noiseFrac)
  }
}

adjustedValuesList <- sapply(fileNames, function(fileName, adjustedValues, allCoordsKey) {
  return(sapply(allCoordsKey, function(coord, fileName, adjustedValues) {
    return(adjustedValues[[fileName]][[coord]])
  }, fileName, adjustedValues, simplify = F))
}, adjustedValues, allCoordsKey, simplify = F)

saveRDS(adjustedValuesList, 'adjustedIntensities.Rdata')



