# sub by species
noID <- function(inData, bioLevel){
  x <- inData[ , bioLevel]
  x <- as.data.frame(table(x))
  
  if("" %in% unique(x$x)){
    numnoID <- x$Freq[1]
  }else{
    numnoID <- 0
  }
  
  return(numnoID)
}


#### process section #### 
generateGrid <- function(occurData, latExt, longExt, resGrid){
  # make ybins (latbins)
  # a sequence of the min y value for each grid
  longBin1 <- seq(from = longExt[1], 
                  to = (longExt[2] - resGrid[1]),
                  by = resGrid[1])
  
  # a sequence of the max y value for each grid
  longBin2 <- seq(from = (longExt[1] + resGrid[1]), 
                  to = longExt[2],
                  by = resGrid[1])
  
  # combine to form a data frame for the min and max y values of the grids
  longBin <- data.frame(minX = longBin1, maxX = longBin2)
  
  # do the same thing for the x
  # make xbins (latbins)
  # a sequence of the min y value for each grid
  latBin1 <- seq(from = latExt[1], 
                 to = (latExt[2] - resGrid[2]),
                 by = resGrid[2])
  
  # a sequence of the max y value for each grid
  latBin2 <- seq(from = (latExt[1] + resGrid[2]), 
                 to = latExt[2],
                 by = resGrid[2])
  
  # combine to form a data frame for the min and max y values of the grids
  latBin <- data.frame(minY = latBin1, maxY = latBin2)
  
  
  # do a nested loop like Cynthia mentioned
  # loop1: iterate through the rows of the longBin
  # loop2: iterate through the rows of the latBin
  
  # before doing the loop; let's make data frame that will
  # hold the results
  # it will be good to keep track of these information
  # 1) miny, 2) maxy, 3) minx, 4) maxx, 5) species richness 6) gridID
  # 7) totalPoints
  # these 7 info will be the columns of the results
  
  # create dataframe for the results
  # this data frame will be filled during the loop
  resultData <- data.frame(minx = numeric(), maxx = numeric(),
                           miny = numeric(), maxy = numeric(),
                           speciesrichness = integer(), gridID= integer(),
                           totalPoints = integer())
  
  # resultData data frame row indicator
  resInd <- 1
  
  # loop1
  # note: long will be used as row index for the longBin
  for(long in 1:nrow(longBin)){
    
    # loop2
    # note: lat will be used as row index for the latBin
    for(lat in 1:nrow(latBin)){
      
      # subset the occurrence data using the following criteria
      # the latitude of the point should be greater than or equal
      # the min lat of the grid AND should be less than or equal
      # the max lat of the grid AND the longitude of the point should
      # be greater than or equal the min long of the grid AND should be
      # less than or equal the max long of the grid
      subTemp <- occurData[ ( (occurData$lon >= longBin[long, "minX"]) &  
                                (occurData$lon <= longBin[long, "maxX"]) ) &
                              ( (occurData$lat >= latBin[lat, "minY"]) &  
                                  (occurData$lat <= latBin[lat, "maxY"]) ), ]
      
      # this shows the number of species: length(unique(subTemp$species))
      # this shows the number of points: nrow(subTemp)
      
      # concatenates the results together and put them in the resultData
      # min x: longBin[long , "minX"]
      # max x: longBin[long , "maxX"]
      # min y: latBin[lat, "minY"]
      # max y: latBin[lat, "maxY"]
      # number of species: length(unique(subTemp$species))
      # a grid ID: paste(long, lat, sep = "_")
      # number of points in grid: nrow(subTemp)
      resultData[resInd, ] <- c(longBin[long , "minX"], longBin[long , "maxX"], 
                                latBin[lat, "minY"], latBin[lat, "maxY"], 
                                length(unique(subTemp$species)), 
                                paste(long, lat, sep = "_"), nrow(subTemp))
      
      # add a count for the resultData row
      resInd <- resInd + 1 
      
    }
    
    
  }
  
  # convert the columns into numeric;
  # will be good to represent each point as the median of the min and max of
  # the coordinates (e.g., (resultData$miny + resultData$maxy)/2))
  resultData[ , c(1:5, 7)] <- apply(resultData[ , c(1:5, 7)], 2, function(x) as.numeric(x))
  
  
  #### convert grids to raster ####
  
  # create a blank raster
  # this uses the extent, and resolution provided above
  blankRaster <- raster::raster(nrows = (latExt[2] - latExt[1])/resGrid[2], 
                                ncols = (longExt[2] - longExt[1])/resGrid[1], 
                                xmn = longExt[1], xmx = longExt[2], 
                                ymn = latExt[1], ymx = latExt[2])
  
  # use the midpoint between min and max as center of grid
  resultData$x <- (resultData$minx + resultData$maxx)/2
  resultData$y <- (resultData$miny + resultData$maxy)/2
  
  # create a spatial points for the result with the resultDate as a dataframe attached
  # to the spatial points
  resPoints <- sp::SpatialPointsDataFrame(cbind(resultData$x, resultData$y), data = resultData)
  
  # generate the raster with 
  rasterData <- raster::rasterize(resPoints, blankRaster, field = resPoints$speciesrichness)
  
  # check how the output looks like
  return(rasterData)
}



gettaxaInfo <- function(genbankFile){
  
  taxonomicInfo <- genbankFile[(grep(x = genbankFile, pattern = "ORGANISM")):((grep(x = genbankFile, pattern = "REFERENCE")[1]) - 1)]
  
  taxonomicInfo <- gsub(x = taxonomicInfo, 
                        pattern = "            |  ORGANISM  |\\."  ,
                        replacement = "")
  
  taxonomicInfo <- paste(taxonomicInfo, collapse = ";")
  
  taxonomicInfo <- gsub(x = taxonomicInfo, 
                        pattern = ";;", 
                        replacement = ";")
  
  taxonomicInfo <- gsub(x = taxonomicInfo, 
                        pattern = " ", 
                        replacement = "")
  
  return(taxonomicInfo)
  
}

getauthorInfo <- function(genbankFile){
  
  authorInfo <- genbankFile[(grep(x = genbankFile, pattern = "^  AUTHORS|^  CONSRTM")[1]):(grep(x = genbankFile, pattern = "TITLE")[1] - 1)]
  
  authorInfo <- gsub(x = authorInfo, pattern = "  AUTHORS   |  CONSRTM   ", replacement = "")
  
  authorInfo <- gsub(x = authorInfo, pattern = "            ", replacement = "")
  
  authorInfo <- paste(authorInfo, collapse = "")
  
  authorInfo <- gsub(x = authorInfo, pattern = " and", replacement = ",")
  
  authorInfo <- gsub(x = authorInfo, pattern = ", ", replacement = ";")
  
  return(authorInfo)
  
}


getinstiInfo <- function(genbankFile){
  
  startLine <- grep(x = genbankFile, pattern = "JOURNAL")
  
  endLine <- grep(x = genbankFile, pattern = "FEATURES|REMARK|^COMMENT")
  
  if(length(startLine) > 1){
    startLine <- startLine[2]
    
    endLine <- endLine[which(endLine > startLine)][1]
    
  }else{
    startLine <- startLine[1]
    
    endLine <- endLine[which(endLine > startLine)][1]
    
  }
  
  instiInfo <- genbankFile[startLine:(endLine - 1)]
  
  instiInfo <- gsub(x = instiInfo,
                    patter = "            |  JOURNAL   |",
                    replacement = "")
  
  instiInfo <- paste(instiInfo, collapse = "")
  
  instiInfo <- gsub(x = instiInfo, pattern = "^.*) ", rep = "")
  
  return(instiInfo)
}

getmetaInfo <- function(genbankFile){
  
  genbankFile <- samp
  
  # tempSource <- genbankFile[(grep(x = genbankFile, pattern = "^     source") + 1):
  #                             (grep(x = genbankFile, pattern = "^     gene|^     rRNA|^     tRNA|^     CDS|^     repeat_region")[1] - 1)]
  
  stLine <- (grep(x = genbankFile, pattern = "^     source") + 1)
  
  if(is.na((grep(x = genbankFile, pattern = "^     [a-z|A-Z]")[2] - 1))){
    
    enLine <- (grep(x = genbankFile, pattern = "^ORIGIN") - 1)
  }else{
    
    enLine <- (grep(x = genbankFile, pattern = "^     [a-z|A-Z]")[2] - 1)
    
  }
  
  tempSource <- genbankFile[stLine:enLine]
  
  metaData <- grep(x = tempSource, pattern = "                    /")
  
  combiInd <- which(!(1:length(tempSource) %in% metaData))
  
  
  if(!(length(combiInd) == 0)){
    
    tempSource[combiInd - 1] <- paste(tempSource[combiInd - 1], tempSource[combiInd], sep = "")
    
    tempSource <- tempSource[-combiInd]
    
  }
  
  tempSource <- gsub(x = tempSource, pattern = "                     ", replacement = "")
  
  headerMeta <- paste(gsub(x = tempSource, pattern = "=.*$|^/", replacement = ""), collapse = ":")
  
  headerMeta <- gsub(x = headerMeta, pattern = " ", replacement = "_")
  
  infoMeta <- paste(gsub(x = tempSource, pattern = '^.*=|\\"', replacement = ""), collapse = "::")
  
  infoMeta <- gsub(x = infoMeta, pattern = " ", replacement = "_")
  
  return(c(headerMeta, infoMeta))
  
  
}

# get sequence

getseqInfo <- function(genbankFile){
  
  temp <- grep(x = genbankFile, pattern = "ORIGIN|//")
  
  seqInfo <- genbankFile[(temp[1] + 1):(temp[2] - 1)]
  
  rm(temp)
  
  seqInfo <- paste(gsub(x = seqInfo, pattern = "[[:digit:]]| ", 
                        replacement = ""), 
                   collapse = "")
  
  return(seqInfo)
  
}

# get gene marker
getgeneInfo <- function(genbankFile){
  
  # temp <- genbankFile[grep("/gene|/product", genbankFile)[1]]
  
  if(length(grep("/gene", genbankFile)) == 0){
    
    if(length(grep("/product", genbankFile))> 1){
      
      temp <- genbankFile[grep('/product.*subunit I"|/product.*COX1|/product.*Cox1|/product.*coi"|/product.*co1|/product.*CO-I|/product.*CO-1|/product.*CO1|/product.*coI|/product.*cox1', genbankFile)[1]]
      # '/product.*subunit I"|COX1|Cox1|coi"|co1|CO-I|CO-1|CO1|coI|cox1'
      # '/gene.*subunit I"|/gene.*COX1|/gene.*Cox1|/gene.*coi$|/gene.*co1|/gene.*CO-I|/gene.*CO-1|/gene.*CO1|/gene.*coI|/gene.*cox1'
    }else{
      
      temp <- genbankFile[grep("/product", genbankFile)[1]]
      
    }
    
    
  }else{
    
    if(length(grep("/gene", genbankFile))> 1){
      
      if(length(grep('/gene.*subunit I"|/gene.*COX1|/gene.*Cox1|/gene.*coi$|/gene.*co1|/gene.*CO-I|/gene.*CO-1|/gene.*CO1|/gene.*coI|/gene.*cox1', genbankFile)) == 0){
        
        temp <- genbankFile[grep("/gene", genbankFile)[1]]
        
      }else{
        
        
        temp <- genbankFile[grep('/gene.*subunit I"|/gene.*COX1|/gene.*Cox1|/gene.*coi$|/gene.*co1|/gene.*CO-I|/gene.*CO-1|/gene.*CO1|/gene.*coI|/gene.*cox1', genbankFile)[1]]  
        
      }
      
      
      
    }else{
      
      temp <- genbankFile[grep("/gene", genbankFile)[1]]
      
    }
    # temp <- genbankFile[grep("/gene", genbankFile)[1]]
    
  }
  
  out <- gsub('/gene=|/product=| |\\"', "", temp)
  
  return(out)
}