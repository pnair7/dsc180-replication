require(R.matlab)
require(analyzeMC)


# directory containing .mat files
###############
### ADJUST ####
###############
dataDir <- "C:/Users/Daniel Tong/Desktop/school stuff/DSC/DSC 180A/BM/_mat/"

# directory to save .rds files
###############
### ADJUST ####
###############
RdataDir <- "C:/Users/Daniel Tong/Desktop/school stuff/DSC/DSC 180A/BM/_rds/inhibitor-data/"

# which data sets are processed?
dataSets <- "BM" # "BND" or "BM"
experiments <- "inhibitor" # "8donor", "inhibitor", or "time"

# list files
dataFiles <- list.files(dataDir)
if(length(grep("info.mat", dataFiles)) > 0){
  dataFiles <- dataFiles[-grep("info.mat", dataFiles)]
}
if(length(grep("inhSampleSizes.mat", dataFiles)) > 0){
  dataFiles <- dataFiles[-grep("inhSampleSizes.mat", dataFiles)]
}

# find unique cell types / define grouping
uniqueCellTypes <- unique(sapply(dataFiles, function(i){

  if(dataSets == "BND"){
    tmp<-strsplit(i, "_")[[1]]
    paste(tmp[-c(length(tmp)-1, length(tmp))], collapse = "_")
  }else if(dataSets == "BM"){
    tmp<-strsplit(i, "_")[[1]]
    paste(tmp[length(tmp)-1], collapse = "_")
  }else{
    stop(paste("dataSets", dataSets, "not supported"))
  }

  }))

nCellTypes <- length(uniqueCellTypes)


#print(readMat('C:/Users/Daniel Tong/Desktop/school stuff/DSC/DSC 180A/BM/_mat/Akti.mat'))


# save one .rds file for each cell type / grouping
for(j in seq_along(uniqueCellTypes)){
  uniqueCellTypeJ <- paste(strsplit(uniqueCellTypes[j], "\\+")[[1]], collapse = "\\+")
  dataFilesCellType <- dataFiles[grep(uniqueCellTypeJ, dataFiles)]
  
  
  p <- ncol(as.matrix(readMat(paste(dataDir, dataFilesCellType[1], sep = ""))$dataset[[1]]))

  
  X <- matrix(nrow=0,ncol=p)
  
  
  ExpInd <- numeric(0)
  remove <- numeric(0)

  
  for(i in seq_along(dataFilesCellType)){
    tmpDat <- readMat(paste(dataDir, dataFilesCellType[i], sep = ""))$dataset
    designMat <- as.matrix(tmpDat[[1]])

    if(nrow(designMat) == 0){
      warning(paste(dataFilesCellType[i], "is empty!"))
      remove <- c(remove, i)
    }else{
      varNames <- unlist(tmpDat[2])
      colnames(designMat) <- varNames
      X <- rbind(X, designMat)
      ExpInd <- c(ExpInd, rep(i, nrow(designMat)))
    }
  }

  
  if(length(remove) > 0) dataFilesCellType <- dataFilesCellType[-remove]
  ExpIndDescription <- data.frame(ExpInd = 1:length(dataFilesCellType),
                                  Description = dataFilesCellType)

  ExpIndDescription <- parseDescription(ExpIndDescription, dataSets, experiments)

  ## save as rds file in _rds folder
  saveRDS(list(X = X, ExpInd = ExpInd, ExpIndDescription = ExpIndDescription),
          file = paste(RdataDir, uniqueCellTypes[j],
                       ifelse(dataSets == "BM", paste("_", experiments), ""),
                       "_", dataSets,".rds", sep = ""))
  print(uniqueCellTypes[j])
}
