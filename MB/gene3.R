gene3<-function (exprMat, scenicOptions,finished = 1, nParts = 10, resumePreviousRun = FALSE, 
          allTFs = getDbTfs(scenicOptions), ...) 
{
  nCores <- getSettings(scenicOptions, "nCores")
  if (is.data.frame(exprMat)) {
    supportedClasses <- paste(gsub("AUCell_buildRankings,", 
                                   "", methods("AUCell_buildRankings")), collapse = ", ")
    supportedClasses <- gsub("-method", "", supportedClasses)
    stop("'exprMat' should be one of the following classes: ", 
         supportedClasses, "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
  }
  if (any(table(rownames(exprMat)) > 1)) 
    stop("The rownames (gene id/name) in the expression matrix should be unique.")
  inputTFs <- allTFs[allTFs %in% rownames(exprMat)]
  percMatched <- length(inputTFs)/length(allTFs)
  if (getSettings(scenicOptions, "verbose")) 
    message("Using ", length(inputTFs), " TFs as potential regulators...")
  if (percMatched < 0.4) 
    warning("Only ", length(inputTFs), " (", round(percMatched * 
                                                     100), "%) of the ", length(allTFs), " TFs in the database were found in the dataset. Do they use the same gene IDs?\n")
  weightMatrices <- list()
  if (resumePreviousRun) {
    if (file.exists(getIntName(scenicOptions, "genie3ll"))) {
      stop("The previous run already finished running (the output already exists: '", 
           getIntName(scenicOptions, "genie3ll"), "'). \nTo re-run GENIE3 set 'resumePreviousRun=FALSE'.")
    }
    else {
      fileNames <- gsub(".Rds$", "_part_", getIntName(scenicOptions, 
                                                      "genie3wm"))
      intDir <- dirname(fileNames)
      fileNames <- list.files(dirname(fileNames), pattern = basename(fileNames))
      if (length(fileNames) > 0) {
        for (fileName in fileNames) {
          i <- gsub(".Rds$", "", strsplit(fileName, 
                                          "_part_")[[1]][[2]])
          weightMatrices[[i]] <- readRDS(file.path(intDir, 
                                                   fileName))
        }
      }
    }
  }
  genesDone <- unname(unlist(lapply(weightMatrices, colnames)))
  if (any(table(genesDone) > 1)) 
    stop("Some genes are in several of the partial runs.")
  genesLeft <- setdiff(sort(rownames(exprMat)), genesDone)
  if (length(genesLeft) > 0) {
    partNames <- as.character(seq_len(nParts))
    partNames <- setdiff(partNames, names(weightMatrices))
    if (length(partNames) == 0) {
      warning("Splitting the ", length(genesLeft), " genes left into ", 
              nParts, " parts.")
      partNames <- as.character(max(as.numeric(names(weightMatrices))) + 
                                  seq_len(nParts))
    }
    genesSplit <- suppressWarnings(split(genesLeft, partNames))
    print( names(genesSplit)[finished:nParts] )
    
    for (i in names(genesSplit)[finished:nParts]) {
      if (getSettings(scenicOptions, "verbose")) 
        message("Running GENIE3 part ", i)
      set.seed(getSettings(scenicOptions, "seed"))
      weightMatrix <- GENIE3::GENIE3(exprMat, regulators = inputTFs, 
                                     nCores = nCores, targets = genesSplit[[i]], 
                                     ...)
      fileName <- gsub(".Rds$", paste0("_part_", i, ".Rds"), 
                       getIntName(scenicOptions, "genie3wm"))
      saveRDS(weightMatrix, file = fileName)
      weightMatrices[[i]] <- weightMatrix
      closeAllConnections()
    }
  }
  linkList_list <- list()
  for (i in names(weightMatrices)) {
    weightMatrix <- weightMatrices[[i]]
    linkList_list[[i]] <- GENIE3::getLinkList(weightMatrix, 
                                              threshold = getSettings(scenicOptions, "modules/weightThreshold"))
    
  }
  rm(weightMatrices)
  linkList <- do.call(rbind, linkList_list)
  colnames(linkList) <- c("TF", "Target", "weight")
  linkList <- linkList[order(linkList[, "weight"], decreasing = TRUE), 
  ]
  linkList <- unique(linkList)
  rownames(linkList) <- NULL
  saveRDS(linkList, file = getIntName(scenicOptions, "genie3ll"))
 
  if (getSettings(scenicOptions, "verbose")) 
    message("Finished running GENIE3.")
  invisible(linkList)
}
