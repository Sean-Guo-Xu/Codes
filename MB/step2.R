step2<-function (scenicOptions, minGenes = 20, coexMethods = NULL, 
          minJakkardInd = 0.8, signifGenesMethod = "aprox", onlyPositiveCorr = TRUE, 
          onlyBestGsPerMotif = TRUE, dbIndexCol = "features") 
{
  nCores <- getSettings(scenicOptions, "nCores")
  tfModules_asDF <- tryCatch(loadInt(scenicOptions, "tfModules_asDF"), 
                             error = function(e) {
                               if (getStatus(scenicOptions, asID = TRUE) < 2) 
                                 e$message <- paste0("It seems the co-expression modules have not been built yet. Please, run runSCENIC_1_coexNetwork2modules() first.\n", 
                                                     e$message)
                               stop(e)
                             })
  if (!is.null(coexMethods)) 
    tfModules_asDF <- tfModules_asDF[which(tfModules_asDF$method %in% 
                                             coexMethods), ]
  if (!is.null(minJakkardInd)) 
    tfModules_asDF <- mergeOverlappingModules(tfModules_asDF, 
                                              minJakkardInd = minJakkardInd)
  if (nrow(tfModules_asDF) == 0) 
    stop("The co-expression modules are empty.")
  if ("BiocParallel" %in% installed.packages() && (nCores > 
                                                   1)) {
    library(BiocParallel)
    register(MulticoreParam(nCores), default = TRUE)
  }
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tStep 2. Identifying regulons")
  if (getSettings(scenicOptions, "verbose")) 
    message(msg)
  library(AUCell)
  library(RcisTarget)
  motifAnnot <- getDbAnnotations(scenicOptions)
  if (is.null(names(getSettings(scenicOptions, "dbs")))) {
    names(scenicOptions@settings$dbs) <- scenicOptions@settings$dbs
    tmp <- sapply(strsplit(getSettings(scenicOptions, "dbs"), 
                           "-", fixed = T), function(x) x[grep("bp|kb", x)])
    if (all(lengths(tmp) > 0)) 
      names(scenicOptions@settings$dbs) <- tmp
  }
  loadAttempt <- sapply(getDatabases(scenicOptions), dbLoadingAttempt, 
                        indexCol = dbIndexCol)
  if (any(!loadAttempt)) 
    stop("It is not possible to load the following databses: \n", 
         paste(dbs[which(!loadAttempt)], collapse = "\n"))
  genesInDb <- unique(unlist(lapply(getDatabases(scenicOptions), 
                                    function(dbFilePath) {
                                      rf <- arrow::ReadableFile$create(dbFilePath)
                                      fr <- arrow::FeatherReader$create(rf)
                                      genesInDb <- names(fr)
                                      rnktype <- "features"
                                      genesInDb <- genesInDb[genesInDb != rnktype]
                                    })))
  featuresWithAnnot <- checkAnnots(scenicOptions, motifAnnot)
  if (any(featuresWithAnnot == 0)) 
    warning("Missing annotations\n", names(which(rankingsInDb == 
                                                   0)))
  tfModules_asDF$TF <- as.character(tfModules_asDF$TF)
  tfModules_asDF$Target <- as.character(tfModules_asDF$Target)
  allTFs <- getDbTfs(scenicOptions)
  tfModules_asDF <- tfModules_asDF[which(tfModules_asDF$TF %in% 
                                           allTFs), ]
  geneInDb <- tfModules_asDF$Target %in% genesInDb
  missingGene <- sort(unique(tfModules_asDF[which(!geneInDb), 
                                            "Target"]))
  if (length(missingGene) > 0) 
    warning(paste0("Genes in co-expression modules not available in RcisTargetDatabases: ", 
                   paste(missingGene, collapse = ", ")))
  tfModules_asDF <- tfModules_asDF[which(geneInDb), ]
  if (all(is.na(tfModules_asDF$corr))) {
    warning("no correlation info available")
    tfModules_Selected <- tfModules_asDF
    tfModules_Selected$geneSetName <- paste(tfModules_Selected$TF, 
                                            tfModules_Selected$method, sep = "_")
  }
  else {
    tfModules_Selected <- tfModules_asDF[which(tfModules_asDF$corr == 
                                                 1), ]
    tfModules_Selected$geneSetName <- paste(tfModules_Selected$TF, 
                                            tfModules_Selected$method, sep = "_")
    if (!onlyPositiveCorr) {
      tfModules_IgnCorr <- tfModules_asDF[which(tfModules_asDF$corr != 
                                                  1), ]
      tfModules_IgnCorr$geneSetName <- paste0(tfModules_IgnCorr$TF, 
                                              "_", tfModules_IgnCorr$method)
      posCorr <- tfModules_Selected[which(tfModules_Selected$geneSetName %in% 
                                            unique(tfModules_IgnCorr$geneSetName)), ]
      tfModules_IgnCorr <- rbind(tfModules_IgnCorr, posCorr)
      tfModules_IgnCorr$geneSetName <- paste0(tfModules_IgnCorr$geneSetName, 
                                              "IgnCorr")
      tfModules_Selected <- rbind(tfModules_Selected, 
                                  tfModules_IgnCorr)
    }
  }
  tfModules_Selected$geneSetName <- factor(as.character(tfModules_Selected$geneSetName))
  allGenes <- unique(tfModules_Selected$Target)
  tfModules <- split(tfModules_Selected$Target, tfModules_Selected$geneSetName)
  tfModules <- setNames(lapply(names(tfModules), function(gsn) {
    tf <- strsplit(gsn, "_")[[1]][1]
    unique(c(tf, tfModules[[gsn]]))
  }), names(tfModules))
  tfModules <- tfModules[which(lengths(tfModules) >= minGenes)]
  saveRDS(tfModules, file = getIntName(scenicOptions, "tfModules_forEnrichment"))
  if (getSettings(scenicOptions, "verbose")) {
    tfModulesSummary <- t(sapply(strsplit(names(tfModules), 
                                          "_"), function(x) x[1:2]))
    message("tfModulesSummary:")
    print(cbind(sort(table(tfModulesSummary[, 2]))))
  }
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Calculating AUC")
  if (getSettings(scenicOptions, "verbose")) 
    message(msg)
  motifs_AUC <- lapply(getDatabases(scenicOptions), function(rnkName) {
    ranking <- importRankings(rnkName, columns = allGenes)
    message("Scoring database: ", ranking@description)
    RcisTarget::calcAUC(tfModules, ranking, aucMaxRank = 0.03 * 
                          getNumColsInDB(ranking), nCores = nCores, verbose = FALSE)
  })
  saveRDS(motifs_AUC, file = getIntName(scenicOptions, "motifs_AUC"))
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Adding motif annotation")
  message(msg)
  motifEnrichment <- lapply(motifs_AUC, function(aucOutput) {
    tf <- sapply(setNames(strsplit(rownames(aucOutput), 
                                   "_"), rownames(aucOutput)), function(x) x[[1]])
    addMotifAnnotation(aucOutput, nesThreshold = 3, digits = 3, 
                       motifAnnot = motifAnnot, motifAnnot_highConfCat = c("directAnnotation", 
                                                                           "inferredBy_Orthology"), motifAnnot_lowConfCat = c("inferredBy_MotifSimilarity", 
                                                                                                                              "inferredBy_MotifSimilarity_n_Orthology"), highlightTFs = tf)
  })
  motifEnrichment <- do.call(rbind, lapply(names(motifEnrichment), 
                                           function(dbName) {
                                             cbind(motifDb = dbName, motifEnrichment[[dbName]])
                                           }))
  saveRDS(motifEnrichment, file = getIntName(scenicOptions, 
                                             "motifEnrichment_full"))
  msg <- paste0("Number of motifs in the initial enrichment: ", 
                nrow(motifEnrichment))
  if (getSettings(scenicOptions, "verbose")) 
    message(msg)
  motifEnrichment_selfMotifs <- motifEnrichment[which(motifEnrichment$TFinDB != 
                                                        ""), , drop = FALSE]
  msg <- paste0("Number of motifs annotated to the matching TF: ", 
                nrow(motifEnrichment_selfMotifs))
  if (getSettings(scenicOptions, "verbose")) 
    message(msg)
  rm(motifEnrichment)
  if (nrow(motifEnrichment_selfMotifs) == 0) 
    stop("None of the co-expression modules present enrichment of the TF motif: There are no regulons.")
  if (onlyBestGsPerMotif) {
    met_byDb <- split(motifEnrichment_selfMotifs, motifEnrichment_selfMotifs$motifDb)
    for (db in names(met_byDb)) {
      met <- met_byDb[[db]]
      met <- split(met, factor(met$highlightedTFs))
      met <- lapply(met, function(x) {
        rbindlist(lapply(split(x, x$motif), function(y) y[which.max(y$NES), 
        ]))
      })
      met_byDb[[db]] <- rbindlist(met)
    }
    motifEnrichment_selfMotifs <- rbindlist(met_byDb)
    rm(met_byDb)
    rm(met)
  }
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Pruning targets")
  if (getSettings(scenicOptions, "verbose")) 
    message(msg)
  dbNames <- getDatabases(scenicOptions)
  motifEnrichment_selfMotifs_wGenes <- lapply(names(dbNames), 
                                              function(motifDbName) {
                                                ranking <- importRankings(dbNames[motifDbName], 
                                                                          columns = allGenes)
                                                addSignificantGenes(resultsTable = motifEnrichment_selfMotifs[motifEnrichment_selfMotifs$motifDb == 
                                                                                                                motifDbName, ], geneSets = tfModules, rankings = ranking, 
                                                                    plotCurve = FALSE, maxRank = 5000, method = signifGenesMethod, 
                                                                    nMean = 100, nCores = 1)
                                                              
                                              })
  suppressPackageStartupMessages(library(data.table))
  motifEnrichment_selfMotifs_wGenes <- rbindlist(motifEnrichment_selfMotifs_wGenes)
  saveRDS(motifEnrichment_selfMotifs_wGenes, file = getIntName(scenicOptions, 
                                                               "motifEnrichment_selfMotifs_wGenes"))
  if (getSettings(scenicOptions, "verbose")) {
    message(format(Sys.time(), "%H:%M"), "\tNumber of motifs that support the regulons: ", 
            nrow(motifEnrichment_selfMotifs_wGenes))
    motifEnrichment_selfMotifs_wGenes[order(motifEnrichment_selfMotifs_wGenes$NES, 
                                            decreasing = TRUE), ][1:5, (1:ncol(motifEnrichment_selfMotifs_wGenes) - 
                                                                          1), with = F]
  }
  if (!file.exists("output")) 
    dir.create("output")
  write.table(motifEnrichment_selfMotifs_wGenes, file = getOutName(scenicOptions, 
                                                                   "s2_motifEnrichment"), sep = "\t", quote = FALSE, row.names = FALSE)
  if ("DT" %in% installed.packages() && nrow(motifEnrichment_selfMotifs_wGenes) > 
      0) {
    nvm <- tryCatch({
      colsToShow <- c("motifDb", "logo", "NES", "geneSet", 
                      "TF_highConf", "TF_lowConf")
      motifEnrichment_2html <- viewMotifs(motifEnrichment_selfMotifs_wGenes, 
                                          colsToShow = colsToShow, options = list(pageLength = 100))
      fileName <- getOutName(scenicOptions, "s2_motifEnrichmentHtml")
      dirName <- dirname(fileName)
      fileName <- basename(fileName)
      suppressWarnings(DT::saveWidget(motifEnrichment_2html, 
                                      fileName))
      file.rename(fileName, file.path(dirName, fileName))
      if (getSettings(scenicOptions, "verbose")) 
        message("\tPreview of motif enrichment saved as: ", 
                file.path(dirName, fileName))
    }, error = function(e) print(e$message))
  }
  motifEnrichment.asIncidList <- apply(motifEnrichment_selfMotifs_wGenes, 
                                       1, function(oneMotifRow) {
                                         genes <- strsplit(oneMotifRow["enrichedGenes"], 
                                                           ";")[[1]]
                                         oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors = FALSE)
                                         data.frame(oneMotifRow[rep(1, length(genes)), c("NES", 
                                                                                         "motif", "highlightedTFs", "TFinDB", "geneSet", 
                                                                                         "motifDb")], genes, stringsAsFactors = FALSE)
                                       })
  motifEnrichment.asIncidList <- rbindlist(motifEnrichment.asIncidList)
  colnames(motifEnrichment.asIncidList)[which(colnames(motifEnrichment.asIncidList) == 
                                                "highlightedTFs")] <- "TF"
  colnames(motifEnrichment.asIncidList)[which(colnames(motifEnrichment.asIncidList) == 
                                                "TFinDB")] <- "annot"
  colnames(motifEnrichment.asIncidList)[which(colnames(motifEnrichment.asIncidList) == 
                                                "genes")] <- "gene"
  motifEnrichment.asIncidList <- data.frame(motifEnrichment.asIncidList, 
                                            stringsAsFactors = FALSE)
  regulonTargetsInfo <- lapply(split(motifEnrichment.asIncidList, 
                                     motifEnrichment.asIncidList$TF), function(tfTargets) {
                                       tfTable <- as.data.frame(do.call(rbind, lapply(split(tfTargets, 
                                                                                            tfTargets$gene), function(enrOneGene) {
                                                                                              highConfAnnot <- "**" %in% enrOneGene$annot
                                                                                              enrOneGeneByAnnot <- enrOneGene
                                                                                              if (highConfAnnot) 
                                                                                                enrOneGeneByAnnot <- enrOneGeneByAnnot[which(enrOneGene$annot == 
                                                                                                                                               "**"), ]
                                                                                              bestMotif <- which.max(enrOneGeneByAnnot$NES)
                                                                                              tf <- unique(enrOneGene$TF)
                                                                                              cbind(TF = tf, gene = unique(enrOneGene$gene), highConfAnnot = highConfAnnot, 
                                                                                                    nMotifs = nrow(enrOneGene), bestMotif = as.character(enrOneGeneByAnnot[bestMotif, 
                                                                                                                                                                           "motif"]), NES = as.numeric(enrOneGeneByAnnot[bestMotif, 
                                                                                                                                                                                                                         "NES"]), motifDb = as.character(enrOneGeneByAnnot[bestMotif, 
                                                                                                                                                                                                                                                                           "motifDb"]), coexModule = gsub(paste0(tf, 
                                                                                                                                                                                                                                                                                                                 "_"), "", as.character(enrOneGeneByAnnot[bestMotif, 
                                                                                                                                                                                                                                                                                                                                                          "geneSet"]), fixed = TRUE))
                                                                                            })), stringsAsFactors = FALSE)
                                       tfTable[order(tfTable$NES, decreasing = TRUE), ]
                                     })
  rm(motifEnrichment.asIncidList)
  regulonTargetsInfo <- rbindlist(regulonTargetsInfo)
  corrMat <- loadInt(scenicOptions, "corrMat", ifNotExists = "null")
  if (!is.null(corrMat)) {
    regulonTargetsInfo$spearCor <- NA_real_
    for (tf in unique(regulonTargetsInfo$TF)) {
      regulonTargetsInfo[which(regulonTargetsInfo$TF == 
                                 tf), "spearCor"] <- corrMat[tf, unlist(regulonTargetsInfo[which(regulonTargetsInfo$TF == 
                                                                                                   tf), "gene"])]
    }
  }
  else warning("It was not possible to add the correlation to the regulonTargetsInfo table.")
  linkList <- loadInt(scenicOptions, "genie3ll", ifNotExists = "null")
  if (!is.null(linkList) & ("weight" %in% colnames(linkList))) {
    if (is.data.table(linkList)) 
      linkList <- as.data.frame(linkList)
    uniquePairs <- nrow(unique(linkList[, c("TF", "Target")]))
    if (uniquePairs == nrow(linkList)) {
      linkList <- linkList[which(linkList$weight >= getSettings(scenicOptions, 
                                                                "modules/weightThreshold")), ]
      rownames(linkList) <- paste(linkList$TF, linkList$Target, 
                                  sep = "__")
      regulonTargetsInfo <- cbind(regulonTargetsInfo, 
                                  CoexWeight = linkList[paste(regulonTargetsInfo$TF, 
                                                              regulonTargetsInfo$gene, sep = "__"), "weight"])
    }
    else {
      warning("There are duplicated regulator-target (gene id/name) pairs in the co-expression link list.", 
              "\nThe co-expression weight was not added to the regulonTargetsInfo table.")
    }
  }
  else warning("It was not possible to add the weight to the regulonTargetsInfo table.")
  saveRDS(regulonTargetsInfo, file = getIntName(scenicOptions, 
                                                "regulonTargetsInfo"))
  write.table(regulonTargetsInfo, file = getOutName(scenicOptions, 
                                                    "s2_regulonTargetsInfo"), sep = "\t", col.names = TRUE, 
              row.names = FALSE, quote = FALSE)
  rm(linkList)
  regulonTargetsInfo_splitByAnnot <- split(regulonTargetsInfo, 
                                           regulonTargetsInfo$highConfAnnot)
  regulons <- NULL
  if (!is.null(regulonTargetsInfo_splitByAnnot[["TRUE"]])) {
    regulons <- lapply(split(regulonTargetsInfo_splitByAnnot[["TRUE"]], 
                             regulonTargetsInfo_splitByAnnot[["TRUE"]][, "TF"]), 
                       function(x) sort(as.character(unlist(x[, "gene"]))))
  }
  regulons_extended <- NULL
  if (!is.null(regulonTargetsInfo_splitByAnnot[["FALSE"]])) {
    regulons_extended <- lapply(split(regulonTargetsInfo_splitByAnnot[["FALSE"]], 
                                      regulonTargetsInfo_splitByAnnot[["FALSE"]][, "TF"]), 
                                function(x) unname(unlist(x[, "gene"])))
    regulons_extended <- setNames(lapply(names(regulons_extended), 
                                         function(tf) sort(unique(c(regulons[[tf]], unlist(regulons_extended[[tf]]))))), 
                                  names(regulons_extended))
    names(regulons_extended) <- paste(names(regulons_extended), 
                                      "_extended", sep = "")
  }
  regulons <- c(regulons, regulons_extended)
  saveRDS(regulons, file = getIntName(scenicOptions, "regulons"))
  incidList <- reshape2::melt(regulons)
  incidMat <- table(incidList[, 2], incidList[, 1])
  saveRDS(incidMat, file = getIntName(scenicOptions, "regulons_incidMat"))
  rm(incidMat)
  if (getSettings(scenicOptions, "verbose")) {
    length(regulons)
    summary(lengths(regulons))
  }
  scenicOptions@status$current <- 2
  invisible(scenicOptions)
}
