ploteqtl<-function (GWAS.df, eQTL.df, Genes.df, LD.df = TRUE, gene, trait, 
                    sigpvalue_GWAS = 5e-08, sigpvalue_eQTL = 0.05, tissue = "all", 
                    range = 200, NESeQTLRange = c(NA, NA), congruence = FALSE, 
                    R2min = 0.2, LDmin = 10, leadSNP = TRUE, LDcolor = "color", 
                    ylima = NA, ylimd = NA, xlimd = NA, genometrackheight = 2, 
                    gbuild = "hg19", res = 300, wi = "wi", CollapseMethod = "min", 
                    getplot = TRUE, saveplot = TRUE, GeneList = FALSE, TissueList = FALSE) 
{
  eQTpLot.gene.list <- function(gene) {
    rangebp <- range * 1000
    startpos <- (min(Genes.df[which(Genes.df$Gene == gene & 
                                      Genes.df$Build == gbuild), ] %>% dplyr::select(Start)) - 
                   rangebp)
    stoppos <- (max(Genes.df[which(Genes.df$Gene == gene & 
                                     Genes.df$Build == gbuild), ] %>% dplyr::select(Stop)) + 
                  rangebp)
    chromosome <- Genes.df$CHR[Genes.df$Gene == gene & Genes.df$Build == 
                                 gbuild]
    gwas.data <- GWAS.df[which(GWAS.df$CHR == chromosome & 
                                 GWAS.df$PHE == trait & GWAS.df$BP >= startpos & 
                                 GWAS.df$BP <= stoppos & !(is.na(GWAS.df$P)) & !(is.na(GWAS.df$BETA))), 
    ]
    if (dim(gwas.data)[1] == 0) 
      stop("Sorry, GWAS.df does not contain data for any SNPs in the ", 
           paste(range), "kb flanking the gene ", paste(gene), 
           " for the trait ", paste(trait))
    if (dim(gwas.data[which(gwas.data$P <= sigpvalue_GWAS), 
    ])[1] == 0) {
      print(paste(sep = "", "WARNING: GWAS.df does not contain any SNPs with p-value < sigpvalue_GWAS within the ", 
                  range, "kb flanking the gene ", gene, " for the trait ", 
                  trait))
    }
    if (GeneList == TRUE & TissueList == FALSE) {
      if (length(tissue) >= 2 & ("all" %in% tissue) == 
          FALSE) {
        eQTL.df <- eQTL.df[which(eQTL.df$Tissue %in% 
                                   tissue), ]
      }
      eqtl.data <- eQTL.df[which(eQTL.df$Gene.Symbol == 
                                   gene & eQTL.df$P.Value <= sigpvalue_eQTL & !(is.na(eQTL.df$NES)) & 
                                   !(is.na(eQTL.df$P.Value))), ]
      if (dim(eqtl.data)[1] == 0) 
        stop("Sorry, eQTL.df does not have any data for the gene ", 
             paste(gene), " meeting your sigpvalue_eQTL threshold")
      if (any(tissue == "all") == TRUE | (length(tissue) >= 
                                          2) == TRUE) {
        PanTissue <- TRUE
      }
      if (("all" %in% tissue) == FALSE & (length(tissue) >= 
                                          2) == FALSE) {
        PanTissue <- FALSE
      }
      if (PanTissue == TRUE & (CollapseMethod == "mean" | 
                               CollapseMethod == "median")) {
        Mean.Ps <- eqtl.data[which(!(is.na(eqtl.data$SNP.Id))), 
        ] %>% dplyr::group_by(SNP.Id) %>% dplyr::summarize(Mean.P = mean(P.Value, 
                                                                         na.rm = TRUE))
        Mean.NESs <- eqtl.data[which(!(is.na(eqtl.data$SNP.Id))), 
        ] %>% dplyr::group_by(SNP.Id) %>% dplyr::summarize(Mean.NES = mean(NES, 
                                                                           na.rm = TRUE))
        Median.Ps <- eqtl.data[which(!(is.na(eqtl.data$SNP.Id))), 
        ] %>% dplyr::group_by(SNP.Id) %>% dplyr::summarize(Median.P = median(P.Value, 
                                                                             na.rm = TRUE))
        Median.NESs <- eqtl.data[which(!(is.na(eqtl.data$SNP.Id))), 
        ] %>% dplyr::group_by(SNP.Id) %>% dplyr::summarize(Median.NES = median(NES, 
                                                                               na.rm = TRUE))
      }
      if (PanTissue == TRUE & CollapseMethod == "meta") {
        if (notrunmeta == TRUE) {
          eqtl.data$N <- 100
        }
        eqtl.data$sign <- ifelse(eqtl.data$NES < 0, 
                                 -1, 1)
        eqtl.data$Z <- (qnorm(eqtl.data$P.Value/2, lower.tail = FALSE)) * 
          (eqtl.data$sign)
        eqtl.data$W <- sqrt(eqtl.data$N)
        eqtl.data$ZW <- (eqtl.data$Z) * (eqtl.data$W)
        eqtl.data$W2 <- (eqtl.data$W) * (eqtl.data$W)
        eqtl.data.sumZW <- eqtl.data %>% dplyr::group_by(SNP.Id) %>% 
          dplyr::summarise(SumZW = sum(ZW, na.rm = TRUE))
        eqtl.data.sumW2 <- eqtl.data %>% dplyr::group_by(SNP.Id) %>% 
          dplyr::summarise(SumW2 = sum(W2, na.rm = TRUE))
        eqtl.data.sumW2$sqrtSumW2 <- sqrt(eqtl.data.sumW2$SumW2)
        eqtl.data.sumW2$SumW2 <- NULL
        eqtl.data.sum.final <- dplyr::left_join(eqtl.data.sumW2, 
                                                eqtl.data.sumZW, by = "SNP.Id")
        eqtl.data.sum.final$Z <- (eqtl.data.sum.final$SumZW)/(eqtl.data.sum.final$sqrtSumW2)
        eqtl.data.sum.final$P <- 2 * (pnorm(-abs(eqtl.data.sum.final$Z)))
        eqtl.data$W <- NULL
        eqtl.data$Z <- NULL
        eqtl.data$ZW <- NULL
        eqtl.data$W2 <- NULL
        eqtl.data.sum.final$sqrtSumW2 <- NULL
        eqtl.data.sum.final$SumZW <- NULL
        eqtl.data.sum.final$Z <- NULL
      }
      if (PanTissue == TRUE) {
        eqtl.data <- eqtl.data[which(!(is.na(eqtl.data$SNP.Id))), 
        ] %>% dplyr::group_by(SNP.Id) %>% dplyr::slice(which.min(P.Value))
        eqtl.data$Tissue <- "PanTissue"
      }
      else {
        eqtl.data <- eqtl.data[which(eqtl.data$Tissue == 
                                       tissue), ]
      }
      if (CollapseMethod == "mean" & PanTissue == TRUE) {
        eqtl.data <- dplyr::left_join(eqtl.data, Mean.Ps, 
                                      by = "SNP.Id")
        eqtl.data <- dplyr::left_join(eqtl.data, Mean.NESs, 
                                      by = "SNP.Id")
        eqtl.data$P.Value <- eqtl.data$Mean.P
        eqtl.data$NES <- eqtl.data$Mean.NES
        eqtl.data$Mean.P <- NULL
        eqtl.data$Mean.NES <- NULL
      }
      if (CollapseMethod == "median" & PanTissue == TRUE) {
        eqtl.data <- dplyr::left_join(eqtl.data, Median.Ps, 
                                      by = "SNP.Id")
        eqtl.data <- dplyr::left_join(eqtl.data, Median.NESs, 
                                      by = "SNP.Id")
        eqtl.data$P.Value <- eqtl.data$Median.P
        eqtl.data$NES <- eqtl.data$Median.NES
        eqtl.data$Median.P <- NULL
        eqtl.data$Median.NES <- NULL
      }
      if (CollapseMethod == "meta" & PanTissue == TRUE) {
        eqtl.data <- dplyr::left_join(eqtl.data, eqtl.data.sum.final, 
                                      by = "SNP.Id")
        eqtl.data$P.Value <- eqtl.data$P
        eqtl.data$NES <- (eqtl.data$sign)/3
        eqtl.data$sign <- NULL
        eqtl.data$N <- NULL
        eqtl.data$P <- NULL
      }
    }
    eqtl.data$P.Value <- ifelse(eqtl.data$P.Value <= 1e-300, 
                                1e-300, eqtl.data$P.Value)
    if (dim(eqtl.data)[1] == 0) 
      stop("Sorry, there are no eQTLs for the tissue", 
           paste(tissue), " with a p-value < sigeQTL")
    eqtl.data <- dplyr::ungroup(eqtl.data)
    gwas.data$SNP <- as.factor(gwas.data$SNP)
    eqtl.data$SNP.Id <- as.factor(eqtl.data$SNP.Id)
    combinedSNPS <- sort(union(levels(gwas.data$SNP), levels(eqtl.data$SNP.Id)))
    Combined.eQTL.GWAS.Data <- dplyr::left_join(dplyr::mutate(gwas.data, 
                                                              SNP = factor(SNP, levels = combinedSNPS)), dplyr::mutate(eqtl.data, 
                                                                                                                       SNP.Id = factor(SNP.Id, levels = combinedSNPS)) %>% 
                                                  dplyr::rename(SNP = SNP.Id), by = "SNP")
    if (dim(Combined.eQTL.GWAS.Data)[1] == 0) {
      stop("Sorry, for the gene ", paste(gene), " and the trait ", 
           paste(trait), " there is no overlap between the SNPs in your GWAS.df and eQTL.df")
    }
    Combined.eQTL.GWAS.Data$DirectionOfEffect_GWAS <- ifelse(Combined.eQTL.GWAS.Data$BETA < 
                                                               0, "Negative", ifelse(Combined.eQTL.GWAS.Data$BETA > 
                                                                                       0, "Positive", NA))
    Combined.eQTL.GWAS.Data$DirectionOfEffect_eQTL <- ifelse(Combined.eQTL.GWAS.Data$NES < 
                                                               0, "DOWN", ifelse(Combined.eQTL.GWAS.Data$NES > 
                                                                                   0, "UP", NA))
    Combined.eQTL.GWAS.Data$Congruence <- (Combined.eQTL.GWAS.Data$BETA * 
                                             Combined.eQTL.GWAS.Data$NES)
    Combined.eQTL.GWAS.Data$Congruence <- ifelse(Combined.eQTL.GWAS.Data$Congruence < 
                                                   0, "Incongruent", ifelse(Combined.eQTL.GWAS.Data$Congruence > 
                                                                              0, "Congruent", NA))
    if (congruence == FALSE) {
      Combined.eQTL.GWAS.Data$Congruence <- ifelse(Combined.eQTL.GWAS.Data$Congruence == 
                                                     "Incongruent", "Congruent", ifelse(Combined.eQTL.GWAS.Data$Congruence == 
                                                                                          "Congruent", "Congruent", NA))
    }
    Combined.eQTL.GWAS.Data$NeglogeQTLpValue <- -(log10(Combined.eQTL.GWAS.Data$P.Value))
    Combined.eQTL.GWAS.Data$Neglog10pvalue_GWAS <- -(log10(Combined.eQTL.GWAS.Data$P))
    Combined.eQTL.GWAS.Data <- Combined.eQTL.GWAS.Data[which(!(is.na(Combined.eQTL.GWAS.Data$P))), 
    ]
    Combined.eQTL.GWAS.Data$significance <- ifelse(Combined.eQTL.GWAS.Data$P >= 
                                                     sigpvalue_GWAS, "Non-significant", "Significant")
    Combined.eQTL.GWAS.Data$Congruence[is.na(Combined.eQTL.GWAS.Data$Congruence)] <- "Non-eQTL"
    Combined.eQTL.GWAS.Data <- Combined.eQTL.GWAS.Data %>% 
      dplyr::mutate(Congruence = factor(Congruence, levels = c("Non-eQTL", 
                                                               "Congruent", "Incongruent"), ordered = TRUE))
    if (dim(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                          Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
    ])[1] == 0) {
      Congruentdata <- FALSE
    }
    else {
      Congruentdata <- TRUE
    }
    if (dim(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                          Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
    ])[1] == 0) {
      Incongruentdata <- FALSE
    }
    else {
      Incongruentdata <- TRUE
    }
    if (Congruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                   Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
    ]) >= 3) {
      pearson.congruent <- suppressWarnings(cor.test(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                                     Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
      ]$NeglogeQTLpValue, Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                          Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
      ]$Neglog10pvalue_GWAS, method = "pearson"))
      pearson.congruent.logic <- TRUE
    }
    else {
      pearson.congruent.logic <- FALSE
      pearson.congruent <- list(NA, NA)
      pearson.congruent$estimate <- 0
      pearson.congruent$p.value <- 1
      print("Not enough data to complete Pearson correlation for Congruent eQTLs")
    }
    if (Incongruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                     Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
    ]) >= 3) {
      pearson.incongruent <- suppressWarnings(cor.test(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                                       Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
      ]$NeglogeQTLpValue, Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                          Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
      ]$Neglog10pvalue_GWAS, method = "pearson"))
      pearson.incongruent.logic <- TRUE
    }
    else {
      if (congruence == TRUE) {
        pearson.incongruent.logic <- FALSE
        pearson.incongruent <- list(NA, NA)
        pearson.incongruent$estimate <- 0
        pearson.incongruent$p.value <- 1
        print("Not enough data to complete Pearson correlation for Incongruent eQTLs")
      }
    }
    if (congruence == FALSE) {
      if (pearson.congruent.logic == TRUE) {
        out1 <- if (pearson.congruent.logic == TRUE) {
          paste(sep = "", "eQTL analysis for gene ", 
                gene, ": Pearson correlation: ", round(pearson.congruent$estimate, 
                                                       3), ", p-value: ", formatC(pearson.congruent$p.value, 
                                                                                  format = "e", digits = 2))
        }
        pvalue1 <- pearson.congruent$p.value
        out <- data.frame(output = out1, pval = pvalue1)
        return(out)
      }
    }
    if (congruence == TRUE) {
      out1 <- if (pearson.congruent.logic == TRUE) {
        paste(sep = "", "Congruent eQTL analysis for gene ", 
              gene, ": Pearson correlation: ", round(pearson.congruent$estimate, 
                                                     3), ", p-value: ", formatC(pearson.congruent$p.value, 
                                                                                format = "e", digits = 2))
      }
      pvalue1 <- pearson.congruent$p.value
      out2 <- if (pearson.incongruent.logic == TRUE) {
        paste(sep = "", "Incongruent eQTL analysis for gene ", 
              gene, ": Pearson correlation: ", round(pearson.incongruent$estimate, 
                                                     3), ", p-value: ", formatC(pearson.incongruent$p.value, 
                                                                                format = "e", digits = 2))
      }
      pvalue2 <- pearson.incongruent$p.value
      out <- data.frame(output = c(out1, out2), pval = c(pvalue1, 
                                                         pvalue2))
      return(out)
    }
  }
  eQTpLot.tissue.list <- function(tissue) {
    rangebp <- range * 1000
    startpos <- (min(Genes.df[which(Genes.df$Gene == gene & 
                                      Genes.df$Build == gbuild), ] %>% dplyr::select(Start)) - 
                   rangebp)
    stoppos <- (max(Genes.df[which(Genes.df$Gene == gene & 
                                     Genes.df$Build == gbuild), ] %>% dplyr::select(Stop)) + 
                  rangebp)
    chromosome <- Genes.df$CHR[Genes.df$Gene == gene & Genes.df$Build == 
                                 gbuild]
    gwas.data <- GWAS.df[which(GWAS.df$CHR == chromosome & 
                                 GWAS.df$PHE == trait & GWAS.df$BP >= startpos & 
                                 GWAS.df$BP <= stoppos & !(is.na(GWAS.df$P)) & !(is.na(GWAS.df$BETA))), 
    ]
    if (dim(gwas.data)[1] == 0) 
      stop("Sorry, GWAS.df does not contain data for any SNPs in the ", 
           paste(range), "kb flanking the gene ", paste(gene), 
           " for the trait ", paste(trait))
    if (dim(gwas.data[which(gwas.data$P <= sigpvalue_GWAS), 
    ])[1] == 0) {
      print(paste(sep = "", "WARNING: GWAS.df does not contain any SNPs with p-value < sigpvalue_GWAS within the ", 
                  range, "kb flanking the gene ", gene, " for the trait ", 
                  trait))
    }
    eqtl.data <- eQTL.df[which(eQTL.df$Gene.Symbol == gene & 
                                 eQTL.df$P.Value <= sigpvalue_eQTL & !(is.na(eQTL.df$NES)) & 
                                 !(is.na(eQTL.df$P.Value))), ]
    if (dim(eqtl.data)[1] == 0) 
      stop("Sorry, eQTL.df does not have any data for the gene ", 
           paste(gene), " meeting your sigpvalue_eQTL threshold")
    eqtl.data <- eqtl.data[which(eqtl.data$Tissue == tissue), 
    ]
    eqtl.data$P.Value <- ifelse(eqtl.data$P.Value <= 1e-300, 
                                1e-300, eqtl.data$P.Value)
    if (dim(eqtl.data)[1] == 0) 
      stop("Sorry, there are no eQTLs for the tissue", 
           paste(tissue), " with a p-value < sigeQTL")
    eqtl.data <- dplyr::ungroup(eqtl.data)
    gwas.data$SNP <- as.factor(gwas.data$SNP)
    eqtl.data$SNP.Id <- as.factor(eqtl.data$SNP.Id)
    combinedSNPS <- sort(union(levels(gwas.data$SNP), levels(eqtl.data$SNP.Id)))
    Combined.eQTL.GWAS.Data <- dplyr::left_join(dplyr::mutate(gwas.data, 
                                                              SNP = factor(SNP, levels = combinedSNPS)), dplyr::mutate(eqtl.data, 
                                                                                                                       SNP.Id = factor(SNP.Id, levels = combinedSNPS)) %>% 
                                                  dplyr::rename(SNP = SNP.Id), by = "SNP")
    if (dim(Combined.eQTL.GWAS.Data)[1] == 0) {
      stop("Sorry, for the gene ", paste(gene), " and the trait ", 
           paste(trait), " there is no overlap between the SNPs in your GWAS.df and eQTL.df")
    }
    Combined.eQTL.GWAS.Data$DirectionOfEffect_GWAS <- ifelse(Combined.eQTL.GWAS.Data$BETA < 
                                                               0, "Negative", ifelse(Combined.eQTL.GWAS.Data$BETA > 
                                                                                       0, "Positive", NA))
    Combined.eQTL.GWAS.Data$DirectionOfEffect_eQTL <- ifelse(Combined.eQTL.GWAS.Data$NES < 
                                                               0, "DOWN", ifelse(Combined.eQTL.GWAS.Data$NES > 
                                                                                   0, "UP", NA))
    Combined.eQTL.GWAS.Data$Congruence <- (Combined.eQTL.GWAS.Data$BETA * 
                                             Combined.eQTL.GWAS.Data$NES)
    Combined.eQTL.GWAS.Data$Congruence <- ifelse(Combined.eQTL.GWAS.Data$Congruence < 
                                                   0, "Incongruent", ifelse(Combined.eQTL.GWAS.Data$Congruence > 
                                                                              0, "Congruent", NA))
    if (congruence == FALSE) {
      Combined.eQTL.GWAS.Data$Congruence <- ifelse(Combined.eQTL.GWAS.Data$Congruence == 
                                                     "Incongruent", "Congruent", ifelse(Combined.eQTL.GWAS.Data$Congruence == 
                                                                                          "Congruent", "Congruent", NA))
    }
    Combined.eQTL.GWAS.Data$NeglogeQTLpValue <- -(log10(Combined.eQTL.GWAS.Data$P.Value))
    Combined.eQTL.GWAS.Data$Neglog10pvalue_GWAS <- -(log10(Combined.eQTL.GWAS.Data$P))
    Combined.eQTL.GWAS.Data <- Combined.eQTL.GWAS.Data[which(!(is.na(Combined.eQTL.GWAS.Data$P))), 
    ]
    Combined.eQTL.GWAS.Data$significance <- ifelse(Combined.eQTL.GWAS.Data$P >= 
                                                     sigpvalue_GWAS, "Non-significant", "Significant")
    Combined.eQTL.GWAS.Data$Congruence[is.na(Combined.eQTL.GWAS.Data$Congruence)] <- "Non-eQTL"
    Combined.eQTL.GWAS.Data <- Combined.eQTL.GWAS.Data %>% 
      dplyr::mutate(Congruence = factor(Congruence, levels = c("Non-eQTL", 
                                                               "Congruent", "Incongruent"), ordered = TRUE))
    if (dim(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                          Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
    ])[1] == 0) {
      Congruentdata <- FALSE
    }
    else {
      Congruentdata <- TRUE
    }
    if (dim(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                          Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
    ])[1] == 0) {
      Incongruentdata <- FALSE
    }
    else {
      Incongruentdata <- TRUE
    }
    if (Congruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                   Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
    ]) >= 3) {
      pearson.congruent <- suppressWarnings(cor.test(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                                     Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
      ]$NeglogeQTLpValue, Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                          Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
      ]$Neglog10pvalue_GWAS, method = "pearson"))
      pearson.congruent.logic <- TRUE
    }
    else {
      pearson.congruent.logic <- FALSE
      pearson.congruent <- list(NA, NA)
      pearson.congruent$estimate <- 0
      pearson.congruent$p.value <- 1
      print(paste("Not enough data to complete Pearson correlation for Congruent eQTLs for tissue", 
                  tissue))
    }
    if (Incongruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                     Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
    ]) >= 3) {
      pearson.incongruent <- suppressWarnings(cor.test(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                                       Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
      ]$NeglogeQTLpValue, Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                          Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
      ]$Neglog10pvalue_GWAS, method = "pearson"))
      pearson.incongruent.logic <- TRUE
    }
    else {
      if (congruence == TRUE) {
        pearson.incongruent.logic <- FALSE
        pearson.incongruent <- list(NA, NA)
        pearson.incongruent$estimate <- 0
        pearson.incongruent$p.value <- 1
        print(paste("Not enough data to complete Pearson correlation for Incongruent eQTLs for tissue", 
                    tissue))
      }
    }
    if (congruence == FALSE) {
      out1 <- if (pearson.congruent.logic == TRUE) {
        paste(sep = "", "eQTL analysis for tissue ", 
              tissue, ": Pearson correlation: ", round(pearson.congruent$estimate, 
                                                       3), ", p-value: ", formatC(pearson.congruent$p.value, 
                                                                                  format = "e", digits = 2))
      }
      pvalue1 <- pearson.congruent$p.value
      out <- data.frame(output = out1, pval = pvalue1)
      return(out)
    }
    if (congruence == TRUE) {
      out1 <- if (pearson.congruent.logic == TRUE) {
        paste(sep = "", "Congruent eQTL analysis for tissue ", 
              tissue, ": Pearson correlation: ", round(pearson.congruent$estimate, 
                                                       3), ", p-value: ", formatC(pearson.congruent$p.value, 
                                                                                  format = "e", digits = 2))
      }
      pvalue1 <- pearson.congruent$p.value
      out2 <- if (pearson.incongruent.logic == TRUE) {
        paste(sep = "", "Incongruent eQTL analysis for tissue ", 
              tissue, ": Pearson correlation: ", round(pearson.incongruent$estimate, 
                                                       3), ", p-value: ", formatC(pearson.incongruent$p.value, 
                                                                                  format = "e", digits = 2))
      }
      pvalue2 <- pearson.incongruent$p.value
      out <- data.frame(output = c(out1, out2), pval = c(pvalue1, 
                                                         pvalue2))
      return(out)
    }
  }
  print("Checking input data...")
  if (GeneList == TRUE & TissueList == TRUE) {
    stop("You can only perform either a GeneList analysis or a TissueList analysis, you cannot perform both at the same time. Please set one of these parameters to false.")
  }
  if (GeneList != TRUE & length(gene) >= 2 & TissueList != 
      TRUE) {
    stop("Please select only a single gene to complete eQTpLot visualization. Current selection is: ", 
         paste(gene, " ", sep = ","))
  }
  if (missing(Genes.df)) {
    Genes.df <- eQTpLot:::genes
  }
  if (all(c("CHR", "BP", "SNP", "BETA", "P") %in% colnames(GWAS.df)) == 
      FALSE) {
    stop("The data supplied to GWAS.df must contain columns 'CHR', 'BP', 'SNP', 'BETA', and 'P'")
  }
  if ("PHE" %in% colnames(GWAS.df)) {
    print(paste(sep = "", "PHE column found in GWAS.df. Analyzing data for phenotype ", 
                trait))
  }
  else {
    print(paste(sep = "", "PHE column not found in GWAS.df. Assuming all data in GWAS.df is for phenotype ", 
                trait))
    GWAS.df$PHE <- trait
  }
  if ((trait %in% GWAS.df$PHE) == FALSE) {
    stop("Sorry, the  phenotype ", paste(trait), " does not exist in the PHE column of the GWAS.df dataframe. Phenotypes included in the data supplied to GWAS.df are:\n", 
         paste("'", as.character(unique(GWAS.df$PHE)), "'", 
               collapse = ", ", sep = ""))
  }
  if (all(c("SNP.Id", "Gene.Symbol", "P.Value", "NES", "Tissue") %in% 
          colnames(eQTL.df)) == FALSE) {
    stop("The data supplied to eQTL.df must contain columns 'SNP.Id', 'Gene.Symbol', 'P.Value', 'NES', and 'Tissue'")
  }
  if (all(c("Gene", "CHR", "Start", "Stop", "Build") %in% 
          colnames(Genes.df)) == FALSE) {
    stop("The data supplied to Genes.df must contain columns 'Gene', 'CHR', 'Start', 'Stop', 'Build'")
  }
  if (is.numeric(eQTL.df$P.Value) == FALSE | is.numeric(eQTL.df$NES) == 
      FALSE) {
    stop("Sorry, the eQTL.df dataframe must contain only numeric data for P.Value and NES")
  }
  if (all("N" %in% colnames(eQTL.df)) == TRUE) {
    if (is.numeric(eQTL.df$N) == FALSE & is.integer(eQTL.df$N) == 
        FALSE) {
      stop("Sorry, the  column N in eQTL.df must contain only numeric values")
    }
  }
  if (is.numeric(GWAS.df$P) == FALSE | is.numeric(GWAS.df$BETA) == 
      FALSE | (is.integer(GWAS.df$BP) == FALSE & is.numeric(GWAS.df$BP) == 
               FALSE) | (is.integer(GWAS.df$CHR) == FALSE & is.numeric(GWAS.df$CHR) == 
                         FALSE)) {
    stop("Sorry, the GWAS.df dataframe must contain only numeric data for CHR, BP, P, and BETA (Note: chromosomes must be coded numerically)")
  }
  if ((is.integer(Genes.df$CHR) == FALSE & is.numeric(Genes.df$CHR) == 
       FALSE) | (is.integer(Genes.df$Start) == FALSE & is.numeric(Genes.df$Start) == 
                 FALSE) | (is.integer(Genes.df$Stop) == FALSE & is.numeric(Genes.df$Stop) == 
                           FALSE)) {
    stop("Sorry, the Genes.df dataframe must contain only integer values for CHR, Start, and Stop (Note: chromosomes must be coded nuemrically)")
  }
  if (length(gene) == 1 & any(gene %in% (Genes.df[which(Genes.df$Build == 
                                                        gbuild), ]$Gene) == FALSE)) {
    stop("Sorry, there is no information for the gene ", 
         paste(gene), " in the Genes.df dataframe for the genomic build ", 
         paste(gbuild), "\nConsider supplying your own Genes.df input file with information on this gene.")
  }
  if (length(gene) >= 2 & any(gene %in% (Genes.df[which(Genes.df$Build == 
                                                        gbuild), ]$Gene) == FALSE)) {
    stop("Sorry, there is no information for at least one of the supplied genes, ", 
         paste(gene, " ", sep = ","), " in the Genes.df dataframe for the genomic build ", 
         paste(gbuild), "\nConsider supplying your own Genes.df input file with information on this gene.")
  }
  if (all(tissue == "all") == FALSE & (all(tissue %in% eQTL.df$Tissue)) == 
      FALSE) {
    stop("Sorry, at least one of the specified tissues, ", 
         paste(tissue, " ", sep = ","), " does not exist in the eQTL.df dataframe. Tissues included in the data supplied to eQTL.df are:\n", 
         paste("'", as.character(unique(eQTL.df$Tissue)), 
               "'", collapse = ", ", sep = ""))
  }
  if (LDcolor != "color" & LDcolor != "black") {
    stop("Sorry, the argument LDcolor must be set to either \"color\" or \"black\"")
  }
  if (isTRUE(LD.df) == FALSE) {
    if (all(c("BP_A", "SNP_A", "BP_B", "SNP_B", "R2") %in% 
            colnames(LD.df)) == FALSE) {
      stop("The data supplied to LD.df must contain columns 'BP_A', 'SNP_A', 'BP_B', 'SNP_B', and 'R2'")
    }
    if ((is.integer(LD.df$BP_A) == FALSE | is.integer(LD.df$BP_B) == 
         FALSE)) {
      stop("Sorry, the LD.df dataframe must contain only integer values for BP_A and BP_B")
    }
    if (is.numeric(LD.df$R2) == FALSE) {
      stop("Sorry, the LD.df dataframe must contain only numeric values for R2")
    }
    if (isTRUE(leadSNP) == FALSE & (leadSNP %in% LD.df$SNP_A) == 
        FALSE & (leadSNP %in% LD.df$SNP_B) == FALSE) {
      stop("Sorry, the specified leadSNP is not present in your LD.df")
    }
    if (isTRUE(leadSNP) == FALSE & (leadSNP %in% GWAS.df$SNP) == 
        FALSE) {
      stop("Sorry, the specified leadSNP is not present in your GWAS.df")
    }
    if (isTRUE(leadSNP) == FALSE & (leadSNP %in% eQTL.df$SNP.Id) == 
        FALSE) {
      stop("Sorry, the specified leadSNP is not present in your eQTL.df")
    }
  }
  if (GeneList == TRUE) {
    if (CollapseMethod == "meta" & "N" %in% colnames(eQTL.df) == 
        FALSE & interactive() == TRUE) {
      notrunmeta <- askYesNo(default = TRUE, msg = paste("To complete a PanTissue or MultiTissue analysis using a meta-analysis approach, eQTL.df must contain a column with header N listing the sample size used for each eQTL calculation\nYour eQTL.df does not have a coumn N.\nDo you want to proceed with eQTL meta-analysis assuming all eQTL sample sizes are equal? CAUTION: this may not yield accurate results."))
    }
    else {
      notrunmeta <- "NA"
    }
    if (notrunmeta == FALSE) {
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      print("Stopping analysis")
      stop()
    }
    else {
      if (notrunmeta == TRUE) {
        print("Proceeding with eQTL meta-analysis assuming all eQTL sample sizes are equal")
        print("CAUTION: this may not yield accurate results")
      }
    }
    if (length(gene) <= 1) {
      stop("Please provide a list of at least two genes, in the format c(\"GENE1\", \"GENE2\" ... ), to complete this function")
    }
    if (length(tissue) <= 1 & ("all" %in% tissue) == FALSE) {
      print(paste(sep = "", "eQTL analysis will be completed for tissue:"))
    }
    if (length(tissue) >= 2) {
      print(paste(sep = "", "MultiTissue eQTL analysis, collapsing by method ", 
                  CollapseMethod, " will be completed across tissues:"))
    }
    if (length(tissue) <= 1 & ("all" %in% tissue) == TRUE) {
      print(paste("PanTissue eQTL analysis, collapsing by method ", 
                  CollapseMethod, " will be completed across all tissues in eQTL.df"))
    }
    if ((length(tissue) <= 1 & ("all" %in% tissue) == FALSE) | 
        (length(tissue) >= 2)) {
      print(paste("'", as.character(unique(tissue)), "'", 
                  collapse = ", ", sep = ""))
    }
    print("For genes:")
    print(paste("'", as.character(unique(gene)), "'", collapse = ", ", 
                sep = ""))
    results <- lapply(gene, eQTpLot.gene.list) %>% dplyr::bind_rows()
    results <- results %>% dplyr::arrange(as.numeric(pval))
    print(results$output)
    return("Complete")
  }
  if (TissueList == TRUE) {
    if (length(gene) >= 2) {
      stop("Please provide only a single gene to perform a TissueList analysis")
    }
    if (length(tissue) <= 1 & tissue != "all") {
      stop("Please provide at least two tissues to perform a TissueList analysis, or set tissue = \"all\" to perform a TissueList analysis for the gene ", 
           paste(gene), " across all tissues in eQTL.df")
    }
    if (tissue == "all") {
      tissue <- as.character(unique(eQTL.df$Tissue))
    }
    print(paste(sep = "", "eQTpLot TissueList analysis will be completed for tissues:"))
    print(paste(as.character(unique(tissue)), "'", collapse = ", ", 
                sep = ""))
    print("For gene:")
    print(paste(gene))
    results <- lapply(tissue, eQTpLot.tissue.list) %>% dplyr::bind_rows()
    results <- results %>% dplyr::arrange(as.numeric(pval))
    print(results$output)
    return("Complete")
  }
  print("Compiling GWAS and eQTL data...")
  rangebp <- range * 1000
  startpos <- (min(Genes.df[which(Genes.df$Gene == gene & 
                                    Genes.df$Build == gbuild), ] %>% dplyr::select(Start)) - 
                 rangebp)
  stoppos <- (max(Genes.df[which(Genes.df$Gene == gene & Genes.df$Build == 
                                   gbuild), ] %>% dplyr::select(Stop)) + rangebp)
  chromosome <- Genes.df$CHR[Genes.df$Gene == gene & Genes.df$Build == 
                               gbuild]
  gwas.data <- GWAS.df[which(GWAS.df$CHR == chromosome & GWAS.df$PHE == 
                               trait & GWAS.df$BP >= startpos & GWAS.df$BP <= stoppos & 
                               !(is.na(GWAS.df$P)) & !(is.na(GWAS.df$BETA))), ]
  if (dim(gwas.data)[1] == 0) 
    stop("Sorry, GWAS.df does not contain data for any SNPs in the ", 
         paste(range), "kb flanking the gene ", paste(gene), 
         " for the trait ", paste(trait))
  if (dim(gwas.data[which(gwas.data$P <= sigpvalue_GWAS), 
  ])[1] == 0) {
    NoFisher <- TRUE
    print(paste(sep = "", "WARNING: GWAS.df does not contain any SNPs with p-value < sigpvalue_GWAS within the ", 
                range, "kb flanking the gene ", gene, " for the trait ", 
                trait, ". eQTL Enrcihment Plot statistics will not be calculated"))
  }
  else {
    NoFisher <- FALSE
  }
  if (length(tissue) >= 2 & ("all" %in% tissue) == FALSE) {
    eQTL.df <- eQTL.df[which(eQTL.df$Tissue %in% tissue), 
    ]
  }
  print(paste(sep = "", "eQTL analysis will be completed for tissues ", 
              paste("'", as.character(unique(tissue)), "'", collapse = ", ", 
                    sep = ""), " and for gene ", paste(gene)))
  eqtl.data <- eQTL.df[which(eQTL.df$Gene.Symbol == gene & 
                               eQTL.df$P.Value <= sigpvalue_eQTL & !(is.na(eQTL.df$NES)) & 
                               !(is.na(eQTL.df$P.Value))), ]
  if (dim(eqtl.data)[1] == 0) 
    stop("Sorry, eQTL.df does not have any data for the gene ", 
         paste(gene), " meeting your sigpvalue_eQTL threshold")
  if (any(tissue == "all") == TRUE | (length(tissue) >= 2) == 
      TRUE) {
    PanTissue <- TRUE
  }
  if (any(tissue == "all") == FALSE & (length(tissue) >= 2) == 
      FALSE) {
    PanTissue <- FALSE
  }
  if (PanTissue == TRUE & (CollapseMethod == "mean" | CollapseMethod == 
                           "median")) {
    Mean.Ps <- eqtl.data[which(!(is.na(eqtl.data$SNP.Id))), 
    ] %>% dplyr::group_by(SNP.Id) %>% dplyr::summarize(Mean.P = mean(P.Value, 
                                                                     na.rm = TRUE))
    Mean.NESs <- eqtl.data[which(!(is.na(eqtl.data$SNP.Id))), 
    ] %>% dplyr::group_by(SNP.Id) %>% dplyr::summarize(Mean.NES = mean(NES, 
                                                                       na.rm = TRUE))
    Median.Ps <- eqtl.data[which(!(is.na(eqtl.data$SNP.Id))), 
    ] %>% dplyr::group_by(SNP.Id) %>% dplyr::summarize(Median.P = median(P.Value, 
                                                                         na.rm = TRUE))
    Median.NESs <- eqtl.data[which(!(is.na(eqtl.data$SNP.Id))), 
    ] %>% dplyr::group_by(SNP.Id) %>% dplyr::summarize(Median.NES = median(NES, 
                                                                           na.rm = TRUE))
  }
  if (PanTissue == TRUE & CollapseMethod == "meta") {
    if ("N" %in% colnames(eQTL.df) == FALSE & interactive() == 
        TRUE) {
      notrunmeta <- askYesNo(default = TRUE, msg = paste("To complete a PanTissue or MultiTissue analysis using a meta-analysis approach, eQTL.df must contain a column with header N listing the sample size used for each eQTL calculation\nYour eQTL.df does not have a coumn N.\nDo you want to proceed with eQTL meta-analysis assuming all eQTL sample sizes are equal? CAUTION: this may not yield accurate results."))
    }
    else {
      notrunmeta <- "NA"
    }
    if (notrunmeta == FALSE) {
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      print("Stopping analysis")
      stop()
    }
    else {
      if (notrunmeta == TRUE) {
        print("Proceeding with eQTL meta-analysis assuming all eQTL sample sizes are equal.")
        print("CAUTION: this may not yield accurate results")
        eqtl.data$N <- 100
      }
    }
    eqtl.data$sign <- ifelse(eqtl.data$NES < 0, -1, 1)
    eqtl.data$Z <- (qnorm(eqtl.data$P.Value/2, lower.tail = FALSE)) * 
      (eqtl.data$sign)
    eqtl.data$W <- sqrt(eqtl.data$N)
    eqtl.data$ZW <- (eqtl.data$Z) * (eqtl.data$W)
    eqtl.data$W2 <- (eqtl.data$W) * (eqtl.data$W)
    eqtl.data.sumZW <- eqtl.data %>% dplyr::group_by(SNP.Id) %>% 
      dplyr::summarise(SumZW = sum(ZW, na.rm = TRUE))
    eqtl.data.sumW2 <- eqtl.data %>% dplyr::group_by(SNP.Id) %>% 
      dplyr::summarise(SumW2 = sum(W2, na.rm = TRUE))
    eqtl.data.sumW2$sqrtSumW2 <- sqrt(eqtl.data.sumW2$SumW2)
    eqtl.data.sumW2$SumW2 <- NULL
    eqtl.data.sum.final <- dplyr::left_join(eqtl.data.sumW2, 
                                            eqtl.data.sumZW, by = "SNP.Id")
    eqtl.data.sum.final$Z <- (eqtl.data.sum.final$SumZW)/(eqtl.data.sum.final$sqrtSumW2)
    eqtl.data.sum.final$P <- 2 * (pnorm(-abs(eqtl.data.sum.final$Z)))
    eqtl.data$W <- NULL
    eqtl.data$Z <- NULL
    eqtl.data$ZW <- NULL
    eqtl.data$W2 <- NULL
    eqtl.data.sum.final$sqrtSumW2 <- NULL
    eqtl.data.sum.final$SumZW <- NULL
    eqtl.data.sum.final$Z <- NULL
  }
  if (PanTissue == TRUE) {
    eqtl.data <- eqtl.data[which(!(is.na(eqtl.data$SNP.Id))), 
    ] %>% dplyr::group_by(SNP.Id) %>% dplyr::slice(which.min(P.Value))
    eqtl.data$Tissue <- "PanTissue"
  }
  else {
    eqtl.data <- eqtl.data[which(eqtl.data$Tissue == tissue), 
    ]
  }
  if (CollapseMethod == "mean" & PanTissue == TRUE) {
    eqtl.data <- dplyr::left_join(eqtl.data, Mean.Ps, by = "SNP.Id")
    eqtl.data <- dplyr::left_join(eqtl.data, Mean.NESs, 
                                  by = "SNP.Id")
    eqtl.data$P.Value <- eqtl.data$Mean.P
    eqtl.data$NES <- eqtl.data$Mean.NES
    eqtl.data$Mean.P <- NULL
    eqtl.data$Mean.NES <- NULL
  }
  if (CollapseMethod == "median" & PanTissue == TRUE) {
    eqtl.data <- dplyr::left_join(eqtl.data, Median.Ps, 
                                  by = "SNP.Id")
    eqtl.data <- dplyr::left_join(eqtl.data, Median.NESs, 
                                  by = "SNP.Id")
    eqtl.data$P.Value <- eqtl.data$Median.P
    eqtl.data$NES <- eqtl.data$Median.NES
    eqtl.data$Median.P <- NULL
    eqtl.data$Median.NES <- NULL
  }
  if (CollapseMethod == "meta" & PanTissue == TRUE) {
    eqtl.data <- dplyr::left_join(eqtl.data, eqtl.data.sum.final, 
                                  by = "SNP.Id")
    eqtl.data$P.Value <- eqtl.data$P
    eqtl.data$NES <- eqtl.data$sign
    eqtl.data$sign <- NULL
    eqtl.data$N <- NULL
    eqtl.data$P <- NULL
  }
  eqtl.data$P.Value <- ifelse(eqtl.data$P.Value <= 1e-300, 
                              1e-300, eqtl.data$P.Value)
  if (dim(eqtl.data)[1] == 0) 
    stop("Sorry, there are no eQTLs for the tissue", paste(tissue), 
         " with a p-value < sigeQTL")
  eqtl.data <- dplyr::ungroup(eqtl.data)
  gwas.data$SNP <- as.factor(gwas.data$SNP)
  eqtl.data$SNP.Id <- as.factor(eqtl.data$SNP.Id)
  combinedSNPS <- sort(union(levels(gwas.data$SNP), levels(eqtl.data$SNP.Id)))
  Combined.eQTL.GWAS.Data <- dplyr::left_join(dplyr::mutate(gwas.data, 
                                                            SNP = factor(SNP, levels = combinedSNPS)), dplyr::mutate(eqtl.data, 
                                                                                                                     SNP.Id = factor(SNP.Id, levels = combinedSNPS)) %>% 
                                                dplyr::rename(SNP = SNP.Id), by = "SNP")
  if (dim(Combined.eQTL.GWAS.Data)[1] == 0) {
    stop("Sorry, for the gene ", paste(gene), " and the trait ", 
         paste(trait), " there is no overlap between the SNPs in your GWAS.df and eQTL.df")
  }
  Combined.eQTL.GWAS.Data$DirectionOfEffect_GWAS <- ifelse(Combined.eQTL.GWAS.Data$BETA < 
                                                             0, "Negative", ifelse(Combined.eQTL.GWAS.Data$BETA > 
                                                                                     0, "Positive", NA))
  Combined.eQTL.GWAS.Data$DirectionOfEffect_eQTL <- ifelse(Combined.eQTL.GWAS.Data$NES < 
                                                             0, "DOWN", ifelse(Combined.eQTL.GWAS.Data$NES > 0, "UP", 
                                                                               NA))
  Combined.eQTL.GWAS.Data$Congruence <- (Combined.eQTL.GWAS.Data$BETA * 
                                           Combined.eQTL.GWAS.Data$NES)
  Combined.eQTL.GWAS.Data$Congruence <- ifelse(Combined.eQTL.GWAS.Data$Congruence < 
                                                 0, "Incongruent", ifelse(Combined.eQTL.GWAS.Data$Congruence > 
                                                                            0, "Congruent", NA))
  if (congruence == FALSE) {
    Combined.eQTL.GWAS.Data$Congruence <- ifelse(Combined.eQTL.GWAS.Data$Congruence == 
                                                   "Incongruent", "Congruent", ifelse(Combined.eQTL.GWAS.Data$Congruence == 
                                                                                        "Congruent", "Congruent", NA))
  }
  Combined.eQTL.GWAS.Data$NeglogeQTLpValue <- -(log10(Combined.eQTL.GWAS.Data$P.Value))
  Combined.eQTL.GWAS.Data$Neglog10pvalue_GWAS <- -(log10(Combined.eQTL.GWAS.Data$P))
  Combined.eQTL.GWAS.Data <- Combined.eQTL.GWAS.Data[which(!(is.na(Combined.eQTL.GWAS.Data$P))), 
  ]
  Combined.eQTL.GWAS.Data$significance <- ifelse(Combined.eQTL.GWAS.Data$P >= 
                                                   sigpvalue_GWAS, "Non-significant", "Significant")
  Combined.eQTL.GWAS.Data$Congruence[is.na(Combined.eQTL.GWAS.Data$Congruence)] <- "Non-eQTL"
  Combined.eQTL.GWAS.Data <- Combined.eQTL.GWAS.Data %>% dplyr::mutate(Congruence = factor(Congruence, 
                                                                                           levels = c("Non-eQTL", "Congruent", "Incongruent"), 
                                                                                           ordered = TRUE))
  if (dim(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                        Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
  ])[1] == 0) {
    Congruentdata <- FALSE
  }
  else {
    Congruentdata <- TRUE
  }
  if (dim(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                        Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
  ])[1] == 0) {
    Incongruentdata <- FALSE
  }
  else {
    Incongruentdata <- TRUE
  }
  if (isTRUE(LD.df) == FALSE) {
    print("Compiling LD information...")
    Combined.eQTL.GWAS.Data$pvaluemult <- Combined.eQTL.GWAS.Data$P.Value * 
      Combined.eQTL.GWAS.Data$P
    LD.df <- LD.df %>% dplyr::filter(SNP_A %in% levels(Combined.eQTL.GWAS.Data$SNP))
    LD.df <- LD.df %>% dplyr::filter(SNP_B %in% levels(Combined.eQTL.GWAS.Data$SNP))
    if (Congruentdata == TRUE) {
      mostsigsnp.cong <- as.character(Combined.eQTL.GWAS.Data %>% 
                                        dplyr::filter(!is.na(pvaluemult)) %>% dplyr::filter(Congruence == 
                                                                                              "Congruent") %>% dplyr::filter(pvaluemult == 
                                                                                                                               min(pvaluemult)) %>% dplyr::pull(SNP))
      mostsigsnp.cong <- sample(mostsigsnp.cong, 1)
      if (isTRUE(leadSNP) == FALSE) {
        if ((as.character(Combined.eQTL.GWAS.Data %>% 
                          dplyr::filter(SNP == leadSNP) %>% dplyr::pull(Congruence)) == 
             "Congruent") == TRUE) {
          mostsigsnp.cong <- leadSNP
        }
      }
      Combined.eQTL.GWAS.Data <- dplyr::left_join(Combined.eQTL.GWAS.Data, 
                                                  LD.df %>% dplyr::filter(SNP_A == mostsigsnp.cong) %>% 
                                                    dplyr::select(c("SNP_B", "R2")), by = c(SNP = "SNP_B"))
      names(Combined.eQTL.GWAS.Data)[names(Combined.eQTL.GWAS.Data) == 
                                       "R2"] <- "R2cong"
      names(Combined.eQTL.GWAS.Data)[names(Combined.eQTL.GWAS.Data) == 
                                       "SNP_B"] <- "SNP_Bcong"
      Combined.eQTL.GWAS.Data.cong <- subset(Combined.eQTL.GWAS.Data, 
                                             SNP == mostsigsnp.cong)
    }
    if (Incongruentdata == TRUE) {
      mostsigsnp.incong <- as.character(Combined.eQTL.GWAS.Data %>% 
                                          dplyr::filter(!is.na(pvaluemult)) %>% dplyr::filter(Congruence == 
                                                                                                "Incongruent") %>% dplyr::filter(pvaluemult == 
                                                                                                                                   min(pvaluemult)) %>% dplyr::pull(SNP))
      mostsigsnp.incong <- sample(mostsigsnp.incong, 1)
      if (isTRUE(leadSNP) == FALSE) {
        if ((as.character(Combined.eQTL.GWAS.Data %>% 
                          dplyr::filter(SNP == leadSNP) %>% dplyr::pull(Congruence)) == 
             "Incongruent") == TRUE) {
          mostsigsnp.incong <- leadSNP
        }
      }
      Combined.eQTL.GWAS.Data <- dplyr::left_join(Combined.eQTL.GWAS.Data, 
                                                  LD.df %>% dplyr::filter(SNP_A == mostsigsnp.incong) %>% 
                                                    dplyr::select(c("SNP_B", "R2")), by = c(SNP = "SNP_B"))
      names(Combined.eQTL.GWAS.Data)[names(Combined.eQTL.GWAS.Data) == 
                                       "R2"] <- "R2incong"
      names(Combined.eQTL.GWAS.Data)[names(Combined.eQTL.GWAS.Data) == 
                                       "SNP_B"] <- "SNP_Bincong"
      Combined.eQTL.GWAS.Data.incong <- subset(Combined.eQTL.GWAS.Data, 
                                               SNP == mostsigsnp.incong)
    }
    LD.df <- LD.df[which(LD.df$SNP_A %in% Combined.eQTL.GWAS.Data$SNP), 
    ]
    LD.df <- LD.df[which(LD.df$SNP_B %in% Combined.eQTL.GWAS.Data$SNP), 
    ]
    LD.df$R2[LD.df$R2 <= R2min] = NA
    LD.df.matrix <- as.data.frame(tidyr::spread(LD.df[(!duplicated(LD.df[, 
                                                                         c("SNP_B", "SNP_A")])), c("SNP_A", "SNP_B", "R2")], 
                                                SNP_A, R2))
    rownames(LD.df.matrix) <- LD.df.matrix$SNP_B
    LD.df.matrix$SNP_B <- NULL
    LD.df.matrix$startpos <- NA
    LD.df.matrix$stoppos <- NA
    dat2 <- data.frame(matrix(nrow = 2, ncol = ncol(LD.df.matrix)))
    rownames(dat2) <- c("startpos", "stoppos")
    colnames(dat2) <- colnames(LD.df.matrix)
    LD.df.matrix <- dplyr::bind_rows(LD.df.matrix, dat2)
    LD.df.matrix[, c("startpos", "stoppos")] <- NA
    matrix <- as.matrix(LD.df.matrix)
    un1 <- unique(sort(c(colnames(LD.df.matrix), rownames(LD.df.matrix))))
    matrix2 <- matrix(NA, length(un1), length(un1), dimnames = list(un1, 
                                                                    un1))
    matrix2[row.names(matrix), colnames(matrix)] <- matrix
    matrix <- t(matrix2)
    LD.df.matrix <- dplyr::coalesce(as.data.frame(matrix), 
                                    as.data.frame(matrix2))
    SNPsWithLDData1 <- rownames(LD.df.matrix[rowSums(!is.na(LD.df.matrix)) >= 
                                               LDmin, ])
    SNPsWithLDData <- c(SNPsWithLDData1, "startpos", "stoppos")
    if (length(SNPsWithLDData) < 4) {
      stop("Sorry, after filtering the LD.df data by the supplied R2min and LDmin thresholds, fewer than 2 SNPs remain that are also present in GWAS.df")
    }
    LD.df.matrix <- LD.df.matrix[rownames(LD.df.matrix) %in% 
                                   SNPsWithLDData, colnames(LD.df.matrix) %in% SNPsWithLDData]
    LD.df.matrix[is.na(LD.df.matrix)] = 0
    SNPPositions <- unique(dplyr::bind_rows(unique(LD.df[which(LD.df$SNP_A %in% 
                                                                 colnames(LD.df.matrix)), c("SNP_A", "BP_A")]), unique(LD.df[which(LD.df$SNP_B %in% 
                                                                                                                                     colnames(LD.df.matrix)), c("SNP_B", "BP_B")] %>% 
                                                                                                                         dplyr::rename(SNP_A = 1, BP_A = 2))))
    SNPPositions2 <- data.frame(matrix(nrow = 2, ncol = 2))
    colnames(SNPPositions2) <- colnames(SNPPositions)
    SNPPositions2$SNP_A <- c("startpos", "stoppos")
    SNPPositions2$BP_A <- c(startpos, stoppos)
    SNPPositions <- dplyr::bind_rows(SNPPositions, SNPPositions2)
    SNPPositions <- SNPPositions[order(SNPPositions$BP_A), 
    ]
    rownames(SNPPositions) <- SNPPositions$SNP_A
    SNPorder <- SNPPositions[order(SNPPositions$BP_A), ]$SNP_A
    positions <- SNPPositions[order(SNPPositions$BP_A), 
    ]$BP_A
    LD.df.matrix <- LD.df.matrix[SNPorder, SNPorder]
  }
  print("Generating main plot...")
  if (is.na(ylima) == TRUE) {
    ylima <- (max(Combined.eQTL.GWAS.Data %>% dplyr::select(Neglog10pvalue_GWAS), 
                  na.rm = TRUE) + 1)
  }
  
  minpos <- min(Combined.eQTL.GWAS.Data$BP, na.rm = TRUE)
  maxpos <- max(Combined.eQTL.GWAS.Data$BP, na.rm = TRUE)
  p1 <- ggplot2::ggplot(data = Combined.eQTL.GWAS.Data, aes(stroke = 0)) + 
    ggplot2::coord_cartesian(xlim = c(minpos, maxpos), expand = FALSE) + 
    ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, 
                                      Congruence == "Non-eQTL"), shape = 15, color = "black", 
                        alpha = 0.2, aes(x = BP, y = Neglog10pvalue_GWAS)) + 
    ggplot2::xlab("") + ggplot2::ylab(bquote(-log10(P[GWAS]))) + 
    ggplot2::scale_y_continuous(limits = c(NA, ylima)) + 
    ggplot2::ggtitle(paste("GWAS of ", trait, ", colored by eQTL data for ", 
                           gene, "\n(Significance thresholds: GWAS, ", sigpvalue_GWAS, 
                           "; eQTL, ", sigpvalue_eQTL, ")", sep = "")) + ggplot2::scale_shape_manual("GWAS Direction\nof Effect", 
                                                                                                     values = c(Negative = 25, Positive = 24), na.value = 22) + 
    ggplot2::guides(alpha = FALSE, size = guide_legend("eQTL Normalized Effect Size", 
                                                       override.aes = list(shape = 24, color = "black", 
                                                                           fill = "grey"), title.position = "top", order = 2), 
                    shape = guide_legend(title.position = "top", direction = "vertical", 
                                         order = 1, override.aes = list(size = 3, fill = "grey"))) + 
    ggplot2::theme(legend.direction = "horizontal", legend.key = element_rect(fill = NA, 
                                                                              colour = NA, size = 0.25)) + ggplot2::geom_hline(yintercept = -log10(sigpvalue_GWAS), 
                                                                                                                               linetype = "solid", color = "red", size = 0.5) + ggplot2::theme(axis.title.x = element_blank(), 
                                                                                                                                                                                               axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    ggplot2::theme(plot.margin = unit(c(0, 1, -0.8, 0), 
                                      "cm"))+
    ggplot2::theme_bw()+
    
    
    if (congruence == TRUE) {
      if (Congruentdata == TRUE) {
        p1 <- p1 + ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, 
                                                     Congruence == "Congruent"), alpha = 1, aes(x = BP, 
                                                                                                y = Neglog10pvalue_GWAS, fill = NeglogeQTLpValue, 
                                                                                                alpha = 1, shape = DirectionOfEffect_GWAS, size = abs(NES))) + 
          ggplot2::scale_fill_gradient(bquote(atop(-log10(P[eQTL]), 
                                                   paste("Congruous SNPs"))), low = "#6BB7CA", 
                                       high = "#9CD1C8", guide = guide_colorbar(title.position = "top"), 
                                       limits = c(min(Combined.eQTL.GWAS.Data %>% 
                                                        dplyr::select(NeglogeQTLpValue), na.rm = TRUE), 
                                                  max(Combined.eQTL.GWAS.Data %>% dplyr::select(NeglogeQTLpValue), 
                                                      na.rm = TRUE)))
      }
      if (Congruentdata == TRUE & Incongruentdata == TRUE) {
        p1 <- p1 + ggnewscale::new_scale_fill()
      }
      if (Incongruentdata == TRUE) {
        p1 <- p1 + ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, 
                                                     Congruence == "Incongruent"), alpha = 1, aes(x = BP, 
                                                                                                  y = Neglog10pvalue_GWAS, fill = NeglogeQTLpValue, 
                                                                                                  alpha = 1, shape = DirectionOfEffect_GWAS, size = abs(NES))) + 
          ggplot2::scale_fill_gradient(bquote(atop(-log10(P[eQTL]), 
                                                   paste("Incongruous SNPs"))), low = "#E07B54", 
                                       high = "#E1C855", guide = guide_colorbar(title.position = "top"), 
                                       limits = c(min(Combined.eQTL.GWAS.Data %>% 
                                                        dplyr::select(NeglogeQTLpValue), na.rm = TRUE), 
                                                  max(Combined.eQTL.GWAS.Data %>% dplyr::select(NeglogeQTLpValue), 
                                                      na.rm = TRUE)))
      }
    }
  if (Congruentdata == TRUE & Incongruentdata == FALSE & congruence != TRUE) {
    p1 <- p1 + ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, 
                                                 Congruence == "Congruent"), alpha = 1, aes(x = BP, 
                                                                                            y = Neglog10pvalue_GWAS, fill = NeglogeQTLpValue, 
                                                                                            alpha = 1, shape = DirectionOfEffect_GWAS, size = abs(NES))) + 
      ggplot2::scale_fill_viridis_c((bquote(-log10(P[eQTL]))), 
                                    option = "C", guide = guide_colorbar(title.position = "top"), 
                                    limits = c(min(Combined.eQTL.GWAS.Data %>% dplyr::select(NeglogeQTLpValue), 
                                                   na.rm = TRUE), max(Combined.eQTL.GWAS.Data %>% 
                                                                        dplyr::select(NeglogeQTLpValue), na.rm = TRUE)))
  }
  p1 <- p1 + ggplot2::scale_size_continuous(limits = NESeQTLRange)
  p1 = p1+ggplot2::scale_fill_steps(low = "#6BB7CA", high = "#E07B54")
  if (CollapseMethod == "meta") {
    p1 <- p1 + ggplot2::guides(size = FALSE)
  }
  if (isTRUE(LD.df) == FALSE & Congruentdata == TRUE) {
    p1 <- p1 + ggrepel::geom_label_repel(aes(x = BP, y = Neglog10pvalue_GWAS, 
                                             label = ifelse(SNP == mostsigsnp.cong, SNP, ""), 
                                             fontface = "bold"), color = ifelse(congruence == 
                                                                                  T, "#51B1B7", "black"), size = 4, data = Combined.eQTL.GWAS.Data, 
                                         max.overlaps = Inf, force = 5, box.padding = 3, 
                                         min.segment.length = unit(0, "lines"))
  }
  if (isTRUE(LD.df) == FALSE & Incongruentdata == TRUE) {
    p1 <- p1 + ggrepel::geom_label_repel(aes(x = BP, y = Neglog10pvalue_GWAS, 
                                             label = ifelse(SNP == mostsigsnp.incong, SNP, ""), 
                                             fontface = "bold"), color = "#E07B54", size = 4, 
                                         data = Combined.eQTL.GWAS.Data, max.overlaps = Inf, 
                                         force = 5, box.padding = 3, min.segment.length = unit(0, 
                                                                                               "lines"))
  }
  print("Generating gene tracks...")
  if (gbuild == "hg19") {
    hostname <- "https://grch37.ensembl.org"
  }
  if (gbuild == "hg38") {
    hostname <- "https://apr2020.archive.ensembl.org"
  }
  bm <- biomaRt::useMart(host = hostname, biomart = "ENSEMBL_MART_ENSEMBL", 
                         dataset = "hsapiens_gene_ensembl")
  biomTrack <- Gviz::BiomartGeneRegionTrack(genome = gbuild, 
                                            chromosome = median(Combined.eQTL.GWAS.Data$CHR, na.rm = TRUE), 
                                            start = minpos, end = maxpos, filter = list(with_refseq_mrna = TRUE), 
                                            name = "ENSEMBL", background.panel = "gray95", biomart = bm, 
                                            margin = c(-3, -3))
  gtrack <- Gviz::GenomeAxisTrack(fontcolor = "#000000", fontsize = 14, 
                                  margin = c(-3, -3))
  genetracks <- patchwork::wrap_elements(panel = (grid::grid.grabExpr(Gviz::plotTracks(list(biomTrack, 
                                                                                            gtrack), collapseTranscripts = "meta", transcriptAnnotation = "symbol", 
                                                                                       chromosome = median(Combined.eQTL.GWAS.Data$CHR, na.rm = TRUE), 
                                                                                       from = min(Combined.eQTL.GWAS.Data$BP, na.rm = TRUE), 
                                                                                       to = max(Combined.eQTL.GWAS.Data$BP, na.rm = TRUE), 
                                                                                       showTitle = FALSE, labelPos = "below", distFromAxis = 10, 
                                                                                       innermargin = 0, maxHeight = (genometrackheight * 10), 
                                                                                       minHeight = (genometrackheight * 10), sizes = c(genometrackheight, 
                                                                                                                                       1), margin = c(-3, -3)))))
  if (isTRUE(LD.df) == FALSE) {
    if (length(SNPsWithLDData) > 1000 & interactive()) {
      notrun <- askYesNo(default = TRUE, msg = paste(sep = "", 
                                                     length(SNPsWithLDData), " variants being used to generate LDHeatMap.\nUsing more than 1000 variants can take a long time to run.\nYou can increase the values for LDmin and R2min to use fewer variants.\nDo you want to continue with ", 
                                                     length(SNPsWithLDData), " variants?"))
    }
    else {
      notrun <- TRUE
    }
    if (notrun == FALSE) {
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      print("Stopping analysis")
      stop()
    }
    else {
      if (notrun == TRUE) {
        print(paste(sep = "", "Generating LDHeatMap with ", 
                    length(SNPsWithLDData), " variants"))
      }
      if (LDcolor == "color") {
        colorscale <- c(viridisLite::viridis(30, option = "C", 
                                             direction = -1), "white")
      }
      else {
        if (LDcolor == "black") {
          colorscale <- c("grey10", "grey20", "grey30", 
                          "grey40", "grey50", "grey60", "grey70", 
                          "grey80", "grey90", "grey100", "white")
        }
      }
      LDmap <- LDheatmap::LDheatmap(as.matrix(LD.df.matrix), 
                                    genetic.distances = positions, color = colorscale, 
                                    flip = TRUE, add.map = TRUE, title = "", geneMapLabelX = NA, 
                                    geneMapLabelY = NA, newpage = FALSE)
      dev.off()
    }
    if (length(SNPsWithLDData) >= 50) {
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$width <- unit(0.69, 
                                                                             "snpc")
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$height <- unit(0.69, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$width <- unit(0.69, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$height <- unit(0.69, 
                                                                               "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$width <- unit(0.69, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$height <- unit(0.69, 
                                                                               "snpc")
    }
    if (length(SNPsWithLDData) < 50 & length(SNPsWithLDData) >= 
        25) {
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$width <- unit(0.7, 
                                                                             "snpc")
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$height <- unit(0.7, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$width <- unit(0.7, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$height <- unit(0.7, 
                                                                               "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$width <- unit(0.7, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$height <- unit(0.7, 
                                                                               "snpc")
    }
    if (length(SNPsWithLDData) < 25 & length(SNPsWithLDData) >= 
        15) {
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$width <- unit(0.725, 
                                                                             "snpc")
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$height <- unit(0.725, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$width <- unit(0.725, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$height <- unit(0.725, 
                                                                               "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$width <- unit(0.725, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$height <- unit(0.725, 
                                                                               "snpc")
    }
    if (length(SNPsWithLDData) < 15 & length(SNPsWithLDData) >= 
        9) {
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$width <- unit(0.75, 
                                                                             "snpc")
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$height <- unit(0.75, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$width <- unit(0.75, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$height <- unit(0.75, 
                                                                               "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$width <- unit(0.75, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$height <- unit(0.75, 
                                                                               "snpc")
    }
    if (length(SNPsWithLDData) < 9 & length(SNPsWithLDData) >= 
        5) {
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$width <- unit(0.8, 
                                                                             "snpc")
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$height <- unit(0.8, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$width <- unit(0.8, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$height <- unit(0.8, 
                                                                               "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$width <- unit(0.8, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$height <- unit(0.8, 
                                                                               "snpc")
    }
    if (length(SNPsWithLDData) < 5) {
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$width <- unit(0.875, 
                                                                             "snpc")
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$height <- unit(0.875, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$width <- unit(0.875, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$height <- unit(0.875, 
                                                                               "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$width <- unit(0.875, 
                                                                              "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$height <- unit(0.875, 
                                                                               "snpc")
    }
    LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$y <- unit(0.875, 
                                                                       "npc")
    LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$y <- unit(0.875, 
                                                                        "npc")
    LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$y <- unit(0.875, 
                                                                        "npc")
    LDmap$LDheatmapGrob$children$Key$vp$y <- unit(0.65, 
                                                  "npc")
    LDmap$LDheatmapGrob$children$Key$vp$x <- unit(0.925, 
                                                  "npc")
    LDmap$LDheatmapGrob$children$Key$vp$width <- unit(0.15, 
                                                      "npc")
    LDmap$LDheatmapGrob$children$Key$vp$height <- unit(0.05, 
                                                       "npc")
    LDmap$flipVP$justification <- c("center", "top")
    g.ld <- ggplotGrob(ggplotify::as.ggplot(LDmap$LDheatmapGrob) + 
                         theme(plot.margin = unit(c(-0.15, 0, -1.8, 0), "npc")))
    g.p1 <- ggplotGrob(p1)
    g.ld$widths <- g.p1$widths
  }
  print("Generating eQTL enrichment plot...")
  Combined.eQTL.GWAS.Data$iseQTL <- Combined.eQTL.GWAS.Data$Congruence
  Combined.eQTL.GWAS.Data$iseQTL <- ifelse(Combined.eQTL.GWAS.Data$iseQTL == 
                                             "Non-eQTL", "Non-eQTL", "eQTL")
  if (nrow(as.table(table(Combined.eQTL.GWAS.Data$iseQTL, 
                          Combined.eQTL.GWAS.Data$significance))) < 2 | ncol(as.table(table(Combined.eQTL.GWAS.Data$iseQTL, 
                                                                                            Combined.eQTL.GWAS.Data$significance))) < 2) {
    NoFisher <- TRUE
    print("Not enough data to compute enrichment significance for Plot C")
  }
  if (NoFisher == FALSE) {
    fisher <- fisher.test(table(Combined.eQTL.GWAS.Data$iseQTL, 
                                Combined.eQTL.GWAS.Data$significance))
  }
  if (NoFisher == FALSE) {
    fpvalue <- fisher$p.value
  }
  p2 <- ggplot2::ggplot(Combined.eQTL.GWAS.Data) + ggplot2::aes(x = significance, 
                                                                y = 1, fill = (Congruence)) + ggplot2::geom_bar(stat = "identity", 
                                                                                                                position = "fill") + ggplot2::ggtitle(paste("Enrichment of eQTLs among\nGWAS-significant SNPs")) + 
    ggplot2::ylab("Proportion of SNPs\nthat are eQTLs") + 
    ggplot2::xlab(paste("GWAS significance\n(threshold p <", 
                        sigpvalue_GWAS, ")")) + ggplot2::ylim(0, 1.2) + 
    if (NoFisher == FALSE) {
      ggpubr::geom_signif(y_position = c(1.1, 1.1), xmin = c("Non-significant"), 
                          xmax = c("Significant"), annotation = (paste("p =", 
                                                                       formatC(fpvalue, format = "e", digits = 2))), 
                          tip_length = 0.05)
    }
  if (Congruentdata == TRUE & Incongruentdata == TRUE) {
    p2 <- p2 + ggplot2::scale_fill_manual(labels = c(Congruent = "Congruent eQTL", 
                                                     Incongruent = "Incongruent eQTL", `Non-eQTL` = "Non-eQTL"), 
                                          values = c(Congruent = "#000099", Incongruent = "#990000", 
                                                     `Non-eQTL` = "#C0C0C0")) + ggplot2::guides(fill = guide_legend(title = NULL))
  }
  if (Congruentdata == TRUE & Incongruentdata == FALSE & congruence == 
      TRUE) {
    p2 <- p2 + ggplot2::scale_fill_manual(labels = c(Congruent = "Congruent eQTL", 
                                                     `Non-eQTL` = "Non-eQTL", ), values = c(Congruent = "#000099", 
                                                                                            `Non-eQTL` = "#C0C0C0")) + ggplot2::guides(fill = guide_legend(title = NULL))
  }
  if (Congruentdata == TRUE & Incongruentdata == FALSE & congruence == 
      FALSE) {
    p2 <- p2 + ggplot2::scale_fill_manual(labels = c(Congruent = "eQTL", 
                                                     `Non-eQTL` = "Non-eQTL"), values = c(Congruent = "#ffee00", 
                                                                                          `Non-eQTL` = "#360f70")) + ggplot2::guides(fill = guide_legend(""))
  }
  if (Congruentdata == FALSE & Incongruentdata == TRUE) {
    p2 <- p2 + ggplot2::scale_fill_manual(labels = c(Incongruent = "Incongruent eQTL", 
                                                     `Non-eQTL` = "Non-eQTL"), values = c(Incongruent = "#990000", 
                                                                                          `Non-eQTL` = "#C0C0C0")) + ggplot2::guides(fill = guide_legend(title = NULL))
  }
  print("Generating P-P plot...")
  p3 <- ggplot2::ggplot(data = Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES)), 
  ]) + ggplot2::guides(color = guide_legend("Direction of Effect")) + 
    ggplot2::xlab((bquote(-log10(P[eQTL])))) + ggplot2::ylab((bquote(-log10(P[GWAS])))) + 
    ggplot2::scale_y_continuous(limits = c(NA, ylimd)) + 
    ggplot2::scale_x_continuous(limits = c(NA, xlimd)) + 
    ggplot2::theme(legend.position = "right") + theme(legend.spacing.y = unit(0.1, 
                                                                              "cm")) + theme(legend.key = element_rect(fill = NA, 
                                                                                                                       colour = NA, size = 0.25))
  if (isTRUE(LD.df) == FALSE & Congruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                                          Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
  ]) >= 2) {
    pearson.congruent <- cor.test(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                  Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
    ]$NeglogeQTLpValue, Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                        Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
    ]$Neglog10pvalue_GWAS, method = "pearson")
    p3.1 <- p3 + ggplot2::geom_point(data = Combined.eQTL.GWAS.Data[which(Combined.eQTL.GWAS.Data$Congruence == 
                                                                            "Congruent" & is.na(Combined.eQTL.GWAS.Data$R2cong) == 
                                                                            TRUE), ], aes(y = Neglog10pvalue_GWAS, x = NeglogeQTLpValue, 
                                                                                          fill = R2cong), shape = 21, size = 3) + ggplot2::geom_point(data = Combined.eQTL.GWAS.Data[which(Combined.eQTL.GWAS.Data$Congruence == 
                                                                                                                                                                                             "Congruent" & is.na(Combined.eQTL.GWAS.Data$R2cong) == 
                                                                                                                                                                                             FALSE), ], aes(y = Neglog10pvalue_GWAS, x = NeglogeQTLpValue, 
                                                                                                                                                                                                            fill = R2cong), shape = 21, size = 3) + ggplot2::geom_smooth(data = Combined.eQTL.GWAS.Data[which(Combined.eQTL.GWAS.Data$Congruence == 
                                                                                                                                                                                                                                                                                                                "Congruent"), ], aes(y = Neglog10pvalue_GWAS, x = NeglogeQTLpValue), 
                                                                                                                                                                                                                                                                         color = "black", method = "lm", formula = (y ~ x))
    p3.1 <- p3.1 + ggplot2::geom_point(data = Combined.eQTL.GWAS.Data.cong, 
                                       aes(y = Neglog10pvalue_GWAS, x = NeglogeQTLpValue, 
                                           shape = SNP), size = 3, fill = "#33FF33") + 
      scale_shape_manual(values = 23) + guides(shape = guide_legend(order = 2, 
                                                                    title = NULL, override.aes = list(size = 3, shape = 23, 
                                                                                                      fill = "#33FF33")))
    p3.1 <- p3.1 + ggplot2::geom_text(x = -Inf, y = Inf, 
                                      label = paste(sep = "", "r = ", round(pearson.congruent$estimate, 
                                                                            2), "\np = ", formatC(pearson.congruent$p.value, 
                                                                                                  format = "e", digits = 2)), hjust = -0.05, vjust = 1.1, 
                                      color = "black")
  }
  else {
    if (isTRUE(LD.df) == FALSE) {
      print("Not enough data to generate P-P plot for Congruent eQTLs")
    }
  }
  if (isTRUE(LD.df) == FALSE & Incongruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                                            Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
  ]) >= 2) {
    pearson.incongruent <- cor.test(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                    Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
    ]$NeglogeQTLpValue, Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                        Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
    ]$Neglog10pvalue_GWAS, method = "pearson")
    p3.2 <- p3 + ggplot2::geom_point(data = Combined.eQTL.GWAS.Data[which(Combined.eQTL.GWAS.Data$Congruence == 
                                                                            "Incongruent" & is.na(Combined.eQTL.GWAS.Data$R2cong) == 
                                                                            TRUE), ], aes(y = Neglog10pvalue_GWAS, x = NeglogeQTLpValue, 
                                                                                          fill = R2incong), shape = 21, size = 3) + ggplot2::geom_point(data = Combined.eQTL.GWAS.Data[which(Combined.eQTL.GWAS.Data$Congruence == 
                                                                                                                                                                                               "Incongruent" & is.na(Combined.eQTL.GWAS.Data$R2cong) == 
                                                                                                                                                                                               FALSE), ], aes(y = Neglog10pvalue_GWAS, x = NeglogeQTLpValue, 
                                                                                                                                                                                                              fill = R2incong), shape = 21, size = 3) + ggplot2::scale_fill_gradient(bquote(R^2 ~ 
                                                                                                                                                                                                                                                                                              "with" ~ {
                                                                                                                                                                                                                                                                                                .(mostsigsnp.incong)
                                                                                                                                                                                                                                                                                              }), limits = c(0, 1), breaks = c(0.2, 0.4, 0.6, 0.8), 
                                                                                                                                                                                                                                                                                     low = "#990000", high = "#FFCC33", na.value = "grey80", 
                                                                                                                                                                                                                                                                                     guide = guide_colorbar(order = 1, direction = "horizontal", 
                                                                                                                                                                                                                                                                                                            title.position = "top", label.position = "bottom")) + 
      ggplot2::geom_smooth(data = Combined.eQTL.GWAS.Data[which(Combined.eQTL.GWAS.Data$Congruence == 
                                                                  "Incongruent"), ], aes(y = Neglog10pvalue_GWAS, 
                                                                                         x = NeglogeQTLpValue), color = "black", method = "lm", 
                           formula = (y ~ x))
    p3.2 <- p3.2 + ggplot2::geom_point(data = Combined.eQTL.GWAS.Data.incong, 
                                       aes(y = Neglog10pvalue_GWAS, x = NeglogeQTLpValue, 
                                           shape = SNP), size = 3, fill = "#33FF33") + 
      scale_shape_manual(values = 23) + guides(shape = guide_legend(order = 2, 
                                                                    title = NULL, override.aes = list(size = 3, shape = 23, 
                                                                                                      fill = "#33FF33"))) + ggplot2::geom_text(x = -Inf, 
                                                                                                                                               y = Inf, label = paste(sep = "", "r = ", round(pearson.incongruent$estimate, 
                                                                                                                                                                                              2), "\np = ", formatC(pearson.incongruent$p.value, 
                                                                                                                                                                                                                    format = "e", digits = 2)), hjust = -0.05, vjust = 1.1, 
                                                                                                                                               color = "black") + ggplot2::ggtitle(paste("P-P plot, incongruent SNPs"))
  }
  else {
    if (isTRUE(LD.df) == FALSE & congruence == TRUE) {
      print("Not enough data to generate P-P plot for Incongruent eQTLs")
    }
  }
  if (isTRUE(LD.df) == FALSE & congruence == FALSE) {
    p3.1 <- p3.1 + ggplot2::scale_fill_viridis_c(bquote(R^2 ~ 
                                                          "with" ~ {
                                                            .(mostsigsnp.cong)
                                                          }), limits = c(0, 1), breaks = c(0.2, 0.4, 0.6, 0.8), 
                                                 option = "C", na.value = "grey80", guide = guide_colorbar(order = 1, 
                                                                                                           direction = "horizontal", title.position = "top", 
                                                                                                           label.position = "bottom", title.vjust = 0)) + 
      ggplot2::ggtitle(paste("P-P plot"))
  }
  else {
    if (isTRUE(LD.df) == FALSE & congruence == TRUE) {
      p3.1 <- p3.1 + ggplot2::scale_fill_gradient(bquote(R^2 ~ 
                                                           "with" ~ {
                                                             .(mostsigsnp.cong)
                                                           }), limits = c(0, 1), breaks = c(0.2, 0.4, 0.6, 
                                                                                            0.8), low = "#000099", high = "#33FFFF", na.value = "grey80", 
                                                  guide = guide_colorbar(order = 1, direction = "horizontal", 
                                                                         title.position = "top", label.position = "bottom", 
                                                                         title.vjust = 0)) + ggplot2::ggtitle(paste("P-P plot, congruent SNPs"))
    }
  }
  if (isTRUE(LD.df) == TRUE) {
    if (Congruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                   Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
    ]) >= 2) {
      pearson.congruent <- cor.test(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                    Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
      ]$NeglogeQTLpValue, Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                          Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
      ]$Neglog10pvalue_GWAS, method = "pearson")
      if (congruence == TRUE) {
        p3 <- p3 + ggplot2::geom_smooth(data = Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                               Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
        ], aes(y = Neglog10pvalue_GWAS, x = NeglogeQTLpValue, 
               color = Congruence), method = "lm", formula = (y ~ 
                                                                x))
      }
      else {
        if (congruence == FALSE) 
          p3 <- p3 + ggplot2::geom_smooth(data = Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                                 Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
          ], aes(y = Neglog10pvalue_GWAS, x = NeglogeQTLpValue, 
                 color = Congruence), method = "lm", formula = (y ~ 
                                                                  x), show.legend = FALSE, color = "#E1C855")
      }
      if (congruence == TRUE) {
        p3 <- p3 + ggplot2::geom_point(data = Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                              Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
        ], aes(y = Neglog10pvalue_GWAS, x = NeglogeQTLpValue, 
               color = Congruence))
      }
      else {
        p3 <- p3 + ggplot2::geom_point(data = Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                              Combined.eQTL.GWAS.Data$Congruence == "Congruent"), 
        ], aes(y = Neglog10pvalue_GWAS, x = NeglogeQTLpValue), 
        color = "#51B1B7")
        p3 <- p3 + ggplot2::geom_text(x = -Inf, y = Inf, 
                                      label = paste(sep = "", "r = ", round(pearson.congruent$estimate, 
                                                                            3), "\np = ", formatC(pearson.congruent$p.value, 
                                                                                                  format = "e", digits = 2)), color = "#E07B54", 
                                      hjust = -0.05, vjust = 1.1)
      }
    }
    else {
      print("Not enough data to generate P-P plot for Congruent eQTLs")
    }
    if (Incongruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                     Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
    ]) >= 2) {
      pearson.incongruent <- cor.test(Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                      Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
      ]$NeglogeQTLpValue, Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                          Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
      ]$Neglog10pvalue_GWAS, method = "pearson")
      p3 <- p3 + ggplot2::geom_point(data = Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                            Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
      ], aes(y = Neglog10pvalue_GWAS, x = NeglogeQTLpValue, 
             color = Congruence)) + ggplot2::geom_smooth(data = Combined.eQTL.GWAS.Data[which(!is.na(Combined.eQTL.GWAS.Data$NES) & 
                                                                                                Combined.eQTL.GWAS.Data$Congruence == "Incongruent"), 
             ], aes(y = Neglog10pvalue_GWAS, x = NeglogeQTLpValue, 
                    color = Congruence), method = "lm", formula = (y ~ 
                                                                     x)) + ggplot2::geom_text(x = -Inf, y = Inf, 
                                                                                              label = paste(sep = "", "\n\nr = ", round(pearson.incongruent$estimate, 
                                                                                                                                        3), "\np = ", formatC(pearson.incongruent$p.value, 
                                                                                                                                                              format = "e", digits = 2)), color = "#990000", 
                                                                                              hjust = -0.05, vjust = 1.1)
      p3 <- p3 + ggplot2::geom_text(x = -Inf, y = Inf, 
                                    label = paste(sep = "", "r = ", round(pearson.congruent$estimate, 
                                                                          3), "\np = ", formatC(pearson.congruent$p.value, 
                                                                                                format = "e", digits = 2)), color = "#000099", 
                                    hjust = -0.05, vjust = 1.1)
    }
    else {
      if (congruence == "TRUE") {
        print("Not enough data to generate P-P Plot for Incongruent eQTLs")
      }
    }
    if (Congruentdata == TRUE & Incongruentdata == TRUE) {
      p3 <- p3 + scale_color_manual(values = c(Congruent = "#000099", 
                                               Incongruent = "#990000"))
    }
    if (Congruentdata == TRUE & Incongruentdata == FALSE) {
      p3 <- p3 + scale_color_manual(values = c(Congruent = "#000099"))
    }
    if (Congruentdata == FALSE & Incongruentdata == TRUE) {
      p3 <- p3 + scale_color_manual(values = c(Incongruent = "#990000"))
    }
    p3 <- p3 + ggplot2::ggtitle(paste("P-P plot"))
    p3 <- p3+ggplot2::theme_bw()+ggplot2::scale_color_manual(values ="#51B1B7" )
  }
  p2<-p2+ggplot2::theme_bw()+ggplot2::scale_fill_manual(values = c("#51B1B7","#E07B54"))
  print("Merging and plotting...")
  CollapseMethod <- ifelse(CollapseMethod == "min", "minimum value", 
                           ifelse(CollapseMethod == "mean", "mean value", ifelse(CollapseMethod == 
                                                                                   "median", "median value", ifelse(CollapseMethod == 
                                                                                                                      "meta", "meta-analysis", NA))))
  if (PanTissue == TRUE & length(tissue) >= 2) {
    tissue <- "MultiTissue"
  }
  if (PanTissue == TRUE & tissue == "all") {
    tissue <- "PanTissue"
  }
  if (PanTissue == TRUE) 
    tissuetitle <- paste(tissue, "analysis, eQTLs collapsed by", 
                         CollapseMethod)
  else tissuetitle <- paste("In", tissue)
  if (isTRUE(LD.df) == FALSE & Incongruentdata == FALSE) {
    p4 <- ((p1 + genetracks + g.ld + patchwork::plot_layout(nrow = 3, 
                                                            byrow = FALSE, heights = c(2, (genometrackheight/5), 
                                                                                       NA))) | (p2/patchwork::plot_spacer()/p3.1/patchwork::plot_spacer()/patchwork::plot_spacer() + 
                                                                                                  patchwork::plot_layout(heights = c(5, 0.5, 5, 0.5, 
                                                                                                                                     5)))) + patchwork::plot_layout(ncol = 2, widths = c(2.5, 
                                                                                                                                                                                         1)) + patchwork::plot_annotation(tag_levels = "A", 
                                                                                                                                                                                                                          tag_suffix = ".") & theme(plot.tag = element_text(size = 18), 
                                                                                                                                                                                                                                                    text = element_text(size = 12), plot.tag.position = c(0, 
                                                                                                                                                                                                                                                                                                          0.995))
  }
  if (isTRUE(LD.df) == FALSE & Congruentdata == FALSE) {
    p4 <- ((p1 + genetracks + g.ld + patchwork::plot_layout(nrow = 3, 
                                                            byrow = FALSE, heights = c(2, (genometrackheight/5), 
                                                                                       NA))) | (p2/patchwork::plot_spacer()/p3.2/patchwork::plot_spacer()/patchwork::plot_spacer() + 
                                                                                                  patchwork::plot_layout(heights = c(5, 0.5, 5, 0.5, 
                                                                                                                                     5)))) + patchwork::plot_layout(ncol = 2, widths = c(2.5, 
                                                                                                                                                                                         1)) + patchwork::plot_annotation(tag_levels = "A", 
                                                                                                                                                                                                                          tag_suffix = ".") & theme(plot.tag = element_text(size = 18), 
                                                                                                                                                                                                                                                    text = element_text(size = 12), plot.tag.position = c(0, 
                                                                                                                                                                                                                                                                                                          0.995))
  }
  if (isTRUE(LD.df) == FALSE & Congruentdata == TRUE & Incongruentdata == 
      TRUE) {
    p4 <- ((p1 + genetracks + g.ld + patchwork::plot_layout(nrow = 3, 
                                                            byrow = FALSE, heights = c(2, (genometrackheight/5), 
                                                                                       NA))) | (p2/patchwork::plot_spacer()/p3.1/patchwork::plot_spacer()/p3.2) + 
             patchwork::plot_layout(heights = c(5, 0.5, 5, 0.5, 
                                                5))) + patchwork::plot_layout(ncol = 2, widths = c(2.5, 
                                                                                                   1)) + patchwork::plot_annotation(tag_levels = "A", 
                                                                                                                                    tag_suffix = ".") & theme(plot.tag = element_text(size = 18), 
                                                                                                                                                              text = element_text(size = 12), plot.tag.position = c(0, 
                                                                                                                                                                                                                    0.995))
  }
  if (isTRUE(LD.df) == TRUE) {
    p4 <- (p1 + genetracks + patchwork::plot_spacer() + 
             (p2 + p3 + patchwork::plot_layout(ncol = 2, widths = c(2, 
                                                                    3))) + patchwork::plot_layout(ncol = 1, height = c(4, 
                                                                                                                       genometrackheight, 0.1, 2)) + patchwork::plot_annotation(tag_levels = "A", 
                                                                                                                                                                                tag_suffix = ".") & theme(plot.tag = element_text(size = 18), 
                                                                                                                                                                                                          text = element_text(size = 12)))
  }
  p4 =  p4+ggplot2::scale_color_manual(values = c( "#6BB7CA"))
  
  
  pfinal <- p4 + patchwork::plot_annotation(title = paste("eQTpLot analysis for ", 
                                                          trait, " and ", gene, "\n", tissuetitle, "\n", sep = ""), 
                                            theme = theme(plot.title = element_text(size = 19, face = "bold")))
  if (isTRUE(LD.df) == FALSE & wi == "wi") {
    wi <- 14
  }
  if (isTRUE(LD.df) == TRUE & wi == "wi") {
    wi <- 12
  }
  if (isTRUE(LD.df) == FALSE) {
    hgt <- ((1.3 * (wi) - 0.01 * (wi)^2 - 5.5))
  }
  if (isTRUE(LD.df) == TRUE) {
    hgt <- wi * 1.1
  }
  if (congruence == TRUE) {
    congruence <- "WithCongruenceData"
  }
  else {
    congruence <- "WithoutCongruenceData"
  }
  if (isTRUE(LD.df) == FALSE) {
    LDinfo <- "WithLinkageData"
  }
  else {
    LDinfo <- "WithoutLinkageData"
  }
  if (saveplot == TRUE) {
    ggsave(pfinal, filename = paste(gene, trait, tissue, 
                                    congruence, LDinfo, "eQTpLot", "png", sep = "."), 
           dpi = res, units = "in", height = hgt, width = wi)
  }
  if (getplot == TRUE) {
    return(pfinal)
  }
  
}
