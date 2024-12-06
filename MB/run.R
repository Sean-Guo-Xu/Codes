run<-function (gene.vec, ori.tbl, sub.tbl, mat, assay.use = "counts", 
          model = c("nb", "zinb", "gaussian", "auto", "qgam"), k = 6, 
          knots = c(0:5/5), fix.weight = TRUE, aicdiff = 10, seed = 123, 
          quant = 0.5, usebam = FALSE, seurat.assay = "RNA", mc.cores = 4, 
          mc.preschedule = TRUE, SIMPLIFY = TRUE) 
{
  set.seed(seed)
  expv.quantile <- gam.fit <- NULL
  BPPARAM <- BiocParallel::bpparam()
  BPPARAM$workers <- mc.cores
  if (!is.null(sub.tbl)) {
    is.subset_all <- sapply(sub.tbl, function(x) {
      all(x$cell %in% ori.tbl$cell)
    })
    if (!all(is.subset_all)) {
      stop("Cells in ", paste0(which(!is.subset_all), 
                               " is not a subset of ori.tbl."))
    }
  }
  res <- BiocParallel::bplapply(gene.vec, function(x, ...) {
    cur_res <- tryCatch( 
                          expr = mypseudotimeDE(gene = x, ori.tbl = ori.tbl, 
                                            sub.tbl = sub.tbl, mat = mat, model = model, 
                                            assay.use = assay.use, seurat.assay = seurat.assay), 
                        error = function(e) {
                          list(fix.pv = NA, emp.pv = NA, para.pv = NA, 
                               ad.pv = NA, rank = NA, test.statistics = NA, 
                               gam.fit = NA, zinf = NA, aic = NA, expv.quantile = NA, 
                               expv.mean = NA, expv.zero = NA)
                        })
    cur_res
  }, assay.use = assay.use, k = k, knots = knots, fix.weight = fix.weight, 
  aicdiff = aicdiff, quant = quant, usebam = usebam, seurat.assay = seurat.assay, 
  BPPARAM =BPPARAM)
  if (SIMPLIFY) {
    print("ending")
    res <- simplify2array(res)
    res <- t(res)
    rownames(res) <- gene.vec
    res <- tibble::as_tibble(res, rownames = "gene")
    res <- tidyr::unnest(res, cols = !(gam.fit | expv.quantile))
  }
  res
}
