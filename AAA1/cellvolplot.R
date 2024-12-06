cellvolplot<-function (diffData = NULL, myMarkers = NULL, order.by = c("avg_log2FC"), 
          log2FC.cutoff = 0.25, pvalue.cutoff = 0.05, adjustP.cutoff = 0.01, 
          topGeneN = 5, col.type = "updown", back.col = "grey93", 
          pSize = 0.75, aesCol = c("#0099CC", "#CC3333"), labelsize=0.2,legend.position = c(0.7, 
                                                                              0.9), base_size = 14, tile.col = jjAnno::useMyCol("paired", 
                                                                                                                                n = 9), cluster.order = NULL, polar = FALSE, expand = c(-1, 
                                                                                                                                                                                        1), flip = FALSE, celltypeSize = 3, ...) {
  diff.marker <- diffData %>% dplyr::filter(abs(avg_log2FC) >= 
                                              log2FC.cutoff & p_val < pvalue.cutoff)
  diff.marker <- diff.marker %>% dplyr::mutate(type = ifelse(avg_log2FC >= 
                                                               log2FC.cutoff, "sigUp", "sigDown")) %>% dplyr::mutate(type2 = ifelse(p_val_adj < 
                                                                                                                                      adjustP.cutoff, paste("adjust Pvalue < ", adjustP.cutoff, 
                                                                                                                                                            sep = ""), paste("adjust Pvalue >= ", adjustP.cutoff, 
                                                                                                                                                                             sep = "")))
  if (!is.null(cluster.order)) {
    diff.marker$cluster <- factor(diff.marker$cluster, levels = cluster.order)
  }
  back.data <- purrr::map_df(unique(diff.marker$cluster), 
                             function(x) {
                               tmp <- diff.marker %>% dplyr::filter(cluster == 
                                                                      x)
                               new.tmp <- data.frame(cluster = x, min = min(tmp$avg_log2FC) - 
                                                       0.2, max = max(tmp$avg_log2FC) + 0.2)
                               return(new.tmp)
                             })
  top.marker.tmp <- diff.marker %>% dplyr::group_by(cluster)
  top.marker.max <- top.marker.tmp %>% dplyr::slice_max(n = topGeneN, 
                                                        order_by = get(order.by))
  top.marker.min <- top.marker.tmp %>% dplyr::slice_min(n = topGeneN, 
                                                        order_by = get(order.by))
  top.marker <- rbind(top.marker.max, top.marker.min)
  if (!is.null(myMarkers)) {
    top.marker <- diff.marker %>% dplyr::filter(gene %in% 
                                                  myMarkers)
  }
  else {
    top.marker <- top.marker
  }
  p1 <- ggplot2::ggplot(diff.marker, ggplot2::aes(x = cluster, 
                                                  y = avg_log2FC)) + ggplot2::geom_col(data = back.data, 
                                                                                       ggplot2::aes(x = cluster, y = min), fill = back.col) + 
    ggplot2::geom_col(data = back.data, ggplot2::aes(x = cluster, 
                                                     y = max), fill = back.col)
  if (col.type == "updown") {
    p2 <- p1 + ggplot2::geom_jitter(ggplot2::aes(color = type), 
                                    size = pSize) + ggplot2::scale_color_manual(values = c(sigDown = aesCol[1], 
                                                                                           sigUp = aesCol[2]))
  }
  else if (col.type == "adjustP") {
    p2 <- p1 + ggplot2::geom_jitter(ggplot2::aes(color = type2), 
                                    size = pSize) + ggplot2::scale_color_manual(values = c(aesCol[2], 
                                                                                           aesCol[1]))
  }
  p3 <- p2 + ggplot2::scale_y_continuous(n.breaks = 6) + ggplot2::theme_classic(base_size = base_size) + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                   legend.position = legend.position, legend.title = ggplot2::element_blank(), 
                   legend.background = ggplot2::element_blank()) + 
    ggplot2::xlab("Clusters") + ggplot2::ylab("Average log2FoldChange") + 
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))
  p4 <- p3 + ggplot2::geom_tile(ggplot2::aes(x = cluster, 
                                             y = 0, fill = cluster), color = "black", height = log2FC.cutoff * 
                                  2, alpha = 0.3, show.legend = FALSE) + ggplot2::scale_fill_manual(values = tile.col) + 
    ggrepel::geom_text_repel(data = top.marker, ggplot2::aes(x = cluster, 
                                                             y = avg_log2FC, label = gene,nudge_y = ifelse(avg_log2FC > 0,log2FC.cutoff*4, -(log2FC.cutoff*4))), max.overlaps = 50, size = labelsize,
                             ...)
  if (polar == TRUE) {
    p5 <- p4 + geomtextpath::geom_textpath(ggplot2::aes(x = cluster, 
                                                        y = 0, label = cluster)) + ggplot2::scale_y_continuous(n.breaks = 6, 
                                                                                                               expand = ggplot2::expansion(mult = expand)) + ggplot2::theme_void(base_size = base_size) + 
      ggplot2::theme(legend.position = legend.position, 
                     legend.title = ggplot2::element_blank()) + ggplot2::coord_polar(clip = "off", 
                                                                                     theta = "x")
  }
  else {
    if (flip == TRUE) {
      p5 <- p4 + ggplot2::scale_y_continuous(n.breaks = 6) + 
        ggplot2::geom_label(ggplot2::aes(x = cluster, 
                                         y = 0, label = cluster), size = celltypeSize) + 
        ggplot2::theme(axis.line.y = ggplot2::element_blank(), 
                       axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) + 
        ggplot2::coord_flip()
    }
    else {
      p5 <- p4 + ggplot2::scale_y_continuous(n.breaks = 6) + 
        ggplot2::geom_text(ggplot2::aes(x = cluster, 
                                        y = 0, label = cluster), size = celltypeSize) + 
        ggplot2::theme(axis.line.x = ggplot2::element_blank(), 
                       axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    }
  }
  return(p5)
}
