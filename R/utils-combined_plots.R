
# R/utils-combined_plots.R

#' combined GMM barchart and count distribution map
#'
#' @param ich_mouse the ICH_Mouse

.combined_GMM_barchat_and_count_distribution_map <- function(ich_mouse) {

  on.exit(gc())

  GMM_barchart <- ICHMousewch:::.plotting_GMM_barchart(ich_mouse = ich_mouse)
  p1 <- ICHMousewch:::.create_count_distribution_map(seu_meta = ich_mouse@seu_metadata_with_cluster_symbol) %>%
    as.ggplot()
  p2 <- as.ggplot(GMM_barchart$`Barchart_Count-Original`) + theme(aspect.ratio = 1/1.414)
  p3 <- as.ggplot(GMM_barchart$`Barchart_Count-Normal`) + theme(aspect.ratio = 1/1.414)
  p4 <- as.ggplot(GMM_barchart$`Barchart_Count-Hematoma`) + theme(aspect.ratio = 1/1.414)

  cp1 <- (plot_spacer() | p1 | plot_spacer() | p2) + plot_layout(ncol = 4,
                                                                 nrow = 1,
                                                                 widths = c(2,97,1,100))

  cp2 <- (p3 | p4) + plot_layout(ncol = 2,
                                 nrow = 1,
                                 widths = c(1,1),
                                 heights = c(1,1))

  p <- (cp1 / cp2) +
    plot_layout(ncol = 1,
                nrow = 2,
                heights = c(1,1)) +
    plot_annotation(theme = theme(plot.background = element_rect(fill = "white",
                                                                 color = "black",
                                                                 linewidth = 1),
                                  plot.margin = margin(t = 5,
                                                       b = 5,
                                                       r = 5,
                                                       l = 5,
                                                       unit = "pt")))

  return(patchworkGrob(p))

}

#' combined volcano plot and venn diagram
#'
#' @param ich_mouse the ICH_Mouse

.combined_volcano_plot_and_venn_diagram <- function(ich_mouse) {

  on.exit(gc())

  volcano_plot <- ICHMousewch:::.plotting_volcano_plot(ich_mouse = ich_mouse)
  venn_diagram <- ICHMousewch:::.create_Venn_Diagram(ich_mouse = ich_mouse)

  p1 <- as.ggplot(venn_diagram$DEG_venn_plot)
  p2 <- as.ggplot(volcano_plot$edge_normal) + theme(aspect.ratio = 1/1.414)
  p3 <- as.ggplot(volcano_plot$center_normal) + theme(aspect.ratio = 1/1.414)
  p4 <- as.ggplot(volcano_plot$edge_center) + theme(aspect.ratio = 1/1.414)

  cp1 <- (p1 | p2) + plot_layout(ncol = 2,
                                 nrow = 1,
                                 widths = c(1,1),
                                 heights = c(1,1))

  cp2 <- (p3 | p4) + plot_layout(ncol = 2,
                                 nrow = 1,
                                 widths = c(1,1),
                                 heights = c(1,1))

  p <- (cp1 / cp2) + plot_layout(ncol = 1,
                                 nrow = 2,
                                 heights = c(1,1)) +
    plot_annotation(theme = theme(plot.background = element_rect(fill = "white",
                                                                 color = "black",
                                                                 linewidth = 1),
                                  plot.margin = margin(t = 5,
                                                       b = 5,
                                                       r = 5,
                                                       l = 5,
                                                       unit = "pt")))

  return(patchworkGrob(p))

}
