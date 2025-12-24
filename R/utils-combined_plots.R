
# R/utils-combined_plots.R

#' combined GMM barchart and count distribution map
#'
#' @param ich_mouse the ICH_Mouse

.combined_GMM_barchat_and_count_distribution_map <- function(ich_mouse)
  {

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

#' combined GO overview plot
#'
#' @param ich_mouse the ICH_Mouse
#' @param group_symbol the group symbol

.combined_GO_overview_plot <- function(ich_mouse,group_symbol) {

  on.exit(gc())

  go_overview_barchart <- ICHMousewch:::.create_GO_result_overview_barchart(ich_mouse = ich_mouse)
  bubble_chart <- ICHMousewch:::.plotting_bubble_chart(ich_mouse = ich_mouse,group_na = group_symbol)
  similarity_heatmap <- ICHMousewch:::.plotting_GO_term_similarity_heatmap(ich_mouse = ich_mouse,
                                                                           cluster_symbol = "total-diff_expr_genes",
                                                                           special_block = TRUE,
                                                                           GO_group_symbol = group_symbol)

  p1 <- as.ggplot(go_overview_barchart$overview_GO_result_barchart)
  p2 <- as.ggplot(similarity_heatmap$GO_heatmap)
  p3 <- as.ggplot(bubble_chart$bubble_chart)

  cp <- (plot_spacer()|p1|p2|plot_spacer()) + plot_layout(ncol = 4,
                                                          nrow = 1,
                                                          widths = c(1,220,200,1))
  p <- (cp/p3) + plot_layout(ncol = 1,
                             nrow = 2,
                             heights = c(3,5)) +
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

#' combine spatial image and umap
#'
#' @param ich_mouse the ICH_Mouse

.combine_spatial_image_and_umap <- function(ich_mouse) {

  on.exit(gc())

  original_hematoma <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(cluster_symbol = "hematoma_symbol",
                                                                               self_definition_color = c("1"="#F5D2A8","2"="#D1352B"),
                                                                               plot_title = "Hematoma",
                                                                               show_plot_title = TRUE,
                                                                               show_legend_label = TRUE,
                                                                               legend_lable = c("1" = "Normal","2" = "Hematoma"),
                                                                               ich_mouse = ich_mouse)

  hematoma_center_edge <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(cluster_symbol = "center_edge_symbol",
                                                                                  self_definition_color = c("1"="#F5D2A8","2"="#D1352B","3"="#3C77AF"),
                                                                                  show_plot_title = TRUE,
                                                                                  show_legend_label = TRUE,
                                                                                  legend_lable = c("1" = "Normal","2" = "Hematoma","3" = "Edge"),
                                                                                  ich_mouse = ich_mouse,
                                                                                  plot_title = "Hematoma Center vs Hematoma Edge")

  umap_plots <- ICHMousewch:::.create_umap_plot(ich_mouse = ich_mouse,
                                                spaceranger_result = list("spacerange_umap_address" = "D:/Pango_Project/ICH_Mouse_Analysis/Original_Data/projection.csv",
                                                                          "spacerange_cluster_address" = "D:/Pango_Project/ICH_Mouse_Analysis/Original_Data/clusters.csv"))

  p1 <- as.ggplot(original_hematoma)
  p2 <- as.ggplot(hematoma_center_edge)
  p3 <- as.ggplot(umap_plots$umap_spaceranger_cluster_spaceranger)
  p4 <- as.ggplot(umap_plots$umap_hematoma_spaceranger)

  cp1 <- (p1|p2) + plot_layout(ncol = 2,
                               nrow = 1,
                               width = c(1,1))
  cp2 <- (p3|p4) + plot_layout(ncol = 2,
                               nrow = 1,
                               width = c(1,1))
  p <- (cp1/cp2) + plot_layout(ncol = 1,
                               nrow = 2,
                               heights = c(5,4)) +
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
