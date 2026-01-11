
#  R/utils-identify_hematoma.R

#' choose hematoma symbol
#'
#' @param hematoma the class of Hematoma

.choose_hematoma_symbol <- function(hematoma)
  {

  on.exit(gc())

  tissue_section <- image_read(hematoma@file_address["background_image_address"])
  tissue_section <- ggdraw() +
    draw_image(tissue_section) +
    ICHMousewch:::.plotting_theme() +
    theme(plot.margin = margin(t = 10,
                               b = 10,
                               r = 10,
                               l = 10,
                               unit = "pt"),
          aspect.ratio = 1)

  seu_meta <- hematoma@seu_metadata_with_cluster_symbol
  seu_meta[,aim_cluster := "0"]
  seu_meta[Louvain_cluster_posi == 1,aim_cluster := "Normal"]

  posi_cluster_num <- seu_meta[Louvain_cluster_posi != 1,Louvain_cluster_posi] %>%
    unique() %>%
    length()

  posi_cluster <- character()
  for (i in 1:posi_cluster_num) {

    seu_meta[Louvain_cluster_posi == i+1,aim_cluster := paste("region",(i+1))]

    posi_cluster <- c(posi_cluster,paste("region",(i+1)))

  }

  posi_cluster <- factor(posi_cluster,levels = posi_cluster)

  posi_cluster_image_ls <- vector("list",length = posi_cluster_num)
  names(posi_cluster_image_ls) <- posi_cluster
  for (i in 1:posi_cluster_num) {

    aim_posi_cluster <- posi_cluster[i]
    aim_posi_cluster_color <- "#D1352B"
    names(aim_posi_cluster_color) <- aim_posi_cluster

    posi_cluster_color_set <- rep("#3C77AF",times = (posi_cluster_num-1))
    names(posi_cluster_color_set) <- posi_cluster[!posi_cluster %in% aim_posi_cluster]

    hematoma@seu_metadata_with_cluster_symbol <- seu_meta
    spatial_image<- ICHMousewch:::.create_spatial_image_with_cluster_symbol(ich_mouse = hematoma,
                                                                            cluster_symbol = "aim_cluster",
                                                                            self_definition_color = c("Normal"="#F5D2A8",
                                                                                                      aim_posi_cluster_color,
                                                                                                      posi_cluster_color_set),
                                                                            show_legend_label = TRUE,
                                                                            show_image = FALSE) %>%
      as.ggplot()

    cbm_plot <- (spatial_image|tissue_section) + plot_layout(ncol = 2,
                                                             nrow = 1,
                                                             widths = c(1,1))

    posi_cluster_image_ls[posi_cluster[i]] <- patchworkGrob(cbm_plot) %>%
      list()

  }

  hematoma_symbol <- numeric()
  for (i in 1:posi_cluster_num) {

    ICHMousewch:::show_image(posi_cluster_image_ls[[posi_cluster[i]]])

    choice_title <- paste("choose Hematoma",as.character(posi_cluster[i]),sep = "-")

    choice <- menu(choices = c("Yes","No"),
                   title = choice_title,
                   graphics = FALSE)

    if (choice == 1) {

      region_number <- posi_cluster[i] %>%
        as.character() %>%
        strsplit(split = " ")
      region_number <- region_number[[1]][2] %>%
        as.numeric()

      hematoma_symbol <- c(hematoma_symbol,region_number)

    }
  }

  return(hematoma_symbol)

}

#' choose center edge symbol
#'
#' @param hematoma the class of Hematoma

.choose_center_edge_symbol <- function(hematoma)
  {

  on.exit(gc())

  tissue_section <- image_read(hematoma@file_address["background_image_address"])
  tissue_section <- ggdraw() +
    draw_image(tissue_section) +
    ICHMousewch:::.plotting_theme() +
    theme(plot.margin = margin(t = 10,
                               b = 10,
                               r = 10,
                               l = 10,
                               unit = "pt"),
          aspect.ratio = 1)

  seu_meta <- hematoma@seu_metadata_with_cluster_symbol
  seu_meta[,aim_cluster := "0"]
  seu_meta[Louvain_cluster_filt_gene == 1,aim_cluster := "Normal"]

  posi_cluster_num <- seu_meta[Louvain_cluster_filt_gene != 1,Louvain_cluster_filt_gene] %>%
    unique() %>%
    length()

  posi_cluster <- character()
  for (i in 1:posi_cluster_num) {

    seu_meta[Louvain_cluster_filt_gene == i+1,aim_cluster := paste("region",(i+1))]

    posi_cluster <- c(posi_cluster,paste("region",(i+1)))

  }

  posi_cluster <- factor(posi_cluster,levels = posi_cluster)

  posi_cluster_image_ls <- vector("list",length = posi_cluster_num)
  names(posi_cluster_image_ls) <- posi_cluster
  for (i in 1:posi_cluster_num) {

    aim_posi_cluster <- posi_cluster[i]
    aim_posi_cluster_color <- "#D1352B"
    names(aim_posi_cluster_color) <- aim_posi_cluster

    posi_cluster_color_set <- rep("#3C77AF",times = (posi_cluster_num-1))
    names(posi_cluster_color_set) <- posi_cluster[!posi_cluster %in% aim_posi_cluster]

    hematoma@seu_metadata_with_cluster_symbol <- seu_meta
    spatial_image<- ICHMousewch:::.create_spatial_image_with_cluster_symbol(ich_mouse = hematoma,
                                                                            cluster_symbol = "aim_cluster",
                                                                            self_definition_color = c("Normal"="#F5D2A8",
                                                                                                      aim_posi_cluster_color,
                                                                                                      posi_cluster_color_set),
                                                                            show_legend_label = TRUE,
                                                                            show_image = FALSE) %>%
      as.ggplot()

    cbm_plot <- (spatial_image|tissue_section) + plot_layout(ncol = 2,
                                                             nrow = 1,
                                                             widths = c(1,1))

    posi_cluster_image_ls[posi_cluster[i]] <- patchworkGrob(cbm_plot) %>%
      list()

  }

  center_edge_symbol <- numeric()
  for (i in 1:posi_cluster_num) {

    ICHMousewch:::show_image(posi_cluster_image_ls[[posi_cluster[i]]])

    choice_title <- paste("choose Hematoma Edge",as.character(posi_cluster[i]),sep = "-")

    choice <- menu(choices = c("Yes","No"),
                   title = choice_title,
                   graphics = FALSE)

    if (choice == 1) {

      region_number <- posi_cluster[i] %>%
        as.character() %>%
        strsplit(split = " ")
      region_number <- region_number[[1]][2] %>%
        as.numeric()

      center_edge_symbol <- c(center_edge_symbol,region_number)

    }
  }

  return(center_edge_symbol)

}

#' show choice plots
#'
#' @param graph the graph

.show_choice_plots <- function(graph)
  {

  on.exit(gc())

  dir_path <- getwd()

  ggsave(filename = ".choice_graph.png",
         path = dir_path,
         device = "png",
         plot = as.ggplot(graph),
         dpi = 600,
         width = 297,
         height = 210,
         unit = "mm")

  img <- image_read(path = paste(dir_path,".choice_graph.png",sep = "/"))

  print(img)

  unlink(dir_path,
         recursive = TRUE)

}
