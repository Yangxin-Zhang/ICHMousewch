
# utils-graph_spatial_image.R

#' create spatial image with cluster symbol
#'
#' @param ich_mouse the ICH_Mouse class
#' @param cluster_symbol the cluster for spatial image
#' @param self_definition_color a set of color defined by user
#' @param plot_title the title of spatial image
#' @param legend_lable the legend label

.create_spatial_image_with_cluster_symbol <- function(ich_mouse,
                                                      cluster_symbol,
                                                      self_definition_color = character(),
                                                      theme_param = list(),
                                                      plot_title = character(),
                                                      legend_lable = character(),
                                                      show_image = TRUE,
                                                      show_legend_label = FALSE,
                                                      show_plot_title = FALSE)
  {

  on.exit(gc())

  in_tissue_metadata <- ich_mouse@seu_metadata_with_cluster_symbol
  in_tissue_metadata[,cell_ID := barcode]

  raw_count_matrix <- ich_mouse@raw_count_matrix

  giotto_instruction <- ich_mouse@giotto_instruction[[1]]

  background_image_address <- ich_mouse@file_address["background_image_address"]

  if (show_legend_label == FALSE) {

    theme_param <- c(theme_param,list(legend.position = "none"))

  } else {

    if (length(legend_lable) != 0) {

      ori_lable <- names(legend_lable)

      for (i in 1:length(ori_lable)) {

        filter_condition <- unlist(in_tissue_metadata[,..cluster_symbol]) %in% ori_lable[i]

        in_tissue_metadata[filter_condition,new_cluster_symbol := legend_lable[ori_lable[i]]]

      }

      in_tissue_metadata <- ICHMousewch::add_new_col_to_data_table(original_data_table = in_tissue_metadata,
                                                                   new_col = in_tissue_metadata[,new_cluster_symbol],
                                                                   new_col_name = cluster_symbol)

      if (length(self_definition_color) != 0) {

        new_self_def_col <- vector("character",length = length(ori_lable))
        names(new_self_def_col) <- legend_lable

        for (i in 1:length(ori_lable)) {

          new_self_def_col[legend_lable[ori_lable[i]]] <- self_definition_color[ori_lable[i]]

        }

        self_definition_color <- new_self_def_col

      }

    }

    theme_param <- c(theme_param,list(legend.position = "right"))

  }

  in_tissue_count_matrix <- raw_count_matrix[,in_tissue_metadata[,barcode]]

  image_obj <- createGiottoImage(mg_object = background_image_address,
                                 name = "background_image",
                                 negative_y = FALSE)

  suppressWarnings(
    giotto_object <- createGiottoObject(expression = in_tissue_count_matrix,
                                        spatial_locs = in_tissue_metadata[,c("imagerow","neg_imagecol")],
                                        cell_metadata = in_tissue_metadata,
                                        instructions = giotto_instruction,
                                        images = list(image_obj))
  )

  color_symbols <- in_tissue_metadata %>%
    subset(select = cluster_symbol) %>%
    unique() %>%
    unlist() %>%
    as.character()

  color_set <- viridis::plasma(length(color_symbols))

  random_color_symbols <- color_symbols[!color_symbols %in% names(self_definition_color)]

  sub_color_set <- color_set[!color_set %in% self_definition_color]

  if(length(random_color_symbols) > 0) {

    select_color <- sample(1:length(sub_color_set),length(random_color_symbols))
    random_colors <- sub_color_set[select_color]
    names(random_colors) <- random_color_symbols

  } else {

    random_colors <- character()

  }

  spatial_image <- spatPlot2D(gobject = giotto_object,
                              cell_color = cluster_symbol,
                              point_size = 0.5,
                              point_shape = "border",
                              point_alpha = 0.5,
                              point_border_stroke = 0,
                              cell_color_code = c(self_definition_color,random_colors),
                              background_color = "#00000000",
                              show_image = show_image,
                              show_legend = show_legend_label,
                              title = plot_title,
                              coord_fix_ratio = 1,
                              return_plot = TRUE,
                              legend_symbol_size = 5,
                              theme_param = c(theme_param,list(ICHMousewch:::.plotting_theme(),
                                                               axis.ticks.x = element_blank(),
                                                               axis.ticks.y = element_blank(),
                                                               axis.text.x = element_blank(),
                                                               axis.text.y = element_blank(),
                                                               axis.title.x = element_blank(),
                                                               axis.title.y = element_blank())))

  return(ggplotGrob(spatial_image))

}

#' create spatial image with single gene
#'
#' @param seu_metadata_with_cluster_symbol the Seurat Object metadata of in tissue barcodes
#' @param gene_ls the gene list for creating spatial image
#' @param raw_count_matrix the matrix of raw count dataset
#' @param background_image_address image for background
#' @param giotto_instruction the instruction of Giotto Object
#' @param show_background_image whether to show background image

.create_spatial_image_with_single_gene <- function(seu_metadata_with_cluster_symbol,
                                                   gene_ls,
                                                   raw_count_matrix,
                                                   background_image_address,
                                                   giotto_instruction,
                                                   show_background_image = TRUE,
                                                   gradient = FALSE)
  {

  on.exit(gc())

  in_tissue_metadata <- seu_metadata_with_cluster_symbol[,cell_ID := barcode]

  in_tissue_count_matrix <- raw_count_matrix[gene_ls,seu_metadata_with_cluster_symbol[,barcode],drop = FALSE]

  if (gradient) {

    edge_barcodes <- in_tissue_metadata[!center_edge_symbol == "1",barcode]
    normal_barcodes <- in_tissue_metadata[center_edge_symbol == "1",barcode]

    logi_zero <- !in_tissue_count_matrix == 0
    logi_edge <-  Matrix::Matrix(rep(colnames(in_tissue_count_matrix) %in% edge_barcodes,times = length(gene_ls)),
                                 nrow = length(gene_ls),
                                 byrow = TRUE,
                                 sparse = TRUE)
    logi_normal <- Matrix::Matrix(rep(colnames(in_tissue_count_matrix) %in% normal_barcodes,times = length(gene_ls)),
                                  nrow = length(gene_ls),
                                  byrow = TRUE,
                                  sparse = TRUE)

    in_tissue_count_matrix[logi_zero & logi_edge] <- 1
    in_tissue_count_matrix[logi_zero & logi_normal] <- 2

  } else {

    in_tissue_count_matrix[!in_tissue_count_matrix == 0] <- 1

  }

  image_obj <- createGiottoImage(mg_object = background_image_address,
                                 name = "background_image",
                                 negative_y = FALSE)

  giotto_object <- createGiottoObject(expression = in_tissue_count_matrix,
                                      spatial_locs = in_tissue_metadata[,c("imagerow","neg_imagecol")],
                                      cell_metadata = in_tissue_metadata,
                                      instructions = giotto_instruction,
                                      images = list(image_obj))

  spatial_image_ls <- list()

  if(gradient) {

    if(length(gene_ls) > 0) {

      for (i in 1:length(gene_ls)) {

        spatial_image <- spatFeatPlot2D(gobject = giotto_object,
                                        expression_values = "raw",
                                        feats = gene_ls[i],
                                        background_color = "white",
                                        point_size = 0.5,
                                        point_alpha = 0.5,
                                        cell_color_gradient = c("#F5D2A8","#D1352B","#3C77AF"),
                                        gradient_midpoint = 1,
                                        gradient_style = "divergent",
                                        show_image = show_background_image,
                                        theme_param = list(plot.background = element_rect(fill = "white",color = NA),
                                                           axis.ticks.x = element_blank(),
                                                           axis.ticks.y = element_blank(),
                                                           axis.text.x = element_blank(),
                                                           axis.text.y = element_blank(),
                                                           axis.title.x = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           legend.position = "none")) %>%
          ggplotGrob()

        spatial_image_ls <- append(spatial_image_ls,list(spatial_image))
        names(spatial_image_ls)[i] <- gene_ls[i]

      }

    }

  } else {

    if(length(gene_ls) > 0) {

      for (i in 1:length(gene_ls)) {

        spatial_image <- spatFeatPlot2D(gobject = giotto_object,
                                        expression_values = "raw",
                                        feats = gene_ls[i],
                                        background_color = "white",
                                        point_size = 0.5,
                                        point_alpha = 0.5,
                                        cell_color_gradient = c("#F5D2A8","#D1352B"),
                                        show_image = show_background_image,
                                        theme_param = list(plot.background = element_rect(fill = "white",color = NA),
                                                           axis.ticks.x = element_blank(),
                                                           axis.ticks.y = element_blank(),
                                                           axis.text.x = element_blank(),
                                                           axis.text.y = element_blank(),
                                                           axis.title.x = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           legend.position = "none")) %>%
          ggplotGrob()

        spatial_image_ls <- append(spatial_image_ls,list(spatial_image))
        names(spatial_image_ls)[i] <- gene_ls[i]

      }

    }

  }


  return(spatial_image_ls)

}

#' create gene distribution map
#'
#' @param ich_mouse the class of ICH_Mouse
#' @param diff_expr_gene_symbol the differential expression genes symbol
#' @param background_genes the background genes for graph
#' @param aim_gene the aim gene for plotting

.create_gene_distribution_map <- function(ich_mouse,
                                          diff_expr_gene_symbol = "edge-normal",
                                          aim_gene = character(),
                                          background_genes = c("Hbb-bt","Hbb-bs","Hba-a2"),
                                          point_size = 0.0001,
                                          normalized_data = FALSE,
                                          scaled_data = TRUE)
  {

  on.exit(gc())

  seu_meta <- ich_mouse@seu_metadata_with_cluster_symbol
  raw_count_matrix <- ich_mouse@raw_count_matrix
  filtered_genes <- ich_mouse@filtered_genes
  filtered_barcodes <- ich_mouse@filtered_barcodes
  diff_expr_gene <- ich_mouse@diff_expr_genes[[diff_expr_gene_symbol]]

  if (scaled_data == normalized_data) {

    scaled_data <- TRUE
    normalized_data <- FALSE

  }

  seu_meta <- seu_meta[,cell_ID := barcode]
  seu_meta <- seu_meta[barcode %in% filtered_barcodes]

  aim_gene <- aim_gene[aim_gene %in% filtered_genes]

  filtered_count_matrix <- raw_count_matrix[filtered_genes,seu_meta[,cell_ID]]

  seu_meta[,plot_row := -row]
  seu_meta[,plot_col := -col]

  if (scaled_data) {

    seu_obj <- CreateSeuratObject(counts = filtered_count_matrix) %>%
      NormalizeData() %>%
      ScaleData(features = diff_expr_gene[avg_log2FC > 1,gene_name])

    gene_expr_matrix <- t(seu_obj@assays$RNA$scale.data) %>%
      as.data.frame()

    color_legend_na <- "Scaled\nExpression"

  }

  if (normalized_data) {

    seu_obj <- CreateSeuratObject(counts = filtered_count_matrix) %>%
      NormalizeData()

    sub_seu_obj <- seu_obj[diff_expr_gene[avg_log2FC > 1,gene_name],]

    gene_expr_matrix <- sub_seu_obj@assays$RNA$data %>%
      as.matrix() %>%
      t() %>%
      as.data.frame()

    color_legend_na <- "Normalized\nExpression"

  }

  gene_expr_matrix <- bind_cols(data.frame(barcode = rownames(gene_expr_matrix)),
                                gene_expr_matrix) %>%
    as.data.table()

  graph_df <- merge(seu_meta,gene_expr_matrix,by = "barcode")

  initialize_color_symbol <- function(background_genes,filtered_count_matrix,graph_df) {

    graph_df[,color_symbol := "0"]

    if (length(background_genes) == 1) {

      background_df <- data.frame(background_genes_name = filtered_count_matrix[background_genes,])

    } else {

      background_df <- filtered_count_matrix[background_genes,] %>%
        as.data.frame() %>%
        t()

    }

    graph_df[,background_count := Matrix::rowSums(background_df != 0)]
    graph_df[background_count != 0,color_symbol := "1"]

    return(graph_df)

  }


  if (length(aim_gene) != 0) {

    distribution_graph_ls <- vector("list",length = length(aim_gene))
    names(distribution_graph_ls) <- aim_gene

    for (i in 1:length(aim_gene)) {

      graph_df <- initialize_color_symbol(background_genes = background_genes,
                                          filtered_count_matrix = filtered_count_matrix,
                                          graph_df = graph_df)

      aim_gene_na <- aim_gene[i]

      aim_gene_df <- data.frame(aim_gene_name = filtered_count_matrix[aim_gene_na,])

      aim_gene_scale_data <- graph_df[,..aim_gene_na]
      graph_df[,plotting_gene := aim_gene_scale_data]

      graph_df[,aim_gene_count := Matrix::rowSums(aim_gene_df != 0)]

      graph_df[aim_gene_count != 0,color_symbol := "2"]

      distribution_graph <- ggplot() +
        geom_point(data = graph_df[color_symbol == "0"],
                   mapping = aes(x = plot_row,y = plot_col,colour = color_symbol),
                   size = point_size,
                   color = "gray90") +
        geom_point(data = graph_df[color_symbol == "1"],
                   mapping = aes(x = plot_row,y = plot_col,colour = color_symbol),
                   size = point_size,
                   color = "white") +
        geom_point(data = graph_df[color_symbol == "2"],
                   mapping = aes(x = plot_row,y = plot_col,colour = plotting_gene),
                   size = point_size) +
        scale_color_gradientn(colors = c("#FEF4E8", "#FED9A6", "#FEB24C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026"),
                              values = scales::rescale(seq(from = min(graph_df[color_symbol == "2",plotting_gene]),
                                                           to = max(graph_df[color_symbol == "2",plotting_gene]),
                                                           length.out = 7),
                                                       to = c(0,1))) +
        coord_fixed() +
        labs(title = aim_gene_na,
             color = color_legend_na) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = "right",
              legend.title = element_text(hjust = 0.5),
              legend.text = element_text(size = 8,
                                         family = "Arial"),
              panel.background = element_rect(fill = "white", color = "black"),
              plot.background = element_rect(fill = "white", color = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 12,
                                        face = "bold",
                                        family = "Arial",
                                        hjust = 0.5,
                                        vjust = 0,
                                        margin = margin(b = 10)))

      distribution_graph_ls[aim_gene_na] <- list(ggplotGrob(distribution_graph))

    }

  } else {

    graph_df <- initialize_color_symbol(background_genes = background_genes,
                                        filtered_count_matrix = filtered_count_matrix,
                                        graph_df = graph_df)

    distribution_graph <- ggplot() +
      geom_point(data = graph_df[color_symbol == "0"],
                 mapping = aes(x = plot_row,y = plot_col,colour = color_symbol),
                 size = 0.01) +
      scale_colour_manual(values = c("0" = "grey")) +
      geom_point(data = graph_df[color_symbol == "1"],
                 mapping = aes(x = plot_row,y = plot_col,colour = color_symbol),
                 size = 0.01,
                 color = "white") +
      labs(title = "Background") +
      coord_fixed() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "right",
            legend.title = element_blank(),
            legend.text = element_text(size = 8,
                                       family = "Arial",
                                       color = "white"),
            legend.key = element_blank(),
            legend.background = element_blank(),
            legend.box.background = element_blank(),
            panel.background = element_rect(fill = "white", color = "black"),
            plot.background = element_rect(fill = "white", color = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 12,
                                      face = "bold",
                                      family = "Arial",
                                      hjust = 0.5,
                                      vjust = 0,
                                      margin = margin(b = 10))) +
      guides(color = guide_legend(override.aes = list(colour = "white")))

    background_na <- background_genes %>%
      paste(collapse = "-") %>%
      paste("background",sep = "_")
    distribution_graph_ls <- list()
    distribution_graph_ls[background_na] <- list(ggplotGrob(distribution_graph))

  }

  return(distribution_graph_ls)

}

#' create count distribution map
#'
#' @param seu_meta the Seurat metadata

.create_count_distribution_map <- function(seu_meta,
                                           point_size = 0.01)
  {

  on.exit(gc())

  seu_meta[,plot_row := -row]
  seu_meta[,plot_col := -col]

  count_distribution_map <- ggplot() +
    geom_point(data = seu_meta,
               mapping = aes(x = plot_row,y = plot_col,color = nCount_log2),
               size = point_size) +
    scale_color_gradientn(colors = c("#0D0887", "#3F007D", "#6A00A8", "#B12A90", "#E16462", "#FCA636", "#F0F921"),
                          limits = range(seu_meta[,nCount_log2])) +
    labs(color = "counts") +
    coord_fixed() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "left",
          legend.text = element_text(size = 8,
                                     family = "Arial"),
          legend.title = element_text(size = 12,
                                      family = "Arial"),
          panel.background = element_rect(fill = "white", color = "black"),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_blank())

  return(count_distribution_map)

}

#' create percent distribution map
#'
#' @param ich_mouse the class of ICH_Mouse
#' @param aim_gene the aim gene list
#' @param distribution_region the region for distribution

.create_percent_distribution_map <- function(ich_mouse,
                                             aim_gene,
                                             diff_expr_gene_symbol = "edge-normal",
                                             distribution_region  = "center_edge_symbol")
  {

  on.exit(gc())

  seu_meta <- ich_mouse@seu_metadata_with_cluster_symbol
  raw_count_matrix <- ich_mouse@raw_count_matrix
  filtered_genes <- ich_mouse@filtered_genes
  filtered_barcodes <- ich_mouse@filtered_barcodes
  diff_expr_gene <- ich_mouse@diff_expr_genes[[diff_expr_gene_symbol]]

  seu_meta[,cell_ID := barcode]
  seu_meta <- seu_meta[cell_ID %in% filtered_barcodes]

  aim_gene <- aim_gene[aim_gene %in% filtered_genes]

  filtered_count_matrix <- raw_count_matrix[filtered_genes,seu_meta[,cell_ID]]

  gene_expr_matrix <- filtered_count_matrix[diff_expr_gene[avg_log2FC > 1,gene_name],] %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()

  gene_expr_matrix <- bind_cols(data.frame(barcode = rownames(gene_expr_matrix)),
                                gene_expr_matrix) %>%
    as.data.table()

  graph_df <- merge(seu_meta,gene_expr_matrix,by = "barcode")

  percent_distribution_map <- function(graph_df,aim_gene,background_symbol,distribution_region = "center_edge_symbol") {

    on.exit(gc())

    if (distribution_region == "center_edge_symbol") {

      if (background_symbol == 1) {

        graph_title <- "Normal"

      }

      if (background_symbol == 2) {

        graph_title <- "Center"

      }

      if (background_symbol == 3) {

        graph_title <- "Edge"

      }

    }

    graph_df[,color_symbol := "0"]

    graph_df[,distribution_logi := graph_df[,..distribution_region]]
    graph_df[distribution_logi == background_symbol,color_symbol := "1"]
    graph_df[,aim_gene_logi := graph_df[,..aim_gene]]
    graph_df[aim_gene_logi != 0 & distribution_logi == background_symbol,color_symbol := "2"]

    per_dist_map <- ggplot() +
      geom_point(data = graph_df[color_symbol == "0"],
                 mapping = aes(x = -row,y = -col,colour = color_symbol),
                 size = 0.01) +
      scale_colour_manual(values = "gray90",
                          labels = "Other") +
      guides(color = guide_none()) +
      new_scale_colour() +
      geom_point(data = graph_df[color_symbol == "1"],
                 mapping = aes(x = -row,y = -col,colour = color_symbol),
                 size = 0.01) +
      scale_colour_manual(values = "#92C5DE",
                          labels = "Background") +
      guides(color = guide_legend(theme = theme(legend.title = element_blank()),
                                  override.aes = list(size = 5))) +
      new_scale_colour() +
      geom_point(data = graph_df[color_symbol == "2"],
                 mapping = aes(x = -row,y = -col,colour = color_symbol),
                 size = 0.01) +
      scale_colour_manual(values = "#E31A1C",
                          labels = "Positive") +
      guides(color = guide_legend(theme = theme(legend.title = element_blank()),
                                  override.aes = list(size = 5)))  +
      coord_fixed() +
      labs(title = paste(aim_gene,graph_title,sep = "-")) +
      ICHMousewch:::.plotting_theme() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank())

    return(ggplotGrob(per_dist_map))

  }

  region_symbol <- graph_df[,..distribution_region] %>%
    unique() %>%
    unlist()
  names(region_symbol) <- NULL
  region_num <- length(region_symbol)
  graph_ls <- list()
  for (i in 1:length(aim_gene)) {

    for (j in 1:region_num) {

      graph_na <- paste(distribution_region,aim_gene[i],sep = "-") %>%
        paste(region_symbol[j],sep = "-")

      graph_ls[graph_na] <- graph_df %>%
        percent_distribution_map(aim_gene = aim_gene[i],
                                 background_symbol = region_symbol[j],
                                 distribution_region = distribution_region) %>%
        list()

    }

  }

  return(graph_ls)

}
