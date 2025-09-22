
# R/utils-graph.R

#' create Giotto Instruction
#'
#' @param python_path the path to a python which the Giotto can use
#' @param results_folder the folder for Giotto save plots

.create_giotto_instruction <- function(python_path, results_folder) {

  on.exit(gc())

  instrs = createGiottoInstructions(save_dir = results_folder,
                                    save_plot = FALSE,
                                    show_plot = FALSE,
                                    python_path = python_path)

  return(instrs)

}

#' create spatial image with cluster symbol
#'
#' @param in_tissue_metadata the Seurat Object metadata of in tissue barcodes
#' @param cluster_symbol the cluster for spatial image
#' @param raw_count_matrix the matrix of raw count dataset
#' @param background_image_address image for background
#' @param color_set a set of color
#' @param self_definition_color a set of color defined by user
#' @param giotto_instruction the instruction of Giotto Object

.create_spatial_image_with_cluster_symbol <- function(in_tissue_metadata,cluster_symbol,raw_count_matrix,background_image_address,color_set,self_definition_color,giotto_instruction) {

  on.exit(gc())

  in_tissue_metadata[,cell_ID := barcode]

  in_tissue_count_matrix <- raw_count_matrix[,in_tissue_metadata[,barcode]]

  giotto_object <- createGiottoObject(expression = in_tissue_count_matrix,
                                      spatial_locs = in_tissue_metadata[,c("imagerow","imagecol")],
                                      cell_metadata = in_tissue_metadata,
                                      instructions = giotto_instruction)

  image_obj <- createGiottoImage(gobject = giotto_object,
                                 mg_object = background_image_address,
                                 name = "background_image",
                                 negative_y = FALSE)

  giotto_object <- addGiottoImage(gobject = giotto_object,images = list(image_obj))

  color_symbols <- in_tissue_metadata %>%
    subset(select = cluster_symbol) %>%
    unique() %>%
    unlist() %>%
    as.character()


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
                              point_alpha = 0.2,
                              cell_color_code = c(self_definition_color,random_colors),
                              background_color = "#00000000",
                              show_image = TRUE)

  return(spatial_image)

}

#' create spatial image with single gene
#'
#' @param seu_metadata_with_cluster_symbol the Seurat Object metadata of in tissue barcodes
#' @param gene_ls the gene list for creating spatial image
#' @param raw_count_matrix the matrix of raw count dataset
#' @param background_image_address image for background
#' @param giotto_instruction the instruction of Giotto Object
#' @param show_background_image whether to show background image

.create_spatial_image_with_single_gene <- function(seu_metadata_with_cluster_symbol,gene_ls,raw_count_matrix,background_image_address,giotto_instruction,show_background_image = TRUE) {

  on.exit(gc())

  seu_metadata_with_cluster_symbol[,cell_ID := barcode]

  in_tissue_count_matrix <- raw_count_matrix[gene_ls,seu_metadata_with_cluster_symbol[,barcode]]

  in_tissue_count_matrix[!in_tissue_count_matrix == 0] <- 1

  giotto_object <- createGiottoObject(expression = in_tissue_count_matrix,
                                      spatial_locs = seu_metadata_with_cluster_symbol[,c("imagerow","imagecol")],
                                      cell_metadata = seu_metadata_with_cluster_symbol,
                                      instructions = giotto_instruction)

  image_obj <- createGiottoImage(gobject = giotto_object,
                                 mg_object = background_image_address,
                                 name = "background_image",
                                 negative_y = FALSE)

  giotto_object <- addGiottoImage(gobject = giotto_object,images = list(image_obj))

  spatial_image_ls <- list()

  if(show_background_image) {

    if(length(gene_ls) > 0) {

      for (i in 1:length(gene_ls)) {

        spatial_image <- spatFeatPlot2D(gobject = giotto_object,
                                        expression_values = "raw",
                                        feats = gene_ls[i],
                                        background_color = "white",
                                        point_size = 0.5,
                                        point_alpha = 0.2,
                                        cell_color_gradient = c("#F5D2A8","#D1352B"),
                                        show_image = TRUE,
                                        show_legend = FALSE) %>%
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
                                        point_alpha = 1,
                                        cell_color_gradient = c("#F5D2A8","#D1352B"),
                                        show_image = FALSE,
                                        show_legend = FALSE) %>%
          ggplotGrob()

        spatial_image_ls <- append(spatial_image_ls,list(spatial_image))
        names(spatial_image_ls)[i] <- gene_ls[i]

      }

    }

  }


  return(spatial_image_ls)

}

