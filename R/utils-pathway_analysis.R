
# R/utils-pathway_analysis.R

#' calculate progeny score
#'
#' @param ich_mouse the class of ICH_Mouse

.calculate_progeny_score <- function(ich_mouse) {

  on.exit(gc())

  raw_count_matrix <- ich_mouse@raw_count_matrix
  filtered_barcode <- ich_mouse@filtered_barcodes
  filtered_genes <- ich_mouse@filtered_genes

  filtered_matrix <- raw_count_matrix[filtered_genes,filtered_barcode]

  seu_obj <- CreateSeuratObject(counts = filtered_matrix) %>%
    NormalizeData(normalization.method = "RC",
                  scale.factor = 1e6)

  barcode_chunks <- Create_Spatial_Chunk(chunk_set_name = "barcode_chunk",
                                         barcodes = filtered_barcode,
                                         bin_scale = 5000)
  barcode_chunks <- barcode_chunks@chunk_barcode
  barcode_chunk_name <- names(barcode_chunks)

  gene_chunks <- Create_Spatial_Chunk(chunk_set_name = "gene_chunk",
                                 barcodes = filtered_genes,
                                 bin_scale = 2000)
  gene_chunks <- gene_chunks@chunk_barcode

  calculate_seperate_progeney_result <- function(barcode_chunk,gene_chunks) {

    on.exit(gc())

    gene_chunk_name <- names(gene_chunks)

    scaled_matrix_ls <- vector("list",length = length(gene_chunks))
    names(scaled_matrix_ls) <- gene_chunk_name
    for (i in 1:length(gene_chunk_name)) {

      scaled_chunk <- gene_chunks[[gene_chunk_name[i]]]

      scaled_matrix <- seu_obj[scaled_chunk,] %>%
        ScaleData() %>%
        GetAssayData(layer = "scale.data")

      scaled_matrix_ls[gene_chunk_name[i]] <- list(scaled_matrix[,barcode_chunk])

      gc()

    }

    combined_matrix <- ICHMousewch::combined_matrix_on_column(matrix_ls = scaled_matrix_ls)

    progeney_result <- progeny(expr = combined_matrix,
                               organism = "Mouse") %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      rename(barcode = rowname) %>%
      as.data.table()

    return(progeney_result)

  }

  progeney_result_ls <- lapply(barcode_chunks,calculate_seperate_progeney_result,gene_chunks = gene_chunks)

  progeney_result <- bind_rows(progeney_result_ls) %>%
    as_tibble() %>%
    column_to_rownames(var = "barcode")
  progeney_result[progeney_result < 0] <- 0
  progeney_result <- progeney_result %>%
    rownames_to_column() %>%
    rename(barcode = rowname) %>%
    as.data.table()

  coord_df <- merge(ich_mouse@seu_metadata_with_cluster_symbol[barcode %in% progeney_result$barcode,],progeney_result,by = "barcode")

  return(coord_df)

}

#' create progeny pathway score map
#'
#' @param progeny_score_df the progeny score dataframe

.create_progeny_pathway_score_map <- function(progeny_score_df) {

  on.exit(gc())

  pathway_symbol_ls <- c("WNT","TNFa","Trail","VEGF",
                         "p53","PI3K","TGFb","JAK-STAT",
                         "MAPK","NFkB","EGFR","Estrogen",
                         "Hypoxia","Androgen")

  create_pathway_graph <- function(pathway_symbol,progeny_score_df) {

    on.exit(gc())

    plotting_col <- progeny_score_df[,..pathway_symbol]
    progeny_score_df[,plotting_pathway := plotting_col]
    progeny_score_df[,point_symbol := "point_symbol"]
    progeny_score_df[center_edge_symbol %in% c("2","3") & plotting_pathway == 0,point_symbol := "0"]
    progeny_score_df[plotting_pathway == 0 & !center_edge_symbol %in% c("2","3"),point_symbol := "1"]
    progeny_score_df[plotting_pathway != 0,point_symbol := "2"]
    coord_df <- progeny_score_df[,c("plotting_pathway","row","col","point_symbol")]

    pathway_graph <- ggplot() +
      geom_point(data = coord_df[point_symbol == "0"],
                 mapping = aes(x = -row,y = -col),
                 color = "white",
                 size = 0.01) +
      geom_point(data = coord_df[point_symbol == "1"],
                 mapping = aes(x = -row,y = -col),
                 color = "grey90",
                 size = 0.01) +
      geom_point(data = coord_df[point_symbol == "2"],
                 mapping = aes(x = -row,y = -col,colour = plotting_pathway),
                 size = 0.01) +
      scale_color_gradientn(colors = c("#FEF4E8", "#FED9A6", "#FEB24C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026")) +
      coord_fixed() +
      labs(title = pathway_symbol) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "right",
            legend.title = element_blank(),
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

  }

  pathway_graph_ls <- vector("list",length = length(pathway_symbol_ls))
  names(pathway_graph_ls) <- pathway_symbol_ls
  for (i in 1:length(pathway_symbol_ls)) {

    pathway_graph_ls[pathway_symbol_ls[i]] <- pathway_symbol_ls[i] %>%
      create_pathway_graph(progeny_score_df = progeny_score_df) %>%
      ggplotGrob() %>%
      list()

  }

  return(pathway_graph_ls)

}
