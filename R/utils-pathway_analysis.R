
# R/utils-pathway_analysis.R

#' calculate progeny score
#'
#' @param ich_mouse the class of ICH_Mouse

.calculate_progeny_score <- function(ich_mouse) {

  on.exit(gc())

  raw_count_matrix <- ich_mouse@raw_count_matrix
  filtered_barcode <- ich_mouse@filtered_genes
  filtered_genes <- ich_mouse@filtered_barcodes

  filtered_matrix <- raw_count_matrix[filtered_barcode,filtered_genes]

  seu_obj <- CreateSeuratObject(counts = filtered_matrix) %>%
    NormalizeData(normalization.method = "RC",
                  scale.factor = 1e6)

  normalized_matrix <- GetAssayData(seu_obj,
                                    layer = "data")

  chunks <- Create_Spatial_Chunk(chunk_set_name = "progeny",
                                 barcodes = ich_mouse@filtered_genes,
                                 bin_scale = 2000)
  chunks <- chunks@chunk_barcode
  chunk_name <- names(chunks)

  progeney_result_ls <- vector("list",length = length(chunks))
  names(progeney_result_ls) <- chunk_name
  for (i in length(chunk_name)) {

    progeney_chunk <- chunks[[chunk_name[i]]]

    progeney_matrix <- normalized_matrix[progeney_chunk,] %>%
      as.matrix()

    progeney_result <- progeny::progeny(expr = progeney_matrix,
                                        organism = "Mouse") %>%
      as.data.frame() %>%
      rownames_to_column()

    progeney_result_ls[chunk_name[i]] <- list(progeney_result)

  }
}
