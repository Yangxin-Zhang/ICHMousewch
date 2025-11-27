
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

}
