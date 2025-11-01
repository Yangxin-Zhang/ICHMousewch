
# R/Gene_Spatial_Pattern.R

#' the Gene_Spatial_Pattern class
#'
#' @slot pattern_name the spatial pattern name
#' @slot gene_list the gene list for classification
#' @slot spatial_chunk the class of Spatial_Chunk
#' @slot characteristic_value_matrix the characteristic value matrix
#' @export

setClass(Class = "Gene_Spatial_Pattern",
         slots = c(pattern_name = "character",
                   gene_list = "character",
                   spatial_chunk = "Spatial_Chunk",
                   characteristic_value_matrix = "data.frame"))

#' initialize class of Gene_Spatial_Pattern
#'
#' @param gene_ls the genes for classfication
#' @param ich_mouse the class of ich_mouse
#' @param pattern_na the spatial pattern name
#' @param bin_scale the scale of the bin

setMethod(f = "initialize",
          signature = signature(.Object = "Gene_Spatial_Pattern"),
          definition = function(.Object,ich_mouse,gene_ls,pattern_na,bin_scale = 5000) {

            on.exit(gc())

            .Object@pattern_name <- pattern_na

            .Object@gene_list <- gene_ls

            .Object@spatial_chunk <- ICHMousewch::Create_Spatial_Chunk(chunk_set_name = .Object@pattern_name,
                                                                       barcodes = ich_mouse@seu_metadata_with_cluster_symbol[,cell_ID],
                                                                       bin_scale = bin_scale)

            .Object@characteristic_value_matrix <- purrr::map(.x = .Object@spatial_chunk@spatial_chunk@chunk_barcode,
                                                              .f = ~ICHMousewch:::.calculate_chunk_characteristic_value(seu_meta = ich_mouse@seu_metadata_with_cluster_symbol,
                                                                                                                        raw_count_matrix = ich_mouse@raw_count_matrix,
                                                                                                                        spatial_chunks = .x,
                                                                                                                        gene_ls = .Object@gene_list)) %>%
              bind_rows()


            validObject(.Object)
            return(.Object)

          })

#' create class of Gene_Spatial_Pattern
#'
#' @param gene_ls the genes for classfication
#' @param ich_mouse the class of ich_mouse
#' @param pattern_na the spatial pattern name
#' @param bin_scale the scale of the bin
#' @export

Create_Gene_Spatial_Pattern <- function(ich_mouse,gene_ls,pattern_na,bin_scale = 5000) {

  on.exit(gc())

  gene_spatial_pattern <- new(Class = "Gene_Spatial_Pattern",
                              ich_mouse = ich_mouse,
                              gene_ls = gene_ls,
                              pattern_na = pattern_na,
                              bin_scale = bin_scale)

  return(gene_spatial_pattern)

}
