
# R/Spatial_Chunk.R

#' the class Spatial_Chunk
#'
#' @slot chunk_bin the bins of the chunk
#' @slot chunk_barcode the barcode of the chunk
#' @slot chunk_parameter the parameter of the chunk
#' @export

setClass(Class = "Spatial_Chunk",
         slots = c(chunk_set_name = "character",
                   chunk_bin = "list",
                   chunk_barcode = "list",
                   chunk_parameter = "list"))

#' initialize the class of Spatial_Chunk
#'
#' @param chunk_set_name the name of the chunk set
#' @param barcodes the barcodes for chunking
#' @param bin_scale the scale of the chunk

setMethod(f = "initialize",
          signature = signature(.Object = "Spatial_Chunk"),
          definition = function(.Object,barcodes,chunk_set_name,bin_scale = 5000) {

            on.exit(gc())

            .Object@chunk_set_name <- chunk_set_name

            .Object@chunk_parameter$total_number <- length(barcodes)
            .Object@chunk_parameter$bin_scale <- bin_scale
            .Object@chunk_parameter$bin_number <- .Object@chunk_parameter$total_number %/% .Object@chunk_parameter$bin_scale

            .Object@chunk_bin <- ICHMousewch:::.create_chunks(total_num = .Object@chunk_parameter$total_number,
                                                              bin_num = .Object@chunk_parameter$bin_number)

            .Object@chunk_barcode <- ICHMousewch:::.genereate_chunk_barcodes(barcodes = barcodes,
                                                                             chunk_bins = .Object@chunk_bin)

            validObject(.Object)
            return(.Object)

          })

#' initialize the class of Spatial_Chunk
#'
#' @param chunk_set_name the name of the chunk set
#' @param barcodes the barcodes for chunking
#' @param bin_scale the scale of the chunk
#' @export

Create_Spatial_Chunk <- function(chunk_set_name,barcodes,bin_scale = 5000) {

  on.exit(gc())

  spatial_chunk <- new(Class = "Spatial_Chunk",
                       barcodes = barcodes,
                       chunk_set_name = chunk_set_name,
                       bin_scale = bin_scale)

  return(spatial_chunk)

}
