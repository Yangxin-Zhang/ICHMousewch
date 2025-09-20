
# R/class-ICH_Mouse.R

#' ICH_Mouse a class for analysis ICH_Mouse
#'
#' @slot diff_expr_genes the differential expression genes
#' @slot symbol_genes the symbol genes

setClass(Class = "ICH_Mouse",
         slots = c(diff_expr_genes = "list",
                   symbol_genes = "list"),
         contains = "Hematoma")

#' Initialize class ICH_Mouse
#'

setMethod(f = "initialize",
          signature = signature(.Object = "ICH_Mouse"),
          definition = function(.Object,analysis_symbol,raw_count_matrix_address,filtered_count_matrix_address,tissue_position_address,background_image_address,giotto_python_path,giotto_results_folder) {

            on.exit(gc())

            .Object <- callNextMethod(.Object,analysis_symbol,raw_count_matrix_address,filtered_count_matrix_address,tissue_position_address,background_image_address,giotto_python_path,giotto_results_folder)

            .Object@symbol_genes$normal_tissue <- ICHMousewch:::.find_symbol_genes(raw_count_matrix = .Object@raw_count_matrix,
                                                                     seu_metadata_with_cluster_symbol = .Object@seu_metadata_with_cluster_symbol,
                                                                     filtered_genes = .Object@filtered_genes,
                                                                     cluster_symbol = 1)

            validObject(.Object)
            return(.Object)

          })

#' ICH_Mouse constructor
#'
#' @param analysis_symbol a symbol of sample
#' @param raw_count_matrix_address address for raw count matrix
#' @param filtered_count_matrix address for filtered count matrix
#' @param tissue_position_address address for tissue position matrix
#' @param background_image_address address for background address
#' @export

Create_ICH_Mouse <- function(analysis_symbol,raw_count_matrix_address,filtered_count_matrix_address,tissue_position_address,background_image_address,giotto_python_path,giotto_results_folder) {

  on.exit(gc())

  ICH_Mouse <- new(Class = "ICH_Mouse",
                   analysis_symbol = analysis_symbol,
                   raw_count_matrix_address = raw_count_matrix_address,
                   filtered_count_matrix_address = filtered_count_matrix_address,
                   tissue_position_address = tissue_position_address,
                   background_image_address = background_image_address,
                   giotto_python_path = giotto_python_path,
                   giotto_results_folder = giotto_results_folder)

  return(ICH_Mouse)

}

















