
# R/class-ICH_Mouse.R

#' ICH_Mouse a class for analysis ICH_Mouse
#'
#' @slot diff_expr_genes the differential expression genes
#' @slot symbol_genes the symbol genes

setClass(Class = "ICH_Mouse",
         slots = c(diff_expr_genes = "list",
                   symbol_genes = "list",
                   spatial_image_with_single_gene = "list"),
         contains = "Hematoma")

#' Initialize class ICH_Mouse
#'
#' @param analysis_symbol a symbol of sample
#' @param raw_count_matrix_address address for raw count matrix
#' @param filtered_count_matrix address for filtered count matrix
#' @param tissue_position_address address for tissue position matrix
#' @param background_image_address address for background address
#' @param giotto_python_path the path to a python which the Giotto can use
#' @param giotto_results_folder the folder for Giotto save plots
#' @param initialization whether start at beginning
#' @param hematoma_symbols the symbols identified by user as the the symbol of hematoma
#' @param center_symbols the symbols identified by user as the the symbol of hematoma center

setMethod(f = "initialize",
          signature = signature(.Object = "ICH_Mouse"),
          definition = function(.Object,analysis_symbol,raw_count_matrix_address,filtered_count_matrix_address,tissue_position_address,background_image_address,giotto_python_path,giotto_results_folder,hematoma_symbols,center_symbols,initialization) {

            on.exit(gc())

            if (initialization) {

              .Object <- callNextMethod(.Object,analysis_symbol,raw_count_matrix_address,filtered_count_matrix_address,tissue_position_address,background_image_address,giotto_python_path,giotto_results_folder,initialization)

              .Object <- ICHMousewch:::identify_hematoma(hematoma = .Object,
                                                         hematoma_symbols = hematoma_symbols)

              .Object <- ICHMousewch:::identify_hematoma_center_and_edge(hematoma = .Object,
                                                                         center_symbols = center_symbols)

              .Object@symbol_genes$normal_tissue <- ICHMousewch:::.find_symbol_genes(raw_count_matrix = .Object@raw_count_matrix,
                                                                                     seu_metadata_with_cluster_symbol = .Object@seu_metadata_with_cluster_symbol,
                                                                                     filtered_genes = .Object@filtered_genes,
                                                                                     cluster_symbol = 1)

              # .Object@diff_expr_genes["edge-normal"]<- ICHMousewch:::.find_differential_expression_genes(raw_count_matrix = .Object@raw_count_matrix,
              #                                                                                            seu_metadata_with_cluster_symbol = .Object@seu_metadata_with_cluster_symbol,
              #                                                                                            filtered_genes = .Object@filtered_genes,
              #                                                                                            cluster_symbol = c(3,1)) %>%
              #   list()
              #
              # .Object@diff_expr_genes["center-normal"]<- ICHMousewch:::.find_differential_expression_genes(raw_count_matrix = .Object@raw_count_matrix,
              #                                                                                              seu_metadata_with_cluster_symbol = .Object@seu_metadata_with_cluster_symbol,
              #                                                                                              filtered_genes = .Object@filtered_genes,
              #                                                                                              cluster_symbol = c(2,1)) %>%
              #   list()
              #
              # .Object@diff_expr_genes["center-edge"]<- ICHMousewch:::.find_differential_expression_genes(raw_count_matrix = .Object@raw_count_matrix,
              #                                                                                            seu_metadata_with_cluster_symbol = .Object@seu_metadata_with_cluster_symbol,
              #                                                                                            filtered_genes = .Object@filtered_genes,
              #                                                                                            cluster_symbol = c(2,3)) %>%
              #   list()

            } else {

              .Object@analysis_symbol = character()
              .Object@file_address = character()
              .Object@color_set = character()
              .Object@tissue_position_matrix = data.table()
              .Object@raw_count_matrix = sparseMatrix(i = integer(0),j = integer(0),x = numeric(0))
              .Object@original_seu_metadata = data.table()
              .Object@seu_metadata_with_cluster_symbol = data.table()
              .Object@spatial_image = list()
              .Object@identification_symbols = list()
              .Object@giotto_instruction = list()
              .Object@filtered_genes = character()
              .Object@diff_expr_genes = list()
              .Object@symbol_genes = list()
              .Object@spatial_image_with_single_gene = list()

            }


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
#' @param initialization whether start at beginning
#' @param hematoma_symbols the symbols identified by user as the the symbol of hematoma
#' @param center_symbols the symbols identified by user as the the symbol of hematoma center
#' @export

Create_ICH_Mouse <- function(analysis_symbol,raw_count_matrix_address,filtered_count_matrix_address,tissue_position_address,background_image_address,giotto_python_path,giotto_results_folder,hematoma_symbols,center_symbols,initialization = TRUE) {

  on.exit(gc())

  ICH_Mouse <- new(Class = "ICH_Mouse",
                   analysis_symbol = analysis_symbol,
                   raw_count_matrix_address = raw_count_matrix_address,
                   filtered_count_matrix_address = filtered_count_matrix_address,
                   tissue_position_address = tissue_position_address,
                   background_image_address = background_image_address,
                   giotto_python_path = giotto_python_path,
                   giotto_results_folder = giotto_results_folder,
                   hematoma_symbols = hematoma_symbols,
                   center_symbols = center_symbols,
                   initialization = initialization)

  return(ICH_Mouse)

}

####
#' find differential expression genes
#'
#' @param ich_mouse the ICH_Mouse class
#' @param cluster_symbol the cluster to contrast

setGeneric(name = "find_differential_expression_genes",
           def = function(ich_mouse,cluster_symbol,gene_set_name) {

             standardGeneric("find_differential_expression_genes")

           })

#' find differential expression genes
#'
#' @param ich_mouse the ICH_Mouse class
#' @param cluster_symbol the cluster to contrast
#' @export

setMethod(f = "find_differential_expression_genes",
          signature = signature(ich_mouse = "ICH_Mouse",cluster_symbol = "numeric",gene_set_name = "character"),
          definition = function(ich_mouse,cluster_symbol,gene_set_name) {

            on.exit(gc())

            diff_expr_gene <- ICHMousewch:::.find_differential_expression_genes(raw_count_matrix = ich_mouse@raw_count_matrix,
                                                                                seu_metadata_with_cluster_symbol = ich_mouse@seu_metadata_with_cluster_symbol,
                                                                                filtered_genes = ich_mouse@filtered_genes,
                                                                                cluster_symbol = cluster_symbol)

            annotation_dt <- ICHMousewch:::.annotate_the_cell_type_based_on_single_gene(ich_mouse = ich_mouse,
                                                                                        gene_ls = diff_expr_gene[,gene_name])

            gene_order <- match(diff_expr_gene[,gene_name],annotation_dt[,gene_name])

            annotation_dt <- annotation_dt[gene_order]

            diff_expr_gene[,cell_type := annotation_dt[,cell_type]]

            ich_mouse@diff_expr_genes[gene_set_name] <- list(diff_expr_gene)

            return(ich_mouse)

          })
####

####
#' save ICH_Mouse class
#'
#' @param ich_mouse the ICH_Mouse class
#' @param saving_path the path to save
#' @export

setGeneric(name = "save_ICH_Mouse",
           def = function(ich_mouse,saving_path) {

             standardGeneric("save_ICH_Mouse")

           })

#' save ICH_Mouse class
#'
#' @param ich_mouse the ICH_Mouse class
#' @param saving_path the path to save

setMethod(f = "save_ICH_Mouse",
          signature = signature(ich_mouse = "ICH_Mouse",saving_path = "character"),
          definition = function(ich_mouse,saving_path) {

            on.exit(gc())

            file_path <- paste(saving_path,"DataBase",sep = "/")

            if(!dir.exists(file_path)) {

              dir.create(file_path,recursive = TRUE)

            } else {

              unlink(file_path,recursive = TRUE)
              dir.create(file_path,recursive = TRUE)

            }

            slot_na <- slotNames(ich_mouse)

            file_na_ls <- list()
            for (i in 1:length(slot_na)) {

              slot_da <- slot(ich_mouse,slot_na[i])

              file_name <- paste(slot_na[i],"rds",sep = ".")

              file_na_ls <- append(file_na_ls,list(file_name))
              names(file_na_ls)[i] <- slot_na[i]

              saveRDS(object = slot_da,
                      file = paste(file_path,file_name,sep = "/"),
                      compress = FALSE)

            }

            file_name <- paste("dataset_name","rds",sep = ".")

            file_na_ls <- unlist(file_na_ls)
            saveRDS(object = file_na_ls,
                    file = paste(file_path,file_name,sep = "/"),
                    compress = FALSE)

          })
####

####
#' load ICH_Mouse from local
#'
#' @param loading_path the path saving ICH_Mouse
#' @export

setGeneric(name = "load_ICH_Mouse",
           def = function(loading_path) {

             standardGeneric("load_ICH_Mouse")

           })

#' load ICH_Mouse from local
#'
#' @param loading_path the path saving ICH_Mouse

setMethod(f = "load_ICH_Mouse",
          signature = signature(loading_path = "character"),
          definition = function(loading_path) {

            on.exit(gc())

            file_path <- paste(loading_path,"DataBase",sep = "/")

            ich_mouse <- new(Class = "ICH_Mouse",
                             analysis_symbol = NULL,
                             raw_count_matrix_address = NULL,
                             filtered_count_matrix_address = NULL,
                             tissue_position_address = NULL,
                             background_image_address = NULL,
                             giotto_python_path = NULL,
                             giotto_results_folder = NULL,
                             hematoma_symbols = NULL,
                             center_symbols = NULL,
                             initialization = FALSE)

            dataset_name <- readRDS(file = paste(file_path,"dataset_name.rds",sep = "/"))

            for (i in 1:length(dataset_name)) {

              slot(object = ich_mouse,name = names(dataset_name)[i]) <- readRDS(file = paste(file_path,dataset_name[i],sep = "/"))

            }

            return(ich_mouse)

          })
####







