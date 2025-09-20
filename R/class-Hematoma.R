
# R/class-Hematoma.R

#' Hematoma: A class for identify the hematoma from spatial trascriptome dataset
#'
#' @slot analysis_symbol a symbol of sample
#' @slot file_address file address of dataset used in analysis
#' @slot tissue_position_matrix save tissue position matrix of dataset
#' @slot color_set a set of code for color
#' @slot raw_count_matrix the raw count matrix from dataset
#' @slot original_seu_metadata the original Seurat Object metadata
#' @slot seu_metadata_with_cluster_symbol Seurat Object metadata with cluster symbol
#' @slot spatial_image the spatial image with cluster
#' @slot identification_symbols the symbols identified by user as the the symbol of hematoma
#' @slot giotto_instruction the instruction of Giotto Object
#' @slot filtered_genes the genes that been filtered

setClass(
  Class = "Hematoma",
  slots = c(
    analysis_symbol = "character",
    file_address = "character",
    color_set = "character",
    tissue_position_matrix = "data.table",
    raw_count_matrix = "dgCMatrix",
    original_seu_metadata = "data.table",
    seu_metadata_with_cluster_symbol = "data.table",
    spatial_image = "list",
    identification_symbols = "list",
    giotto_instruction = "list",
    filtered_genes = "character")
)

#' Initialize class Hematoma
#'
#' @param analysis_symbol a symbol of sample
#' @param raw_count_matrix_address address for raw count matrix
#' @param filtered_count_matrix address for filtered count matrix
#' @param tissue_position_address address for tissue position matrix
#' @param background_image_address address for background address
#' @param giotto_python_path the path to a python which the Giotto can use
#' @param giotto_results_folder the folder for Giotto save plots
#' @param initialization whether start at beginning

setMethod(f = "initialize",
          signature = signature(.Object = "Hematoma"),
          definition = function(.Object,analysis_symbol,raw_count_matrix_address,filtered_count_matrix_address,tissue_position_address,background_image_address,giotto_python_path,giotto_results_folder,initialization) {

            on.exit(gc())

            if(initialization) {

              # add color set
              .Object@color_set <- c("#F5D2A8","#3C77AF","#7DBFA7","#EE934E","#9B5B33","#B383B9","#FCED82","#BBDD78","#F5D2A8","#D1352B","#8FA4AE","#F5CFE4","#D2EBC8","#F3AE63","#E69F84","#AA3538","#5891BF","#89558D","#79B99D","#AFC2D9","#D0AFC4","#C6307C","#E9E55A")
              # add analysis symbol
              .Object@analysis_symbol <- analysis_symbol
              # add file address
              .Object@file_address <- ICHMousewch ::: .integrate_file_address(raw_count_matrix_address = raw_count_matrix_address,
                                                                              filtered_count_matrix_address = filtered_count_matrix_address,
                                                                              tissue_position_address = tissue_position_address,
                                                                              background_image_address = background_image_address)

              # add Giotto instruction
              .Object@giotto_instruction <- ICHMousewch:::.create_giotto_instruction(python_path = giotto_python_path,
                                                                                     results_folder = giotto_results_folder) %>%
                list()

              # load tissue position matrix
              .Object@tissue_position_matrix <- ICHMousewch:::.load_tissue_position_matrix(tissue_position_address = tissue_position_address)

              # load raw count matrix
              .Object@raw_count_matrix <- ICHMousewch:::.load_raw_count_matrix(raw_count_matrix_address = raw_count_matrix_address)

              # generate original Seurat Object metadata
              .Object@original_seu_metadata <- ICHMousewch:::.generate_original_seu_metadata(raw_count_matrix = .Object@raw_count_matrix,
                                                                                             tissue_position_matrix = .Object@tissue_position_matrix)

              # filter genes
              .Object@filtered_genes <- ICHMousewch:::.find_filtered_genes(raw_count_matrix = .Object@raw_count_matrix,
                                                                           original_seu_metadata = .Object@original_seu_metadata)

              # generate Seurat Object with cluster symbol
              .Object@seu_metadata_with_cluster_symbol <- ICHMousewch:::.generate_seu_metadata_with_cluster_symbol(original_seu_metadata = .Object@original_seu_metadata,
                                                                                                                   raw_count_matrix = .Object@raw_count_matrix)

              # create spatial image plot
              .Object@spatial_image$GMM_cluster <-ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = .Object@seu_metadata_with_cluster_symbol,
                                                                                                          cluster_symbol = "GMM_cluster",
                                                                                                          raw_count_matrix = .Object@raw_count_matrix,
                                                                                                          background_image_address = .Object@file_address["background_image_address"],
                                                                                                          color_set = .Object@color_set,
                                                                                                          self_definition_color = c("1"="#F5D2A8","2"="#D1352B"),
                                                                                                          giotto_instruction = .Object@giotto_instruction[[1]]) %>%
                ggplotGrob()

              .Object@spatial_image$Louvain_cluster_posi <-ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = .Object@seu_metadata_with_cluster_symbol,
                                                                                                                   cluster_symbol = "Louvain_cluster_posi",
                                                                                                                   raw_count_matrix = .Object@raw_count_matrix,
                                                                                                                   background_image_address = .Object@file_address["background_image_address"],
                                                                                                                   color_set = .Object@color_set,
                                                                                                                   self_definition_color = c("1"="#F5D2A8"),
                                                                                                                   giotto_instruction = .Object@giotto_instruction[[1]]) %>%
                ggplotGrob()

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

            }

            validObject(.Object)
            return(.Object)

          })

#' Hematoma constructor
#'
#' @param analysis_symbol a symbol of sample
#' @param raw_count_matrix_address address for raw count matrix
#' @param filtered_count_matrix address for filtered count matrix
#' @param tissue_position_address address for tissue position matrix
#' @param background_image_address address for background address
#' @param initialization whether start at beginning
#' @export

Create_Hematoma <- function(analysis_symbol,raw_count_matrix_address,filtered_count_matrix_address,tissue_position_address,background_image_address,giotto_python_path,giotto_results_folder,initialization = TRUE) {

  on.exit(gc())

  Hematoma <- new(Class = "Hematoma",
                  analysis_symbol = analysis_symbol,
                  raw_count_matrix_address = raw_count_matrix_address,
                  filtered_count_matrix_address = filtered_count_matrix_address,
                  tissue_position_address = tissue_position_address,
                  background_image_address = background_image_address,
                  giotto_python_path = giotto_python_path,
                  giotto_results_folder = giotto_results_folder,
                  initialization = initialization)

  return(Hematoma)

}

####
#' the generic function of identify_hematoma
#'
#' @param hematoma the class of HematomaICHMousewch
#' @param hematoma_symbols the symbols identified by user as the the symbol of hematoma
#' @export

setGeneric(name = "identify_hematoma",
           def = function(hematoma,hematoma_symbols) {

             standardGeneric("identify_hematoma")

           })

#'  identify the exact hematoma location
#'
#'  @param hematoma the class of HematomaICHMousewch
#'  @param hematoma_symbols the symbols identified by user as the the symbol of hematoma

setMethod(f = "identify_hematoma",
          signature = signature(hematoma = "Hematoma",hematoma_symbols = "numeric"),
          definition = function(hematoma,hematoma_symbols) {

            on.exit(gc())

            hematoma@identification_symbols$hematoma_symbols <- hematoma_symbols

            hematoma@seu_metadata_with_cluster_symbol <- ICHMousewch:::.generate_seu_metadata_with_hematoma_symbol(hematoma = hematoma)

            hematoma@spatial_image$hematoma <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = hematoma@seu_metadata_with_cluster_symbol,
                                                                                                       cluster_symbol = "hematoma_symbol",
                                                                                                       raw_count_matrix = hematoma@raw_count_matrix,
                                                                                                       background_image_address = hematoma@file_address["background_image_address"],
                                                                                                       color_set = hematoma@color_set,
                                                                                                       self_definition_color = c("1"="#F5D2A8","2"="#D1352B"),
                                                                                                       giotto_instruction = hematoma@giotto_instruction[[1]]) %>%
              ggplotGrob()

            hematoma@spatial_image$Louvain_cluster_filt_gene <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = hematoma@seu_metadata_with_cluster_symbol,
                                                                                                                        cluster_symbol = "Louvain_cluster_filt_gene",
                                                                                                                        raw_count_matrix = hematoma@raw_count_matrix,
                                                                                                                        background_image_address = hematoma@file_address["background_image_address"],
                                                                                                                        color_set = hematoma@color_set,
                                                                                                                        self_definition_color = c("1"="#F5D2A8"),
                                                                                                                        giotto_instruction = hematoma@giotto_instruction[[1]]) %>%
              ggplotGrob()

            return(hematoma)

          })
####

####
#' identify hematoma center and edge
#'
#' @param hematoma the class of HematomaICHMousewch
#' @param center_symbols the symbols identified by user as the the symbol of hematoma center
#' @export

setGeneric(name = "identify_hematoma_center_and_edge",
           def = function(hematoma,center_symbols) {

             standardGeneric("identify_hematoma_center_and_edge")

           })


#'  identify the exact hematoma center and edge
#'
#'  @param hematoma the class of HematomaICHMousewch
#'  @param center_symbols the symbols identified by user as the the symbol of hematoma center

setMethod(f = "identify_hematoma_center_and_edge",
          signature = signature(hematoma = "Hematoma",center_symbols = "numeric"),
          definition = function(hematoma,center_symbols) {

            on.exit(gc())

            hematoma@identification_symbols$center_symbols <- center_symbols

            Louvain_cluster_filt_gene <- hematoma@seu_metadata_with_cluster_symbol[,Louvain_cluster_filt_gene] %>%
              unique() %>%
              unlist()

            hematoma@identification_symbols$edge_symbols <- Louvain_cluster_filt_gene[!Louvain_cluster_filt_gene %in% c(center_symbols,1)]

            hematoma@seu_metadata_with_cluster_symbol <- ICHMousewch:::.generate_seu_metadata_with_hematoma_center_and_edge_symbol(hematoma = hematoma)

            hematoma@spatial_image$hematoma_center_edge <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = hematoma@seu_metadata_with_cluster_symbol,
                                                                                                                   cluster_symbol = "center_edge_symbol",
                                                                                                                   raw_count_matrix = hematoma@raw_count_matrix,
                                                                                                                   background_image_address = hematoma@file_address["background_image_address"],
                                                                                                                   color_set = hematoma@color_set,
                                                                                                                   self_definition_color = c("1"="#F5D2A8","2"="#D1352B","3"="#3C77AF"),
                                                                                                                   giotto_instruction = hematoma@giotto_instruction[[1]]) %>%
              ggplotGrob()

            return(hematoma)

          })
####

####
#' generic function for create spatial image with highlighted clusters
#'
#' @param hematoma the class of HematomaICHMousewch
#' @param cluster_symbol the cluster symbol for plotiing
#' @param cluster_ls the cluster list for plot
#' @export

setGeneric(name = "create_spatial_image_with_highlighted_clusters",
           def = function(hematoma,cluster_symbol,cluster_ls) {

             standardGeneric("create_spatial_image_with_highlighted_clusters")

           })

#' create spatial image with highlighted clusters
#' @param hematoma the Hematoma Class
#' @param cluster_symbol the cluster symbol for plotiing
#' @param cluster_ls the cluster list for plot

setMethod(f = "create_spatial_image_with_highlighted_clusters",
          signature = signature(hematoma = "Hematoma",cluster_symbol = "character",cluster_ls = "numeric"),
          definition = function(hematoma,cluster_symbol,cluster_ls) {

            on.exit(gc())

            highlighted_color_set <- rep("black",length(cluster_ls))
            names(highlighted_color_set) <- as.character(cluster_ls)

            highlighted_spatial_image <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = hematoma@seu_metadata_with_cluster_symbol,
                                                                                                 cluster_symbol = cluster_symbol,
                                                                                                 raw_count_matrix = hematoma@raw_count_matrix,
                                                                                                 background_image_address = hematoma@file_address["background_image_address"],
                                                                                                 color_set = hematoma@color_set,
                                                                                                 self_definition_color = c("1"="#F5D2A8",highlighted_color_set),
                                                                                                 giotto_instruction = hematoma@giotto_instruction[[1]]) %>%
              ggplotGrob()

            grid.draw(highlighted_spatial_image)

          })
####

####
#' save Hematoma class
#'
#' @param hematoma the Hematoma class
#' @param saving_path the path to save
#' @export

setGeneric(name = "save_Hematoma",
           def = function(hematoma,saving_path) {

             standardGeneric("save_Hematoma")

           })

#' save Hematoma class
#'
#' @param hematoma the Hematoma class
#' @param saving_path the path to save

setMethod(f = "save_Hematoma",
          signature = signature(hematoma = "Hematoma",saving_path = "character"),
          definition = function(hematoma,saving_path) {

            on.exit(gc())

            file_path <- paste(saving_path,"DataBase",sep = "/")

            if(!dir.exists(file_path)) {

              dir.create(file_path)

            }

            slot_na <- slotNames(hematoma)

            file_na_ls <- list()
            for (i in 1:length(slot_na)) {

              slot_da <- slot(hematoma,slot_na[i])

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
#' load Hematoma from local
#'
#' @param loading_path the path saving Hematoma
#' @export

setGeneric(name = "load_Hematoma",
           def = function(loading_path) {

             standardGeneric("load_Hematoma")

           })

#' load Hematoma from local
#'
#' @param loading_path the path saving Hematoma

setMethod(f = "load_Hematoma",
          signature = signature(loading_path = "character"),
          definition = function(loading_path) {

            on.exit(gc())

            file_path <- paste(loading_path,"DataBase",sep = "/")

            hematoma <- new(Class = "Hematoma",
                            analysis_symbol = NULL,
                            raw_count_matrix_address = NULL,
                            filtered_count_matrix_address = NULL,
                            tissue_position_address = NULL,
                            background_image_address = NULL,
                            giotto_python_path = NULL,
                            giotto_results_folder = NULL,
                            initialization = FALSE)

            dataset_name <- readRDS(file = paste(file_path,"dataset_name.rds",sep = "/"))

            for (i in 1:length(dataset_name)) {

              slot(object = hematoma,name = names(dataset_name)[i]) <- readRDS(file = paste(file_path,dataset_name[i],sep = "/"))

            }

            return(hematoma)

            })
####

####

#' show ggplotGrob image
#'
#' @param ggplot_image a ggplot image in the form of gtable

setGeneric(name = "show_image",
           def = function(ggplot_image) {

             standardGeneric("show_image")

           })

#' show ggplotGrob image
#'
#' @param ggplot_image a ggplot image in the form of gtable
#' @export

setMethod(f = "show_image",
          definition = function(ggplot_image) {

            on.exit(gc())

            grid.draw(ggplot_image)

          })
####
