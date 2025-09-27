
# R/class-Spatial_Image.R

#' the class fot save Giottto Spatial Image
#'
#' @slot image_set_name the name of image set
#' @slot gene_name a gene list for plotting the spatial image
#' @slot spatial_image a list of spatial image

setClass(Class = "Spatial_Image",
         slots = c(image_set_name = "character",
                   gene_name = "character",
                   spatial_image = "list"))

#' Initialize the class of Spatial_Image
#'
#' @param spatial_image_ls a list contain a set of spatial image
#' @param image_set_name the name of image set

setMethod(f = "initialize",
          signature = signature(.Object = "Spatial_Image"),
          definition = function(.Object,image_set_name,spatial_image_ls) {

            .Object@image_set_name <- image_set_name

            .Object@gene_name <- names(spatial_image_ls)

            .Object@spatial_image <- spatial_image_ls

            validObject(.Object)
            return(.Object)

          })

#' Spatial Image constructor
#'
#' @param spatial_image_ls a list contain a set of spatial image
#' @param image_set_name the name of image set

.Create_Spatial_Image <- function(image_set_name,spatial_image_ls) {

  on.exit(gc())

  spatial_image <- new(Class = "Spatial_Image",
                       image_set_name = image_set_name,
                       spatial_image_ls = spatial_image_ls)

  return(spatial_image)

}

####
#' save spatial image
#'
#' @param saving_path the path for saving
#' @param spatial_image the class of Spatial Image

setGeneric(name = "save_spatial_image",
           def = function(spatial_image,saving_path,height,width) {

             standardGeneric("save_spatial_image")

           })

#' save spatial image
#'
#' @param saving_path the path for saving
#' @param spatial_image the class of Spatial Image
#' @export

setMethod(f = "save_spatial_image",
          signature = signature(spatial_image = "Spatial_Image",saving_path = "character"),
          definition = function(spatial_image,saving_path,height = 8,width = 8) {

            on.exit(gc())

            ICHMousewch:::.save_spatial_image(spatial_image = spatial_image,
                                              saving_path = saving_path,
                                              width = width,
                                              height = height)

          })
####

####
#' create single gene spatial image
#'
#' @param ich_mouse the class of ICH_Mouse
#' @param gene_ls a gene list
#' @param show_background_image whether to show background image
#' @param image_set_name the name of the image set

setGeneric(name = "create_single_gene_spatial_image",
           def = function(ich_mouse,gene_ls,image_set_name,show_background_image = TRUE) {

             standardGeneric("create_single_gene_spatial_image")

           })

#' create single gene spatial image
#'
#' @param ICH_Mouse the class of ICH_Mouse
#' @param gene_ls a gene list
#' @param show_background_image whether to show background image
#' @param image_set_name the name of the image set
#' @export

setMethod(f = "create_single_gene_spatial_image",
          signature = signature(ich_mouse = "ICH_Mouse",gene_ls = "character",image_set_name = "character"),
          definition = function(ich_mouse,gene_ls,image_set_name,show_background_image = TRUE) {

            on.exit(gc())

            spatial_image_ls <- ICHMousewch:::.create_spatial_image_with_single_gene(seu_metadata_with_cluster_symbol = ich_mouse@seu_metadata_with_cluster_symbol,
                                                                                     gene_ls = gene_ls,
                                                                                     raw_count_matrix = ich_mouse@raw_count_matrix,
                                                                                     background_image_address = ich_mouse@file_address["background_image_address"],
                                                                                     giotto_instruction = ich_mouse@giotto_instruction[[1]],
                                                                                     show_background_image = show_background_image)

            spatial_image_set <- ICHMousewch:::.Create_Spatial_Image(image_set_name = image_set_name,
                                                                     spatial_image_ls = spatial_image_ls)

            return(spatial_image_set)

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

            show(as.ggplot(ggplot_image))

          })
####

####
#' create single gene contrasting spatial image
#'
#' @param ich_mouse1 the class of ICH_Mouse
#' @param ich_mouse2 the class of ICH_Mouse
#' @param gene_ls a gene list
#' @param show_background_image whether to show background image
#' @param image_set_name the name of the image set

setGeneric(name = "create_single_gene_contrasting_spatial_image",
           def = function(ich_mouse1,ich_mouse2,gene_ls,image_set_name,show_background_image = TRUE) {

             standardGeneric("create_single_gene_contrasting_spatial_image")

           })

#' create single contrasting gene spatial image
#'
#' @param ich_mouse1 the class of ICH_Mouse
#' @param ich_mouse2 the class of ICH_Mouse
#' @param gene_ls a gene list
#' @param show_background_image whether to show background image
#' @param image_set_name the name of the image set
#' @export

setMethod(f = "create_single_gene_contrasting_spatial_image",
          signature = signature(ich_mouse1 = "ICH_Mouse",ich_mouse2 = "ICH_Mouse",gene_ls = "character",image_set_name = "character"),
          definition = function(ich_mouse1,ich_mouse2,gene_ls,image_set_name,show_background_image = TRUE) {

            on.exit(gc())

            spatial_image_ls1 <- ICHMousewch:::.create_spatial_image_with_single_gene(seu_metadata_with_cluster_symbol = ich_mouse1@seu_metadata_with_cluster_symbol,
                                                                                      gene_ls = gene_ls,
                                                                                      raw_count_matrix = ich_mouse1@raw_count_matrix,
                                                                                      background_image_address = ich_mouse1@file_address["background_image_address"],
                                                                                      giotto_instruction = ich_mouse1@giotto_instruction[[1]],
                                                                                      show_background_image = show_background_image)

            spatial_image_ls2 <- ICHMousewch:::.create_spatial_image_with_single_gene(seu_metadata_with_cluster_symbol = ich_mouse2@seu_metadata_with_cluster_symbol,
                                                                                      gene_ls = gene_ls,
                                                                                      raw_count_matrix = ich_mouse2@raw_count_matrix,
                                                                                      background_image_address = ich_mouse2@file_address["background_image_address"],
                                                                                      giotto_instruction = ich_mouse2@giotto_instruction[[1]],
                                                                                      show_background_image = show_background_image)

            spatial_image_ls <- vector("list",length = length(gene_ls))
            names(spatial_image_ls) <- gene_ls
            for (i in 1:length(gene_ls)) {

              spatial_image <- wrap_plots(as.ggplot(spatial_image_ls1[[gene_ls[i]]]),as.ggplot(spatial_image_ls2[[gene_ls[i]]])) %>%
                patchworkGrob()

              spatial_image_ls[gene_ls[i]] <- list(spatial_image)

            }

            spatial_image_set <- ICHMousewch:::.Create_Spatial_Image(image_set_name = image_set_name,
                                                                     spatial_image_ls = spatial_image_ls)

            return(spatial_image_set)

          })
####






