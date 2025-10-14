
# R/class-Plotting_Dataset.R

#' the class for Plotting
#'
#' @slot spatial_image the spatial image
#' @slot spatial_single_gene_image the spatial single gene image
#' @slot violin_plot the violin plot
#' @slot volcano_plot the volcano plot
#' @slot heatmap the heatmap
#' @slot bubble_chart the bubble chart
#' @slot barchart the bar chart

setClass(Class = "Plotting_Dataset",
         slots = c(spatial_image = "list",
                   spatial_single_gene_image = "list",
                   violin_plot = "list",
                   volcano_plot = "list",
                   heatmap = "list",
                   bubble_chart = "list",
                   barchart = "list"))

#' initialize Plotting_Dataset
#'
#' @param ich_mouse the class of ICH_Mouse

setMethod(f = "initialize",
          signature = signature(.Object = "Plotting_Dataset"),
          definition = function(.Object,ich_mouse,initialization) {

            on.exit(gc())

            if (initialization) {

            # create spatial image
              .Object@spatial_image <- ICHMousewch:::.plotting_spatial_image(ich_mouse = ich_mouse)

            # create volcano plot
              .Object@volcano_plot <- ICHMousewch:::.plotting_volcano_plot(ich_mouse = ich_mouse)

            # create heatmap
              .Object@heatmap <- ICHMousewch:::.plotting_GO_term_similarity_heatmap(ich_mouse = ich_mouse)

            # create barchart
              .Object@barchart <- ICHMousewch:::.plotting_GMM_barchart(ich_mouse = ich_mouse)

            # create violin plot
              .Object@violin_plot <- ICHMousewch:::.plotting_violin_plot(ich_mouse =ich_mouse)

            # create bubble chart
              .Object@bubble_chart <- ICHMousewch:::.plotting_bubble_chart(ich_mouse = ich_mouse)

            } else {

              .Object@spatial_image <- list()
              .Object@spatial_single_gene_image <- list()
              .Object@violin_plot <- list()
              .Object@volcano_plot <- list()
              .Object@heatmap <- list()
              .Object@bubble_chart <- list()
              .Object@barchart <- list()

            }



            validObject(.Object)
            return(.Object)

          })

#' constructor of Plotting_Dataset
#'
#' @param ich_mouse the class of ICH_Mouse
#' @export

Create_Plotting_Dataset <- function(ich_mouse,initialization = TRUE) {

  on.exit(gc())

  plotting_dataset <- new(Class = "Plotting_Dataset",
                          ich_mouse = ich_mouse,
                          initialization = initialization)

  return(plotting_dataset)

}

####
#' save Plotting_Dataset class
#'
#' @param plotting_dataset the Plotting_Dataset class
#' @param saving_path the path to save
#' @export

setGeneric(name = "save_Plotting_Dataset",
           def = function(plotting_dataset,saving_path) {

             standardGeneric("save_Plotting_Dataset")

           })

#' save Plotting_Dataset class
#'
#' @param plotting_dataset the Plotting_Dataset class
#' @param saving_path the path to save

setMethod(f = "save_Plotting_Dataset",
          signature = signature(plotting_dataset = "Plotting_Dataset",saving_path = "character"),
          definition = function(plotting_dataset,saving_path) {

            on.exit(gc())

            file_path <- paste(saving_path,"DataBase-plotting_dataset",sep = "/")

            if(!dir.exists(file_path)) {

              dir.create(file_path,recursive = TRUE)

            } else {

              unlink(file_path,recursive = TRUE)
              dir.create(file_path,recursive = TRUE)

            }

            slot_na <- slotNames(plotting_dataset)

            file_na_ls <- list()
            for (i in 1:length(slot_na)) {

              slot_da <- slot(plotting_dataset,slot_na[i])

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
#' load Plotting_Dataset from local
#'
#' @param loading_path the path saving Plotting_Dataset
#' @export

setGeneric(name = "load_Plotting_Dataset",
           def = function(loading_path) {

             standardGeneric("load_Plotting_Dataset")

           })

#' load Plotting_Dataset from local
#'
#' @param loading_path the path saving Plotting_Dataset

setMethod(f = "load_Plotting_Dataset",
          signature = signature(loading_path = "character"),
          definition = function(loading_path) {

            on.exit(gc())

            file_path <- paste(loading_path,"DataBase-plotting_dataset",sep = "/")

            plotting_dataset <- new(Class = "Plotting_Dataset",
                                    ich_mouse = NULL,
                                    initialization = FALSE)

            dataset_name <- readRDS(file = paste(file_path,"dataset_name.rds",sep = "/"))

            for (i in 1:length(dataset_name)) {

              slot(object = plotting_dataset,name = names(dataset_name)[i]) <- readRDS(file = paste(file_path,dataset_name[i],sep = "/"))

            }

            return(plotting_dataset)

          })
####

####
#' export plotting dataset
#'
#' @param plotting_dataset the Plotting_Dataset class
#' @param saving_path the path for saving plotting_dataset

setGeneric(name = "export_plotting_dataset",
           def = function(plotting_dataset,saving_path) {

             standardGeneric("export_plotting_dataset")

           })

#' export plotting dataset
#'
#' @param plotting_dataset the Plotting_Dataset class
#' @param saving_path the path for saving plotting_dataset
#' @export

setMethod(f = "export_plotting_dataset",
          signature = signature(plotting_dataset = "Plotting_Dataset",saving_path = "character"),
          definition = function(plotting_dataset,saving_path) {

            on.exit(gc())

            saving_directory <- paste(saving_path,"ICH_Mouse_plots",sep = "/")

            if(!dir.exists(saving_directory)) {

              dir.create(saving_directory,recursive = TRUE)

            } else {

              unlink(saving_directory,recursive = TRUE)
              dir.create(saving_directory,recursive = TRUE)

            }

            slot_na <- slotNames(plotting_dataset)

            for (i in 1:length(slot_na)) {

              slot_da <- slot(plotting_dataset,slot_na[i])

              if (length(slot_da) != 0) {

                plot_na <- names(slot_da)

                slot_da_path <- paste(saving_directory,slot_na[i],sep = "/")

                dir.create(slot_da_path,recursive = FALSE)

                for (j in 1:length(plot_na)) {

                  saving_plotting <- as.ggplot(slot_da[[plot_na[j]]])

                  ggsave(filename = paste(plot_na[j],"png",sep = "."),
                         plot = saving_plotting,
                         device = "png",
                         path = slot_da_path,
                         dpi = 600,
                         width = (297/25.4),
                         height = (210/25.4),
                         unit = "in")

                }

              }

            }

          })

####

