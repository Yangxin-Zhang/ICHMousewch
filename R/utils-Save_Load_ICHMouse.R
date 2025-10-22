
# R/Save_ICHMouse.R

####
#' save ICHMouse data
#' @param Class_ICHMouse the class of ICHMouse
#' @param saving_path the path for saving

setGeneric(name = "Save_ICHMouse",
           def = function(Class_ICHMouse,data_symbol = "ich_mouse",saving_path = NULL) {

             standardGeneric("Save_ICHMouse")

           })

#' load ICHMouse data
#'
#' @param loading_path the path saving data
#' @param data_symbol the data symbol

setGeneric(name = "Load_ICHMouse",
           def = function(data_symbol = "ich_mouse",loading_path = NULL) {

             standardGeneric("Load_ICHMouse")

           })
####

####
#' save Hematoma
#' @param Class_ICHMouse the class of ICHMouse
#' @param saving_path the path for saving
#' @param data_symbol the saving symbol
#' @export

setMethod(f = "Save_ICHMouse",
          signature = signature(Class_ICHMouse = "Hematoma"),
          definition = function(Class_ICHMouse,data_symbol = "ich_mouse",saving_path = NULL) {

            on.exit(gc())

            if (is.null(saving_path)) {

              saving_path <- getwd()

            }

            hematoma <- Class_ICHMouse

            file_path <- paste(saving_path,"ICHMouse_Database",sep = "/") %>%
              paste(data_symbol,sep = "/") %>%
              paste("Class-Hematoma",sep = "/")

            if(!dir.exists(file_path)) {

              dir.create(file_path,recursive = TRUE)

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

#' load Hematoma
#'
#' @param loading_path the path loading data
#' @param data_symbol the data symbol
#' @export

setMethod(f = "Load_ICHMouse",
          definition = function(data_symbol = "ich_mouse",loading_path = NULL) {

            on.exit(gc())

            if (is.null(loading_path)) {

              loading_path <- getwd()

            }

            file_path <- paste(loading_path,"ICHMouse_Database",sep = "/") %>%
              paste(data_symbol,sep = "/") %>%
              paste("Class-Hematoma",sep = "/")

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
#' save ICH_Mouse
#' @param Class_ICHMouse the class of ICHMouse
#' @param saving_path the path for saving
#' @param data_symbol the saving symbol
#' @export

setMethod(f = "Save_ICHMouse",
          signature = signature(Class_ICHMouse = "ICH_Mouse"),
          definition = function(Class_ICHMouse,data_symbol = "ich_mouse",saving_path = NULL) {

            on.exit(gc())

            if (is.null(saving_path)) {

              saving_path <- getwd()

            }

            file_path <- paste(saving_path,"ICHMouse_Database",sep = "/") %>%
              paste(data_symbol,sep = "/") %>%
              paste("Class-ICH_Mouse",sep = "/")

            ich_mouse <- Class_ICHMouse

            if(!dir.exists(file_path)) {

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

#' load ICH_Mouse
#'
#' @param loading_path the path loading data
#' @param data_symbol the data symbol
#' @export

setMethod(f = "Load_ICHMouse",
          definition = function(data_symbol = "ich_mouse",loading_path = NULL) {

            on.exit(gc())

            if (is.null(loading_path)) {

              loading_path <- getwd()

            }

            file_path <- paste(loading_path,"ICHMouse_Database",sep = "/") %>%
              paste(data_symbol,sep = "/") %>%
              paste("Class-ICH_Mouse",sep = "/")

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

####
#' export plotting dataset
#'
#' @param plotting_dataset the Plotting_Dataset class
#' @param saving_path the path for saving plotting_dataset

setGeneric(name = "export_plotting_dataset",
           def = function(plotting_dataset,saving_path) {

             standardGeneric("export_plotting_dataset")

           })


###
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
###

###
#' export spatial image dataset
#'
#' @param plotting_dataset the Spatial_Image class
#' @param saving_path the path for saving Spatial_Image
#' @export

setMethod(f = "export_plotting_dataset",
          signature = signature(plotting_dataset = "Spatial_Image",saving_path = "character"),
          definition = function(plotting_dataset,saving_path) {

            on.exit(gc())

            saving_directory <- paste(saving_path,"single_gene_spatial_image",sep = "/")

            if(!dir.exists(saving_directory)) {

              dir.create(saving_directory,recursive = TRUE)

            }

            slot_na <- "spatial_image"
            slot_da <- slot(plotting_dataset,slot_na)

            if (length(slot_da) != 0) {

              plot_na <- names(slot_da)

              slot_da_path <- paste(saving_directory,plotting_dataset@image_set_name,sep = "/")

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

          })
###

###
####
