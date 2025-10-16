
# R/Save_ICHMouse.R

#' save ICHMouse data
#' @param Class_ICHMouse the class of ICHMouse
#' @param saving_path the path for saving

setGeneric(name = "Save_ICHMouse",
           def = function(Class_ICHMouse,saving_path = getwd()) {

             standardGeneric("Save_ICHMouse")

           })

#' save Hematoma
#' @param Class_ICHMouse the class of ICHMouse
#' @param saving_path the path for saving
#' @export

setMethod(f = "Save_ICHMouse",
          signature = signature(Class_ICHMouse = "Hematoma",saving_path = "character"),
          definition = function(Class_ICHMouse,saving_path) {

            on.exit(gc())

            file_path <- paste(saving_path,"DataBase",sep = "/")

            if(!dir.exists(file_path)) {

              dir.create(file_path,recursive = TRUE)

            } else {

              unlink(file_path,recursive = TRUE)
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
