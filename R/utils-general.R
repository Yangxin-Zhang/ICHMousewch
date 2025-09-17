
# R/utils-general.R

#' integrate file address
#'
#' @param raw_count_matrix_address address for raw count matrix
#' @param filtered_count_matrix address for filtered count matrix
#' @param tissue_position_address address for tissue position matrix
#' @param background_image_address address for background address

.integrate_file_address <- function(raw_count_matrix_address,filtered_count_matrix_address,tissue_position_address,background_image_address) {

  on.exit(gc())

  file_address <- list("raw_count_matrix_address" = raw_count_matrix_address,
                       "filtered_count_matrix_address" = filtered_count_matrix_address,
                       "tissue_position_address" = tissue_position_address,
                       "background_image_address" = background_image_address)

  return(file_address)

}

#' load tissue position matrix
#'
#' @param tissue_position_address the address of tissue position matrix

.load_tissue_position_matrix <- function(tissue_position_address) {

  on.exit(gc())

  tissue_position_matrix <- Read10X_Coordinates(tissue_position_address,
                                                filter.matrix = FALSE)

  tissue_position_matrix$barcode <- rownames(tissue_position_matrix)

  return(tissue_position_matrix)

}

#' load raw count matrix
#'
#' @param raw_count_matrix_address the file address of raw count matrix

.load_raw_count_matrix <- function(raw_count_matrix_address) {

  on.exit(gc())

  raw_count_matrix <- Read10X_h5(raw_count_matrix_address)

  return(raw_count_matrix)

}

#' save barcodes
#'
#' @param tissue_position_matrix the matrix of tissue position
#' @param in_tissue the location of barcodes value: TRUE or FALSE

.save_barcodes <- function(tissue_position_matrix,in_tissue) {

  on.exit(gc())

  if(in_tissue) {

    barcodes <- tissue_position_matrix[tissue_position_matrix$tissue == 1,]$barcode

  } else {

    barcodes <- tissue_position_matrix[tissue_position_matrix$tissue == 0,]$barcode

  }

  return(barcodes)

}

#' find filtered genes
#'
#' @param raw_count_matrix the matrix of raw count dataset
#' @param in_tissue_barcodes the barcodes in tissue

.find_filtered_genes <- function(raw_count_matrix,in_tissue_barcodes) {

  on.exit(gc())

  filtered_count_matrix <- raw_count_matrix[,in_tissue_barcodes]

  num_gene <- filtered_count_matrix %>%
    rownames() %>%
    length()

  bins <- seq(1,num_gene,5000)

  per_ls <- list()
  for (i in 1:length(bins)) {

    if (i != length(bins)) {

      per <- Matrix::rowSums(filtered_count_matrix[bins[i]:(bins[i]+4999),] == 0)/length(colnames(filtered_count_matrix))
      per_ls <- append(per_ls,list(per))
      names(per_ls)[i] <- paste(bins[i],(bins[i]+4999),sep = ":")

    } else {

      per <- Matrix::rowSums(filtered_count_matrix[bins[i]:num_gene,] == 0)/length(colnames(filtered_count_matrix))
      per_ls <- append(per_ls,list(per))
      names(per_ls)[i] <- paste(bins[i],num_gene,sep = ":")

    }

  }

  filtered_gene <- list()
  for (i in 1:length(per_ls)) {

    logi <- per_ls[[i]] >= 0.99
    filtered_gene <- append(filtered_gene,per_ls[[i]][!logi])

  }

  return(names(filtered_gene))

}


