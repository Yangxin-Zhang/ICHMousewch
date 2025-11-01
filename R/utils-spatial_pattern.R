
# R/utils-spatial_pattern.R

#' generate the distance matrix of barcodes
#'
#' @param seu_meta the Seurat Object metadata
#' @param spatial_chunks the chunks for calculation
#' @param kernel_function choose the kernel function

.genetate_distance_matrix <- function(seu_meta,spatial_chunks,kernel_function = "Gaussian") {

  on.exit(gc())

  # calculate standard distance

  sub_seu_meta <- seu_meta[cell_ID %in% spatial_chunks]

  sub_seu_meta[,st_x := as.numeric(scale(imagerow))]
  sub_seu_meta[,st_y := as.numeric(scale(imagecol))]

  # import kernel function

  if (kernel_function == "Gaussian") {

    ker_fun <- rbfdot(sigma = 0.5)

  }

  dt_mt <- sub_seu_meta[,c("st_x","st_y")] %>%
    as.matrix()

  kernel_matrix <- kernelMatrix(ker_fun,dt_mt)

  rownames(kernel_matrix) <- spatial_chunks
  colnames(kernel_matrix) <- spatial_chunks

  return(kernel_matrix)

}

#' create chunks
#'
#' @param total_num the total number for creating chunks
#' @param bin_num the number of chunks

.create_chunks <- function(total_num,bin_num) {

  on.exit(gc())

  chunk_scale <- total_num %/% bin_num
  rest_value <- total_num %% bin_num

  chunk_edge <- seq(from = 1,
                    to = total_num,
                    by = chunk_scale)

  chunk_ls <- vector("list",length = bin_num)
  names(chunk_ls) <- paste("chunk",c(1:bin_num),sep = "_")
  for (i in seq(from = 1,to = bin_num,by = 1)) {

    if (i == bin_num) {

      sub_chunk <- seq(from = chunk_edge[i],
                       to = total_num,
                       by = 1) %>%
        list()

    } else {

      sub_chunk <- seq(from = chunk_edge[i],
                       to = (chunk_edge[i+1]-1),
                       by = 1) %>%
        list()

    }

    chunk_ls[i] <- sub_chunk

  }

  return(chunk_ls)

}

#' generate chunk barcodes
#'
#' @param barcodes the barcode for chunk
#' @param chunk_bins the bins of chunk set

.genereate_chunk_barcodes <- function(barcodes,chunk_bins) {

  on.exit(gc())

  set.seed(2025)

  barcodes_sample <- barcodes

  chunk_name <- names(chunk_bins)

  barcodes_ls <- vector("list",length = length(chunk_name))
  names(barcodes_ls) <- chunk_name
  for (i in 1:length(chunk_name)) {

    selected_barcodes <- sample(barcodes_sample,length(chunk_bins[[chunk_name[i]]]))

    barcodes_ls[chunk_name[i]] <- list(selected_barcodes)

    barcodes_sample <- barcodes_sample[!barcodes_sample %in% selected_barcodes]

  }

  return(barcodes_ls)

}

#' generate expression symbol matrix
#'
#' @param barcodes the barcodes of a chunk
#' @param raw_count_matrix the raw count matrix
#' @param gene_ls the gene list from expression symbol matrix

.generate_expression_symbol_matrix <- function(barcodes,raw_count_matrix,gene_ls) {

  on.exit(gc())

  expr_sym_matrix <- raw_count_matrix[gene_ls,barcodes]

  expr_sym_matrix[expr_sym_matrix != 0] <- 1
  expr_sym_matrix[expr_sym_matrix == 0] <- -1

  return(expr_sym_matrix)

}

#' calculate chunk characteristic value
#'
#' @param seu_meta the Seurat Object metadata
#' @param spatial_chunks the chunks for calculation
#' @param kernel_function choose the kernel function
#' @param raw_count_matrix the raw count matrix
#' @param gene_ls the gene list from expression symbol matrix

.calculate_chunk_characteristic_value <- function(seu_meta,raw_count_matrix,spatial_chunks,gene_ls,kernel_function = "Gaussian") {

  on.exit(gc())

  dis_mat <- ICHMousewch:::.genetate_distance_matrix(seu_meta = seu_meta,
                                                     spatial_chunks = spatial_chunks,
                                                     kernel_function = kernel_function)

  dis_mat <- Matrix::Matrix(dis_mat@.Data,sparse = TRUE)

  expr_sym_mat <- ICHMousewch:::.generate_expression_symbol_matrix(barcodes = spatial_chunks,
                                                                   raw_count_matrix = raw_count_matrix,
                                                                   gene_ls = gene_ls)

  cha_val_mat_ls <- vector("list",length = length(gene_ls))
  names(cha_val_mat_ls) <- gene_ls

  for (i in 1:length(gene_ls)) {

    cha_val_mat <- dis_mat*expr_sym_mat[gene_ls[i],]

    cha_val_mat_ls[gene_ls[i]] <- cha_val_mat %>%
      Matrix::colSums() %>%
      scale() %>%
      list()

  }

  cha_val_mat_ls <- cha_val_mat_ls %>%
    as.data.frame()

  return(cha_val_mat_ls)

}

