
# R/utils-general.R

#' integrate file address
#'
#' @param raw_count_matrix_address address for raw count matrix
#' @param filtered_count_matrix address for filtered count matrix
#' @param tissue_position_address address for tissue position matrix
#' @param background_image_address address for background address

.integrate_file_address <- function(raw_count_matrix_address,filtered_count_matrix_address,tissue_position_address,background_image_address) {

  on.exit(gc())

  file_address <- c("raw_count_matrix_address" = raw_count_matrix_address,
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

  return(as.data.table(tissue_position_matrix))

}

#' load raw count matrix
#'
#' @param raw_count_matrix_address the file address of raw count matrix

.load_raw_count_matrix <- function(raw_count_matrix_address) {

  on.exit(gc())

  raw_count_matrix <- Read10X_h5(raw_count_matrix_address)

  return(raw_count_matrix)

}

#' find filtered genes
#'
#' @param raw_count_matrix the matrix of raw count dataset
#' @param original_seu_metadata the original Seurat Object metadata

.find_filtered_genes <- function(raw_count_matrix,original_seu_metadata) {

  on.exit(gc())

  filtered_count_matrix <- raw_count_matrix[,original_seu_metadata[tissue == 1,barcode]]

  num_gene <- rownames(filtered_count_matrix) %>%
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

#' conduct statistic test on rate
#'
#' @param sample_dt the sample for analysis

.conduct_statistic_test_on_rate <- function(sample_dt) {

  on.exit(gc())

  sample_dt[,log_diff_pct := log2(pct.1/pct.2)]

  p_value_ls <- list()
  conf_up_ls <- list()
  conf_down_ls <- list()
  for (i in 1:length(rownames(sample_dt))) {

    test_result <- prop.test(x = sample_dt[i,c(posi.1,posi.2)],
                             n = sample_dt[i,c(total.1,total.2)])

    p_value_ls <- append(p_value_ls,list(test_result$p.value))
    conf_up_ls <- append(conf_up_ls,list(test_result$conf.int[1]))
    conf_down_ls <- append(conf_down_ls,list(test_result$conf.int[2]))

  }

  sample_dt[,p_value_rate := unlist(p_value_ls)]
  sample_dt[,conf_up_rate := unlist(conf_up_ls)]
  sample_dt[,conf_down_rate := unlist(conf_down_ls)]

  return(sample_dt)

}

#' save spatial image
#'
#' @param saving_path the path for saving
#' @param spatial_image the class of Spatial Image

.save_spatial_image <- function(saving_path,spatial_image,width,height) {

  on.exit(gc())

  file_path <- paste(saving_path,spatial_image@image_set_name,sep = "/")

  if(!dir.exists(file_path)) {

    dir.create(file_path,recursive = TRUE)

  } else {

    unlink(file_path,recursive = TRUE)
    dir.create(file_path,recursive = TRUE)

  }

num_image <- length(spatial_image@spatial_image)

  for (i in 1:num_image) {

    file_name <- paste(names(spatial_image@spatial_image)[i],"png",sep = ".")

    image_plot <- as.ggplot(spatial_image@spatial_image[[i]])

    ggsave(filename = file_name,
           plot = image_plot,
           path = file_path,
           width = width,
           height = height,
           units = "in",
           dpi = 300)

  }

}

#' filter differential expression genes
#'
#' @param diff_expr_genes the differential expression genes

.filter_diff_expr_genes <- function(diff_expr_genes) {

  on.exit(gc())

  diff_expr_genes <- diff_expr_genes[abs(avg_log2FC) > 2]

  setorder(diff_expr_genes,-log2pct)

  return(diff_expr_genes)

}

#' annotate the cell type based on single gene
#'
#' @param ich_mouse the ICH_Mouse class
#' @param gene_ls the gene ls used to annotate

.annotate_the_cell_type_based_on_single_gene <- function(ich_mouse,gene_ls) {

  on.exit(gc())

  ref_ds <- ICHMousewch::mouse_cell_variable_genes

  cell_type <- names(ref_ds)

  annotation_ls <- vector("list",length = length(gene_ls))
  names(annotation_ls) <- gene_ls
  for (i in 1:length(gene_ls)) {

    single_gene <- gene_ls[i]

    cell_type_ls <- list()
    for (j in 1:length(cell_type)) {

      same_gene <- single_gene  %in% ref_ds[[cell_type[j]]]

      if(same_gene) {

        cell_type_ls <- append(cell_type_ls,cell_type[j])

      }

    }

    if(length(cell_type_ls) == 0) {

      annotation_ls[single_gene] <- "no value"

    } else {

      annotation_ls[single_gene] <- list(unlist(cell_type_ls))

    }

  }

  annotation_dt <- data.table()
  annotation_dt[,gene_name := names(annotation_ls)]
  annotation_dt[,cell_type := annotation_ls]

  return(annotation_dt)

}

#' refresh internal dataset
#'
#' @export

refresh_internal_dataset <- function() {

  on.exit(gc())

  ICHMousewch:::.download_internal_dataset()

}











