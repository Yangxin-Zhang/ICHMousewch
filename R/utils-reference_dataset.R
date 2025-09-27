
# R/utils-reference_dataset.R

#' integrate mouse RNA seq dataset based on cell type
#'

.integrate_mouse_RNA_seq_dataset <- function() {

  on.exit(gc())

  ref_da <- celldex::MouseRNAseqData()

  anno_data <- as.data.table(colData(ref_da))
  anno_data[,assay_label := colData(ref_da)@rownames]

  cell_type <- anno_data[,label.fine] %>%
    unique()

  integrated_expr_matrix <- data.table()
  integrated_expr_matrix[,gene_name := rownames(assay(ref_da))]
  for (i in 1:length(cell_type)) {

    ass_lb <- anno_data[label.fine == cell_type[i],assay_label]

    if (length(ass_lb) == 1) {

      integrated_expr_matrix[,cell_type[i] := assay(ref_da)[,ass_lb]]

    } else {

      expr_matrix <- rowMeans(assay(ref_da)[,ass_lb])

      integrated_expr_matrix[,cell_type[i] := expr_matrix]

    }
  }

  return(integrated_expr_matrix)

}

#' integrate immune mouse RNA seq dataset based on cell type
#'

.integrate_immune_mouse_RNA_seq_dataset <- function() {

  on.exit(gc())

  ref_da <- celldex::ImmGenData()

  anno_data <- as.data.table(colData(ref_da))
  anno_data[,assay_label := colData(ref_da)@rownames]

  cell_type <- anno_data[,label.fine] %>%
    unique()

  integrated_expr_matrix <- data.table()
  integrated_expr_matrix[,gene_name := rownames(assay(ref_da))]
  for (i in 1:length(cell_type)) {

    ass_lb <- anno_data[label.fine == cell_type[i],assay_label]

    if (length(ass_lb) == 1) {

      integrated_expr_matrix[,cell_type[i] := assay(ref_da)[,ass_lb]]

    } else {

      expr_matrix <- rowMeans(assay(ref_da)[,ass_lb])

      integrated_expr_matrix[,cell_type[i] := expr_matrix]

    }
  }

  return(integrated_expr_matrix)

}

#' transfer reference dataset to matrix
#'
#' @param ref_dt the reference dataset in the form of data.table

.transfer_ref_dt_to_mat <- function(ref_dt) {

  on.exit(gc())

  cell_type <- colnames(ref_dt)[!(colnames(ref_dt) == "gene_name")]

  ref_mat <- as.matrix(ref_dt[,..cell_type])
  rownames(ref_mat) <- ref_dt[,gene_name]

  return(ref_mat)

}

#' calculate the log2FC one type to other type
#'
#' @param reference_dataset the reference dataset

.the_log2FC_one_type_to_other_type <- function(reference_dataset) {

  on.exit(gc())

  ref_dt <- reference_dataset

  cell_type <- colnames(ref_dt)[!(colnames(ref_dt) == "gene_name")]

  result_dt <- data.table()
  result_dt[,gene_name := ref_dt[,gene_name]]
  for (i in 1:length(cell_type)) {

    one_type <- cell_type[i]
    one_type_dt <- ref_dt[,..one_type] %>%
      unlist() %>%
      as.numeric()

    other_type <- cell_type[!(cell_type == cell_type[i])]
    other_type_dt <- ref_dt[,..other_type]

    other_type_dt <- rowMeans(other_type_dt)

    result_dt[,cell_type[i] := log2(one_type_dt/other_type_dt)]

  }

  return(result_dt)

}

#' calculate the log2FC one type to one type
#'
#' @param reference_dataset the reference dataset

.the_log2FC_one_type_to_one_type <- function(reference_dataset) {

  on.exit(gc())

  ref_dt <- reference_dataset

  cell_type <- colnames(ref_dt)[!(colnames(ref_dt) == "gene_name")]

  result_ls <- vector("list",length = length(cell_type))
  names(result_ls) <- cell_type

  for (i in 1:length(cell_type)) {

    result_dt <- data.table()
    result_dt[,gene_name := ref_dt[,gene_name]]

    col_i <- cell_type[i]
    data_i <- ref_dt[,..col_i] %>%
      unlist() %>%
      as.numeric()

    for (j in 1:length(cell_type)) {

      if (cell_type[j] != cell_type[i]) {

        col_j <- cell_type[j]
        data_j <- ref_dt[,..col_j] %>%
          unlist() %>%
          as.numeric()

        result_dt[,paste(cell_type[i],cell_type[j],sep = "-") := log2(data_i/data_j)]

      } else {

        result_dt[,paste(cell_type[i],cell_type[j],sep = "-") := 0]

      }

    }

    result_ls[cell_type[i]] <- list(result_dt)

  }

  return(result_ls)

}

#' get differential exoression genes from log2FC data.table
#' @param log2FC_dt the log2FC data.table

.get_de_genes_from_log2FC_datatable <- function(log2FC_dt) {

  on.exit(gc())

  cell_type <- colnames(log2FC_dt)[!(colnames(log2FC_dt) == "gene_name")]

  de_gene_ls <- vector("list",length = length(cell_type))
  names(de_gene_ls) <- cell_type
  for (i in 1:length(cell_type)) {

    cell_sym <- cell_type[i]
    dt_col <- c(cell_sym,"gene_name")
    dt <- log2FC_dt[,..dt_col]

    utils_col <- dt[,..cell_sym] %>%
      unlist() %>%
      as.numeric()

    dt[,utils_col := utils_col]

    setorderv(dt,cols = cell_type[i],order = -1)

    de_gene <- dt[utils_col > 1,gene_name]

    de_gene_ls[cell_type[i]] <- list(de_gene)

  }

  return(de_gene_ls)

}

#' differential exoression genes one type to other type according to mouse RNA seq

.diff_expr_genes_cell_type_one_to_other <- function() {

  on.exit(gc())

  mouse_RNA_seq_dataset <- ICHMousewch::integrated_mouse_RNA_seq_dataset

  log2FC_dt <- ICHMousewch:::.the_log2FC_one_type_to_other_type(reference_dataset = mouse_RNA_seq_dataset)

  results <- ICHMousewch:::.get_de_genes_from_log2FC_datatable(log2FC_dt = log2FC_dt)

  return(results)

}

#' differential exoression genes one type to other type according to immune mouse RNA seq

.diff_expr_genes_immune_cell_type_one_to_other <- function() {

  on.exit(gc())

  mouse_RNA_seq_dataset <- ICHMousewch::integrated_immune_mouse_RNA_seq_dataset

  log2FC_dt <- ICHMousewch:::.the_log2FC_one_type_to_other_type(reference_dataset = mouse_RNA_seq_dataset)

  results <- ICHMousewch:::.get_de_genes_from_log2FC_datatable(log2FC_dt = log2FC_dt)

  return(results)

}

#' differential exoression genes one type to one type

.diff_expr_genes_one_type_to_one_type <- function() {

  on.exit(gc())

  log2FC_ls <- ICHMousewch:::.the_log2FC_one_type_to_one_type()

  cell_group <- names(log2FC_ls)

  results <- vector("list",length = length(cell_group))
  names(results) <- cell_group
  for (i in 1:length(cell_group)) {

    log2FC_dt <- log2FC_ls[[cell_group[i]]]

    de_gene_ls <- ICHMousewch:::.get_de_genes_from_log2FC_datatable(log2FC_dt = log2FC_dt)

    results[cell_group[i]] <- list(de_gene_ls)

  }

  return(results)

}

#' download gene id information

.download_gene_id_information <- function() {

  on.exit(gc())

  ensembl_mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

  mouse_genes <- union(celldex::MouseRNAseqData()@NAMES,celldex::ImmGenData()@NAMES)

  gene_id_information <- getBM(attributes = c("mgi_symbol", "ensembl_gene_id","external_gene_name","gene_biotype","entrezgene_id"),
                               filters = "mgi_symbol",
                               values = mouse_genes,
                               mart = ensembl_mouse)

  gene_id_information <- as.data.table(gene_id_information)

  return(gene_id_information)

}
#' download internal dataset

.download_internal_dataset <- function() {

  on.exit(gc())

  gene_id_information <- ICHMousewch:::.download_gene_id_information()
  integrated_mouse_RNA_seq_dataset <- ICHMousewch:::.integrate_mouse_RNA_seq_dataset()
  integrated_immune_mouse_RNA_seq_dataset <- ICHMousewch:::.integrate_immune_mouse_RNA_seq_dataset()

  de_genes_cell_type_one_to_other <- ICHMousewch:::.diff_expr_genes_cell_type_one_to_other()
  de_genes_immune_cell_type_one_to_other <- ICHMousewch:::.diff_expr_genes_immune_cell_type_one_to_other()

  usethis::use_data(gene_id_information,
                    integrated_mouse_RNA_seq_dataset,
                    integrated_immune_mouse_RNA_seq_dataset,
                    de_genes_cell_type_one_to_other,
                    de_genes_immune_cell_type_one_to_other,
                    overwrite = TRUE,
                    internal = FALSE)

}

