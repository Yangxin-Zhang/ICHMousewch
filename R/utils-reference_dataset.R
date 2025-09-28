
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

#' find variable gene
#' @param reference_dataset the reference dataset

.variable_genes <- function(reference_dataset) {

  on.exit(gc())

  cell_type <- colnames(reference_dataset)[!colnames(reference_dataset) %in% "gene_name"]

  working_path <- "/home/youngxin/Desktop/saving_directory"

  if(!dir.exists(working_path)) {

    dir.create(working_path,recursive = TRUE)

  }

##
  file_name <- paste(working_path,"var_gene_ls.rds",sep = "/")
  if (file.exists(file_name)) {

    var_gene_ls <- readRDS(file = file_name)

    sub_cell_type <- list()
    for (i in 1:length(var_gene_ls)) {

      if (is.null(var_gene_ls[[i]])) {

        sub_cell_type <- append(sub_cell_type,list(names(var_gene_ls)[i]))

      }

    }

    sub_cell_type <- unlist(sub_cell_type)

    for (i in 1:length(sub_cell_type)) {

      contrasting_cell_type <- sub_cell_type[!sub_cell_type %in% sub_cell_type[i]]

      gene_ls <- vector("list",length = length(contrasting_cell_type))
      names(gene_ls) <- contrasting_cell_type
      for (j in 1:length(contrasting_cell_type)) {

        aimming_cell_type <- c(sub_cell_type[i],contrasting_cell_type[j],"gene_name")
        ref_dt <- reference_dataset[,..aimming_cell_type]

        cell_type1 <- sub_cell_type[i]
        cell_type1 <- ref_dt[,..cell_type1] %>%
          unlist() %>%
          as.numeric()

        cell_type2 <- contrasting_cell_type[j]
        cell_type2 <- ref_dt[,..cell_type2] %>%
          unlist() %>%
          as.numeric()

        ref_dt[,cell_type1 := cell_type1]
        ref_dt[,cell_type2 := cell_type2]

        ref_dt[,avglogFC := cell_type1 - cell_type2]

        var_genes <- ref_dt[abs(avglogFC) > 2,gene_name]
        gene_ls[contrasting_cell_type[j]] <- list(var_genes)

      }

      var_gene_ls[sub_cell_type[i]] <- list(gene_ls)

      saveRDS(var_gene_ls,file = file_name)

    }

  } else {

    var_gene_ls <- vector("list",length = length(cell_type))
    names(var_gene_ls) <- cell_type
    for (i in 1:length(cell_type)) {

      contrasting_cell_type <- cell_type[!cell_type %in% cell_type[i]]

      gene_ls <- vector("list",length = length(contrasting_cell_type))
      names(gene_ls) <- contrasting_cell_type
      for (j in 1:length(contrasting_cell_type)) {

        aimming_cell_type <- c(cell_type[i],contrasting_cell_type[j],"gene_name")
        ref_dt <- reference_dataset[,..aimming_cell_type]

        cell_type1 <- cell_type[i]
        cell_type1 <- ref_dt[,..cell_type1] %>%
          unlist() %>%
          as.numeric()

        cell_type2 <- contrasting_cell_type[j]
        cell_type2 <- ref_dt[,..cell_type2] %>%
          unlist() %>%
          as.numeric()

        ref_dt[,cell_type1 := cell_type1]
        ref_dt[,cell_type2 := cell_type2]

        ref_dt[,avglogFC := cell_type1 - cell_type2]

        var_genes <- ref_dt[abs(avglogFC) > 2,gene_name]
        gene_ls[contrasting_cell_type[j]] <- list(var_genes)

      }

      var_gene_ls[cell_type[i]] <- list(gene_ls)

      saveRDS(var_gene_ls,file = file_name)

    }

  }
##

##
  file_name <- paste(working_path,"symbol_genes.rds",sep = "/")
  if (file.exists(file_name)) {

    symbol_genes <- readRDS(file = file_name)

    sub_cell_type <- list()
    for (i in 1:length(symbol_genes)) {

      if (is.null(symbol_genes[[i]])) {

        sub_cell_type <- append(sub_cell_type,list(names(symbol_genes)[i]))

      }

    }

    sub_cell_type <- unlist(sub_cell_type)

    for (i in 1:length(sub_cell_type)) {

      genes <- Reduce(union,var_gene_ls[[sub_cell_type[i]]])
      totoal_gene <- unlist(var_gene_ls[[sub_cell_type[i]]])

      rate_ls <- vector("list",length = length(genes))
      names(rate_ls) <- genes
      for (j in 1:length(genes)) {

        num <- sum(totoal_gene %in% genes[j])
        rate <- num/length(var_gene_ls[[sub_cell_type[i]]])
        rate_ls[genes[j]] <- list(rate)

      }


      results <- data.table(present_rate = unlist(rate_ls),
                            gene_name = names(rate_ls))
      setorder(results,-present_rate)
      results <- results[present_rate > 0.8,gene_name]

      symbol_genes[sub_cell_type[i]] <- list(results)

      saveRDS(symbol_genes,file = file_name,compress = FALSE)

    }

  } else {

    symbol_genes <- vector("list",length = length(cell_type))
    names(symbol_genes) <- cell_type
    for (i in 1:length(cell_type)) {

      genes <- Reduce(union,var_gene_ls[[cell_type[i]]])
      totoal_gene <- unlist(var_gene_ls[[cell_type[i]]])

      rate_ls <- vector("list",length = length(genes))
      names(rate_ls) <- genes
      for (j in 1:length(genes)) {

        num <- sum(totoal_gene %in% genes[j])
        rate <- num/length(var_gene_ls[[cell_type[i]]])
        rate_ls[genes[j]] <- list(rate)
        }

      results <- data.table(present_rate = unlist(rate_ls),
                            gene_name = names(rate_ls))
      setorder(results,-present_rate)
      results <- results[present_rate > 0.8,gene_name]

      symbol_genes[cell_type[i]] <- list(results)

      saveRDS(symbol_genes,file = file_name,compress = FALSE)
      }

  }
##

  return(symbol_genes)

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

  gene_ls <- reference_dataset[,gene_name]
  contrasting_genes <- gene_ls[gene_ls %in% ICHMousewch::integrated_mouse_RNA_seq_dataset[,gene_name]]

  ref_dt <- reference_dataset[gene_name %in% contrasting_genes]

  cell_type <- colnames(ref_dt)[!(colnames(ref_dt) == "gene_name")]

  other_type_dt <- ICHMousewch::integrated_mouse_RNA_seq_dataset
  other_type <- colnames(other_type_dt)
  other_type <- other_type[!other_type %in% "gene_name"]

  other_type_dt <- other_type_dt[gene_name %in% contrasting_genes,..other_type]

  other_type_dt <- rowMeans(other_type_dt)

  result_dt <- data.table()
  result_dt[,gene_name := ref_dt[,gene_name]]
  for (i in 1:length(cell_type)) {

    one_type <- cell_type[i]
    one_type_dt <- ref_dt[,..one_type] %>%
      unlist() %>%
      as.numeric()

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

  mouse_cell_variable_genes <- ICHMousewch:::.variable_genes(reference_dataset = integrated_mouse_RNA_seq_dataset)
  mouse_immune_cell_variable_genes <- ICHMousewch:::.variable_genes(reference_dataset = integrated_immune_mouse_RNA_seq_dataset)

  usethis::use_data(gene_id_information,
                    integrated_mouse_RNA_seq_dataset,
                    integrated_immune_mouse_RNA_seq_dataset,
                    mouse_cell_variable_genes,
                    mouse_immune_cell_variable_genes,
                    overwrite = TRUE,
                    internal = FALSE)

}

