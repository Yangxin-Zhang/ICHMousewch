
# R/utils-ANOVA.R

#' find symbol genes according to cluster
#'
#' @param raw_count_matrix the matrix of raw count dataset
#' @param seu_metadata_with_cluster_symbol the Seurat Object metadata with cluster symbol
#' @param filtered_genes the gene that been filtered
#' @param cluster_symbol the symbol of cluster

.find_symbol_genes <- function(raw_count_matrix,seu_metadata_with_cluster_symbol,filtered_genes,cluster_symbol) {

  on.exit(gc())

  seu_obj <- CreateSeuratObject(raw_count_matrix[filtered_genes,seu_metadata_with_cluster_symbol[center_edge_symbol %in% cluster_symbol,barcode]]) %>%
    NormalizeData(scale.factor = 1000000) %>%
    FindVariableFeatures(selection.method = "vst")

  assay_metadata <- as.data.table(seu_obj@assays$RNA@meta.data)

  assay_metadata[,gene_name := rownames(seu_obj)]

  symbol_genes <- assay_metadata[abs(vf_vst_counts_variance.standardized) >= 2,gene_name]

  return(symbol_genes)

}

#' find differential expression genes
#'
#' @param raw_count_matrix the matrix of raw count dataset
#' @param seu_metadata_with_cluster_symbol the Seurat Object metadata with cluster symbol
#' @param filtered_genes the gene that been filtered
#' @param cluster_symbol the symbol of cluster

.find_differential_expression_genes <- function(raw_count_matrix,seu_metadata_with_cluster_symbol,filtered_genes,cluster_symbol) {

  on.exit(gc())

  seu_obj <- CreateSeuratObject(raw_count_matrix[filtered_genes,seu_metadata_with_cluster_symbol[,barcode]]) %>%
    NormalizeData(scale.factor = 1000000)

  seu_obj@meta.data$barcode <- rownames(seu_obj@meta.data)
  order_metadata <- match(seu_obj@meta.data$barcode,seu_metadata_with_cluster_symbol[,barcode])
  seu_obj@meta.data$center_edge_symbol <- seu_metadata_with_cluster_symbol[order_metadata,center_edge_symbol]

  diff_expr_genes_dt <- seu_obj %>%
    FindMarkers(ident.1 = cluster_symbol[1],
                ident.2 = cluster_symbol[2],
                group.by = "center_edge_symbol")

  diff_expr_genes_dt$gene_name <- rownames(diff_expr_genes_dt)

  diff_expr_genes_dt <- as.data.table(diff_expr_genes_dt) %>%
    setorder(-avg_log2FC)

  count_matrix_1 <- raw_count_matrix[diff_expr_genes_dt[,gene_name],seu_metadata_with_cluster_symbol[center_edge_symbol == cluster_symbol[1],barcode]]
  count_matrix_2 <- raw_count_matrix[diff_expr_genes_dt[,gene_name],seu_metadata_with_cluster_symbol[center_edge_symbol == cluster_symbol[2],barcode]]

  diff_expr_genes_dt[,total.1 := ncol(count_matrix_1)]
  diff_expr_genes_dt[,total.2 := ncol(count_matrix_2)]

  bins <- seq(1,nrow(diff_expr_genes_dt),2000)
  num_posi.1_ls <- list()
  num_posi.2_ls <- list()
  for (i in 1:length(bins)) {

    if (i != length(bins)) {

      num_posi.1 <- Matrix::rowSums(!count_matrix_1[bins[i]:(bins[i]+1999),] == 0)
      num_posi.2 <- Matrix::rowSums(!count_matrix_2[bins[i]:(bins[i]+1999),] == 0)

    } else {

      num_posi.1 <- Matrix::rowSums(!count_matrix_1[bins[i]:nrow(diff_expr_genes_dt),] == 0)
      num_posi.2 <- Matrix::rowSums(!count_matrix_2[bins[i]:nrow(diff_expr_genes_dt),] == 0)

    }

    num_posi.1_ls <- append(num_posi.1_ls,list(num_posi.1))
    num_posi.2_ls <- append(num_posi.2_ls,list(num_posi.2))

    gc()

  }

  diff_expr_genes_dt[,posi.1 := unlist(num_posi.1_ls)]
  diff_expr_genes_dt[,posi.2 := unlist(num_posi.2_ls)]

  diff_expr_genes_dt <- ICHMousewch:::.conduct_statistic_test_on_rate(sample_dt = diff_expr_genes_dt)

  return(diff_expr_genes_dt)

}
