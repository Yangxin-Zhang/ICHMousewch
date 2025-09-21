
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

}
