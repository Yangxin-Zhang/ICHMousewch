
# R/utils-ANOVA.R

#' find symbol genes according to cluster
#'
#' @param raw_count_matrix the matrix of raw count dataset
#' @param seu_metadata_with_cluster_symbol the Seurat Object metadata with cluster symbol
#' @param filtered_genes the gene that been filtered
#' @param cluster_symbol the symbol of cluster

.find_symbol_genes <- function(raw_count_matrix,seu_metadata_with_cluster_symbol,filtered_genes,cluster_symbol) {

  on.exit(gc())

  seu_obj <- CreateSeuratObject(raw_count_matrix[filtered_genes,seu_metadata_with_cluster_symbol[center_edge_symbol %in% cluster_symbol,center_edge_symbol]]) %>%
    NormalizeData(scale.factor = 1000000) %>%
    FindVariableFeatures(selection.method = "vst")

  assay_metadata <- as.data.table(seu_obj@assay$RNA@meta.data)

  assay_metadata[,gene_name := rownames(seu_obj)]

  symbol_genes <- assay_metadata[abs(vf_vst_counts_variance.standardized) >= 2,gene_name]

  return(symbol_genes)

}
