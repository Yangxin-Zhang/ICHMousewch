
# R/utils-import.R

#' @import org.Mm.eg.db
#' @importFrom Seurat Read10X_Coordinates
#' @importFrom Seurat Read10X_h5
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat ScaleData
#' @importFrom Seurat RunPCA
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat FindMarkers
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom Giotto createGiottoInstructions
#' @importFrom Giotto createGiottoObject
#' @importFrom Giotto createGiottoImage
#' @importFrom Giotto addGiottoImage
#' @importFrom Giotto spatPlot2D
#' @importFrom Giotto spatFeatPlot2D
#' @importFrom ggplot2 ggplotGrob
#' @importFrom ggplot2 ggsave
#' @importFrom ggplotify as.ggplot
#' @importFrom patchwork patchworkGrob
#' @importFrom patchwork wrap_plots
#' @importFrom dplyr bind_rows
#' @importFrom dplyr union
#' @importFrom dplyr intersect
#' @importFrom Matrix rowSums
#' @importFrom Matrix sparseMatrix
#' @importFrom data.table as.data.table
#' @importFrom data.table rbindlist
#' @importFrom data.table data.table
#' @importFrom data.table setorder
#' @importFrom data.table setorderv
#' @importFrom celldex MouseRNAseqData
#' @importFrom celldex ImmGenData
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assay
#' @importFrom biomaRt useEnsembl
#' @importFrom biomaRt getBM
#' @importFrom usethis use_data
#' @importFrom clusterProfiler enrichGO
#' @importFrom clusterProfiler gseKEGG
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom clusterProfiler simplify
#' @importFrom simplifyEnrichment GO_similarity

NULL
