
# R/utils-seu_obj.R

#' get the Seurat Object metadata
#'
#' @param seu_obj A Seurat Object

.get_seu_obj_metadata <- function(seu_obj) {

  on.exit(gc())

  return(seu_obj@meta.data)

}

#' generate original Seurat Object metadata
#'
#' @param raw_count_matrix the count  matrix of raw dataset
#' @param tissue_position_matrix the matrix of tissue position

.generate_original_seu_metadata <- function(raw_count_matrix,tissue_position_matrix) {

  on.exit(gc())

  seu_obj <- CreateSeuratObject(raw_count_matrix)

  # add barcode to Seurat Object metadata
  seu_obj@meta.data$barcode <- rownames(seu_obj@meta.data)

  # calculate mito and ribo percent from Seurat Object
  seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj,pattern = "^mt-")
  seu_obj[["percent.ribo"]] <- PercentageFeatureSet(seu_obj,pattern = "Rp[sl]")

  # extract Seurat Object metadata
  seu_obj_metadata <- as.data.table(seu_obj@meta.data)

  # calculate log2 value of nCount and nFeature
  seu_obj_metadata[,nCount_log2 := log2(seu_obj_metadata[,nCount_RNA]+1)]

  seu_obj_metadata[,nFeature_log2 := log2(seu_obj_metadata[,nFeature_RNA]+1)]

  # add tissue position information to Seurat Object metadata
  adj_tissue_position_matrix <- tissue_position_matrix[tissue_position_matrix[,barcode] %in% seu_obj_metadata[,barcode]]
  barcode_order <- match(seu_obj_metadata[,barcode],adj_tissue_position_matrix[,barcode])
  adj_tissue_position_matrix <- adj_tissue_position_matrix[barcode_order]

  seu_obj_metadata[,imagerow := adj_tissue_position_matrix[,imagerow]]
  seu_obj_metadata[,imagecol := adj_tissue_position_matrix[,imagecol]]
  seu_obj_metadata[,tissue := adj_tissue_position_matrix[,tissue]]
  seu_obj_metadata[,row := adj_tissue_position_matrix[,row]]
  seu_obj_metadata[,col := adj_tissue_position_matrix[,col]]

  return(seu_obj_metadata)

}

#' generate Seurat Object metadata with cluster symbol
#'
#' @param original_seu_metadata original Seurat Object metadata
#' @param raw_count_matrix the matrix of raw count dataset

.generate_seu_metadata_with_cluster_symbol <- function(original_seu_metadata,in_tissue_barcode,raw_count_matrix) {

  on.exit(gc())

  metadata_with_cluster_symbol <- ICHMousewch:::.conduct_Louvain_cluster_based_on_tissue_position(original_seu_metadata = original_seu_metadata,
                                                                                                  raw_count_matrix = raw_count_matrix)

  return(metadata_with_cluster_symbol)

}

#' add hematoma symbol to Seurat Object metadata
#'
#' @param hematoma the Hematoma class

.generate_seu_metadata_with_hematoma_symbol <- function(hematoma) {

  on.exit(gc())

  # add hematoma symbol
  hematoma_symbols <- hematoma@identification_symbols$hematoma_symbols

  seu_metadata <- hematoma@seu_metadata_with_cluster_symbol

  Louvain_cluster_posi <- seu_metadata$Louvain_cluster_posi %>%
    unique() %>%
    unlist()

  non_hematoma_symbols <- Louvain_cluster_posi[!Louvain_cluster_posi %in% hematoma_symbols]

  seu_metadata$hematoma_symbol <- seu_metadata$Louvain_cluster_posi

  for (i in 1:length(hematoma_symbols)) {

    seu_metadata[seu_metadata$Louvain_cluster_posi == hematoma_symbols[i],]$hematoma_symbol <- 2

  }

  for (i in 1:length(non_hematoma_symbols)) {

    seu_metadata[seu_metadata$Louvain_cluster_posi == non_hematoma_symbols[i],]$hematoma_symbol <- 1

  }

  # add hematoma center and edge symbol
  filtered_genes <- hematoma@filtered_genes
  raw_count_matrix <- hematoma@raw_count_matrix

  metadata_hematoma <- seu_metadata[seu_metadata$hematoma_symbol == 2,]
  metadata_normal_tissue <- seu_metadata[seu_metadata$hematoma_symbol == 1,]

  seu_obj_hematoma <- CreateSeuratObject(raw_count_matrix[filtered_genes,metadata_hematoma$barcode]) %>%
    NormalizeData(scale.factor = 1000000) %>%
    ScaleData() %>%
    RunPCA(features = filtered_genes,
           npcs = 100) %>%
    FindNeighbors(reduction = "pca",
                  dims = 1:10) %>%
    FindClusters(algorithm = 1,
                 resolution = 0.5,
                 random.seed = 2025,
                 cluster.name = "Louvain_cluster_filt_gene")

  metadata_normal_tissue$Louvain_cluster_filt_gene <- 1

  diff_value <- seu_obj_hematoma@meta.data$Louvain_cluster_filt_gene %>%
    as.numeric() %>%
    unique() %>%
    unlist() %>%
    min()

  metadata_hematoma$Louvain_cluster_filt_gene <- as.numeric(seu_obj_hematoma@meta.data$Louvain_cluster_filt_gene) + 2 - diff_value
  metadata_hematoma <- metadata_hematoma[,colnames(metadata_normal_tissue)]

  metadata_hematoma_filt_gene <- bind_rows(metadata_hematoma,metadata_normal_tissue)
  metadata_order <- match(seu_metadata$barcode,metadata_hematoma_filt_gene$barcode)

  return(metadata_hematoma_filt_gene[metadata_order,])

}

#' add hematoma center and edge symbols to Seurat Object metadata
#'
#' @param hematoma the Hematoma class

.generate_seu_metadata_with_hematoma_center_and_edge_symbol <- function(hematoma) {

  on.exit(gc())

  seu_metadata <- hematoma@seu_metadata_with_cluster_symbol

  center_symbols <- hematoma@identification_symbols$center_symbols
  edge_symbols <- hematoma@identification_symbols$edge_symbols

  seu_metadata$center_edge_symbol <- seu_metadata$Louvain_cluster_filt_gene

  for (i in 1:length(center_symbols)) {

    seu_metadata[seu_metadata$Louvain_cluster_filt_gene == center_symbols[i],]$center_edge_symbol <- 2

  }

  for (i in 1:length(edge_symbols)) {

    seu_metadata[seu_metadata$Louvain_cluster_filt_gene == edge_symbols[i],]$center_edge_symbol <- 3

  }

  return(seu_metadata)

}

