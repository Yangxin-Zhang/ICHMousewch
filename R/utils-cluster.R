
# R/utils-cluster.R

#' conduct GMM classification
#'
#' @param original_seu_metadata original Seurat Object metadata

.conduct_GMM_classification <- function(original_seu_metadata,hematoma_larger_than_normal_tissue = TRUE) {

  on.exit(gc())

  set.seed(2025)

  in_tissue_metadata <- original_seu_metadata[tissue == 1]

  GMM_result <- Mclust(data = in_tissue_metadata[,c("nCount_log2","nFeature_log2")],
                       G = 2)

  in_tissue_metadata[,GMM_cluster := GMM_result$classification]

  num_1 <- sum(in_tissue_metadata[GMM_cluster == 1,GMM_cluster])
  num_2 <- sum(in_tissue_metadata[GMM_cluster == 2,GMM_cluster])

  if(num_1 < num_2 && hematoma_larger_than_normal_tissue) {

    in_tissue_metadata[GMM_cluster == 1,GMM_cluster := 0]
    in_tissue_metadata[GMM_cluster == 2,GMM_cluster := 1]
    in_tissue_metadata[GMM_cluster == 0,GMM_cluster := 2]

  }

  return(in_tissue_metadata)

}

#' conduct Louvain cluster based on tissue position matrix
#'
#' @param original_seu_metadata original Seurat Object metadata
#' @param raw_count_matrix the matrix of raw count dataset

.conduct_Louvain_cluster_based_on_tissue_position <- function(original_seu_metadata,raw_count_matrix) {

  on.exit(gc())

  metadata_with_GMM <- ICHMousewch:::.conduct_GMM_classification(original_seu_metadata = original_seu_metadata)

  metadata_with_GMM_1 <- metadata_with_GMM[GMM_cluster == 1]

  metadata_with_GMM_2 <- metadata_with_GMM[GMM_cluster == 2]

  # create Seurat Object according to GMM cluster
  raw_count_matrix_GMM_2 <- raw_count_matrix[,metadata_with_GMM_2[,barcode]]
  seu_obj_GMM_2 <- CreateSeuratObject(raw_count_matrix_GMM_2)

  # conduct Louvain algorithm
  metadata_for_KNN <- metadata_with_GMM_2[,c("imagerow","imagecol")] %>%
    as.matrix()
  rownames(metadata_for_KNN) <- metadata_with_GMM_2[,barcode]

  KNN_graph <- metadata_for_KNN %>%
    FindNeighbors(k.param = 100)

  seu_obj_GMM_2@graphs$position_nn <- KNN_graph[["nn"]]
  seu_obj_GMM_2@graphs$position_snn <- KNN_graph[["snn"]]

  seu_obj_GMM_2 <- FindClusters(object = seu_obj_GMM_2,
                                resolution = 0.5,
                                algorithm = 1,
                                random.seed = 2025,
                                graph.name = "position_snn",
                                cluster.name = "Louvain_cluster_posi")

  seu_obj_GMM_2@meta.data$barcode <- rownames(seu_obj_GMM_2@meta.data)
  seu_metadata <- as.data.table(seu_obj_GMM_2@meta.data)

  # add cluster symbol
  metadata_with_GMM_1[,Louvain_cluster_posi := 1]

  diff_value <- seu_metadata[,Louvain_cluster_posi] %>%
    as.numeric() %>%
    unique() %>%
    unlist() %>%
    min()

  metadata_with_GMM_2[,Louvain_cluster_posi := (as.numeric(seu_metadata[,Louvain_cluster_posi]) + 2 - diff_value)]

  metadata <- rbindlist(list(metadata_with_GMM_1,metadata_with_GMM_2))

  metadata_order <- match(metadata_with_GMM[,barcode],metadata[,barcode])

  return(metadata[metadata_order,])

}

#' conduct Louvain cluster on custom features
#'
#' @param seu_meta the Seurat Object metadata
#' @param raw_count_matrix the row count matrix
#' @param custom_features the custom features

.conduct_Louvain_cluster_on_custom_features <- function(seu_meta,raw_count_matrix,custom_features) {

  on.exit(gc())

  barcodes <- seu_meta[,cell_ID]

  filtered_count_matrix <- raw_count_matrix[custom_features,barcodes]

  dataset_for_KNN <- filtered_count_matrix %>%
    as.matrix() %>%
    t()

  rownames(dataset_for_KNN) <- barcodes

  KNN_graph <- dataset_for_KNN %>%
    FindNeighbors(k.param = 100)

  seu_obj <- CreateSeuratObject(filtered_count_matrix)

  seu_obj@graphs$custom_nn <- KNN_graph[["nn"]]
  seu_obj@graphs$custom_snn <- KNN_graph[["snn"]]

  seu_obj <- FindClusters(object = seu_obj,
                          resolution = 0.5,
                          algorithm = 1,
                          random.seed = 2025,
                          graph.name = "custom_snn",
                          cluster.name = "Louvain_cluster_custom_feature")


  return(seu_obj@meta.data)

}

#' conduct HMRF cluster
#'
#' @param giotto_instruction the instruction of Giotto Object
#' @param seu_meta the Seurat Object metadata
#' @param raw_count_matrix the raw count matrix
#' @param symbol_genes the symbol genes

.conduct_HMRF_cluster <- function(giotto_instruction,seu_meta,raw_count_matrix,symbol_genes) {

  on.exit(gc())

  filtered_matrix <- raw_count_matrix[symbol_genes,seu_meta[,cell_ID]]

  gio_obj <- createGiottoObject(expression = filtered_matrix,
                                spatial_locs = seu_meta[,c("imagerow","neg_imagecol")],
                                instructions = giotto_instruction) %>%
    normalizeGiotto(scalefactor = 1e6) %>%
    addStatistics() %>%
    createSpatialDelaunayNetwork() %>%
    doHMRF(k = 2)

  return(gio_obj)

}




