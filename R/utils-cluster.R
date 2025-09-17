
# R/utils-cluster.R

#' conduct GMM classification
#'
#' @param original_seu_metadata original Seurat Object metadata
#' @param in_tissue_barcode the barcodes in tissue

.conduct_GMM_classification <- function(original_seu_metadata,in_tissue_barcode) {

  on.exit(gc())

  set.seed(2025)

  in_tissue_metadata <- original_seu_metadata[in_tissue_barcode,]

  in_tissue_metadata$GMM_cluster <- Mclust(data = in_tissue_metadata[,c("nCount_log2","nFeature_log2")],G = 2)[["classification"]]

  num_1 <- sum(in_tissue_metadata$GMM_cluster == 1)
  num_2 <- sum(in_tissue_metadata$GMM_cluster == 2)

  if(num_1 < num_2) {

    in_tissue_metadata$GMM_cluster[metadata$GMM_cluster == 1] <- 0
    in_tissue_metadata$GMM_cluster[metadata$GMM_cluster == 2] <- 1
    in_tissue_metadata$GMM_cluster[metadata$GMM_cluster == 0] <- 2

  }

  rownames(in_tissue_metadata) <- in_tissue_metadata$barcode

  return(in_tissue_metadata)

}

#' conduct Louvain cluster based on tissue position matrix
#'
#' @param original_seu_metadata original Seurat Object metadata
#' @param in_tissue_barcode the barcodes in tissue
#' @param raw_count_matrix the matrix of raw count dataset

.conduct_Louvain_cluster_based_on_tissue_position <- function(original_seu_metadata,in_tissue_barcode,raw_count_matrix) {

  on.exit(gc())

  metadata_with_GMM <- ICHMousewch:::.conduct_GMM_classification(original_seu_metadata = original_seu_metadata,
                                                                 in_tissue_barcode = in_tissue_barcode)

  metadata_with_GMM_1 <- metadata_with_GMM[metadata_with_GMM$GMM_cluster == 1,]

  metadata_with_GMM_2 <- metadata_with_GMM[metadata_with_GMM$GMM_cluster == 2,]

  # create Seurat Object according to GMM cluster
  raw_count_matrix_GMM_2 <- raw_count_matrix[,metadata_with_GMM_2$barcode]
  seu_obj_GMM_2 <- CreateSeuratObject(raw_count_matrix_GMM_2)

  # conduct Louvain algorithm
  KNN_graph <- metadata_with_GMM_2[,c("imagerow","imagecol")] %>%
    as.matrix() %>%
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

  # add cluster symbol
  metadata_with_GMM_1$Louvain_cluster_posi <- 1

  diff_value <- seu_obj_GMM_2@meta.data$Louvain_cluster_posi %>%
    as.numeric() %>%
    unique() %>%
    unlist() %>%
    min()

  metadata_with_GMM_2$Louvain_cluster_posi <- as.numeric(seu_obj_GMM_2@meta.data$Louvain_cluster_posi) + 2 - diff_value
  metadata_with_GMM_2 <- metadata_with_GMM_2[,colnames(metadata_with_GMM_1)]

  metadata <- bind_rows(metadata_with_GMM_1,metadata_with_GMM_2)

  metadata_order <- match(metadata_with_GMM$barcode,metadata$barcode)

  return(metadata[metadata_order,])

}







