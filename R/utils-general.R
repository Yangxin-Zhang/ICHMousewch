
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

  diff_expr_genes <- diff_expr_genes[abs(avg_log2FC) > 1]

  setorder(diff_expr_genes,-log2pct)

  return(diff_expr_genes)

}

#' annotate the cell type based on single gene
#'
#' @param ich_mouse the ICH_Mouse class
#' @param gene_ls the gene ls used to annotate
#' @param reference_dataset the reference_dataset for annotation

.annotate_the_cell_type_based_on_single_gene <- function(ich_mouse,gene_ls,reference_dataset) {

  on.exit(gc())

  ref_ds <- reference_dataset

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

#' export data.table as excel
#'
#' @param data.table_obj a data.table object or list
#' @param saving_path the saving path
#' @param file_name the workbook name
#' @export

export_data.table_as_excel <- function(data.table_obj,saving_path,file_name) {

  on.exit(gc())

  wb <- createWorkbook()

  if (is.list(data.table_obj)) {

    sheet_na <- names(data.table_obj)

    for (i in 1:length(sheet_na)) {

      addWorksheet(wb = wb,
                   sheetName = sheet_na[i])

      writeData(wb = wb,
                sheet = sheet_na[i],
                x = data.table_obj[[sheet_na[i]]])

    }

  }

  if (is.data.frame(data.table_obj)) {

    sheet_na <- file_name

    addWorksheet(wb = wb,
                 sheetName = sheet_na)

    writeData(wb = wb,
              sheet = sheet_na,
              x = data.table_obj)

  }

  saving_name <- paste(file_name,"xlsx",sep = ".")
  file_path <- paste(saving_path,saving_name,sep = "/")
  saveWorkbook(wb = wb,
               file = file_path,
               overwrite = TRUE)

}

#' plotting spatial image
#'
#' @param ich_mouse the class of ICH_Mouse

.plotting_spatial_image <- function(ich_mouse) {

  on.exit(gc())

  GMM_cluster <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = ich_mouse@seu_metadata_with_cluster_symbol,
                                                                         cluster_symbol = "GMM_cluster",
                                                                         raw_count_matrix = ich_mouse@raw_count_matrix,
                                                                         background_image_address = ich_mouse@file_address["background_image_address"],
                                                                         self_definition_color = c("1"="#F5D2A8","2"="#D1352B"),
                                                                         giotto_instruction = ich_mouse@giotto_instruction[[1]],
                                                                         plot_title = "GMM_Cluster") %>%
    ggplotGrob()

  Louvain_cluster_posi <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = ich_mouse@seu_metadata_with_cluster_symbol,
                                                                                  cluster_symbol = "Louvain_cluster_posi",
                                                                                  raw_count_matrix = ich_mouse@raw_count_matrix,
                                                                                  background_image_address = ich_mouse@file_address["background_image_address"],
                                                                                  self_definition_color = c("1"="#F5D2A8"),
                                                                                  giotto_instruction = ich_mouse@giotto_instruction[[1]],
                                                                                  plot_title = "Louvain_Cluster_Posi") %>%
    ggplotGrob()

  hematoma <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = ich_mouse@seu_metadata_with_cluster_symbol,
                                                                      cluster_symbol = "hematoma_symbol",
                                                                      raw_count_matrix = ich_mouse@raw_count_matrix,
                                                                      background_image_address = ich_mouse@file_address["background_image_address"],
                                                                      self_definition_color = c("1"="#F5D2A8","2"="#D1352B"),
                                                                      giotto_instruction = ich_mouse@giotto_instruction[[1]],
                                                                      plot_title = "Hematoma")

  Louvain_cluster_filt_gene <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = ich_mouse@seu_metadata_with_cluster_symbol,
                                                                                       cluster_symbol = "Louvain_cluster_filt_gene",
                                                                                       raw_count_matrix = ich_mouse@raw_count_matrix,
                                                                                       background_image_address = ich_mouse@file_address["background_image_address"],
                                                                                       self_definition_color = c("1"="#F5D2A8","2"="#D1352B"),
                                                                                       giotto_instruction = ich_mouse@giotto_instruction[[1]],
                                                                                       plot_title = "Louvain_Cluster_Filt_Gene")

  hematoma_center_edge <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = ich_mouse@seu_metadata_with_cluster_symbol,
                                                                                  cluster_symbol = "center_edge_symbol",
                                                                                  raw_count_matrix = ich_mouse@raw_count_matrix,
                                                                                  background_image_address = ich_mouse@file_address["background_image_address"],
                                                                                  self_definition_color = c("1"="#F5D2A8","2"="#D1352B","3"="#3C77AF"),
                                                                                  giotto_instruction = ich_mouse@giotto_instruction[[1]],
                                                                                  plot_title = "Hematoma_Center_Edge")

  plotting_list <- list("GMM_cluster" = GMM_cluster,
                        "Louvain_cluster_posi" = Louvain_cluster_posi,
                        "hematoma" = hematoma,
                        "Louvain_cluster_filt_gene" = Louvain_cluster_filt_gene,
                        "hematoma_center_edge" = hematoma_center_edge)

  return(plotting_list)

}

#' get another dataset for plotting
#'
#' @param ich_mouse the class of ich_mouse

.get_another_dataset_for_plotting <- function(ich_mouse) {

  on.exit(gc())

  another_dataset_na <- c("raw_count_matrix","background_image_address","color_set","giotto_instruction")

  another_dataset <- vector("list",length = length(another_dataset_na))

  another_dataset["raw_count_matrix"] <- list(ich_mouse@raw_count_matrix)
  another_dataset["background_image_address"] <- list(ich_mouse@file_address$background_image_address)
  another_dataset["color_set"] <- list(ich_mouse@color_set)
  another_dataset["giotto_instruction"] <- list(ich_mouse@giotto_instruction)

  return(another_dataset)

}

#' plotting volcano plot
#'
#' @param ich_mouse the ICH_Mouse class

.plotting_volcano_plot <- function(ich_mouse) {

  on.exit(gc())

  diff_expr_mat <- ich_mouse@diff_expr_genes$`edge-normal`

  diff_expr_mat[,Threshold := rep("No",nrow(diff_expr_mat))]

  diff_expr_mat[avg_log2FC >= 1 & p_val_adj < 0.01,Threshold := "Up"]

  diff_expr_mat[,log10p_val_adj := -log10(p_val_adj)]

  diff_expr_mat[log10p_val_adj > 300,log10p_val_adj := 300]
  diff_expr_mat[avg_log2FC > 10, avg_log2FC := 10]
  diff_expr_mat[avg_log2FC < -10, avg_log2FC := -10]

  volcano_plot <- ggplot(data = diff_expr_mat,mapping = aes(x = avg_log2FC, y = log10p_val_adj, color = Threshold)) +
    labs(title = "Volcano_Plot") +
    geom_point(size = 1,alpha = 0.7) +
    coord_cartesian(xlim = c(-10, 10), ylim = c(-1, 300)) +
    scale_x_continuous(breaks = c(-10,-5,-1,0,1,5,10)) +
    scale_y_continuous(breaks = c(0,100,200,300)) +
    scale_color_manual(values = c("No" = "#999999", "Up" = "#E41A1C")) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
    theme(plot.title = element_text(family = "Arial",size = 12,color = "black",face = "bold",hjust = 0.5,vjust = 0.5,margin = margin(b = 10, t = 10)),
          axis.text = element_text(family = "Arial",size = 12,color = "black",hjust = 0.5,vjust = 0.5),
          axis.title.x = element_text(family = "Arial",size = 12,color = "black",hjust = 0.5,vjust = 0.5,margin = margin(b = 10, t = 10)),
          axis.title.y = element_text(family = "Arial",size = 12,color = "black",hjust = 0.5,vjust = 0.5,margin = margin(r = 10, l = 10)),
          legend.text = element_text(family = "Arial",size = 12,color = "black",hjust = 0.5,vjust = 0.5),
          legend.title = element_text(family = "Arial",size = 12,color = "black",hjust = 0.5,vjust = 0.5,margin = margin(b = 10)),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(fill = "white", color = NA),
          legend.position = "right",
          panel.background = element_rect(fill = "white", color = "black"),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(keywidth = unit(0.8, "cm"),keyheight = unit(0.8, "cm"),
                                override.aes = list(size = 3,alpha = 1)))

  volcano_plot_ls <- list("diff_expr_gene" = ggplotGrob(volcano_plot))

  return(volcano_plot_ls)

}









