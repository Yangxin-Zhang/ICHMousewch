
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

#' find filtered barcodes
#'
#' @param raw_count_matrix the matrix of raw count dataset
#' @param original_seu_metadata the original Seurat Object metadata

.find_filtered_barcodes <- function(raw_count_matrix,original_seu_metadata) {

  on.exit(gc())

  up_nCount <- quantile(original_seu_metadata[tissue == 1,nCount_RNA],probs = 0.99)
  down_nCount <- 0
  up_nFeature <- quantile(original_seu_metadata[tissue == 1,nFeature_RNA],probs = 0.99)
  down_nFeature <- 0
  up_per.mt <- quantile(original_seu_metadata[tissue == 1,percent.mt],probs = 0.99,na.rm = TRUE)
  down_per.mt <- 0

  filtered_barcodes <- original_seu_metadata[tissue == 1 &
                                               nCount_RNA > down_nCount &
                                               nCount_RNA <= up_nCount &
                                               nFeature_RNA > down_nFeature &
                                               nFeature_RNA <= up_nFeature &
                                               percent.mt >= down_per.mt &
                                               percent.mt <= up_per.mt,barcode]

  return(filtered_barcodes)

}

#' find filtered genes
#'
#' @param raw_count_matrix the matrix of raw count dataset
#' @param original_seu_metadata the original Seurat Object metadata

.find_filtered_genes <- function(raw_count_matrix,original_seu_metadata) {

  on.exit(gc())

  filtered_barcodes <- ICHMousewch:::.find_filtered_barcodes(raw_count_matrix = raw_count_matrix,
                                                             original_seu_metadata = original_seu_metadata)

  filtered_count_matrix <- raw_count_matrix[,filtered_barcodes]

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
                                                                         giotto_instruction = ich_mouse@giotto_instruction[[1]])

  Louvain_cluster_posi <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = ich_mouse@seu_metadata_with_cluster_symbol,
                                                                                  cluster_symbol = "Louvain_cluster_posi",
                                                                                  raw_count_matrix = ich_mouse@raw_count_matrix,
                                                                                  background_image_address = ich_mouse@file_address["background_image_address"],
                                                                                  self_definition_color = c("1"="#F5D2A8"),
                                                                                  giotto_instruction = ich_mouse@giotto_instruction[[1]])

  original_hematoma <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = ich_mouse@seu_metadata_with_cluster_symbol,
                                                                      cluster_symbol = "hematoma_symbol",
                                                                      raw_count_matrix = ich_mouse@raw_count_matrix,
                                                                      background_image_address = ich_mouse@file_address["background_image_address"],
                                                                      self_definition_color = c("1"="#F5D2A8","2"="#D1352B"),
                                                                      giotto_instruction = ich_mouse@giotto_instruction[[1]],
                                                                      plot_title = "Hematoma",
                                                                      show_plot_title = TRUE)

  Louvain_cluster_filt_gene <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = ich_mouse@seu_metadata_with_cluster_symbol,
                                                                                       cluster_symbol = "Louvain_cluster_filt_gene",
                                                                                       raw_count_matrix = ich_mouse@raw_count_matrix,
                                                                                       background_image_address = ich_mouse@file_address["background_image_address"],
                                                                                       self_definition_color = c("1"="#F5D2A8","2"="#D1352B"),
                                                                                       giotto_instruction = ich_mouse@giotto_instruction[[1]])

  hematoma_center_edge <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(in_tissue_metadata = ich_mouse@seu_metadata_with_cluster_symbol,
                                                                                  cluster_symbol = "center_edge_symbol",
                                                                                  raw_count_matrix = ich_mouse@raw_count_matrix,
                                                                                  background_image_address = ich_mouse@file_address["background_image_address"],
                                                                                  self_definition_color = c("1"="#F5D2A8","2"="#D1352B","3"="#3C77AF"),
                                                                                  giotto_instruction = ich_mouse@giotto_instruction[[1]],
                                                                                  plot_title = "Center__Edge",
                                                                                  show_plot_title = TRUE)

  log2Count <- ICHMousewch:::.create_count_distribution_map(seu_meta = ich_mouse@seu_metadata_with_cluster_symbol)

  plotting_list <- list("log2Count" = ggplotGrob(log2Count),
                        "GMM_cluster" = ggplotGrob(GMM_cluster),
                        "Louvain_cluster_posi" = ggplotGrob(Louvain_cluster_posi),
                        "original_hematoma" = ggplotGrob(original_hematoma),
                        "Louvain_cluster_filt_gene" = ggplotGrob(Louvain_cluster_filt_gene),
                        "hematoma_center_edge" = ggplotGrob(hematoma_center_edge))

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
    ICHMousewch:::.plotting_theme() +
    guides(color = guide_legend(keywidth = unit(0.8, "cm"),keyheight = unit(0.8, "cm"),
                                override.aes = list(size = 3,alpha = 1)))

  volcano_plot_ls <- list("diff_expr_gene" = ggplotGrob(volcano_plot))

  return(volcano_plot_ls)

}

#' plotting bubble chart
#'
#' @param ich_mouse the ICH_Mouse class

.plotting_bubble_chart <- function(ich_mouse,gene_num = 10) {

  on.exit(gc())

  GO_term_set_ls <- ich_mouse@GO_ID_group

  enri_set <- ICHMousewch::Create_Enrichment_Set(ich_mouse = ich_mouse)
  enri_set <- ICHMousewch::add_GO_term_set(enrichment_set = enri_set,
                                           GO_term_set_ls = GO_term_set_ls,
                                           ich_mouse = ich_mouse)

  gene_info <- enri_set@gene_information
  GO_info <- enri_set@GO_set

  plotting_dataset <- ICHMousewch:::.generate_bubble_dataset(gene_info = gene_info,
                                                             GO_info = GO_info,
                                                             gene_num = gene_num)

  posi_y_da <- plotting_dataset[y_symbol == "posi_y"] %>%
    as.data.frame() %>%
    mutate( GO_ID = factor(GO_ID,levels = unique(plotting_dataset[,GO_ID])),
            GO_group = factor(GO_group,levels = unique(plotting_dataset[,GO_group])),
            text_y = rep(0,times = length(plotting_dataset[y_symbol == "posi_y",gene_name])))

  negt_y_da <- plotting_dataset[y_symbol == "negt_y"] %>%
    as.data.frame() %>%
    mutate(GO_ID = factor(GO_ID,levels = unique(plotting_dataset[,GO_ID])),
           GO_group = factor(GO_group,levels = unique(plotting_dataset[,GO_group])))

  y_tick <- ICHMousewch:::.generate_bubble_chart_y_tick(plotting_dataset = plotting_dataset)

  bubble_color <- plotting_dataset[y_symbol == "negt_y",bubble_color] %>%
    unique()
  names(bubble_color) <- bubble_color

  bubble_chart <- ggplot(data = plotting_dataset) +
    geom_point(data = negt_y_da,
               mapping = aes(x = GO_ID, y = bubble_y,color = bubble_color),
               size = negt_y_da$avg_log2FC,
               alpha = 0.5) +
    scale_color_manual(values = bubble_color) +
    geom_bar(data = posi_y_da,
             mapping = aes(x = GO_ID, y = bubble_y,fill = GO_group),
             stat = "identity",
             width = 0.5) +
    geom_hline(yintercept = 0,
               color = "#000000FF",
               linewidth = 0.5) +
    geom_text(data = posi_y_da,
              mapping = aes(x = GO_ID,y = text_y,label = GO_ID),
              vjust = 1.5,
              size = 0,
              size.unit = "pt",
              family = "Arial") +
    scale_y_continuous(breaks = y_tick,
                       labels = names(y_tick)) +
    ICHMousewch:::.plotting_theme() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    guides(fill = guide_legend(position = "top",
                               nrow = 2,
                               theme = theme(legend.text = element_text(size = 12,
                                                                        family = "Arial",
                                                                        vjust = 0.5,
                                                                        hjust = 0),
                                             legend.title = element_blank(),
                                             legend.key.size = unit(12,"pt"))),
           color = guide_none())

  bubble_chart_ls <- list("bubble_chart" = ggplotGrob(bubble_chart))

  return(bubble_chart_ls)

}

#' plotting GMM barchart
#'
#' @slot ich_mouse the class of ICH_Mouse

.plotting_GMM_barchart <- function(ich_mouse) {

  on.exit(gc())

  seu_metadata <- ich_mouse@seu_metadata_with_cluster_symbol
  seu_metadata_normal <- ich_mouse@seu_metadata_with_cluster_symbol[hematoma_symbol == 1]
  seu_metadata_hematoma <- ich_mouse@seu_metadata_with_cluster_symbol[hematoma_symbol == 2]

  barchart_ori <- ICHMousewch:::.create_ncount_nfeature_histogram(seu_metadata = seu_metadata,
                                                                  plotting_title = "Origin")

  barchart_normal <- ICHMousewch:::.create_ncount_nfeature_histogram(seu_metadata = seu_metadata_normal,
                                                                     plotting_title = "Normal")
  barchart_hematoma <- ICHMousewch:::.create_ncount_nfeature_histogram(seu_metadata = seu_metadata_hematoma,
                                                                       plotting_title = "Hematoma",
                                                                       reso = 0.4)

  barchart <- c(barchart_ori,barchart_normal,barchart_hematoma)

  return(barchart)

}

#' choose bins for histogram
#'
#' @param his_dataset histogram dataset

.choose_bins_for_histogram <- function(his_dataset,reso = 0.05,bins = 0) {

  on.exit(gc())

  his_dataset <- data.table(his_data = his_dataset)

  if (bins == 0) {

    bins <- (length(unique(his_dataset[,his_data]))*reso) %>%
      round()

  }

  count_da <- ggplot_build(ggplot(his_dataset,mapping = aes(x=his_data)) +
                             geom_histogram(bins = bins))
  count_da <- count_da$data[[1]] %>%
    as.data.table()

  count_max <- max(count_da[,y])

  high <- count_max %>%
    as.character() %>%
    strsplit(split = "")
  high <- high[[1]]

  high_two <- as.integer(high[2])

  if (length(high) < 3) {
    if (length(high) == 1) {

      high = 10

    } else {

      high_one <- (as.integer(high[1]) + 1) %>%
        as.character()
      high[1] <- high_one
      high[2:length(high)] <- rep("0",times = (length(high)-1))

    }

  } else {

    if (high_two == 9) {

      high_one <- (as.integer(high[1]) + 1) %>%
        as.character()
      high[1] <- high_one
      high[2:length(high)] <- rep("0",times = (length(high)-1))

      high <- paste(high,collapse = "") %>%
        as.integer()

    } else {

      if (high_two == 5) {

        high_two <- high_two + 2
        high[2] <- high_two
        high[3:length(high)] <- rep("0",times = (length(high)-2))

        high <- paste(high,collapse = "") %>%
          as.integer()

      } else {

        high_two <- high_two + 1
        high[2] <- high_two
        high[3:length(high)] <- rep("0",times = (length(high)-2))

        high <- paste(high,collapse = "") %>%
          as.integer()

      }



    }

  }

  y_tick <- high %>%
    as.character() %>%
    strsplit(split = "")

  y_tick <- y_tick[[1]]

  inter <- character(length(y_tick)-1)
  inter[1] <- "5"
  inter[2:length(inter)] <- "0"
  inter <- paste(inter,collapse = "") %>%
    as.integer()

  y_tick <- seq(from = 0,to = high,by = inter)

  if (y_tick[length(y_tick)] != high) {

    y_tick <- c(y_tick,high)

  }

  setorder(count_da,x)
  mo_avg <- stats::filter(count_da[,y],rep(1/3,3),sides = 2) %>%
    as.numeric()
  mo_avg[is.na(mo_avg)] <- 0
  sel_y <- mo_avg  == 0
  count_zero <- count_da[sel_y]
  x_count_max <- count_da[y == count_max,x]

  x_min <- max(count_zero[xmax < x_count_max,xmax]) %>%
    floor()

  x_max <- min(count_zero[xmin > x_count_max,xmin]) %>%
    ceiling()

  x_tick <- seq(from = x_min,to = x_max, length.out = 5)

  x_width <- c(x_min,x_max)

  results <- list("bins" = bins,
                  "y_high" = c(0,high),
                  "y_tick" = y_tick,
                  "x_width" = c(x_min,x_max),
                  "x_tick" = x_tick)

  return(results)

}

#' plotting violin plot
#'
#' @param ich_mouse the class of ICH_Mouse

.plotting_violin_plot <- function(ich_mouse) {

  on.exit(gc())

  ncount_violin <- ICHMousewch:::.create_violin_plot_of_quality_control(seu_metadata = ich_mouse@original_seu_metadata,
                                                                        filtered_barcodes = ich_mouse@filtered_barcodes,
                                                                        plotting_var = "nCount_RNA")
  nfeature_violin <- ICHMousewch:::.create_violin_plot_of_quality_control(seu_metadata = ich_mouse@original_seu_metadata,
                                                                          filtered_barcodes = ich_mouse@filtered_barcodes,
                                                                          plotting_var = "nFeature_RNA")
  percent.mt_violin <- ICHMousewch:::.create_violin_plot_of_quality_control(seu_metadata = ich_mouse@original_seu_metadata,
                                                                        filtered_barcodes = ich_mouse@filtered_barcodes,
                                                                        plotting_var = "percent.mt")

  violin_ls <- list("ncount_violin" = ncount_violin,
                    "nfeature_violin" = nfeature_violin,
                    "percent.mt_violin" = percent.mt_violin)

  return(violin_ls)

}

#' plotting umap plots
#'
#' @param ich_mouse the class of ICH_Mouse

.plotting_umap_plot <- function(ich_mouse,diff_expr_genes_symbol = "edge-normal") {

  on.exit(gc())

  raw_count_matrix <- ich_mouse@raw_count_matrix
  seu_meta <- ich_mouse@seu_metadata_with_cluster_symbol
  filter_genes <- ich_mouse@filtered_genes
  diff_expr_genes <- ich_mouse@diff_expr_genes[[diff_expr_genes_symbol]]

  umap_center_edge <- ICHMousewch:::.create_umap_plot(raw_count_matrix = raw_count_matrix,
                                                      seu_meta = seu_meta,
                                                      filter_genes = filter_genes,
                                                      diff_expr_genes = diff_expr_genes,
                                                      group_symbol = "center_edge_symbol")

  umap_Louvain_posi <- ICHMousewch:::.create_umap_plot(raw_count_matrix = raw_count_matrix,
                                                       seu_meta = seu_meta,
                                                       filter_genes = filter_genes,
                                                       diff_expr_genes = diff_expr_genes,
                                                       group_symbol = "Louvain_cluster_posi")

  umap_plot_ls <- list("umap_center_edge" = umap_center_edge,
                       "umap_Louvain_posi" = umap_Louvain_posi)

  return(umap_plot_ls)

}

#' generate GO term gene contained matrix
#'
#' @param GO_ID_group the GO terms
#' @param GO_enrichment the result of GO enrichment

.generate_GO_term_gene_contained_matrix <- function(GO_ID_group,GO_enrichment) {

  on.exit(gc())

  sub_GO_enri <- GO_enrichment[ID %in% GO_ID_group]

  gene_ls <- vector("list",length = length(GO_ID_group))
  names(gene_ls) <- GO_ID_group
  for (i in 1:length(GO_ID_group)) {

    genes <- sub_GO_enri[ID == GO_ID_group[i],geneID] %>%
      strsplit(split = "/")
    gene_ls[GO_ID_group[i]] <- genes

  }

  contained_matrix <- vector("list",length = length(GO_ID_group))
  names(contained_matrix) <- GO_ID_group
  for (i in 1:length(GO_ID_group)) {

    con_vec <- vector("integer",length = length(GO_ID_group))
    names(con_vec) <- GO_ID_group
    for (j in 1:length(GO_ID_group)) {

      num <- sum(!gene_ls[[GO_ID_group[i]]] %in% gene_ls[[GO_ID_group[j]]])

      con_vec[GO_ID_group[j]] <- num

    }

    contained_matrix[GO_ID_group[i]] <- list(con_vec)

  }

  contained_matrix <- as.data.frame(contained_matrix)

  return(contained_matrix)

}

#' filter GO term based on contained matrix
#'
#' @param GO_ID_group the GO terms
#' @param GO_enrichment the result of GO enrichment

.filter_GO_term_based_on_contained_matrix <- function(GO_ID_group,GO_enrichment) {

  on.exit(gc())

  contained_matrix <- ICHMousewch:::.generate_GO_term_gene_contained_matrix(GO_ID_group = GO_ID_group,
                                                                            GO_enrichment = GO_enrichment)

  zero_symbol <- which(contained_matrix == 0,arr.ind = TRUE) %>%
    as.data.table()

  zero_symbol <- zero_symbol[!row == col]

  de <- GO_ID_group[unique(zero_symbol[,col])]

  sa <- GO_ID_group[!GO_ID_group %in% de]

  return(GO_enrichment[ID %in% sa])

}
