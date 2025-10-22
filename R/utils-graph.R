
# R/utils-graph.R

#' create Giotto Instruction
#'
#' @param python_path the path to a python which the Giotto can use
#' @param results_folder the folder for Giotto save plots

.create_giotto_instruction <- function(python_path, results_folder) {

  on.exit(gc())

  instrs = createGiottoInstructions(save_dir = results_folder,
                                    save_plot = FALSE,
                                    show_plot = FALSE,
                                    python_path = python_path)

  return(instrs)

}

#' create spatial image with cluster symbol
#'
#' @param in_tissue_metadata the Seurat Object metadata of in tissue barcodes
#' @param cluster_symbol the cluster for spatial image
#' @param raw_count_matrix the matrix of raw count dataset
#' @param background_image_address image for background
#' @param self_definition_color a set of color defined by user
#' @param giotto_instruction the instruction of Giotto Object

.create_spatial_image_with_cluster_symbol <- function(in_tissue_metadata,cluster_symbol,raw_count_matrix,background_image_address,giotto_instruction,self_definition_color = character(),theme_param = list(),plot_title = character(),legend_lable = character(),show_image = TRUE,show_legend_label = FALSE,show_plot_title = FALSE) {

  on.exit(gc())

  in_tissue_metadata[,cell_ID := barcode]

  if (show_plot_title == FALSE) {

    theme_param <- c(theme_param,list(plot.title = element_blank()))

  } else {

    if (length(plot_title) != 0) {

      aim_col <- in_tissue_metadata[,..cluster_symbol]

      in_tissue_metadata[,(plot_title) := aim_col]

      cluster_symbol <- plot_title

    }

    theme_param <- c(theme_param,list(plot.title = element_text(family = "Arial",
                                                                size = 16,
                                                                color = "black",
                                                                face = "bold",
                                                                hjust = 0.5,
                                                                vjust = 0.5,
                                                                margin = margin(b = 10, t = 10))))

  }

  if (show_legend_label == FALSE) {

    theme_param <- c(theme_param,list(legend.position = "none"))

  } else {

    if (length(legend_lable) != 0) {

      ori_lable <- names(legend_lable)

      for (i in 1:length(ori_lable)) {

        filter_condition <- as.character(in_tissue_metadata[,..cluster_symbol]) %in% ori_lable[i]

        in_tissue_metadata[filter_condition,..cluster_symbol := legend_lable[ori_lable[i]]]

      }

      if (length(self_definition_color) != 0) {

        for (i in 1:length(legend_lable)) {

          new_self_def_col <- vector("character",length = length(ori_lable))
          names(new_self_def_col) <- legend_lable

          new_self_def_col[legend_lable[i]] <- self_definition_color[names(legend_lable[i])]

        }

        self_definition_color <- new_self_def_col

      }

    }

    theme_param <- c(theme_param,list(legend.position = "top",
                                      legend.text = element_text(family = "Arial",
                                                                 size = 12,
                                                                 color = "black",
                                                                 face = "bold",
                                                                 hjust = 0,
                                                                 vjust = 0.5),
                                      legend.key.size = unit(12,"pt")))

  }

  in_tissue_count_matrix <- raw_count_matrix[,in_tissue_metadata[,barcode]]

  image_obj <- createGiottoImage(mg_object = background_image_address,
                                 name = "background_image",
                                 negative_y = FALSE)

  giotto_object <- createGiottoObject(expression = in_tissue_count_matrix,
                                      spatial_locs = in_tissue_metadata[,c("imagerow","neg_imagecol")],
                                      cell_metadata = in_tissue_metadata,
                                      instructions = giotto_instruction,
                                      images = list(image_obj))

  color_symbols <- in_tissue_metadata %>%
    subset(select = cluster_symbol) %>%
    unique() %>%
    unlist() %>%
    as.character()

  color_set <- viridis::plasma(length(color_symbols))

  random_color_symbols <- color_symbols[!color_symbols %in% names(self_definition_color)]

  sub_color_set <- color_set[!color_set %in% self_definition_color]

  if(length(random_color_symbols) > 0) {

    select_color <- sample(1:length(sub_color_set),length(random_color_symbols))
    random_colors <- sub_color_set[select_color]
    names(random_colors) <- random_color_symbols

  } else {

    random_colors <- character()

  }

  spatial_image <- spatPlot2D(gobject = giotto_object,
                              cell_color = cluster_symbol,
                              point_size = 0.5,
                              point_shape = "border",
                              point_alpha = 0.5,
                              point_border_stroke = 0,
                              cell_color_code = c(self_definition_color,random_colors),
                              background_color = "#00000000",
                              show_image = TRUE,
                              theme_param = c(theme_param,list(plot.background = element_rect(fill = "white",color = NA),
                                                               axis.ticks.x = element_blank(),
                                                               axis.ticks.y = element_blank(),
                                                               axis.text.x = element_blank(),
                                                               axis.text.y = element_blank(),
                                                               axis.title.x = element_blank(),
                                                               axis.title.y = element_blank())))

  return(spatial_image)

}

#' create spatial image with single gene
#'
#' @param seu_metadata_with_cluster_symbol the Seurat Object metadata of in tissue barcodes
#' @param gene_ls the gene list for creating spatial image
#' @param raw_count_matrix the matrix of raw count dataset
#' @param background_image_address image for background
#' @param giotto_instruction the instruction of Giotto Object
#' @param show_background_image whether to show background image

.create_spatial_image_with_single_gene <- function(seu_metadata_with_cluster_symbol,gene_ls,raw_count_matrix,background_image_address,giotto_instruction,show_background_image = TRUE,gradient = FALSE) {

  on.exit(gc())

  in_tissue_metadata <- seu_metadata_with_cluster_symbol[,cell_ID := barcode]

  in_tissue_count_matrix <- raw_count_matrix[gene_ls,seu_metadata_with_cluster_symbol[,barcode]]

  if (gradient) {

    edge_barcodes <- in_tissue_metadata[center_edge_symbol == "3",barcode]
    normal_barcodes <- in_tissue_metadata[center_edge_symbol == "1",barcode]

    in_tissue_count_matrix[!in_tissue_count_matrix == 0 & colnames(in_tissue_count_matrix) %in% edge_barcodes] <- 1
    in_tissue_count_matrix[!in_tissue_count_matrix == 0 & colnames(in_tissue_count_matrix) %in% normal_barcodes] <- 2

  } else {

    in_tissue_count_matrix[!in_tissue_count_matrix == 0] <- 1

  }

  image_obj <- createGiottoImage(mg_object = background_image_address,
                                 name = "background_image",
                                 negative_y = FALSE)

  giotto_object <- createGiottoObject(expression = in_tissue_count_matrix,
                                      spatial_locs = in_tissue_metadata[,c("imagerow","neg_imagecol")],
                                      cell_metadata = in_tissue_metadata,
                                      instructions = giotto_instruction,
                                      images = list(image_obj))

  spatial_image_ls <- list()

  if(gradient) {

    if(length(gene_ls) > 0) {

      for (i in 1:length(gene_ls)) {

        spatial_image <- spatFeatPlot2D(gobject = giotto_object,
                                        expression_values = "raw",
                                        feats = gene_ls[i],
                                        background_color = "white",
                                        point_size = 0.5,
                                        point_alpha = 0.5,
                                        cell_color_gradient = c("#F5D2A8","#D1352B","#3C77AF"),
                                        gradient_midpoint = 1,
                                        gradient_style = "divergent",
                                        show_image = show_background_image,
                                        theme_param = list(plot.background = element_rect(fill = "white",color = NA),
                                                           axis.ticks.x = element_blank(),
                                                           axis.ticks.y = element_blank(),
                                                           axis.text.x = element_blank(),
                                                           axis.text.y = element_blank(),
                                                           axis.title.x = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           legend.position = "none")) %>%
          ggplotGrob()

        spatial_image_ls <- append(spatial_image_ls,list(spatial_image))
        names(spatial_image_ls)[i] <- gene_ls[i]

      }

    }

  } else {

    if(length(gene_ls) > 0) {

      for (i in 1:length(gene_ls)) {

        spatial_image <- spatFeatPlot2D(gobject = giotto_object,
                                        expression_values = "raw",
                                        feats = gene_ls[i],
                                        background_color = "white",
                                        point_size = 0.5,
                                        point_alpha = 0.5,
                                        cell_color_gradient = c("#F5D2A8","#D1352B"),
                                        show_image = show_background_image,
                                        theme_param = list(plot.background = element_rect(fill = "white",color = NA),
                                                           axis.ticks.x = element_blank(),
                                                           axis.ticks.y = element_blank(),
                                                           axis.text.x = element_blank(),
                                                           axis.text.y = element_blank(),
                                                           axis.title.x = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           legend.position = "none")) %>%
          ggplotGrob()

        spatial_image_ls <- append(spatial_image_ls,list(spatial_image))
        names(spatial_image_ls)[i] <- gene_ls[i]

      }

    }

  }


  return(spatial_image_ls)

}

#' create Categorical Network Plot of GO term
#'
#' @param GO_results the results of GO enrichment

.create_categorical_network_plot <- function(GO_results) {

  on.exit(gc())

  GO_terms <- GO_results[,Description]

  cnet_ls <- vector("list",length = length(GO_terms))
  names(cnet_ls) <- GO_terms
  for (i in 1:length(GO_terms)) {

    gene_na <- GO_results[Description == GO_terms[i],geneID] %>%
      strsplit(split = "/") %>%
      unlist()

    cnet_ls[GO_terms[i]] <- list(gene_na)

  }

  cnet_plot <- cnetplot(cnet_ls,
                        showCategory = 10) %>%
    ggplotGrob()

  return(cnet_plot)

}

#' plotting cnetplot
#'
#' @param GO_results the results of GO enrichment

.plotting_cnetplot <- function(GO_results) {

  on.exit(gc())

  cluster_na <- names(GO_results)

  cnet_plot_ls <- vector("list",length = length(GO_results))
  names(cnet_plot_ls) <- cluster_na
  for (i in 1:length(cluster_na)) {

    cnet_plot <- ICHMousewch:::.create_categorical_network_plot(GO_results = GO_results[[i]])
    cnet_plot_ls[cluster_na[i]] <- list(cnet_plot)

  }

  return(cnet_plot_ls)

}

#' save gtable plot
#'
#' @param gtable_plot a plot in the form of gtable
#' @param saving_path the path for saving plot
#' @param file_name the file name
#' @export

save_gtable_plot <- function(gtable_plot,saving_path,file_name) {

  on.exit(gc())

  file_name <- paste(file_name,"png",sep = ".")

  ggsave(filename = paste(saving_path,file_name,sep = "/"),
         plot = as.ggplot(gtable_plot),
         dpi = 600,
         width = 297,height = 210,
         units = "mm",
         device = "png")

}

#' plotting GO term similarity heatmap
#'
#' @param ich_mouse the ICH_Mouse class

.plotting_GO_term_similarity_heatmap <- function(ich_mouse) {

  on.exit(gc())

  GO_results <- rbindlist(ich_mouse@GO_cluster$`edge-normal`)

  sim_mat <- GO_similarity(go_id = GO_results[,ID],
                           ont = "BP",
                           db = "org.Mm.eg.db",
                           measure = "Sim_Resnik_1999")

  GO_id <- colnames(sim_mat)
  GO_id_x <- rep(colnames(sim_mat),each = ncol(sim_mat))
  GO_id_y <- rep(rownames(sim_mat),times = nrow(sim_mat))

  sim_dt <- data.table(similarity = numeric(ncol(sim_mat)*nrow(sim_mat)),
                       GO_id_x = GO_id_x,
                       GO_id_y = GO_id_y)

  sim_data <- c()
  for(i in 1:length(GO_id)) {

    sim_data <- c(sim_data,sim_mat[,GO_id[i]])

  }

  sim_dt[,similarity := sim_data]
  sim_dt[similarity < 0.2,similarity := 0]
  sim_dt[similarity > 0.8,similarity := 1]

  sim_dt <- as.data.frame(sim_dt) %>%
    mutate(GO_id_x = factor(GO_id_x,levels = colnames(sim_mat)),
           GO_id_y = factor(GO_id_y,levels = rownames(sim_mat)))

  sim_heatmap <- ggplot(data = sim_dt,mapping = aes(x = GO_id_x,y = GO_id_y,fill = similarity)) +
    geom_tile() +
    scale_fill_gradientn(colours = c("#FFFFFF", "#FFA500", "#E6550D"),
                         values = c(0,0.75,1),
                         limits = c(0,1)) +
    ICHMousewch:::.plotting_theme() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")

  heatmap_ls <- list("GO_heatmap" = ggplotGrob(sim_heatmap))

  return(heatmap_ls)

}

#' create ncount and nfeature histogram
#'
#' @param seu_metadata the Seurat Object metadata
#' @param plotting_title the title of plotting

.create_ncount_nfeature_histogram <- function(seu_metadata,plotting_title,reso = 0.05) {

  on.exit(gc())

  para_ncount <- ICHMousewch:::.choose_bins_for_histogram(his_dataset = seu_metadata[,nCount_log2],reso = reso)
  barchart_ncount <- ggplot(data = seu_metadata,mapping = aes(x = nCount_log2)) +
    geom_histogram(bins = para_ncount[["bins"]],
                   linewidth = 0.9,
                   fill = viridis::plasma(1),
                   color = "white") +
    scale_y_continuous(breaks = para_ncount[["y_tick"]]) +
    scale_x_continuous(breaks = para_ncount[["x_tick"]]) +
    coord_cartesian(xlim = para_ncount[["x_width"]],
                    ylim = para_ncount[["y_high"]]) +
    ICHMousewch:::.plotting_theme()

  para_nfeature <- ICHMousewch:::.choose_bins_for_histogram(his_dataset = seu_metadata[,nFeature_log2],bins = para_ncount[["bins"]],reso = reso)
  barchart_nfeature <- ggplot(data = seu_metadata,mapping = aes(x = nFeature_log2))+
    geom_histogram(bins = para_nfeature[["bins"]],
                   linewidth = 0.9,
                   fill = viridis::plasma(1),
                   color = "white") +
    scale_y_continuous(breaks = para_nfeature[["y_tick"]]) +
    scale_x_continuous(breaks = para_nfeature[["x_tick"]]) +
    coord_cartesian(xlim = para_nfeature[["x_width"]],
                    ylim = para_nfeature[["y_high"]]) +
    ICHMousewch:::.plotting_theme()

  barchart <- barchart_ncount + barchart_nfeature +
    plot_layout(ncol = 2) +
    plot_annotation(title = plotting_title,
                    theme = ICHMousewch:::.plotting_theme())

  return(patchworkGrob(barchart))

}

#' plotting theme

.plotting_theme <- function() {

  theme(plot.title = element_text(family = "Arial",size = 12,color = "black",face = "bold",hjust = 0.5,vjust = 0.5,margin = margin(b = 10, t = 10)),
        axis.text.x = element_text(family = "Arial",size = 8,color = "black",hjust = 0.5,vjust = 0.5),
        axis.text.y = element_text(family = "Arial",size = 8,color = "black",hjust = 0.5,vjust = 0.5),
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
        panel.grid.minor = element_blank())

}

#' create violin plot of quality control
#'
#' @param seu_metadata the Seurat Object metadata
#' @param filtered_barcodes the filtered barcodes
#' @param plotting_var the variable for plotting

.create_violin_plot_of_quality_control <- function(seu_metadata,filtered_barcodes,plotting_var) {

  on.exit(gc())

  unfiltered_seu_metadata <- seu_metadata[tissue == 1]
  filtered_seu_metadata <- seu_metadata[barcode %in% filtered_barcodes & tissue == 1]

  unfiltered_violin_da <- unfiltered_seu_metadata[,..plotting_var]
  unfiltered_seu_metadata[,violin_var := unfiltered_violin_da]

  filtered_violin_da <- filtered_seu_metadata[,..plotting_var]
  filtered_seu_metadata[,violin_var := filtered_violin_da]

  unfiltered_violin <- ggplot(data = unfiltered_seu_metadata,mapping = aes(x = orig.ident, y = violin_var)) +
    geom_jitter(alpha = 0.3,size = 0.5,color = "#35B779",width = 0.5) +
    geom_violin(alpha = 0,width = 1) +
    scale_x_discrete(expand = expansion(0)) +
    scale_y_continuous(expand = expansion(0)) +
    ICHMousewch:::.plotting_theme() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())

  filtered_violin <- ggplot(data = filtered_seu_metadata,mapping = aes(x = orig.ident, y = violin_var)) +
    geom_jitter(alpha = 0.3,size = 0.5,color = "#35B779",width = 0.5) +
    geom_violin(alpha = 0,width = 1) +
    scale_x_discrete(expand = expansion(0)) +
    scale_y_continuous(expand = expansion(0)) +
    ICHMousewch:::.plotting_theme() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())

  co_violin <- unfiltered_violin + filtered_violin +
    plot_layout(ncol = 2) +
    plot_annotation(title = plotting_var,
                    theme = ICHMousewch:::.plotting_theme())

  return(patchworkGrob(co_violin))

}

#' generate bubble dataset
#'
#' @param gene_info the gene information
#' @param GO_info the GO term information

.generate_bubble_dataset <- function(gene_info,GO_info,gene_num) {

  on.exit(gc())

  for (i in 1:length(GO_info)) {

    setorder(GO_info[[i]],-Count)

  }

  GO_results <- rbindlist(GO_info)

  setorder(gene_info,-GO_term_Count)
  gene_ls <- gene_info[c(1:gene_num),gene_name]

  diff_da <- gene_info[gene_name %in% gene_ls]

  diff_order <- match(gene_ls,diff_da[,gene_name])

  diff_da <- diff_da[diff_order,avg_log2FC]

  GO_id <- GO_results[,ID]
  GO_num <- length(GO_id)

  posi_y_da <- data.table(GO_ID = GO_id,
                          GO_group = character(GO_num),
                          bubble_y = integer(GO_num),
                          y_symbol = rep("posi_y",times = GO_num),
                          avg_log2FC = rep(NA,times = GO_num),
                          gene_name = rep(NA,times = GO_num),
                          bubble_color = rep(NA,times = GO_num))
  GO_order <- match(GO_id,GO_results[,ID])
  posi_y_da[,bubble_y := GO_results[GO_order,Count]]

  negt_y_da <- data.table(GO_ID = rep(GO_id,times = gene_num),
                          GO_group = character(GO_num*gene_num),
                          bubble_y = rep(seq(from = -5, to = -max(GO_results[,Count]), length.out = gene_num),each = GO_num),
                          y_symbol = rep("negt_y",times = GO_num*gene_num),
                          avg_log2FC = rep(diff_da,each = GO_num),
                          gene_name = rep(gene_ls,each = GO_num),
                          bubble_color = character(GO_num))

  bubble_dataset <- rbindlist(list(posi_y_da,negt_y_da))

  group_na <- names(GO_info)

  k <- 0
  for (i in 1:length(group_na)) {

    group_GO_id <- GO_info[[group_na[i]]][,ID]

    for (j in 1:length(group_GO_id)) {

      k <- k+1

      gene_na_ls <- GO_info[[group_na[i]]][ID %in% group_GO_id[j],gene][[1]]

      bubble_dataset[GO_ID %in% group_GO_id[j],GO_group := group_na[i]]

      selected_gene_na <- bubble_dataset[!is.na(gene_name) & gene_name %in% gene_na_ls,gene_name]
      unselected_gene_na <- bubble_dataset[!gene_name %in% selected_gene_na,gene_name]

      bubble_dataset[GO_ID %in% group_GO_id[j] & !is.na(gene_name) & gene_name %in% selected_gene_na,bubble_color := viridis::plasma(length(GO_id))[k]]
      bubble_dataset[GO_ID %in% group_GO_id[j] & !is.na(gene_name) & gene_name %in% unselected_gene_na,bubble_color := "#000000FF"]

    }

  }

  return(bubble_dataset)

}

#' generate bubble chart y tick
#'
#' @param plotting_dataset the plotting dataset

.generate_bubble_chart_y_tick <- function(plotting_dataset) {

  on.exit(gc())

  gene_ls <- plotting_dataset[!is.na(gene_name),gene_name] %>%
    unique()

  negt_y_tick <- plotting_dataset[y_symbol == "negt_y",bubble_y] %>%
    unique()
  names(negt_y_tick) <- gene_ls

  max_count <- max(plotting_dataset[y_symbol == "posi_y",bubble_y]) %>%
    as.character()

  max_count_dig <- strsplit(max_count,split = "")[[1]]

  if (length(max_count_dig) >= 2) {

    tick_step <- rep("0",times = length(max_count_dig)-1)

    tick_step[1] <- "5"

    tick_step <- paste(tick_step,collapse = "") %>%
      as.numeric()

  }

  posi_y_tick <- seq(from = 0 ,to = max_count, by = tick_step)
  names(posi_y_tick) <- as.character(posi_y_tick)

  y_tick <- c(posi_y_tick,negt_y_tick)

  return(y_tick)

}

#' plotting GO term cluster heatmap
#'
#' @param ich_mouse the class of ICH_Mouse

.plotting_GO_term_cluster_heatmap <- function(ich_mouse,diff_expr_symbol = "edge-normal") {

  on.exit(gc())

  GO_id_group <- ich_mouse@GO_ID_group
  group_na <- names(GO_id_group)
  GO_enrichment_result <- ich_mouse@GO_enrichment[[diff_expr_symbol]]

  GO_enrichment_result[,log2_fold_enrichment := log2(FoldEnrichment)]

  heatmap_ls <- vector("list",length = length(group_na))
  names(heatmap_ls) <- group_na

  for (i in 1:length(group_na)) {

    GO_id_ls <- GO_id_group[[group_na[i]]]

    sub_GO_enrichment_result <- GO_enrichment_result[ID %in% GO_id_ls]

    gene_ls <- sub_GO_enrichment_result[,geneID] %>%
      unlist() %>%
      paste(collapse = "/") %>%
      strsplit(split = "/")
    gene_ls <- unique(gene_ls[[1]])

    id_order <- match(GO_id_ls,sub_GO_enrichment_result[,ID])

    heatmap_dataset <- data.table(GO_ID = rep(GO_id_ls,times = length(gene_ls)),
                                  gene_name = rep(gene_ls,each = length(GO_id_ls)),
                                  log2fold_enrichment = sub_GO_enrichment_result[id_order,log2_fold_enrichment],
                                  relation = character(length(GO_id_ls)*length(gene_ls)))

    for (j in 1:length(GO_id_ls)) {

      genes <- GO_enrichment_result[ID %in% GO_id_ls[j],geneID] %>%
        strsplit(split = "/")
      genes <- genes[[1]]

      heatmap_dataset[GO_ID %in% GO_id_ls[j] & gene_name %in% genes,relation := "yes"]
      heatmap_dataset[GO_ID %in% GO_id_ls[j] & !gene_name %in% genes,relation := "no"]

    }

    heatmap_dataset[relation == "no",log2fold_enrichment := 0]

    GO_heatmap <- ggplot() +
      geom_tile(data = heatmap_dataset,
                mapping = aes(x = GO_ID,y = gene_name,fill = log2fold_enrichment)) +
      scale_fill_gradient(low = "white",high = "red")

    heatmap_ls[group_na] <- list(GO_heatmap)

  }

  return(heatmap_ls)

}

