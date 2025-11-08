
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

  in_tissue_count_matrix <- raw_count_matrix[gene_ls,seu_metadata_with_cluster_symbol[,barcode],drop = FALSE]

  if (gradient) {

    edge_barcodes <- in_tissue_metadata[!center_edge_symbol == "1",barcode]
    normal_barcodes <- in_tissue_metadata[center_edge_symbol == "1",barcode]

    logi_zero <- !in_tissue_count_matrix == 0
    logi_edge <-  Matrix::Matrix(rep(colnames(in_tissue_count_matrix) %in% edge_barcodes,times = length(gene_ls)),
                                 nrow = length(gene_ls),
                                 byrow = TRUE,
                                 sparse = TRUE)
    logi_normal <- Matrix::Matrix(rep(colnames(in_tissue_count_matrix) %in% normal_barcodes,times = length(gene_ls)),
                                  nrow = length(gene_ls),
                                  byrow = TRUE,
                                  sparse = TRUE)

    in_tissue_count_matrix[logi_zero & logi_edge] <- 1
    in_tissue_count_matrix[logi_zero & logi_normal] <- 2

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

  GO_ID_group <- ich_mouse@GO_ID_group

  group_na <- names(GO_ID_group)

  GO_ID_pair <- vector("list",length = length(GO_ID_group))
  names(GO_ID_pair) <- group_na

  for (i in 1:length(group_na)) {

    id_pair <- combn(GO_ID_group[[group_na[i]]],2) %>%
      as.data.table()

    pair_da_1 <- data.table(GO_id_x = unlist(id_pair[1]),
                            GO_id_y = unlist(id_pair[2]),
                            GO_id_group = rep(group_na[i],times = ncol(id_pair)))

    pair_da_2 <- data.table(GO_id_x = unlist(id_pair[2]),
                            GO_id_y = unlist(id_pair[1]),
                            GO_id_group = rep(group_na[i],times = ncol(id_pair)))

    GO_ID_pair[group_na[i]] <- list(pair_da_1,pair_da_2) %>%
      rbindlist() %>%
      list()

  }

  special_blocks <- rbindlist(GO_ID_pair) %>%
    as.data.frame() %>%
    mutate(GO_id_x = factor(GO_id_x,levels = unique(GO_id_x)),
           GO_id_y = factor(GO_id_y,levels = unique(GO_id_y)),
           GO_id_group = factor(GO_id_group,levels = unique(GO_id_group)))

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

  sim_heatmap <- ggplot() +
    geom_tile(data = sim_dt,mapping = aes(x = GO_id_x,y = GO_id_y,fill = similarity)) +
    geom_tile(data = special_blocks,
              mapping = aes(x = GO_id_x,y = GO_id_y,colour = GO_id_group),
              fill = "white",
              alpha = 0,
              linetype = 1) +
    scale_fill_gradientn(colours = c("#FFFFFF", "#FFA500", "#E6550D"),
                         values = c(0,0.75,1),
                         limits = c(0,1)) +
    ICHMousewch:::.plotting_theme() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()) +
    guides(color = guide_legend(position = "top",
                                nrow = 2,
                                theme = theme(legend.text = element_text(size = 12,
                                                                         family = "Arial",
                                                                         vjust = 0.5,
                                                                         hjust = 0),
                                              legend.title = element_blank(),
                                              legend.key.size = unit(12,"pt"))),
           fill = guide_none())

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
    labs(title = plotting_title) +
    ICHMousewch:::.plotting_theme() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())

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
    labs(title = plotting_title) +
    ICHMousewch:::.plotting_theme() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())

  ncount_na <- paste("Barchart_Count",plotting_title,sep = "-")
  nfeature_na <- paste("Barchart_Feature",plotting_title,sep = "-")

  barchart_ls <- list()
  barchart_ls[ncount_na] <- list(ggplotGrob(barchart_ncount))
  barchart_ls[nfeature_na] <- list(ggplotGrob(barchart_nfeature))

  return(barchart_ls)

}

#' plotting theme

.plotting_theme <- function() {

  theme(plot.title = element_text(family = "Arial",size = 12,color = "black",face = "bold",hjust = 0.5,vjust = 0.5,margin = margin(b = 10, t = 10)),
        axis.text.x = element_text(family = "Arial",size = 8,color = "black",hjust = 0.5),
        axis.text.y = element_text(family = "Arial",size = 8,color = "black",hjust = 0.5),
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

#' create umap plot
#'
#' @param raw_count_matrix the raw count matrix
#' @param seu_meta the seurat metadata
#' @param umap_feature the features for umap
#' @param diff_expr_genes the different expression genes
#' @param filter_genes the filter genes
#' @param group_symbol the group symbol

.create_umap_plot <- function(raw_count_matrix,seu_meta,diff_expr_genes,filter_genes,group_symbol = "center_edge_symbol",umap_feature = "variable_genes") {

  on.exit(gc())

  fil_count_matrix <- raw_count_matrix[filter_genes,seu_meta[,barcode]]

  seu_obj <- CreateSeuratObject(fil_count_matrix) %>%
    NormalizeData()

  if (group_symbol == "center_edge_symbol") {

    seu_meta[center_edge_symbol == "1", cluster := "normal"]
    seu_meta[center_edge_symbol == "2", cluster := "center"]
    seu_meta[center_edge_symbol == "3", cluster := "edge"]

  }

  if (group_symbol == "hematoma_symbol") {

    seu_meta[hematoma_symbol == "1", cluster := "normal"]
    seu_meta[hematoma_symbol == "2", cluster := "hematoma"]

  }

  if (group_symbol == "GMM_cluster") {

    seu_meta[,cluster := GMM_cluster]

  }

  if (group_symbol == "Louvain_cluster_posi") {

    seu_meta[,cluster := Louvain_cluster_posi]

  }

  if (group_symbol == "Louvain_cluster_filt_gene") {

    seu_meta[,cluster := Louvain_cluster_filt_gene]

  }

  if (umap_feature == "diff_expr_genes") {

    fil_genes <- diff_expr_genes

    seu_obj <- seu_obj[fil_genes]

  }

  if (umap_feature == "filter_genes") {

    fil_genes <- filter_genes
    seu_obj <- seu_obj[fil_genes]

  }

  if (umap_feature == "variable_genes") {

    variable_genes <- seu_obj %>%
      FindVariableFeatures() %>%
      VariableFeatures()

    fil_genes <- variable_genes

    seu_obj <- seu_obj[fil_genes]

  }

  seu_obj <- seu_obj %>%
    NormalizeData() %>%
    ScaleData() %>%
    RunPCA(features = fil_genes,
           npcs = 10) %>%
    RunUMAP(dims = 1:10)

  cluster_order <- match(rownames(seu_obj@meta.data),seu_meta[,barcode])

  seu_obj@meta.data$cluster <- seu_meta[cluster_order,cluster]

  umap_plot <- DimPlot(seu_obj,
                       group.by = "cluster") %>%
    patchworkGrob()

  return(umap_plot)

}

#' create gene distribution map
#'
#' @param seu_meta the Seurat metadata
#' @param raw_count_matrix the raw count matrix
#' @param filtered_genes the filtered genes
#' @param filtered_barcodes the filtered barcodes
#' @param diff_expr_gene the differential expression genes
#' @param background_genes the background genes for graph
#' @param aim_gene the aim gene for plotting

.create_gene_distribution_map <- function(seu_meta,raw_count_matrix,filtered_genes,filtered_barcodes,diff_expr_gene,aim_gene = character(),background_genes = c("Hbb-bt","Hbb-bs","Hba-a2"),point_size = 0.0001) {

  on.exit(gc())

  seu_meta <- seu_meta[cell_ID %in% filtered_barcodes]

  filtered_count_matrix <- raw_count_matrix[filtered_genes,seu_meta[,cell_ID]]

  seu_meta[,plot_row := -row]
  seu_meta[,plot_col := -col]

  seu_obj <- CreateSeuratObject(counts = filtered_count_matrix) %>%
    NormalizeData() %>%
    ScaleData(features = diff_expr_gene[avg_log2FC > 1,gene_name])

  gene_expr_matrix <- t(seu_obj@assays$RNA$scale.data) %>%
    as.data.frame()
  gene_expr_matrix <- bind_cols(data.frame(barcode = rownames(gene_expr_matrix)),
                                gene_expr_matrix) %>%
    as.data.table()

  graph_df <- merge(seu_meta,gene_expr_matrix,by = "barcode")

  initialize_color_symbol <- function(background_genes,filtered_count_matrix,graph_df) {

    graph_df[,color_symbol := "0"]

    if (length(background_genes) == 1) {

      background_df <- data.frame(background_genes_name = filtered_count_matrix[background_genes,])

    } else {

      background_df <- filtered_count_matrix[background_genes,] %>%
        as.data.frame() %>%
        t()

    }

    graph_df[,background_count := Matrix::rowSums(background_df != 0)]
    graph_df[background_count != 0,color_symbol := "1"]

    return(graph_df)

  }


  if (length(aim_gene) != 0) {

    distribution_graph_ls <- vector("list",length = length(aim_gene))
    names(distribution_graph_ls) <- aim_gene

    for (i in 1:length(aim_gene)) {

      graph_df <- initialize_color_symbol(background_genes = background_genes,
                                          filtered_count_matrix = filtered_count_matrix,
                                          graph_df = graph_df)

      aim_gene_na <- aim_gene[i]

      aim_gene_df <- data.frame(aim_gene_name = filtered_count_matrix[aim_gene_na,])

      aim_gene_scale_data <- graph_df[,..aim_gene_na]
      graph_df[,plotting_gene := aim_gene_scale_data]

      graph_df[,aim_gene_count := Matrix::rowSums(aim_gene_df != 0)]

      graph_df[aim_gene_count != 0,color_symbol := "2"]

      distribution_graph <- ggplot() +
        geom_point(data = graph_df[color_symbol == "0"],
                   mapping = aes(x = plot_row,y = plot_col,colour = color_symbol),
                   size = point_size,
                   color = "gray90") +
        geom_point(data = graph_df[color_symbol == "1"],
                   mapping = aes(x = plot_row,y = plot_col,colour = color_symbol),
                   size = point_size,
                   color = "white") +
        geom_point(data = graph_df[color_symbol == "2"],
                   mapping = aes(x = plot_row,y = plot_col,colour = plotting_gene),
                   size = point_size) +
        scale_color_gradientn(colors = c("#FEF4E8", "#FED9A6", "#FEB24C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026"),
                              limits = range(graph_df[color_symbol == "2",plotting_gene])) +
        coord_fixed() +
        labs(title = aim_gene_na) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = "right",
              legend.title = element_blank(),
              legend.text = element_text(size = 8,
                                         family = "Arial"),
              panel.background = element_rect(fill = "white", color = "black"),
              plot.background = element_rect(fill = "white", color = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 12,
                                        face = "bold",
                                        family = "Arial",
                                        hjust = 0.5,
                                        vjust = 0,
                                        margin = margin(b = 10)))

      distribution_graph_ls[aim_gene_na] <- list(ggplotGrob(distribution_graph))

    }

  } else {

    graph_df <- initialize_color_symbol(background_genes = background_genes,
                                        filtered_count_matrix = filtered_count_matrix,
                                        graph_df = graph_df)

    distribution_graph <- ggplot() +
      geom_point(data = graph_df[color_symbol == "0"],
                 mapping = aes(x = plot_row,y = plot_col,colour = color_symbol),
                 size = 0.01) +
      scale_colour_manual(values = c("0" = "grey")) +
      geom_point(data = graph_df[color_symbol == "1"],
                 mapping = aes(x = plot_row,y = plot_col,colour = color_symbol),
                 size = 0.01,
                 color = "white") +
      labs(title = "Background") +
      coord_fixed() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "right",
            legend.title = element_blank(),
            legend.text = element_text(size = 8,
                                       family = "Arial",
                                       color = "white"),
            legend.key = element_blank(),
            legend.background = element_blank(),
            legend.box.background = element_blank(),
            panel.background = element_rect(fill = "white", color = "black"),
            plot.background = element_rect(fill = "white", color = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 12,
                                      face = "bold",
                                      family = "Arial",
                                      hjust = 0.5,
                                      vjust = 0,
                                      margin = margin(b = 10))) +
      guides(color = guide_legend(override.aes = list(colour = "white")))

    background_na <- background_genes %>%
      paste(collapse = "-") %>%
      paste("background",sep = "_")
    distribution_graph_ls <- list()
    distribution_graph_ls[background_na] <- list(ggplotGrob(distribution_graph))

  }

  return(distribution_graph_ls)

}

#' create GO cluster graph combination set
#'
#' @param GO_ID_group the GO ID group
#' @param GO_enrichment the result of GO enrichment
#' @param seu_meta the Seurat metadata
#' @param raw_count_matrix the raw count matrix
#' @param filtered_genes the filtered genes
#' @param filtered_barcodes the filtered barcodes
#' @param diff_expr_gene the differential expression genes

.create_GO_cluster_graph_combination_set <- function(GO_ID_group,GO_enrichment,seu_meta,raw_count_matrix,filtered_genes,filtered_barcodes,diff_expr_gene) {

  on.exit(gc())

  sub_GO_enri <- GO_enrichment[ID %in% GO_ID_group]

  gene_ls <- vector("list",length = length(GO_ID_group))
  names(gene_ls) <- GO_ID_group

  for (i in 1:length(GO_ID_group)) {

    genes <- sub_GO_enri[ID == GO_ID_group[i],geneID] %>%
      strsplit(split = "/")
    gene_ls[GO_ID_group[i]] <- genes

  }

  aim_genes <- gene_ls %>%
    unlist() %>%
    unique()

  gene_distribution_map <- ICHMousewch:::.create_gene_distribution_map(seu_meta = seu_meta,
                                                                       raw_count_matrix = raw_count_matrix,
                                                                       filtered_genes = filtered_genes,
                                                                       filtered_barcodes = filtered_barcodes,
                                                                       diff_expr_gene = diff_expr_gene,
                                                                       aim_gene = aim_genes)

}

#' create Venn Diagram
#'
#' @param data_list the data list for Venn Diagram

.create_Venn_Diagram <- function(data_list) {

  on.exit(gc())

  venn_plot <- venn.diagram(x = data_list,
                            filename = NULL)

  return(venn_plot)

}

#' create count distribution map
#'
#' @param seu_meta the Seurat metadata

.create_count_distribution_map <- function(seu_meta,point_size = 0.01) {

  on.exit(gc())

  seu_meta[,plot_row := -row]
  seu_meta[,plot_col := -col]

  count_distribution_map <- ggplot() +
    geom_point(data = seu_meta,
               mapping = aes(x = plot_row,y = plot_col,color = nCount_log2),
               size = point_size) +
    scale_color_gradientn(colors = c("#0D0887", "#3F007D", "#6A00A8", "#B12A90", "#E16462", "#FCA636", "#F0F921"),
                          limits = range(seu_meta[,nCount_log2])) +
    labs(title = "log2Count") +
    coord_fixed() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 8,
                                     family = "Arial"),
          panel.background = element_rect(fill = "white", color = "black"),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 12,
                                    face = "bold",
                                    family = "Arial",
                                    hjust = 0.5,
                                    vjust = 0,
                                    margin = margin(b = 10)))

  return(count_distribution_map)

}

#' create the barchart of the mouse cell reference dataset
#'
#' @param gene_ls the gene list for the barchart

.create_the_barchart_of_the_mouse_cell_reference_dataset <- function(gene_ls) {

  on.exit(gc())

  reference_dataset <- ICHMousewch::mouse_cell_plotting_dataset

  gene_ls <- gene_ls[gene_ls %in% colnames(reference_dataset)]

  cell_type_color <- qualitative_hcl(length(unique(reference_dataset[,cell_type])), palette = "Dark 3")
  names(cell_type_color) <- unique(reference_dataset[,cell_type])

  barchart_ls <- vector("list",length = length(gene_ls))
  names(barchart_ls) <- paste(gene_ls,"mouse_cell",sep = "-")
  for (i in 1:length(gene_ls)) {

    aim_gene_na <- gene_ls[i]
    avg_aim_gene_na <- paste("avg",aim_gene_na,sep = "-")

    orig_da_col <- c("cell_type",aim_gene_na,avg_aim_gene_na)
    orig_da <- reference_dataset[,..orig_da_col]

    aim_gene_da <- orig_da[,..aim_gene_na] %>%
      unlist()
    avg_aim_gene_da <- orig_da[,..avg_aim_gene_na] %>%
      unlist()

    adj_da <- data.table(cell_type = orig_da[,cell_type],
                         aim_gene = aim_gene_da,
                         avg_aim_gene = avg_aim_gene_da)

    plot_da <- adj_da[,c("cell_type","avg_aim_gene")] %>%
      unique() %>%
      setorder(-avg_aim_gene) %>%
      as.data.frame() %>%
      mutate(cell_type = factor(cell_type,levels = cell_type))

    cell_type_num <- seq(from = 1,
                         to = 20,
                         by = 1)

    plotting_cell_type <- plot_da$cell_type[cell_type_num]

    barchart <- ggplot() +
      geom_col(data = plot_da[cell_type_num,],
               mapping = aes(x = cell_type,y = avg_aim_gene,
                             colour = cell_type,
                             fill = cell_type),
               width = 0.8,
               alpha = 0.2) +
      geom_point(data = adj_da[cell_type %in% plotting_cell_type],
                 mapping = aes(x = cell_type,
                               y = aim_gene,
                               colour = cell_type),
                 position = position_jitter(width = 0.3,
                                            seed = 2025),
                 size = 1,
                 shape = 24,
                 fill = NA)  +
      geom_boxplot(data = adj_da[cell_type %in% plotting_cell_type],
                   mapping = aes(x = cell_type,
                                 y = aim_gene,
                                 colour = cell_type),
                   width = 0.3,
                   fill = NA,
                   outlier.shape = NA,
                   linewidth = 0.8) +
      scale_colour_manual(values = cell_type_color) +
      scale_fill_manual(values = cell_type_color) +
      labs(title = gene_ls[i]) +
      ICHMousewch:::.plotting_theme() +
      theme(panel.background = element_rect(fill = "white", color = "black"),
            plot.background = element_rect(fill = "white", color = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(family = "Arial",
                                       size = 8,
                                       color = "black",
                                       hjust = 1),
            axis.text.x = element_text(family = "Arial",
                                       size = 8,
                                       color = "black",
                                       hjust = 1,
                                       angle = 30),
            plot.margin = margin(l = 1,
                                 unit = "cm"))

    barchart_na <- paste(gene_ls[i],"mouse_cell",sep = "-")
    barchart_ls[barchart_na] <- list(ggplotGrob(barchart))

  }

  return(barchart_ls)

}
