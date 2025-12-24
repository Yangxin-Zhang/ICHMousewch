
# R/utils-graph.R

#' create Giotto Instruction
#'
#' @param python_path the path to a python which the Giotto can use
#' @param results_folder the folder for Giotto save plots

.create_giotto_instruction <- function(python_path,
                                       results_folder) {

  on.exit(gc())

  instrs = createGiottoInstructions(save_dir = results_folder,
                                    save_plot = FALSE,
                                    show_plot = FALSE,
                                    python_path = python_path)

  return(instrs)

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
#' @param cluster_symbol the GO cluster symbol
#' @param special_block logic
#' @param GO_group_symbol the GO group symbol

.plotting_GO_term_similarity_heatmap <- function(ich_mouse,cluster_symbol = "edge-normal",special_block = FALSE,GO_group_symbol = list()) {

  on.exit(gc())

  group_color <- scales::hue_pal()(length(GO_group_symbol))
  names(group_color) <- GO_group_symbol

  if (special_block) {

    if (length(GO_group_symbol) == 0) {

      GO_cluster <- ich_mouse@GO_cluster[[cluster_symbol]]
      clu_na <- names(GO_cluster)
      GO_ID_group <- vector("list",length = length(GO_cluster))
      names(GO_ID_group) <- clu_na

      for (i in 1:length(clu_na)) {

        GO_ID_group[clu_na[i]] <- list(GO_cluster[[clu_na[i]]][,ID])

      }

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

        pair_da_self <- data.table(GO_id_x = GO_ID_group[[group_na[i]]],
                                   GO_id_y = GO_ID_group[[group_na[i]]],
                                   GO_id_group = rep(group_na[i],times = length(GO_ID_group[[group_na[i]]])))

        GO_ID_pair[group_na[i]] <- list(pair_da_1,pair_da_2,pair_da_self) %>%
          rbindlist() %>%
          list()

      }

      special_blocks <- rbindlist(GO_ID_pair) %>%
        as.data.frame() %>%
        mutate(GO_id_x = factor(GO_id_x,levels = unique(GO_id_x)),
               GO_id_y = factor(GO_id_y,levels = unique(GO_id_y)),
               GO_id_group = factor(GO_id_group,levels = unique(GO_id_group)))

    } else {

      GO_ID_group <- ich_mouse@GO_ID_group[GO_group_symbol]

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

        pair_da_self <- data.table(GO_id_x = GO_ID_group[[group_na[i]]],
                                   GO_id_y = GO_ID_group[[group_na[i]]],
                                   GO_id_group = rep(group_na[i],times = length(GO_ID_group[[group_na[i]]])))

        GO_ID_pair[group_na[i]] <- list(pair_da_1,pair_da_2,pair_da_self) %>%
          rbindlist() %>%
          list()

      }

      special_blocks <- rbindlist(GO_ID_pair) %>%
        as.data.frame() %>%
        mutate(GO_id_x = factor(GO_id_x,levels = unique(GO_id_x)),
               GO_id_y = factor(GO_id_y,levels = unique(GO_id_y)),
               GO_id_group = factor(GO_id_group,levels = unique(GO_id_group)))

    }

  }

  if (length(GO_group_symbol) != 0) {

    GO_results <- ich_mouse@GO_cluster[[cluster_symbol]] %>%
      ICHMousewch:::.adjust_special_block_GO_group_order(GO_group = GO_ID_group) %>%
      rbindlist()

  } else {

    GO_results <- rbindlist(ich_mouse@GO_cluster[[cluster_symbol]])

  }

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

  color_set <- c("#92C5DE","#D1E5F0","#F7F7F7",
                 "#FFF7BC","#FEE391", "#FEC44F",
                 "#FB9A29", "#EC7014", "#CC4C02",
                 "#993404", "#662506")
  color_break <- seq(from = min(sim_data),
                     to = max(sim_data),
                     length.out = length(color_set))

  sim_dt <- as.data.frame(sim_dt) %>%
    mutate(GO_id_x = factor(GO_id_x,levels = colnames(sim_mat)),
           GO_id_y = factor(GO_id_y,levels = rownames(sim_mat)))

  if (special_block) {

    sim_heatmap <- ggplot() +
      geom_tile(data = sim_dt,mapping = aes(x = GO_id_x,y = GO_id_y,fill = similarity)) +
      scale_fill_gradientn(colours = color_set,
                           values = scales::rescale(color_break,to = c(0,1))) +
      geom_tile(data = special_blocks,
                mapping = aes(x = GO_id_x,y = GO_id_y,color = GO_id_group),
                fill = "white",
                alpha = 0) +
      scale_color_manual(values = group_color) +
      labs(title = "GO Similarity Matrix",
           x = "GO IDs",
           y = "GO IDs",
           color = "GO Group",
           fill = "Similarity") +
      coord_fixed() +
      ICHMousewch:::.plotting_theme() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            legend.title = element_text(size = 12)) +
      guides(fill = guide_colorbar(position = "left",
                                   theme = theme(legend.key.size = unit(13,"pt"),
                                                 legend.text = element_text(size = 10))),
             color = guide_legend(position = "right",
                                  theme = theme(legend.key.size = unit(10,"pt")),
                                  override.aes = list(fill = "white",
                                                      linewidth = 1)))

  } else {

    sim_heatmap <- ggplot() +
      geom_tile(data = sim_dt,mapping = aes(x = GO_id_x,y = GO_id_y,fill = similarity)) +
      scale_fill_gradientn(colours = color_set,
                           values = scales::rescale(color_break,to = c(0,1))) +
      labs(title = "GO Similarity Matrix",
           x = "GO IDs",
           y = "GO IDs",
           fill = "Similarity") +
      coord_fixed(ratio = 1) +
      ICHMousewch:::.plotting_theme() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank()) +
      guides(fill = guide_colorbar(position = "left",
                                   theme = theme(legend.key.size = unit(12,"pt"))))


  }


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
    labs(title = plotting_title,
         x = "counts",
         y = "n_spots") +
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
    labs(title = plotting_title,
         x = "n_features",
         y = "n_spots") +
    ICHMousewch:::.plotting_theme()

  ncount_na <- paste("Barchart_Count",plotting_title,sep = "-")
  nfeature_na <- paste("Barchart_Feature",plotting_title,sep = "-")

  barchart_ls <- list()
  barchart_ls[ncount_na] <- list(ggplotGrob(barchart_ncount))
  barchart_ls[nfeature_na] <- list(ggplotGrob(barchart_nfeature))

  return(barchart_ls)

}

#' plotting theme

.plotting_theme <- function()
  {

  theme(plot.title = element_text(family = "Arial",
                                  size = 12,
                                  color = "black",
                                  face = "bold",
                                  hjust = 0.5,
                                  vjust = 0,
                                  margin = margin(t = 5,
                                                  b = 5,
                                                  unit = "pt")),
        plot.margin = margin(r = 10,
                             l = 10,
                             b = 10,
                             t = 10,
                             unit = "pt"),
        axis.text.x = element_text(family = "Arial",size = 8,color = "black",hjust = 0.5),
        axis.text.y = element_text(family = "Arial",size = 8,color = "black",hjust = 0.5),
        axis.title.x = element_text(family = "Arial",
                                    size = 12,
                                    color = "black",
                                    hjust = 0.5,
                                    vjust = 0.5,
                                    margin = margin(b = 8,
                                                    t = 3,
                                                    unit = "pt")),
        axis.title.y = element_text(family = "Arial",
                                    size = 12,
                                    color = "black",
                                    hjust = 0.5,
                                    vjust = 0.5,
                                    margin = margin(r = 5,
                                                    l = 5,
                                                    unit = "pt")),
        legend.text = element_text(family = "Arial",
                                   size = 12,
                                   color = "black",
                                   hjust = 0.5,
                                   vjust = 0.5,
                                   margin = margin(r = 5,
                                                   l = 5,
                                                   unit = "pt")),
        legend.title = element_text(family = "Arial",
                                    size = 12,
                                    color = "black",
                                    hjust = 0.5,
                                    vjust = 0.5,
                                    margin = margin(b = 10,
                                                    unit = "pt")),
        legend.background = element_rect(fill = "white",
                                         color = NA),
        legend.key = element_rect(fill = "white",
                                  color = NA),
        legend.position = "right",
        panel.background = element_rect(fill = "white",
                                        color = "black"),
        plot.background = element_rect(fill = "white",
                                       color = NA),
        panel.grid.major = element_line(color = "gray90",
                                        linewidth = 0.2),
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
#' @param ich_mouse the class of ICH_Mouse
#' @param spaceranger_result the file address of spaceranger_result

.create_umap_plot <- function(ich_mouse,spaceranger_result,DEGs_dims = FALSE) {

  on.exit(gc())

  create_umap_plot <- function(umap_coord,cluster_symbol) {

    on.exit(gc())

    if (cluster_symbol == "hematoma") {

      umap_coord$cluster <-umap_coord$cluster_hematoma

      umap_title <- "Hematoma vs Normal Tissue"

      color_set <- c("center" = "#F8766D",
                     "edge" = "#00BA38",
                     "normal" = "#619CFF")

      umap_labels <- c("center","edge","normal")

    }

    if (cluster_symbol == "spaceranger") {

      umap_coord$cluster <-umap_coord$cluster_spaceranger

      umap_title <- "SpaceRanger Cluster"

      n_colors <- umap_coord$cluster %>%
        unique() %>%
        length()

      color_set <- scales::hue_pal()(n_colors)

      clu_sym <- seq(from = 1,
                     to = n_colors) %>%
        as.character()

      names(color_set) <- clu_sym

      umap_labels <- paste("cluster",clu_sym,sep = " ")

    }

    umap_plot <- ggplot() +
      geom_point(data = umap_coord,
                 mapping = aes(x = UMAP_1, y = UMAP_2, color = cluster),
                 size = 0.01) +
      scale_color_manual(values = color_set,
                         labels = umap_labels) +
      labs(title = umap_title,
           x = "UMAP-1",
           y = "UMAP-2") +
      ICHMousewch:::.plotting_theme() +
      theme(legend.position = "right",
            legend.title = element_blank(),
            panel.grid.major = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank()) +
      guides(color = guide_legend(theme = theme(legend.key.size = unit(10,"pt")),
                                  override.aes = list(size = 5)))

    return(umap_plot)

  }

  raw_count_matrix <- ich_mouse@raw_count_matrix
  seu_meta <- ich_mouse@seu_metadata_with_cluster_symbol
  diff_expr_genes <- ich_mouse@diff_expr_genes
  filter_genes <- ich_mouse@filtered_genes

  spaceranger_cluster <- read_csv(spaceranger_result[["spacerange_cluster_address"]]) %>%
    as_tibble()
  spaceranger_cluster <- tibble(barcode = spaceranger_cluster$Barcode,
                                cluster_spaceranger = as.character(spaceranger_cluster$Cluster))

  spaceranger_umap <- read_csv(spaceranger_result[["spacerange_umap_address"]]) %>%
    as_tibble()
  spaceranger_umap <- tibble(barcode = spaceranger_umap$Barcode,
                             UMAP_1 = spaceranger_umap$`UMAP-1`,
                             UMAP_2 = spaceranger_umap$`UMAP-2`)

  barcodes <- seu_meta[barcode %in% spaceranger_umap$barcode,barcode]


  seu_meta[center_edge_symbol == "1", cluster := "normal"]
  seu_meta[center_edge_symbol == "2", cluster := "center"]
  seu_meta[center_edge_symbol == "3", cluster := "edge"]

  hematoma_cluster <- seu_meta[barcode %in% barcodes,c("barcode","cluster")]
  hematoma_cluster <- tibble(barcode = hematoma_cluster$barcode,
                             cluster_hematoma = hematoma_cluster$cluster)

  spaceranger_umap <- merge(spaceranger_umap,hematoma_cluster,by = "barcode") %>%
    merge(spaceranger_cluster,by = "barcode") %>%
    column_to_rownames(var = "barcode")

  umap_hematoma_spaceranger <- create_umap_plot(umap_coord = spaceranger_umap,
                                                cluster_symbol = "hematoma")
  umap_spaceranger_cluster_spaceranger <- create_umap_plot(umap_coord = spaceranger_umap,
                                                           cluster_symbol = "spaceranger")

  umap_plot_ls <- list("umap_hematoma_spaceranger" = ggplotGrob(umap_hematoma_spaceranger),
                       "umap_spaceranger_cluster_spaceranger" = ggplotGrob(umap_spaceranger_cluster_spaceranger))

  if (DEGs_dims) {

    fil_count_matrix <- raw_count_matrix[filter_genes,barcodes]

    seu_obj <- CreateSeuratObject(fil_count_matrix) %>%
      NormalizeData()

    DEG_genes <- c(diff_expr_genes$`edge-normal`[avg_log2FC > 1 & p_val_adj < 0.01,gene_name],
                   diff_expr_genes$`center-edge`[abs(avg_log2FC) > 1 & p_val_adj < 0.01,gene_name],
                   diff_expr_genes$`center-normal`[avg_log2FC > 1 & p_val_adj < 0.01,gene_name]) %>%
      unique()

    seu_obj <- seu_obj[DEG_genes,]

    seu_obj <- ScaleData(seu_obj) %>%
      RunPCA(seed.use = 2025,
             npcs = 10,
             features = DEG_genes) %>%
      RunUMAP(dims = 1:10,
              min.dist = 1,
              seed.use = 2025)

    umap_coord_DEG <- Embeddings(seu_obj,"umap") %>%
      as.data.frame() %>%
      rownames_to_column()
    umap_coord_DEG <- tibble(barcode = umap_coord_DEG$rowname,
                             UMAP_1 = umap_coord_DEG$umap_1,
                             UMAP_2 = umap_coord_DEG$umap_2)

    umap_coord_DEG <- merge(umap_coord_DEG,hematoma_cluster,by = "barcode") %>%
      merge(spaceranger_cluster,by = "barcode") %>%
      column_to_rownames(var = "barcode")

    umap_hematoma_ichmouse <- create_umap_plot(umap_coord = umap_coord_DEG,
                                               cluster_symbol = "hematoma")

    umap_spaceranger_cluster_ichmouse <- create_umap_plot(umap_coord = umap_coord_DEG,
                                                          cluster_symbol = "spaceranger")

    umap_plot_ls <- append(umap_plot_ls,
                           list("umap_hematoma_ichmouse" = ggplotGrob(umap_hematoma_ichmouse),
                                "umap_spaceranger_cluster_ichmouse" = ggplotGrob(umap_spaceranger_cluster_ichmouse)))

  }

  return(umap_plot_ls)

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
#' @param ich_mouse the class of ICH_Mouse

.create_Venn_Diagram <- function(ich_mouse) {

  on.exit(gc())

  data_ls <- list("edge vs normal" = ich_mouse@diff_expr_genes$`edge-normal`[avg_log2FC > 1 & p_val_adj < 0.01,gene_name],
                  "center vs normal" = ich_mouse@diff_expr_genes$`center-normal`[avg_log2FC > 1 & p_val_adj < 0.01,gene_name],
                  "center vs edge" = ich_mouse@diff_expr_genes$`center-edge`[abs(avg_log2FC) > 1 & p_val_adj < 0.01,gene_name])

  venn_plot <- ggvenn::ggvenn(data_ls,
                              fill_color = c("#1F77B4","#FF7F0E","#2CA02C"),
                              fill_alpha = 0.8,
                              text_size = 3.5,
                              set_name_size = 4.5) +
    theme(text = element_text(family = "Arial"))

  venn_plot_ls <- list("DEG_venn_plot" = ggplotGrob(venn_plot))

  return(venn_plot_ls)

}

#' create the barchart of the mouse cell reference dataset
#'
#' @param gene_ls the gene list for the barchart
#' @param reference_data_symbol the reference data symbol

.create_the_barchart_of_the_mouse_cell_reference_dataset <- function(gene_ls,reference_data_symbol,cell_type_label = "main") {

  on.exit(gc())

  if (reference_data_symbol == "mouse_cell") {

    reference_dataset <- ICHMousewch::mouse_cell_plotting_dataset
    gene_ls <- gene_ls[gene_ls %in% colnames(reference_dataset)]

  }

  if (reference_data_symbol == "immune_cell") {

    reference_dataset <- ICHMousewch::mouse_immune_cell_plotting_dataset
    gene_ls <- gene_ls[gene_ls %in% colnames(reference_dataset)]

    gene_ls_na <- paste("avg",gene_ls,sep = "-")

    if (cell_type_label == "main") {

      reference_dataset[,cell_type := cell_type_main]

      gene_ls_na_main <- paste(gene_ls_na,"main",sep = "_")

      avg_da <- reference_dataset[,..gene_ls_na_main] %>%
        as.matrix()

      rownames(avg_da) <- reference_dataset[,dataset_symbol]
      colnames(avg_da) <- gene_ls_na

      avg_da <- bind_cols(data.frame(dataset_symbol = rownames(avg_da)),
                          avg_da) %>%
        as.data.table()

      reference_dataset <- merge(reference_dataset,avg_da,by = "dataset_symbol")

    }

    if (cell_type_label == "fine") {

      reference_dataset[,cell_type := cell_type_fine]

      gene_ls_na_fine <- paste(gene_ls_na,"fine",sep = "_")

      avg_da <- reference_dataset[,..gene_ls_na_fine] %>%
        as.matrix()

      rownames(avg_da) <- reference_dataset[,dataset_symbol]
      colnames(avg_da) <- gene_ls_na

      avg_da <- bind_cols(data.frame(dataset_symbol = rownames(avg_da)),
                          avg_da) %>%
        as.data.table()

      reference_dataset <- merge(reference_dataset,avg_da,by = "dataset_symbol")

    }

  }

  cell_type_color <- qualitative_hcl(length(unique(reference_dataset[,cell_type])), palette = "Dark 3")
  names(cell_type_color) <- unique(reference_dataset[,cell_type])

  if (reference_data_symbol == "immune_cell") {

    barchart_na <- paste(gene_ls,"mouse_immune_cell",sep = "-")

  }

  if (reference_data_symbol == "mouse_cell") {

    barchart_na <- paste(gene_ls,"mouse_cell",sep = "-")

  }

  barchart_ls <- vector("list",length = length(gene_ls))
  names(barchart_ls) <- barchart_na
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

    barchart_ls[barchart_na] <- list(ggplotGrob(barchart))

  }

  return(barchart_ls)

}

#' create heatmap of mouse cell reference dataset
#'
#' @param gene_ls the gene list for heatmap
#' @param reference_data_symbol the reference data symbol

.create_heatmap_of_mouse_cell_reference_dataset <- function(gene_ls,reference_data_symbol,cell_type_label = "main") {

  on.exit(gc())

  if (reference_data_symbol == "mouse_cell") {

    reference_dataset <- ICHMousewch::mouse_cell_plotting_dataset
    gene_ls <- gene_ls[gene_ls %in% colnames(reference_dataset)]

  }

  if (reference_data_symbol == "immune_cell") {

    reference_dataset <- ICHMousewch::mouse_immune_cell_plotting_dataset
    gene_ls <- gene_ls[gene_ls %in% colnames(reference_dataset)]

    gene_ls_na <- paste("avg",gene_ls,sep = "-")

    if (cell_type_label == "main") {

      reference_dataset[,cell_type := cell_type_main]

      gene_ls_na_main <- paste(gene_ls_na,"main",sep = "_")

      avg_da <- reference_dataset[,..gene_ls_na_main] %>%
        as.matrix()

      rownames(avg_da) <- reference_dataset[,dataset_symbol]
      colnames(avg_da) <- gene_ls_na

      avg_da <- bind_cols(data.frame(dataset_symbol = rownames(avg_da)),
                          avg_da) %>%
        as.data.table()

      reference_dataset <- merge(reference_dataset,avg_da,by = "dataset_symbol")

    }

    if (cell_type_label == "fine") {

      reference_dataset[,cell_type := cell_type_fine]

      gene_ls_na_fine <- paste(gene_ls_na,"fine",sep = "_")

      avg_da <- reference_dataset[,..gene_ls_na_fine] %>%
        as.matrix()

      rownames(avg_da) <- reference_dataset[,dataset_symbol]
      colnames(avg_da) <- gene_ls_na

      avg_da <- bind_cols(data.frame(dataset_symbol = rownames(avg_da)),
                          avg_da) %>%
        as.data.table()

      reference_dataset <- merge(reference_dataset,avg_da,by = "dataset_symbol")

    }

  }

  avg_gene_ls <- paste("avg",gene_ls,sep = "-")
  avg_gene_ls <- c("cell_type",avg_gene_ls)

  avg_dataset <- unique(reference_dataset[,..avg_gene_ls])

  plotting_dataset <- data.table(cell_type = rep(avg_dataset[,cell_type],times = length(gene_ls)),
                                 gene_name = rep(gene_ls,each = length(avg_dataset[,cell_type])),
                                 avg_expression = numeric(length(gene_ls)*length(avg_dataset[,cell_type])))

  for (i in 1:length(gene_ls)) {

    gene_na <- paste("avg",gene_ls[i],sep = "-")

    ce_order <- match(plotting_dataset[gene_name == gene_ls[i],cell_type],avg_dataset[,cell_type])

    plotting_dataset[gene_name == gene_ls[i],avg_expression := avg_dataset[ce_order,..gene_na]]

  }

  reference_heatmap <- ggplot() +
    geom_tile(data = plotting_dataset,
              mapping = aes(x = cell_type,
                            y = gene_name,
                            fill = avg_expression)) +
    scale_fill_gradientn(colors = c("#FEF4E8", "#FED9A6", "#FEB24C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026"),
                          limits = c(0,20)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(family = "Arial",
                                     size = 8,
                                     color = "black",
                                     hjust = 1,
                                     angle = 30),
          axis.text.y = element_text(family = "Arial",
                                     size = 8,
                                     color = "black",
                                     hjust = 1),
          panel.background = element_rect(fill = "white", color = "black"),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")

  return(reference_heatmap)

}

#' create gene group barchart
#'
#' @param ich_mouse the class of ICH_Mouse
#' @param gene_name the gene name

.create_gene_group_barchart <- function(ich_mouse,gene_name) {

  on.exit(gc())

  GO_result <- ich_mouse@GO_enrichment$`filtered_edge-normal`

  group_GO_result <- ICHMousewch:::.extract_group_GO_result(GO_result = GO_result,
                                             gene_set = gene_name)

  bubble_chart_ls <- list()
  for (i in 1:length(gene_name)) {

    if (length(group_GO_result[[gene_name[i]]]) != 0) {

      sub_GO_result <- GO_result[ID %in% group_GO_result[[gene_name[i]]]]
      sub_GO_result[,log10p := -log10(p.adjust)]
      setorder(sub_GO_result,-log10p)

      if (length(group_GO_result[[gene_name[i]]]) <= 15) {

        plotting_dataset <- sub_GO_result[,c("Description","log10p")] %>%
          setorder(log10p) %>%
          as_tibble() %>%
          mutate(Description = factor(Description,levels = Description))

      } else {

        plotting_dataset <- sub_GO_result[1:15,c("Description","log10p")] %>%
          setorder(log10p) %>%
          as_tibble() %>%
          mutate(Description = factor(Description,levels = Description))

      }

      bubble_chart <- ggplot() +
        geom_col(data = plotting_dataset,
                 mapping = aes(x = log10p,y = Description),
                 fill = "grey90",
                 color = "black") +
        ICHMousewch:::.plotting_theme() +
        theme(panel.grid.major = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())

      bubble_chart_ls[gene_name[i]] <- list(ggplotGrob(bubble_chart))

    }

  }


  return(bubble_chart_ls)

}

#' create GO result overview barchart
#'
#' @param ich_mouse the class of ICH_Mouse

.create_GO_result_overview_barchart <- function(ich_mouse) {

  on.exit(gc())

  filtered_GO <- ich_mouse@GO_enrichment$`filtered_total-diff_expr_genes`
  filtered_GO[,log10p := -log10(p.adjust)]
  setorder(filtered_GO,-FoldEnrichment)

  plotting_dataset <- filtered_GO[1:20] %>%
    setorder(FoldEnrichment) %>%
    as_tibble() %>%
    mutate(Description = factor(Description,levels = Description))

  bubble_chart <- ggplot() +
    geom_col(data = plotting_dataset,
             mapping = aes(x = FoldEnrichment,y = Description),
             fill = "grey90",
             color = "black",
             width = 0.7) +
    labs(x = "fold enrichment",
         title = "GO Enrichment Overview") +
    ICHMousewch:::.plotting_theme() +
    theme(panel.grid.major = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 8,
                                     hjust = 1),
          aspect.ratio = ((1/30)*nrow(plotting_dataset)))

  return(list("overview_GO_result_barchart" = ggplotGrob(bubble_chart)))

}

#' create GO result custom barchart
#'
#' @param ich_mouse the class of ICH_Mouse
#' @param GO_term_ls the list of GO term

.create_GO_result_custom_barchart <- function(ich_mouse,GO_term_ls) {

  on.exit(gc())

  plot_barchart <- function(ich_mouse,GO_term_set,GO_term_set_name) {

    on.exit(gc())

    filtered_GO <- ich_mouse@GO_enrichment$`filtered_total-diff_expr_genes`[ID %in% GO_term_set]
    filtered_GO[,log10p := -log10(p.adjust)]
    setorder(filtered_GO,-FoldEnrichment)

    plotting_dataset <- filtered_GO %>%
      setorder(FoldEnrichment) %>%
      as_tibble() %>%
      mutate(Description = factor(Description,levels = Description))

    bar_chart <- ggplot() +
      geom_col(data = plotting_dataset,
               mapping = aes(x = FoldEnrichment,y = Description),
               fill = "grey90",
               color = "black",
               width = 0.7) +
      labs(title = GO_term_set_name,
           x = "fold enrichment") +
      ICHMousewch:::.plotting_theme() +
      theme(panel.grid.major = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_text(size = 8,
                                       hjust = 1),
            aspect.ratio = ((1/30)*nrow(plotting_dataset)))

    return(bar_chart)

  }

  if (length(GO_term_ls) != 0) {

    barchart_name <- names(GO_term_ls)
    barchart_ls <- vector("list",length = length(barchart_name))
    names(barchart_ls) <- barchart_name

    for (i in 1:length(barchart_name)) {

      barchart_ls[barchart_name[i]] <- plot_barchart(ich_mouse = ich_mouse,
                                                     GO_term_set = GO_term_ls[[barchart_name[i]]],
                                                     GO_term_set_name = barchart_name[i]) %>%
        ggplotGrob() %>%
        list()

    }

  }

  return(barchart_ls)

}

#' create single gene information heatmap
#'
#' @param ich_mouse the class of ICH_Mouse
#' @param go_term_ls the GO term list
#' @param GO_symbol the GO analysis symbol

.create_single_gene_infotmation_heatmap <- function(ich_mouse,go_term_ls,GO_symbol = "filtered_total-diff_expr_genes") {

  on.exit(gc())

  go_results <- ich_mouse@GO_enrichment[[GO_symbol]]

  aim_go_results <- go_results[ID %in% go_term_ls] %>%
    ICHMousewch:::.split_GO_result_genes()

  gene_ls <- aim_go_results[,split_genes] %>%
    unlist() %>%
    unique()

  gene_info_df <- ICHMousewch:::.calculate_single_gene_average_expression(ich_mouse = ich_mouse,
                                                                          gene_ls = gene_ls) %>%
    ICHMousewch:::.add_go_term_condition_information(go_term_list = go_term_ls,
                                                     go_result = go_results) %>%
    setorder(-log_pt_edge_FC)

  plotting_df_log_pt_normal_FC <- gene_info_df[,c("gene_name","log_pt_normal_FC")]
  plotting_df_log_pt_normal_FC <- plotting_df_log_pt_normal_FC[,x_tick := "percent_normal"] %>%
    as_tibble() %>%
    rename(scale_dt = log_pt_normal_FC,
           y_tick = gene_name)

  plotting_df_log_pt_edge_FC <- gene_info_df[,c("gene_name","log_pt_edge_FC")]
  plotting_df_log_pt_edge_FC <- plotting_df_log_pt_edge_FC[,x_tick := "percent_edge"] %>%
    as_tibble() %>%
    rename(scale_dt = log_pt_edge_FC,
           y_tick = gene_name)

  plotting_df_log_pt_center_FC <- gene_info_df[,c("gene_name","log_pt_center_FC")]
  plotting_df_log_pt_center_FC <- plotting_df_log_pt_center_FC[,x_tick := "percent_center"] %>%
    as_tibble() %>%
    rename(scale_dt = log_pt_center_FC,
           y_tick = gene_name)

  plotting_df_log_expr_normal <- gene_info_df[,c("gene_name","log_expr_normal")]
  plotting_df_log_expr_normal <- plotting_df_log_expr_normal[,x_tick := "expression_normal"] %>%
    as_tibble() %>%
    rename(scale_dt = log_expr_normal,
           y_tick = gene_name)

  plotting_df_log_expr_center <- gene_info_df[,c("gene_name","log_expr_center")]
  plotting_df_log_expr_center <- plotting_df_log_expr_center[,x_tick := "expression_center"] %>%
    as_tibble() %>%
    rename(scale_dt = log_expr_center,
           y_tick = gene_name)

  plotting_df_log_expr_edge <- gene_info_df[,c("gene_name","log_expr_edge")]
  plotting_df_log_expr_edge <- plotting_df_log_expr_edge[,x_tick := "expression_edge"] %>%
    as_tibble() %>%
    rename(scale_dt = log_expr_edge,
           y_tick = gene_name)

  plotting_df_go_term_ls <- vector("list",length = length(go_term_ls))
  names(plotting_df_go_term_ls) <- go_term_ls
  for (i in 1:length(go_term_ls)) {

    go_description <- ich_mouse@GO_enrichment[[GO_symbol]][ID %in% go_term_ls[i],Description]

    col_chars <- c("gene_name","size_data")
    tool_df <- gene_info_df[,..col_chars]
    x_tick_chars <- rep(go_description,times = nrow(tool_df))
    tool_df <- tool_df[,x_tick := x_tick_chars]

    go_term_chars <- go_term_ls[i]
    color_dt_df <- gene_info_df[,..col_chars]
    color_dt_df <- color_dt_df[,color_dt_logi := gene_info_df[,..go_term_chars]]
    color_dt_df <- color_dt_df[,color_dt := "0"]
    color_dt_df <- color_dt_df[color_dt_logi == TRUE,color_dt := "1"]

    tool_df <- merge(tool_df,
                     color_dt_df[,c("gene_name","color_dt")],
                     by = "gene_name")
    plotting_df_go_term_ls[go_term_ls[i]] <- tool_df %>%
      as_tibble() %>%
      rename(y_tick = gene_name) %>%
      list()

  }

  plotting_df_go_term <- bind_rows(plotting_df_go_term_ls) %>%
    as.data.table()

  plotting_df <- bind_rows(list(plotting_df_log_pt_normal_FC,
                                plotting_df_log_pt_edge_FC,
                                plotting_df_log_pt_center_FC,
                                plotting_df_log_expr_normal,
                                plotting_df_log_expr_center,
                                plotting_df_log_expr_edge)) %>%
    as.data.table()

  plotting_df[,size_data := NA]
  plotting_df[,color_dt := NA]
  plotting_df[,graph_type := 0]
  plotting_df_go_term[,scale_dt := NA]
  plotting_df_go_term[,graph_type := 1]

  plotting_df <- bind_rows(plotting_df,
                           plotting_df_go_term) %>%
    as_tibble() %>%
    mutate(x_tick = factor(x_tick,levels = c("percent_normal",
                                             "expression_normal",
                                             "percent_edge",
                                             "expression_edge",
                                             "percent_center",
                                             "expression_center",
                                             ich_mouse@GO_enrichment[[GO_symbol]][ID %in% go_term_ls,Description])),
           y_tick = factor(y_tick,levels = gene_info_df[,gene_name]))

  plotting_df_tile <- plotting_df[plotting_df[,"graph_type"] == 0,]
  plotting_df_point <- plotting_df[plotting_df[,"graph_type"] == 1,]
  single_gene_info_heatmap <- ggplot() +
    geom_tile(data = plotting_df_tile,
              mapping = aes(x = x_tick,
                            y = y_tick,
                            fill = scale_dt)) +
    scale_fill_gradientn(colors = c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#F4A582","#D6604D","#B2182B"),
                         values = scales::rescale(c(-6,-5,-4,-3,-2,0,1,2,3),
                                                  to = c(0,1),
                                                  from = c(-6,3)),
                         name = "Percent(log2)") +
    geom_tile(data = plotting_df_point[plotting_df_point["color_dt"] == "0",],
              mapping = aes(x = x_tick,
                            y = y_tick),
              fill = "grey90",
              color = "white") +
    new_scale_fill() +
    geom_tile(data = plotting_df_point[plotting_df_point["color_dt"] == "1",],
              mapping = aes(x = x_tick,
                            y = y_tick,
                            fill = size_data),
              color = "white") +
    scale_fill_gradient2(low = "#2166AC",
                         mid = "#F7F7F7",
                         high = "#D6604D",
                         name = "Expression") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_text(family = "Arial",
                                      size = 8,
                                      color = "black",
                                      hjust = 1,
                                      angle = 45),
          axis.text.y = element_text(family = "Arial",
                                      size = 8,
                                      color = "black"))

  return(ggplotGrob(single_gene_info_heatmap))

}
