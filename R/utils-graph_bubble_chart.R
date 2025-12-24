
# R/utils-graph_bubble_chart.R

#' plotting bubble chart
#'
#' @param ich_mouse the ICH_Mouse class
#' @param group_na the group name

.plotting_bubble_chart <- function(ich_mouse,
                                   gene_num = 10,
                                   group_na = character())
  {

  on.exit(gc())

  if (length(group_na) == 0) {

    GO_term_set_ls <- ich_mouse@GO_ID_group

  } else {

    GO_term_set_ls <- ich_mouse@GO_ID_group[group_na]

  }

  num_go_term <- unlist(GO_term_set_ls) %>%
    length()

  enri_set <- ICHMousewch::Create_Enrichment_Set(ich_mouse = ich_mouse)
  enri_set <- ICHMousewch::add_GO_term_set(enrichment_set = enri_set,
                                           GO_term_set_ls = GO_term_set_ls,
                                           ich_mouse = ich_mouse)

  gene_info <- enri_set@gene_information
  GO_info <- enri_set@GO_set
  total_GO_info <- rbindlist(GO_info) %>%
    unique(by = "ID")

  plotting_dataset <- ICHMousewch:::.generate_bubble_dataset(gene_info = gene_info,
                                                             GO_info = GO_info,
                                                             gene_num = gene_num)

  posi_y_da <- plotting_dataset[y_symbol == "posi_y"] %>%
    as.data.frame() %>%
    mutate( GO_Description = factor(GO_Description,levels = unique(plotting_dataset[,GO_Description])),
            GO_group = factor(GO_group,levels = unique(plotting_dataset[,GO_group])),
            text_y = (plotting_dataset[y_symbol == "posi_y",bubble_y] + 5))

  negt_y_da <- plotting_dataset[y_symbol == "negt_y"] %>%
    as.data.frame() %>%
    mutate(GO_Description = factor(GO_Description,levels = unique(plotting_dataset[,GO_Description])),
           GO_group = factor(GO_group,levels = unique(plotting_dataset[,GO_group])))

  y_tick <- ICHMousewch:::.generate_bubble_chart_y_tick(plotting_dataset = plotting_dataset)

  bubble_color_value <- plotting_dataset[,bubble_color] %>%
    unique()
  names(bubble_color_value) <- plotting_dataset[,color_symbol] %>%
    unique()

  color_order <- names(bubble_color_value)
  color_order <- c(color_order[!color_order %in% "not included"],"not included")

  plotting_size <- plotting_dataset[y_symbol == "negt_y",plotting_size_data] %>%
    unique()

  size_breaks <- seq(from = min(plotting_size),
                     to = max(plotting_size),
                     length.out = 4)

  bubble_chart <- ggplot(data = plotting_dataset) +
    geom_point(data = negt_y_da,
               mapping = aes(x = GO_Description, y = bubble_y,color = color_symbol,size = plotting_size_data),
               alpha = 0.5) +
    scale_color_manual(values = bubble_color_value,
                       limits = color_order) +
    scale_size_identity(guide = "legend",
                        breaks = size_breaks) +
    geom_hline(yintercept = 0,
               color = "#000000FF",
               linewidth = 0.5) +
    geom_bar(data = posi_y_da,
             mapping = aes(x = GO_Description, y = bubble_y,fill = color_symbol),
             stat = "identity",
             width = 0.5) +
    scale_fill_manual(values = bubble_color_value) +
    scale_y_continuous(breaks = y_tick,
                       labels = names(y_tick)) +
    labs(title = "GO Cluster Overview",
         size = "Expression",
         color = "GO Group",
         y = "n_genes") +
    ICHMousewch:::.plotting_theme() +
    theme(axis.title.y = element_text(hjust = 0.8),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30,
                                     hjust = 1),
          aspect.ratio = (15/(num_go_term))) +
    guides(color = guide_legend(position = "left",
                                theme = theme(legend.text = element_text(size = 12,
                                                                         family = "Arial",
                                                                         vjust = 0.5,
                                                                         hjust = 0),
                                              legend.title = element_text(size = 12,
                                                                          family = "Arial",
                                                                          hjust = 0.5),
                                              legend.key.size = unit(12,"pt")),
                                override.aes = list(size = 5)),
           size = guide_legend(position = "right",
                               theme = theme(legend.text = element_blank(),
                                             legend.title = element_text(size = 12,
                                                                         family = "Arial",
                                                                         hjust = 0.5),
                                             legend.key.size = unit(12,"pt")),
                               override.aes = list(size = scales::rescale(size_breaks,to = c(2,5)))),
           fill = guide_none())

  bubble_chart_ls <- list("bubble_chart" = ggplotGrob(bubble_chart))

  return(bubble_chart_ls)

}

#' generate bubble dataset
#'
#' @param gene_info the gene information
#' @param GO_info the GO term information

.generate_bubble_dataset <- function(gene_info,
                                     GO_info,
                                     gene_num,
                                     size_param = 0.5)
  {

  on.exit(gc())

  for (i in 1:length(GO_info)) {

    setorder(GO_info[[i]],-Count)

  }

  GO_results <- rbindlist(GO_info) %>%
    unique(by = "ID")

  setorder(gene_info,-GO_term_Count)
  gene_ls <- gene_info[c(1:gene_num),gene_name]

  diff_da <- gene_info[gene_name %in% gene_ls]

  diff_order <- match(gene_ls,diff_da[,gene_name])

  diff_da <- diff_da[diff_order,size_data]

  GO_id <- GO_results[,ID]
  GO_description <- GO_results[,Description]
  GO_num <- length(GO_id)

  posi_y_da <- data.table(GO_ID = GO_id,
                          GO_Description = GO_description,
                          GO_group = character(GO_num),
                          bubble_y = integer(GO_num),
                          y_symbol = rep("posi_y",times = GO_num),
                          size_data = rep(NA,times = GO_num),
                          plotting_size_data = rep(NA,times = GO_num),
                          gene_name = rep(NA,times = GO_num),
                          bubble_color = rep(NA,times = GO_num))
  GO_order <- match(GO_id,GO_results[,ID])
  posi_y_da[,bubble_y := GO_results[GO_order,Count]]

  negt_y_da <- data.table(GO_ID = rep(GO_id,times = gene_num),
                          GO_Description = rep(GO_description,times = gene_num),
                          GO_group = character(GO_num*gene_num),
                          bubble_y = rep(seq(from = -5, to = -max(GO_results[,Count]), length.out = gene_num),each = GO_num),
                          y_symbol = rep("negt_y",times = GO_num*gene_num),
                          size_data = rep(diff_da,each = GO_num),
                          plotting_size_data = rep(diff_da*size_param,each = GO_num),
                          gene_name = rep(gene_ls,each = GO_num),
                          bubble_color = character(GO_num))

  bubble_dataset <- rbindlist(list(posi_y_da,negt_y_da))

  group_na <- names(GO_info)

  group_color <- scales::hue_pal()(length(group_na))
  names(group_color) <- group_na

  for (i in 1:length(group_na)) {

    group_GO_id <- GO_info[[group_na[i]]][,ID]

    for (j in 1:length(group_GO_id)) {

      gene_na_ls <- GO_info[[group_na[i]]][ID %in% group_GO_id[j],gene][[1]]

      bubble_dataset[GO_ID %in% group_GO_id[j],GO_group := group_na[i]]

      selected_gene_na <- bubble_dataset[!is.na(gene_name) & gene_name %in% gene_na_ls,gene_name]
      unselected_gene_na <- bubble_dataset[!gene_name %in% selected_gene_na,gene_name]

      bubble_dataset[GO_ID %in% group_GO_id[j] & !is.na(gene_name) & gene_name %in% selected_gene_na,bubble_color := group_color[group_na[i]]]
      bubble_dataset[GO_ID %in% group_GO_id[j] & !is.na(gene_name) & gene_name %in% unselected_gene_na,bubble_color := "#000000FF"]

      bubble_dataset[GO_ID %in% group_GO_id[j] & is.na(gene_name),bubble_color := group_color[group_na[i]]]
    }

  }

  bubble_dataset[,color_symbol := bubble_dataset[,GO_group]]
  bubble_dataset[bubble_color == "#000000FF",color_symbol := "not included"]

  return(bubble_dataset)

}

#' generate bubble chart y tick
#'
#' @param plotting_dataset the plotting dataset

.generate_bubble_chart_y_tick <- function(plotting_dataset)
  {

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
