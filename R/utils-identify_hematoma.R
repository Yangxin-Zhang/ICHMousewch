
#  R/utils-identify_hematoma.R

#' choose hematoma symbol
#'
#' @param hematoma the class of Hematoma

.choose_hematoma_symbol <- function(hematoma)
  {

  on.exit(gc())

  seu_meta <- hematoma@seu_metadata_with_cluster_symbol
  seu_meta[,aim_cluster := "0"]
  seu_meta[Louvain_cluster_posi == 1,aim_cluster := "Normal"]

  posi_cluster_num <- seu_meta[Louvain_cluster_posi != 1,Louvain_cluster_posi] %>%
    unique() %>%
    length()

  for (i in 1:posi_cluster_num) {

    seu_meta[Louvain_cluster_posi == i+1,aim_cluster := paste("region",(i+1))]

  }

  posi_cluster <- seu_meta[aim_cluster != "Normal",aim_cluster] %>%
    unique()

  posi_cluster_image_ls <- vector("list",length = posi_cluster_num)
  names(posi_cluster_image_ls) <- posi_cluster
  for (i in 1:posi_cluster_num) {

    aim_posi_cluster <- posi_cluster[i]
    aim_posi_cluster_color <- "#D1352B"
    names(aim_posi_cluster_color) <- aim_posi_cluster

    posi_cluster_color_set <- rep("#3C77AF",times = (posi_cluster_num-1))
    names(posi_cluster_color_set) <- posi_cluster[!posi_cluster %in% aim_posi_cluster]

    posi_cluster_image_ls[posi_cluster[i]] <- ICHMousewch:::.create_spatial_image_with_cluster_symbol(ich_mouse = hematoma,
                                                                                                      cluster_symbol = "Louvain_cluster_posi",
                                                                                                      self_definition_color = c("1"="#F5D2A8",
                                                                                                                                aim_posi_cluster_color,
                                                                                                                                posi_cluster_color_set),
                                                                                                      show_legend_label = TRUE,
                                                                                                      show_image = FALSE) %>%
      list()

  }

  hematoma_symbol <- numeric()
  for (i in 1:posi_cluster_num) {

    ICHMousewch::show_image(posi_cluster_image_ls[posi_cluster[i]])

    choice <- menu(choices = c("Yes","No"),
                   title = posi_cluster[i],
                   graphics = FALSE)

    if (choice == 1) {

      hematoma_symbol <- c(hematoma_symbol,i)

    }
  }

  return(hematoma_symbol)

}
