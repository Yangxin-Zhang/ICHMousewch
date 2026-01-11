
# R/utils-identify_gene_distribution_pattern.R

#' identify gene distribution pattern
#'
#' @param ich_mouse the class of ICH_Mouse
#' @param threshold the avg_log2FC threshold

.identify_gene_distribution_pattern <- function(ich_mouse,
                                                threshold)
  {

  on.exit(gc())

  gene_distribution_map <- ICHMousewch:::.create_gene_distribution_map(ich_mouse = ich_mouse,
                                                                       aim_gene = ich_mouse@diff_expr_genes$`edge-normal`[avg_log2FC >= threshold,gene_name])

  gene_na_ls <- names(gene_distribution_map)

  choice_normal_edge <- character()
  choice_hematoma_edge <- character()
  choice_hematoma_center <- character()
  choice_no_significance <- character()

  for (i in 1:length(gene_na_ls)) {

    ICHMousewch:::show_image(gene_distribution_map[[gene_na_ls[i]]])

    choice_title <- paste("Choose Distribution Pattern","Gene Name:",sep = "\n") %>%
      paste(gene_na_ls[i],sep = " ") %>%
      paste("Order Number:",sep = "\n") %>%
      paste(length(gene_na_ls),sep = " ") %>%
      paste(i,sep = "-")
    choice <- menu(choices = c("Normal-Edge","Hematoma-Edge","Hematoma-Center","No-Significance"),
                   title = choice_title,
                   graphics = FALSE)

    if (choice == 1) {

      choice_normal_edge <- c(choice_normal_edge,gene_na_ls[i])

    }

    if (choice == 2) {

      choice_hematoma_edge <- c(choice_hematoma_edge,gene_na_ls[i])

    }

    if (choice == 3) {

      choice_hematoma_center <- c(choice_hematoma_center,gene_na_ls[i])

    }

    if (choice == 4) {

      choice_no_significance <- c(choice_no_significance,gene_na_ls[i])

    }

  }

  choice_results <- list("Normal-Edge" = choice_normal_edge,
                         "Hematoma-Edge" = choice_hematoma_edge,
                         "Hematoma-Center" = choice_hematoma_center,
                         "No-Significance" = choice_no_significance)

  return(choice_results)

}
