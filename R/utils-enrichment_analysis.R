
# R/utils-GO_analysis.R

#' conduct GO enrichment
#'
#' @param gene_ls a gene list
#' @param filtered_genes the filtered genes

.conduct_GO_enrichment <- function(gene_ls,filtered_genes) {

  on.exit(gc())

  gene_id <- ICHMousewch::gene_id_information[mgi_symbol %in% gene_ls]
  universe_id <- ICHMousewch::gene_id_information[mgi_symbol %in% filtered_genes]

  GO_results <- enrichGO(gene = as.character(gene_id[,entrezgene_id]),
                         universe = as.character(universe_id[,entrezgene_id]),
                         OrgDb = org.Mm.eg.db,
                         readable = TRUE,
                         ont = "ALL",
                         pvalueCutoff = 0.01) %>%
    simplify(cutoff = 0.7,
             measure = "Wang")

  GO_results <- as.data.table(GO_results@result)
  setorder(GO_results,-FoldEnrichment)

  return(GO_results)

}

#' cluster the GO term by dynamicTreeCut
#'
#' @param GO_results the results of GO enrichment

.cluster_GO_terms <- function(GO_results) {

  on.exit(gc())

  similarity_matrix <- GO_similarity(go_id = GO_results[,ID],
                                     ont = "BP",
                                     db = "org.Mm.eg.db",
                                     measure = "Sim_Resnik_1999")

  simplify_result <- simplifyGO(mat = similarity_matrix,
                                method = "dynamicTreeCut",
                                draw_word_cloud = FALSE,
                                plot = TRUE,
                                control = list(minClusterSize =3)) %>%
    as.data.table()

  cluster_symbol <- unique(simplify_result[,cluster])

  avg_similarity <- data.table(cluster = cluster_symbol,
                               avg_similarity = numeric(length(cluster_symbol)),
                               GO_num = numeric(length(cluster_symbol)))

  for (i in 1:length(cluster_symbol)) {

    GO_id <- simplify_result[cluster == cluster_symbol[i],id]
    sub_sim_mat <- similarity_matrix[GO_id,GO_id]

    if (length(GO_id) != 1) {

      avg_similarity[cluster == cluster_symbol[i],avg_similarity := mean(rowMeans(sub_sim_mat))]
      avg_similarity[cluster == cluster_symbol[i],GO_num := length(GO_id)]

    } else {

      avg_similarity[cluster == cluster_symbol[i],avg_similarity := sub_sim_mat]
      avg_similarity[cluster == cluster_symbol[i],GO_num := length(GO_id)]

    }

  }

  if (length(avg_similarity[GO_num == 1,cluster]) == 0) {

    filtered_cluster <- avg_similarity[avg_similarity > 0.5,cluster]
    unfiltered_cluster <- avg_similarity[!cluster %in% filtered_cluster,cluster]

    if (length(filtered_cluster) != 0) {

      clustered_GO_results <- vector("list",length = length(filtered_cluster))
      for (i in 1:length(filtered_cluster)) {

        clustered_GO_results[i] <- list(GO_results[ID %in% simplify_result[cluster == filtered_cluster[i],id],])

      }

    } else {

      filtered_cluster <- avg_similarity[avg_similarity == max(avg_similarity),cluster]
      unfiltered_cluster <- avg_similarity[!cluster %in% filtered_cluster,cluster]

      clustered_GO_results <- vector("list",length = length(filtered_cluster))
      for (i in 1:length(filtered_cluster)) {

        clustered_GO_results[i] <- list(GO_results[ID %in% simplify_result[cluster == filtered_cluster[i],id],])

      }

    }

    clustered_GO_results["unclustered"] <- list(GO_results[ID %in% simplify_result[cluster %in% unfiltered_cluster,id],])

  } else {

    clustered_GO_results <- list(GO_results)
    clustered_GO_results["unclustered"] <- NULL

  }

  return(clustered_GO_results)

}

#' conduct iteration cluster on GO terms
#'
#' @param GO_results the results of GO enrichment

.conduct_iteration_cluster_on_GO_terms <- function(GO_results) {

  on.exit(gc())

  GO_results <- list(unclustered = GO_results)

  condition <- TRUE
  while (condition) {

    clustered_GO_term <- ICHMousewch:::.cluster_GO_terms(GO_results = GO_results[["unclustered"]])

    GO_results <- c(GO_results,clustered_GO_term[!names(clustered_GO_term) %in% "unclustered"])
    GO_results["unclustered"] <- clustered_GO_term["unclustered"]

    if (is.null(GO_results[["unclustered"]])) {

      GO_results <- GO_results[!names(GO_results) %in% "unclustered"]
      condition <- FALSE

    }

  }

  for (i in 1:length(GO_results)) {

    names(GO_results)[i] <- paste("cluster",i,sep = ".")

  }

  return(GO_results)

}

#' integrate iteration cluster results
#'
#' @param GO_results the results of GO enrichment

.integrate_iteration_cluster_results <- function(GO_results) {

  on.exit(gc())

  similarity_matrix <- GO_similarity(go_id = GO_results[,ID],
                                     ont = "BP",
                                     db = "org.Mm.eg.db",
                                     measure = "Sim_Resnik_1999")

  iteration_cluster_results <- ICHMousewch:::.conduct_iteration_cluster_on_GO_terms(GO_results = GO_results)

  cluster_na <- names(iteration_cluster_results)

  avg_similarity_matrix <- vector("list",length = length(iteration_cluster_results))
  names(avg_similarity_matrix) <- cluster_na
  for (i in 1:length(cluster_na)) {

    avg_sim_ls <- vector("list",length = length(iteration_cluster_results))
    names(avg_sim_ls) <- cluster_na
    for (j in 1:length(cluster_na)) {

      go_id <- c(iteration_cluster_results[[cluster_na[i]]][,ID],iteration_cluster_results[[cluster_na[j]]][,ID]) %>%
        unique()

      sub_sim_mat <- similarity_matrix[go_id,go_id]

      avg_sim <- mean(rowMeans(sub_sim_mat))
      avg_sim_ls[cluster_na[j]] <- list(avg_sim)

    }

    avg_similarity_matrix[cluster_na[i]] <- list(unlist(avg_sim_ls))

  }

  avg_similarity_matrix <- as.data.frame(avg_similarity_matrix)
  rownames(avg_similarity_matrix) <- colnames(avg_similarity_matrix)
  avg_similarity_matrix <- as.matrix(avg_similarity_matrix)
  max_similarity <- max(diag(avg_similarity_matrix))

  avg_similarity_matrix <- avg_similarity_matrix - max_similarity + 1

  return(avg_similarity_matrix)

}

#' conduct KEGG enrichment
#'
#' @param gene_ls a gene list
#' @param filtered_genes the filtered genes
#'
.conduct_KEGG_enrichment <- function(gene_ls,filtered_genes) {

  on.exit(gc())

  gene_id <- ICHMousewch::gene_id_information[mgi_symbol %in% gene_ls]
  universe_id <- ICHMousewch::gene_id_information[mgi_symbol %in% filtered_genes]

  KEGG_results <- enrichKEGG(gene = as.character(gene_id[,entrezgene_id]),
                             universe = as.character(universe_id[,entrezgene_id]),
                             organism = "mmu",
                             keyType = "kegg",
                             pvalueCutoff = 0.05,
                             use_internal_data = FALSE)

  KEGG_results <- as.data.table(KEGG_results@result)
  setorder(KEGG_results,-FoldEnrichment)

  return(KEGG_results)

}

#' conduct GSEA
#'
#' @param gene_ls a gene list

.conduct_GSEA <- function(gene_ls) {

  on.exit(gc())

  gene_id <- ICHMousewch::gene_id_information[mgi_symbol %in% gene_ls]

  GSEA_results <- gseKEGG(geneList = gene_id[,entrezgene_id],
                          organism = "mmu",
                          keyType = "kegg",
                          pvalueCutoff = 0.05)

  return(GSEA_results)

}
















