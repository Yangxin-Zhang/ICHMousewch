
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
    simplify(cutoff = 0.8,
             measure = "Wang")

  GO_results <- as.data.table(GO_results@result)

  GO_results <- GO_results[Count > 3]
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
                                plot = FALSE) %>%
    as.data.table()

  cluster_symbol <- unique(simplify_result[,cluster])

  avg_similarity <- data.table(cluster = cluster_symbol,
                               avg_similarity = numeric(length(cluster_symbol)),
                               GO_num = numeric(length(cluster_symbol)),
                               minimal_value = numeric(length(cluster_symbol)))

  for (i in 1:length(cluster_symbol)) {

    GO_id <- simplify_result[cluster == cluster_symbol[i],id]
    sub_sim_mat <- similarity_matrix[GO_id,GO_id]

    if (length(GO_id) != 1) {

      avg_similarity[cluster == cluster_symbol[i],avg_similarity := mean(rowMeans(sub_sim_mat))]
      avg_similarity[cluster == cluster_symbol[i],GO_num := length(GO_id)]
      avg_similarity[cluster == cluster_symbol[i],minimal_value := sum(sub_sim_mat < 0.2)/(length(GO_id)*length(GO_id))]

    } else {

      avg_similarity[cluster == cluster_symbol[i],avg_similarity := sub_sim_mat]
      avg_similarity[cluster == cluster_symbol[i],GO_num := length(GO_id)]
      avg_similarity[cluster == cluster_symbol[i],minimal_value := sum(sub_sim_mat < 0.2)]

    }

  }

  if (length(avg_similarity[,cluster]) != 1) {

    if (length(avg_similarity[GO_num == 1,cluster]) == 0) {

      if (nrow(avg_similarity[minimal_value < 0.2]) != 0) {

        filtered_avg_similarity <- avg_similarity[minimal_value < 0.2]
        max_sim <- max(filtered_avg_similarity[,avg_similarity])
        filtered_cluster <- filtered_avg_similarity[avg_similarity == max_sim,cluster]
        unfiltered_cluster <- avg_similarity[!cluster %in% filtered_cluster,cluster]

      } else {

        max_sim <- max(avg_similarity[,avg_similarity])
        filtered_cluster <- avg_similarity[avg_similarity == max_sim,cluster]
        unfiltered_cluster <- avg_similarity[!cluster %in% filtered_cluster,cluster]

      }

      clustered_GO_results <- vector("list",length = length(filtered_cluster))
      for (i in 1:length(filtered_cluster)) {

        clustered_GO_results[i] <- list(GO_results[ID %in% simplify_result[cluster == filtered_cluster[i],id],])

      }

      clustered_GO_results["unclustered"] <- list(GO_results[ID %in% simplify_result[cluster %in% unfiltered_cluster,id],])

    } else {

      if (max(avg_similarity[,GO_num]) <= 3) {

        clustered_GO_results <- list(GO_results)
        clustered_GO_results["unclustered"] <- list(NULL)

      } else {

        max_sim <- max(avg_similarity[GO_num > 3,avg_similarity])
        filtered_cluster <- avg_similarity[avg_similarity == max_sim,cluster]
        unfiltered_cluster <- avg_similarity[!cluster %in% filtered_cluster,cluster]

        clustered_GO_results <- vector("list",length = length(filtered_cluster))
        for (i in 1:length(filtered_cluster)) {

          clustered_GO_results[i] <- list(GO_results[ID %in% simplify_result[cluster == filtered_cluster[i],id],])

        }

        clustered_GO_results["unclustered"] <- list(GO_results[ID %in% simplify_result[cluster %in% unfiltered_cluster,id],])

      }

    }

  } else {

    clustered_GO_results <- list(GO_results)
    clustered_GO_results["unclustered"] <- list(NULL)

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

    if (is.null(GO_results[["unclustered"]])) {

      GO_results <- GO_results[!names(GO_results) %in% "unclustered"]
      condition <- FALSE

    } else {

      clustered_GO_term <- ICHMousewch:::.cluster_GO_terms(GO_results = GO_results[["unclustered"]])
      GO_results["unclustered"] <- clustered_GO_term["unclustered"]
      GO_results <- c(GO_results,clustered_GO_term[!names(clustered_GO_term) %in% "unclustered"])

    }

  }

  for (i in 1:length(GO_results)) {

    names(GO_results)[i] <- paste("cluster",i,sep = ".")

  }

  return(GO_results)

}

#' adjust iteration cluster results
#'
#' @param GO_results the results of GO enrichment

.adjust_iteration_cluster_results <- function(GO_results) {

  on.exit(gc())

  iteration_cluster_results <- ICHMousewch:::.conduct_iteration_cluster_on_GO_terms(GO_results = GO_results)

  similarity_matrix <- GO_similarity(go_id = rbindlist(iteration_cluster_results)[,ID],
                                     ont = "BP",
                                     db = "org.Mm.eg.db",
                                     measure = "Sim_Resnik_1999")

  cluster_na <- names(iteration_cluster_results)
  GO_id_ls <- vector("list",length = length(cluster_na))
  names(GO_id_ls) <- cluster_na
  for (i in 1:length(cluster_na)) {

    GO_id_ls[cluster_na[i]] <- list(iteration_cluster_results[[cluster_na[i]]][,ID])

  }

  pair_ls <- list()
  for (i in 1:length(cluster_na)) {

    pair_mat <- combn(GO_id_ls[[cluster_na[i]]],2)

    for (j in 1:ncol(pair_mat)) {

      pair_ls <- append(pair_ls,list(pair_mat[,j]))
      names(pair_ls)[length(pair_ls)] <- paste(cluster_na[i],j,sep = "_")

    }

  }

  pair_na <- names(pair_ls)
  pair_overlap <- data.table(pair_name = pair_na,
                             overlap_num = numeric(length(pair_na)),
                             GO_cluster = character(length(pair_na)))
  for (i in 1:length(pair_na)) {

    clu_na <- strsplit(pair_na[i],split = "_")[[1]][1]
    GO_id_1 <- pair_ls[[pair_na[i]]][1]
    GO_id_2 <- pair_ls[[pair_na[i]]][2]
    pair_sim_1 <- similarity_matrix[!rownames(similarity_matrix) %in% GO_id_ls[[clu_na]],GO_id_1]
    pair_sim_2 <- similarity_matrix[!rownames(similarity_matrix) %in% GO_id_ls[[clu_na]],GO_id_2]

    overlap_na_1 <- names(pair_sim_1)[pair_sim_1 > 0.2]
    overlap_na_2 <- names(pair_sim_2)[pair_sim_2 > 0.2]

    pair_overlap[pair_name == pair_na[i],overlap_num := sum(overlap_na_1 %in% overlap_na_2)]
    pair_overlap[pair_name == pair_na[i],GO_cluster := clu_na]

  }

  for (i in 1:length(cluster_na)) {

    sub_pair_overlap <- pair_overlap[GO_cluster %in% cluster_na[i]]

    setorder(sub_pair_overlap,-overlap_num)

    ordered_GO_id <- pair_ls[sub_pair_overlap[,pair_name]] %>%
      unlist() %>%
      unique()

    GO_id_order <- match(ordered_GO_id,iteration_cluster_results[[cluster_na[i]]][,ID])

    iteration_cluster_results[cluster_na[i]] <- list(iteration_cluster_results[[cluster_na[i]]][GO_id_order])

  }

  pair_cluster <- combn(cluster_na,2)
  pair_cluster_ls <- vector("list", length = ncol(pair_cluster))
  for (i in 1:length(pair_cluster_ls)) {

    pair_cluster_ls[i] <- list(pair_cluster[,i])
    names(pair_cluster_ls)[i] <- paste(pair_cluster_ls[[i]][1],pair_cluster_ls[[i]][2],sep = "_")

  }

  pair_cluster_na <- names(pair_cluster_ls)

  avg_similarity_dt <- data.table(pair_cluster = pair_cluster_na,
                                  avg_similarity_diff = numeric(length(pair_cluster_na)))
  for (i in 1:length(pair_cluster_na)) {

    sub_GO_results <- rbindlist(iteration_cluster_results[pair_cluster_ls[[pair_cluster_na[i]]]])
    sub_sim_mat <- similarity_matrix[sub_GO_results[,ID],sub_GO_results[,ID]]

    sub_GO_results_1 <- rbindlist(iteration_cluster_results[pair_cluster_ls[[pair_cluster_na[i]]][1]])
    sub_GO_results_2 <- rbindlist(iteration_cluster_results[pair_cluster_ls[[pair_cluster_na[i]]][2]])
    sub_sim_mat_1 <- similarity_matrix[sub_GO_results_1[,ID],sub_GO_results_1[,ID]]
    sub_sim_mat_2 <- similarity_matrix[sub_GO_results_2[,ID],sub_GO_results_2[,ID]]

    avg_sim <- sub_sim_mat %>%
      rowMeans() %>%
      mean()

    if (ncol(sub_sim_mat_1) != 1) {

      avg_sim_1 <- sub_sim_mat_1 %>%
        rowMeans() %>%
        mean()

    } else {

      avg_sim_1 <- sub_sim_mat_1[1,1]

    }

    if (ncol(sub_sim_mat_2) != 1) {

      avg_sim_2 <- sub_sim_mat_2 %>%
        rowMeans() %>%
        mean()

    } else {

      avg_sim_2 <- sub_sim_mat_2[1,1]

    }

    avg_similarity_dt[pair_cluster == pair_cluster_na[i],avg_similarity_diff := (avg_sim - mean(avg_sim_1,avg_sim_2))/mean(avg_sim_1,avg_sim_2)]

  }

  setorder(avg_similarity_dt,-avg_similarity_diff)

  pair_cluster_vec <- pair_cluster_ls[avg_similarity_dt[,pair_cluster]] %>%
    unlist() %>%
    unique()

  pair_cluster_order <- match(pair_cluster_vec,names(iteration_cluster_results))

  iteration_cluster_results <- iteration_cluster_results[pair_cluster_order]

  names(iteration_cluster_results) <- cluster_na

  return(iteration_cluster_results)

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

#' combine GO cluster based on cluster annotation
#'
#' @param GO_results GO enrichment results
#' @param cluster_annotation the cluster annotation

.combine_GO_cluster_based_on_cluster_annotation <- function(GO_results,cluster_annotation) {

  on.exit(gc())

  cluster_na <- names(GO_results)

  annotation_na <- unique(cluster_annotation)

  inte_GO_results <- vector("list", length = length(annotation_na))
  names(inte_GO_results) <- annotation_na
  for (i in 1:length(annotation_na)) {

    inte_cluster <- names(cluster_annotation)[cluster_annotation %in% annotation_na[i]]

    inte_GO_results[annotation_na[i]] <- list(rbindlist(GO_results[inte_cluster]))

  }

  return(inte_GO_results)

}

#' combine a GO enrichment heatmap
#'
#' @param GO_results_ls the GO results list

.combine_GO_enrichment_heatmap <- function(GO_results_ls) {

  on.exit(gc())

  GO_results <- rbindlist(GO_results_ls)

  sim_mat <- GO_similarity(go_id = GO_results[,ID],
                           ont = "BP",
                           db = "org.Mm.eg.db",
                           measure = "Sim_Resnik_1999")

  return()

}












