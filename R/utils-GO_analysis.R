
# R/utils-GO_analysis.R

#' conduct GO enrichment
#'
#' @param gene_ls a gene list

.conduct_GO_enrichment <- function(gene_ls) {

  on.exit(gc())

  gene_id <- ICHMousewch::gene_id_information[mgi_symbol %in% gene_ls]

  GO_results <- enrichGO(gene = gene_id[],
                         OrgDb = org.Mm.eg.db,
                         readable = TRUE,
                         ont = "ALL",
                         pvalueCutoff = 0.05)

  return(GO_results)

}

#' conduct KEEG enrichment
#'
#' @param gene_ls a gene list

.conduct_KEEG_enrichment <- function(gene_ls) {

  on.exit(gc())

  gene_id <- ICHMousewch::gene_id_information[mgi_symbol %in% gene_ls]

  KEGG_results <- enrichKEGG(gene = gene_id[,],
                             organism = "mmu",
                             keyType = "kegg",
                             pvalueCutoff = 0.05)

  return(KEGG_results)

}

#' conduct GSEA
#'
#' @param gene_ls a gene list

.conduct_GSEA <- function(gene_ls) {

  on.exit(gc())

  gene_id <- ICHMousewch::gene_id_information[mgi_symbol %in% gene_ls]

  GSEA_results <- gseKEGG(geneList = gene_id[,],
                          organism = "mmu",
                          keyType = "kegg",
                          pvalueCutoff = 0.05)

  return(GSEA_results)

}
















