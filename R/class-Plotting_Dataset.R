
# R/class-Plotting_Dataset.R

#' the class for Plotting
#'
#' @slot seurat_metadata_dataset the Seurat Object metadata dataset for plotting
#' @slot differential_expression_genes_dataset the differential expression genes dataset
#' @slot GO_enrichment_dataset the dataset of GO enrichment
#' @slot plots the plots created from dataset

setClass(Class = "Plotting_Dataset",
         slots = c(seurat_metadata_dataset = "data.table",
                   differential_expression_genes_dataset = "data.table",
                   GO_enrichment_dataset = "data.table",
                   plots = "list"))

#' initialize Plotting_Dataset
#'
#' @param ich_mouse the class of ICH_Mouse
#' @param diff_expr_gene_symbol the differential expression genes symbol

setMethod(f = "initialize",
          signature = signature(.Object = "Plotting_Dataset"),
          definition = function(.Object,ich_mouse,diff_expr_gene_symbol) {

            on.exit(gc())

            .Object@seurat_metadata_dataset <- ich_mouse@seu_metadata_with_cluster_symbol

            .Object@differential_expression_genes_dataset <- ich_mouse@diff_expr_genes[[diff_expr_gene_symbol]]

            .Object@GO_enrichment_dataset <- ICHMousewch:::.conduct_GO_enrichment(gene_ls = ich_mouse@diff_expr_genes[["edge-normal"]][avg_log2FC > 1,gene_name],
                                                                                  filtered_genes = ich_mouse@filtered_genes)


            validObject(.Object)
            return(.Object)

          })

#' constructor of Plotting_Dataset
#'
#' @param ich_mouse the class of ICH_Mouse
#' @param diff_expr_gene_symbol the differential expression genes symbol

Create_Plotting_Dataset <- function(ich_mouse,diff_expr_gene_symbol = "edge-normal") {

  on.exit(gc())

  plotting_dataset <- new(Class = "Plotting_Dataset",
                          ich_mouse = ich_mouse,
                          diff_expr_gene_symbol = diff_expr_gene_symbol)

  return(plotting_dataset)

}
