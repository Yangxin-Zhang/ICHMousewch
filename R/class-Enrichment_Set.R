
# R/class-Enrichment_Set.R

#' a class for storing enrichment result
#'
#' @slot GO_enrich the result of GO enrich
#' @slot KEGG_enrich the result of KEGG enrich

setClass(Class = "Enrichment_Set",
         slots = c(GO_enrich = "data.table",
                   KEGG_enrich = "data.table"))

#' initialize Enrichment_Set
#'
#' @param gene_list_symbol the gene list to conduct enrichment analysis
#' @param ich_mouse the ICH_Mouse class

setMethod(f = "initialize",
          signature = signature(.Object = "Enrichment_Set"),
          definition = function(.Object,gene_list_symbol,ich_mouse) {

            gene_list <- ich_mouse@diff_expr_genes[[gene_list_symbol]]
            gene_list <- gene_list[avg_log2FC > 1,gene_name]

            .Object@GO_enrich <- ICHMousewch:::.conduct_GO_enrichment(gene_ls = gene_list,
                                                                      filtered_genes = ich_mouse@filtered_genes)

            .Object@KEGG_enrich <- ICHMousewch:::.conduct_KEGG_enrichment(gene_ls = gene_list,
                                                                          filtered_genes = ich_mouse@filtered_genes)

            validObject(.Object)
            return(.Object)

          })

#' constructor of Enrichment_Set
#'
#' @param gene_list_symbol the gene list to conduct enrichment analysis
#' @param ich_mouse the ICH_Mouse class
#' @export

Create_Enrichment_Set <- function(gene_list_symbol,ich_mouse) {

  on.exit(gc())

  enrichment_set <- new(Class = "Enrichment_Set",
                        gene_list_symbol = gene_list_symbol,
                        ich_mouse = ich_mouse)

}








