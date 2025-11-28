
# R/class-Enrichment_Set.R

#' a class for storing enrichment result
#'
#' @slot GO_set the GO term set for analysis
#' @slot GO_enrich the result of GO enrich
#' @slot KEGG_enrich the result of KEGG enrich
#' @slot gene_information the gene information

setClass(Class = "Enrichment_Set",
         slots = c(GO_set = "list",
                   KEGG_enrich = "list",
                   GO_enrich = "list",
                   gene_information = "data.table"))

#' initialize Enrichment_Set
#'
#' @param ich_mouse the ICH_Mouse class

setMethod(f = "initialize",
          signature = signature(.Object = "Enrichment_Set"),
          definition = function(.Object,ich_mouse,initialization) {

            if (initialization) {

              .Object@GO_enrich <- ich_mouse@GO_enrichment
              .Object@gene_information <- data.table(gene_name = character(0),
                                                     nbarcodes = integer(0),
                                                     avg_expr = numeric(0),
                                                     avg_log2FC = numeric(0))

            } else {

              .Object@GO_set <- list()
              .Object@GO_enrich <- list()
              .Object@KEGG_enrich <- list()
              .Object@gene_information <- data.table()

            }

            validObject(.Object)
            return(.Object)

          })

#' constructor of Enrichment_Set
#'
#' @param ich_mouse the ICH_Mouse class
#' @export

Create_Enrichment_Set <- function(ich_mouse,initialization = TRUE) {

  on.exit(gc())

  enrichment_set <- new(Class = "Enrichment_Set",
                        initialization = initialization,
                        ich_mouse = ich_mouse)

}

####
#' add GO term set
#'
#' @param enrichment_set the class of Enrichment_Set
#' @param GO_term_set_ls the set of GO ID
#' @param ich_mouse the class of ICH_Mouse

setGeneric(name = "add_GO_term_set",
           def = function(enrichment_set,GO_term_set_ls,ich_mouse,diff_symbol = "edge-normal") {

             standardGeneric("add_GO_term_set")

           })

#' add GO term set
#'
#' @param enrichment_set the class of Enrichment_Set
#' @param GO_term_set_ls the set of GO ID
#' @param ich_mouse the class of ICH_Mouse
#' @export

setMethod(f = "add_GO_term_set",
          signature = signature(enrichment_set = "Enrichment_Set",GO_term_set_ls = "list"),
          definition = function(enrichment_set,GO_term_set_ls,ich_mouse) {

            on.exit(gc())

            enrichment_set <- enrichment_set %>%
              ICHMousewch:::.generate_single_gene_GO_infomation(ich_mouse = ich_mouse,
                                                                GO_term_set_ls = GO_term_set_ls) %>%
              ICHMousewch:::.calculate_GO_term_Count()

            return(enrichment_set)

          })

####






