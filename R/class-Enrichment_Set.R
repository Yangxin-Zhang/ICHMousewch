
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
              .Object@gene_information <- data.table()

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
#' @param GO_term_set the set of GO ID
#' @param GO_set_name the name of the GO set
#' @param ich_mouse the class of ICH_Mouse

setGeneric(name = "add_GO_term_set",
           def = function(enrichment_set,GO_term_set,GO_set_name,ich_mouse,diff_symbol = "edge-normal") {

             standardGeneric("add_GO_term_set")

           })

#' add GO term set
#'
#' @param enrichment_set the class of Enrichment_Set
#' @param GO_term_set the set of GO ID
#' @param GO_set_name the name of the GO set
#' @param ich_mouse the class of ICH_Mouse
#' @export

setMethod(f = "add_GO_term_set",
          signature = signature(enrichment_set = "Enrichment_Set",GO_term_set = "character",GO_set_name = "character"),
          definition = function(enrichment_set,GO_term_set,GO_set_name,ich_mouse,diff_symbol) {

            on.exit(gc())

            GO_results <- enrichment_set@GO_enrich[[diff_symbol]]

            sub_GO_results <- GO_results[ID %in% GO_term_set]
            sub_GO_results[,gene := vector("list",length = length(GO_term_set))]

            gene_set <- character()
            for (i in 1:length(GO_term_set)) {

              genes <- sub_GO_results[ID %in% GO_term_set[i],geneID] %>%
                strsplit(split = "/")

              genes <- genes[[1]]

              gene_set <- c(gene_set,genes)

              sub_GO_results[ID %in% GO_term_set[i],gene := list(genes)]

            }

            enrichment_set@GO_set[GO_set_name] <- list(sub_GO_results)

            gene_information <- data.table(gene_name = gene_set,
                                           nbarcodes = integer(length = length(gene_set)),
                                           avg_expr = numeric(length = length(gene_set)),
                                           avg_log2FC = numeric(length = length(gene_set)))

            sub_count_mat <- ich_mouse@raw_count_matrix[ich_mouse@filtered_genes,ich_mouse@filtered_barcodes]

            sub_count_mat_gene_set <- sub_count_mat[gene_set,]

            sub_count_mat_order <- match(gene_set,rownames(sub_count_mat_gene_set))

            gene_information[,nbarcodes := Matrix::rowSums(sub_count_mat_gene_set[sub_count_mat_order,] != 0)]

            au_seu_obj <- CreateSeuratObject(counts = sub_count_mat) %>%
              NormalizeData(normalization.method = "LogNormalize",
                            scale.factor = 1e6)

            cpm_count_mat <- au_seu_obj@assays$RNA$data[gene_set,]

            cpm_count_mat_order <- match(gene_set,rownames(cpm_count_mat))

            avg_cpm_expr <- Matrix::rowMeans(cpm_count_mat[cpm_count_mat_order,])

            gene_information[,avg_expr := avg_cpm_expr]

            diff_expr <- ich_mouse@diff_expr_genes[[diff_symbol]][gene_name %in% gene_set]

            diff_order <- match(gene_set,diff_expr[,gene_name])

            gene_information[,avg_log2FC := diff_expr[diff_order,avg_log2FC]]

            enrichment_set@gene_information <- gene_information

            return(enrichment_set)

          })

####






