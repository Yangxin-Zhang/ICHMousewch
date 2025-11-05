
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
          definition = function(enrichment_set,GO_term_set_ls,ich_mouse,diff_symbol = "edge-normal") {

            on.exit(gc())

            GO_results <- enrichment_set@GO_enrich[[diff_symbol]]

            GO_set_name_ls <- names(GO_term_set_ls)

            GO_names <- unlist(GO_term_set_ls)
            gene_set_ls <- vector("list",length = length(GO_names))
            names(gene_set_ls) <- GO_names
            for (i in 1:length(GO_names)) {

              GO_genes <- GO_results[ID %in% GO_names[i],geneID] %>%
                strsplit(split = "/")

              GO_genes <- GO_genes[[1]]

              GO_genes <- GO_genes[GO_genes %in% ich_mouse@filtered_genes]

              gene_set_ls[GO_names[i]] <- list(GO_genes)

            }

            sub_GO_results <- vector("list",length = length(GO_set_name_ls))
            names(sub_GO_results) <- GO_set_name_ls
            for (i in 1:length(GO_set_name_ls)) {

              go_re <- GO_results[ID %in% GO_term_set_ls[[GO_set_name_ls[i]]]]

              for (j in 1:length(GO_term_set_ls[[GO_set_name_ls[i]]])) {

                go_re[ID == GO_term_set_ls[[GO_set_name_ls[i]]][j], gene := gene_set_ls[GO_term_set_ls[[GO_set_name_ls[i]]][j]]]

              }

              sub_GO_results[GO_set_name_ls[i]] <- list(go_re)

            }

            enrichment_set@GO_set <- sub_GO_results

            gene_set_whole <- gene_set_ls %>%
              unlist() %>%
              unique()

            gene_information <- data.table(gene_name = gene_set_whole,
                                           nbarcodes = integer(length = length(gene_set_whole)),
                                           avg_expr = numeric(length = length(gene_set_whole)),
                                           avg_log2FC = numeric(length = length(gene_set_whole)))

            sub_count_mat <- ich_mouse@raw_count_matrix[ich_mouse@filtered_genes,ich_mouse@filtered_barcodes]

            sub_count_mat_gene_set <- sub_count_mat[gene_set_whole,]

            sub_count_mat_order <- match(gene_set_whole,rownames(sub_count_mat_gene_set))

            gene_information[,nbarcodes := Matrix::rowSums(sub_count_mat_gene_set[sub_count_mat_order,] != 0)]

            au_seu_obj <- CreateSeuratObject(counts = sub_count_mat) %>%
              NormalizeData(normalization.method = "LogNormalize",
                            scale.factor = 1e6)

            cpm_count_mat <- au_seu_obj@assays$RNA$data[gene_set_whole,]

            cpm_count_mat_order <- match(gene_set_whole,rownames(cpm_count_mat))

            avg_cpm_expr <- Matrix::rowMeans(cpm_count_mat[cpm_count_mat_order,])

            gene_information[,avg_expr := avg_cpm_expr]

            diff_expr <- ich_mouse@diff_expr_genes[[diff_symbol]][gene_name %in% gene_set_whole]

            diff_order <- match(gene_set_whole,diff_expr[,gene_name])

            gene_information[,avg_log2FC := diff_expr[diff_order,avg_log2FC]]

            gene_information <- rbindlist(list(enrichment_set@gene_information,gene_information))

            enrichment_set@gene_information <- gene_information

            enrichment_set <- ICHMousewch:::.calculate_GO_term_Count(enrichment_set = enrichment_set)

            return(enrichment_set)

          })

####






