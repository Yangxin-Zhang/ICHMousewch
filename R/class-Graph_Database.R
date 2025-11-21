
# R/class-Graph_Database

#' the class of Graph_Database
#'
#' @slot graph_database the graph database
#' @slot combined_graph the combined graph
#' @export

setClass(Class = "Graph_Database",
         slots = c(graph_database = "list",
                   combined_graph = "list"))

#' create empty Graph_Database
#'
#'@export

Graph_Database <- function(){

  graph_database <- new(Class = "Graph_Database")

  return(graph_database)

}

#' combine plots
#'
#' @param graph_database the class of graph database
#' @param plot_name_ls the plot for combination
#' @param combination_num the number to combination
#' @param graph_name the graph name

setGeneric(name = "combine_plots",
           def = function(plot_name_ls,graph_database = ICHMousewch::Graph_Database(),combination_num = 2,graph_name = character()) {

             standardGeneric("combine_plots")

           })

#' combine plots
#'
#' @param graph_database the class of graph database
#' @param plot_name_ls the plot for combination
#' @param combination_num the number to combination
#' @param graph_name the graph name

setMethod(f = "combine_plots",
          signature = signature(graph_database = "Graph_Database"),
          definition = function(plot_name_ls,graph_database = ICHMousewch::Graph_Database(),combination_num = 2,graph_name = character()) {

            on.exit()

            if (sum(!plot_name_ls %in% names(graph_database@graph_database)) != 0) {

              graph_database <- ICHMousewch::load_graph(graph_name = plot_name_ls,
                                                        graph_database = graph_database)

            }

            if (combination_num == 2) {

              plot_1 <- graph_database@graph_database[[plot_name_ls[1]]]
              plot_2 <- graph_database@graph_database[[plot_name_ls[2]]]

              plot_name_com <- paste(plot_name_ls,collapse = "-")

              com_plot <- plot_1 + plot_2

            }

            if (combination_num == 3) {

              plot_1 <- graph_database@graph_database[[plot_name_ls[1]]]
              plot_2 <- graph_database@graph_database[[plot_name_ls[2]]]
              plot_3 <- graph_database@graph_database[[plot_name_ls[3]]]

              plot_name_com <- paste(plot_name_ls,collapse = "-")

              com_plot <- (plot_1 | plot_2) / (plot_3)

            }

            if (combination_num == 4) {

              plot_1 <- graph_database@graph_database[[plot_name_ls[1]]]
              plot_2 <- graph_database@graph_database[[plot_name_ls[2]]]
              plot_3 <- graph_database@graph_database[[plot_name_ls[3]]]
              plot_4 <- graph_database@graph_database[[plot_name_ls[4]]]

              plot_name_com <- paste(plot_name_ls,collapse = "-")

              com_plot <- (plot_1 | plot_2) / (plot_3 | plot_4)

            }

            if (length(graph_name) != 0) {

              com_plot <- com_plot +
                plot_annotation(title = graph_name,
                                theme = theme(plot.title = element_text(size = 16,
                                                                        face = "bold",
                                                                        family = "Arial",
                                                                        hjust = 0.5
                                                                        ,vjust = 0.5,
                                                                        margin = margin(b = 1,t = 30))))

            }

            new_com_plot <- com_plot %>%
              patchworkGrob() %>%
              list()
            names(new_com_plot) <- plot_name_com
            graph_database@combined_graph <- append(graph_database@combined_graph,new_com_plot)

            return(graph_database)

          })

