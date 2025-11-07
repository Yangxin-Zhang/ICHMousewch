
# R/class-Graph_Database

#' the class of Graph_Database
#'
#' @slot graph_database the graph database
#' @slot combined_graph the combined graph
#' @export

setClass(Class = "Graph_Database",
         slots = c(graph_database = "list",
                   combined_graph = "list"))

#' combine plots
#'
#' @param graph_database the class of graph database
#' @param plot_name_ls the plot for combination
#' @param combination_num the number to combination

setGeneric(name = "combine_plots",
           def = function(graph_database,plot_name_ls,combination_num = 2) {

             standardGeneric("combine_plots")

           })

#' combine plots
#'
#' @param graph_database the class of graph database
#' @param plot_name_ls the plot for combination
#' @param combination_num the number to combination

setMethod(f = "combine_plots",
          signature = signature(graph_database = "Graph_Database"),
          definition = function(graph_database,plot_name_ls,combination_num = 2) {

            on.exit()

            if (combination_num == 2) {

              plot_1 <- graph_database@graph_database[[plot_name_ls[1]]]
              plot_2 <- graph_database@graph_database[[plot_name_ls[2]]]

              plot_name_com <- paste(plot_name_ls,collapse = "-")

              com_plot <- plot_1 + plot_2

              graph_database@combined_graph[plot_name_com] <- com_plot %>%
                patchworkGrob() %>%
                list()

            }

            if (combination_num == 4) {

              plot_1 <- graph_database@graph_database[[plot_name_ls[1]]]
              plot_2 <- graph_database@graph_database[[plot_name_ls[2]]]
              plot_3 <- graph_database@graph_database[[plot_name_ls[3]]]
              plot_4 <- graph_database@graph_database[[plot_name_ls[4]]]

              plot_name_com <- paste(plot_name_ls,collapse = "-")

              com_plot <- (plot_1 | plot_2) / (plot_3 | plot_4)

              graph_database@combined_graph[plot_name_com] <- com_plot %>%
                patchworkGrob() %>%
                list()

            }

            return(graph_database)

          })

