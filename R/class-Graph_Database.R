
# R/class-Graph_Database

#' the class of Graph_Database
#'
#' @slot graph_database the graph database
#' @export

setClass(Class = "Graph_Database",
         slots = c(graph_database = "list"))

