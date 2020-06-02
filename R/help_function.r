#--------------------------------------------------- #
# Functions to create bootstrap community detection  #
#--------------------------------------------------- #
# rm(list = ls())
# Libraries
extract_years <- function(x) {
  graph_attr(x, "year")
}

find_range <- function(x)
{
  c(min(x), max(x))
}

identify_years_function <- function(input, interval)
{
  year_data <- map(input, function(x){graph_attr(x, interval)}) 
  return(year_data)
}

identify_weighted_edges <- function(input)
{
  edge_names <- map(input, edge_attr_names)
}

zero_one_mat <- function(x){
  return(
    (x -  min(x, na.rm = TRUE))/(
      max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  )
}

