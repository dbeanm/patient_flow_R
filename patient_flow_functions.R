library(lubridate)

#weight variability centrality for one node
#this score is a proportion of the maximum possible, the scale is linear
#note this score is only defined for nodes with more than one input or output edge
#output NA means no edges, NaN means only one edge, either way that's the correct result.
# w = input or output weights for the node
weight_variability_centrality <- function(w){
  weights <- w[which(w > 0)]
  n <- length(weights)
  if(n == 0){
    return(NA)
  }
  balanced_weight <- mean(weights)
  diff <- sum(abs(balanced_weight - weights))
  score <- diff/balanced_weight
  
  #score for max possible deviation
  a <- (n-1) * (balanced_weight - 1)
  b <- sum(weights) - balanced_weight - n + 1
  max_score <- (a + b)/balanced_weight
  
  prop_max <- score/max_score
  
  return(prop_max)
  
} 

#calculate weight variability centrality for all nodes
#in an adjacency matrix
# adj_matrix = output of adjacency()
# type = "in" or "out" to calculate for input or output edges respectively
weight_variability_arr <- function(adj, type){
  if(type == 'out'){
    axis = 1
  } else {
    axis = 2
  }
  wvc <- apply(adj, axis, weight_variability_centrality)
  return(wvc)
}

## weight variability centrality applied to all nodes over time
# adj_matrix = output of adjacency_over_time
# type = "in" or "out" to calculate for input or output edges respectively
# returns a matrix with one row per time period
weight_variability_over_time <- function(adj_matrix, type){
  wvc_time <- apply(adj_matrix, 3, weight_variability_arr, type)
  return(t(wvc_time))
}

## apply any function to each timepoint in the time-adjacency matrix
# returns a matrix with one row per time period
# any additional ... arugments are passed to fn when applied
#e.g. 
#weight_variability_over_time(adj, 'in')
#is exactly equivalent to:
#apply_over_time(adj, weight_variability_arr, 'in')
apply_over_time <- function(adj_matrix, fn, ...){
  res <- apply(adj_matrix, 3, fn, ...)
  return(t(res))
}

## weight variability centrality for igraph


## make time-adjacency matrix
#source node in rows, target nodes in columns
#node out degree is the number of nonzero elements in the corresponding ROW
#node in degree is the number of nonzero elements in the corresponding COLUMN
# e.g.
# e <- data.frame('Source' = c('a','a','b'), 'Target' = c('b','b','c'))
# w <- c('a','b','c')
# adjacency(w,e)
#   a b c
# a 0 2 0
# b 0 0 1
# c 0 0 0
adjacency <- function(wards, edges){
  weights <- matrix(nrow=length(wards), ncol=length(wards))
  rownames(weights) <- wards
  colnames(weights) <- wards
  
  for(i in 1:length(wards)){
    tmp <- edges[edges$Source == wards[i],]
    for(j in 1:length(wards)){
      tmp2 <- tmp[tmp$Target == wards[j],]
      weights[i,j] <- dim(tmp2)[1]
      #remove self loops
      if(i == j){
        weights[i,j] <- 0
      }
    }
  }
  return(weights)
}

adjacency_over_time <- function(transfer_data, groupby){
  wards <- unique(c(transfer_data$Source, transfer_data$Target))
  periods <- sort(unique(transfer_data[[groupby]]))
  weights_array <- array(dim=c(length(wards), length(wards), length(periods)))
  rownames(weights_array) <- wards
  colnames(weights_array) <- wards
  for(i in c(1:length(periods))){
    p <- periods[i]
    contained <- transfer_data[[groupby]] == p
    period_data <- transfer_data[contained, ]
    weights_array[,,i] <- adjacency(wards, period_data)
  }
  return(weights_array)
}


#calculate degree from an adjacency matrix
# type = "in" or "out" to calculate for input or output edges respectively
#combine with apply_over_time to get degree of each node over time
#e.g. apply_over_time(adj, degree_adj, 'in')
degree_adj <- function(adj, type){
  if(type == 'out'){
    axis = 1
  } else {
    axis = 2
  }
  deg <- apply(adj > 0, axis, sum)
  return(deg)
}

######################## Misc
# misc functions mainly to make the released flow data usable with 
# the same functing as for raw flow data

#convert from year-week to year-month
calculate_start_of_week = function(week, year) {
  date <- ymd(paste(year, 1, 1, sep="-"))
  week(date) = week
  return(date)
}

p_to_str <- function(p){
  if(p >= 0.05){
    s <- "n.s."
  } else if(p < 0.05){
    if(p >= 0.01){
      s <- "*"
    } else if (p >= 0.001) {
      s <- "**"
    } else {
      s <- "***"
    }
  }
  return(s)
}
