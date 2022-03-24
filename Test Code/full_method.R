###############################################################
# This is the main source for all function needed to perform
# Clustering on cells based on density estimations

###############################################################



#### Function for estimating density

calc_rho <- function(df, k = 7, estimator){
  ## A function to estimate densities with 3 methods
  
  # :Param df: An (N by d) matrix with coordinates of datapoints 
  # :Param k: Number of neighbors when calculating the Adjacency Matrix.
  # :Param estimator: 1 for 1/distance, 2 for stationary distribution
  
  suppressWarnings(suppressMessages(library(dbscan)))
  #Weighted Adjaceny matrix
  W <- kNNdist(df, k, all = TRUE)
  #Density Estimations
  if (estimator == 1){
    rho = 1/rowSums(W)
  }
  
  else {
    rho = 1/(1 + rowSums(W^2))
  }
  return(rho)
}



#### Calculating Delta
calc_delta <- function(df, rho, track){
  
  delta <- matrix(NA, ncol = 2, nrow = nrow(df))
  data <- cbind(df, rho, 1:nrow(df))[order(rho, decreasing = TRUE),]
  point1 <- data[1, 1:2]
  distance1 <- apply(data[,1:2], 1, function(x) sqrt(sum((x - point1)^2)))
  delta[1, ] <- c(max(distance1), data[1, 4]) #First two points does not work in the Loop but are easy to do manually
  delta[2, ] <- c(sqrt(sum((data[1, 1:2] - data[2, 1:2])^2)), data[2, 4])
  for (i in 3:nrow(df)){
      point <- data[i, 1:2]
      distance <- apply(data[1:i-1, 1:2], 1, function(x) sqrt(sum((x - point)^2)))
      delta[i, ] <- c(min(distance), data[i, 4])
      
      if (track == TRUE & i %% 500 == 0){
        print(i)
      }

  }
  return(delta[order(delta[,2])])
}



#### Data for Clusters vs Constant plot
findStats <- function(df, k = 7, estimator = 2, track = FALSE){
  
  # By default k = 7 is used for the k-nearest neighbor distances
  # and the stationary distribution is used for density estimation
  
  rho = calc_rho(df, k, estimator)
  delta = calc_delta(df, rho, track)
  #Find slope of ln(rho) + a*ln(delta)
  
  ## This is for outliers but not sure if we should use that yet
  ##non_outliers <- which(rho/delta > quantile(rho/delta, 0.05))
  ##delta_trimmed <- delta[non_outliers]
  ##rho_trimmed <- rho[non_outliers]
  
  a <- lm(log(delta)~log(rho))$coefficient[[2]]
  #Find which constant to use for the decision line ln(rho) + a*ln(delta) = const
  threshold <-  seq(0, 1, by = 0.001)
  clusters <- c()
  
  for (i in 1:length(threshold)){
    clusters[i] <- sum(delta > (threshold[i]*rho)^a)
  }
  stats <- list(df = df, rho = rho, delta = delta, coefs = cbind(threshold, clusters), slope = a)
  class(stats) <- 'fuzzyDensClust'
  stats
}

# Decision Graph
plot.fuzzyDensClust <- function(x, ...){
  plot(x$coefs[,1], x$coefs[,2], type = "l", main = "Decision Graph", xlab = "Threshold", ylab = "Clusters")
}


#### Returns the ID of potential cluster centers
clustering <- function(stats, threshold){
  suppressWarnings(suppressMessages(library(Rfast)))
  delta = stats$delta
  rho = stats$rho
  a = stats$slope
  center_id = which(delta > (threshold*rho)^a)
  cent_to_point_dist <- 1/(1+dista(stats$df[center_id,], stats$df))
  assignments <- sweep(cent_to_point_dist, 2, colSums(cent_to_point_dist), "/")
  
  #Set all cluster centers to 100% chance of being to its own cluster
  
  for (i in 1:length(center_id)) {
    j <- which.max(assignments[,center_id[i]])
    assignments[,center_id[i]] = c(rep(0, i-1), 1, rep(0, length(center_id)-i))
  }
  
  return(list(centerID = center_id, assignments = assignments))
}

