###############################################################
# This is the main source for all function needed to perform
# Clustering on cells based on density estimations

###############################################################
library(dbscan)


#### Function for estimating density

calc_rho <- function(df, k = 7, estimator){
  ## A function to estimate densities with 3 methods
  
  # :Param df: An (N by d) matrix with coordinates of datapoints 
  # :Param k: Number of neighbors when calculating the Adjacency Matrix.
  # :Param estimator: 1 for 1/distance, 2 for stationary distribution
  
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
  for (i in 1:nrow(df)){
    data <- cbind(df, rho, 1:nrow(df))[order(rho, decreasing = TRUE),]
    
    if (i == 1) {
      point <- data[1, 1:2]
      distance <- apply(data[,1:2], 1, function(x) sqrt(sum((x - point)^2)))
      delta[1, ] <- c(max(distance), data[1, 4])
    }
    
    else if (i == 2) {
      delta[2, ] <- c(sqrt(sum((data[1, 1:2] - data[2, 1:2])^2)), data[2, 4])
    }
    
    else {
      point <- data[i, 1:2]
      distance <- apply(data[1:i-1, 1:2], 1, function(x) sqrt(sum((x - point)^2)))
      delta[i, ] <- c(min(distance), data[i, 4])
      
      if (track == TRUE & i %% 200 == 0){
        print(i)
      }
    }
  }
  return(delta[order(delta[,2])])
}



#### Data for Clusters vs Constant plot
clust_vs_const <- function(df, k = 7, estimator = 2, track = FALSE){
  
  # By default k = 7 is used for the k-nearest neighbor distances
  # and the stationary distribution is used for density estimation
  
  rho = calc_rho(df, k, estimator)
  delta = calc_delta(df, rho, track)
  #Find slope of ln(rho) + a*ln(delta)
  non_extreme <- which(delta > quantile(delta, 0.05) & delta < quantile(delta, 0.95))
  delta_trimmed <- delta[non_extreme]
  rho_trimmed <- rho[non_extreme]
  a <- lm(log(delta_trimmed)~log(rho_trimmed))$coefficient[[2]]
  #Find which constant to use for the decision line ln(rho) + a*ln(delta) = const
  const <-  seq(0, 1, by = 0.001)
  clusters <- c()
  
  for (i in 1:length(const)){
    clusters[i] <- sum(delta > (const[i]*rho)^a)
  }
  return(list(rho = rho, delta = delta, coefs = cbind(const, clusters), slope = a))
}


#### Returns the ID of potential cluster centers
cluster_centers <- function(delta, rho, treshold){
  #Find slope of ln(rho) + a*ln(delta)
  non_extreme <- which(delta > quantile(delta, 0.05) & delta < quantile(delta, 0.95))
  delta_trimmed <- delta[non_extreme]
  rho_trimmed <- rho[non_extreme]
  a <- lm(log(delta_trimmed)~log(rho_trimmed))$coefficient[[2]]
  center_id = which(delta > (treshold*rho)^a)
  return(centers = center_id)
}

