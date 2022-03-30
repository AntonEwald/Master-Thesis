#' This is the main source for all function needed to perform
#' Clustering on cells based on density estimations


#' Function for estimating pointdensity
#' 
#' @param df: An (N by d) matrix with coordinates of datapoints 
#' @param k: Number of neighbors when calculating the Adjacency Matrix.
#' @param estimator: 1 for 1/distance, 2 for stationary distribution


calc_rho <- function(df, k = 7, estimator){
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


#'  Automatically decides on cluster peaks from decision graph
#'  This method takes in an Nx3 matrix with point ID, delta and rho estimations.
#'  It then fits a line through log(delta) as a function of log(rho) by
#'  monotonic regression. This regression line is then shifted from the bottom
#'  of the decision graph and up by alterning a threshold constructing a table
#'  with columns threshold and nr of clusters chosen. The number of clusters chosen
#'  is the most frequently occurence. 
#'  Possible improvment: Make decision based on max/min values instead of count
#'  
#'  @param sub_stats: An Nx3 matrix containing datapoint ID, Delta, Rho

auto_peak_finder <- function(sub_stats, avgPeakDist){
  L = 3*avgPeakDist
  delta = sub_stats[,2] + 0.0001 #R Crash if you input log(0) in y vector into isoreg()
  rho = sub_stats[,3]
  indicator = c(rep(1, floor(length(delta)*0.1)), rep(0, length(delta) - floor(length(delta)*0.1)))
  iso <- as.stepfun(isoreg(x = log(rho), y = -log(delta))) #monotonic regression line (minus on y cause we have decreasing correlation)
  a = -iso(log(rho)) # Minus again to revert back the minus above
  threshold <-  seq(exp(- min(a)), exp(log(L/4) - a[which.max(delta*indicator)]), length.out = 30) 
  clusters <- c()
  # Decides number of clusters
  for (l in 1:length(threshold)){
    clusters[l] <- sum((delta) > (exp(a)*threshold[l]))
  }
  if (length(unique(clusters)) == 1){
    nr_clusters = clusters[1]
  }
  else {
    nr_clusters <- sort(table(clusters), decreasing = TRUE) %>% 
      as.data.frame() %>% 
      mutate(clusters = as.numeric(as.character(clusters))) %>% 
      filter(Freq == max(Freq)) %>% .[,1] %>% max() #Cluster with most occurence
  }
  # Add derivative automation instead
  
  coefs <- cbind(threshold, clusters) 
  suggested_threshold <- median(coefs[coefs[,2] == nr_clusters, 1]) #Take the threshold in the middle of the peak and valley
  peak_id <- sub_stats[which(delta > (exp(a)*suggested_threshold)), 1] #Datapoint ID for cluster centers
  return(peak_id)
}



#' @param clusterStats: Object from findclusterPeaks function
clusterAssignments <- function(clusterStats){
  #Set cluster ID to cluster peaks
  cluster_id <- rep(0, length(clusterStats$rho))
  cluster_id[clusterStats$ClusterPeaks] = 1:length(clusterStats$ClusterPeaks)
  
  stats <- cbind(1:nrow(clusterStats$delta), clusterStats$delta, clusterStats$rho, cluster_id) %>% 
    .[order(.[,4], decreasing = TRUE), ] 
  
  #Give each point a cluster ID based on method from paper "Clustering by Fast Search"
  for (i in 1:nrow(stats)){
    if (stats[i, 5] == 0){ #Only look at non cluster centers
      if (is.na(stats[i, 3])){ #Points with nearest neighbor with higher density outside of our border
        stats[i, 5] = length(clusterStats$ClusterPeaks) + 1
      }
      else if (stats[i,3] != 0){ #If there is a neighbor with higher density
        stats[i, 5] = stats[stats[,1] == stats[i, 3], 5] #Give same cluster ID as nearest neighbor with higher density
      }
      else {
        stats[i, 5] = 0
      }
    }
    if (i %% 100 == 0){
      print(i)
    }
  }
  
  assignments <- as.factor(stats[order(stats[,1]),][,5])
  
}



#' Function for estimating Delta by grids instead of over entire data. 
#' Since the computational time is of order o(N^2), when the dataset is too large 
#' the computation time is not feasible. Hence we introduce a method for estimating
#' delta by grids, decreasing the computaional time significantly. We use a 9x9 grid
#' with total size as L \times L where L is 3 times the average distance between clusters.
#' We then calculate the delta for everypoint in the center part of the 9x9 grid by
#' only using data within the 9x9 grid. By doing so we do not calculate delta for the
#' points at the border grids. These will get delta = NA
#' 
#' @param df: An (N by d) matrix with coordinates of datapoints 
#' @param avgPeakDist The average distance between the center of clusters.
#' @param k: Number of neighbors when calculating the Adjacency Matrix.
#' @param estimator: 1 for 1/distance, 2 for stationary distribution

findClusterPeaks <- function(df, avgPeakDist, k = 7, estimator = 2, clustering = FALSE){
  L = avgPeakDist*3
  rho = calc_rho(df, k=k, estimator = estimator)
  #Order by rho
  sorted_rho <- cbind(df, rho, 1:length(rho)) %>% 
    .[order(.[,3], decreasing = TRUE), ]
  sorted <- cbind(1:length(rho), sorted_rho)
  xmin <- min(df[,1])
  xmax <- max(df[,1])
  ymin <- min(df[,2])
  ymax <- max(df[,2])
  xlength <- xmax - xmin
  ylength <- ymax - ymin
  grid_size <- avgPeakDist
  xgrids <- ceiling(xlength/grid_size)
  ygrids <- ceiling(ylength/grid_size)
  deltas <- matrix(NA, ncol = 2, nrow = nrow(df))
  cluster_centers <- c()
  for (i in 2:(ygrids-1)){
    for (j in 2:(xgrids-1)){
      
      large_grid <- sorted[sorted[, 2] < xmin + (j+1)*grid_size &
                             sorted[, 2] >= xmin + (j-2)*grid_size & 
                             sorted[, 3] < ymin + (i+1)*grid_size &
                             sorted[, 3] >= ymin + (i-2)*grid_size, ]
      
      small_grid <- sorted[sorted[, 2] < xmin + j*grid_size &
                             sorted[, 2] >= xmin + (j-1)*grid_size & 
                             sorted[, 3] < ymin + i*grid_size &
                             sorted[, 3] >= ymin + (i-1)*grid_size, ]
      
      if(length(small_grid) > 2*5){ #more than 6 points
        #Atleast 6 points in the small grid
        sub_stats <- matrix(NA, ncol = 3, nrow = nrow(small_grid))
        #Calculate our deltas for the small grid
        for (k in 1:nrow(small_grid)){ #For every point in the small grid
          
          if (sum(large_grid[,1] < small_grid[k, 1]) > 1){
            dists <- apply(large_grid[large_grid[,1] < small_grid[k,1], ][,2:3], 1, function(x) sqrt(sum((x-small_grid[k,2:3])^2)))
            delta = min(dists)
            neighbor_id = large_grid[which.min(dists), 5]
            deltas[small_grid[k, 5], 1] <- delta
            deltas[small_grid[k, 5], 2] <- neighbor_id
            sub_stats[k, 1] <- small_grid[k, 5]
            sub_stats[k, 2] <- delta
            sub_stats[k, 3] <- small_grid[k, 4]
          }
          else if (sum(large_grid[,1] < small_grid[k, 1]) == 0){ #There is no point with larger density
            deltas[small_grid[k, 5], 1] = L #Give max delta
            deltas[small_grid[k, 5], 2] = 0
            sub_stats[k, 1] <- small_grid[k, 5] 
            sub_stats[k, 2] <- L 
            sub_stats[k, 3] <- small_grid[k, 4]
          }
          
          else { # Only one point with higher density
            delta = sqrt((large_grid[large_grid[,1] < small_grid[k,1], ][2] - small_grid[k,2])^2 + (large_grid[large_grid[,1] < small_grid[k,1], ][3] - small_grid[k,3])^2)
            deltas[small_grid[k,5], 1] <- delta
            deltas[small_grid[k,5], 2] <- large_grid[1, 5]
            sub_stats[k, 1] <-small_grid[k, 5]
            sub_stats[k, 2] <- delta
            sub_stats[k, 3] <- small_grid[k, 4]
            
          }
        }
        #Here we should check for density peaks
        cluster_centers <- c(cluster_centers, auto_peak_finder(sub_stats, avgPeakDist))
      }
    }
    print(i)
    
  }
  
  cluster_info <- list(ClusterPeaks = cluster_centers, rho = rho, delta = deltas)
  
  if (clustering == TRUE){
    assignments = clusterAssignments(cluster_info)
  }
  return(append(cluster_info, list(assignments = assignments)))
}