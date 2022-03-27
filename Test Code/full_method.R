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

#Function to find the center peaks in each of the loops
findCenters <- function(sub_stats){
  delta = sub_stats[,2]
  rho = sub_stats[,3]
  a <- lm(log(delta)~log(rho))$coefficient[[2]]
  #Find which constant to use for the decision line ln(rho) + a*ln(delta) = const
  threshold <-  seq(0.001, 1, by = 0.001)
  clusters <- c()
  
  for (i in 1:length(threshold)){
    clusters[i] <- sum(delta > (threshold[i]*rho)^a)
  }
  library(pracma)
  #Project decisiongraph on regression line to find peaks and valleys
  b <- lm(clusters ~ threshold)$coefficient[[2]]
  new_coefs <- cbind(as.numeric(threshold), as.numeric(clusters - threshold*b))
  if (sum(new_coefs[,2]) == 0){ #There are no cluster centers
    center_id = c()
  }
  else{
    #Find all peaks and valleys
    p <- findpeaks(new_coefs[,2])
    v <- findpeaks(-new_coefs[,2])
    #Find the best candidate for threshold
    #If we start with a valley remove it for easier computation (only interested in peak vs next valley)
    #Also if we never have a peak, there are no clusters
    if(length(p) == 0){
      center_id = c()
    }
    
    else{
      if (p[1,2] > v[1,2]) { 
        p = rbind(c(0, 1, 0, 0), p)
      }
    #Add a minimum at the end if the curve was still going down
    
      if (nrow(p) > nrow(v)){
        v = rbind(v, c(0, length(threshold), 0, 0))
      }
      
      else{
      ind_pv <- which.max(v[, 2] - p[, 2])
      ind_bin <- (p[ind_pv, 2] + v[ind_pv,2])/2
      Threshold <- stats$coefs[[ind_bin, 1]]
      #Threshold found
      }
      center_id = sub_stats[which(delta > (Threshold*rho)^a),1]
      
    }
  return(center_id = center_id)
  }
}

#### Calculating Delta ###
#This one works perfect
find_delta_rho <- function(df, avgCellDist, k = 7, estimator = 2, track = FALSE){
  
  #Input df must be of size Nx2, with x and y columns
  
  L <- 3*avgCellDist #Large Grid Size
  rho = calc_rho(df, k = 7, estimator = estimator)
  #Order by rho
  sorted_rho <- cbind(df, rho, 1:length(rho)) %>% 
    .[order(.[,3], decreasing = TRUE), ]
  sorted <- cbind(1:length(rho), sorted_rho)
  
  #Statistics needed to construct grid
  xmin <- min(df[,1])
  xmax <- max(df[,1])
  ymin <- min(df[,2])
  ymax <- max(df[,2])
  xlength <- xmax - xmin
  ylength <- ymax - ymin
  grid_size <- L/3
  xgrids <- ceiling(xlength/grid_size)
  ygrids <- ceiling(ylength/grid_size)
  deltas <- c()
  
  for (i in 2:(ygrids-1)){ #For each row
    for (j in 2:(xgrids-1)){ #Check each grid
      large_grid <- sorted[sorted[, 2] < xmin + (j+1)*grid_size &
                             sorted[, 2] >= xmin + (j-2)*grid_size & 
                             sorted[, 3] < ymin + (i+1)*grid_size &
                             sorted[, 3] >= ymin + (i-2)*grid_size, ]
      
      small_grid <- sorted[sorted[, 2] < xmin + j*grid_size &
                             sorted[, 2] >= xmin + (j-1)*grid_size & 
                             sorted[, 3] < ymin + i*grid_size &
                             sorted[, 3] >= ymin + (i-1)*grid_size, ]
      
      if (is.vector(small_grid) == TRUE){ #Only one point within our smaller grid
        # It should be an outlier most likely so we don't need to consider it
      }
      else if (nrow(small_grid) < 6){
        # Dont add anything
      }
      #Atleast 6 points in the small grid
      else {
        sub_stats <- matrix(NA, ncol = 3, nrow = nrow(small_grid))
        #Calculate our deltas for the small grid
        for (k in 1:nrow(small_grid)){ #For every point in the small grid
          if (sum(large_grid[,1] < small_grid[k, 1]) == 0){ #There is no point with larger density
            deltas[small_grid[k, 5]] = L #Give max delta
            
            sub_stats[k, 1] <- small_grid[k, 5] 
            sub_stats[k, 2] <- L 
            sub_stats[k, 3] <- small_grid[k, 4]
          }
          
          else if (sum(large_grid[,1] < small_grid[k,1]) == 1){
            delta = sqrt((large_grid[large_grid[,1] < small_grid[k,1], ][2] - small_grid[k,2])^2 + (large_grid[large_grid[,1] < small_grid[k,1], ][3] - small_grid[k,3])^2)
            deltas[small_grid[k,5]] <- delta
            sub_stats[k, 1] <-small_grid[k, 5]
            sub_stats[k, 2] <- delta
            sub_stats[k, 3] <- small_grid[k, 4]
            
          }
          else{
            delta = min(apply(large_grid[large_grid[,1] < small_grid[k,1], ][,2:3], 1, function(x) sqrt(sum((x-small_grid[k,2:3])^2))))
            deltas[small_grid[k, 5]] <- delta
            sub_stats[k, 1] <- small_grid[k, 5]
            sub_stats[k, 2] <- delta
            sub_stats[k, 3] <- small_grid[k, 4]
          }
        }
        
      }
    }
    if (track == TRUE){
    print(i)
    }
  }
  return(list(delta = deltas, rho = rho))
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

#Find the threshold to use when deciding on nr of clusters
findThreshold <- function(stats){
  library(pracma)
  #Project decisiongraph on regression line to find peaks and valleys
  b <- lm(stats$coefs[,2] ~ stats$coefs[,1])$coefficient[[2]]
  new_coefs <- cbind(as.numeric(stats$coefs[,1]), as.numeric(stats$coefs[,2] - stats$coefs[,1]*b))
  #Find all peaks and valleys
  p <- findpeaks(new_coefs[,2])
  v <- findpeaks(-new_coefs[,2])
  #Find the best candidate for threshold
  #If we start with a valley remove it for easier computation (only interested in peak vs next valley)
  if (p[1,2] > v[1,2]) {
    v = v[-1, ]
  }
  #Remove last peak/valley if we have one more of peaks/valleys for correct dimensions.
  if (nrow(p) != nrow(v)){
    ind_pv <- which.max(v[-max(nrow(v), nrow(p)),2] - p[-max(nrow(v), nrow(p)),2])
  }
  else {
    ind_pv <- which.max(v[, 2] - p[, 2])
  }
  ind_bin <- (p[ind_pv, 2] + v[ind_pv,2])/2
  Threshold <- stats$coefs[[ind_bin, 1]]
  return(list(Threshold = Threshold, b = b, new_coefs = new_coefs))
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

