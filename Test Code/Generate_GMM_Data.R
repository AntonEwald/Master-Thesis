

#### Simulate Data ####

#' A function to create a dataframe consisting of a GMM sample
#' comparable to some in-situ sample
#' @param clusters: Decides the number of cells
#' @param avg_genes: Average number of genes per cell
#' @param cellwidth: Used as component variance in covariance matrix (cell width in um)
#' @param tissue_x: Size of the tissue in mikrometers in x axis
#' @param tissue_y: Size of the tissue in mikrometers in y axis
#' @param Mode: 3 Modes (1 = fixt cellsize, 2 = only change size, 3 = change size and shape)
Generate_GMM_Data <- function(clusters, tissue_x, tissue_y, Cellwidth, Mode){
  library(MASS)
  library(emdbook)
  library(tidyverse)
  library(mvtnorm)

  mu <-cbind(runif(clusters, 0, tissue_x), runif(clusters, 0, tissue_y))
  variance <- matrix(NA, nrow = nrow(mu), ncol = 2)
  # Decide which gaussian to sample from (same probabilities)
  # Create Sample from GMM
  gene_coord <- matrix(NA, ncol = 2, nrow = 0)
  true_dens <- matrix(NA, ncol = 1, nrow = 0)
  samples <- c()
  for (i in 1:clusters){
    mu1 <- mu[i, ]
    eps1 <- runif(1, 3, Cellwidth)
    eps2 <- runif(1, 3, Cellwidth)
    # Varying covariance matrix.
    if (Mode == 1){
      Cov = matrix(c(sqrt(Cellwidth), 0, 0, sqrt(Cellwidth)), nrow = 2, ncol = 2)
      area = qnorm(0.95, 0, Cellwidth)*10 #area of the circle
    } else if (Mode == 2){
      Cov = matrix(c(eps1, 0, 0, eps1), nrow = 2, ncol = 2)
      area = (qnorm(0.95, 0, eps1)/2)^2*pi #area of the circle
    } else {
      Cov = matrix(c(eps1, 0, 0, eps2), nrow = 2, ncol = 2)
      area = qnorm(0.95, 0, eps1)/2 * qnorm(0.95, 0, eps2)/2 * pi
      
    }
    n <- max(30, round(area*2 + rnorm(1, 0, 10))) #Number of genes per cell is based on its area
    point <- mvrnorm(n = n, mu = mu1, Sigma = Cov)
    gene_coord <- rbind(gene_coord, point)
    true_dens <- c(true_dens, emdbook::dmvnorm(point, mu = mu1, Sigma = Cov))
    variance[i, 1] <- eps1
    variance[i, 2] <- eps2
    samples <- c(samples, rep(i, n))
  }
  
  return(list(Coordinates = gene_coord, Density = true_dens, Mu = mu, Samples = as.factor(samples), Variance = variance))
}
