#### Simulate Data ####

generate_data_varyingShape <- function(n, clusters, tissue_x, tissue_y){
  
  ## A function to create a dataframe consisting of a GMM sample
  ## comparable to some in-situ sample
  
  # :Param clusters: Decides the number of cells
  # :Param avg_genes: Average number of genes per cell
  # :Param cellwidth: Used as component variance in covariance matrix (cell width in um)
  # :Param tissue_x: Size of the tissue in mikrometers in x axis
  # :Param tissue_y: Size of the tissue in mikrometers in y axis
  # :Clust: different mean vectors (different gaussians in GMM)
  mu <-cbind(runif(clusters, 0, tissue_x), runif(clusters, 0, tissue_y))
  # Decide which gaussian to sample from (same probabilities)
  GMM_sample <- sample.int(clusters, n, replace = TRUE) %>% 
    table()
  # Create Sample from GMM
  gene_coord <- matrix(NA, ncol = 2, nrow = 0)
  true_dens <- matrix(NA, ncol = 1, nrow = 0)
  for (i in 1:length(GMM_sample)){
    k <- as.numeric(dimnames(GMM_sample)[[1]][i])
    mu1 <- mu[k, ]
    eps1 <- runif(1, 1, 6)
    eps2 <- runif(1, 1, 6)
    # Varying covariance matrix.
    Cov <- matrix(c(eps1, 0, 0, eps2), nrow = 2, ncol = 2)
    point <- mvrnorm(n = GMM_sample[i], mu = mu1, Sigma = Cov)
    gene_coord <- rbind(gene_coord, point)
    true_dens <- c(true_dens, dmvnorm(point, mean = mu1, sigma = Cov))
    variance <- c(variance, eps)
  }
  return(list(gene_coord, true_dens, mu, GMM_sample, variance))
}
