run_pca <- function(X){
  #Generating principal component analysis using SVD
  scaled_x =  scale(X, center = TRUE, scale = TRUE)
  k = min(ncol(X), nrow(X))
  x_svd = svd(scaled_x) #running base r's svd to get corresponding matrices
  #nicer way to get eigenvectors of X^TX without messy computation
  pc_matrix = scaled_x %*% x_svd$v[, 1:k, drop = FALSE]
  loading_matrix <- x_svd$v[, 1:k, drop = FALSE]
  sdevs <- x_svd$d[1:k] / sqrt(nrow(X) - 1)
  
  colnames(pc_matrix) = paste0("PC", 1:k)
  colnames(loading_matrix) = paste0("PC", 1:k)
  rownames(loading_matrix) = colnames(X)
  pc_output = list(sdev = sdevs, rotation = loading_matrix, x = pc_matrix)
  pc_output
}

run_lloyds <- function(X, k) {
  #pseudo code:
  # pick random k points as centroids
  # while not converged :
  # # assign each observation to the cluster with the nearest centroid
  # # compute new centroid as averages of the updated group
  
  dimensions = ncol(X)
  #calculating distance using |a-b|^2 = |a|^2 + |b|^2 - 2a dot b
  centroids <- matrix(runif(k * dimensions), nrow = k, ncol = dimensions)
  X_norms = rowSums(X^2)
  centroid_norms = rowSums(centroids^2)
  distances = outer(X_norms, centroid_norms, "+") - 2 * (X %*% t(centroids))
  min_cols <- apply(distances, 1, which.min) #getting clusters
  cluster_matrix = cbind(X, min_cols)
  total_diff = 1e100
  while(total_diff < 10){
    for(cluster in 1:dimensions){
      centroids[cluster,] = colMeans(cluster_matrix[cluster_matrix[, dimensions + 1] == cluster, ][,1:dimensions])
    }
    
  }
}
X = matrix(runif(10), nrow = 5)
run_lloyds(X, 2)
vect1 = c(0.1697951, 0.754010689)
vect2 = c(0.62922483, 0.4816428)
sum((vect1 - vect2)^2)
