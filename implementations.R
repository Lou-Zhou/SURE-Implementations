library(ggplot2)
#PCA

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

#K-Means

run_lloyds <- function(X, k, random_function = complete_random) {
  #pseudo code:
  # pick random k points as centroids
  # while not converged :
  # # assign each observation to the cluster with the nearest centroid
  # # compute new centroid as averages of the updated group
  scaled_x =  scale(X, center = TRUE, scale = TRUE)
  dimensions = ncol(scaled_x)
  centroids <- random_function(scaled_x, k)
  cluster_matrix <- generate_clusters(scaled_x, centroids)
  same_clusters = FALSE
  while(!same_clusters){
    for(cluster in 1:k){
      #generating new centroids
      cluster_rows <- cluster_matrix[, "min_cols"] == cluster
      if(any(cluster_rows)){#recalculate centroid if there are members of the cluster
        centroids[cluster, ] <- colMeans(cluster_matrix[cluster_rows, 1:dimensions, drop = FALSE])
      }
    }
    new_clusters <- generate_clusters(scaled_x,  centroids)
    same_clusters = all(cluster_matrix[,"min_cols"] == new_clusters[,"min_cols"])
    cluster_matrix = new_clusters
  }
  cluster_matrix
}
generate_clusters <- function(X_scaled, centroids) {
  #calculating distance using |a-b|^2 = |a|^2 + |b|^2 - 2a dot b
  distances <- get_pairwise_distances(X_scaled, centroids)
  min_cols <- apply(distances, 1, which.min)
  cluster_matrix = cbind(X_scaled, min_cols)
  cluster_matrix
}
get_pairwise_distances <- function(X, Y){
  #Euclidean Distance Matrix between X and Y
  X_norms <- rowSums(X^2)
  Y_norms <- rowSums(Y^2)
  outer(X_norms, Y_norms, "+") - 2 * (X %*% t(Y))
}
complete_random <- function(X, k) {
  #randomly generates centroids within bounds of each dimension
  mins <- apply(X, 2, min)
  maxs <- apply(X, 2, max)
  matrix(runif(k * ncol(X), min = rep(mins, each = k), max = rep(maxs, each = k)),
         nrow = k, byrow = FALSE)
}
kpp <- function(X, k){
  #Using k++ method of initializing centroids
  centroids = matrix(nrow = 1, ncol = ncol(X))
  centroids[1,] <- X[sample(nrow(X),size=1,replace=TRUE),]
  for(centroid in 2:k + 1){
    distances <- get_pairwise_distances(X, centroids)
    min_cols <- apply(distances, 1, min)
    probabilities <- min_cols^2 / sum(min_cols^2)
    random_point  <- sample(1:nrow(X), size = 1, prob = probabilities)
    centroids <- rbind(centroids, X[random_point,])
  }
  centroids
}

#Hierarchical

get_cluster_distance <- function(cluster1, cluster2, distance_matrix, method = "average") {
  dists <- distance_matrix[unlist(cluster1), unlist(cluster2), drop = FALSE]
  if (method == "average") {
    return(mean(dists))
  } else if (method == "single") {
    return(min(dists))
  } else if (method == "complete") {
    return(max(dists))
  } else{
    print("ERROR: Invalid Method")
  }
}
check_id <- function(cluster) {
  #if singleton
  if((length(cluster)) == 1){
    return(-cluster[1])
  }
  return(attr(cluster, 'cluster_step'))
}
hierach_cluster <- function(X, method = "average") {
  #Runs hierarchical clustering from base R
  #TODO: Test this!
  n = nrow(X)
  merge_mat = matrix(NA, ncol = 2, nrow = n - 1)
  heights = seq(-1, n - 1)
  distances = get_pairwise_distances(X, X)
  clusters = lapply(1:n, function(i) i )
  #cluster of just singular points
  for(step in 1:(n-1)){
    min_dist = Inf
    num_clusters = length(clusters)
    cluster_merge = c(NA, NA)
    #find minimum distance by iterating through all
    for(idx_1 in 1:(num_clusters - 1)) {
      for(idx_2 in (idx_1 + 1):(num_clusters)) {
        distance <- get_cluster_distance(clusters[[idx_1]], clusters[[idx_2]], distances, method)
        if(distance < min_dist) {
          min_dist = distance
          clusters_merge = c(idx_1, idx_2)
        }
      }
    }
    #log merge
    idx_1 <- clusters_merge[1]
    idx_2 <- clusters_merge[2]
    merge_mat[step,] <- c(check_id(clusters[[idx_1]]), check_id(clusters[[idx_2]]))
    heights[step] <- min_dist
    
    #merge the two clusters
    new_cluster = c(clusters[[idx_1]], clusters[[idx_2]])
    attr(new_cluster, 'cluster_step') = step
    clusters <- clusters[-c(idx_1, idx_2)]
    
    #append new cluster
    clusters[[length(clusters) + 1]] <- new_cluster
  }
  list(merge = merge_mat, heights = heights)
}

# KDE

sqrt_matrix <- function(X){
  e_vectors <- eigen(X)$vectors
  e_values <- eigen(X)$values
  D_sqrt <- if (length(e_values) == 1) {
    matrix(sqrt(e_values), 1, 1)
  } else {
    diag(sqrt(e_values))
  }
  return(e_vectors %*% D_sqrt %*% solve(e_vectors))
}

kde <- function(x, data, H, K) {
  H_inv = solve(H)
  differences <- sweep(data, 2, x, FUN = "-")
  summation <- sum(apply(differences, 1, function(row)  K(solve(sqrt_matrix(H)) %*% row)))
  return(1/(nrow(data) * sqrt(det(H))) * summation)
}
gaussian_kernel <- function(z) {
  d <- length(z)
  
  const <- (2 * pi)^(-d/2)
  exponent <- (-1/2) *(t(z) %*% z)
  
  return(const * exp(exponent))
}
