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
#ground truth: summary(prcomp(X, , center = TRUE, scale. = TRUE))$x
