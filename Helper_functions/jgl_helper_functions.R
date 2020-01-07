# Helper function internal to JGL package
make.adj.matrix <- function(theta, separate = FALSE) {
  K = length(theta)
  adj = list()
  if (separate) {
    for(k in 1:K)
    {
      adj[[k]] = (abs(theta[[k]]) > 1e-5) * 1
    }
  }
  if(!separate) {
    adj = 0 * theta[[1]]
    for(k in 1:K) {
      adj = adj + (abs(theta[[k]]) > 1e-5)*2^(k-1)
    }
  }
  return(adj)
}

# Create function to calculate AIC
get_aic <- function(lambda1, lambda2, Y, penalty) {
  K <- length(Y)
  AIC <- 0  
  
  jgl_res <- JGL(Y = Y, penalty = penalty, lambda1 = lambda1, lambda2 = lambda2,
                 return.whole.theta = TRUE)
  theta <- jgl_res$theta
  
  for (k in 1:K) {
    # Get data for current sample group
    cur_Y <- Y[[k]]
    n_k <- nrow(cur_Y)
    E_k <- length(net.edges(theta)[[k]])
    S <- 1 / n_k * (t(cur_Y) %*% cur_Y)
    AIC <- AIC + n_k * sum(diag(S %*% theta[[k]])) - n_k * log(det(theta[[k]])) + 2 * E_k
  }
  
  AIC
}

# Function to do grid search to find optimal params using AIC - this is slow!
get_optimal_params <- function(Y, penalty, lambda1_opts, lambda2_opts) {
  AIC_matrix <- matrix(NA, length(lambda1_opts), length(lambda2_opts))
  for (i in 1:nrow(AIC_matrix)) {
    for (j in 1:ncol(AIC_matrix)) {
      AIC_matrix[i, j] <- get_aic(lambda1_opts[i], lambda2_opts[j], Y, penalty)
    }
  }
  
  # What if the min is not unique? Favor largest lambda1, then largest lambda2 since
  # they recommend search over lamba1 first, then lambda2 when necessary
  indices <- which(AIC_matrix == min(AIC_matrix), arr.ind = TRUE)
  
  lambda1_index <- max(indices[ , "row"])
  lambda2_index <- max(which(AIC_matrix[lambda1_index, ] == min(AIC_matrix)))
  
  lambda1 <- lambda1_opts[lambda1_index]
  lambda2 <- lambda2_opts[lambda2_index]
  
  # Check to make sure we don't need to expand grid
  if (lambda1 == max(lambda1_opts) || lambda2 == max(lambda2_opts)) {
    print("Max lambda value selected: need to expand grid")
  }
  
  c(lambda1, lambda2)
}