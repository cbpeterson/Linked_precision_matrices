# Number of simulated data sets
niter <- 25

# Number of methods compared
nmethod <- 5

# Number of performance metrics computed
nmetric <- 9

# Number of groups in simulated data sets
K <- 3

# Number of variables in simulated data sets
p <- 100

# Table 1 = table summarizing performance of edge selection and precision matrix estimatation
# FL = Frobenius loss
table1 <- data.frame(TPR = rep(NA, nmethod * 2), FPR = rep(NA, nmethod * 2),
                     MCC = rep(NA, nmethod * 2), AUC = rep(NA, nmethod * 2),
                     FL = rep(NA, nmethod * 2),
                     num_edges = rep(NA, nmethod * 2))

# Table 2 = table summarizing performance of differential edge selection
table2 <- data.frame(TPR_diff = rep(NA, nmethod * 2), FPR_diff = rep(NA, nmethod * 2),
                     MCC_diff = rep(NA, nmethod * 2), AUC_diff = rep(NA, nmethod * 2))

for (method in 1:nmethod) {
  if (method == 1) {
    # Fused graphical lasso
    dir <- "Simulation/Output_files/fgl/"
    bayesian <- FALSE
    method_name <- 'fgl'
  } else if (method == 2) {
    # Group graphical lasso
    dir <- "Simulation/Output_files/ggl/"
    bayesian <- FALSE
    method_name <- 'ggl'
  } else if (method == 3) {
    # Separate Bayesian estimation
    dir <- "Simulation/Output_files/sep/"
    bayesian <- TRUE
    method_name <- 'sep'
  } else if (method == 4) {
    # Joint Bayesian estimation
    dir <- "Simulation/Output_files/joint/"
    bayesian <- TRUE
    method_name <- 'joint'
  } else {
    # Linked precision matrix approach
    dir <- "Simulation/Output_files/lpm/"
    bayesian <- TRUE
    method_name <- 'lpm'
  }

  perf_across_iters <- matrix(NA, nrow = niter, ncol = nmetric)
  num_edges <- rep(NA, niter)
  for (iter in 1:niter) {
    # Get performance summary
    cur_filename <- paste(dir, "perf_iter", iter, "_n150_", method_name, "_sim_p100.csv", sep = "")
    cur_results <- as.numeric(read.csv(cur_filename, header = FALSE))
    perf_across_iters[iter, ] <- cur_results
    
    # Get number of edges included (either from estimated precision matrices for JGL methods,
    # or by thresholding the PPI matrices for the Bayesian methods)
    if (method == 1 || method == 2) {
      cur_filename <- paste(dir, "./Theta_", toupper(method_name), "_iter", iter, "_n150_sim_p100.csv", sep = "")
      cur_results <- as.matrix(read.csv(cur_filename, header = FALSE))
      # Count nonzero entries as edges
      cur_adj <- abs(cur_results) > 1e-6
    } else {
      cur_filename <- paste(dir, "PPI_iter", iter, "_n150_", method_name, "_sim_p100.csv", sep = "")
      cur_results <- as.matrix(read.csv(cur_filename, header = FALSE))
      # Count edges as entries with PPI > 0.5
      cur_adj <- cur_results > 0.5
    }
    
    # Avg number of edges over all 3 graphs = (sum(cur_adj) - 3 * p) / 2 / K
    num_edges[iter] <- (sum(cur_adj) - 3 * p) / 2 / K 
  }
  mean_perf <- colMeans(perf_across_iters)
  se_perf <- apply(perf_across_iters, 2, sd)  / sqrt(niter)
  mean_edges <- mean(num_edges)
  se_edges <- sd(num_edges) / sqrt(niter)
  
  table1[(method - 1) * 2 + 1, 1:(nmetric - 4)] <- mean_perf[-c(5:8)]
  table1[(method - 1) * 2 + 2, 1:(nmetric - 4)] <- se_perf[-c(5:8)]
  table1[(method - 1) * 2 + 1, 6] <- mean_edges
  table1[(method - 1) * 2 + 2, 6] <- se_edges
  table2[(method - 1) * 2 + 1, ] <- mean_perf[c(5:8)]
  table2[(method - 1) * 2 + 2, ] <- se_perf[c(5:8)]
}

# Get estimated mean of Phi
mean_Phi <- matrix(0, K, K)
for (iter in 1:niter) {
  cur_filename <- paste(dir, "Phi_", iter, "_n150_lpm_sim_p100.csv", sep = "")
  cur_results <- as.matrix(read.csv(cur_filename, header = FALSE))
  mean_Phi <- mean_Phi + cur_results / niter
}

# Simulation performance summary for edge selection and Frobenius loss of precision matrix estimation
row.names(table1) <- c("fgl", "fgl_se", "ggl", "ggl_se", "sep", "sep_se", "joint", "joint_se", "lpm", "lpm_se")
round(table1, 4)

# Simulation performance summary for differential edge selection
row.names(table2) <- c("fgl", "fgl_se", "ggl", "ggl_se", "sep", "sep_se", "joint", "joint_se", "lpm", "lpm_se")
round(table2, 4)

# Posterior estimated value of Phi
round(mean_Phi, 2)