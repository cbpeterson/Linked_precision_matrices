library(JGL)

# Source functions
source("./Helper_functions/jgl_helper_functions.R")
source("./Helper_functions/get_perf_jgl.r")

# Run fused and group lasso on simulated data, keeping track of performance and lambda value
K <- 3
niter <- 25
n <- 50

# Read in true precision matrices
A1 <- read.csv('./Simulation/Input_files/A1_sim_p100.csv', header = FALSE)
A2 <- read.csv('./Simulation/Input_files/A2_sim_p100.csv', header = FALSE)
A3 <- read.csv('./Simulation/Input_files/A3_sim_p100.csv', header = FALSE)
true_C <- list(A1, A2, A3)

# Compute true adjacency matrices based on these
adjTrue1 <-  1 * (A1 != 0);
adjTrue2 <-  1 * (A2 != 0);
adjTrue3 <-  1 * (A3 != 0);
true_adj <- list(adjTrue1, adjTrue2, adjTrue3)

# Lambda1 options (tuning parameter for graphical lasso penalty)
lambda1_opts <- seq(0.1, 0.5, by = 0.02)

# Lambda2 options (tuning parameter for fused or group lasso penalty)
lambda2_opts_fgl <- c(0.0004, seq(0.0008, 0.02, by = 0.0008))
lambda2_opts_ggl <- c(0, 0.0001, 0.0002, 0.0004, seq(0.0008, 0.016, by = 0.0008))

# Range of parameter values for fused and group graphical lasso to use when
# computing the auc.
lambda2s <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.225, 0.25, 0.5) 
omega2s <- c(0.5, 0.7, 1, 1.3, 1.6, 2, 2.5, 5, 7.5)

for (iter in 1:niter) {
  set.seed(13643 + iter)
  
  print(paste('Current iter =', iter))
  
  cur_filename <- paste('./Simulation/Input_files/Y1_iter', iter, '_n', n, '_lpm_p100.csv', sep = '')
  Y1 <- as.matrix(read.csv(cur_filename, header = FALSE))
  
  cur_filename <- paste('./Simulation/Input_files/Y2_iter', iter, '_n', n, '_lpm_p100.csv', sep = '')
  Y2 <- as.matrix(read.csv(cur_filename, header = FALSE))
  
  cur_filename <- paste('./Simulation/Input_files/Y3_iter', iter, '_n', n, '_lpm_p100.csv', sep = '')
  Y3 <- as.matrix(read.csv(cur_filename, header = FALSE))
  
  Y <- list(Y1, Y2, Y3)
  
  # Run fused joint graphical lasso
  lambda <- get_optimal_params(Y, "fused", lambda1_opts, lambda2_opts_fgl)
  fused_res <- JGL(Y = Y, penalty = "fused", lambda1 = lambda[1], lambda2 = lambda[2],
                   return.whole.theta = TRUE)
  
  # Write theta matrix to file (i.e. estimated precision matrices)
  cur_theta_filename <- paste('Theta_FGL_iter', iter, '_n', n, '_sim_p100.csv', sep = '')
  write.table(fused_res$theta, cur_theta_filename, quote = FALSE, row.names = FALSE,
              col.names = FALSE, sep = ",")
  
  # Get estimated adjacency matrix
  est_adj <- make.adj.matrix(fused_res$theta, separate = TRUE)
  
  # Get performance (tpr, fpr, mcc, tpr_diff, fpr_diff, mcc_diff, fl1, fl2)
  perf_summary <- get_perf_jgl(true_adj, est_adj, true_C, fused_res$theta, "fused")
  
  # Run fused graphical lasso for different values of lambda2 to get range of AUC values
  auc_fgl <- rep(NA, length(lambda2s))
  auc_diff_fgl <- rep(NA, length(lambda2s))
  for (l in 1:length(lambda2s)) {
    auc_fused_res <- get_auc_jgl(param2 = lambda2s[l], true_adj, penalty = "fused")
    auc_fgl[l] <- auc_fused_res$auc
    auc_diff_fgl[l] <- auc_fused_res$auc_diff
  }
  
  # Take best AUC as final result + write performance summary to disk
  cur_perf_filename <- paste('perf_iter', iter, '_n', n, '_fgl_sim_p100.csv', sep = '')
  write.table(t(c(perf_summary[1], perf_summary[2], perf_summary[3], max(auc_fgl),
                perf_summary[4], perf_summary[5], perf_summary[6], max(auc_diff_fgl),
                perf_summary[7])),
              cur_perf_filename, quote = FALSE, row.names = FALSE,
              col.names = FALSE, sep = ",")

  # Run group graphical lasso
  lambda <- get_optimal_params(Y, "group", lambda1_opts, lambda2_opts_ggl)
  group_res <- JGL(Y = Y, penalty = "group", lambda1 = lambda[1], lambda2 = lambda[2],
                   return.whole.theta = TRUE)
  
  # Write theta matrix to file (i.e. estimated precision matrices)
  cur_theta_filename <- paste('Theta_GGL_iter', iter, '_n', n, '_sim_p100.csv', sep = '')
  write.table(group_res$theta, cur_theta_filename, quote = FALSE, row.names = FALSE,
              col.names = FALSE, sep = ",")
  
  # Get adjacency matrix from this
  est_adj <- make.adj.matrix(group_res$theta, separate = TRUE)
  
  # Get performance (tpr, fpr, mcc, tpr_diff, fpr_diff, mcc_diff, fl1, fl2)
  perf_summary <- get_perf_jgl(true_adj, est_adj, true_C, group_res$theta, "group")
  
  # Run group graphical lasso with different values of omega2
  auc_ggl <- rep(NA, length(omega2s))
  auc_diff_ggl <- rep(NA, length(omega2s))
  for (l in 1:length(omega2s)) {
    print(l)
    auc_group_res <- get_auc_jgl(param2 = omega2s[l], true_adj, penalty = "group")
    auc_ggl[l] <- auc_group_res$auc
    auc_diff_ggl[l] <- auc_group_res$auc_diff
  }
  
  # Take best AUC as final result + write performance summary to disk
  cur_perf_filename <- paste('perf_iter', iter, '_n', n, '_ggl_sim_p100.csv', sep = '')
  write.table(t(c(perf_summary[1], perf_summary[2], perf_summary[3], max(auc_ggl),
                perf_summary[4], perf_summary[5], perf_summary[6], max(auc_diff_ggl),
                perf_summary[7])),
              cur_perf_filename, quote = FALSE, row.names = FALSE,
              col.names = FALSE, sep = ",")
}