# Function to calculate true positive and false positive rate from true adjacency and estimated
# adjacency matrices (as a list of length K of binary matrices)
get_tpr_fpr_mcc <- function(truth, results) {
  K <- length(truth)
  
  # Get indices of upper triangular entries i.e. unique edges
  upperind <- which(upper.tri(truth[[1]]))
  
  # Get number of fp, tp, fn and tn
  fp <- 0
  tp <- 0
  fn <- 0
  tn <- 0
  for (k in 1:K) {
    cur_truth <- truth[[k]]
    true_up <- cur_truth[upperind]
    cur_res <- results[[k]]
    res_up <- cur_res[upperind]
    tp <- tp + sum(true_up & res_up)
    fp <- fp + sum(!true_up & res_up)
    tn <- tn + sum(!true_up & !res_up)
    fn <- fn + sum(true_up & !res_up)
  }
  
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  
  c(tpr, fpr, mcc)
}

# Function to calculate true positive and false positive rate for detection of differential edges
# using either fused graphical lasso or group graphical lasso
get_tpr_fpr_diff_jgl <- function(truth, true_C, method) {
  
  K <- length(truth)
  p <- nrow(truth[[1]])
  
  # For fused lasso, entries should be equal - we assume this is difference less than 10^-6
  if (method == "fused") {
    threshold <- 10^-6
  } else {
    # For group lasso, the criteria for difference in Danaher is abs diff > 10^-2
    threshold <- 10^-2
  }
  
  # Get number of fp, tp, fn and tn
  fp <- 0
  tp <- 0
  fn <- 0
  tn <- 0
  for (k1 in 1:(K - 1)) {
    for (k2 in (k1 + 1):K) {
      cur_truth1 <- truth[[k1]]
      cur_truth2 <- truth[[k2]]
      cur_theta1 <- true_C[[k1]]
      cur_theta2 <- true_C[[k2]] 
      
      for (i in 1:(p - 1)) {
        for (j in (i + 1): p) {
          true_diff <- abs(cur_truth1[i, j] - cur_truth2[i, j])
          if (method == "adj") {
            # Using adjacency matrix
            res_diff <- abs(cur_theta1[i, j] - cur_theta2[i, j])
          } else {
            theta_diff <- abs(cur_theta1[i, j] - cur_theta2[i, j])
            res_diff <- as.numeric(theta_diff > threshold)
          }
          
          tp <- tp + sum(true_diff & res_diff)
          fp <- fp + sum(!true_diff & res_diff)
          tn <- tn + sum(!true_diff & !res_diff)
          fn <- fn + sum(true_diff & !res_diff)
        }
      }
    }
  }
  
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  
  c(tpr, fpr, mcc)
}

# Function to compute Frobenius loss given true precision matrix + estimated precision matrix
get_fl <- function(true_C, est_C) {
  K <- length(true_C)
  fl <- 0
  
  for (k in 1:K) {
    fl <- fl + 
      norm(as.matrix(true_C[[k]]) - as.matrix(est_C[[k]]), type = 'F')^2 / 
      norm(as.matrix(true_C[[k]]), type = 'F')^2 / K
  }
  fl
}

get_perf_jgl <- function(true_adj, est_adj, true_C, est_C, method) {
  tpr_fpr_mcc_res <- get_tpr_fpr_mcc(true_adj, est_adj)
  tpr_fpr_mcc_diff_res <- get_tpr_fpr_diff_jgl(true_adj, est_C, method)
  fl <- get_fl(true_C, est_C)
  c(tpr_fpr_mcc_res, tpr_fpr_mcc_diff_res, fl)
}

get_auc_jgl <- function(param2, true_adj, penalty) {
  # Get AUC by increasing sparsity penalty until solution is all zeros
  
  if (penalty == "fused") {    
    # Possible values for lambda1s
    lambda1s <- (0:350) / 100
    
    # lmabda2 is fixed
    lambda2s <- rep(param2, length(lambda1s))
  } else {
    # Translate range of omega1s to lambda1s and lambda2s based on omega
    omega1s <- (0:250) / 100
    omega2 <- param2
    
    lambda1s <- omega1s * (1 - omega2)
    lambda2s <- sqrt(2) * omega1s * omega2    
  }
  
  # Keep track of TPR and FPR, and TPR and FPR for diff
  tpr_roc <- rep(NA, length(lambda1s))
  fpr_roc <- rep(NA, length(lambda1s))
  tpr_roc_diff <- rep(NA, length(lambda1s))
  fpr_roc_diff <- rep(NA, length(lambda1s))
  
  # Keep track of tpr and fpr at each iteration for regular edge detection and detection of difference
  for (i in 1:length(lambda1s)) { 
    jgl_res <- JGL(Y = Y, penalty = penalty, lambda1 = lambda1s[i], lambda2 = lambda2s[i],
                   return.whole.theta = TRUE)   
    
    # Get adjacency matrix from this
    est_adj <- make.adj.matrix(jgl_res$theta, separate = TRUE)
    
    # Get performance of edge learning
    perf_summary <- get_tpr_fpr_mcc(true_adj, est_adj)
    tpr_roc[i] <- perf_summary[1]
    fpr_roc[i] <- perf_summary[2]
    
    # Get performance for difference learning USING ADJACENCY MATRIX
    perf_summary <- get_tpr_fpr_diff_jgl(true_adj, est_adj, "adj")
    tpr_roc_diff[i] <- perf_summary[1]
    fpr_roc_diff[i] <- perf_summary[2]
  }
  
  auc <- sum((fpr_roc[-length(fpr_roc)] - fpr_roc[-1]) *
               (tpr_roc[-1] + tpr_roc[-length(tpr_roc)]) / 2)
  
  # Keep upper boundary of diff points - keep whichever portion of curve turns out to be better
  stop_ind <- min(which(fpr_roc_diff == max(fpr_roc_diff)))
  tpr_roc_diff1 <- tpr_roc_diff[1:stop_ind]
  fpr_roc_diff1 <- fpr_roc_diff[1:stop_ind]
  
  # Diff results need sorting
  fpr_tpr_diff1 <- data.frame(fpr = fpr_roc_diff1, tpr = tpr_roc_diff1)
  fpr_tpr_diff1 <- fpr_tpr_diff1[order(-fpr_tpr_diff1$fpr, -fpr_tpr_diff1$tpr), ]
  auc_diff1 <- sum((fpr_tpr_diff1$fpr[-nrow(fpr_tpr_diff1)] - fpr_tpr_diff1$fpr[-1]) *
                     (fpr_tpr_diff1$tpr[-1] + fpr_tpr_diff1$tpr[-nrow(fpr_tpr_diff1)]) / 2)
  
  start_ind <- min(which(fpr_roc_diff == max(fpr_roc_diff)))
  tpr_roc_diff2 <- tpr_roc_diff[start_ind:length(tpr_roc_diff)]
  fpr_roc_diff2 <- fpr_roc_diff[start_ind:length(tpr_roc_diff)]
  
  # Diff results need sorting
  fpr_tpr_diff2 <- data.frame(fpr = fpr_roc_diff2, tpr = tpr_roc_diff2)
  fpr_tpr_diff2 <- fpr_tpr_diff2[order(-fpr_tpr_diff2$fpr, -fpr_tpr_diff2$tpr), ]
  auc_diff2 <- sum((fpr_tpr_diff2$fpr[-nrow(fpr_tpr_diff2)] - fpr_tpr_diff2$fpr[-1]) *
                     (fpr_tpr_diff2$tpr[-1] + fpr_tpr_diff2$tpr[-nrow(fpr_tpr_diff2)]) / 2)
  
  if (auc_diff1 > auc_diff2) {
    auc_diff <- auc_diff1
    fpr_tpr_diff <- fpr_tpr_diff1
  } else {
    auc_diff <- auc_diff2
    fpr_tpr_diff <- fpr_tpr_diff2
  }
  
  # Compute and return AUC and AUC of differential edge detection, as well as full tpr and fpr results
  return(list(auc = auc, auc_diff = auc_diff,
              tpr_roc = tpr_roc, fpr_roc = fpr_roc,
              tpr_roc_diff = fpr_tpr_diff$tpr, fpr_roc_diff = fpr_tpr_diff$fpr))
}

