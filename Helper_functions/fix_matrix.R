fix_matrix <- function(A, denom_factor = 1) {
  # Fix to ensure positive definiteness from Danaher et al
  # Divide each off-diagonal element by sum of absolute values of
  # off-diagonal elements in its row
  p <- nrow(A)
  
  for (cur_row in 1:p) {
    cur_sum <- sum(abs(A[cur_row, ])) - 1
    if (cur_sum != 0) {
      A[cur_row, ] <- A[cur_row, ] / (denom_factor * cur_sum)
    }
    
    # Make sure diagonal entries are still 1
    A[cur_row, cur_row] <- 1
  }
  
  # Final matrix is average of matrix with its transpose
  A <- (A + t(A)) / 2
}
