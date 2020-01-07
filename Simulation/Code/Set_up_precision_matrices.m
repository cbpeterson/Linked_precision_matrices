% Load required package
library(igraph)

% Set seed for random number generator
set.seed(62145)

% Add path to directory containing helper functions
addpath('./Helper_functions')

p <- 100

# Simulate a network with five communities, each a scale-free network on 20 nodes
# a1 = adjacency matrix for group 1
a1 <- matrix(0, nrow = p, ncol = p)
n_lobes <- 5
start_ind <- 1
stop_ind <- p / 5
for (cur_lobe in 1:n_lobes) {
  g1_cur_lobe <- barabasi.game(p / 5, power = 1.1, directed = FALSE, m = 2)
  a1_cur_lobe <- as.matrix(as_adjacency_matrix(g1_cur_lobe))
  a1[start_ind:stop_ind, start_ind:stop_ind] <- a1_cur_lobe
  start_ind <- start_ind + p / 5
  stop_ind <- stop_ind + p / 5
}

# Generate a symmetric precision (concentration) matrix c1 corresponding to a1
c1 <- matrix(NA, nrow = p, ncol = p)
for (i in 1:p) {
  for (j in i:p) {
    if (i == j) {
      c1[i, j] <- 1
    } else if (a1[i, j] == 1) {
      c1[i, j] <- runif(1, 0.4, 0.6) * sample(c(1,-1), 1)
      c1[j, i] <- c1[i, j]
    } else {
      c1[i, j] <- 0
      c1[j, i] <- 0
    }
  }
}
isSymmetric(c1)

# Add/remove 5 edges at random to get second graph
n_add_rem <- 5
zero_inds <- which(c1 == 0, arr.ind = TRUE)
nonzero_inds <- which((c1 - diag(1, p)) != 0, arr.ind = TRUE)
inds_to_add <- zero_inds[sample(nrow(zero_inds), n_add_rem), ]
inds_to_remove <- nonzero_inds[sample(nrow(nonzero_inds), n_add_rem), ]
c2 <- c1
for (i in 1:n_add_rem) {
  c2[inds_to_add[i, 1], inds_to_add[i, 2]] <- c1[inds_to_remove[i, 1], inds_to_remove[i, 2]]
  c2[inds_to_add[i, 2], inds_to_add[i, 1]] <- c1[inds_to_remove[i, 1], inds_to_remove[i, 2]]
  c2[inds_to_remove[i, 1], inds_to_remove[i, 2]] <- 0
  c2[inds_to_remove[i, 2], inds_to_remove[i, 1]] <- 0
}
isSymmetric(c2)

# Now adjust the first two precision matrices following Danaher et al approach
c1 <- fix_matrix(c1, 1.5)
c2 <- fix_matrix(c2, 1.5)

# To get third graph, remove 20 edges from graph 2
n_rem <- 20
nonzero_inds <- which((c2 - diag(1, p)) != 0, arr.ind = TRUE)
inds_to_remove <- nonzero_inds[sample(nrow(nonzero_inds), n_rem), ]
c3 <- c2
for (i in 1:n_rem) {
  c3[inds_to_remove[i, 1], inds_to_remove[i, 2]] <- 0
  c3[inds_to_remove[i, 2], inds_to_remove[i, 1]] <- 0
}
isSymmetric(c3)

# Check that all eigenvalues are positive i.e. the matrices are positive definite
min(eigen(c1)$values)
min(eigen(c2)$values)
min(eigen(c3)$values)

# Write these to disk
write.table(c1, 'A1_sim_p100.csv', row.names = FALSE, quote = FALSE, col.names = FALSE, sep = ",")
write.table(c2, 'A2_sim_p100.csv', row.names = FALSE, quote = FALSE, col.names = FALSE, sep = ",")
write.table(c3, 'A3_sim_p100.csv', row.names = FALSE, quote = FALSE, col.names = FALSE, sep = ",")





