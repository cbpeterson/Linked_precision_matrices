library(JGL)
source("./Helper_functions/jgl_helper_functions.R")

hpHC <- read.table('./Case_Study/Input/fake_hpHC.csv')
HC <- read.table('./Case_Study/Input/fake_HC.csv')
MCI <- read.table('./Case_Study/Input/fake_MCI.csv')
AD <- read.table('./Case_Study/Input/fake_AD.csv')

AD.data <- list(as.matrix(hpHC),as.matrix(HC),as.matrix(MCI),as.matrix(AD))


max_lambda1 <- 0.02
lambda1_opts <- seq(0, max_lambda1, by = 0.0001)

# Lambda2 options (tuning parameter for fused or group lasso penalty)
# Orginal version: lambda2_opts <- seq(0, 0.025, by = 0.0005)
# Mostly 0 values were getting chosen, so try using grid with values closer to 0
max_lambda2 <- 0.0001
lambda2_opts <- seq(0, max_lambda2, by = 0.00001)

lambda <- get_optimal_params(AD.data, "fused", lambda1_opts, lambda2_opts)
fused <- JGL(Y = AD.data, penalty = "fused", lambda1 = lambda[1], lambda2 = lambda[2],
                 return.whole.theta = TRUE)
theta <- fused$theta
adj <- make.adj.matrix(fused$theta,separate = T)

write.table(adj[[1]],'./Case_Study/Output/Fused_hpHC.csv',row.names = F,col.names = F)
write.table(adj[[2]],'./Case_Study/Output/Fused_HC.csv',row.names = F,col.names = F)
write.table(adj[[3]],'./Case_Study/Output/Fused_MCI.csv',row.names = F,col.names = F)
write.table(adj[[4]],'./Case_Study/Output/Fused_AD.csv',row.names = F,col.names = F)
