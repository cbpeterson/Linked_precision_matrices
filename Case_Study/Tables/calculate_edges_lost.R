Frontal.Lobe <- seq(1:44)
Temporal.Lobe <- seq(45,60)
Parietal.Lobe <- seq(61,72)
Occipital.Lobe <- seq(73,88)
Limbic.Cortex <- seq(89,100)

ppi.hpHC <- read.table("./Case_Study/Output/hpHC.csv",sep = ",")
ppi.HC <- read.table("./Case_Study/Output/HC.csv",sep = ",")
ppi.MCI <- read.table("./Case_Study/Output/MCI.csv",sep = ",")
ppi.AD <- read.table("./Case_Study/Output/AD.csv",sep = ",")

median_model.hpHC <- ppi.hpHC > .5
median_model.HC <- ppi.HC > .5
median_model.MCI <- ppi.MCI > .5
median_model.AD <- ppi.AD > .5

diag(median_model.hpHC) <- rep(0,100)
diag(median_model.HC) <- rep(0,100)
diag(median_model.MCI) <- rep(0,100)
diag(median_model.AD) <- rep(0,100)

median_models <- array(0,c(100,100,4))
median_models[,,1] <- median_model.hpHC
median_models[,,2] <- median_model.HC
median_models[,,3] <- median_model.MCI
median_models[,,4] <- median_model.AD

median_model_colSums <- apply(median_models,3,colSums)
table_8 <- t(sapply(1:50,function(x) apply(median_model_colSums[c(2*x - 1,2* x),],2,sum)))

table_8


