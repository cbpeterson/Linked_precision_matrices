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

median_models <- array(0,c(100,100,4))
median_models[,,1] <- median_model.hpHC
median_models[,,2] <- median_model.HC
median_models[,,3] <- median_model.MCI
median_models[,,4] <- median_model.AD

test <- list()
test[[1]] <- c(1,0,0,0)
test[[2]] <- c(0,1,0,0)
test[[3]] <- c(0,0,1,0)
test[[4]] <- c(0,0,0,1)

lobes <- list(regions,Frontal.Lobe,Temporal.Lobe,Parietal.Lobe,Occipital.Lobe,Limbic.Cortex)
table_2 <- matrix(0,ncol=4,nrow=6)
tables <- list()

for(i in 1:4) {
  for(j in 1:6) {
    l <- lobes[[j]]
    median_models.sum <- apply(median_models[l,l,],c(1,2),function(x) sum(all.equal(x,test[[i]]) == TRUE))
    table_2[j,i] <- sum(median_models.sum[lower.tri(median_models.sum)])
  }
}

table_2
