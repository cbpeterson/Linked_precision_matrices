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

table_1 <- matrix(0,nrow = 30, ncol = 5)
colnames(table_1) <- c("","hpHC", "HC" , "MCI" , "AD")
regions <- seq(1:100)
lobes <- list(regions,Frontal.Lobe,Temporal.Lobe,Parietal.Lobe,Occipital.Lobe,Limbic.Cortex)
lobe_names <- c("All ROIs","Frontal","Temporal","Parietal","Occipital","Limbic")

test <- c(1,1,1,1)

for(lobe in 0:5) {
  ind <- lobe * 5 + 1
  lobe_name <- lobe_names[lobe + 1]
  table_1[ind,] <- c(lobe_name,"hpHC", "HC" , "MCI" , "AD")
  table_1[(ind + 1):(ind + 4),1] <- c("hpHC", "HC" , "MCI" , "AD")
  l <- lobes[[lobe + 1]]
  for(i in 1:4) {
    for(j in 1:4) {
      median_models.sum <- apply(median_models[l,l,c(i,j)],c(1,2),function(x) sum(all.equal(x,test[c(i,j)]) == TRUE))
      table_1[i + ind,j + 1] <- sum(median_models.sum[lower.tri(median_models.sum)])
    }
  }
}

table_1