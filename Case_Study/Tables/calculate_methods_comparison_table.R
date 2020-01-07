Frontal.Lobe <- seq(1:44)
Temporal.Lobe <- seq(45,60)
Parietal.Lobe <- seq(61,72)
Occipital.Lobe <- seq(73,88)
Limbic.Cortex <- seq(89,100)
regions <- seq(1:100)

## load linked precision results

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

## load single graph estimate results

ppi_SG.hpHC <- read.table("./Case_Study/Output/SingleGraph_hpHC.csv",sep = ",")
ppi_SG.HC <- read.table("./Case_Study/Output/SingleGraph_HC.csv",sep = ",")
ppi_SG.MCI <- read.table("./Case_Study/Output/SingleGraph_MCI.csv",sep = ",")
ppi_SG.AD <- read.table("./Case_Study/Output/SingleGraph_AD.csv",sep = ",")

median_model_SG.hpHC <- ppi_SG.hpHC > .5
median_model_SG.HC <- ppi_SG.HC > .5
median_model_SG.MCI <- ppi_SG.MCI > .5
median_model_SG.AD <- ppi_SG.AD > .5

diag(median_model_SG.hpHC) <- rep(0,100)
diag(median_model_SG.HC) <- rep(0,100)
diag(median_model_SG.MCI) <- rep(0,100)
diag(median_model_SG.AD) <- rep(0,100)

## load fused graphical lasso results

ppi_fused.hpHC <- read.table("./Case_Study/Output/Fused_hpHC.csv",sep = ",")
ppi_fused.HC <- read.table("./Case_Study/Output/Fused_HC.csv",sep = ",")
ppi_fused.MCI <- read.table("./Case_Study/Output/Fused_MCI.csv",sep = ",")
ppi_fused.AD <- read.table("./Case_Study/Output/Fused_AD.csv",sep = ",")

median_model_fused.hpHC <- ppi_fused.hpHC > .5
median_model_fused.HC <- ppi_fused.HC > .5
median_model_fused.MCI <- ppi_fused.MCI > .5
median_model_fused.AD <- ppi_fused.AD > .5

diag(median_model_fused.hpHC) <- rep(0,100)
diag(median_model_fused.HC) <- rep(0,100)
diag(median_model_fused.MCI) <- rep(0,100)
diag(median_model_fused.AD) <- rep(0,100)

## load joint graph estimate results

ppi_joint.hpHC <- read.table("./Case_Study/Output/hpHC_joint.csv",sep = ",")
ppi_joint.HC <- read.table("./Case_Study/Output/HC_joint.csv",sep = ",")
ppi_joint.MCI <- read.table("./Case_Study/Output/MCI_joint.csv",sep = ",")
ppi_joint.AD <- read.table("./Case_Study/Output/AD_joint.csv",sep = ",")

median_model_joint.hpHC <- ppi_joint.hpHC > .5
median_model_joint.HC <- ppi_joint.HC > .5
median_model_joint.MCI <- ppi_joint.MCI > .5
median_model_joint.AD <- ppi_joint.AD > .5

diag(median_model_joint.hpHC) <- rep(0,100)
diag(median_model_joint.HC) <- rep(0,100)
diag(median_model_joint.MCI) <- rep(0,100)
diag(median_model_joint.AD) <- rep(0,100)


###

med.mod.hpHC <- med.mod.HC <- med.mod.MCI <- med.mod.AD <- array(0,dim = c(100,100,4))
med.mod.hpHC[,,4] <- median_model.hpHC
med.mod.hpHC[,,2] <- median_model_SG.hpHC
med.mod.hpHC[,,1] <- median_model_fused.hpHC
med.mod.hpHC[,,3] <- median_model_joint.hpHC

med.mod.HC[,,4] <- median_model.HC
med.mod.HC[,,2] <- median_model_SG.HC
med.mod.HC[,,1] <- median_model_fused.HC
med.mod.HC[,,3] <- median_model_joint.HC

med.mod.MCI[,,4] <- median_model.MCI
med.mod.MCI[,,2] <- median_model_SG.MCI
med.mod.MCI[,,1] <- median_model_fused.MCI
med.mod.MCI[,,3] <- median_model_joint.MCI

med.mod.AD[,,4] <- median_model.AD
med.mod.AD[,,2] <- median_model_SG.AD
med.mod.AD[,,1] <- median_model_fused.AD
med.mod.AD[,,3] <- median_model_joint.AD

med.mods <- array(0,c(100,100,4,4))
med.mods[,,,1] <- med.mod.hpHC
med.mods[,,,2] <- med.mod.HC
med.mods[,,,3] <- med.mod.MCI
med.mods[,,,4] <- med.mod.AD

med.mods <- med.mods * 1

lobes <- list(regions,Frontal.Lobe,Temporal.Lobe,Parietal.Lobe,Occipital.Lobe,Limbic.Cortex)
lobe_names <- c("All ROIs","Frontal","Temporal","Parietal","Occipital","Limbic")

test <- c(1,1,1,1)

table_4 <- matrix(0,nrow = 25,ncol = 18)
colnames(table_4) <- c("Lobe","Method",rep("hpHC",4),rep("HC",4),rep("MCI",4),rep("AD",4))
table_4[1,] <- c("","",rep(c("Linked","Single","Fused","Shaddox"),4))
table_4[,2] <- c("",rep(c("Linked","Single","Fused","Shaddox"),6))

for(group in 1:4) {
  col.ind <- (group - 1)*4 + 2
  median_models <- med.mods[,,,group]
  print(sum(median_models[,,4]))
  for(lobe in 1:6) {
    l <- lobes[[lobe]]
    row.ind <- (lobe - 1)*4 + 1
    # table_4[row.ind,] <- row.fill
    for(i in 1:4) {
      for(j in 1:4) {
        if(i >= j) {
          median_models.sum <- apply(median_models[l,l,c(i,j)],c(1,2),function(x) sum(all.equal(x,test[c(i,j)]) == TRUE))
          table_4[row.ind + i,col.ind + j] <- sum(median_models.sum[lower.tri(median_models.sum)])
          table_4[row.ind + i,1] <- lobe_names[lobe]
        }
        else{
          table_4[row.ind + i,col.ind + j] <- ""
        }
      }
    }
  }
}

table_4

