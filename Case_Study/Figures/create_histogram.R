library(plotrix)

ppi.hpHC <- read.table("./Case_Study/Output/hpHC.csv",sep = ",")
ppi.HC <- read.table("./Case_Study/Output/HC.csv",sep = ",")
ppi.MCI <- read.table("./Case_Study/Output/MCI.csv",sep = ",")
ppi.AD <- read.table("./Case_Study/Output/AD.csv",sep = ",")

lower_tri <- lower.tri(ppi.hpHC)
ppi.hpHC <- ppi.hpHC[lower_tri]
ppi.HC <- ppi.HC[lower_tri]
ppi.MCI <- ppi.MCI[lower_tri]
ppi.AD <- ppi.AD[lower_tri]

# list of PPI's for each group
ppis <- list(ppi.hpHC,ppi.HC,ppi.MCI,ppi.AD)

groups = c("hpHC" , "HC", "MCI" , "AD")
breaks = seq(0,1,.025)
ylim.breaks = seq(0,2250,250)

dev.off()
par(mfrow = c(4,4))
for(i in 1:4) {
  for(j in 1:4) {
    if(i == j) {
      breaks = seq(0,1,.05)
      h <- hist(ppis[[i]],breaks = breaks, plot = F)
      h$counts[1] <- h$counts[1] - 2000
      plot(h,ylim = c(0,2500), axes = F , ylab = "Percentage" , xlab = "PPI" , main =  paste("Histogram of PPI for ",groups[i], sep = "") , cex.main = .8)
      axis(side=1, at=breaks, labels=breaks)
      axis(side=2, at=c(0,500,1500,2000,2500), labels=c(0,500,3500,4000,4500)/5000)
      axis.break(axis = 2, breakpos = 1000,brw = .05)
    }
    else if(i > j) {
      q1 <- sum(ppis[[j]] > 0.5 & ppis[[i]] > 0.5)
      q2 <- sum(ppis[[j]] < 0.5 & ppis[[i]] > 0.5)
      q3 <- sum(ppis[[j]] < 0.5 & ppis[[i]] < 0.5)
      q4 <- sum(ppis[[j]] > 0.5 & ppis[[i]] < 0.5)
      plot(0,0,cex=0,ylim = c(0,1),xlim=c(0,1),axes = T,xlab = paste("PPIs for",groups[j]) , ylab = paste("PPIs for",groups[i]))
      abline(h=.5,v=.5)
      text(c(.75,.25,.25,.75),c(.75,.75,.25,.25),paste(round((c(q1,q2,q3,q4)/sum(q1,q2,q4,q3))*100,2),"%",sep=""),cex=1.5)
    }
    else{
      plot(ppis[[i]],ppis[[j]],xlab = paste("PPIs for",groups[i]) , ylab = paste("PPIs for",groups[j]))
      abline(v=.5,h=.5)
    }
  }
}
