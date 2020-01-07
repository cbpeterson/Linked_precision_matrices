library(ggplot2)
library(reshape2)
library(grid)
library(plyr)

ppi.hpHC <- read.table("./Case_Study/Output/hpHC.csv",sep=',')
ppi.HC <- read.table("./Case_Study/Output/HC.csv",sep=',')
ppi.MCI <- read.table("./Case_Study/Output/MCI.csv",sep=',')
ppi.AD <- read.table("./Case_Study/Output/AD.csv",sep=',')

mm.hpHC <- ppi.hpHC > .5
mm.HC <- ppi.HC > .5
mm.MCI <- ppi.MCI > .5
mm.AD <- ppi.AD > .5

# create a vector with number of ROIs in each lobe
lobes.size <- c(44,16,12,16,12)
v.lines <- cumsum(lobes.size)

#create matrix of ppi's for each lobe, only including ppis > 0.50
ppi.hpHC <- ppi.hpHC * mm.hpHC
ppi.HC <- ppi.HC * mm.HC
ppi.MCI <- ppi.MCI * mm.MCI
ppi.AD <- ppi.AD * mm.AD

colnames(ppi.hpHC) <- seq(1:100)

ppi.hpHC.hm <- melt(as.matrix(ppi.hpHC))
colnames(ppi.hpHC.hm) <- c("Column","Row","PPI")

ggplot(data = ppi.hpHC.hm, aes(x=Column, y=Row, fill=PPI)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "white", high = "black", midpoint = .5) + 
  geom_vline(xintercept = v.lines+.38, colour = 'red') + 
  geom_hline(yintercept = v.lines+.38, colour = 'red') + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("High Performing Healthy Control")

grid.text("Frontal \n Lobe",x = unit(.033, "npc"),y = unit(0.25, "npc"),gp=gpar(fontsize=7))
grid.text("Temporal \n Lobe",x = unit(.033, "npc"),y = unit(0.50, "npc"),gp=gpar(fontsize=7))
grid.text("Parietal \n Lobe",x = unit(.033, "npc"),y = unit(0.61, "npc"),gp=gpar(fontsize=7))
grid.text("Occipital \n Lobe",x = unit(.033, "npc"),y = unit(0.72, "npc"),gp=gpar(fontsize=7))
grid.text("Limbic \n Cortex",x = unit(.033, "npc"),y = unit(0.84, "npc"),gp=gpar(fontsize=7))

ppi.HC.hm <- melt(as.matrix(ppi.HC))
colnames(ppi.HC.hm) <- c("Column","Row","PPI")

ggplot(data = ppi.HC.hm, aes(x=Column, y=Row, fill=PPI)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "white", high = "black", midpoint = .5) + 
  geom_vline(xintercept = v.lines+.38, colour = 'red') + 
  geom_hline(yintercept = v.lines+.38, colour = 'red') + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Healthy Control")

grid.text("Frontal \n Lobe",x = unit(.033, "npc"),y = unit(0.25, "npc"),gp=gpar(fontsize=7))
grid.text("Temporal \n Lobe",x = unit(.033, "npc"),y = unit(0.50, "npc"),gp=gpar(fontsize=7))
grid.text("Parietal \n Lobe",x = unit(.033, "npc"),y = unit(0.61, "npc"),gp=gpar(fontsize=7))
grid.text("Occipital \n Lobe",x = unit(.033, "npc"),y = unit(0.72, "npc"),gp=gpar(fontsize=7))
grid.text("Limbic \n Cortex",x = unit(.033, "npc"),y = unit(0.84, "npc"),gp=gpar(fontsize=7))

ppi.MCI.hm <- melt(as.matrix(ppi.MCI))
colnames(ppi.MCI.hm) <- c("Column","Row","PPI")

ggplot(data = ppi.MCI.hm, aes(x=Column, y=Row, fill=PPI)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "white", high = "black", midpoint = .5) + 
  geom_vline(xintercept = v.lines+.38, colour = 'red') + 
  geom_hline(yintercept = v.lines+.38, colour = 'red') + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Mild Cognitive Impairment")

grid.text("Frontal \n Lobe",x = unit(.033, "npc"),y = unit(0.25, "npc"),gp=gpar(fontsize=7))
grid.text("Temporal \n Lobe",x = unit(.033, "npc"),y = unit(0.50, "npc"),gp=gpar(fontsize=7))
grid.text("Parietal \n Lobe",x = unit(.033, "npc"),y = unit(0.61, "npc"),gp=gpar(fontsize=7))
grid.text("Occipital \n Lobe",x = unit(.033, "npc"),y = unit(0.72, "npc"),gp=gpar(fontsize=7))
grid.text("Limbic \n Cortex",x = unit(.033, "npc"),y = unit(0.84, "npc"),gp=gpar(fontsize=7))

ppi.AD.hm <- melt(as.matrix(ppi.AD))
colnames(ppi.AD.hm) <- c("Column","Row","PPI")

ggplot(data = ppi.AD.hm, aes(x=Column, y=Row, fill=PPI)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "white", high = "black", midpoint = .5) + 
  geom_vline(xintercept = v.lines+.38, colour = 'red') + 
  geom_hline(yintercept = v.lines+.38, colour = 'red') + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Alzheimer's Disease")

grid.text("Frontal \n Lobe",x = unit(.033, "npc"),y = unit(0.25, "npc"),gp=gpar(fontsize=7))
grid.text("Temporal \n Lobe",x = unit(.033, "npc"),y = unit(0.50, "npc"),gp=gpar(fontsize=7))
grid.text("Parietal \n Lobe",x = unit(.033, "npc"),y = unit(0.61, "npc"),gp=gpar(fontsize=7))
grid.text("Occipital \n Lobe",x = unit(.033, "npc"),y = unit(0.72, "npc"),gp=gpar(fontsize=7))
grid.text("Limbic \n Cortex",x = unit(.033, "npc"),y = unit(0.84, "npc"),gp=gpar(fontsize=7))
