library(brainGraph)
library(igraph)
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

## create graph objects

g1 <- graph_from_adjacency_matrix(median_model.hpHC,mode = "undirected")
g2 <- graph_from_adjacency_matrix(median_model.HC,mode = "undirected")
g3 <- graph_from_adjacency_matrix(median_model.MCI,mode = "undirected")
g4 <- graph_from_adjacency_matrix(median_model.AD,mode = "undirected")

## remove nodes that are disconnected from main network

g1 <- induced_subgraph(g1,components(g1)$membership == which.max(components(g1)$csize))
g2 <- induced_subgraph(g2,components(g2)$membership == which.max(components(g2)$csize))
g3 <- induced_subgraph(g3,components(g3)$membership == which.max(components(g3)$csize))
g4 <- induced_subgraph(g4,components(g4)$membership == which.max(components(g4)$csize))

g1$Cp <- transitivity(g1,c('average'))
g2$Cp <- transitivity(g2,c('average'))
g3$Cp <- transitivity(g3,c('average'))
g4$Cp <- transitivity(g4,c('average'))

g1$Lp <- mean_distance(g1)
g2$Lp <- mean_distance(g2)
g3$Lp <- mean_distance(g3)
g4$Lp <- mean_distance(g4)

g1$density <- g2$density <- g3$density <- g4$density <- NA

g1.rand <- sim.rand.graph.par(g1,N=1000,clustering = F)
g2.rand <- sim.rand.graph.par(g2,N=1000,clustering = F)
g3.rand <- sim.rand.graph.par(g3,N=1000,clustering = F)
g4.rand <- sim.rand.graph.par(g4,N=1000,clustering = F)

sw.hpHC <- small.world(g1,g1.rand)
sw.HC <- small.world(g2,g2.rand)
sw.MCI <- small.world(g3,g3.rand)
sw.AD <- small.world(g4,g4.rand)

sw.list <- list(sw.hpHC,sw.HC,sw.MCI,sw.AD)

table_3 <- matrix(0,nrow = 3, ncol = 4)
for(i in 1:4) {
  table_3[1,i] <- sw.list[[i]]$Lp.norm
  table_3[2,i] <- sw.list[[i]]$Cp.norm
  table_3[3,i] <- sw.list[[i]]$sigma
}

table_3
