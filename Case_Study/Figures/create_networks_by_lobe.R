library(igraph)

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

# tally the number of groups that have edge (i,j)
wts <- apply(median_models,c(1,2),sum)

median_model.hpHC <- median_model.hpHC * wts
median_model.HC <- median_model.HC * wts
median_model.MCI <- median_model.MCI * wts
median_model.AD <- median_model.AD * wts

###
# Frontal_Lobe
Frontal.Lobe <- seq(1:44)
Frontal <- c(Frontal.Lobe[which(Frontal.Lobe %% 2 == 0)],rev(Frontal.Lobe[which(Frontal.Lobe %% 2 == 1)]))

sample.graph <- graph_from_adjacency_matrix(median_model.hpHC[Frontal,Frontal],mode = "undirected",diag = F,weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.95*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "hpHC",vertex.label = c(seq(2,44,2),rev(seq(1,44,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.HC[Frontal,Frontal],mode = "undirected",diag = F,weighted = T)
l <- layout_in_circle(sample.graph)
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 2] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
right <- 2.95*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "HC",vertex.label = c(seq(2,44,2),rev(seq(1,44,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.MCI[Frontal,Frontal],mode = "undirected",diag = F,weighted = T)
l <- layout_in_circle(sample.graph)
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
right <- 2.95*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "MCI",vertex.label = c(seq(2,44,2),rev(seq(1,44,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.AD[Frontal,Frontal],mode = "undirected",diag = F,weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.95*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "AD",vertex.label = c(seq(2,44,2),rev(seq(1,44,2))))

# Temporal_Lobe
Temporal.Lobe <- seq(45,60)
Temporal <- c(Temporal.Lobe[which(Temporal.Lobe %% 2 == 0)],rev(Temporal.Lobe[which(Temporal.Lobe %% 2 == 1)]))

sample.graph <- graph_from_adjacency_matrix(median_model.hpHC[Temporal,Temporal],mode = "undirected",diag = F,weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.88*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "hpHC",vertex.label = c(seq(46,60,2),rev(seq(45,60,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.HC[Temporal,Temporal],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.88*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "HC",vertex.label = c(seq(46,60,2),rev(seq(45,60,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.MCI[Temporal,Temporal],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.88*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "MCI",vertex.label = c(seq(46,60,2),rev(seq(45,60,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.AD[Temporal,Temporal],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.88*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "AD",vertex.label = c(seq(46,60,2),rev(seq(45,60,2))))

# Parietal_Lobe 
Parietal.Lobe <- seq(61,72)
Parietal <- c(Parietal.Lobe[which(Parietal.Lobe %% 2 == 0)],rev(Parietal.Lobe[which(Parietal.Lobe %% 2 == 1)]))

sample.graph <- graph_from_adjacency_matrix(median_model.hpHC[Parietal,Parietal],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.83*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "hpHC",vertex.label = c(seq(62,72,2),rev(seq(61,72,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.HC[Parietal,Parietal],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.83*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "HC",vertex.label = c(seq(62,72,2),rev(seq(61,72,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.MCI[Parietal,Parietal],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.83*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "MCI",vertex.label = c(seq(62,72,2),rev(seq(61,72,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.AD[Parietal,Parietal],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.83*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "AD",vertex.label = c(seq(62,72,2),rev(seq(61,72,2))))

# Occipital_Lobe
Occipital.Lobe <- seq(73,88)
Occipital <- c(Occipital.Lobe[which(Occipital.Lobe %% 2 == 0)],rev(Occipital.Lobe[which(Occipital.Lobe %% 2 == 1)]))

sample.graph <- graph_from_adjacency_matrix(median_model.hpHC[Occipital,Occipital],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.87*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "hpHC",vertex.label = c(seq(74,88,2),rev(seq(73,88,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.HC[Occipital,Occipital],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.87*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "HC",vertex.label = c(seq(74,88,2),rev(seq(73,88,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.MCI[Occipital,Occipital],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.87*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "MCI",vertex.label = c(seq(74,88,2),rev(seq(73,88,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.AD[Occipital,Occipital],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.87*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "AD",vertex.label = c(seq(74,88,2),rev(seq(73,88,2))))

# Limbic 
Limbic.Cortex <- seq(89,100)
Limbic <- c(Limbic.Cortex[which(Limbic.Cortex %% 2 == 0)],rev(Limbic.Cortex[which(Limbic.Cortex %% 2 == 1)]))

sample.graph <- graph_from_adjacency_matrix(median_model.hpHC[Limbic,Limbic],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.84*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "hpHC",vertex.label = c(seq(90,100,2),rev(seq(89,100,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.HC[Limbic,Limbic],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.84*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "HC",vertex.label = c(seq(90,100,2),rev(seq(89,100,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.MCI[Limbic,Limbic],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.84*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "MCI",vertex.label = c(seq(90,100,2),rev(seq(89,100,2))))

sample.graph <- graph_from_adjacency_matrix(median_model.AD[Limbic,Limbic],mode = "undirected",diag = F, weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.84*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=10, vertex.color = 'white', 
     main = "AD",vertex.label = c(seq(90,100,2),rev(seq(89,100,2))))
