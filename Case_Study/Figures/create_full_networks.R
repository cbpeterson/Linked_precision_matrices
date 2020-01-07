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

wts <- apply(median_models,c(1,2),sum)

median_model.hpHC <- median_model.hpHC * wts
median_model.HC <- median_model.HC * wts
median_model.MCI <- median_model.MCI * wts
median_model.AD <- median_model.AD * wts

#
regions <- seq(1,100)
region <- c(regions[which(regions %% 2 == 0)],rev(regions[which(regions %% 2 == 1)]))

sample.graph <- graph_from_adjacency_matrix(median_model.hpHC[region,region],mode = "undirected",diag = F, weighted = TRUE)
l <- layout_in_circle(sample.graph)
right <- 2.98*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
E(sample.graph)$width <- 1
E(sample.graph)$width[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$width[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=5, vertex.color = 'white', 
     main = "hpHC",vertex.label = c(seq(2,100,2),rev(seq(1,100,2))),
     vertex.label.cex = .5)

sample.graph <- graph_from_adjacency_matrix(median_model.HC[region,region],mode = "undirected",diag = F,weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.98*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
E(sample.graph)$width <- 1
E(sample.graph)$width[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$width[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=5, vertex.color = 'white', 
     main = "HC",vertex.label = c(seq(2,100,2),rev(seq(1,100,2))),
     vertex.label.cex = .5)

sample.graph <- graph_from_adjacency_matrix(median_model.MCI[region,region],mode = "undirected",diag = F,weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.98*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
E(sample.graph)$width <- 1
E(sample.graph)$width[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$width[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=5, vertex.color = 'white', 
     main = "MCI",vertex.label = c(seq(2,100,2),rev(seq(1,100,2))),
     vertex.label.cex = .5)

sample.graph <- graph_from_adjacency_matrix(median_model.AD[region,region],mode = "undirected",diag = F,weighted = T)
l <- layout_in_circle(sample.graph)
right <- 2.98*pi/2
rotation <- matrix(c(cos(right),sin(right),-sin(right),cos(right)),ncol=2)
l <- l %*% rotation
E(sample.graph)$color <- "lightblue"
E(sample.graph)$color[E(sample.graph)$weight == 1] <- "red"
E(sample.graph)$color[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- "black"
E(sample.graph)$lty <- 5
E(sample.graph)$lty[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$lty[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
E(sample.graph)$width <- 1
E(sample.graph)$width[E(sample.graph)$weight == 1] <- 3
E(sample.graph)$width[E(sample.graph)$weight == 2 | E(sample.graph)$weight == 3] <- 1
plot(sample.graph, layout=l, vertex.size=5, vertex.color = 'white', 
     main = "AD",vertex.label = c(seq(2,100,2),rev(seq(1,100,2))),
     vertex.label.cex = .5)
