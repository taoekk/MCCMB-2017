ppi_data <- read_delim("G:/guru13/Yandex disk/R_Projects/Applied Graph Analysis of Arabidopsis data/Data/ppi_data.txt", 
                       +     "\t", escape_double = FALSE, col_types = cols(Date_updated = col_character()), 
                       +     trim_ws = TRUE)

go_data <- read_csv("G:/guru13/Yandex disk/R_Projects/Applied Graph Analysis of Arabidopsis data/Data/go_data.csv")

library(dplyr)
library(igraph)
library(intergraph)
library(subgraphMining)
library(tidyr)
library(RColorBrewer)
source("https://bioconductor.org/biocLite.R")
biocLite("org.At.tair.db")


graph_ppi <- as.matrix(select(ppi_data, Locus_name, InteractorLocus_name)) %>%  
  graph_from_edgelist() %>%
  simplify(remove.multiple = TRUE, remove.loops = TRUE)
graph_ppi

#Vertex properties
V(graph_ppi)$cat <- go_data$cat
V(graph_ppi)$GO_term <- go_data$`GO term (GO ID)`
V(graph_ppi)$code <- go_data$code

#Extract components
is.connected(graph_ppi)
comps <- decompose(graph_ppi, mode = "weak")
table(sapply(comps, vcount))
# (or clusters <- clusters(graph_ppi), table(clusters$csize))
max_comp <- decompose.graph(graph_ppi, mode = "weak", max.comps = NA, min.vertices = 700)[[1]]
is.connected(max_comp)

#Plot giant component
gl <- layout.kamada.kawai(max_comp)
plot(max_comp, layout=gl, vertex.size=2, vertex.color = "royalblue", edge.color = "tan2",  vertex.label=NA, edge.arrow.size=0.1)


# Properties of giant component
average.path.length(max_comp)
diameter(max_comp)
reciprocity(max_comp)
transitivity(max_comp, type="global")
transitivity(max_comp, type="local")
edge_density(max_comp, loops = F)
triad_census(max_comp)

# Degree properties
deg1 <- degree(max_comp, mode="all")
hist(deg1, breaks=1:vcount(max_comp)-1, main="Histogram of node degree")
table(deg1)
deg.dist <- degree_distribution(max_comp, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

deg2 <- degree(max_comp, mode="in")
hist(deg2, breaks=1:vcount(max_comp)-1, main="Histogram of node degree")
table(deg2)
deg.dist2 <- degree_distribution(max_comp, cumulative=T, mode="in")
plot( x=0:max(deg2), y=1-deg.dist2, pch=19, cex=1.2, col="orange",
      xlab="in-Degree", ylab="Cumulative Frequency")


deg3 <- degree(max_comp, mode="out")
hist(deg3, breaks=1:vcount(max_comp)-1, main="Histogram of node degree")
table(deg3)
deg.dist3 <- degree_distribution(max_comp, cumulative=T, mode="out")
plot( x=0:max(deg3), y=1-deg.dist3, pch=19, cex=1.2, col="orange",
      xlab="in-Degree", ylab="Cumulative Frequency")

# Degree distribution log plotting
deg.dist.cum1 <- degree.distribution(max_comp, mode = "all")
d1 <- 1:max(deg1)-1
ind <- (deg.dist.cum1 != 0)
plot(d1[ind], deg.dist.cum1[ind], log="xy", col="royalblue",xlab=c("Log all-Degree"), ylab=c("Log-Intensity"),
     main="Log-Log Degree Distribution")

deg.dist.cum2 <- degree.distribution(max_comp, mode = "in")
d2 <- 1:max(deg2)-1
ind <- (deg.dist.cum2 != 0)
plot(d2[ind], deg.dist.cum2[ind], log="xy", col="royalblue",xlab=c("Log in-Degree"), ylab=c("Log-Intensity"),
     main="Log-Log Degree Distribution")

deg.dist.cum3 <- degree.distribution(max_comp, mode = "out")
d3 <- 1:max(deg3)-1
ind <- (deg.dist.cum3 != 0)
plot(d3[ind], deg.dist.cum3[ind], log="xy", col="royalblue",xlab=c("Log out-Degree"), ylab=c("Log-Intensity"),
     main="Log-Log Degree Distribution")

a.nn.deg.at <- graph.knn(max_comp,V(max_comp))$knn
plot(deg, a.nn.deg.at, log="xy",col="royalblue", xlab=c("Log Vertex Degree"),ylab=c("Log Average Neighbor Degree"))


# Vertex centrality
l <- layout.kamada.kawai(max_comp) 
hs <- hub_score(max_comp, weights=NA)$vector
plot(max_comp, layout = l, edge.arrow.size = .1, vertex.size=rescale(hs, 1, 7), main="Hubs", 
     vertex.color = "royalblue", edge.color = "tan2", vertex.label=NA, edge.arrow.size=0.1)

as <- authority_score(max_comp, weights=NA)$vector
plot(max_comp, layout = l, edge.arrow.size = .1, vertex.size=rescale(as, 1, 7), main="Authority", 
     vertex.color = "royalblue", edge.color = "tan2", vertex.label=NA, edge.arrow.size=0.1)

# PageRank
pg <- page_rank(max_comp, vids = V(max_comp), directed = TRUE, damping = 0.5)
pg <- as.vector(pg[[1]])
plot(max_comp, layout = l, edge.arrow.size = .1, vertex.size=rescale(as.vector(pg[[1]]), 1, 7), main="PageRank", 
     vertex.color = "royalblue", edge.color = "tan2", vertex.label=NA, edge.arrow.size=0.1)

# Edge betweeness centrality
eb <- edge.betweenness(max_comp) 
E(max_comp)[order(eb, decreasing=T)[1:10]]
plot(max_comp, layout = l, vertex.size=2, edge.arrow.size = .1, edge.width=rescale(as.vector(eb), 1, 20), main="Edge betweeness centrality", 
     vertex.color = "royalblue", edge.color = "tan2", vertex.label=NA)

# Cliques and K-cores analysis
net.sym <- as.undirected(max_comp, mode= "collapse")
cliques(net.sym)
table(sapply(cliques(net.sym), length))
largest_cliques(net.sym)

coreness <- graph.coreness(max_comp)
table(coreness)
maxCoreness <- max(coreness)
coreness_net <- max_comp
V(coreness_net)$cor <- coreness

# Articulation points
ap <- articulation_points(coreness_net)
ap1 <- V(coreness_net)[ap]
V(coreness_net)$shape = "circle"
V(coreness_net)[ap]$shape = "square"
plot(max_comp, edge.arrow.size = .1, vertex.color = V(max_comp)$color, vertex.size = 2.3, 
     vertex.label = ifelse(degree(max_comp) > 20, V(max_comp)$name, NA), vertex.label.cex = 0.5, edge.color = "black", layout = layout_with_kk)

display.brewer.pal(8, "YlOrRd")
pal <- brewer.pal(8, "YlOrRd")
plot(coreness_net, layout = l, edge.arrow.size = .1, vertex.size = 3.5, 
     vertex.color = ifelse(V(coreness_net)$cor  == 1, pal[1], 
                           ifelse(V(coreness_net)$cor  == 2, pal[2],
                                  ifelse(V(coreness_net)$cor  == 3, pal[3],
                                         ifelse(V(coreness_net)$cor  == 4, pal[4],
                                                ifelse(V(coreness_net)$cor  == 5, pal[5],
                                                       ifelse(V(coreness_net)$cor  == 6, pal[6],
                                                              ifelse(V(coreness_net)$cor  == 7, pal[7],
                                                                     ifelse(V(coreness_net)$cor  == 8, pal[8], "white")))))))),
     vertex.label.cex = 0.4, vertex.label=NA, edge.width=rescale(as.vector(eb), 1, 20),
     edge.color = "tan4", vertex.shape = V(coreness_net)$shape, main = "all k-cores" )














