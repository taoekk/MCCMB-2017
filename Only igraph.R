ppi_data <- read_delim("G:/guru13/Yandex disk/R_Projects/Applied Graph Analysis of Arabidopsis data/Data/ppi_data.txt", 
                       +     "\t", escape_double = FALSE, col_types = cols(Date_updated = col_character()), 
                       +     trim_ws = TRUE)

go_data <- read_csv("G:/guru13/Yandex disk/R_Projects/Applied Graph Analysis of Arabidopsis data/Data/go_data.csv")

# GO anthologies, under construction
#full_go <- read_delim("G:/guru13/Yandex disk/R_Projects/Applied Graph Analysis of Arabidopsis data/Data/ATH_GO_GOSLIM.txt", 
#                      +     "\t", escape_double = FALSE, col_names = FALSE, 
#                      +     trim_ws = TRUE)
#colnames(full_go) <- c("locus_name", "TAIR_accession", "object_name", "relationship_type", "GO_term", "GO_ID", "TAIR_Keyword_ID", "Aspect", "GOslim_term", "Evidence_code", "Evidence_description", "Evidence_with", "Reference", "Annotator", "Date_annotated")
#
#spread_go <- unique(full_go[, c(1, 4, 5)]) %>% 
#                    group_by(locus_name, relationship_type) 

#tmp <- paste(spread_go$Class_a, dd_1$Class_b, sep='-')

#%>%
#                    dcast(locus_name ~ relationship_type)                    
#  group_by(locus_name, relationship_type)
#    spread(relationship_type, GO_term, GOslim_term )


#group_by(spread_go, locus_name, relationship_type)                              

#123 <-paste(as.vector(spread_go$GO_term), as.vector(spread_go$GO_ID), sep = ";")  
#spread_go %>%
#  group_by(locus_name, relationship_type) %>% 
#  summarise(GO_term=toString(unique(GO_term))) %>% 
#  spread(relationship_type, GO_term, fill='')

library(dplyr)
library(igraph)
library(intergraph)
library(subgraphMining)
library(tidyr)
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

# Fullgraph image
l <- layout.kamada.kawai(graph_ppi)
plot(graph_ppi, layout=l, vertex.size=2, vertex.color = "royalblue",  vertex.label=NA, edge.arrow.size=0.1)

#Extract components
full_comps <- decompose(graph_ppi, mode = "weak", max.comps = NA, min.vertices = 10)
max_comp <- decompose(graph_ppi, mode = "weak", max.comps = NA, min.vertices = 100)
clusters <- clusters(graph_ppi)
table(clusters$csize) # Clusters enumeration

#
max_comp <- max_comp[[1]]
edge_density(max_comp, loops = F)
diameter(max_comp)
reciprocity(max_comp)
transitivity(max_comp, type="global")
transitivity(max_comp, type="local")
triad_census(max_comp)

deg <- degree(max_comp, mode="all")
hist(deg, breaks=1:vcount(max_comp)-1, main="Histogram of node degree")
table(deg)
deg.dist <- degree_distribution(max_comp, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

# Degree distribution log plotting
deg.dist.cum <- degree.distribution(max_comp)
d <- 1:max(deg)-1
ind <- (deg.dist.cum != 0)
plot(d[ind], deg.dist.cum[ind], log="xy", col="blue",xlab=c("Log-Degree"), ylab=c("Log-Intensity"),
     main="Log-Log Degree Distribution")

a.nn.deg.at <- graph.knn(max_comp,V(max_comp))$knn
plot(deg, a.nn.deg.at, log="xy",col="goldenrod", xlab=c("Log Vertex Degree"),ylab=c("Log Average Neighbor Degree"))

degree(max_comp, mode="in")
centr_degree(max_comp, mode="in", normalized=T)

closeness(max_comp, mode="all", weights=NA)
centr_clo(max_comp, mode="all", normalized=T)

centr_eigen(max_comp, directed=T, normalized=T)

betweenness(max_comp, directed=T, weights=NA)
edge_betweenness(max_comp, directed=T, weights=NA)
centr_betw(max_comp, directed=T, normalized=T)

# Hubs and Authorities

hs <- hub_score(max_comp, weights=NA)$vector
as <- authority_score(max_comp, weights=NA)$vector
plot(max_comp, edge.arrow.size = .1, vertex.size=sqrt(hs*50 +1), main="Hubs", vertex.label = NA)
plot(max_comp, edge.arrow.size = .1, vertex.size=sqrt(as*50 +1), main="authority", vertex.label = NA, layout = layout_with_kk)

# PageRank
pg <- page_rank(max_comp, vids = V(max_comp), directed = TRUE, damping = 0.5)
pg <- as.vector(pg[[1]])

# HITS
as <- authority.score(max_comp, scale = TRUE) %>% as.data.frame()
hs <- hub.score(max_comp, scale = TRUE) %>% as.data.frame()

# All data frame
df_cent <- data.frame(degree = degree(max_comp))
df_cent$closeness = closeness(max_comp)
df_cent$betweenness = betweenness(max_comp)
df_cent$authority = as
df_cent$hub = hs
df_cent$page_rank = pg
df_cent$row.names <- V(max_comp)$name
df_cent$cat <- V(max_comp)$cat
df_cent <- arrange(df_cent, desc(degree, closeness, betweenness, authority, hub, page_rank))

# Proximity Measure


# Frequent subgraph mining #not Working wit any graph
# gSpan
database <- array(dim = 1)
database[1] <- list(max_comp)
gSpan <- gspan(database, support = "60%") 

# SUBDUE
subdue <- subdue(max_comp)

# Graph cluster analysis
# Community structure detection based on edge betweenness
csd <- cluster_edge_betweenness(max_comp, weights = NULL, directed = TRUE,
                         edge.betweenness = TRUE, merges = TRUE, bridges = TRUE,
                         modularity = TRUE, membership = TRUE)
dendPlot(csd, mode="hclust")
plot(csd, max_comp, vertex.label = NA, vertex.size = 2.5, vertex.size = rescale(hs, 1, 8), edge.width = 0.2, edge.arrow.size = 0.2,
     edge.arrow.width = 0.2, layout = layout_with_kk)
length(csd)
modularity(csd)


# Finding communities in graphs based on statistical meachanics
smc <- cluster_spinglass(max_comp, weights = NULL, vertex = NULL, spins = 25,
                         parupdate = FALSE, start.temp = 1, stop.temp = 0.01, cool.fact = 0.99,
                         update.rule = c("config"), gamma = 1,
                         implementation = c("orig", "neg"), gamma.minus = 1) 
#Cliques
net.sym <- as.undirected(max_comp, mode= "collapse")
cliques(net.sym)
table(sapply(cliques(net.sym), length))
largest_cliques(net.sym)

# Spectral partitioning
at.lap <- graph.laplacian(max_comp)
eig.at <- eigen(at.lap)
plot(eig.at$values, col="blue",ylab="Eigenvalues of Graph Laplacian")

# Articulation points
ap <- articulation_points(max_comp)
ap1 <- V(max_comp)[ap]
V(max_comp)$color = "red"
V(max_comp)[ap]$color = "green"
plot(max_comp, edge.arrow.size = .1, vertex.color = V(max_comp)$color, vertex.size = 2.3, 
     vertex.label = ifelse(degree(max_comp) > 20, V(max_comp)$name, NA), vertex.label.cex = 0.5, edge.color = "black", layout = layout_with_kk)

# Subgroups
clique.number(max_comp)
largest.cliques(max_comp)

# k-Cores
coreness <- graph.coreness(max_comp)
table(coreness)
maxCoreness <- max(coreness)

V(max_comp)$color <- coreness
op <- par(mar = rep(0, 4))
plot(max_comp, edge.arrow.size = .1, vertex.size = 3, vertex.label.cex = 0.2)
par(op)

colors <- rainbow(maxCoreness)
op <- par(mar = rep(0, 4))
plot(max_comp, edge.arrow.size = .1, vertex.size = 3, vertex.label.cex = 0.3, 
     vertex.label = coreness, vertex.color = colors[coreness])
par(op)

coreness_net <- max_comp
V(coreness_net)$name <- coreness
V(coreness_net)$color <- colors[coreness]
coreness_net2_8 <- induced.subgraph(coreness_net, 
                                    vids = which(coreness > 1), impl = "copy_and_delete")
coreness_net3_8 <- induced.subgraph(coreness_net, 
                                    vids = which(coreness > 2), impl = "copy_and_delete")
coreness_net4_8 <- induced.subgraph(coreness_net, 
                                    vids = which(coreness > 3), impl = "copy_and_delete")
coreness_net5_8 <- induced.subgraph(coreness_net, 
                                    vids = which(coreness > 4), impl = "copy_and_delete")
coreness_net6_8 <- induced.subgraph(coreness_net, 
                                    vids = which(coreness > 5), impl = "copy_and_delete")
coreness_net7_8 <- induced.subgraph(coreness_net, 
                                    vids = which(coreness > 6), impl = "copy_and_delete")
coreness_net8_8 <- induced.subgraph(coreness_net, 
                                    vids = which(coreness > 7), impl = "copy_and_delete")
lay <- layout.kamada.kawai(coreness_net)
op <- par(mfrow = c(1,1), mar = rep(1, 4))
plot(coreness_net, edge.arrow.size = .1, 
     vertex.size = 3.5, layout = lay, vertex.label.cex = 0.4, main = "All k-cores" )
plot(coreness_net2_8, edge.arrow.size = .1, 
     vertex.size = 3.5, layout = lay, vertex.label.cex = 0.4, main = "2-8 k-cores" )
plot(coreness_net3_8, edge.arrow.size = .1, 
     vertex.size = 3.5, layout = lay, vertex.label.cex = 0.4, main = "3-8 k-cores" )
plot(coreness_net4_8, edge.arrow.size = .1, 
     vertex.size = 3.5, layout = lay, vertex.label.cex = 0.4, main = "4-8 k-cores" )
plot(coreness_net5_8, edge.arrow.size = .1, 
     vertex.size = 3.5, layout = lay, vertex.label.cex = 0.4, main = "5-8 k-cores" )
plot(coreness_net6_8, edge.arrow.size = .1, 
     vertex.size = 3.5, layout = lay, vertex.label.cex = 0.4, main = "6-8 k-cores" )
plot(coreness_net7_8, edge.arrow.size = .1, 
     vertex.size = 3.5, layout = lay, vertex.label.cex = 0.4, main = "7-8 k-cores" )
plot(coreness_net8_8, edge.arrow.size = .1, 
     vertex.size = 3.5, layout = lay, vertex.label.cex = 0.4, main = "8-8 k-cores" )
par(op)

V(coreness_net)$cor <- coreness
plot(coreness_net, layout = lay, edge.arrow.size = .1, vertex.size = 3.5, 
     vertex.color = ifelse(V(coreness_net)$cor > 1, "yellowgreen", "white"),
     vertex.label.cex = 0.4, main = "2-8 k-cores" )
plot(coreness_net, layout = lay, edge.arrow.size = .1, vertex.size = 3.5, 
     vertex.color = ifelse(V(coreness_net)$cor > 2, "yellowgreen", "white"),
     vertex.label.cex = 0.4, main = "3-8 k-cores" )
plot(coreness_net, layout = lay, edge.arrow.size = .1, vertex.size = 3.5, 
     vertex.color = ifelse(V(coreness_net)$cor > 3, "yellowgreen", "white"),
     vertex.label.cex = 0.4, main = "4-8 k-cores" )
plot(coreness_net, layout = lay, edge.arrow.size = .1, vertex.size = 3.5, 
     vertex.color = ifelse(V(coreness_net)$cor > 4, "yellowgreen", "white"),
     vertex.label.cex = 0.4, main = "5-8 k-cores" )
plot(coreness_net, layout = lay, edge.arrow.size = .1, vertex.size = 3.5, 
     vertex.color = ifelse(V(coreness_net)$cor > 5, "yellowgreen", "white"),
     vertex.label.cex = 0.4, main = "6-8 k-cores" )
plot(coreness_net, layout = lay, edge.arrow.size = .1, vertex.size = 3.5, 
     vertex.color = ifelse(V(coreness_net)$cor > 6, "yellowgreen", "white"),
     vertex.label.cex = 0.4, main = "7-8 k-cores" )
plot(coreness_net, layout = lay, edge.arrow.size = .1, vertex.size = 3.5, 
     vertex.color = ifelse(V(coreness_net)$cor > 7, "yellowgreen", "white"),
     vertex.label.cex = 0.4, main = "8-8 k-cores" )

plot(coreness_net, layout = lay, edge.arrow.size = .1, vertex.size = 3.5, 
     vertex.color = ifelse(V(coreness_net)$cor  == 1, "yellow", 
                           ifelse(V(coreness_net)$cor  == 2,"green" ,
                           ifelse(V(coreness_net)$cor  == 3,"blue" ,
                           ifelse(V(coreness_net)$cor  == 4,"red" ,
                           ifelse(V(coreness_net)$cor  == 5,"orange" ,
                           ifelse(V(coreness_net)$cor  == 6,"violet" ,
                           ifelse(V(coreness_net)$cor  == 7,"black" ,
                           ifelse(V(coreness_net)$cor  == 8,"grey" , "white")))))))),
     vertex.label.cex = 0.4, main = "all k-cores" )





# Community detection algorithms

table(V(max_comp)$code)
grp_num <- as.numeric(factor(V(max_comp)$code))
modularity(max_comp, grp_num)

