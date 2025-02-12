
setwd("/home/bkopper/projects/empirical_shifts/")

library(ggtree)
library(ape)
library(treeio)
library(dplyr)

#tr <- read.beast.newick("output/empirical/newick/Actinopterygii_Rabosky2018.tre")
#tr <- read.beast.newick("output/empirical/newick/Mammalia_AlvarezCarretero2022.tre")
tr <- read.beast.newick("/home/bkopper/fishtree.tre")
tr <- read.beast.newick("/home/bkopper/amphibians.tre")
tr <- read.beast.newick("/home/bkopper/mammals.tre")

#is_signif <- tr@data$delta_netdiv > 0.03
is_signif <- tr@data$shift_bf > 10


node_indices <- as.numeric(tr@data[is_signif,][["node"]])
treeheight <- max(node.depth.edgelength(tr@phylo))
node_time <- treeheight - node.depth.edgelength(tr@phylo)[node_indices]

ggtree(tr, aes(color = mean_netdiv))
hist(node_time, breaks = 10)
plot(node_time, df2$delta_netdiv, xlim = c(0.0, treeheight))

df2 <- tr@data %>%
  dplyr::filter(shift_bf > 10)

df2

df2

get_tip_labels <- function(phy, node_index){
  if (node_index  <= length(phy$tip.label)){
    return(phy$tip.label[node_index])
  }else{
    descendant_edges <- which(phy$edge[,1] == node_index)
    descendant_nodes <- phy$edge[descendant_edges,2]
    
    r <- lapply(descendant_nodes, function(n) get_tip_labels(phy, n))
    tl <- do.call(c, r) 
  }
  return(tl)
}

#signif_nodes#

#tips <- get_tip_labels(tr@phylo, 8019)
#plot(keep.tip(tr@phylo, tips))



tips


df2 <- tr@data %>%
  dplyr::filter(shift_bf > 100)

df2



