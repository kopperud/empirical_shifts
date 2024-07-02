library(ape)

fishtree <- read.tree("data/empirical/Actinopterygii_Rabosky2018.tree")

ntip <- length(fishtree$tip.label)
internal_node_indices <- (ntip+2):max(fishtree$edge)


set.seed(123)
#nodes <- sample(internal_node_indices, 1000)
nodes <- internal_node_indices


trees <- list()
i <- 1L
used_nodes <- numeric()
for (node in nodes){
  tree <- extract.clade(fishtree, node)
  if (length(tree$tip.label) > 30){
    trees[[i]] <- tree
    used_nodes <- append(used_nodes, node)
    i <- i + 1L
  }
  if (i > 12000){
    break
  }
}


for (j in seq_along(trees)){
  fpath <- paste0("data/empirical_fish_subtrees/Actinopterygii_node_", used_nodes[j], ".tre")
  write.tree(trees[[j]], fpath)
}


