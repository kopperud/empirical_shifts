library(TESS)

## find all active lineages at time t
draw_random_node <- function(time1, phy){
  nd <- node.depth.edgelength(phy)
  #nd <- max(nd) - nd
  is_active <- apply(phy$edge, 1, function(x) (nd[x[1]] < time1) && (nd[x[2]] > time1))
  active_branch_indices <- which(is_active)

  edge_index <- sample(active_branch_indices, 1)
  node_index <- phy$edge[edge_index,1]
  return(node_index)
}

## graft a young tree onto an older tree
graft_tree <- function(old_tree, young_tree, break_time){
  node_index <- draw_random_node(break_time, old_tree)
 
  nd <- node.depth.edgelength(old_tree)
  old_height <- max(nd)
  
  ## add some extra tip label
  young_tree$tip.label <- paste0(young_tree$tip.label, "x")
  
  ## ensure that the grafted tree stays ultrametric
  time_node_to_present <- old_height - nd[node_index] 
  young_height <- max(node.depth.edgelength(young_tree))
  delta_time <- time_node_to_present - young_height
  young_tree$root.edge <- delta_time
  
  hybrid_tree <- bind.tree(old_tree, young_tree, where = node_index)
   
  ## find the polytomy
  tbl <- table(hybrid_tree$edge[,1])
  polytomy_node <- as.numeric(names(tbl)[which.max(tbl)])
  
  ## find the tips in one of the subtrees descending from node_index
  left_branch <- which(hybrid_tree$edge[,1] == polytomy_node)[2]
  left_node <- hybrid_tree$edge[left_branch,2]
  
  if (left_node > length(hybrid_tree$tip.label)){
    tips_to_drop <- extract.clade(hybrid_tree, left_node)$tip.label
  }else{
    tips_to_drop <- hybrid_tree$tip.label[[left_node]]
  }
 
  ## drop those tips from the old tree
  if (any(grepl("x", tips_to_drop))){
    stop("error")  
  }
  hybrid_tree <- drop.tip(hybrid_tree, tips_to_drop, collapse.singles = TRUE) 

  return(hybrid_tree)
}

## up-shift trees
age <- 60
set.seed(123)

tr <- tess.sim.age(n = 500, age = age, lambda = 0.3, mu = 0.21)

backbone_trees <- list()
j = 1; for (i in seq_along(tr)){
  if (length(tr[[i]]$tip.label) > 99){
    backbone_trees[[j]] = tr[[i]]
    j = j + 1
  }
}

## trees to be grafted on
break_time <- 2*age/3

high_div_trees <- tess.sim.age(n = 500, age = 15, lambda = 0.5, mu = 0.21)
low_div_trees <- tess.sim.age(n = 500, age = 40, lambda = 0.1, mu = 0.21)


upshift_trees <- list()
downshift_trees <- list()
for (i in 1:350){
  print(i)
  upshift_trees[[i]] <- graft_tree(backbone_trees[[i]], high_div_trees[[i]], 45)
  downshift_trees[[i]] <- graft_tree(backbone_trees[[i]], low_div_trees[[i]], 20)
}

plot(upshift_trees[[3]])
plot(downshift_trees[[21]])

ntip1 <- sapply(backbone_trees, function(x) length(x$tip.label))
ntip2 <- sapply(upshift_trees, function(x) length(x$tip.label))
ntip3 <- sapply(downshift_trees, function(x) length(x$tip.label))

xmax <- max(ntip1, ntip2, ntip3)
xmin <- min(ntip1, ntip2, ntip3)
xlim = c(xmin, xmax)

n = 20
par(mfrow = c(3,1)); hist(ntip1, breaks = n, xlim = xlim); hist(ntip2, breaks = n, xlim = xlim); hist(ntip3, breaks = n, xlim = xlim); par(mfrow=c(1,1))

par(mfrow = c(1,2)); 
hist(ntip2 - ntip1[1:350], main = "upshift - backbone", breaks = 20);
hist(ntip3 - ntip1[1:350], main = "downshift - backbone", breaks = 20);
par(mfrow=c(1,1))



## save to file
for (i in seq_along(upshift_trees)){
  fpath <- paste0("data/simulations/single_shift_grafts/backbone/", i, ".tre")
  write.tree(backbone_trees[[i]], fpath) 
  
  fpath <- paste0("data/simulations/single_shift_grafts/upshift/", i, ".tre")
  write.tree(upshift_trees[[i]], fpath) 
  
  fpath <- paste0("data/simulations/single_shift_grafts/downshift/", i, ".tre")
  write.tree(downshift_trees[[i]], fpath) 
}


## some quick plots
fpaths <- Sys.glob("data/simulations/single_shift_grafts/*/119.tre")
tr <- lapply(fpaths, read.tree)
par(mfrow=c(1,3))
plot.phylo(tr[[1]], show.tip.label = F, main = paste0("backbone (n=", length(tr[[1]]$tip.label), ")"))
plot.phylo(tr[[2]], show.tip.label = F, main = paste0("downshift (n=", length(tr[[2]]$tip.label), ")"))
plot.phylo(tr[[3]], show.tip.label = F, main = paste0("upshift (n=", length(tr[[3]]$tip.label), ")"))
