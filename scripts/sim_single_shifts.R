library(TESS)

## find all active lineages at time t
## and draw a node index at random (uniform)
draw_random_node <- function(time1, phy){
  nd <- node.depth.edgelength(phy)
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

## simulate backbone trees trees
age <- 60
set.seed(123)

#tr <- tess.sim.age(n = 500, age = age, lambda = 0.3, mu = 0.21)
tr <- tess.sim.age(n = 500, age = age, lambda = 0.325, mu = 0.235)

backbone_trees <- list()
j = 1; for (i in seq_along(tr)){
  if (length(tr[[i]]$tip.label) > 99){ ## discard too small trees
    backbone_trees[[j]] = tr[[i]]
    j = j + 1
  }
}

## simulate high-diversity trees for a short time 
#high_div_lambda_trees <- tess.sim.age(n = 500, age = 15, lambda = 0.5, mu = 0.21)
#high_div_mu_trees <- tess.sim.age(n = 500, age = 15, lambda = 0.3, mu = 0.01)
high_div_lambda_trees <- tess.sim.age(n = 500, age = 15, lambda = 0.52, mu = 0.235)
high_div_mu_trees <- tess.sim.age(n = 500, age = 15, lambda = 0.325, mu = 0.0)
## simulate low-diversity trees for a longer time (so we might have a chance of recovering them)
low_div_trees <- tess.sim.age(n = 500, age = 40, lambda = 0.13, mu = 0.235)

expected_number_of_tips <- function(lambda, mu, age){
  res <- 2 * (1 + (lambda/(lambda-mu))*(exp(age*(lambda-mu))-1))
  return(res)
}
expected_number_of_tips(0.325, 0.235, 60)
expected_number_of_tips(0.52, 0.235, 15)
expected_number_of_tips(0.325, 0.0, 15)
expected_number_of_tips(0.13, 0.235, 40)


upshift_lambda_trees <- list()
upshift_mu_trees <- list()
downshift_trees <- list()

for (i in 1:350){
  upshift_lambda_trees[[i]] <- graft_tree(backbone_trees[[i]], high_div_lambda_trees[[i]], 45)
  upshift_mu_trees[[i]] <- graft_tree(backbone_trees[[i]], high_div_mu_trees[[i]], 45)
  downshift_trees[[i]] <- graft_tree(backbone_trees[[i]], low_div_trees[[i]], 20)
}

par(mfrow=c(1,2))
plot.phylo(upshift_lambda_trees[[3]], show.tip.label = FALSE)
plot(upshift_mu_trees[[3]], show.tip.label = FALSE)

ntip1 <- sapply(backbone_trees, function(x) length(x$tip.label))
ntip2 <- sapply(upshift_lambda_trees, function(x) length(x$tip.label))
ntip3 <- sapply(upshift_mu_trees, function(x) length(x$tip.label))
ntip4 <- sapply(downshift_trees, function(x) length(x$tip.label))

xmax <- max(ntip1, ntip2, ntip3, ntip4)
xmin <- min(ntip1, ntip2, ntip3, ntip4)
xlim = c(xmin, xmax)

n = 20
par(mfrow = c(3,1)); hist(ntip1, breaks = n, xlim = xlim); hist(ntip2, breaks = n, xlim = xlim); hist(ntip3, breaks = n, xlim = xlim); par(mfrow=c(1,1))

par(mfrow = c(1,2)); 
hist(ntip2 - ntip1[1:350], main = "upshift - backbone", breaks = 20);
hist(ntip3 - ntip1[1:350], main = "downshift - backbone", breaks = 20);
par(mfrow=c(1,1))



## save to file
#for (i in seq_along(upshift_trees)){
for (i in 1:350){
  fpath <- paste0("data/simulations/grafts/backbone/", i, ".tre")
  write.tree(backbone_trees[[i]], fpath) 
  
  fpath <- paste0("data/simulations/grafts/upshift_lambda/", i, ".tre")
  write.tree(upshift_lambda_trees[[i]], fpath) 
  
  fpath <- paste0("data/simulations/grafts/upshift_mu/", i, ".tre")
  write.tree(upshift_mu_trees[[i]], fpath) 
  
  fpath <- paste0("data/simulations/grafts/downshift/", i, ".tre")
  write.tree(downshift_trees[[i]], fpath) 
}


## some quick plots
fpaths <- Sys.glob("data/simulations/grafts/*/119.tre")
tr <- lapply(fpaths, read.tree)
par(mfrow=c(1,3))
plot.phylo(tr[[1]], show.tip.label = F, main = paste0("backbone (n=", length(tr[[1]]$tip.label), ")"))
plot.phylo(tr[[2]], show.tip.label = F, main = paste0("downshift (n=", length(tr[[2]]$tip.label), ")"))
plot.phylo(tr[[3]], show.tip.label = F, main = paste0("upshift (n=", length(tr[[3]]$tip.label), ")"))


## plot some more
par(mfrow = c(1,1))

hist(ntip2 - ntip3)

par(mfrow = c(2,1))
upshift_tips1 <- sapply(high_div_lambda_trees, function(x) length(x$tip.label))
upshift_tips2 <- sapply(high_div_mu_trees, function(x) length(x$tip.label))
hist(upshift_tips1, main = paste0("model w/ high lambda, mean ntip = ",mean(upshift_tips1)), breaks = 30, xlim = c(0.0, 1500))
hist(upshift_tips2, main = paste0("model w/ low lambda, mean ntip = ",mean(upshift_tips2)), breaks = 30, xlim = c(0.0, 1500))

var(upshift_tips1)
var(upshift_tips2)

mean(upshift_tips1)
mean(upshift_tips2)



hist(ntip3, main = "upshift (mu)", breaks = 30)
mean(ntip2)
mean(ntip3)

