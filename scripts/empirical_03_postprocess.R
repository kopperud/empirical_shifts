library(rhdf5)
library(ape)
library(treeio)
library(dplyr)
library(ggtree)
library(ggplot2)
library(patchwork)
library(readr)


setwd("~/projects/empirical_shifts")


readNumberOfShifts <- function(name, subdir = "empirical"){
  fpath <- paste0("output/", subdir, "/jld2/", name, ".jld2")
  mu <- h5read(fpath, "mu")
  lambda <- h5read(fpath, "lambda")
  lambdaml <- h5read(fpath, "lambdaml")
  muml <- h5read(fpath, "muml")
  etaml <- h5read(fpath, "etaml")
  Ns <- h5read(fpath, "Nsum")
  N <- h5read(fpath, "N")
  
  
  Ntotal <- sum(N)
  
  fpath <- paste0("output/", subdir, "/rates/", name, ".csv")
  df10 <- read.csv(fpath) |> as_tibble()
  df11 <- df10 %>%
    #filter(shift_bf > 10, nshift > 0.5)
    filter(shift_bf > 100)
  number_of_supported <- nrow(df11)
  N_supported <- N[df11$edge,,,drop=FALSE]
  
  phypath <- paste0("output/", subdir, "/newick/", name, ".tre")
  tree <- read.beast.newick(phypath)
  edge_df <- tree@data
  edge_df <- arrange(edge_df, edge)
  
  bls <- tree@phylo$edge.length
  tree_netdiv <- sum(edge_df$mean_netdiv * bls) / sum(bls)
  tree_lambda <- sum(edge_df$mean_lambda * bls) / sum(bls)
  tree_mu <- sum(edge_df$mean_mu * bls) / sum(bls)
  tree_relext <- sum(edge_df$mean_relext * bls) / sum(bls)
  
  phy <- tree@phylo
  height <- max(node.depth.edgelength(phy))
  tl <- sum(phy$edge.length)
  
  df1 <- tibble(
    "name" = name,
    "height" = height,
    "N_total" = Ntotal,
    "treelength" = tl,
    "muml" = muml,
    "etaml" = etaml,
    "lambdaml" = lambdaml,
    "tree_netdiv" = tree_netdiv,
    "tree_lambda" = tree_lambda,
    "tree_mu" = tree_mu,
    "tree_relext" = tree_relext,
    "inference" = subdir,
    "number_of_supported" = number_of_supported,
    "type" = "pooled",
  )
  
  df2 <- tibble(
    "name" = name,
    "height" = height,
    "N_total" = sum(N_supported),
    "treelength" = tl,
    "muml" = muml,
    "etaml" = etaml,
    "lambdaml" = lambdaml,
    "tree_netdiv" = tree_netdiv,
    "tree_lambda" = tree_lambda,
    "tree_mu" = tree_mu,
    "tree_relext" = tree_relext,
    "inference" = subdir,
    "number_of_supported" = number_of_supported,
    "type" = "strong support",
  )
  
  df <- bind_rows(df1, df2)
  return(df)
}


fpaths <- Sys.glob("output/empirical/jld2/*.jld2")

dfs <- list()
ix <- length(fpaths)
pb <- txtProgressBar(min = 1, max = ix, initial = 1) 
for (i in 1:ix){
  fpath <- fpaths[i]
  name <- strsplit(basename(fpath), "\\.")[[1]][[1]]
  dfs[[i]] <- readNumberOfShifts(name, "empirical")
  setTxtProgressBar(pb,i)
};close(pb)

df <- bind_rows(dfs)

###############
##
## write to file
##
#############

df[["N_per_time"]] <- df$N_total / df$treelength
write.csv(df, "output/empirical_munged.csv")


