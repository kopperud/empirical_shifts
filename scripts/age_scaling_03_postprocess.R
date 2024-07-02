library(rhdf5)
library(ape)
library(treeio)
library(dplyr)
library(ggtree)
library(ggplot2)
library(patchwork)
library(readr)


setwd("~/projects/pesto_empirical_clean/")


readNumberOfShifts <- function(dir = "age_scaling_effect", fpath){
  hi <- strsplit(
    strsplit(fpath, "/")[[1]][5],
    "_")[[1]] |>
    readr::parse_number()
  height <- hi[1]
  tree_index <- hi[2]
  
  
  fpath <- paste0("output/simulations/", dir, "/jld2/h", height, "_", tree_index, ".jld2")
  mu <- h5read(fpath, "mu")
  lambda <- h5read(fpath, "lambda")
  lambdaml <- h5read(fpath, "lambdaml")
  muml <- h5read(fpath, "muml")
  etaml <- h5read(fpath, "etaml")
  
  
  
  
  fpath <- paste0("output/simulations/", dir, "/rates/h", height, "_", tree_index, ".csv")
  df10 <- read.csv(fpath) |> as_tibble()
  df11 <- df10 %>%
    filter(shift_bf > 10.0)
  Ntotal <- sum(df10$nshift)
  N_supported <- df11$nshift
  
  phypath <- paste0("output/simulations/", dir, "/newick/h", height, "_", tree_index, ".tre")
  phy <- read.tree(phypath)
  tl <- sum(phy$edge.length)
  ntip <- length(phy$tip.label)
  
  df1 <- tibble(
    "tree_index" = tree_index,
    "ntip" = ntip,
    "height" = height,
    "N_total" = Ntotal,
    "treelength" = tl,
    "muml" = muml,
    "etaml" = etaml,
    "lambdaml" = lambdaml,
    "type" = "pooled",
    "inference" = dir,
    "how_many_supported" = nrow(df11)
  )
  
  df2 <- tibble(
    "tree_index" = tree_index,
    "ntip" = ntip,
    "height" = height,
    "N_total" = sum(N_supported),
    "treelength" = tl,
    "muml" = muml,
    "etaml" = etaml,
    "lambdaml" = lambdaml,
    "type" = "strong support",
    "inference" = dir,
    "how_many_supported" = nrow(df11)
  )
  
  df <- bind_rows(df1, df2)
  return(df)
}

heights <- seq(30, 100, length.out = 8)
heights



fpaths<- Sys.glob("output/simulations/age_scaling_effect/jld2/*.jld2")

dfs1 <- list()
ix <- length(fpaths)
pb <- txtProgressBar(min = 1, max = ix, initial = 1) 
for (i in 1:ix){
  fpath <- fpaths[i]
  dfs1[[i]] <- readNumberOfShifts(dir = "age_scaling_effect", fpath)
  setTxtProgressBar(pb,i)
};close(pb)
df1 <- bind_rows(dfs1)


df <- df1
df[["N_per_time"]] <- df$N_total / df$treelength
df[["support_per_time"]] <- df$how_many_supported / df$treelength

write.csv(df, "output/age_scaling_effect_munged.csv")












  




