library(ggtree)
library(treeio)
library(ggplot2)
library(patchwork)
library(deeptime)
library(tidytree)


setwd("~/projects/empirical_shifts/")

################################################################
##
##  read tree files and calculate supported nodes
##
################################################################


tree1 <- read.beast.newick("output/empirical/newick/Mammalia_AlvarezCarretero2022.tre")
tree2 <- read.beast.newick("output/empirical/newick/Anura_Portik2023.tre")

mammal_supported_upshifts <- as.numeric(tree1@data[tree1@data$shift_bf > 10 & tree1@data$delta_netdiv > 0,][["node"]])
mammal_supported_downshifts <- as.numeric(tree1@data[tree1@data$shift_bf > 10 & tree1@data$delta_netdiv < 0,][["node"]])
#mammal_supported_nodes <- as.numeric(tree1@data[tree1@data$shift_bf > 10,][["node"]])
frog_supported_upshifts <- as.numeric(tree2@data[tree2@data$shift_bf > 10 & tree2@data$delta_netdiv > 0,][["node"]])
frog_supported_downshifts <- as.numeric(tree2@data[tree2@data$shift_bf > 10 & tree2@data$delta_netdiv < 0,][["node"]])


################################################################
##
## find the taxonomic groups that are interesting
## in terms of understanding the clades, and some highlights
## where the clades contain one or more rate shifts
##
################################################################

mammal_clade_highlights <- read.table("data/mammal_clade_labels.tsv", sep = "\t", header = TRUE)
mammal_clade_highlights$fake_label <- ""
#mammal_clade_highlights <-
  #mammal_clade_highlights %>%
  #filter(Skip == 0)


frog_clade_highlights <- read.table("data/frog_clade_labels.tsv", sep = "\t", header = TRUE)
frog_clade_highlights$fake_label <- ""
#frog_clade_highlights <-
#  frog_clade_highlights %>%
#  filter(Skip == 0)

################################################################
##
## plot the main figure 1, mammals and frogs together
##
################################################################

#dot_color <- "#f56b5b"
dot_color <- "#fc9272"

offset <- 15
offset.text <- 35
p1 <- revts(ggtree(tree1, aes(color = mean_netdiv), linewidth = 0.2)) +
  coord_radial(theta = "y", start = -0.5 * pi, end = 1.35 * pi) +
  scale_x_continuous(breaks = seq(-60, 0, 20), labels = abs(seq(-60, 0, 20)),
                     expand = expansion(mult = c(0.05, 0)),
                     guide = guide_axis_stack(guide_geo("periods", neg = TRUE,
                                                        rot = -90, size = "auto",
                                                        height = unit(1, "line")),
                                              guide_axis(),
                                              spacing = unit(0, "line"))) +
  scale_y_continuous(guide = NULL, expand = expansion(mult = c(0.01, 0.01))) + 
  geom_point2(aes(subset=(node %in% mammal_supported_upshifts)), size = 1.5, color = "black", fill = dot_color, shape = 21, stroke = 0.5) +
  geom_point2(aes(subset=(node %in% mammal_supported_downshifts)), size = 1.5, color = "black", fill = "lightgray", shape = 21, stroke = 0.5) +
  ggtitle("Mammalia") +
  scale_colour_gradient(
    name = "", 
    #low = "#132B43",
    #high = "#3496eb",
    low = "#08306b",
    high = "#c6dbef",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  ) +
  #geom_cladelab(data = mammal_clade_highlights, mapping = aes(node = NodeNumber, label = CladeName), offset = 15, offset.text = 35, align = FALSE, angle = 0)
  geom_cladelab(data = mammal_clade_highlights, mapping = aes(node = NodeNumber, label = fake_label), offset = 15, offset.text = 35, align = FALSE, angle = 0)
p1

p2 <- revts(ggtree(tree2, aes(color = mean_netdiv), layout = "circular", linewidth = 0.2)) +
  coord_radial(theta = "y", start = -0.5 * pi, end = 1.35 * pi) +
  scale_x_continuous(breaks = seq(-60, 0, 20), labels = abs(seq(-60, 0, 20)),
                     expand = expansion(mult = c(0.05, 0)),
                     guide = guide_axis_stack(guide_geo("periods", neg = TRUE,
                                                        rot = -90, size = "auto",
                                                        height = unit(1, "line")),
                                              guide_axis(),
                                              spacing = unit(0, "line"))) +
  scale_y_continuous(guide = NULL, expand = expansion(mult = c(0.01, 0.01))) + 
  geom_point2(aes(subset=(node %in% frog_supported_upshifts)), size = 1.5, color = "black", fill = dot_color, shape = 21, stroke = 0.5) +
  geom_point2(aes(subset=(node %in% frog_supported_downshifts)), size = 1.5, color = "black", fill = "lightgray", shape = 21, stroke = 0.5) +
  ggtitle("Anura") +
  scale_colour_gradient(
    name = "", 
    #low = "#132B43",
    #high = "#4a990e",
    low = "#00441b",
    high = "#c7e9c0",
    breaks = c(0.05, 0.10, 0.15),
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  ) +
  #geom_cladelab(data = frog_clade_highlights, mapping = aes(node = NodeNumber, label = CladeName), offset = 15, offset.text = 35, align = FALSE, angle = 0)
  geom_cladelab(data = frog_clade_highlights, mapping = aes(node = NodeNumber, label = fake_label), offset = 15, offset.text = 35, align = FALSE, angle = 0)
p2

p <- p1 + p2 + plot_layout(ncol = 2, guides = "collect") &
  theme(plot.title = element_text(hjust = 0.5, size = 28),
        legend.position = "bottom",
        theme(plot.margin = unit(c(0,0,0,0), "mm")),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        ) 

p
ggsave("figures/two_circular_trees.pdf", p, width = 300, height = 150, units = "mm")

################################################################
##
## plot very large versions with node labels and tip labels
##
################################################################
p1a <- ggtree(tree1, aes(color = mean_netdiv), layout = "circular") +
  geom_tiplab(size = 1) +
  geom_nodelab(size = 1, aes(filter=shift_bf > 10, label = node)) +
  geom_point2(aes(subset=(node %in% mammal_supported_nodes)), size = 3, color = "red")
ggsave("figures/mammal_inspect.png", p1a, width = 1400, height = 1400, units = "mm", limitsize = FALSE)

p2a <- ggtree(tree2, aes(color = mean_netdiv), layout = "circular") +
  geom_tiplab(size = 1) +
  geom_nodelab(size = 1, aes(filter=shift_bf > 10, label = node)) +
  geom_point2(aes(subset=(node %in% frog_supported_nodes)), size = 3, color = "red")
ggsave("figures/frog_inspect.png", p2a, width = 1400, height = 1400, units = "mm", limitsize = FALSE)





################################################################
##
## plot all of the phylogenies in the dataset (quite slow to read files)
##
################################################################
fpaths <- Sys.glob("output/empirical/newick/*.tre")

for (i in seq_along(fpaths)){
  cat(".")
  
  tree <- read.beast.newick(fpaths[i])
  title_name <- strsplit(fpaths[i], "/")[[1]][4]
  
  p <- ggtree(tree, aes(color = mean_netdiv), layout = "circular") +
    scale_colour_gradient(
      name = "net-div", 
      low = "#132B43",
      high = "green",
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "colour"
    ) +
    ggtitle(title_name)
  ggsave(paste0("figures/empirical/", title_name, ".pdf"), p, width = 400, height = 400, units = "mm")
}; cat("\n")

