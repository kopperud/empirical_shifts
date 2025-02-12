library(treeio)
library(ggtree)
library(ggplot2)

setwd("~/projects/empirical_shifts/")

fpaths <- Sys.glob("output/simulations/grafts/upshift_lambda/newick/*.tre")[1:10]

scalexformat <- function(x) sprintf("%.0f Ma", abs(th - round(x, 0)))
th <- 60

trees <- lapply(fpaths, read.beast.newick)

library(patchwork)

mu_low <- min(sapply(trees, function(x) min(x@data$mean_mu)))
mu_high <- max(sapply(trees, function(x) max(x@data$mean_mu)))

lambda_low <- min(sapply(trees, function(x) min(x@data$mean_lambda)))
lambda_high <- max(sapply(trees, function(x) max(x@data$mean_lambda)))

netdiv_low <- min(sapply(trees, function(x) min(x@data$mean_netdiv)))
netdiv_high <- max(sapply(trees, function(x) max(x@data$mean_netdiv)))

print(paste0("range of speciation rates: ", lambda_low, ", ", lambda_high))
print(paste0("range of exinction rates: ", mu_low, ", ", mu_high))
print(paste0("range of netdiv rates: ", netdiv_low, ", ", netdiv_high))

ps_left <- list()
ps_middle <- list()
ps_right <- list()
for (i in seq_along(fpaths)){
  ps_left[[i]] <- ggtree(trees[[i]], aes(col = mean_lambda)) +
  scale_colour_gradient(
    low = "black",
    high = "cyan",
    limits = c(lambda_low, lambda_high),
    name = "speciation") +
   theme(
     legend.key.size = unit(15, "mm"),
     legend.title = element_text(size=20),
     legend.text = element_text(size=15),
     plot.title = element_text(size=22, hjust = 0.5),
     axis.title = element_text(size=20),
         ) +
    labs(y = paste0("tree ", i))
  
 ps_middle[[i]] <- ggtree(trees[[i]], aes(col = mean_mu)) +
  scale_colour_gradient(
    low = "black",
    high = "red",
    limits = c(mu_low, mu_high),
    name = "extinction") + 
   theme(
     legend.key.size = unit(15, "mm"),
     legend.title = element_text(size=20),
     legend.text = element_text(size=15),
     plot.title = element_text(size=22, hjust = 0.5),
         )
 
 ps_right[[i]] <- ggtree(trees[[i]], aes(col = mean_netdiv)) +
  scale_colour_gradient(
    low = "black",
    high = "green",
    limits = c(netdiv_low, netdiv_high),
    name = "netdiv") +
   theme(
     legend.key.size = unit(15, "mm"),
     legend.title = element_text(size=20),
     legend.text = element_text(size=15),
     plot.title = element_text(size=22, hjust = 0.5),
         ) 
}
ps_left[[1]] <- ps_left[[1]] + ggtitle("true speciation = [0.325, 0.52]")
ps_middle[[1]] <- ps_middle[[1]] + ggtitle("true extinction = 0.235")
ps_right[[1]] <- ps_right[[1]] + ggtitle("true netdiv = [0.09, 0.285]")

left <- Reduce("+", ps_left) + plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")
middle <- Reduce("+", ps_middle) + plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")
right <- Reduce("+", ps_right) + plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

p <- left | middle | right

ggsave("figures/ten_trees_upshift_lambda.pdf", p, width = 450, height = 600*1.0, units = "mm", limitsize=FALSE)


## upshift mu

