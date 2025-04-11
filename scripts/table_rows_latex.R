library(tibble)
library(dplyr)

df_analyses <- read.csv("output/empirical_munged.csv") |> 
    as_tibble() |>
    filter(inference == "empirical") |>
    filter(type == "pooled")

df_metadata <- read.csv("data/empirical/metadata.csv") |> 
    filter(skip != 1) |>
    as_tibble()

names1 <- gsub("\\.tree", "", df_metadata$Filename)
names1 <- gsub("\\.tre", "", names1)
df_metadata$name <- names1

df <- inner_join(df_metadata, df_analyses, by = "name") |>
    arrange(Clade)

for (i in 1:nrow(df)){
    cite_label <- paste0("\\cite{", df$bibtex_key[i], "}")
    taxon_sampl_cite_label <- paste0("\\cite{", df$bibtex_sampling[i], "}")
    
    cat(df$Clade[i])
    cat("  \t &  ")
    cat(df$Taxonomic.rank[i])
    cat("  \t &  ")
    cat(cite_label)
    if (df$Clade[i] == "Aves"){
      #cat("$^*$")
      cat("\footnote{We obtained the phylogeny tree file from \cite{quintero2022macroevolutionary}, however the analysis of inferring the phylogeny was conducted by \cite{Jetz2012}.}")
    }
    cat("  \t &  ")
    cat(format(df$Root.age[i], digits = 1))
    cat("  \t &  ")
    cat(format(df$NTips[i], digits = 1))
    cat("  \t &  ")
    cat(format(df$tree_netdiv[i], digits = 3, nsmall = 3))
    cat("  \t &  ")
    #cat(format(df$N_per_time[i], digits = 3, nsmall = 3))
    cat(format(df$number_of_supported[i], digits = 1, nsmall = 1))
    cat("  \t &  ")
    cat(format(df$P.extant.sampling[i], digits = 2, nsmall = 2))
    cat("  \t &  ")
    cat(taxon_sampl_cite_label)
    cat("  \t \\\\ ")

    cat("\n")
}
