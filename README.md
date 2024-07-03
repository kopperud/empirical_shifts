## Empirical analyses

Diversification analyses for large empirical phylogenies done with [Pesto](https://github.com/kopperud/Pesto.jl)

## Contents

This repository contains several scripts that perform simulations of reconstructed phylogenetic trees, infer branch-specific diversification rates and number of rate shifts using `Pesto`, perform post-processing and produce pdf figures.

## Organization 

We have the following directories

* `data`: this is where the data (empirical and simulated trees) go, as well as the metadata table
* `figures`: this is where the pdf figures go. The contents are not committed to the git repo
* `output`: this is where the results of the analyses go
* `output/*.csv`: the `csv` files in this directory are post-processed (or munged) results that are more easily plotted than the raw result files

In the `output` directory, there are several subdirectories for each of type of content, i.e. empirical, validation study etc. Within each of these we have
* `output/*/newick/*.tre`: these are extended newick strings with metadata for each node in the tree, including estimates of branch-specific diversification rates, the **change** (i.e the delta) in diversification rates, number of rate shifts, Bayes factors for the support of rate shifts, etc.
* `output/*/rates/*.csv`: these are `csv` files that essentially have the same information as the newick strings, except that the estimates are given in a comma-separated format (for easier parsing) and the information of the phylogeny is not shown except for the node and edge indices
* `output/*/jld2/*.jld2`: these are binary files that take up a lot of space (gigabytes). These include some useful information for each analysis, including the speciation and extinction rate categories, the log-likelihood, as well as the large three-dimensional array `N`, where the first index is the branch index, and the second and third indices represent the arrival and departure indices for the rate categories.
* `output/*/shift_rate_through_time/*.csv`: these is an intermediate storage for estimates instantaneous shift rates within a phylogeny, given as a time series 

## Scripts

Most of the scripts are written in the `julia` programming language, however there are also a few R scripts and some bash files for scheduling the analyses on a high-performance computer.

The prefix for the script name will tell you what analysis the script belongs to. For example, the `empirical_01_inference.jl` file is the first script one should run in the analyses of the empirical datasets. The second script one should run is `empirical_02_shift_through_time.jl`, etc. In an overview, these are the main scripts:

* `empirical*`: these scripts are for the analyses of the empirical phylogenies
* `age_scaling*`: these scripts are for the validation study using variable tree heights
* `up_vs_down*`: these scripts are for the validation study where we assess whether and how precise upwards vs downward shifts can be inferred
* `rate_shift_type*`: these scripts are for the validation study of the rate shift type, i.e. can we infer whether it was a shift in the speciation rate, the extinction rate, or whether it was a simultaneous shift in both rates that led to the divergence in species richness
* `fishtree*`: these scripts pertain to the supplemental materials where we subsampled many smaller clades in the actinopterygii tree from Rabosky (2018), and ran the inference procedure in each subtree.
* `shifts_through_time*`: these scripts are setting up the figure for the instantaneous shift rate within phylogenies, averaged across the number of active lineages and plotted as a function of time
* `table_rows_latex.R`: this script formats and prints the contents of the data table in the supplementary 

## Dependencies

The scripts require several dependencies (either `julia` modules or `R` packages) in order to run. The dependencies are usually listed in the very top of the script. 

I recommend to first start a separate julia environment for this project, to avoid mixing dependencies with the general environment. Open the REPL in the in the root directory, and type in

```julia
import Pkg

Pkg.activate(".")
```

The two most important dependencies, for simulating trees and analyzing the rates, can be installed in `julia` as follows.

```julia
Pkg.add(url="https://github.com/kopperud/BirthDeathSimulation.jl")
Pkg.add(url="https://github.com/kopperud/Pesto.jl")
```

The other dependencies can be added from the main package registry, by typing in 
```julia
Pkg.add("Distributions")
Pkg.add("CairoMakie")
Pkg.add("CSV")
Pkg.add("ProgressMeter")
Pkg.add("LaTeXStrings")
Pkg.add("DataFrames")
Pkg.add("Glob")
Pkg.add("JLD2")
Pkg.add("Revise")
Pkg.add("RCall")
Pkg.add("Printf")
Pkg.add("Measures")
Pkg.add("Statistics")
```

Although not all of these dependencies will be necessary for any one script. The R code dependencies can be installed in the usual way with `install.packages("...")`, except for the `treeio` and `ggtree` packages which need to be installed from github using 

```R
install.packages("remotes")
remotes::install_github("YuLab-SMU/treeio")
remotes::install_github("YuLab-SMU/ggtree")
```

## How to run the scripts

The scripts can be run from the command line as `julia scripts/script_name.jl` or `Rscript scripts/script_name.R`, however the postprocessing and figure making scripts are more intended for interactive use. I recommend to run these in an interactive way, either by pasting in the lines in the REPL, or by using the `julia` extension in `visual studio code`, which behaves similarly to how `R Studio` works, with different panels for the source code, the REPL and a preview window for the figures.

The inference scripts, in particular for the validation study, can be quite resource intensive due to the number of phylogenies (up to 4000 in one set of analyses). These scripts will not be feasible to run unless you reduce the number of phylogenies, or use very powerful computer. The likelihood calculation (and other parts) in `Pesto` are multithreaded, so using a HPC with many cores will speed it up greatly. I used a HPC setup where I scheduled to use 80-100 CPU cores, which worked much faster than using a personal PC with the typical 4-8 CPU cores. Note however that the CPU cores need to be located on the same physical computer. Parallel execution across diferent machines is not supported. As the multithreading depends on the tree structure, it will only be possible to use two CPU cores when the tree traversal is close to the root. When the tree traversal is closer to the tips, where there are hundreds of branches, then many more threads can be utilized. Averaged over time this means that scheduling many (say 80 cores) will result in some cores being idle for some of the time. In other words, as the problem set is not trivially parallelizeable, the CPU utilization will not be perfectly optimized.


