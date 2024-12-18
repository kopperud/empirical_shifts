using Distributions
using Glob, DataFrames, CSV, RCall, ProgressMeter
using LaTeXStrings, Measures
using JLD2
using Printf
using CairoMakie
using Pesto

########################
##
##   some helper functions
##
##########################

function ols_regression(x, y)
    X = hcat([1 for _ in 1:length(x)], x)
    n, p = size(X)

    ## OLS
    β = (X' * X) \ X' * y
    yhat = X * β
    sigma_squared = (1 / (n - p - 1)) * sum((y .- yhat).^2) ## MLE for sigma^2
    #s_squared = (y .- yhat)' * (y .- yhat) ./ (n - p) ## OLS for sigma^2
    Varβ = inv(X' * X) .* sigma_squared
    yVar = x -> Varβ[1,1] + 2*x*Varβ[1,2] + (x^2)*Varβ[2,2]
    ySE = x -> sqrt.(yVar(x))
    return(β, Varβ, ySE)
end

function magnitude(model, N)
    r = model.λ .- model.μ
    Δr = N .* (r .- r')
    mean_magnitude = sum(Δr) / sum(N)
    return(mean_magnitude)
end

########################
##
##         read files
##
##########################


inference = "age_scaling_effect"

df = CSV.read("output/age_scaling_effect_munged.csv", DataFrame)
#df = df[df[!,:inference] .== inference,:]

########################
##
##           N tensor
## Dimensions:
## 1: edge index (M)
## 2: row in the N matrix (index i), arrival state
## 3: column in the N matrix (index j), departure state
##
##########################


d = Dict()
#names = unique(df[!,:name])
fpaths = Glob.glob("output/simulations/age_scaling_effect/jld2/*.jld2")
@showprogress for fpath in fpaths
    fname = split(Base.basename(fpath), ".")[1]
    height, i = split(fname, "_")
    
    height = parse(Int64, height[2:end])
    i = parse(Int64, i)

    x = JLD2.load(fpath)
    d[height, i] = x
end

models = Dict()
heights = [30, 40, 50, 60, 70, 80, 90, 100]

for fpath in fpaths
    fname = split(Base.basename(fpath), ".")[1]
    height, i = split(fname, "_")

    height = parse(Int64, height[2:end])
    i = parse(Int64, i)

    λ = d[height, i]["lambda"]
    μ = d[height, i]["mu"]
    η = d[height, i]["etaml"]
    models[height, i] = SSEconstant(λ, μ, η)
end


datasets = Dict()
treepaths = Glob.glob("output/simulations/age_scaling_effect/newick/h*.tre")

@showprogress for fpath in treepaths
    tree = readtree(fpath)
    fname = split(Base.basename(fpath), ".")[1]
    height, i = split(fname, "_")

    height = parse(Int64, height[2:end])
    i = parse(Int64, i)

    datasets[height,i] = SSEdata(tree, 1.0) ## complete taxon sampling
end


rates = Dict()
ratepaths = Glob.glob("output/simulations/age_scaling_effect/rates/h*.csv")
@showprogress for fpath in ratepaths
    #tree = readtree(fpath)
    fname = split(Base.basename(fpath), ".")[1]
    height, i = split(fname, "_")

    height = parse(Int64, height[2:end])
    i = parse(Int64, i)

    rates_df = CSV.read(fpath, DataFrame)
    sort!(rates_df, :edge)
    rates[height, i] = rates_df[2:end,:]
end

##################
##
##     calculate magnitudes
##
##################

n_iters = 500
iters = 1:n_iters
magnitudes_matrix = zeros(length(heights), n_iters)
magnitudes_matrix[:,:] .= NaN
ntips_matrix = similar(magnitudes_matrix) 
ntips_matrix[:,:] .= NaN

@showprogress for (height_index, height) in enumerate(heights)
    for i in iters

        if (height, i) in keys(rates)
            model = models[height,i]
            
            Nsum = sum(d[height,i]["N"], dims = 1)[1,:,:]
            m = magnitude(model, Nsum)
            magnitudes_matrix[height_index, i] = m

            ntips_matrix[height_index, i] = length(datasets[height, i].tiplab)

        end
    end 
end

plotdf = DataFrame(
    "magnitudes" => vcat(magnitudes_matrix[:,:]...),
    "heights" => repeat(heights, 500),
    "replicate" => vcat([repeat([i], 8) for i in 1:500]...),
    "ntip" => vcat(ntips_matrix...),
    #"name" => collect(keys(d))
)
filter!(:magnitudes => x -> !isnan(x), plotdf)

##################
##
##   set up the makie figure
##
###############


fig2 = Figure(size = (500, 300), fontsize = 14, 
                figure_padding = (5,8,1,1))

magnitudes = plotdf[!,:magnitudes]
heights_vec = plotdf[!,:heights]

xt = collect(range(extrema(magnitudes)...; length = 5))
xtl = [@sprintf("%.2f", x) for x in xt]

ax_hist = Axis(fig2[1,1], 
            ylabel = L"\text{frequency}", 
            xlabel = L"\text{magnitude }(\Delta r)",
            xgridvisible = false, 
            ygridvisible = false,
            xticks = (xt, xtl),
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelsize = 9,
            yticklabelsize = 9)

xt = collect(range(extrema(magnitudes)...; length = 5))
xtl = [@sprintf("%.2f", x) for x in xt]


CairoMakie.hist!(ax_hist, magnitudes, bins = 10, color = "gray")


xt = 10 .^ (collect(range(extrema(log10.(plotdf[!,:heights]))...; length = 5)))
xtl = [@sprintf("%.1f", x) for x in xt]


yr = collect(range(extrema(plotdf[!,:magnitudes])...; length = 5))
yt = yr
ytl = [@sprintf("%.2f", y) for y in yt]

ax_scatter = Axis(fig2[2,1], 
        ylabel = L"\text{magnitude } (\Delta r)", 
        xgridvisible = false, 
        ygridvisible = false,
        xscale = log10,
        xticks = (xt, xtl),
        yticks = (yt, ytl),
        xlabel = L"\text{tree height (Ma)}",
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelsize = 9,
        yticklabelsize = 9)


#for (magnitudes, heights, ax) in zip([magnitudes, magnitudes_support], [heights, heights_support], [ax_scatter_pooled, ax_scatter_support])

    β, Varβ, ySE = ols_regression(log10.(heights_vec), magnitudes)

    x = collect(Pesto.lrange(30.0, 100.0, 20))
    linefit(x) = β[1] + β[2]*log10.(x)

    yarithmetic = linefit.(x)

    yVar = Varβ[1,1] .+ 2 .* log10.(x) .* Varβ[1,2] .+ (log10.(x) .^ 2) .* Varβ[2,2]

    yupper = yarithmetic .+ 2*sqrt.(yVar)
    ylower = yarithmetic .- 2*sqrt.(yVar)
    y = linefit.(x)

    CairoMakie.band!(ax_scatter, x, ylower, yupper, color = "#e0e0e0")
    CairoMakie.lines!(ax_scatter, x, y; label = "OLS", markersize = 7, color = "gray", linestyle = :dash)

    CairoMakie.scatter!(ax_scatter, 
                        heights_vec, 
                        magnitudes,
                        color = "black",
                        markersize = 7)
    
#end

colgap!(fig2.layout, 5)
rowgap!(fig2.layout, 5)
fig2

##CairoMakie.save("figures/magnitude_simulated_trees.pdf", fig2)

ntips = plotdf[!,:ntip]

ax = Axis(fig2[1:2,2], 
    xlabel = L"\text{number of tips}", 
    ylabel = L"\text{magnitude }(\Delta r)",
    xgridvisible = false, 
    ygridvisible = false,
    xscale = log10,
    topspinevisible = false,
    rightspinevisible = false)
CairoMakie.scatter!(ax, ntips, magnitudes, color = "black")
CairoMakie.lines!(ax, [2, 30_000], [0.0, 0.0], color = "red", linestyle = :dash)


fig2

CairoMakie.save("figures/magnitude_simulated_trees.pdf", fig2)


fig3
ntips

## estimation error
## calculate true magnitudes

simulated_tree_paths = Glob.glob("data/simulations/age_scaling_effect/*.tre")
@rput simulated_tree_paths
R""" ## this is kind of slow and not necessary, I didn't save the N matrices before
library(treeio)

items <- list()
ntrees <- length(simulated_tree_paths)
tree_indices <- list()
tree_heights_R <- list()

for (i in 1:ntrees){
    bn <- strsplit(basename(simulated_tree_paths[i]), "\\.")[[1]][[1]]
    height <- as.numeric(gsub("h", "", strsplit(bn, "_")[[1]][1]))
    tree_index <- as.numeric(strsplit(bn, "_")[[1]][2])

    tr <- read.beast.newick(simulated_tree_paths[i])
    l <- lapply(tr@data$N, function(x) matrix(x, 3, 3))
    N <- Reduce("+", l)
    items[[i]] <- N
    tree_indices[[i]] <- tree_index
    tree_heights_R[[i]] <- height
}
"""
@rget items
@rget tree_indices
@rget tree_heights_R

tree_indices = Int64.(tree_indices)
tree_heights_R = Int64.(tree_heights_R)

N_true = Dict()
for i in eachindex(items)
    N_true[tree_heights_R[i], tree_indices[i]] = items[i]
end




r0 = 0.04
r = [r0, 0.07, 0.10]
ϵ = 2/3
λ = r ./ (1 - ϵ)
μ = λ .- r

Δr = r .- r'
n_heights = length(heights)
Ntrue = zeros(Int64, n_heights, n_iters, 3, 3)
magnitudes_true = zeros(Float64, n_heights, n_iters)
for i in 1:4000
    height = tree_heights_R[i]
    height_index = argmax(height .== heights)

    tree_index = tree_indices[i]
    Ntrue[height_index, tree_index,:,:] .= items[i]
    mag = sum(Ntrue[height_index, tree_index,:,:] .* Δr) / sum(Ntrue[height_index, tree_index,:,:])
    magnitudes_true[height_index, tree_index] = mag
end


magnitudes_estimated = magnitudes[:,:,1]
error = magnitudes_estimated .- magnitudes_true

function rm_na(x)
    res = x[.!isnan.(x)]
    return(res)
end


fig6 = Figure()

ax1 = Axis(fig6[1,1],
    xticks = heights,
    #xlabel = L"\text{tree height (Ma)}",
    ylabel = L"\text{magnitude}",
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false)

#xt = collect(range(30,100; length = 8))
ax3 = Axis(fig6[2,1],
    xticks = heights,
    xlabel = L"\text{tree height (Ma)}",
    ylabel = L"\text{error (}\text{mag}_\text{estimated} - \text{mag}_\text{true})",
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false)


for (i, height) in enumerate(heights)
    these_true = rm_na(magnitudes_true[i,:])
    these_estimated = rm_na(magnitudes_estimated[i,:])
    
    hist!(ax1, these_true, offset = height, scale_to = 8.0, direction = :x, bins = 25, color = (:black, 0.6), label = "true")
    hist!(ax1, these_estimated, offset = height, scale_to = 8.0, direction = :x, bins = 25, color = (:orange, 0.6), label = "estimated")

    these_errors = rm_na(error[i,:])
    hist!(ax3, these_errors, offset = height, 
            scale_to = 8.0, direction = :x, bins = 25, color = (:black, 0.6), label = "simulated trees")
end
lines!(ax3, [25, 110], [0.0, 0.0], linestyle = :dash, color = :red, label = "zero error")
axislegend(ax1, position = :lt, unique = true)
axislegend(ax3, position = :lt, unique = true)

fig6
[length(rm_na(error[i,:])) for i in 1:8]

save("figures/magnitude_vs_age_estimation_error.pdf", fig6)










