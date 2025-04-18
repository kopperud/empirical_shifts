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

function compute_ratios(name, model, rates)
    N = d[name]["N"]

    is_supported = rates[!,:is_signif]
    N = N[is_supported,:,:]

    Nmatrix = sum(N, dims = 1)[1,:,:]
    nbins = 14

    filters = ["extinction", "speciation", ""]
    #limits = [-1.2, 1.2]
    r = model.λ .- model.μ
    Δr = r .- r'
    netdiv_extrema = extrema(Δr)

    dfs = []
    for (i, filter) in enumerate(filters)
        mids, bins = makebins(Nmatrix, model, netdiv_extrema...; filter = filter, nbins = nbins)
        df = DataFrame(
            "Δr" => bins[:,3],
            "mids" => mids,
            "subset" => [i for _ in 1:nbins]
        )
        append!(dfs, [df])
    end
    df = DataFrame(
        "Δr" => dfs[3][!, "Δr"] .- dfs[1][!, "Δr"] .- dfs[2][!, "Δr"],
        "mids" => dfs[3][!, "mids"],
        "subset" => 3
    )
    dfs[3] = df

    Nλ = sum(dfs[1][!,:Δr])
    Nμ = sum(dfs[2][!,:Δr])
    Njoint = sum(dfs[3][!,:Δr])
    Nall = sum([Nλ, Nμ, Njoint])
    rs = [Nλ, Nμ, Njoint] ./ Nall
    return(rs)
end

inference = "empirical"

df = CSV.read("output/empirical_munged.csv", DataFrame)
df = df[df[!,:inference] .== inference,:]
df = df[df[!,:type] .== "strong support",:]

## N tensor
#
# Dimensions:
# 1: edge index (M)
# 2: row in the N matrix (index i), arrival state
# 3: column in the N matrix (index j), departure state
d = Dict()
names = unique(df[!,:name])
fpaths = Glob.glob("output/" * inference * "/jld2/*.jld2")
#fpaths = Glob.glob("output/empirical/jld2/*.jld2")
@showprogress for fpath in fpaths
    name = split(Base.basename(fpath), ".")[1]
    x = JLD2.load(fpath)
    d[name] = x
end

models = Dict{String, SSEconstant}()
for name in names
    λ = d[name]["lambda"]
    μ = d[name]["mu"]
    η = d[name]["etaml"]
    models[name] = SSEconstant(λ, μ, η)
end



#fpaths = glob("*.tree", "data/empirical")
fpaths = ["data/empirical/" * name * ".tree" for name in names]

df = CSV.read("data/empirical/metadata.csv", DataFrame)
ρs = Dict()
for row in eachrow(df)
    fn = row["Filename"]
    ρ = row["P extant sampling"]
    ρs[fn] = ρ
end

trees = Dict()
datasets = Dict()
for fpath in fpaths
    println(fpath)
    tree = readtree(fpath)
    
    if all(tree.edge_length .> 0)
        bn = Base.Filesystem.basename(fpath)
        trees[bn] = tree
        datasets[bn] = SSEdata(trees[bn], ρs[bn])
    end
end

rates_fpaths = Glob.glob("output/" * inference * "/rates/*.csv")
d_rates = Dict{String, DataFrame}()
@showprogress for fpath in rates_fpaths
    name = split(Base.basename(fpath), ".")[1]
    df = CSV.read(fpath, DataFrame)
    sort!(df, :edge)
    filter!(:edge => x -> x != 0, df) ## remove the root edge
    is_signif = df[!,:shift_bf] .> 100
    df[!,:is_signif] = is_signif

    d_rates[name] = df
end



#heights = [maximum(d.node_depth) for (key, d) in datasets]




####################################
##
##   Fig. 1: The rate shift type (speciation vs extinction)
##
####################################
fig = Figure(size = (220, 220), fontsize = 14, 
            figure_padding = (1,1,3,1));

colors = [:steelblue, "orange", "gray"]
labels = [
    L"\text{Shift in \lambda}", 
    L"\text{Shift in \mu}", 
    L"\text{Shift in both}"
    ]


##############
#
#   VIOLIN FIGURE 
#
############

ratios1 = zeros(length(datasets), 3)

for (dataset_index, name) in enumerate(names)
    model = models[name]

    rs1 = compute_ratios(name, 
            models[name],
            d_rates[name])

    ratios1[dataset_index,:] .= rs1
end


ax = Axis(fig[1,1], 
        ylabel = L"\hat{N}_\text{rate} / \hat{N}_\text{all}", 
        xlabel = L"\text{type of rate shift}",
        xgridvisible = false, 
        ygridvisible = false,
        xticks = (1:3, [L"\lambda", L"\mu", L"\text{both}"]),
        yticks = range(0.0,1.0; length = 5),
        titlesize = 9,
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelsize = 9,
        yticklabelsize = 9)

ylims!(ax, -0.05, 1.05)

## hard code the priors
priors = [1/11, 1/11, 9/11]
prior_ratios = zeros(length(datasets), 3)
for i in 1:size(prior_ratios)[1]
    prior_ratios[i,:] .= priors
end

is_larger1 = ratios1 .> prior_ratios
sum(is_larger1, dims = 1)

for (i, prior) in enumerate(priors)
    x = [i-0.5, i+0.5]
    y = [prior, prior]
    CairoMakie.lines!(ax, x, y, color = colors[i], linestyle = :solid, alpha = 0.5, label = "prior")
end

ratios2 = ratios1[.!isnan.(ratios1[:,1]),:]

xs = vcat([repeat([i], size(ratios2)[1]) for i in 1:3]...)
cs = vcat([repeat([i], size(ratios2)[1]) for i in colors]...)


CairoMakie.rainclouds!(ax, xs .- 0.25, vcat(ratios2...), 
                     color = cs, label = "Posterior",
                     clouds=nothing, boxplot_width=0.5,
                     jitter_width=0.35,
                     side_nudge = 0.5,
                     markersize=4)

#axislegend(ax)
#set_theme!(fig, figure_padding = 0)
fig

CairoMakie.save("figures/empirical_rate_shift_type_significant.pdf", fig)
