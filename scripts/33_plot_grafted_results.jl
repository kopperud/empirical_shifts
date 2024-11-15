using Revise
using CairoMakie
using Glob
using JLD2
using Statistics
using LaTeXStrings

function netdiv_mag(d::Dict{String,Any})
    N = d["N"] 
    λ = d["lambda"]
    μ = d["mu"]

    r = λ .- μ 
    Δr = r .- r'

    nom = sum(Δr .* N)
    denom = sum(N)
    m = nom/denom

    return(denom)
end

function speciation_mag(d::Dict{String,Any})
    N = d["N"] 
    λ = d["lambda"]
    μ = d["mu"]

    Δλ = λ .- λ'

    nom = sum(Δλ .* N)
    denom = sum(N)
    m = nom/denom
    return(m)
end


fpaths = [
    glob("output/simulations/single_shift_grafts/backbone/jld2/*.jld2"),
    glob("output/simulations/single_shift_grafts/upshift/jld2/*.jld2"),
    glob("output/simulations/single_shift_grafts/downshift/jld2/*.jld2"),
]

xs = [[load(fpath) for fpath in fpaths[i]] for i in 1:3]

mags = [[netdiv_mag(x) for x in xs[i]] for i in 1:3]

etas = [[load(fpath, "etaml") for fpath in fpaths[i]] for i in 1:3]
eta_times_tl = [[load(fpath, "etaml")*load(fpath, "treelength") for fpath in fpaths[i]] for i in 1:3]


fig = Figure(size = (500, 450));
axs = [Axis(fig[i,1]) for i in 1:3]
axs2 = [Axis(fig[i,2], xticks = ([-4, -2, 0, 2], [L"10^{-4}", L"10^{-2}", L"10^0", L"10^2"])) for i in 1:3]
labels = ["backbone", "upshift", "downshift"]
hs = []
hs2 = []
colors = [:black, :orange, :green]
for (i, label, color) in zip(1:3, labels, colors)
    h = hist!(axs[i], mags[i], label = label, color = color, bins = 30)
    push!(hs, h)

    #h2 = hist!(axs2[i], log10.(etas[i]), label = label, color = color, bins = 30)
    h2 = hist!(axs2[i], log10.(eta_times_tl[i]), label = label, color = color, bins = 30)
    #axislegend(axs2[i])
    push!(hs2, h2)
end

lines!(axs[1], [0.0, 0.0], [0.0, 100], linewidth = 2, linestyle = :dash, color = :red, label = "true")
lines!(axs[2], [0.2, 0.2], [0.0, 45], linewidth = 2, linestyle = :dash, color = :red, label = "true")
lines!(axs[3], [-0.2, -0.2], [0.0, 65], linewidth = 2, linestyle = :dash, color = :red, label = "true")
for i in 1:3
    axislegend(axs[i])
end

lines!(axs2[2], [0.0, 0.0], [0.0, 55], linewidth = 2, linestyle = :dash, color = :red)
lines!(axs2[3], [0.0, 0.0], [0.0, 30], linewidth = 2, linestyle = :dash, color = :red)

linkxaxes!(axs...)
linkxaxes!(axs2...)
ylab = Label(fig[1:3,0], "number of simulated trees", rotation = π/2, tellheight = false, tellwidth = true)
lab = Label(fig[4,1], "magnitude (Δr)", tellwidth = false, tellheight = true)
lab = Label(fig[4,2], L"\text{E}[N] = \eta \times \text{treelength}", tellwidth = false, tellheight = true)
fig

#save("figures/grafted_trees_results.pdf", fig)

[mean(mags[i]) for i in 1:3]
## squared errors
[(m-(-0.2))]

fpaths

x = load(fpaths[2][1])

1 / x["treelength"]

#############################################
##
##      compute mean(netdiv_tips) - netdiv_root
##
#############################################

fpaths2 = [
    glob("output/simulations/single_shift_grafts2/backbone/jld2/*.jld2"),
    glob("output/simulations/single_shift_grafts2/upshift/jld2/*.jld2"),
    glob("output/simulations/single_shift_grafts2/downshift/jld2/*.jld2"),
]

using RCall

R"""
library(ape)
fp_up <- Sys.glob("data/simulations/single_shift_grafts/upshift/*.tre")
tr_up <- lapply(fp, read.tree)

fp_down <- Sys.glob("data/simulations/single_shift_grafts/downshift/*.tre")
tr_down <- lapply(fp, read.tree)

netdiv_root <- 0.09
netdiv_upshift <- 0.29
netdiv_downshift <- -0.11

foo_up <- function(tr){
    ntip <- length(tr$tip.label)
    n_shifted <- sum(grepl("x", tr$tip.label))
    netdiv_tip <- (n_shifted * netdiv_upshift + (ntip - n_shifted) * netdiv_root) / ntip
    true_delta_netdiv <- netdiv_tip - netdiv_root
    return(true_delta_netdiv)
}

foo_down <- function(tr){
    ntip <- length(tr$tip.label)
    n_shifted <- sum(grepl("x", tr$tip.label))
    netdiv_tip <- (n_shifted * netdiv_downshift + (ntip - n_shifted) * netdiv_root) / ntip
    true_delta_netdiv <- netdiv_tip - netdiv_root
    return(true_delta_netdiv)
}

delta_netdivs_upshift <- sapply(tr_up, foo_up)
delta_netdivs_downshift <- sapply(tr_down, foo_down)
"""
@rget delta_netdivs_upshift 
@rget delta_netdivs_downshift 

xs = [[load(fpath) for fpath in fpaths2[i]] for i in 1:3]

xs[1][1]

delta_netdiv(x::Dict{String,Any}) = mean(x["netdiv_tips"]) - x["netdiv_root"]
Δnetdivs = [[delta_netdiv(x) for x in xs[i]] for i in 1:3]

using Printf

fig2 = Figure(size = (400, 500))
mean_backbone = round(mean(Δnetdivs[1]); digits = 4)
mean_upshift = round(mean(Δnetdivs[2]); digits = 4)
mean_downshift = round(mean(Δnetdivs[3]); digits = 4)
ax1 = Axis(fig2[1,1], title = "backbone (mean = $mean_backbone)")
ax2 = Axis(fig2[2,1], title = "upshift (mean = $mean_upshift)")
ax3 = Axis(fig2[3,1], title = "downshift (mean = $mean_downshift)", 
    xlabel = L"\frac{1}{n}\left (\sum_{m=1}^n r_m(t_0) \right ) - r_\text{root}")
axs = [ax1, ax2, ax3]
for (i, ax) in enumerate(axs)
    Δnetdiv = Δnetdivs[i] 
    hist!(ax, Δnetdiv, bins = 25, label = "inferred", color = :black, normalization = :pdf)
end
hist!(axs[2], delta_netdivs_upshift, bins = 25, label = "true", color = (:orange, 0.7), normalization = :pdf)
hist!(axs[3], delta_netdivs_downshift, bins = 25, label = "true", color = (:orange, 0.7), normalization = :pdf)
linkxaxes!(axs...)
axislegend(axs[3])

fig2

save("figures/")







using CSV
using DataFrames

df = CSV.read(fpaths_rates[2][1], DataFrame)

dn = df[!,:delta_netdiv]
dn = dn[.!isnan.(dn)]

sum(dn)

netdiv_mag(xs[1][1])


