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


fpaths_rates = [
    glob("output/simulations/single_shift_grafts/backbone/rates/*.csv"),
    glob("output/simulations/single_shift_grafts/upshift/rates/*.csv"),
    glob("output/simulations/single_shift_grafts/downshift/rates/*.csv"),
]

using CSV
using DataFrames

df = CSV.read(fpaths_rates[2][1], DataFrame)

dn = df[!,:delta_netdiv]
dn = dn[.!isnan.(dn)]

sum(dn)

netdiv_mag(xs[1][1])


