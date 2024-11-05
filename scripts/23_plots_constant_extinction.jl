using Revise
using Pesto
using JLD2
using Glob
using DataFrames

fpaths = Glob.glob("output/simulations/constant_extinction/jld2/*.jld2")


N = sum(N1, dims = 1)[1,:,:]


function compute_ratios(N, model)

    Nmatrix = sum(N, dims = 1)[1,:,:]
    nbins = 14

    Ns = zeros(3,nbins)
    filters = ["extinction", "speciation", ""]
    limits = [-1.2, 1.2]

    dfs = []
    for (i, filter) in enumerate(filters)
        mids, bins = makebins(Nmatrix, model, limits...; filter = filter, nbins = nbins)
        Ns[i,:] = bins[:,3]
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

function make_model(d::Dict)
    λ = d["lambda"]
    μ = d["mu"]
    η = d["etaml"]
    BDSconstant(λ, μ, η)
end


ratios = zeros(length(fpaths), 3)

using ProgressMeter

ntips = Int64[]
@showprogress for (i, fpath) in enumerate(fpaths)
    x = load(fpath)
    model = make_model(x)

    r = compute_ratios(x["N"], model, "")
    ratios[i,:] = r
    push!(ntips, x["ntip"])
end

using CairoMakie
using Statistics



[mean(ratios[:,i]) for i in 1:3]

fig = Figure(size = (400, 400));

titles = [
    round(mean(ratios[:,i]); digits =  4) for i in 1:3
]


ax1 = Axis(fig[1,1], xlabel = "number of tips", xticklabelrotation = π/2, title = "$(length(fpaths)) trees")
ax2 = Axis(fig[1,2], xlabel = "ratio µ", title = "r = $(titles[1])")
ax3 = Axis(fig[2,1], xlabel = "ratio λ", title = "r = $(titles[2])")
ax4 = Axis(fig[2,2], xlabel = "ratio µ+λ", title = "r = $(titles[3])")


hist!(ax1, ntips)  
hist!(ax2, ratios[:,1], bins = 10, color = "orange")
hist!(ax3, ratios[:,2], bins = 10, color = "blue")
hist!(ax4, ratios[:,3], bins = 10, color = "gray")

linkxaxes!(ax2, ax3, ax4)

fig

save("figures/constant_extinction_preliminary.pdf", fig)