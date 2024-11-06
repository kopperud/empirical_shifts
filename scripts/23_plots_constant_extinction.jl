using Revise
using Pesto
using JLD2
using Glob
using DataFrames

fpaths = []


push!(fpaths, Glob.glob("output/simulations/constant_extinction/jld2/*.jld2"))
push!(fpaths, Glob.glob("output/simulations/constant_extinction2/jld2/*.jld2"))
push!(fpaths, Glob.glob("output/simulations/constant_speciation/jld2/*.jld2"))
push!(fpaths, Glob.glob("output/simulations/constant_speciation2/jld2/*.jld2"))

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


ratios = zeros(4, length(fpaths), 3)

using ProgressMeter

ntips = Vector{Int64}[]

prog = Progress(reduce(+, map(length, fpaths)));

dfs = DataFrame[]

for j in 1:4
    paths = fpaths[j]
    for path in paths
        x = load(path)
        iter = split(Base.basename(path), ".")[1]
        model = make_model(x)
        
        r = compute_ratios(x["N"], model)
        df = DataFrame(
            "ratio μ" => r[1],
            "ratio λ" => r[2],
            "ratio μ+λ" => r[3],
            "study" => j,
            "ntip" => x["ntip"],
            "iter" => iter,
            )
        push!(dfs, df)

        next!(prog)
    end
end
finish!(prog)

dfx = vcat(dfs...)


using CairoMakie
using Statistics




fig = Figure(size = (700, 700));

titles = [
    round(mean(ratios[:,i]); digits =  4) for i in 1:3
]

## iter over models
axs = []
for j in 1:4
    ## iter over rows
    df = filter(:study => x -> x == j, dfx)
    ntrees = size(df)[1]
     
    ax1 = Axis(fig[1,j], xlabel = "number of tips", xticklabelrotation = π/2, title = "$(ntrees) trees")
    ax2 = Axis(fig[2,j], xlabel = L"\hat{N}_\mu/\hat{N}", xticklabelrotation = π/2)#, title = "r = $(titles[1])")
    ax3 = Axis(fig[3,j], xlabel = L"\hat{N}_\lambda/\hat{N}", xticklabelrotation = π/2)#, title = "r = $(titles[2])")
    ax4 = Axis(fig[4,j], xlabel = L"\hat{N}_{\lambda+\mu}/\hat{N}")#", title = "r = $(titles[3])")

    if ntrees > 0

        hist!(ax1, df[!,:ntip])  
        hist!(ax2, df[!,"ratio μ"], bins = 10, color = "orange")
        hist!(ax3, df[!,"ratio λ"], bins = 10, color = "blue")
        hist!(ax4, df[!,"ratio μ+λ"], bins = 10, color = "gray")
        for ax in (ax2, ax3, ax4)
            xlims!(ax, (0.0, 1.0))
            #push!(axs, ax)
        end
    end
end

Label(fig[-1,1:2], "constant extinction (true model)")
Label(fig[-1,3:4], "constant speciation (true model)")

Label(fig[0,1], "inference sd=0.587")
Label(fig[0,2], "inference sd=0.2")
Label(fig[0,3], "inference sd=0.587")
Label(fig[0,4], "inference sd=0.2")

for j in 1:4
    colsize!(fig.layout, j, Relative(1/4))
end

fig


save("figures/partial_constant_models_ratios.pdf", fig)