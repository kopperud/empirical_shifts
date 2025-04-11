using JLD2
using CairoMakie
using Glob
using DataFrames

function delta_netdiv(d::Dict{String,Any})
    N = d["N"] 
    Ns = sum(N, dims = 1)[1,:,:]

    λ = d["lambda"]
    μ = d["mu"]

    r = λ .- μ 
    Δr = r .- r'

    delta = sum(Δr .* Ns)
    return(delta)
end

function delta_lambda(d::Dict{String,Any})
    N = d["N"] 
    Ns = sum(N, dims = 1)[1,:,:]

    λ = d["lambda"]
    Δλ = λ .- λ'

    delta = sum(Δλ .* Ns)
    return(delta)
end

function delta_mu(d::Dict{String,Any})
    N = d["N"] 
    Ns = sum(N, dims = 1)[1,:,:]

    μ  = d["mu"]
    Δμ = μ .- μ'

    delta = sum(Δμ .* Ns)
    return(delta)
end


backbone_fpaths = Glob.glob("output/simulations/grafts/backbone/jld2/*.jld2")
upshift_lambda_fpaths = Glob.glob("output/simulations/grafts/upshift_lambda/jld2/*.jld2")
upshift_mu_fpaths = Glob.glob("output/simulations/grafts/upshift_mu/jld2/*.jld2")
downshift_fpaths = Glob.glob("output/simulations/grafts/downshift/jld2/*.jld2")

fpaths = [
    backbone_fpaths,
    upshift_lambda_fpaths,
    upshift_mu_fpaths,
    downshift_fpaths,
]

data = [
    [load(fpath) for fpath in x] for x in fpaths
];

ntips = [[x["ntip"] for x in xs] for xs in data]
dlambdas = [[delta_lambda(x) for x in xs] for xs in data]
dmus = [[delta_mu(x) for x in xs] for xs in data]
dnetdivs = [[delta_netdiv(x) for x in xs] for xs in data]

fig = Figure(size = (550, 600))

side_titles = [
    L"\text{constant}",
    L"\text{upshift }(\lambda)",
    L"\text{upshift }(\mu)",
    L"\text{downshift }(\lambda)",
]

vertical_lines = [
    [0.0, 0.0, 0.0],
    [0.2, 0.0, 0.2],
    [0.0, -0.2, 0.2],
    [-0.2, 0.0, -0.2],
]

top_titles = [
    L"\text{speciation }(\Delta \lambda)", 
    L"\text{extinction }(\Delta \mu)", 
    L"\text{netdiv }(\Delta r)", 
]

xlabels = [
    L"\sum_{i,j}(\lambda_i - \lambda_j) \hat{N}_{ij}",
    L"\sum_{i,j}(\mu_i - \mu_j) \hat{N}_{ij}",
    L"\sum_{i,j}(r_i - r_j) \hat{N}_{ij}",
]

all_axs = []
for (row_index, (title, vl)) in enumerate(zip(titles, vertical_lines))
    ntip = ntips[row_index]
    dlambda = dlambdas[row_index]
    dmu = dmus[row_index]
    dnetdiv = dnetdivs[row_index]

    if row_index == 1
        title1 = top_titles[1]
        title2 = top_titles[2]
        title3 = top_titles[3]
    else
        title1 = ""
        title2 = ""
        title3 = ""
    end

    if row_index == 4
        xlab1 = xlabels[1]
        xlab2 = xlabels[2]
        xlab3 = xlabels[3]
    else
        xlab1 = ""
        xlab2 = ""
        xlab3 = ""
    end

    ax1 = Axis(fig[row_index,1],
        topspinevisible = false,
        rightspinevisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xticklabelrotation = π/2,
        title = title1,
        xlabel = xlab1,
        )
    ax2 = Axis(fig[row_index,2],
        topspinevisible = false,
        rightspinevisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xticklabelrotation = π/2,
        title = title2,
        xlabel = xlab2,
        )
    ax3 = Axis(fig[row_index,3],
        topspinevisible = false,
        rightspinevisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xticklabelrotation = π/2,
        title = title3,
        xlabel = xlab3,
        )

    hideydecorations!(ax2, ticks = false)
    hideydecorations!(ax3, ticks = false)

    if row_index < 4
         hidexdecorations!(ax1, ticks = false)
         hidexdecorations!(ax2, ticks = false)
         hidexdecorations!(ax3, ticks = false)
    end

    axs = [ax1, ax2, ax3]
    append!(all_axs, axs)

    ex = extrema(vcat(dlambda, dmu, dnetdiv))

    linkyaxes!(axs...)

    hist!(ax1, dlambda, label = L"\Delta \lambda", color = :steelblue, normalization = :none, bins = 25) #, strokecolor = :black, strokewidth=1)
    hist!(ax2, dmu, label = L"\Delta \mu", color = :red, bins = 25)
    hist!(ax3, dnetdiv, label = L"\Delta r", color = :green, bins = 25)

    for (ax, v) in zip(axs, vl)
        vlines!(ax, [v], color = :black, linestyle = :dash, label = L"\text{true}")
    end

    fig
end

for (row_index, label) in enumerate(side_titles)
    Label(fig[row_index,4], label, rotation = π/2)
end

## make the y-axes linked per row
linkxaxes!([all_axs[j] for j in [1,4,7,10]]...)
linkxaxes!([all_axs[j] for j in [2,5,8,11]]...)
linkxaxes!([all_axs[j] for j in [3,6,9,12]]...)

ylabel = Label(fig[1:4,0], L"\text{number of trees}", rotation = π/2)

elem_1 = [
    PolyElement(color = :steelblue, strokecolor = :black),
    PolyElement(color = :red, strokecolor = :black),
    PolyElement(color = :black, strokecolor = :black),
]
elem_2 = LineElement(color = :black, linestyle = :dash)
#elem_2 = PolyElement(color = :red, strokecolor = :black)
#elem_3 = PolyElement(color = :green, strokecolor = :black)
elem_4 = LineElement(color = :black, linestyle = :dash)

Legend(fig[5,1:3],
    [elem_1, elem_2],
    [
        L"\text{inferred change}", 
        L"\text{true change}",
        ],
    orientation = :vertical,
)

#xlabel = Label(fig[5,1:3], L"\text{inferred change} =  \sum_j \sum_i (r_i - r_j) \hat{N}_{ij}", rotation = 0)

for j in 1:4
    rowsize!(fig.layout, j, Relative(0.23))
end
#rowsize!(fig.layout, 5, Relative(0.05))
rowsize!(fig.layout, 5, Relative(0.07))

fig

save("figures/significant_branches_rate_change.pdf", fig)


fig2 = Figure()
ax = Axis(fig2[1,1],
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    xlabel = L"\text{change in rate}",
)
density!(ax, dlambda, color = (:teal, 0.1), strokecolor = :teal, strokearound = false, strokewidth = 3)
density!(ax, dmu, color = (:red, 0.1), strokecolor = :red, strokearound = false, strokewidth = 3)
density!(ax, dnetdiv, color = (:green, 0.1), strokecolor = :green, strokearound = false, strokewidth = 3)

Legend(fig2[1,2],
    [elem_1, elem_2, elem_3, elem_4],
    [L"\Delta \lambda", L"\Delta \mu", L"\Delta r", L"\text{true}"],
)

fig2

