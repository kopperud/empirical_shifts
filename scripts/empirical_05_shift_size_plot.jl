using Distributions
using Glob, DataFrames, CSV, RCall, ProgressMeter
using LaTeXStrings, Measures
using JLD2
using Printf
using CairoMakie
using Pesto

########################
##
##   read the data files
##
##########################

#inference = "empirical"
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
@showprogress for fpath in fpaths
    name = split(Base.basename(fpath), ".")[1]
    x = JLD2.load(fpath)
    d[name] = x

end

rates_fpaths = Glob.glob("output/" * inference * "/rates/*.csv")
d_rates = Dict{String, DataFrame}()
@showprogress for fpath in rates_fpaths
    name = split(Base.basename(fpath), ".")[1]
    df = CSV.read(fpath, DataFrame)
    sort!(df, :edge)
    filter!(:edge => x -> x != 0, df) ## remove the root edge
    is_signif = df[!,:shift_bf] .> 10
    df[!,:is_signif] = is_signif

    d_rates[name] = df
end


models = Dict{String, SSEconstant}()
for name in names
    λ = d[name]["lambda"]
    μ = d[name]["mu"]
    η = d[name]["etaml"]
    models[name] = SSEconstant(λ, μ, η)
end


fpaths = ["data/empirical/" * name * ".tree" for name in names]

meta = CSV.read("data/empirical/metadata.csv", DataFrame)
ρs = Dict()
for row in eachrow(meta)
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

#heights = [maximum(d.node_depth) for (key, d) in datasets]

####################################
##
##   Fig. 2: The empirical dataset exploration (not for any of the hypotheses in particular)
##
####################################
fig = Figure(size = (480, 280), fontsize = 14, 
            figure_padding = (1,1,1,1));

colors = [:steelblue, "orange", "gray"]
labels = [
    L"\text{Shift in \lambda}", 
    L"\text{Shift in \mu}", 
    L"\text{Shift in both}"
    ]

g = fig[1,1] = GridLayout()

name_subset = [
    "Actinopterygii_Rabosky2018",
    "Mammalia_AlvarezCarretero2022",
    "Rosidae_Sun2020",
    "Chondrichthyes_Stein2018",
    "Squamata_Zheng2016",
    "Asteraceae_Palazzesi2022",
    "Agaricomycetes_Varga2019",
    "Anura_Portik2023",
    "Aves_Quintero2022",
    "Polypodiophyta_Nitta2022",
    "Lecanoromycetes_Nelsen2020",
]

nbins1 = [
    20, 20, 10,
    4, 6, 20, 18,
    6, 14, 4, 20
]

titles = []
dnames = collect(keys(datasets))
for dname in dnames
    append!(titles, [split(dname, "_")[1]])
end

axs = []
q = 1
for i in 1:3, j in 1:4
    if (i == 1) & (j == 4)
        #ax = Axis(fig[i,j], scene = false)
        #hidedecorations!(ax)
    else
        if i < 3
            xlabel = ""
        else
            xlabel = L"\Delta r"
        end

        title = split(name_subset[q], "_")[1]

        ax = Axis(g[i,j], 
        xgridvisible = false, 
        ygridvisible = false,
        title = title,
        titlesize = 9,
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelrotation = π/2,
        xticklabelsize = 9,
        yticklabelsize = 9)

        if i < 3
            hidexdecorations!(ax, ticks = false)
        end
    
        #=if j > 1
            hideydecorations!(ax, ticks = false, ticklabels = false)
        end=#

        append!(axs, [ax])
        q += 1
    end
end

for (q, name) in enumerate(name_subset)
    model = models[name]
    println(q, "\t", name)

    N_signif = d[name]["N"][d_rates[name][!,:is_signif],:,:]
    Nmatrix = sum(N_signif, dims = 1)[1,:,:]
    nbins = nbins1[q]

    Ns = zeros(3,nbins)
    filters = ["extinction", "speciation", ""]

    if true
        r = model.λ .- model.μ
        Δr = r .- r'
        netdiv_extrema = extrema(Δr)
    else 
        netdiv_extrema = [-1.2, 1.2]
    end

    dfs = []
    for (i, filter) in enumerate(filters)
        mids, bins = makebins(Nmatrix, model, netdiv_extrema...; filter = filter, nbins = nbins)
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

    bar_df = vcat(dfs...)

    barplot!(axs[q], 
        bar_df[!, "mids"], 
        bar_df[!, "Δr"], 
        stack = bar_df[!, "subset"],
        color = [colors[x] for x in bar_df[!, "subset"]],
        )
end

elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
leg = Legend(g[1,4], elements, labels, labelsize=9, 
            patchsize = (10.0f0, 10.0f0),
            framevisible = false)

leg.tellheight = true
leg.tellwidth = true

linkxaxes!(axs[1], axs[4], axs[8])
linkxaxes!(axs[2], axs[5], axs[9])
linkxaxes!(axs[3], axs[6], axs[10])
linkxaxes!(axs[7], axs[11])

ylabel = Label(g[1:3, 0], L"\text{number of rate shifts }(\hat{N})", rotation = π/2)
xlabel = Label(g[4, 1:4], L"\text{shift size in net diversification }(\Delta r)")




for col in 1:3
    colsize!(g, col, Relative(0.23))
end
colsize!(g, 0, Relative(0.05))

for i in 1:3
    rowsize!(g, i, Relative(0.30))
end
rowsize!(g, 4, Relative(0.1))

colgap!(g, 5)
rowgap!(g, 3)

fig
#set_theme!(fig, figure_padding = 0)
#CairoMakie.save("figures/fig1_empirical_significant.pdf", fig)

################################
##
## plot each phylogeny individually
##
################################


for name in keys(models)  
    model = models[name]

    #Nmatrix = N[dataset_index,:,:,1]
    Nmatrix = sum(d[name]["N"], dims = 1)[1,:,:]
    nbins = 20

    Ns = zeros(3,nbins)
    filters = ["extinction", "speciation", ""]

    dfs = []
    #netdiv_extrema = [-1.2, 1.2]
    r = model.λ .- model.μ
    Δr = r .- r'
    netdiv_extrema = extrema(Δr)

    for (i, filter) in enumerate(filters)
        mids, bins = makebins(Nmatrix, model, netdiv_extrema...; filter = filter, nbins = nbins)
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

    bar_df = vcat(dfs...)

    fig = Figure(size = (200,200))

    ax = Axis(fig[1,1], 
            xgridvisible = false, 
            ygridvisible = false,
            title = name,
            titlesize = 9,
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelrotation = π/2,
            xticklabelsize = 9,
            xlabel = "shift size (netdiv)",
            ylabel = "number of shifts",
            yticklabelsize = 9)

    barplot!(ax, 
        bar_df[!, "mids"], 
        bar_df[!, "Δr"], 
        stack = bar_df[!, "subset"],
        color = [colors[x] for x in bar_df[!, "subset"]],
        )

    
    save("/tmp/empirical_shiftfigs/$name.pdf", fig)    

end


###
## plot them in a big plot
##

fig3 = Figure(size = (800, 800));

colors = [:steelblue, "orange", "gray"]

global q = 1
qnames = [keys(models)...]

titles2 = []
dnames = collect(keys(models))
for dname in dnames
    item = [split(dname, "_")[1]] 
    item = replace(item, "Vascular" => "Vascular Plants")
    append!(titles2, item)
end


axs = []
for i in 1:5
    for j in 1:5
        if q <= length(keys(models))
            ax = Axis(fig3[i,j], 
                xgridvisible = false, 
                ygridvisible = false,
                titlesize = 7,
                title = titles2[q],
                topspinevisible = false,
                rightspinevisible = false,
                xticklabelrotation = π/2,
                xticklabelsize = 9,
                yticklabelsize = 9)
            push!(axs, ax)
            q += 1
        end
    end
end

i = 1
for name in keys(models)
    rates = deepcopy(d_rates[name])
    filter!(row -> row[:shift_bf] > 10, rates)

    hist!(axs[i], rates[!,:delta_netdiv])
   
    #=
    model = models[name]
    Nmatrix = sum(d[name]["N"], dims = 1)[1,:,:]
    nbins = 20

    Ns = zeros(3,nbins)
    filters = ["extinction", "speciation", ""]

    dfs = []
    if false
        netdiv_extrema = [-1.5, 1.5]
    else
        r = model.λ .- model.μ
        Δr = r .- r'
        netdiv_extrema = extrema(Δr)
    end

    for (i, filter) in enumerate(filters)
        mids, bins = makebins(Nmatrix, model, netdiv_extrema...; filter = filter, nbins = nbins)
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

    bar_df = vcat(dfs...)


    ax = axs[i]
    barplot!(ax, 
        bar_df[!, "mids"], 
        bar_df[!, "Δr"], 
        stack = bar_df[!, "subset"],
        color = [colors[x] for x in bar_df[!, "subset"]],
        )
       
        =#
    i += 1
end

## add x-axis label and y-axis label
ylabel = Label(fig3[1:5, 0], L"\text{number of rate shifts }(\hat{N})", rotation = π/2)
xlabel = Label(fig3[6, 1:5], L"\text{shift size in net diversification }(\Delta r)")

colgap!(fig3.layout, 4)
rowgap!(fig3.layout, 2)

fig3

CairoMakie.save("figures/shiftsize_free_xlimit.pdf", fig3)
#CairoMakie.save("figures/shiftsize_fixed_xlimit.pdf", fig3)

########################
##
## Plot the number of significant rate shift events
## as a function of the groups
##
########################


kingdom = Dict{String, String}()

kingdom["Actinopterygii_Rabosky2018"] = "Animalia"
kingdom["Agaricomycetes_Varga2019"] = "Fungi"
kingdom["Anthophila_HenriquezPiskulich"] = "Animalia"
kingdom["Anura_Portik2023"] = "Animalia"
kingdom["Aristolochiaceae_Allio2021"] = "Plantae"
kingdom["Asteraceae_Palazzesi2022"] = "Plantae"
kingdom["Aves_Quintero2022"] = "Animalia"
kingdom["Caryophyllales_Smithetal2017_align_w_shtest_treepl_3samp"] = "Plantae"
kingdom["Chondrichthyes_Stein2018"] = "Animalia"
kingdom["Cornales_Rose2018"] = "Plantae"
kingdom["Cycadaceae_Liu2021"] = "Plantae"
kingdom["Ericales_Rose2018"] = "Plantae"
kingdom["Gesnerioideae_SerranoSerrano2017"] = "Plantae"
kingdom["Lecanoromycetes_Nelsen2020"] = "Fungi"
kingdom["Mammalia_AlvarezCarretero2022"] = "Animalia"
kingdom["Orobanchaceae_Mortimer2022"] = "Plantae"
kingdom["Phasmatodea_Bank2022"] = "Animalia"
kingdom["Pinophyta_Leslie2018"] = "Plantae"
kingdom["Poaceae_Spriggs2015"] = "Plantae"
kingdom["Polypodiophyta_Nitta2022"] = "Plantae"
kingdom["Rhopalocera_Kawahara2023"] = "Animalia"
kingdom["Rosidae_Sun2020"] = "Plantae"
kingdom["Salvia_Kriebel2019"] = "Plantae"
kingdom["Squamata_Zheng2016"] = "Animalia"




rainclouds_datasets = reverse([
    "Actinopterygii_Rabosky2018",
    "Agaricomycetes_Varga2019",
    "Anthophila_HenriquezPiskulich",
    "Anura_Portik2023",
    "Aristolochiaceae_Allio2021",
    "Asteraceae_Palazzesi2022",
    "Aves_Quintero2022",
    "Caryophyllales_Smithetal2017_align_w_shtest_treepl_3samp",
    "Chondrichthyes_Stein2018",
    "Cornales_Rose2018",
    "Cycadaceae_Liu2021",
    "Ericales_Rose2018",
    "Gesnerioideae_SerranoSerrano2017",
    "Lecanoromycetes_Nelsen2020",
    "Mammalia_AlvarezCarretero2022",
    "Orobanchaceae_Mortimer2022",
    "Phasmatodea_Bank2022",
    "Pinophyta_Leslie2018",
    "Poaceae_Spriggs2015",
    "Polypodiophyta_Nitta2022",
    "Rhopalocera_Kawahara2023",
    "Rosidae_Sun2020",
    "Salvia_Kriebel2019",
    "Squamata_Zheng2016",
])

colormap = Dict("Animalia" => :gray, "Plantae" => (:green, 0.8), "Fungi" => (:orange, 0.8))
color = [colormap[kingdom[name]] for name in rainclouds_datasets]

rainclouds_labels = reverse([
    L"\text{Actinopterygii}",
    L"\text{Agaricomycetes}",
    L"\text{Anthophila}",
    L"\text{Anura}",
    L"\text{Aristolochiaceae}",
    L"\text{Asteraceae}",
    L"\text{Aves}",
    L"\text{Caryophyllales}",
    L"\text{Chondrichthyes}",
    L"\text{Cornales}",
    L"\text{Cycadaceae}",
    L"\text{Ericales}",
    L"\text{Gesnerioideae}",
    L"\text{Lecanoromycetes}",
    L"\text{Mammalia}",
    L"\text{Orobanchaceae}",
    L"\text{Phasmatodea}",
    L"\text{Pinophyta}",
    L"\text{Poaceae}",
    L"\text{Polypodiophyta}",
    L"\text{Rhopalocera}",
    L"\text{Rosidae}",
    L"\text{Salvia}",
    L"\text{Squamata}",
])

fig2 = Figure(size = (500, 350))

yt = collect(1:length(rainclouds_datasets)/2)
ytl1 = rainclouds_labels[1:12]
ytl2 = rainclouds_labels[13:24]

ax1 = Axis(fig2[1,1], 
    #xlabel = L"\text{shift size in net diversification }(\Delta r)", 
    xgridvisible = false, ygridvisible = false,
    yticks = (yt, ytl2),
    topspinevisible = false,
    rightspinevisible = false,
    xticklabelrotation = π/2,
    xticklabelsize = 9,
    yticklabelsize = 9)

ax2 = Axis(fig2[1,2], 
    #xlabel = L"\text{shift size in net diversification }(\Delta r)", 
    xgridvisible = false, ygridvisible = false,
    yticks = (yt, ytl1),
    topspinevisible = false,
    rightspinevisible = false,
    xticklabelrotation = π/2,
    xticklabelsize = 9,
    yticklabelsize = 9)

    
Label(fig2[2, 1:2], L"\text{shift size in net diversification }(\Delta r)", rotation = 0)


shifts = Dict{String, Vector{Float64}}()
for (name, rates) in d_rates
    rates_signif = rates[rates[!,:is_signif],:]
    shifts[name] = rates_signif[!,:delta_netdiv]
end

xs1 = Float64[]
ys1 = Float64[]

xs2 = Float64[]
ys2 = Float64[]
wc = Makie.wong_colors()
#colors = reverse([:gray, :gray, :gray, :gray, :gray, :gray, wc[3], wc[3], wc[3], wc[3], wc[4], wc[4]])
cs1 = []
cs2 = []
for (i, name) in enumerate(rainclouds_datasets[1:12])
    for item in shifts[name]
        push!(xs1, i)
        push!(ys1, item)
        push!(cs1, color[i])
    end
end

for (i, name) in enumerate(rainclouds_datasets[13:24])
    for item in shifts[name]
        push!(xs2, i)
        push!(ys2, item)
        push!(cs2, color[13:24][i])
        #push!(cs, colors[i])
    end
end


rainclouds!(ax2, xs1 .+ 0.25, ys1, 
    clouds = nothing, 
    #color = :gray,
    color = cs1,
    boxplot_width = 0.5, cloud_width = 0.5,
    orientation = :horizontal,
    jitter_width=0.35,
    side_nudge = 0.55,
    markersize=3,
    )

rainclouds!(ax1, xs2 .+ 0.25, ys2, 
    clouds = nothing, 
    #color = :gray,
    color = cs2,
    boxplot_width = 0.5, cloud_width = 0.5,
    orientation = :horizontal,
    jitter_width=0.35,
    side_nudge = 0.55,
    markersize=3,
    )

linkxaxes!(ax1, ax2)


Legend(fig2[1,3], 
    [
        MarkerElement(color = :gray, marker = :circle, markersize = 9),
        MarkerElement(color = :green, marker = :circle, markersize = 9),
        MarkerElement(color = :orange, marker = :circle, markersize = 9),
    ],
    ["Animalia", "Plantae", "Fungi"],
    patchsize = (8, 8), 
    rowgap = 5,
)

for i in 1:2
    colsize!(fig2.layout, i, Relative(0.4))
end
colsize!(fig2.layout, 3, Relative(0.2))

fig2


save("figures/all_datasets_shiftsizes.pdf", fig2)

########################
##
## Plot the number of significant rate shift events
## as a function of the clades
##
########################


name_subset2 = [
    "Actinopterygii_Rabosky2018",
    "Mammalia_AlvarezCarretero2022",
    "Rosidae_Sun2020",
    "Polypodiophyta_Nitta2022",
    "Chondrichthyes_Stein2018",
    "Squamata_Zheng2016",
    "Asteraceae_Palazzesi2022",
    "Agaricomycetes_Varga2019",
    "Anura_Portik2023",
    "Aves_Quintero2022",
    "Poaceae_Spriggs2015",
    "Lecanoromycetes_Nelsen2020",
]



fig3 = Figure(size = (480, 280), fontsize = 14, 
            figure_padding = (1,1,1,1));

nbins1 = [
    20, 20, 10, 10,
    4, 6, 20, 18,
    6, 14, 4, 20
]

axs = []
q = 1
for i in 1:3, j in 1:4
    if i < 3
        xlabel = ""
    else
        xlabel = L"\Delta r"
    end

    title = split(name_subset2[q], "_")[1]

    ax = Axis(fig3[i,j], 
    xgridvisible = false, 
    ygridvisible = false,
    title = title,
    titlesize = 9,
    topspinevisible = false,
    rightspinevisible = false,
    xticklabelrotation = π/2,
    xticklabelsize = 9,
    yticklabelsize = 9)

    if i < 3
        hidexdecorations!(ax, ticks = false)
    end

    append!(axs, [ax])
    q += 1
end

nbins3 = [
    20, 20, 14, 20,
    2, 6, 20, 4,
    7, 14, 10, 10
]

for (i, name) in enumerate(name_subset2)
    ax = axs[i]
    hist!(ax, shifts[name], color = :gray, bins = nbins3[i])
end

linkxaxes!(axs[1], axs[5], axs[9])
linkxaxes!(axs[2], axs[6], axs[10])
linkxaxes!(axs[3], axs[7], axs[11])
linkxaxes!(axs[4], axs[8], axs[12])

ylabel = Label(fig3[1:3, 0], L"\text{no. significant rate shift events}", rotation = π/2)
xlabel = Label(fig3[4, 1:4], L"\text{shift size in net diversification }(\Delta r)")

for col in 1:3
    colsize!(fig3.layout, col, Relative(0.23))
end
colsize!(fig3.layout, 0, Relative(0.05))

for i in 1:3
    rowsize!(fig3.layout, i, Relative(0.30))
end
rowsize!(fig3.layout, 4, Relative(0.1))

colgap!(fig3.layout, 5)
rowgap!(fig3.layout, 3)

fig3

save("figures/fig1_empirical_atomized.pdf", fig3)



## compute the variance
kingdom = Dict{String, String}()
for row in eachrow(meta)
    fn = row[:Filename]
    if ismissing(fn)
        continue
    end
    name = replace(fn, ".tree" => "")
    kingdom[name] = row[:Kingdom]
end



name = "Sigmodontinae_VallejosGarrido2023"


mags = Float64[]
vars = Float64[]
heights = Float64[]
names = String[]
kingdoms = String[]
netdivs = Float64[]
sp_rates = Float64[]

@showprogress for name in keys(models)

    model = models[name]
    N = sum(d[name]["N"], dims = 1)[1,:,:]

    Δr = model.λ .- model.μ

    mag = sum(Δr .* N) / sum(N)
    var = sum((Δr .* N .- mag) .^2) #/ sum(N)
    push!(mags, mag)
    push!(vars, var)
    height = maximum(datasets[name * ".tree"].node_depth)
    push!(heights, height)
    push!(names, name)
    push!(kingdoms, kingdom[name])

    fname = "output/" * inference * "/newick/" * name * ".tre"
    @rput fname
    R"""
    library(treeio)
    tr <- read.beast.newick(fname)
    xdf <- tr@data
    xdf <- xdf[order(xdf$edge),]
    netdiv <- sum(xdf$mean_netdiv * tr@phylo$edge.length) / sum(tr@phylo$edge.length)
    speciationrate <- sum(xdf$mean_lambda * tr@phylo$edge.length) / sum(tr@phylo$edge.length)
    """
    @rget netdiv
    @rget speciationrate

    push!(netdivs, netdiv)
    push!(sp_rates, speciationrate)
end

newdf = DataFrame(
    "magnitude" => mags,
    "variance" => vars,
    "height" => heights,
    "name" => names,
    "kingdom" => kingdoms,
    "mean_netdiv" => netdivs,
    "mean_speciation" => sp_rates,
)

plants = filter(:kingdom => x -> x == "Plantae", newdf)
animals = filter(:kingdom => x -> x == "Animalia", newdf)
fungi = filter(:kingdom => x -> x == "Fungi", newdf)

fig2 = Figure();
ax1 = Axis(fig2[1,1], 
        xscale = log10, 
        yscale = log10,
        xgridvisible = false, 
        ygridvisible = false,
        titlesize = 9,
        topspinevisible = false,
        rightspinevisible = false,
        xlabel = L"\text{tree height (Ma)}",
        ylabel = L"\text{variation in shift size}",
        xticklabelrotation = π/2,
        xticklabelsize = 9,
        yticklabelsize = 9);
ax2 = Axis(fig2[1,2], 
        xscale = log10, 
        yscale = log10,
        xgridvisible = false, 
        ygridvisible = false,
        titlesize = 9,
        topspinevisible = false,
        rightspinevisible = false,
        xlabel = L"\text{mean sp rate}",
        ylabel = L"\text{variation in shift size}",
        xticklabelrotation = π/2,
        xticklabelsize = 9,
        yticklabelsize = 9);
#scatter!(ax, heights, sqrt.(vars))
for (kng, label) in zip((plants, animals, fungi), ("Plantae", "Animalia", "Fungi"))
    scatter!(ax1, kng[!,:height], kng[!,:variance], label = label)
    scatter!(ax2, kng[!,:mean_speciation], kng[!,:variance], label = label)
end
axislegend(ax1; position = :lt);
#axislegend(ax2; position = :ct);
fig2




x = Float64[]
y = Float64[]

vars = Dict(key => var(val[!,:mean_netdiv]) for (key, val) in d_rates)

for row in eachrow(df)
    name = row[:name]
    v = vars[name]
    println(name)

    push!(x, row[:height])
    push!(y, v)
    println(typeof(v))
end

f = Figure()
ax = Axis(f[1,1],
xscale = log10,
yscale = log10)
scatter!(ax, x, y)
f

y

