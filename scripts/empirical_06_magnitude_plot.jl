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

function shift_size(model, N)
    r = model.λ .- model.μ
    Δr = N .* (r .- r')
    mean_magnitude = sum(Δr) 
    return(mean_magnitude)
end

function direction(model, N)
    r = model.λ .- model.μ

    Δr = r .- r'

    Δr_sign = sign.(Δr)

    mean_direction = sum(Δr_sign .* N) / sum(N)
    return(mean_direction)
end

########################
##
##         read files
##
##########################

inference = "empirical"

df = CSV.read("output/empirical_munged.csv", DataFrame)
df = df[df[!,:inference] .== inference,:]
df = df[df[!,:type] .== "pooled",:]

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
names = unique(df[!,:name])
fpaths = Glob.glob("output/" * inference * "/jld2/*.jld2")
@showprogress for fpath in fpaths
    name = split(Base.basename(fpath), ".")[1]
    x = JLD2.load(fpath)
    d[name] = x
end

models = Dict{String, BDSconstant}()
for name in names
    λ = d[name]["lambda"]
    μ = d[name]["mu"]
    η = d[name]["etaml"]
    models[name] = SSEconstant(λ, μ, η)
end

fpaths = ["data/empirical/" * name * ".tree" for name in names]
meta = CSV.read("data/empirical/Phylogenies for DeepILS project.csv", DataFrame)
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

fpaths = Glob.glob("output/" * inference * "/rates/*.csv")
rates = Dict{String, DataFrame}()
for fpath in fpaths
    name = split(Base.basename(fpath), ".")[1]
    rates_df = CSV.read(fpath, DataFrame)
    sort!(rates_df, :edge)
    rates[name] = rates_df[2:end,:]
end

##################
##
##     calculate magnitudes and direction
##
###############

#n_datasets = length(datasets)
n_datasets = length(d)
magnitudes = zeros(n_datasets)
directions = zeros(n_datasets)
heights = zeros(n_datasets)
treelengths = zeros(n_datasets)
N_per_time = zeros(n_datasets)

for (i, name) in enumerate(keys(d))
    println(name)
    model = models[name]
    shift_bf = rates[name][!,:shift_bf]

    heights[i] = maximum(datasets[name .* ".tree"].node_depth)
    treelengths[i] = sum(datasets[name * ".tree"].branch_lengths)
   
    is_signif = findall(shift_bf .> 10)

    Nsum = sum(d[name]["N"][is_signif,:,:], dims = 1)[1,:,:]
    m = shift_size(model, Nsum)
    dir = direction(model, Nsum)
    magnitudes[i] = m    
    directions[i] = dir
    N_per_time[i] = df[df[:,:name] .== name,:N_per_time][1]

end

#filter(:name => x -> x == "Vascular_Plants_Zanne2014", df)[1,:]

## write to a file

plotdf = DataFrame(
    "magnitude" => magnitudes,
    "direction" => directions,
    #"height" => heights,
    "name" => collect(keys(d)),
    #"N_per_time" => N_per_time,
)

meta2 = filter(:Filename => x -> !ismissing(x), meta)
meta2[:,:name] = map(x -> split(x, ".")[1], meta2.Filename)
dfx = innerjoin(plotdf, meta2, on = :name)
dfx2 = innerjoin(dfx, df, on = :name)
CSV.write("output/munged_magnitude.csv", dfx2)

##################
##
##   set up the makie figure
##
###############

fig2 = Figure(size = (600, 180), fontsize = 14, 
                figure_padding = (5,8,1,1));

#xt = collect(range(extrema(magnitudes_pooled)...; length = 5))
#xt = [-0.60, -0.30, 0.0, 0.30, 0.60, 0.90, 1.20]
xt = [-1.0, -0.50, 0.0, 0.50, 1.0]
xtl = [@sprintf("%.2f", x) for x in xt]

ax_direction = Axis(fig2[1,1], 
            ylabel = L"\text{frequency}", 
            xlabel = L"\text{shift direction}",
            title = L"\text{a) direction}",
            xgridvisible = false, 
            ygridvisible = false,
            xticks = (xt, xtl),
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelsize = 9,
            yticklabelsize = 9)
xlims!(ax_direction, -1.1, 1.1)

CairoMakie.hist!(ax_direction, 
    directions, bins = 15, color = "gray")


#xt = [-0.25, 0.0, 0.25, 0.50, 0.75, 1.0, 1.25]
xt = [-1.0, -0.50, 0.0, 0.50, 1.0]
xtl = [@sprintf("%.2f", x) for x in xt]

ax_magnitude = Axis(fig2[1,2], 
            #ylabel = L"\text{frequency}", 
            xlabel = L"\text{magnitude }(\Delta r)",
            title = L"\text{b) magnitude}",
            xgridvisible = false, 
            ygridvisible = false,
            xticks = (xt, xtl),
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelsize = 9,
            yticklabelsize = 9)
xlims!(ax_magnitude, -1.35, 1.35)

CairoMakie.hist!(ax_magnitude, 
    magnitudes, bins = 15, color = "gray")

    

xt = 10 .^ (collect(range(extrema(log10.(plotdf[!,:heights]))...; length = 5)))
xtl = [@sprintf("%.1f", x) for x in xt]


yr = collect(range(extrema(plotdf[!,:magnitudes])...; length = 5))
#yr = [-0.25, 0.0, 0.25, 0.50, 0.75, 1.0, 1.25]
yr = [-1.0, -0.50, 0.0, 0.50, 1.0]
yt = yr
ytl = [@sprintf("%.2f", y) for y in yt]

ax_scatter = Axis(fig2[1,3], 
        ylabel = L"\text{magnitude } (\Delta r)", 
        title = L"\text{c) time scaling}",
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

ylims!(ax_scatter, -1.35, 1.35)


yr = collect(range(extrema(magnitudes)...; length = 5))
yt = yr
ytl = [@sprintf("%.2f", y) for y in yt]


β, Varβ, ySE = ols_regression(log10.(heights), magnitudes)

x = collect(Pesto.lrange(extrema(heights)..., 20))
linefit(x) = β[1] + β[2]*log10.(x)

yarithmetic = linefit.(x)

#yVar = Varβ[1,1] .+ 2 .* x .* Varβ[1,2] .+ (x .^ 2) .* Varβ[2,2]
yVar = Varβ[1,1] .+ 2 .* log10.(x) .* Varβ[1,2] .+ (log10.(x) .^ 2) .* Varβ[2,2]
yupper = yarithmetic .+ 2*sqrt.(yVar)
ylower = yarithmetic .- 2*sqrt.(yVar)
y = linefit.(x)


CairoMakie.band!(ax_scatter, x, ylower, yupper, color = "#e0e0e0")
CairoMakie.lines!(ax_scatter, x, y; label = "OLS", markersize = 7, color = "gray", linestyle = :dash)

CairoMakie.scatter!(ax_scatter, 
                    heights, 
                    magnitudes,
                    color = "black",
                    markersize = 7)


colgap!(fig2.layout, 5)
rowgap!(fig2.layout, 5)
fig2

###################
n_positive = sum(magnitudes .> 0)
n_negative = sum(magnitudes .< 0)

println("$n_positive magnitudes are positive, while $n_negative are negative")


#CairoMakie.save("figures/magnitude_empirical.pdf", fig2)
#CairoMakie.save("figures/direction_magnitude_empirical.pdf",fig2)

