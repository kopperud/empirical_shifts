using Revise
using Pesto
using JLD2
using ProgressMeter
using Glob
using DataFrames
using CSV

fpaths = Glob.glob("data/simulations/age_scaling_effect/h90_*.tre")

datasets = Dict{Tuple{Int64, Int64}, SSEdata}()

@showprogress for (i, fpath) in enumerate(fpaths)
    phy = readtree(fpath)
    ρ = 1.0
    dataset = SSEdata(phy, ρ)
    datasets[(90, i)] = dataset
end


## load the fitted model
models_fitted = Dict{Tuple{Int64, Int64}, SSEconstant}()
models_true = Dict{Tuple{Int64, Int64}, SSEconstant}()
@showprogress for i in 1:500
    fpath = "output/simulations/age_scaling_effect/jld2/h90_$i.jld2"

    if isfile(fpath)

        λ = JLD2.load(fpath, "lambda")
        μ = JLD2.load(fpath, "mu")
        η = JLD2.load(fpath, "etaml")
        model_fitted = SSEconstant(λ, μ, η)
        models_fitted[(90, i)] = model_fitted
    
        r0 = 0.04
        r = [r0, 0.07, 0.10]
        ϵ = 2/3
        λ = r ./ (1 - ϵ)
        μ = λ .- r
        η = r0 / 50

        model_true = SSEconstant(λ, μ, η)
        models_true[(90, i)] = model_true
    end
end

times = Dict{Tuple{Int64,Int64},Vector{Float64}}()
mean_shift_rate_fitted = Dict{Tuple{Int64,Int64},Vector{Float64}}()
mean_shift_rate_true = Dict{Tuple{Int64,Int64},Vector{Float64}}()

@showprogress for (h, i) in keys(models_fitted)
    x, sr_fitted = Pesto.shift_rate_through_time(models_fitted[90,i], datasets[90,i]; condition = [])
    times[90,i] = x
    mean_shift_rate_fitted[90,i] = sr_fitted

    x, sr_true = Pesto.shift_rate_through_time(models_true[90,i], datasets[90,i]; condition = [])
    times[90,i] = x
    mean_shift_rate_true[90,i] = sr_true
end


for (h, i) in keys(models_fitted)
    df = DataFrame(
        "times" => times[90,i], 
        "nshift_fitted" => mean_shift_rate_fitted[90,i], 
        "nshift_true" => mean_shift_rate_true[90,i])
    fpath = "output/simulations/age_scaling_effect/shift_rate_through_time/h90_$i.csv"
    CSV.write(fpath, df) 
end
## quick and dirty plots

using LaTeXStrings
using CairoMakie

figs = []

for (h, i) in keys(models_true)
    ntip = length(datasets[90,i].tiplab)
    t = times[90,i]
    y_fitted = mean_shift_rate_fitted[90,i]
    η_fitted = models_fitted[90,i].η

    fig1 = Figure()

    ax = Axis(fig1[1,1],
        xreversed = true, 
        title = "idx = $i, ntip = $ntip",
        #yscale = :log10,
        xgridvisible = false,
        ygridvisible = false,
        )

    lines!(ax, [extrema(t)...],  [η_fitted, η_fitted], color = :orange, linestyle = :dash)

    y_true = mean_shift_rate_true[90,i]
    η_true = 0.04 / 50.0
    lines!(ax, t, y_true, color = :black)
    lines!(ax, [extrema(t)...],  [η_true, η_true], color = :black, linestyle = :dash)

    #yboth = vcat(y_fitted, y_true)
    #ylims!(ax, (minimum(yboth)*0.9, maximum(yboth)*1.1))

    push!(figs, fig1)
end

figs[10]

for j in 1:10
    r = ((j-1)*25+1):(25+(j-1)*25)
    px = plot(ps[r]..., layout = (5,5), size = (1200, 1200))
    savefig(px, "figures/simulated_shift_through_time$j.pdf")
end








recent_mrca_slope = [(x[end] - x[1])/h for ((h,i), x) in mean_shift_rate_true]
ph1 = histogram(recent_mrca_slope, bins = 150, xlabel = "slope from oldest and youngest")
plot!(ph1, [0.0, 0.0], [0.0, 200], linestyle = :dash, label = "zero")

mean_diffs = [Statistics.mean(diff(x)) for ((h,i), x) in mean_shift_rate_true]
xlab = L"\text{mean}(\frac{dN}{dt}(t_i) - \frac{dN}{dt}(t_{i+1})∀ i)"
ph2 = histogram(mean_diffs, bins = 150, xlabel = xlab)
plot!(ph2, [0.0, 0.0], [0.0, 200], linestyle = :dash, label = "zero")


fig = Figure(size = (400, 250));

ax1 = Axis(fig[1,1], 
        ylabel = L"\text{number of datasets}", 
        xgridvisible = false, 
        ygridvisible = false,
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelsize = 9,
        yticklabelsize = 9)
hist1 = CairoMakie.hist!(ax1, mean_diffs, bins = 150, label = L"\text{simulated trees}", color = :gray)
lines!(ax1, [0.0, 0.0], [0.0, 481], label = L"\text{zero line}", color = :red, linestyle = :dash)
axislegend(ax1; position = :lt)


ax2 = Axis(fig[1,2], 
        xgridvisible = false, 
        ygridvisible = false,
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelsize = 9,
        yticklabelsize = 9)
hist2 = CairoMakie.hist!(ax2, mean_diffs, bins = 500, label = L"\text{simulated trees}", color = :gray)
CairoMakie.xlims!(ax2, -0.000001, 0.000002)
lines!(ax2, [0.0, 0.0], [0.0, 481], label = L"\text{zero line}", color = :red, linestyle = :dash)

xlabel = Label(fig[2,1:2], L"\text{mean}(\frac{dN}{dt}(t_i) - \frac{dN}{dt}(t_{i+1})\text{ for all } i)")

fig


branch_index = 2

Ds, Fs = backwards_forwards_pass(models_true[90,1], datasets[90,1]; condition = []);

branch_index = datasets[90,1].po[end]


branch_index = 2339

fig5 = Figure()
ax1 = Axis(fig5[1,1], title = "category 1")
ax2 = Axis(fig5[1,2], title = "category 2")
ax3 = Axis(fig5[1,3], title = "category 3")

for (i, ax) in enumerate((ax1, ax2, ax3))
    x = Fs[branch_index].t
    y = [F[i,1] for F in Fs[branch_index].u]
    plot!(ax, x, y)
    lines!(ax, x, y, linestyle = :dash, label  = "E(t) preorder")

    x = Ds[branch_index].t
    y = [D[i,1] for D in Ds[branch_index].u]
    #plot!(ax, x, y, alpha = 0.5)
    lines!(ax, x, y, label = "E(t) postorder")


end

#axislegend(ax1, position = :rb)
axislegend(ax2, position = :rb)

linkaxes!(ax1, ax2, ax3)

fig5
