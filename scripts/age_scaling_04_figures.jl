using Makie
using CSV
using DataFrames
using CairoMakie
using LaTeXStrings
using Printf


function ols_regression(x, y)
    X = hcat([1 for _ in 1:length(x)], x)
    n, p = size(X)

    ## OLS
    β = (X' * X) \ X' * y
    yhat = X * β
    #sigma_squared = (1 / (n - p - 1)) * sum((y .- yhat).^2) ## MLE for sigma^2
    s_squared = (y .- yhat)' * (y .- yhat) ./ (n - p) ## OLS for sigma^2
    Varβ = inv(X' * X) .* s_squared
    yVar = x -> Varβ[1,1] + (x^2)*Varβ[2,2] + 2*x*Varβ[1,2]
    ySE = x -> sqrt.(yVar(x))
    return(β, Varβ, ySE)
end

ms = ["NA", "NAN", "NULL", "NaN"]
df = CSV.read("output/age_scaling_effect_munged.csv", DataFrame; missingstring = ms)

#df = df[df[!,:type] .== "pooled",:]
df = df[df[!,:type] .== "pooled",:]
#df = df[df[!,:inference] .== "age_scaling_effect",:]
df = df[df[!,:inference] .== "age_scaling_effect",:]

## filter for at least one strongly supported shift
#df = df[df[!,:how_many_supported ] .> 0,:]

df[:,["N_total", "N_per_time", "how_many_supported", "support_per_time"]]

df_true = CSV.read("output/age_scaling_effect_true_nshift.csv", DataFrame)
df_true[!,:height] = Int64.(round.(df_true[!,:height]))


fig = Makie.Figure(size = (450, 220), fontsize = 14, figure_padding = (1,1,1,1))

xt = range(30, 100; length = 8) |> collect

kwargs = (; 
topspinevisible = false,
rightspinevisible = false,
xgridvisible = false,
ygridvisible = false,
titlealign = :center,
xticks = (xt, [@sprintf("%.0f", x) for x in xt])
)

ax1 = Axis(fig[1,1]; 
    ylabel = L"N/t",
    xlabel = "tree height (Ma)",
    #title = "a) all pooled (n = $a)",       
    #ylabel = L"N_\text{rate} / N_\text{all}",
    kwargs...)

#xs = 30

colors = [:steelblue, "orange", "gray"]
#cs = vcat([repeat([i], size(df)[1]) for i in colors]...)

rainclouds!(ax1, 
    df[!,:height] .- 2,
    df[!,:N_per_time]; clouds = nothing,
    plot_boxplots = true,
    #color = cs,
    boxplot_width = 4,
    side_nudge = 4,
    cloud_width = 4,
    markersize = 4,
    label = "label",
    jitter_width = 4)
fig


fig2 = Makie.Figure(size = (450, 450), fontsize = 14, figure_padding = (1,1,1,1))

n = size(df)[1]
ax1 = Axis(fig2[1,1]; 
    ylabel = L"\hat{N}/t \text{ (all branches)}",
    #xlabel = "tree height (Ma)",
    title = L"\text{simulated trees, filtered for }\geq 1\text{ supported branch }(n = %$n)",
    titlealign = :center,
    kwargs...)

ax2 = Axis(fig2[2,1]; 
    ylabel = L"\text{no. supported branches}/t",
    xlabel = L"\text{tree height (Ma)}",
    kwargs...)

for h in unique(df[!,:height])
    y = df[df[!,:height] .== h,:N_per_time]
    hist!(ax1, y, offset = h, scale_to = -5, direction = :x, bins = 50)
    boxplot!(ax1, repeat([h+2.5], length(y)), y, show_outliers = true, width = 5)

    y = df[df[!,:height] .== h,:support_per_time]
    hist!(ax2, y, offset = h, scale_to = -5, direction = :x, bins = 50)
    boxplot!(ax2, repeat([h+2.5], length(y)), y, show_outliers = true, width = 5)
    #boxplot!(repeat([h+2.5], 1000), y, show_outliers = true, width = 5)
end

#ylims!(ax2, -0.001, 0.03)
rowgap!(fig2.layout, 1.0)
fig2

save("figures/age_scaling_effect.pdf", fig2)


#### TODO
## need to replicate the scatter plot regression, 
## and calculate the slope for the simulated data


fig3 = Figure(size = (650, 300))

heights = [30, 40, 50, 60, 70, 80, 90, 100]

xt = heights
yt = [0.03, 0.02, 0.01, 0.005]

ax1 = Axis(
    fig3[1,1],
    xticks = xt,
    #yticks = yt,
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    yscale = log10,
    #yscale = Makie.Symlog10(-0.005, 0.005),
    xlabel = L"\text{tree height (Ma)}",
    ylabel = L"\text{no. rate shift events per time}",
)

ylims!(ax1, 1e-6*0.8, 0.3)

fig3


colors = [:black, :orange]
labels = ["estimated shifts", "true shifts"]

for (i, this_df) in enumerate([df, df_true])
    for height in heights
        df1 = filter(:height => h -> h == height, this_df)

        nrows = size(df1)[1]
        #y = df1[!,:N_per_time]
        if hasproperty(df1, :support_per_time)
            y = df1[!,:etaml]
        else
            y = df1[!,:N_per_time]
        end

        x = df1[!,:height]
        w = 3.0
        x_jitter = x .+ (rand(nrows).-0.5) .* w
        #x_jitter 
        if i == 1
            x_jitter .-= w/2.0
        else
            x_jitter .+= w/2.0
        end

        scatter!(ax1, x_jitter, y, color = colors[i], label = labels[i], markersize = 4)
    end
end

η_true = 0.0008
lines!(ax1, [25, 105], [η_true, η_true], color = :red, linestyle = :dash)

elem_1 = MarkerElement(color = :orange, marker = :circle, markersize = 15)
elem_2 = MarkerElement(color = :black, marker = :circle, markersize = 15)
elem_3 = LineElement(color = :red, linestyle = :dash)

Legend(fig3[1,2], 
    [elem_1, elem_2, elem_3],
    [
        L"\text{Realized shift rate }(N_\text{true}/t)",
        L"\text{Estimated shift rate }(\eta)", 
        L"\text{True shift rate }(\eta = 0.0008)",
        ],
    patchsize = (35, 35), rowgap = 5,
    framevisible = false
    )
fig3

save("figures/age_scaling_effect2.pdf", fig3)




## estimation error

df_true[!,:N_per_time_true] = df_true[!,:N_per_time]
df_true[!,:N_total_true] = df_true[!,:N]

dfx = innerjoin(df, df_true, on = [:tree_index,:height]; makeunique=true)

error = dfx[!,:N_per_time] .- dfx[!,:N_per_time_true]


fig4 = Figure(size = (500, 500))

heights = [30, 40, 50, 60, 70, 80, 90, 100]
xt = heights

ax2 = Axis(
    fig4[1,1],
    xticks = xt,
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    #yscale = log10,
    title = L"\text{a) error in shift rate}",
    #ylabel = L"\text{estimation error (}\eta - N_\text{true}/t)",
)
    
ax3 = Axis(
    fig4[2,1],
    xticks = xt,
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    #yscale = log10,
    title = L"\text{b) error in supported shift events}",
    xlabel = L"\text{tree height (Ma)}",
    #ylabel = L"\text{estimation error (}N^*/t - N_\text{true}/t)",
)

Label(fig4[1,0], L"\text{estimation error (}\eta - N_\text{true}/t)", rotation = π/2)
Label(fig4[2,0], L"\text{estimation error (}N^*/t - N_\text{true}/t)", rotation = π/2)


for (i, height) in enumerate(heights)
    this_df = filter(:height => h -> h == height, dfx)

    error = this_df[!,:etaml] .- this_df[!,:N_per_time_true]
    hist!(ax2, error, color = (:black, 0.5), direction = :x, scale_to = 8, offset = height, bins = 20, label = "simulated trees")

    error_normalized = this_df[!,:support_per_time] .- this_df[!,:N_per_time_true]
    hist!(ax3, error_normalized, color = (:black, 0.5), direction = :x, scale_to = 8, offset = height, bins = 20, label = "simulated trees")
end
for ax in (ax2, ax3)
    lines!(ax, [28, 108], [0.0, 0.0], linestyle = :dash, color = :red, label = "zero error")
end
axislegend(ax3, position = :rb, unique = true)
fig4

save("figures/age_scaling_effect_estimation_error.pdf", fig4)


this_df = filter(:height => h -> h == 100, dfx)
error = this_df[!,:N_total] .- this_df[!,:N_total_true]
#error = (this_df[!,:N_total] .- this_df[!,:N_total_true]) ./ this_df[!,:N_total_true]
hist(error)

y = this_df[!,:N_total]
x = this_df[!,:N_total_true]
scatter(x, y)

error


error = dfx[!,:N_per_time] .- dfx[!,:N_per_time_true]
save("figures/age_scaling_effect_estimation_error_ntip.pdf", fig5)



fig5 = Figure(size = (500, 500))

xt = [1, 10, 100, 1000, 10_000]

ax1 = Axis(
    fig5[1,1],
    xscale = log10,
    xticks = xt,
    #xlabel = L"\text{number of tips}",
    title = L"\text{a) error in shift rate}",
    xgridvisible = false,
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
    )

ax2 = Axis(
    fig5[2,1],
    xscale = log10,
    xticks = xt,
    xlabel = L"\text{number of tips}",
    title = L"\text{b) error in supported shift events}",
    xgridvisible = false,
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
    )

    
Label(fig5[1,0], L"\text{estimation error (}\eta - N_\text{true}/t)", rotation = π/2)
Label(fig5[2,0], L"\text{estimation error (}N^*/t - N_\text{true}/t)", rotation = π/2)

scatter!(
    ax1,
    dfx.ntip, 
    dfx.etaml .- dfx.N_per_time_true,
    color = (:black, 0.3))

scatter!(
    ax2,
    dfx.ntip, 
    dfx.support_per_time .- dfx.N_per_time_true,
    color = (:black, 0.3),
    label = "simulated trees")


for ax in (ax1, ax2)
    lines!(ax, [extrema(dfx.ntip)...], [0.0, 0.0], linestyle = :dash, color = :red, label = "zero error")
end

for i in 1:2
    rowsize!(fig5.layout, i, Relative(0.5))
end

axislegend(ax2, position = :rb, unique = true)

fig5
save("figures/age_scaling_effect_estimation_error_ntip.pdf", fig5)

