using Distributions
using DataFrames, CSV, RCall, ProgressMeter
using LaTeXStrings, Measures
using JLD2
using Printf
using CairoMakie
using Glob


df = CSV.read("output/empirical_munged.csv", DataFrame)
metadata = CSV.read("data/empirical/metadata.csv", DataFrame)

#df = df[df[!,:inference] .== "empirical",:]
df[!,:type] = String.(df[!,:type]) ## ensure it's a string
df2 = df[df[!,:type] .== "strong support",:]

skip = [
    "Caenophidia_Zaher2019",
    "Charadriiformes_Cerny2022",
    "Cichlidae_McGee2020",
    "Ericaceae_Schwery2015",
    "Hyloidea_Hutter2017",
    "Liolaemidae_Esquerre2018",
    "Melastomataceae_Reginato2020",
    "Mimosa_Vasconcelos2020",
    "Papilionidae_Allio2021",
    #"Poaceae_Spriggs2015",
    "Monocots_Howard2019",
    "Polemoniaceae_Landis2018",
    "Primates_Springer2012",
    "Sigmodontinae_VallejosGarrido2023",
    "Viperidae_Alencar2016",
] 

keep = [!(name in skip) for name in df2[!,:name]]
df2 = df2[keep,:]

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
#kingdom["Squamata_Zheng2016"] = "Animalia"
kingdom["Squamata_Title2024"] = "Animalia"
kingdom["Coenagrionoidea_Willink2024"] = "Animalia"
kingdom["Anisoptera_Letsch2016"] = "Animalia"

shift_df = DataFrame(
    "log_height" => log10.(df2[!,:height]),
    "height" => df2[!,:height],
    "log_eta" => log10.(df2[!,:etaml]),
    "log_eta_times_tl" => log10.(df2[!,:etaml] .* df2[!,:treelength]),
    "log_N" => log10.(df2[!,:N_total]),
    "log_N_by_t" => log10.(df2[!,:N_per_time]),
    "log_netdiv" => log10.(df2[!,:tree_netdiv]),
    "log_mu" => log10.(df2[!,:tree_mu]),
    "log_lambda" => log10.(df2[!,:tree_lambda]),
    "number_of_supported" => df2[!,:number_of_supported],
    "number_of_supported_per_time" => df2[!,:number_of_supported] ./ df2[!,:treelength],
    "log_number_of_supported" => log10.(df2[!,:number_of_supported]),
    "log_number_of_supported_per_time" => log10.(df2[!,:number_of_supported] ./ df2[!,:treelength]),
    "kingdom" => [kingdom[name] for name in df2[!,:name]],
)


##################################
##
## some helper functions for the plots
##
##################################

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

## lrange
function lrange3(from::Float64, to::Float64, length::Int64 = 6)
    exp.(collect(range(log(from), log(to); length = length)))
end

function plot_panel!(ax, xvar, yvar, ymin, ymax, color)
    ## band
    β, Varβ, ySE = ols_regression(xvar, yvar)
    linefit(x) = β[1] + β[2]*x
    #x = collect(lrange3(extrema(xvar)..., 20))
    x = collect(range(extrema(xvar)..., 20))
    y = linefit.(x)

    yupper = y .+ 2 .* ySE.(x)
    ylower = y .- 2 .* ySE.(x)

    #ϵ = 0.7 * minimum(yvar)
    #ϵ = minimum(log10.((10 .^ yvar) .* 0.7))
    ϵ = minimum(yvar) - 0.15
    ϵ = ymin - 0.15
    ylower[ylower .< ϵ] .= ϵ

    α = maximum(yvar) + 0.15
    α = ymax + 0.15
    yupper[yupper .> α] .= α

    CairoMakie.band!(ax, x, ylower, yupper, color = "#e0e0e0")


    ## regression line
    y = linefit.(x)
    CairoMakie.lines!(ax, x, y; 
                label = "OLS", #markersize = 7, 
                color = "gray", linestyle = :dash)

    ## the points            
    scatter!(ax, xvar, yvar, color = color)

    nothing
end

xvar = log10.(shift_df[!,:height])
yvar = log10.(shift_df[!,:number_of_supported])


function make_tick_labels(x, nticks, ndigits)
    ticks = collect(range(extrema(x)...; length = nticks)) 
    fmt = Printf.Format("%.$(ndigits)f")
    ticklabels = [Printf.format(fmt, 10^tick) for tick in ticks]

    return (ticks, ticklabels)
end

############################################# 
##
##  Figure 3 -- with clade age
##
############################################# 

fig1 = Figure(size= (520, 410), figure_padding = 0)

colormap = Dict("Animalia" => :black, "Plantae" => :green, "Fungi" => :orange)
color = [colormap[k] for k in shift_df.kingdom]

###
xt, xtl = make_tick_labels(shift_df[!,:log_height], 5, 1)
yt = log10.([1, 10, 100, 1000])
ytl = ["1", "10", "100", "1000"]
xtl = ["" for _ in 1:5]
ax1 = Axis(fig1[2,2], 
    title = L"\text{a)}",
    titlealign = :left,
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    xticks = (xt, xtl),
    yticks = (yt, ytl),
    xticklabelrotation = π/2,
    )

ymin = minimum([minimum(shift_df[!,:log_number_of_supported]), minimum(shift_df[!,:log_eta_times_tl])])
ymax = maximum([maximum(shift_df[!,:log_number_of_supported]), maximum(shift_df[!,:log_eta_times_tl])])


plot_panel!(ax1, shift_df[!,:log_height], shift_df[!,:log_number_of_supported], ymin, ymax, color)

##
xt, xtl = make_tick_labels(shift_df[!,:log_height], 5, 2)
xtl = ["" for _ in 1:5]
ax2 = Axis(fig1[2,3], 
    title = L"\text{b)}",
    titlealign = :left,
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    xticks = (xt, xtl),
    yticks = (yt, ytl),
    xticklabelrotation = π/2,
    )

    
#ylims!(ax1, (log10(0.7), log10(1800)))
#ylims!(ax2, (log10(0.7), log10(1800)))

plot_panel!(ax2, shift_df[!,:log_height], shift_df[!,:log_eta_times_tl], ymin, ymax, color)

yt2, ytl2 = make_tick_labels(
    vcat(shift_df[!,:log_number_of_supported_per_time],
         shift_df[!,:log_eta]),
    5, 4)

##
ax3 = Axis(fig1[3,2], 
    xlabel = L"\text{clade age (Ma)}",
    #title = L"\text{c) slope}=-0.3\pm0.1,\text{ R}^2=33%",
    title = L"\text{c)}",
    titlealign = :left,
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    xticks = make_tick_labels(shift_df[!,:log_height], 5, 1),
    yticks = (yt2, ytl2),
    xticklabelrotation = π/2,
    )



fig1




ymin = minimum([minimum(shift_df[!,:log_number_of_supported_per_time]), minimum(shift_df[!,:log_eta])])
ymax = maximum([maximum(shift_df[!,:log_number_of_supported_per_time]), maximum(shift_df[!,:log_eta])])

plot_panel!(ax3, shift_df[!,:log_height], shift_df[!,:log_number_of_supported_per_time], ymin, ymax, color)

#
ax4 = Axis(fig1[3,3], 
    xlabel = L"\text{clade age (Ma)}",
    title = L"\text{d)}",
    titlealign = :left,
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    xticks = make_tick_labels(shift_df[!,:log_height], 5, 1),
    yticks = (yt2, ytl2),
    xticklabelrotation = π/2,
    )

plot_panel!(ax4, shift_df[!,:log_height], shift_df[!,:log_eta], ymin, ymax, color)

colgap!(fig1.layout, 12)
rowgap!(fig1.layout, 8)

linkxaxes!(ax1, ax2, ax3, ax4)
linkyaxes!(ax1, ax2)
linkyaxes!(ax3, ax4)

for ax in (ax2, ax4)
    hideydecorations!(ax,ticks = false)
end


left_title = Label(fig1[1,2], L"\text{supported branches }(N^*)")
right_title = Label(fig1[1,3], L"\text{all branches }(\eta)")

row1 = Label(fig1[2,1], L"\text{rate shift events}", rotation = π/2)
row2 = Label(fig1[3,1], L"\text{rate shift events/myr}", rotation = π/2)

Legend(fig1[2:3,4], 
    [
        MarkerElement(color = :black, marker = :circle, markersize = 9),
        MarkerElement(color = :green, marker = :circle, markersize = 9),
        MarkerElement(color = :orange, marker = :circle, markersize = 9),
    ],
    ["Animalia", "Plantae", "Fungi"],
    patchsize = (10, 10), 
    rowgap = 5,
)

for i in 2:3
    rowsize!(fig1.layout, i, Relative(0.35))
    colsize!(fig1.layout, i, Relative(0.35))
end

rowsize!(fig1.layout, 1, Relative(0.05))
colsize!(fig1.layout, 1, Relative(0.05))
colsize!(fig1.layout, 4, Relative(0.21))

## add coefficients
positions = [
    log10.((12.5, 500)),
    log10.((12.5, 1.0)),
    log10.((12.5, 0.0425)),
    log10.((12.5, 0.0001)),
]

labels = [
    L"\text{}s=0.65\pm0.41,\text{ R}^2=10%",
    L"\text{}s=0.67\pm0.29,\text{ R}^2=18%",
    L"\text{}s=-0.79\pm0.30,\text{ R}^2=22%",
    L"\text{}s=-0.77\pm0.15,\text{ R}^2=52%",
]

for (ax, pos, label) in zip([ax1, ax2, ax3, ax4], positions, labels)
    text!(ax, pos[1], pos[2],
        text = label, 
        align = (:left, :bottom),
        fontsize = 9,
    )
end

fig1


CairoMakie.save("figures/scatter_clade_age.pdf", fig1)

############################################# 
##
##  Figure 4 -- with net-div rates
##
############################################# 

fig2 = Figure(size= (520, 410), figure_padding = 0)

###
xt, xtl = make_tick_labels(shift_df[!,:log_netdiv], 5, 1)
yt = log10.([1, 10, 100, 1000])
ytl = ["1", "10", "100", "1000"]
xtl = ["" for _ in 1:5]
ax1 = Axis(fig2[2,2], 
    title = L"\text{a)}",
    titlealign = :left,
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    xticks = (xt, xtl),
    yticks = (yt, ytl),
    xticklabelrotation = π/2,
    )

ymin = minimum([minimum(shift_df[!,:log_number_of_supported]), minimum(shift_df[!,:log_eta_times_tl])])
ymax = maximum([maximum(shift_df[!,:log_number_of_supported]), maximum(shift_df[!,:log_eta_times_tl])])

plot_panel!(ax1, shift_df[!,:log_netdiv], shift_df[!,:log_number_of_supported], ymin, ymax, color)

##
xt, xtl = make_tick_labels(shift_df[!,:log_netdiv], 5, 2)
xtl = ["" for _ in 1:5]
ax2 = Axis(fig2[2,3], 
    title = L"\text{b)}",
    titlealign = :left,
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    xticks = (xt, xtl),
    yticks = (yt, ytl),
    xticklabelrotation = π/2,
    )

    

plot_panel!(ax2, shift_df[!,:log_netdiv], shift_df[!,:log_eta_times_tl], ymin, ymax, color)

yt2, ytl2 = make_tick_labels(
    vcat(shift_df[!,:log_number_of_supported_per_time],
         shift_df[!,:log_eta]),
    5, 4)

##
ax3 = Axis(fig2[3,2], 
    xlabel = L"\text{netdiv }(\lambda - \mu)",
    title = L"\text{c)}",
    titlealign = :left,
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    xticks = make_tick_labels(shift_df[!,:log_netdiv], 5, 2),
    yticks = (yt2, ytl2),
    xticklabelrotation = π/2,
    )


ymin = minimum([minimum(shift_df[!,:log_number_of_supported_per_time]), minimum(shift_df[!,:log_eta])])
ymax = maximum([maximum(shift_df[!,:log_number_of_supported_per_time]), maximum(shift_df[!,:log_eta])])
    
plot_panel!(ax3, shift_df[!,:log_netdiv], shift_df[!,:log_number_of_supported_per_time], ymin, ymax, color)

#
ax4 = Axis(fig2[3,3], 
    xlabel = L"\text{netdiv }(\lambda - \mu)",
    title = L"\text{d)}",
    titlealign = :left,
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    xticks = make_tick_labels(shift_df[!,:log_netdiv], 5, 2),
    yticks = (yt2, ytl2),
    xticklabelrotation = π/2,
    )

plot_panel!(ax4, shift_df[!,:log_netdiv], shift_df[!,:log_eta], ymin, ymax, color)

colgap!(fig2.layout, 12)
rowgap!(fig2.layout, 8)

linkxaxes!(ax1, ax2, ax3, ax4)
linkyaxes!(ax1, ax2)
linkyaxes!(ax3, ax4)

for ax in (ax2, ax4)
    hideydecorations!(ax,ticks = false)
end


left_title = Label(fig2[1,2], L"\text{supported branches }(N^*)")
right_title = Label(fig2[1,3], L"\text{all branches }(\eta)")

row1 = Label(fig2[2,1], L"\text{rate shift events}", rotation = π/2)
row2 = Label(fig2[3,1], L"\text{rate shift events/myr}", rotation = π/2)

Legend(fig2[2:3,4], 
    [
        MarkerElement(color = :black, marker = :circle, markersize = 9),
        MarkerElement(color = :green, marker = :circle, markersize = 9),
        MarkerElement(color = :orange, marker = :circle, markersize = 9),
    ],
    ["Animalia", "Plantae", "Fungi"],
    patchsize = (10, 10), 
    rowgap = 5,
)

for i in 2:3
    rowsize!(fig2.layout, i, Relative(0.35))
    colsize!(fig2.layout, i, Relative(0.35))
end

rowsize!(fig2.layout, 1, Relative(0.05))
colsize!(fig2.layout, 1, Relative(0.05))
colsize!(fig2.layout, 4, Relative(0.21))

## add coefficients
positions = [
    log10.((0.03, 570)),
    log10.((0.03, 1.0)),
    log10.((0.03, 0.0225)),
    log10.((0.03, 0.0001)),
]

labels = [
    L"\text{}s=0.35\pm0.47,\text{ R}^2=3%",
    L"\text{}s=-0.28\pm0.35,\text{ R}^2=3%",
    L"\text{}s=1.35\pm0.26,\text{ R}^2=53%",
    L"\text{}s=0.72\pm0.19,\text{ R}^2=37%",
]

for (ax, pos, label) in zip([ax1, ax2, ax3, ax4], positions, labels)
    text!(ax, pos[1], pos[2],
        text = label, 
        align = (:left, :bottom),
        fontsize = 9,
    )
end

fig2


CairoMakie.save("figures/scatter_netdiv.pdf", fig2)

## 


using GLM

function aic(m; smallsample = true)
    n, k = size(m.mm.m)
    res = 2*k - 2*GLM.loglikelihood(m)
    if smallsample
        res += (2*k^2 + 2*k)/(n - k -1)
    end
    return(res)
end

m0 = lm(@formula(log_eta ~ log_height), shift_df)
m1 = lm(@formula(log_eta_times_tl ~ log_height), shift_df)
#m2 = lm(@formula(log_eta ~ log_height + log_netdiv), shift_df)
m2 = lm(@formula(log_number_of_supported ~ log_height), shift_df)
m3 = lm(@formula(log_number_of_supported_per_time ~ log_height), shift_df)

aics = map(aic, [m0, m1, m2])
delta_aics = [x - minimum(aics) for x in aics]

m0a = lm(@formula(log_eta ~ log_netdiv), shift_df)
m1a = lm(@formula(log_eta_times_tl ~ log_netdiv), shift_df)
m2a = lm(@formula(log_number_of_supported ~ log_netdiv), shift_df)
m3a = lm(@formula(log_number_of_supported_per_time ~ log_netdiv), shift_df)



lm(@formula(log10(number_of_supported_per_time) ~ log_height), shift_df)
lm(@formula(log_eta ~ log_height), shift_df)


