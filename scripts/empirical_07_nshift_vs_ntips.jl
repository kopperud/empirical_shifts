##################################
##
## now do it in Makie
## very annoying to do but a lot of possibility for finer control
##
##################################

using CSV
using CairoMakie

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

function lrange3(from, to, length::Int64 = 6)
    exp.(collect(range(log(from), log(to); length = length)))
end




#df = CSV.read("output/munged_magnitude.csv", DataFrame)
df = CSV.read("output/empirical_munged.csv", DataFrame)
filter!(x -> x[:type] == "strong support", df)

ntips = df[!,:ntips]
shifts_per_time = df[!,:number_of_supported] ./ df[!,:treelength]
shifts = df[!,:number_of_supported]
shifts_per_tips = df[!,:number_of_supported] ./ df[!,:ntips]
tips_per_shift = df[!,:ntips] ./ df[!,:number_of_supported] 

df[!,:number_of_supported_per_time] = df[!,:number_of_supported] ./ df[!,:treelength]

sort(df, :number_of_supported_per_time, rev = true)[!,:number_of_supported_per_time]


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


df[!,:kingdom] = [kingdom[name] for name in df[!,:name]]
colormap = Dict("Animalia" => :black, "Plantae" => :green, "Fungi" => :orange)
color = [colormap[k] for k in df.kingdom]


xt = collect(lrange3(extrema(Float64.(ntips))..., 5))
yt = collect(lrange3(extrema(shifts)..., 5))
yt2 = collect(lrange3(extrema(shifts_per_time)..., 5))
yt3 = collect(lrange3(extrema(tips_per_shift)..., 5))

fmt = Printf.Format("%.0f")
xtl = [Printf.format(fmt, x) for x in xt]

fmt = Printf.Format("%.0f")
ytl = [Printf.format(fmt, y) for y in yt]

fmt2 = Printf.Format("%.4f")
ytl2 = [Printf.format(fmt2, y) for y in yt2]

fmt3 = Printf.Format("%.1f")
ytl3 = [Printf.format(fmt3, y) for y in yt3]

fig = Figure(size = (650, 210), fontsize = 14);

## number of rate shifts (not per time)
ax1 = Axis(fig[1,1], 
            ylabel = L"N^*", 
            #xlabel = L"\text{number of tips}",
            title = L"\text{a) shifts}",
            xgridvisible = false, 
            ygridvisible = false,
            yscale = log10, xscale = log10,
            xticks = (xt, xtl),
            yticks = (yt, ytl),
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelrotation = π/2,
            xticklabelsize = 9,
            yticklabelsize = 9)

CairoMakie.ylims!(ax1, minimum(yt)*0.7, maximum(yt)*1.3)

β, Varβ, ySE = ols_regression(log10.(ntips), log10.(shifts)) 
linefit(x) = β[1] + β[2]*x
x = collect(lrange3(Float64.(extrema(ntips))..., 20))
y = linefit.(log10.(x))


yupper = 10 .^ (y .+ 2 .* ySE.(log10.(x)))
ylower = 10 .^ (y .- 2 .* ySE.(log10.(x)))
ϵ = 0.7 * minimum(log10.(shifts))
ylower[ylower .< ϵ] .= ϵ

CairoMakie.band!(ax1, x, ylower, yupper, color = "#e0e0e0")
CairoMakie.scatter!(ax1, ntips, shifts; 
                    label = "asd", markersize = 7, color = color)

x = [extrema(Float64.(ntips))...]
y = 10 .^(linefit.(log10.(x)))
CairoMakie.lines!(ax1, x, y; 
                label = "OLS", color = "gray", linestyle = :dash)


## number of rate shifts per time
ax2 = Axis(fig[1,2], 
            ylabel = L"N^*/t", 
            #xlabel = L"\text{number of tips}",
            title = L"\text{b) shifts per time}",
            xgridvisible = false, 
            ygridvisible = false,
            yscale = log10, xscale = log10,
            xticks = (xt, xtl),
            yticks = (yt2, ytl2),
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelrotation = π/2,
            xticklabelsize = 9,
            yticklabelsize = 9)

CairoMakie.ylims!(ax2, minimum(yt2)*0.7, maximum(yt2)*1.3)

β, Varβ, ySE = ols_regression(log10.(ntips), log10.(shifts_per_time)) 
linefit(x) = β[1] + β[2]*x
x = collect(lrange3(Float64.(extrema(ntips))..., 20))
y = linefit.(log10.(x))


yupper = 10 .^ (y .+ 2 .* ySE.(log10.(x)))
ylower = 10 .^ (y .- 2 .* ySE.(log10.(x)))
ϵ = 0.7 * minimum(log10.(shifts_per_time))
ylower[ylower .< ϵ] .= ϵ

CairoMakie.band!(ax2, x, ylower, yupper, color = "#e0e0e0")
CairoMakie.scatter!(ax2, ntips, shifts_per_time; 
                    label = "asd", markersize = 7, color = color)

x = [extrema(Float64.(ntips))...]
y = 10 .^(linefit.(log10.(x)))
CairoMakie.lines!(ax2, x, y; 
                label = "OLS", color = "gray", linestyle = :dash)

               


ax3 = Axis(fig[1,3], 
            yticks = (log10.(yt3), ytl3), 
            xgridvisible = false, 
            #xlabel = L"\text{frequency}",
            ygridvisible = false,
            topspinevisible = false,
            rightspinevisible = false,
            title = L"\text{c) tips per shift}",
            xticklabelrotation = π/2,
            xticklabelsize = 9,
            yticklabelsize = 9,
        ylabel = L"\text{no. tips}/N^*")

Label(fig[2,1], L"\text{number of tips}")
Label(fig[2,2], L"\text{number of tips}")
Label(fig[2,3], L"\text{frequency}")

hist!(ax3, log10.(tips_per_shift), direction = :x, color = :gray, bins = 20)

Legend(fig[1,4], 
    [
        MarkerElement(color = :black, marker = :circle, markersize = 9),
        MarkerElement(color = :green, marker = :circle, markersize = 9),
        MarkerElement(color = :orange, marker = :circle, markersize = 9),
    ],
    ["Animalia", "Plantae", "Fungi"],
    patchsize = (10, 10), 
    rowgap = 5,
)

rowsize!(fig.layout, 1, Relative(0.85))
rowsize!(fig.layout, 2, Relative(0.15))
rowgap!(fig.layout, 0)

for i in 1:3
    colsize!(fig.layout, i, Relative(0.27))
end
colsize!(fig.layout,4,Relative(0.19))

fig

CairoMakie.save("figures/shifts_vs_ntips.pdf", fig)

