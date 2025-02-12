using Distributions
using DataFrames, CSV, RCall, ProgressMeter
using LaTeXStrings, Measures
using JLD2
using Printf
using CairoMakie
using Glob


df = CSV.read("output/empirical_munged.csv", DataFrame)
inference = "empirical"

df = df[df[!,:inference] .== inference,:]

df[!,:type] = String.(df[!,:type]) ## ensure it's a string

df2 = df[df[!,:type] .== "strong support",:]


heights = df2[!,:height]


shift_df = DataFrame(
    "log_height" => log10.(df2[!,:height]),
    "log_eta" => log10.(df2[!,:etaml]),
    "log_eta_times_tl" => log10.(df2[!,:etaml] .* df2[!,:treelength]),
    "log_N" => log10.(df2[!,:N_total]),
    "log_N_by_t" => log10.(df2[!,:N_per_time]),
    #"log_netdiv" => log10.(df2[!,:lambdaml] .- df2[!,:muml])
    "log_netdiv" => log10.(df2[!,:tree_netdiv]),
    "log_mu" => log10.(df2[!,:tree_mu]),
    "log_lambda" => log10.(df2[!,:tree_lambda])
)


##################################
##
## now do it in Makie
## very annoying to do but a lot of possibility for finer control
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

function foobar!(fig, xvars, xlabs, shift_df, xdigits = [2,1], ydigits = [1,4])
    xt = [
        10 .^ (collect(range(extrema(xvars[i])...; length = 5))) for i in eachindex(xvars)
    ]

    offset = length(xvars)

    xtl = []
    for i in eachindex(xvars)
        fmt = Printf.Format("%.$(xdigits[i])f")
        s = [Printf.format(fmt, x) for x in xt[i]]
        append!(xtl, [s])
    end
    yt = 10 .^ (collect(range(extrema(shift_df[!,:log_eta_times_tl])...; length = 5)))
    
    fmt = Printf.Format("%.$(ydigits[1])f")
    ytl = [Printf.format(fmt, y) for y in yt]

    yt2 = 10 .^ (collect(range(extrema(shift_df[!,:log_eta])...; length = 5)))
    fmt = Printf.Format("%.$(ydigits[2])f")
    ytl2 = [Printf.format(fmt, y) for y in yt2]
    

    axs = []
    for i in eachindex(xvars)
        ax = Axis(fig[1, i], ylabel = L"\text{shifts }(\eta \times t)", xlabel = xlabs[i],
        xgridvisible = false, ygridvisible = false,
        yscale = log10, xscale = log10,
        xticks = (xt[i], xtl[i]),
        yticks = (yt, ytl),
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelrotation = π/2,
        xticklabelsize = 9,
        yticklabelsize = 9)

        CairoMakie.ylims!(ax, minimum(yt)*0.7, maximum(yt)*1.3)
        append!(axs, [ax])
    end

    for i in eachindex(xvars)
        ax = Axis(fig[2, i], ylabel = L"\text{shift rate }(\eta)", 
        xlabel = xlabs[i],
        xgridvisible = false, ygridvisible = false,
        yscale = log10, xscale = log10,
        xticks = (xt[i], xtl[i]),
        yticks = (yt2, ytl2),
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelrotation = π/2,
        xticklabelsize = 9,
        yticklabelsize = 9)
        CairoMakie.ylims!(ax, minimum(yt2)*0.7, maximum(yt2)*1.3)
        append!(axs, [ax])
    end

    ylowers = zeros(5, 20)
    yuppers = zeros(5, 20)
    for i in eachindex(xvars)
        β, Varβ, ySE = ols_regression(xvars[i], shift_df[!,:log_eta_times_tl])
        linefit(x) = β[1] + β[2]*x
        x = collect(lrange3(extrema(10 .^ (xvars[i]))..., 20))
        y = linefit.(log10.(x))

        yupper = 10 .^ (y .+ 2 .* ySE.(log10.(x)))
        ylower = 10 .^ (y .- 2 .* ySE.(log10.(x)))

        #ϵ = 1e-6
        ϵ = 0.7 * minimum(shift_df[!,:log_N])
        ylower[ylower .< ϵ] .= ϵ

        CairoMakie.band!(axs[i], x, ylower, yupper, color = "#e0e0e0")
        CairoMakie.scatter!(axs[i], 10 .^ xvars[i], 10 .^ shift_df[!,:log_eta_times_tl]; 
                        label = "asd", markersize = 7, color = "black")

        x = [extrema(10 .^ (xvars[i]))...]
        y = 10 .^(linefit.(log10.(x)))
        CairoMakie.lines!(axs[i], x, y; 
                    label = "OLS", markersize = 7, color = "gray", linestyle = :dash)

    end

    for i in eachindex(xvars)
        β, Varβ, ySE = ols_regression(xvars[i], shift_df[!,:log_eta])
        linefit(x) = β[1] + β[2]*x

        x = collect(lrange3(extrema(10 .^ xvars[i])..., 20))
        y = linefit.(log10.(x))
        yupper = 10 .^ (y .+ 2 .* ySE.(log10.(x)))
        ylower = 10 .^ (y .- 2 .* ySE.(log10.(x)))
        #ϵ = 1e-6
        ϵ = 0.7 * minimum(10 .^ shift_df[!,:log_eta])
        ylower[ylower .< ϵ] .= ϵ
        CairoMakie.band!(axs[i+offset], x, ylower, yupper, color = "#e0e0e0")
        CairoMakie.scatter!(axs[i+offset], 10 .^ xvars[i], 10 .^ shift_df[!,:log_eta]; 
                            label = "asd", markersize = 7, color = "black")
        x = [extrema(10 .^ xvars[i])...]
        y = 10 .^(linefit.(log10.(x)))
        CairoMakie.lines!(axs[i+offset], x, y; 
                    label = "OLS", markersize = 7, color = "gray", linestyle = :dash)

        intercept = @sprintf "%.2f" β[1]
        slope = @sprintf "%.2f" β[2]
        tag = string(
            "log(y) = ",
            intercept,
            " + ",
            slope,
            " log(x)"
            )

        if false ## print the intercept and slope
            tag = replace(tag, "+ -"=> "- ")
            tag = replace(tag, "-"=> "- ")
            tag = replace(tag, "log"=> "\\log")
            tag = LaTeXString(string("\$", tag, "\$"))
            text!(axs[i+offset], minimum(xt[i])*1.1, maximum(yt2)*0.65, text = tag, fontsize = 7)
        end
    end


    for i in 2:offset
        hideydecorations!(axs[i], ticks = false)
        hideydecorations!(axs[i+offset], ticks = false)
    end
    for i in 1:offset
        hidexdecorations!(axs[i], ticks = false)
    end

    colgap!(fig.layout, 8)
    rowgap!(fig.layout, 8)

    return()
end



fig1 = Figure(size=(350, 300), fontsize = 14);
xlabs = [
    #L"\lambda",
    #L"\mu",
    "", #L"\lambda - \mu",
    "", ##L"\text{tree height (Ma)}",
]
xvars = [
    #shift_df[!,:log_lambda],
    #shift_df[!,:log_mu],
    shift_df[!,:log_netdiv],
    shift_df[!,:log_height],
]
xdigits = [2, 1]
ydigits = [1, 4]

foobar!(fig1, xvars, xlabs, shift_df, xdigits, ydigits)
xlabel1 = Label(fig1[3,1], L"\lambda - \mu")
xlabel2 = Label(fig1[3,2], L"\text{tree height (Ma)}")
fig1

rowgap!(fig1.layout, 9)

for i in 1:2
    colsize!(fig1.layout, i, Relative(0.5))
    rowsize!(fig1.layout, i, Relative(0.495))
end
rowsize!(fig1.layout, 3, Relative(0.01))
fig1


fig1

#CairoMakie.save("figures/scatter1.pdf", fig1)
CairoMakie.save("figures/scatter1_empirical.pdf", fig1)

β, Varβ, ySE = ols_regression(xvars[2], shift_df[!,:log_eta])

print(β[2], " ± ", sqrt(Varβ[2,2]))

######################
######################
######################
######################


using GLM

function aic(m; smallsample = true)
    k = size(m.mm.m)[2]
    res = 2*k - 2*GLM.loglikelihood(m)
    if smallsample
        res += (2*k^2 + 2*k)/(n )
    end
    return(res)
end

m0 = lm(@formula(log_N_by_t ~ log_height), shift_df)
m1 = lm(@formula(log_N_by_t ~ log_netdiv), shift_df)
m2 = lm(@formula(log_N_by_t ~ log_height + log_netdiv), shift_df)

aics = map(aic, [m0, m1, m2])
delta_aics = [x - minimum(aics) for x in aics]

s = summary(m0)
lm(@formula(log_eta ~ log_height), shift_df)
lm(@formula(log_eta ~ log_netdiv), shift_df)
lm(@formula(log_eta ~ log_height + log_netdiv), shift_df)


