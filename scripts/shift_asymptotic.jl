using CairoMakie


fig = Figure(size = (700, 450))

axs = Axis[]

prefixes = ["a) ", "b) ", "c) "]
for i in 1:3
    ax = Axis(fig[1,(i*2-1):(2*i)],
        title = prefixes[i] * "$(i+4) rate classes",
        topspinevisible = false,
        rightspinevisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xticks = (0:8, ["" for _ in 0:8]),
        yticks = (0:8, ["" for _ in 0:8]),
    )
    hidexdecorations!(ax, ticks = false)
    hideydecorations!(ax, ticks = false)

    push!(axs, ax)
end
         

ax1, ax2, ax3 = axs



Z = rand(0:2, (10, 10));

X1 = [
    3 3 2 3 3
    3 3 2 3 3
    1 1 0 1 1
    3 3 2 3 3
    3 3 2 3 3
]

X2 = [
    3 3 2 3 3 3
    3 3 2 3 3 3
    1 1 0 1 1 1
    3 3 2 3 3 3
    3 3 2 3 3 3
    3 3 2 3 3 3
]

X3 = [
    3 3 2 3 3 3 3
    3 3 2 3 3 3 3
    1 1 0 1 1 1 1
    3 3 2 3 3 3 3
    3 3 2 3 3 3 3
    3 3 2 3 3 3 3
    3 3 2 3 3 3 3
]

#cgrad(:Blues, 4; categorical = true)
colors = [:white, :orange, :teal, :lightgray]
x = [1,2,3,4,5,6]
#hm = heatmap!(ax, [1,2,3,4,5], [1,2,3,4,5], X'; colormap=cgrad(:Blues, 4; categorical=true), label = "as")

hm = heatmap!(ax1, 1:5, 1:5, X1'; colormap=colors, label = "as")
hm = heatmap!(ax2, 1:6, 1:6, X2'; colormap=colors, label = "as")
hm = heatmap!(ax3, 1:7, 1:7, X3'; colormap=colors, label = "as")
#Colorbar(fig[1, 2], hm; ticks=0:2)
#colsize!(fig.layout, 1, Aspect(1, 1.0))
#resize_to_layout!(f)

xlabel = Label(fig[2,1:2], "speciation rate classes")
xlabel = Label(fig[2,3:4], "speciation rate classes")
xlabel = Label(fig[2,5:6], "speciation rate classes")
ylabel = Label(fig[1,0], "extinction rate classes", rotation = π/2)


linkaxes!(axs...)

elem_1 = MarkerElement(color = :blue, marker = 'π', markersize = 15,
        points = Point2f[(0.2, 0.2), (0.5, 0.8), (0.8, 0.2)])

elem_1 = [PolyElement(color = :white, strokecolor = :black, strokewidth = 1)]
elem_2 = [PolyElement(color = :orange, strokecolor = :black, strokewidth = 1)]
elem_3 = [PolyElement(color = :teal, strokecolor = :black, strokewidth = 1)]
elem_4 = [PolyElement(color = :lightgray, strokecolor = :black, strokewidth = 1)]

elem_5 = PolyElement(color = :green, strokecolor = :black, strokewidth = 2,
        points = Point2f[(0, 0), (1, 0), (0, 1)])

Legend(fig[3, 1:2],
    [elem_1, elem_2, elem_3, elem_4],
    ["no change", "change in speciation rate", "change in extinction rate", "change in speciation & extinction rate"],
    patchsize = (20, 20), rowgap = 5)

## new figure
using LaTeXStrings

ax4 = Axis(fig[3,4:6],
    xlabel = "number of rate classes",
    #ylabel = L"n_\text{rate}/n_\text{total}",
    ylabel = "n(rate) / n(total)",
    title = "d) ratio to total",
    xgridvisible = false,
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
)

x = collect(range(3, 25))

n_ex = [(xi -1) / (xi^2) for xi in x]
n_sp = [(xi -1) / (xi^2) for xi in x]
n_both = [(xi^2 - 2*xi -1) / (xi^2) for xi in x]

lines!(ax4, x, n_ex, color = :orange, linewidth = 4, label = L"(n-1)/(n^2-1)")
lines!(ax4, x, n_ex, color = :teal, linewidth = 4, linestyle = :dash, label = L"(n-1)/(n^2-1)")
lines!(ax4, x, n_both, color = :lightgray, linewidth = 4, label = L"(n^2-2n-1)/(n^2-1)")

axislegend(ax4, position = :rc)

colsize!(fig.layout, 0, Relative(0.05))
for i in 1:6
    colsize!(fig.layout, i, Relative(1.0/6))
end
rowsize!(fig.layout, 1, Relative(0.5))
rowsize!(fig.layout, 2, Relative(0.05))
rowsize!(fig.layout, 3, Relative(0.45))

fig

CairoMakie.resize_to_layout!(fig)

CairoMakie.save("figures/shift_ratio_asymptotic.pdf", fig)

