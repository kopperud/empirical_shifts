using Revise
using Pesto
using Glob
using JLD2


fpaths = glob("output/simulations/single_shift_grafts/backbone/jld2/*.jld2")

for (i, fpath) in enumerate(fpaths)
    phy = readtree("data/simulations/single_shift_grafts/backbone/$i.tre")
    tl = sum(phy.edge_length)

    x = load(fpath)
    x["treelength"] = tl
    save(fpath, x)
end

fpaths = glob("output/simulations/single_shift_grafts/upshift/jld2/*.jld2")

for (i, fpath) in enumerate(fpaths)
    phy = readtree("data/simulations/single_shift_grafts/upshift/$i.tre")
    tl = sum(phy.edge_length)

    x = load(fpath)
    x["treelength"] = tl
    save(fpath, x)
end


fpaths = glob("output/simulations/single_shift_grafts/downshift/jld2/*.jld2")

for (i, fpath) in enumerate(fpaths)
    phy = readtree("data/simulations/single_shift_grafts/downshift/$i.tre")
    tl = sum(phy.edge_length)

    x = load(fpath)
    x["treelength"] = tl
    save(fpath, x)
end


