using Revise
using Pesto
#using Distributions
using CSV
using JLD2
using ProgressMeter
using Glob

scratch = "/sto/nfsscratch/grp_shoehna/empirical_shifts/"

fpaths = Glob.glob("data/simulations/grafts/backbone/*.tre")


n_iters = length(fpaths)

io = open("output/prog_backbone.txt", "w")

completed_jobs = [
                  split(Base.basename(x), ".")[1] for x in Glob.glob("output/simulations/grafts/backbone/jld2/*.jld2", scratch)
 ]

prog = Progress(n_iters; desc = "Inference (backbone): ", output= io);

for fpath in fpaths

    phy = readtree(fpath)
    sampling_probability = 1.0
    data = SSEdata(phy, sampling_probability)
    ntip = length(data.tiplab)
    tl = sum(data.branch_lengths)

    name = split(Base.basename(fpath), ".")[1]

    if name in completed_jobs
        continue ## skip this job if already completed
    end

    try
        optres, model, n_attempts = optimize_hyperparameters(data; n = 10, sd = 0.578, n_attempts = 20)

        λ = model.λ
        μ = model.μ
        ηml = model.η

        Ds, Fs = backwards_forwards_pass(model, data);
        Ss = ancestral_state_probabilities(data, Ds, Fs);

        rates = birth_death_shift(model, data);
        shift_bf = rates[1:end-1,:shift_bf]
        is_significant = findall(shift_bf .> 100.0)

        N = state_shifts(model, data, Ds, Fs);
        N = N[is_significant,:,:];

        ## calculate S root
        root_index = length(data.tiplab)+1 
        root_children = findall(data.edges[:,1] .== root_index)
        left, right = root_children
        D_root = Ds[left].u[end] .* Ds[left].u[end] .* λ
        S_root = D_root ./ sum(D_root)

        netdiv_root = sum((λ .- μ) .* S_root)
        netdiv_tips = tip_rates(model, data, Ds, Fs)[!,:netdiv]

        ## save data
        fpath = string(scratch, "output/simulations/grafts/backbone/newick/", name, ".tre")
        writenewick(fpath, data, rates)

        fpath = string(scratch, "output/simulations/grafts/backbone/rates/", name, ".csv")
        CSV.write(fpath, rates)

        fpath = string(scratch, "output/simulations/grafts/backbone/jld2/", name, ".jld2")

        save(fpath, 
            "N", N,
            "lambda", λ,
            "ntip", ntip,
            "mu", μ,
            "etaml", ηml,
            "treelength", tl,
            "netdiv_root", netdiv_root,
            "netdiv_tips", netdiv_tips,
           )

    catch e
        if isa(e, Pesto.ConvergenceException)
            continue
        else
            rethrow(e)
        end
    end

    next!(prog)
end
finish!(prog)
close(io)


