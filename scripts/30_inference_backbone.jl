using Revise
using Pesto
#using Distributions
using CSV
using JLD2
using ProgressMeter
using Glob

scratch = "/sto/nfsscratch/grp_shoehna/empirical_shifts/"

fpaths = Glob.glob("data/simulations/single_shift_grafts/backbone/*.tre")


n_iters = length(fpaths)

io = open("output/prog_backbone.txt", "w")

completed_jobs = [
                  split(Base.basename(x), ".")[1] for x in Glob.glob("output/simulations/single_shift_grafts/backbone/*.jld2", scratch)
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

        rates = tree_rates(data, model, Fs, Ss);
        N = state_shifts(model, data, Ds, Fs);

        nshift = sum(N, dims = (2,3))[:,1,1];
        append!(nshift, 0.0)
        rates[!,"nshift"] = nshift

        Nsum = sum(N, dims = 1)[1,:,:]

        mag = magnitude(model, data, N);

        ## save data
        fpath = string(scratch, "output/simulations/single_shift_grafts/backbone/newick/", name, ".tre")
        writenewick(fpath, data, rates)

        fpath = string(scratch, "output/simulations/single_shift_grafts/backbone/rates/", name, ".csv")
        CSV.write(fpath, rates)

        fpath = string(scratch, "output/simulations/single_shift_grafts/backbone/jld2/", name, ".jld2")

        save(fpath, 
            "N", Nsum,
            "lambda", λ,
            "ntip", ntip,
            "mu", μ,
            "etaml", ηml,
            "magnitude", mag,
            "treelength", tl,
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


