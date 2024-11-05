using Revise
using Pesto
#using Distributions
using CSV
using JLD2
using ProgressMeter
using Glob

scratch = "/sto/nfsscratch/grp_shoehna/empirical_shifts/"

fpaths = Glob.glob("data/simulations/constant_speciation/*.tre")


n_iters = length(fpaths)

io = open("output/prog_constant_speciation.jl", "w")

completed_jobs = [
                  split(Base.basename(x), ".")[1] for x in Glob.glob("output/simulations/constant_speciation/jld2/*.jld2", scratch)
 ]

prog = Progress(n_iters; desc = "Inference (constant speciation): ", output= io);

for fpath in fpaths

    phy = readtree(fpath)
    sampling_probability = 1.0
    data = SSEdata(phy, sampling_probability)

    name = split(Base.basename(fpath), ".")[1]

    if name in completed_jobs
        continue ## skip this job if already completed
    end
    
    optres, model, n_attempts = optimize_hyperparameters(data; n = 10, n_attempts = 100)
    ntip = length(data.tiplab)

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


    #bf = posterior_prior_shift_odds(model,data)
    #append!(bf, NaN)
    #rates[!,"shift_bf"] = bf
    #rates[!,"shift_bf_log"] = log10.(bf)

    ## save data
    # scratch space
    # /sto/nfsscratch/grp_shoehna/empirical_shifts/output/simulations/constant_speciation
    fpath = string(scratch, "output/simulations/constant_speciation/newick/", name, ".tre")
    writenewick(fpath, data, rates)

    fpath = string(scratch, "output/simulations/constant_speciation/rates/", name, ".csv")
    CSV.write(fpath, rates)

    fpath = string(scratch, "output/simulations/constant_speciation/jld2/", name, ".jld2")

    save(fpath, 
        "N", N,
        "lambda", λ,
        "ntip", ntip,
        "mu", μ,
        "etaml", ηml)
    next!(prog)
end
finish!(prog)
close(io)


