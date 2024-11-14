using Revise
using Pesto
#using Distributions
using CSV
using JLD2
using ProgressMeter
using Glob

scratch = "/sto/nfsscratch/grp_shoehna/empirical_shifts/"

fpaths = Glob.glob("data/simulations/single_shift_grafts/upshift/*.tre")

using Distributions

rhat = 0.1
dr = LogNormal(log(rhat), 0.587)
r = Pesto.make_quantiles(dr, 36)
rel_ext = 0.67 ## 

λ = r ./ (1 .- rel_ext)
μ = λ .- r


α = 0.05; β = 0.01
Q = zeros(36, 36); 
for i in 1:36
    for j in 1:36
        if i > j
            Q[i,j] = α
        elseif i < j
            Q[i,j] = β
        end
    end
end
for i in 1:36
    Q[i,i] = -sum(Q[:,i])
end

E = [0.1 for _ in 1:36]
dE1 = zeros(36)
for i in 1:36
    for j in 1:36
        dE[i] += Q[j,i] * E[j]
    end
end
dE2 = Q * E

dE1
dE2
#dE1 = deepcopy(dE)

dE2

model = Pesto.BDSconstantQ(λ, μ, Q)

phy = readtree("data/simulations/single_shift_grafts/upshift/5.tre")
data = SSEdata(phy, 1.0)

E = extinction_probability(model, data)
logL_root(model, data)


αs = collect(lrange(0.000001, 0.04, 20))
βs = collect(lrange(0.000001, 0.04, 21))

logLs = zeros(20, 21)
for (αidx, α) in enumerate(αs)
    for (βidx, β) in enumerate(βs)
        Q = zeros(36, 36); 
        for i in 1:36
            for j in 1:36
                if i > j
                    Q[i,j] = α / (36-1)
                elseif i < j
                    Q[i,j] = β / (36-1)
                end
            end
        end
        for i in 1:36
            Q[i,i] = -sum(Q[:,i])
        end
        model = Pesto.BDSconstantQ(λ, μ, Q)

        logLs[αidx,βidx] = logL_root(model, data) 
    end
end

Q

extrema(logLs)
fig = Figure();
ax = Axis3(fig[1,1], xlabel = "α: rate of down shifts", ylabel = "β: rate of up shifts", zlabel ="log lik")
surface!(ax, αs, βs, logLs)
fig

argmax(logLs)

a = hcat([αs for _ in 1:length(βs)]...)[argmax(logLs)]
b = transpose(hcat([βs for _ in 1:length(αs)]...))[argmax(logLs)]

println("α = $a")
println("β = $b")



plot(E)


using LoopVectorization
using BenchmarkTools

dE = zeros(36); E = zeros(36);
p = (36, λ, μ, Q)

@benchmark foo(dE, E, p, 0.0)

foo(dE, E, p, 0.0)

dE1 = deepcopy(dE)
dE2 = deepcopy(dE)

dE1 .- dE2

function foo(dE, E, p, t)
    K, λ, μ, Q = p

    sumE = sum(E)
    #dE[:] .= μ .- (λ .+ μ .+ Q[1,1]) .* E  .+ λ .* E .* E .+ (Q[1,1]/(K-1)) .* (sumE .- E) 
    
    
    LoopVectorization.@turbo warn_check_args=false for i in eachindex(E) 
        dE[i] = μ[i] - (λ[i] + μ[i]) * E[i] + λ[i] * E[i] * E[i]
        for j in eachindex(E)
            dE[i] += Q[i,j] * E[j]
        end
    end
    
    nothing
end



n_iters = length(fpaths)

io = open("output/prog_upshift.txt", "w")

completed_jobs = [
                  split(Base.basename(x), ".")[1] for x in Glob.glob("output/simulations/single_shift_grafts/upshift/*.jld2", scratch)
 ]

prog = Progress(n_iters; desc = "Inference (upshift): ", output= io);

for fpath in fpaths

    phy = readtree(fpath)
    sampling_probability = 1.0
    data = SSEdata(phy, sampling_probability)
    ntip = length(data.tiplab)

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
        fpath = string(scratch, "output/simulations/single_shift_grafts/upshift/newick/", name, ".tre")
        writenewick(fpath, data, rates)

        fpath = string(scratch, "output/simulations/single_shift_grafts/upshift/rates/", name, ".csv")
        CSV.write(fpath, rates)

        fpath = string(scratch, "output/simulations/single_shift_grafts/upshift/jld2/", name, ".jld2")

        save(fpath, 
            "N", Nsum,
            "lambda", λ,
            "ntip", ntip,
            "mu", μ,
            "etaml", ηml,
            "magnitude", mag)

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


