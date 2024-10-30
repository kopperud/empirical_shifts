using Revise
using BirthDeathSimulation
using ProgressMeter
using JLD2
using Distributions
using Random

## simulate trees conditional on survival until time t
Random.seed!(1234)

## a model with moderate rate variation
r0 = 0.07

λ1 = 0.30
μ = [λ1 - r0 for _ in 1:3]
λ = [λ1 - 0.03, λ1, λ1 + 0.03]

η = r0 / 50

model = bdsmodel(λ, μ, η)

starting = [2]
tree_height = 70.0

## simulate trees conditional on survival until time t
Random.seed!(1234)

maxtaxa = 50000
n_trees = 500
n_states = length(model.λ)

##
trees = Array{Tree, 1}(undef, n_trees)
N = Array{Int64, }(undef, n_trees, n_states, n_states)

begin 
    prog = ProgressUnknown("Complete trees sim:")
    i = 1
    while i <= n_trees
        maxtime = tree_height
        #maxtime = tree_heights[j]
        state = 2
        tree = sim_bdshift(model, maxtime, maxtaxa, state)
        if length(tree.Leaves) < maxtaxa ## reject complete trees that termined when too many taxa
            if length(tree.Leaves) > 5 ## Reject completete trees where all taxa went extinct
                prune_extinct!(tree)
                if abs(treeheight(tree) - maxtime) < 0.001 ## reject complete trees where both root children did not survive
                    N0 = +([branch.N for (idx, branch) in tree.Branches]...)


                    if sum(N0) > 0 ## reject trees without any shifts
                        N[i,:,:] .= N0
                        trees[i] = tree
                        i += 1
                     end
                end
            end
        end
        ProgressMeter.next!(prog)
    end
    ProgressMeter.finish!(prog)
end


length(trees)

ntaxa = [length(tree.Leaves) for tree in trees]
hist(ntaxa)



prog = Progress(n_trees, "Writing trees:")
for j in 1:n_trees
    fpath = string(
        "data/simulations/constant_extinction/", string(j), ".tre"
        )
    writenewick(fpath, trees[j], model)
    next!(prog)
end

## do the same but for constant speciation
r0 = 0.07

λ1 = 0.30
μ1 = λ1 - r0

λ = [λ1 for _ in 1:3]
μ = [μ1 - 0.03, μ1, μ1 + 0.03]

η = r0 / 50

model2 = bdsmodel(λ, μ, η)



##
trees2 = Array{Tree, 1}(undef, n_trees)
N2 = Array{Int64, }(undef, n_trees, n_states, n_states)

begin 
    prog = ProgressUnknown("Complete trees sim:")
    i = 1
    while i <= n_trees
        maxtime = tree_height
        #maxtime = tree_heights[j]
        state = 2
        tree = sim_bdshift(model2, maxtime, maxtaxa, state)
        if length(tree.Leaves) < maxtaxa ## reject complete trees that termined when too many taxa
            if length(tree.Leaves) > 5 ## Reject completete trees where all taxa went extinct
                prune_extinct!(tree)
                if abs(treeheight(tree) - maxtime) < 0.001 ## reject complete trees where both root children did not survive
                    N0 = +([branch.N for (idx, branch) in tree.Branches]...)


                    if sum(N0) > 0 ## reject trees without any shifts
                        N2[i,:,:] .= N0
                        trees2[i] = tree
                        i += 1
                     end
                end
            end
        end
        ProgressMeter.next!(prog)
    end
    ProgressMeter.finish!(prog)
end



prog = Progress(n_trees, "Writing trees:")
for j in 1:n_trees
    fpath = string(
        "data/simulations/constant_speciation/", string(j), ".tre"
        )
    writenewick(fpath, trees2[j], model)
    next!(prog)
end


ntaxa2 = [length(tree.Leaves) for tree in trees2]
hist(ntaxa2)

sum(N2, dims = (2,3))[:,1,1] |> extrema