using PhyloNetworks
using EcologicalNetworks
using DelimitedFiles
using Plots
using StatsBase
using Statistics
using SparseArrays
using EcologicalNetworksPlots
using CSV
using DataFrames
using GBIF
using Base.Threads
using StatsBase
using StatsModels
using Distributions

# Get ancillary functions
include("lib/pwar.jl")

# Read the tree
tree_net = readTopology(joinpath("data", "mammals.newick"));

# Read the names
namelist = DataFrame(CSV.File(joinpath("artifacts", "names_metaweb_tree_gbif.csv")))

# Prepare the metaweb
eurometa = readdlm(joinpath("artifacts", "europeanmetaweb.csv"), ',', String)
mwspecies = unique(eurometa)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(eurometa, 1)
    M[eurometa[i, :]...] = true
end

# Partition the metaweb into two latent subspaces
L, R = rdpg(M, 5)

treeleaves = tipLabels(tree_net)

canada = DataFrame(CSV.File(joinpath("artifacts", "iucn_gbif_names.csv")))
cancodes = replace.(unique(filter(!ismissing, canada.gbifname)), " " => "_")
tree_cleanup = DataFrame(CSV.File(joinpath("artifacts", "upham_gbif_names.csv")))

csp = dropmissing!(DataFrame(; gbifname=canada.gbifname))
csp = dropmissing(leftjoin(csp, tree_cleanup; on=:gbifname))

canmammals = unique(csp.code)

## Here I diverge from Tim, who does k=3 phylogenetic kmeans

pool = DataFrame(; tipNames=canmammals)

metawab_usnames = replace.(M.S, " " => "_")
metawab_usnames ∩ canmammals

#### I am missing some species from the metaweb
filter((x) -> (x) ∉ (metawab_usnames ∩ treeleaves), metawab_usnames)

traitframe = DataFrame(; tipNames=treeleaves)

matching_tree_reconstruction = DataFrame(;
    tipNames=treeleaves, nodeNumber=range(1, length(treeleaves); step=1)
)

# Automatically name the traits with L or R prefixes

leftnames = "L".*string.(1:size(L,2))
traits_L = DataFrame(L, Symbol.(leftnames))
traits_L[!, "tipNames"] = metawab_usnames

rightnames = "R".*string.(1:size(R,1))
traits_R = DataFrame(R', rightnames)
traits_R[!, "tipNames"] = metawab_usnames

# Use a single dataframe for the traits
traits = leftjoin(traitframe, traits_L; on=:tipNames)
traits = leftjoin(traits, traits_R; on=:tipNames)

imputedtraits = DataFrame(;
    tipNames=[treeleaves; fill(missing, tree_net.numNodes - tree_net.numTaxa)]
);

Threads.@threads for coord in 1:size(L,2)
    for prefix in ["L", "R"]
        lower, upper, mean_trait = leaf_traits_reconstruction(
            traits[!, ["$(prefix)$(coord)", "tipNames"]],
            tree_net
        )
        imputedtraits[!, "$(prefix)$(coord)_low"] = lower
        imputedtraits[!, "$(prefix)$(coord)_up"] = upper
        imputedtraits[!, "$(prefix)$(coord)_mean"] = mean_trait
    end
end

# We save the reconstructed values
canadian_rec = innerjoin(dropmissing(imputedtraits), pool; on=:tipNames)

l = Array(canadian_rec_L[!, ["x1mean", "x2mean", "x3mean", "x4mean", "x5mean"]])

r = transpose(Array(canadian_rec_R[!, ["x1mean", "x2mean", "x3mean", "x4mean", "x5mean"]]))

lₗ = Array(canadian_rec_L[!, ["x1low", "x2low", "x3low", "x4low", "x5low"]])

lᵤ = Array(canadian_rec_L[!, ["x1up", "x2up", "x3up", "x4up", "x5up"]])

rₗ = transpose(Array(canadian_rec_R[!, ["x1low", "x2low", "x3low", "x4low", "x5low"]]))

rᵤ = transpose(Array(canadian_rec_R[!, ["x1up", "x2up", "x3up", "x4up", "x5up"]]))

ld = Matrix{Uniform}(undef, size(l))
for i in eachindex(ld)
    ld[i] = Uniform(lₗ[i], lᵤ[i])
end

rd = Matrix{Uniform}(undef, size(r))
for i in eachindex(rd)
    rd[i] = Uniform(rₗ[i], rᵤ[i])
end

draws = 20_000

Ld = [rand.(ld) for i in 1:draws]
Rd = [rand.(rd) for i in 1:draws]

Ns = [(Ld[i] * Rd[i]) .> 0.11 for i in 1:length(Ld)]
P = UnipartiteProbabilisticNetwork(
    reduce(.+, Ns) ./ draws, replace.(canadian_rec_L.tipNames, "_" => " ")
)

sort(interactions(P); by=(x) -> x.probability, rev=true)

histogram([x.probability for x in interactions(P)])
heatmap(adjacency(P))

omn = omnivory.(rand(P, 100))
O = Dict{String,Float64}([sp => mean([o[sp] for o in omn]) for sp in species(P)])
