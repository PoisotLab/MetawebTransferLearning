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
L, R = rdpg(M, 8)

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
)

for coord in 1:size(L,2)
    for prefix in ["L", "R"]
        @info "Reconstructing $(prefix)$(coord)"
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

# TODO add the reference values from the European metaweb if they are known

# Get the left and right subspaces (and the lower and upper values)
𝓁 = Array(canadian_rec[!, leftnames.*"_mean"])
𝓇 = transpose(Array(canadian_rec[!, rightnames.*"_mean"]))

𝓁ₗ= Array(canadian_rec[!, leftnames.*"_low"])
𝓇ₗ = transpose(Array(canadian_rec[!, rightnames.*"_low"]))

𝓁ᵤ = Array(canadian_rec[!, leftnames.*"_up"])
𝓇ᵤ = transpose(Array(canadian_rec[!, rightnames.*"_up"]))

ℒ = Matrix{Uniform}(undef, size(𝓁))
for i in eachindex(ℒ)
    ℒ[i] = Uniform(𝓁ₗ[i], 𝓁ᵤ[i])
end

ℛ = Matrix{Uniform}(undef, size(𝓇))
for i in eachindex(ℛ)
    ℛ[i] = Uniform(𝓇ₗ[i], 𝓇ᵤ[i])
end

draws = 20_000

𝐋 = [rand.(ℒ) for i in 1:draws]
𝐑 = [rand.(ℛ) for i in 1:draws]

# TODO set the correct threshold
Ns = [(𝐋[i] * 𝐑[i]) .> 0.22 for i in 1:length(𝐋)]
P = UnipartiteProbabilisticNetwork(
    reduce(.+, Ns) ./ draws, replace.(canadian_rec.tipNames, "_" => " ")
)

sort(interactions(P); by=(x) -> x.probability, rev=true)

histogram([x.probability for x in interactions(P)])

heatmap(adjacency(P))

omn = omnivory.(rand(P, 100))

O = Dict{String,Float64}([sp => mean([o[sp] for o in omn]) for sp in species(P)])
