1.1

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
using StatsBase
using StatsModels

# Read the tree
tree_net = readTopology(joinpath("data", "mammals.newick"));

# Read the names
namelist = DataFrame(CSV.File(joinpath("artifacts", "names_metaweb_tree_gbif.csv")))

# Prepare the metaweb
mwspecies = unique(namelist.name)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)

# Read the metaweb
speciescodes = readdlm(joinpath("data", "Spp_Id.txt"))[2:end,:]
speciesdict = Dict([speciescodes[i,1] => speciescodes[i,2] for i in 1:size(speciescodes,1)])
mwlines = readlines(joinpath("data", "Metaweb_adults.csv"));
mwhead = [speciesdict[sp] for sp in replace.(split(mwlines[1], ","), '"' => "")[2:end]]
for row in mwlines[2:end]
    splitrow = replace.(split(row, ","), '"' => "")
    from = speciesdict[splitrow[1]]
    # Real name?
    realname = namelist[isequal(from).(namelist.metaweb),:name]
    if length(realname) == 0
        continue
    else
        sp_from = only(realname)
        int = findall(isequal("1"), splitrow[2:end])
        if !isempty(int)
            to = mwhead[int]
            to_names = namelist[map(n -> n in to, namelist.metaweb),:name]
            for t in to_names
                M[sp_from,t] = true
            end
        end
    end
end

simplify!(M)
# Partition the metaweb into two latent subspaces
L, R = rdpg(M, 5)

treeleaves = tipLabels(tree_net)
canada = DataFrame(CSV.File(joinpath("data", "canada.csv")))
canada = canada[canada.taxonRank .== "SPECIES", :]
canada = canada[canada.numberOfOccurrences .> 1 , :]
canada = canada[canada.order .!== "Cetacea" , :]
cancodes = replace.(unique(filter(!ismissing, canada.species)), " " => "_")
tree_cleanup = DataFrame(CSV.File(joinpath("artifacts", "upham_gbif_names.csv")))


csp = dropmissing!(DataFrame(gbifname = canada.species))
csp = dropmissing(leftjoin(csp, tree_cleanup, on=:gbifname))

canmammals = unique(csp.code)

## Here I diverge from Tim, who does k=3 phylogenetic kmeans

pool = DataFrame(tipNames = canmammals)

metawab_usnames =  replace.(M.S, " " => "_")
metawab_usnames ∩ canmammals

#### I am missing some species from the metaweb
filter((x) -> (x) ∉ (metawab_usnames ∩ treeleaves), metawab_usnames)

traitframe = DataFrame(tipNames = treeleaves)
matching_tree_reconstruction  = DataFrame(
    tipNames = treeleaves,
    nodeNumber = range(1,length(treeleaves), step = 1)
)

traits_L = DataFrame(L)
traits_L[!,"tipNames"] = metawab_usnames
traits_R = DataFrame(R')
traits_R[!,"tipNames"] = metawab_usnames

traits_L = leftjoin(traitframe,traits_L, on = :tipNames)
traits_R = leftjoin(traitframe,traits_R, on = :tipNames)

# phyloNetworklm(@formula(x1 ~ 1), traits_L, tree_net, model="lambda")

Results_L = DataFrame(tipNames = [treeleaves; fill(missing,tree_net.numNodes - tree_net.numTaxa)]);
Results_R = DataFrame(tipNames = [treeleaves; fill(missing,tree_net.numNodes - tree_net.numTaxa)]);


for coord in 1:size(L)[2]
    
    lower,upper,mean_trait = leaf_traits_reconstruction(traits_L[!,["x$(coord)", "tipNames"]], tree_net);

    Results_L[!,"x$(coord)low"] = lower;
    Results_L[!,"x$(coord)up"] = upper;
    Results_L[!,"x$(coord)mean"] =  mean_trait;

end

for coord in 1:size(R')[2]
    lower,upper,mean_trait = leaf_traits_reconstruction(traits_R[!,["x$(coord)", "tipNames"]], tree_net);

    Results_R[!,"x$(coord)low"] = lower;
    Results_R[!,"x$(coord)up"] = upper;
    Results_R[!,"x$(coord)mean"] = mean_trait;
end


canadian_rec_L = innerjoin(dropmissing(Results_L), pool, on = :tipNames)
canadian_rec_R = innerjoin(dropmissing(Results_R), pool, on = :tipNames)

l = canadian_rec_L[!,
    ["x1mean", "x2mean", "x3mean", "x4mean", "x5mean"]] |>
    Array

r = canadian_rec_R[!,
    ["x1mean", "x2mean", "x3mean", "x4mean", "x5mean"]] |>
    Array |>
    transpose

#P = simplify(UnipartiteQuantitativeNetwork(, pool))
l*r |> heatmap

# we compare it with Timothée work
using Phylo
tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]
k = 3
# I need to change this to my ordering of the beasts
pool = canadian_rec_R.tipNames
poolsize = length(pool)
lₜ = zeros(Float64, poolsize, size(L, 2))
rₜ = permutedims(lₜ)
Threads.@threads for i in 1:length(pool)
    target = pool[i]
    D = [distance(tree, target, first(namelist[isequal(n).(namelist.name),:upham])) for n in species(M)]
    p = sortperm(D)
    if iszero(minimum(D))
        rₜ[:,i] = R[:,p][:,1]
        lₜ[i,:] = L[p,:][1,:]
    else
        rₜ[:,i] = mean(R[:,p][:,1:k]; dims=2)
        lₜ[i,:] = mean(L[p,:][1:k,:]; dims=1)
    end
end

lₜ*rₜ |> heatmap

# It's all over the place
# but not thaaaat much given the huge uncertainty
# we have in the data itself 
scatter(vec(l*r),vec(lₜ*rₜ),
    lab = "",
    markersize = 1
)
title!("Raw interaction scores from two traits inference methods:")
xaxis!("Phylogenetic Brownian motion")
yaxis!("Phylogenetic kmeans")

savefig("figures/Phylobbm_vs_phylo_kmeans.png")