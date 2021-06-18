using Phylo
using LinearAlgebra
using EcologicalNetworks
using DelimitedFiles
using Plots
using StatsBase
using Statistics
using EcologicalNetworksPlots
using CSV
using DataFrames
using GBIF
using ProgressMeter

theme(:mute)
default(; frame=:box)

# Read the metaweb
eurometa = readdlm(joinpath("artifacts", "europeanmetaweb.csv"), ',', String)
mwspecies = unique(eurometa)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(eurometa, 1)
    M[eurometa[i, :]...] = true
end

# Get the correct rank to cut the matrix at
rnk = rank(Array(M.edges))
eig = svd(M).S[1:rnk]
eig = eig./sum(eig)

𝒻 = (d, μ, σ²) -> (1.0 / (2π * σ²)) * exp(-((d - μ)^2.0) / (2σ²))
ℒ = zeros(Float64, rnk)

for p in 2:(rnk-2)
    θ₁ = eig[1:p]
    θ₂ = eig[(p+1):end]
    μ₁ = mean(θ₁)
    μ₂ = mean(θ₂)
    σ² = (1.0/(rnk-2))*((p-1)*var(θ₁) + (rnk-p-1)*var(θ₂))
    𝒻₁ = (d) -> 𝒻(d, μ₁, σ²)
    𝒻₂ = (d) -> 𝒻(d, μ₂, σ²)
    𝓁ₚ = sum(log.(𝒻₁.(θ₁))) + sum(log.(𝒻₂.(θ₂)))
    ℒ[p] = 𝓁ₚ 
    @info "$(p) → $(𝓁ₚ)"
end

# Partition the metaweb into two latent subspaces
L, R = rdpg(M, 21)

A = adjacency(M)
thresholds = LinRange(extrema(L * R)..., 100)
tp = zeros(Float64, length(thresholds))
tn = similar(tp)
fp = similar(tp)
fn = similar(tp)
for (i, t) in enumerate(thresholds)
    PN = (L * R) .>= t
    tp[i] = sum((A) .& (PN)) / sum(A)
    tn[i] = sum((.!(A)) .& (.!(PN))) / sum(.!(A))
    fp[i] = sum((.!(A)) .& (PN)) / sum(.!(A))
    fn[i] = sum((A) .& (.!(PN))) / sum(A)
end

Y = tp .+ tn .- 1
threshold = thresholds[findmax(Y)[2]]

plot(fp, tp; lab="", aspectratio=1, xlim=(0, 1), fill=(0, 0.2), ylim=(0, 1), dpi=600)
xaxis!("False positives")
yaxis!("True positives")
savefig("figures/roc.png")

heatmap(L * R)
heatmap((L * R) .> threshold)

# Get a random species pool from the tree
tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]
treenodes = [n.name for n in tree.nodes if !startswith(n.name, "Node ")]
canada = DataFrame(CSV.File(joinpath("artifacts", "iucn_gbif_names.csv")))
cancodes = replace.(unique(filter(!ismissing, canada.gbifname)), " " => "_")
tree_cleanup = DataFrame(CSV.File(joinpath("artifacts", "upham_gbif_names.csv")))

csp = dropmissing!(DataFrame(; gbifname=canada.gbifname))
csp = dropmissing(leftjoin(csp, tree_cleanup; on=:gbifname))

canmammals = unique(csp.code)
namelist = DataFrame(CSV.File(joinpath("artifacts", "names_metaweb_tree_gbif.csv")))

k = 6
pool = canmammals
poolsize = length(pool)
l = zeros(Float64, poolsize, size(L, 2))
r = permutedims(l)
prog = Progress(length(pool))
Threads.@threads for i in 1:length(pool)
    target = pool[i]
    D = [
        distance(tree, target, first(namelist[isequal(n).(namelist.name), :upham])) for
        n in species(M)
    ]
    p = sortperm(D)
    if iszero(minimum(D))
        r[:, i] = R[:, p][:, 1]
        l[i, :] = L[p, :][1, :]
    else
        r[:, i] = mean(R[:, p][:, 1:k]; dims=2)
        l[i, :] = mean(L[p, :][1:k, :]; dims=1)
    end
    next!(prog)
end

P = simplify(UnipartiteQuantitativeNetwork(l * r, pool))

s_int = sort(
    filter(x -> x.strength >= threshold, interactions(P)); by=(x) -> x.strength, rev=true
)
predint = DataFrame(; predator=String[], prey=String[], evidence=Float64[])
canpool = unique(csp.gbifname)
N = UnipartiteNetwork(zeros(Bool, length(canpool), length(canpool)), canpool)
for i in s_int
    fr = first(unique(csp[csp.code .== i.from, :gbifname]))
    to = first(unique(csp[csp.code .== i.to, :gbifname]))
    push!(predint, (fr, to, i.strength))
    N[fr, to] = true
end

simplify!(N)

open(joinpath("artifacts", "canadianmetaweb.csv"), "w") do caio
    for int in sort(interactions(N); by=x -> x.from)
        println(caio, "$(int.from),$(int.to)")
    end
end

