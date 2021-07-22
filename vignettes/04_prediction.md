# Step 5 - the actual prediction

This is the largest step in the entire pipeline, but not necessarily a complex
one. In short, this will generate the entire paper, minus some of the data
inflation and post-processing steps. As such, there are a few dependencies in
play.

````julia
using PhyloNetworks
using ProgressMeter
using LinearAlgebra
using EcologicalNetworks
using DelimitedFiles
using Random
using StatsPlots
using StatsBase
using Statistics
using SparseArrays
using CSV
using DataFrames
using GBIF
using Base.Threads
using StatsBase
using StatsModels
using Distributions
````

These lines are ensuring that the figures all look uniform, and we also set a
seed for reproducibility.

````julia
theme(:mute)
default(; frame=:box)
Random.seed!(01189998819991197253)
````

Some functions that are required for phylogenetic imputation are stored in
another file, and we will load them here to be done with that.

````julia
include("lib/pwar.jl")
````

## Reading the data pieces

Reading the tree is done with `PhyloNetworks` this time -- this is where we
will get the Brownian motion code from:

````julia
tree_net = readTopology(joinpath("data", "mammals.newick"));
````

We will grab back the corrected Europe/tree/GBIF names, as they will be used
quite a lot to match simulations to species names:

````julia
namelist = DataFrame(CSV.File(joinpath("artifacts", "names_metaweb_tree_gbif.csv")))
````

Finally, we grab the European metaweb *not* from the original file, but from
our edgelist artifact.

````julia
eurometa = readdlm(joinpath("artifacts", "europeanmetaweb.csv"), ',', String)
mwspecies = unique(eurometa)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(eurometa, 1)
    M[eurometa[i, :]...] = true
end
````

## Finding the rank for truncation

The ideal way to find a cutoff for the rank at which the matrix should be cut
would involve the profile likelihood, or approximating the maximum curvature
point on the screeplot using the central difference approximation for the
second order partial derivative. Sadly all of these approaches will conclude
that only the first dimension is required (there is a biological reason for
this, to which we will return when we dig into the biological meaning of the
latent variables).

The only eigenvalues we care about at the ones that are up to the rank of the
adjacency matrix - we will extract them, and range them so that they sum to
unity. It makes no difference in the results, and allows to read the
proportion of variance explained directly.

````julia
rnk = rank(adjacency(M))
eig = svd(M).S[1:rnk]
neig = eig ./ sum(eig)
````

The first piece of information is a screeplot of the eigenvalues. As a rule,
the document rendered from the literate files will not include the figures, as
we do not want to go into the analysis (there is, after all, a whole paper for
that).

````julia
scatter(eig; lab="", dpi=600, size=(500, 500))
xaxis!("Dimension", (1, 38))
yaxis!("Eigenvalue", (0, 40))
savefig("figures/screeplot.png")
````

The second diagnosis plot is the proportion of variance explained. We added
lines marking the points that explain between 50% and 90% of the total
variance, every 10%.

````julia
scatter(cumsum(neig); lab="", dpi=600, size=(500, 500))
for ve in [0.5, 0.6, 0.7, 0.8, 0.9]
    i = findfirst(cumsum(neig) .> ve)
    x = cumsum(neig)[i]
    plot!([0, i, NaN, i, i], [x, x, NaN, x, 0]; lab="", c=:grey, ls=:dash)
end
xaxis!("Dimension", (1, 38))
yaxis!("Cumulative variance explained", (0, 1.0))
savefig("figures/varexplained.png")
````

At this point, we need to make an arbitrary call, which is to say to decide on
a proportion of variance explained. We decided on 60%, which means that the
first 12 ranks will be used.

## Embedding the European metaweb

````julia
L, R = rdpg(M, 12)
````

The `L` and `R` variables are the left and right subspaces from the random dot
product graph (which in turn are the left and right subspace from the t-SVD
mutliplied by the square root of the diagonal matrix containing the
eigenvalues). Multiplying `L` and `R` (using a matrix multiplication!) will
give a rank-12 approximation of the European metaweb.

## Thresholding the embedded network

The reconstucted network (`L*R`) will not have Boolean values. Essentially,
the multiplication will approximate something with values in ${0,1}$ by
something with values in ℝ, and so we need to find a cutoff to separate
interactions from non-interactions. This threshold is important because it
holds for any reconstruction made using these latent variables - that is to
say, when we infer the Canadian interactions using reconstrusted left and
right subspaces, we will apply the same threshold.

Our basic approach to threshdolding is to examine 500 cutoffs points, get the
confusion matrix, and retain as a threshold the value that maximizes
Informedness:

````julia
A = adjacency(M)
thresholds = LinRange(extrema(L * R)..., 500)
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
````

We can calculate Youden's J at every cutoff, and the threshold is going to be
the threshold associated to the maximum Y:

````julia
Y = tp .+ tn .- 1
maxY, posY = findmax(Y)
threshold = thresholds[posY]
````

We can plot the optimal cutoff - this is not presented to the manuscript, but
useful in case you need to satisfy your curiosity. The best cutoff gives an
informedness essentially = 1, which is expected.

````julia
plot(thresholds, Y; lw=0.0, fill=(0, 0.2), lab="", dpi=600, size=(500, 500))
scatter!([threshold], [maxY]; lab="", c=:darkgrey, msw=0.0)
yaxis!("Youden's J", (0, 1))
xaxis!("Cutoff", extrema(L * R))
savefig("figures/optimalcutoff.png")
````

## Visual examination of the subspaces

The left and right subspaces *do* hold ecological information, and

````julia
plot(
    heatmap(L; c=:GnBu, frame=:none, cbar=false),
    heatmap(R'; c=:BuPu, frame=:none, cbar=false);
    dpi=600,
    size=(500, 400),
)
savefig("figures/subspaces.png")
````

Get the species names

````julia
treeleaves = tipLabels(tree_net)

canada = DataFrame(CSV.File(joinpath("artifacts", "iucn_gbif_names.csv")))
cancodes = replace.(unique(filter(!ismissing, canada.gbifname)), " " => "_")
tree_cleanup = DataFrame(CSV.File(joinpath("artifacts", "upham_gbif_names.csv")))

csp = dropmissing!(DataFrame(; gbifname=canada.gbifname))
csp = dropmissing(leftjoin(csp, tree_cleanup; on=:gbifname))

canmammals = unique(csp.code)
````

We start getting a species pool to infer the traits

````julia
pool = DataFrame(; tipNames=canmammals)

metaweb_can_names = replace.(M.S, " " => "_")
metaweb_can_names ∩ canmammals

#### I am missing some species from the metaweb
filter((x) -> (x) ∉ (metaweb_can_names ∩ treeleaves), metaweb_can_names)

traitframe = DataFrame(; tipNames=treeleaves)

matching_tree_reconstruction = DataFrame(;
    tipNames=treeleaves, nodeNumber=range(1, length(treeleaves); step=1)
)
````

Automatically name the traits with L or R prefixes

````julia
leftnames = "L" .* string.(1:size(L, 2))
traits_L = DataFrame(L, Symbol.(leftnames))
traits_L[!, "tipNames"] = metaweb_can_names

rightnames = "R" .* string.(1:size(R, 1))
traits_R = DataFrame(R', rightnames)
traits_R[!, "tipNames"] = metaweb_can_names
````

Use a single dataframe for the traits

````julia
traits = leftjoin(traitframe, traits_L; on=:tipNames)
traits = leftjoin(traits, traits_R; on=:tipNames)
````

Imputed traits dataframe

````julia
imputedtraits = DataFrame(;
    tipNames=[treeleaves; fill(missing, tree_net.numNodes - tree_net.numTaxa)]
)
````

Imputation loop!

````julia
for coord in 1:size(L, 2)
    for prefix in ["L", "R"]
        @info "Reconstructing $(prefix)$(coord)"
        lower, upper, mean_trait = leaf_traits_reconstruction(
            traits[!, ["$(prefix)$(coord)", "tipNames"]], tree_net
        )
        imputedtraits[!, "$(prefix)$(coord)_low"] = lower
        imputedtraits[!, "$(prefix)$(coord)_up"] = upper
        imputedtraits[!, "$(prefix)$(coord)_mean"] = mean_trait
    end
end
````

We save the reconstructed values

````julia
canadian_rec = innerjoin(dropmissing(imputedtraits), pool; on=:tipNames)
````

Get the left and right subspaces (and the lower and upper values)

````julia
𝓁 = Array(canadian_rec[!, leftnames .* "_mean"])
𝓇 = transpose(Array(canadian_rec[!, rightnames .* "_mean"]))

plot(
    heatmap(𝓁; c=:GnBu, frame=:none, cbar=false),
    heatmap(𝓇'; c=:BuPu, frame=:none, cbar=false);
    dpi=600,
    size=(500, 400),
)
savefig("figures/imputed-subspaces.png")

𝓁ₗ = Array(canadian_rec[!, leftnames .* "_low"])
𝓇ₗ = transpose(Array(canadian_rec[!, rightnames .* "_low"]))

𝓁ᵤ = Array(canadian_rec[!, leftnames .* "_up"])
𝓇ᵤ = transpose(Array(canadian_rec[!, rightnames .* "_up"]))

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
````

We get the thresholded networks here

````julia
Ns = [(𝐋[i] * 𝐑[i]) .> threshold for i in 1:length(𝐋)]
P = UnipartiteProbabilisticNetwork(
    reduce(.+, Ns) ./ draws, replace.(canadian_rec.tipNames, "_" => " ")
)
````

Deterministic version

````julia
N = UnipartiteNetwork(𝓁 * 𝓇 .>= threshold, replace.(canadian_rec.tipNames, "_" => " "))

sort(interactions(P); by=(x) -> x.probability, rev=true)

histogram([x.probability for x in interactions(P)])

sporder = sortperm(vec(sum(adjacency(P); dims=2)))
h1 = heatmap(
    adjacency(P)[sporder, sporder];
    c=:Greys,
    frame=:none,
    cbar=false,
    dpi=600,
    size=(500, 500),
    aspectratio=1,
)

sporder = sortperm(vec(sum(adjacency(M); dims=2)))
h2 = heatmap(
    adjacency(M)[sporder, sporder];
    c=:Greys,
    frame=:none,
    cbar=false,
    dpi=600,
    size=(500, 500),
    aspectratio=1,
)

plot(h2, h1; size=(1000, 500))
savefig("figures/adjacencymatrices.png")

output = DataFrame(; from=String[], to=String[], score=Float64[], pair=Bool[], int=Bool[])
for int in interactions(P)
    pair = (int.from in species(M)) & (int.to in species(M))
    mint = pair ? M[int.from, int.to] : false
    push!(output, (int.from, int.to, int.probability, pair, mint))
end
sort!(output, [:score, :from, :to]; rev=[true, false, false])
````

Save the basic network (no corrections)

````julia
CSV.write("artifacts/canadian_uncorrected.csv", output)
````

Exploration of the relationship between subspaces and network properties

````julia
kout = degree(P; dims=1)
kin = degree(P; dims=2)
ordered_sp = replace.(canadian_rec.tipNames, "_" => " ")

scatter(
    𝓁[:, 1], [kout[s] / richness(P) for s in ordered_sp]; dpi=600, size=(500, 500), lab=""
)
xaxis!("Position in the left subspace", extrema(vcat(𝓇', 𝓁)))
yaxis!("Probabilistic generality", (0, 1))
savefig("figures/left-gen.png")

scatter(
    𝓇'[:, 1],
    [kin[s] / richness(P) for s in ordered_sp];
    dpi=600,
    size=(500, 500),
    lab="",
    legend=:bottomright,
)
xaxis!("Position in the right subspace", extrema(vcat(𝓇', 𝓁)))
yaxis!("Probabilistic vulnerability")
savefig("figures/right-vuln.png")
````

Corrections assuming
- if species don't interact in Europe, no interaction in Canada
- if species interact in Europe, interaction in Canada

````julia
N = copy(P)
shared_species = filter(s -> s in species(M), species(P))
for s1 in shared_species
    for s2 in shared_species
        N[s1, s2] = M[s1, s2] ? 1.0 : 0.0
    end
end
SparseArrays.dropzeros!(N.edges)
simplify!(N)
````

Final metaweb

````julia
final = DataFrame(; from=String[], to=String[], score=Float64[])
for int in interactions(N)
    push!(final, (int.from, int.to, int.probability))
end
sort!(final, [:score, :from, :to]; rev=[true, false, false])
````

Save the corrected network

````julia
CSV.write("artifacts/canadian_corrected.csv", final)

#%% Plot
l = @layout [
    a{0.3w} b
    c{0.7h} d
]

sporder = sortperm(vec(sum(adjacency(N); dims=2)))

plot(
    plot(; legend=false, axes=false, frame=:none),
    heatmap(𝓇[:, sporder]; frame=:none, legend=false, c=:BrBG, clim=(-1, 1)),
    heatmap(𝓁[sporder, :]; frame=:none, legend=false, c=:PRGn, clim=(-1, 1)),
    heatmap(adjacency(N)[sporder, sporder]; c=:Greys, frame=:none, legend=false);
    layout=l,
    size=(500, 500),
    dpi=600,
)

savefig("figures/combined-prediction.png")

sporder = sortperm(vec(sum(adjacency(M); dims=2)))

l = @layout [
    a{0.3w} b
    c{0.7h} d
]

plot(
    plot(; legend=false, axes=false, frame=:none),
    heatmap(R[:, sporder]; frame=:none, legend=false, c=:BrBG, clim=(-1, 1)),
    heatmap(L[sporder, :]; frame=:none, legend=false, c=:PRGn, clim=(-1, 1)),
    heatmap(adjacency(M)[sporder, sporder]; c=:Greys, frame=:none, legend=false);
    layout=l,
    size=(500, 500),
    dpi=600,
)

savefig("figures/combined-empirical.png")
````

MaxEnt configuration model

````julia
a = zeros(Float64, size(adjacency(P)))
C = UnipartiteProbabilisticNetwork(a, EcologicalNetworks._species_objects(P)...)

for s1 in species(P; dims=1), s2 in species(P; dims=2)
    C[s1, s2] = 0.5(kout[s1] / richness(P) + kin[s2] / richness(P))
end

rec = 𝓁[:, 1] * hcat(𝓇'[:, 1]...)

Ceff = adjacency(C)
Cinf = rec

Zeff = (Ceff .- mean(Ceff)) ./ std(Ceff)
Zinf = (Cinf .- mean(Cinf)) ./ std(Cinf)

density(sqrt.(vec(Zeff .- Zinf) .^ 2.0); size=(500, 500), dpi=600, fill=(0, 0.2), lab="")
xaxis!("Mean squared error", (0, 5))
yaxis!("Density", (0, 3))
savefig("figures/distance-configuration.png")

sporder = sortperm(vec(sum(adjacency(C); dims=2)))

ΔZ = Zeff .- Zinf

heatmap(
    ΔZ[sporder, sporder];
    clim=(-3, 3),
    c=:PuOr,
    size=(500, 500),
    dpi=600,
    frame=:none,
    aspectratio=1,
)
savefig("figures/heatmap-configuration.png")
````

