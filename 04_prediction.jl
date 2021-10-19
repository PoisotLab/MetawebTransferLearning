# # Step 5 - the actual prediction

# This is the largest step in the entire pipeline, but not necessarily a complex
# one. In short, this will generate the entire paper, minus some of the data
# inflation and post-processing steps. As such, there are a few dependencies in
# play.

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

# These lines are ensuring that the figures all look uniform, and we also set a
# seed for reproducibility.

theme(:mute)
default(; frame=:box)
Random.seed!(01189998819991197253)

# This function will be used for leaf trait reconstruction:

function leaf_traits_reconstruction(traits,tree)
    ancestral_rec = ancestralStateReconstruction(traits, tree);
    predint_rec = predint(ancestral_rec, level = 0.3)
    lower = predint_rec[:,1]
    upper = predint_rec[:,2]
    mean_trait = mean(predint_rec,dims = 2)[:,1]
    return [lower,upper,mean_trait]
end

# ## Reading the data pieces

# Reading the tree is done with `PhyloNetworks` this time -- this is where we
# will get the Brownian motion code from:

tree_net = readTopology(joinpath("data", "mammals.newick"));

# We will grab back the corrected Europe/tree/GBIF names, as they will be used
# quite a lot to match simulations to species names:

namelist = DataFrame(CSV.File(joinpath("artifacts", "names_metaweb_tree_gbif.csv")))

# Finally, we grab the European metaweb *not* from the original file, but from
# our edgelist artifact.

eurometa = readdlm(joinpath("artifacts", "europeanmetaweb.csv"), ',', String)
mwspecies = unique(eurometa)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(eurometa, 1)
    M[eurometa[i, :]...] = true
end

# ## Finding the rank for truncation

# The ideal way to find a cutoff for the rank at which the matrix should be cut
# would involve the profile likelihood, or approximating the maximum curvature
# point on the screeplot using the central difference approximation for the
# second order partial derivative. Sadly all of these approaches will conclude
# that only the first dimension is required (there is a biological reason for
# this, to which we will return when we dig into the biological meaning of the
# latent variables).

# The only eigenvalues we care about are the ones that are up to the rank of the
# adjacency matrix - we will extract them, and range them so that they sum to
# unity. It makes no difference in the results, and allows us to read the
# proportion of variance explained directly.

rnk = rank(adjacency(M))
eig = svd(M).S[1:rnk]
neig = eig ./ sum(eig)

# The first piece of information is a screeplot of the eigenvalues. As a rule,
# the document rendered from the literate files will not include the figures, as
# we do not want to go into the analysis (there is, after all, a whole paper for
# that).

scatter(eig; lab="", dpi=600, size=(500, 500))
xaxis!("Dimension", (1, 38))
yaxis!("Eigenvalue", (0, 40))
savefig("figures/screeplot.png")

# The second diagnosis plot is the proportion of variance explained. We added
# lines marking the points that explain between 50% and 90% of the total
# variance, every 10%.

scatter(cumsum(neig); lab="", dpi=600, size=(500, 500))
for ve in [0.5, 0.6, 0.7, 0.8, 0.9]
    i = findfirst(cumsum(neig) .> ve)
    x = cumsum(neig)[i]
    plot!([0, i, NaN, i, i], [x, x, NaN, x, 0]; lab="", c=:grey, ls=:dash)
end
xaxis!("Dimension", (1, 38))
yaxis!("Cumulative variance explained", (0, 1.0))
savefig("figures/varexplained.png")

# At this point, we need to make an arbitrary call, which is to say to decide on
# a proportion of variance explained. We decided on 60%, which means that the
# first 12 ranks will be used.

# ## Embedding the European metaweb

L, R = rdpg(M, 12)

# The `L` and `R` variables are the left and right subspaces from the random dot
# product graph (which in turn are the left and right subspace from the t-SVD
# multiplied by the square root of the diagonal matrix containing the
# eigenvalues). Multiplying `L` and `R` (using a matrix multiplication!) will
# give a rank-12 approximation of the European metaweb.

# ## Thresholding the embedded network

# The reconstucted network (`L*R`) will not have Boolean values. Essentially,
# the multiplication will approximate something with values in {0,1} by
# something with values in ℝ, and so we need to find a cutoff to separate
# interactions from non-interactions. This threshold is important because it
# holds for any reconstruction made using these latent variables - that is to
# say, when we infer the Canadian interactions using reconstrusted left and
# right subspaces, we will apply the same threshold.

# Our basic approach to threshdolding is to examine 500 cutoffs points, get the
# confusion matrix, and retain as a threshold the value that maximizes
# Informedness:

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

# We can calculate Youden's J at every cutoff, and the threshold is going to be
# the threshold associated to the maximum Y:

Y = tp .+ tn .- 1
maxY, posY = findmax(Y)
threshold = thresholds[posY]

# We can plot the optimal cutoff - this is not presented in the manuscript, but
# useful in case you need to satisfy your curiosity. The best cutoff gives an
# informedness essentially = 1, which is expected.

plot(thresholds, Y; lw=0.0, fill=(0, 0.2), lab="", dpi=600, size=(500, 500))
scatter!([threshold], [maxY]; lab="", c=:darkgrey, msw=0.0)
yaxis!("Youden's J", (0, 1))
xaxis!("Cutoff", extrema(L * R))
savefig("figures/optimalcutoff.png")

# ## Visual examination of the subspaces

# The left and right subspaces *do* hold ecological information, and so it is a
# good idea to check them out visually. One striking result is that species with
# a value of 0 in the left subspace also have no preys: this is a strong clue
# that the left subspace is associated to generality (in the sense of Schoener
# 1989), a fact we will epxloit later on.

plot(
    heatmap(L; c=:GnBu, frame=:none, cbar=false),
    heatmap(R'; c=:BuPu, frame=:none, cbar=false);
    dpi=600,
    size=(500, 400),
)
savefig("figures/subspaces.png")

# ## Preparing the tree for inference

# We are now ready to move on to the next step: infering the values of the left
# and right subspaces for the species that are not in the European metaweb, but
# are in the tree. To do so, we will need the tip names:

treeleaves = tipLabels(tree_net)

# We are now ready to match the tree to the Canadian species pool names:

canada = DataFrame(CSV.File(joinpath("artifacts", "iucn_gbif_names.csv")))
cancodes = replace.(unique(filter(!ismissing, canada.gbifname)), " " => "_")
tree_cleanup = DataFrame(CSV.File(joinpath("artifacts", "upham_gbif_names.csv")))

csp = dropmissing!(DataFrame(; gbifname=canada.gbifname))
csp = dropmissing(leftjoin(csp, tree_cleanup; on=:gbifname))

# We can store the cleaned Canadian species name into a data frame - this is
# mandated by the `PhyloNetworks` interface:

canmammals = unique(csp.code)
pool = DataFrame(; tipNames=canmammals)

# We finally convert the network names to their underscored versions:

metaweb_can_names = replace.(M.S, " " => "_")
filter((x) -> (x) ∉ (metaweb_can_names ∩ treeleaves), metaweb_can_names)

# This prepare a data frame for the trait values:

traitframe = DataFrame(; tipNames=treeleaves)

# ## Infering the subspaces from the phylogeny

# We will prepare the data frame required by `PhyloNetworks` to store the
# reconstructed traits:

matching_tree_reconstruction = DataFrame(;
    tipNames=treeleaves, nodeNumber=range(1, length(treeleaves); step=1)
)

# To avoid writing long code, we prepare a series of columns with L and R as a
# prefix, and the number (from 1 to 12) of the corresponding dimension as a
# suffix:

leftnames = "L" .* string.(1:size(L, 2))
traits_L = DataFrame(L, Symbol.(leftnames))
traits_L[!, "tipNames"] = metaweb_can_names

rightnames = "R" .* string.(1:size(R, 1))
traits_R = DataFrame(R', rightnames)
traits_R[!, "tipNames"] = metaweb_can_names

# We can now merge these dataframes, and have something fully ready to be filled
# by the phylogenetic simulation:

traits = leftjoin(traitframe, traits_L; on=:tipNames)
traits = leftjoin(traits, traits_R; on=:tipNames)

imputedtraits = DataFrame(;
    tipNames=[treeleaves; fill(missing, tree_net.numNodes - tree_net.numTaxa)]
)

# The last step is to reconstruct each trait (L1 to L12, R1 to R12) for the
# entire tree. This is by far the longest part of the script, but it is not
# terribly long. This could probably be made thread parallel fairly easily but,
# this would also take more time to do than it takes to run.

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

# When the loop is done, we extract the values for the species in the Canadian
# species pool:

canadian_rec = innerjoin(dropmissing(imputedtraits), pool; on=:tipNames)

# ## Extracting the left and right subspaces from canada

# This part of the script uses the same characters as the equations in the paper
# - if you do not have a typeface with good unicode mathematical support
# installed, you might not get the full effect.

# Recall that a RDGP approximation is a matrix multiplication - we will
# therefore get the reconstructed average values after the Brownian motion
# model, and have a little look at them:

𝓁 = Array(canadian_rec[!, leftnames .* "_mean"])
𝓇 = transpose(Array(canadian_rec[!, rightnames .* "_mean"]))

plot(
    heatmap(𝓁; c=:GnBu, frame=:none, cbar=false),
    heatmap(𝓇'; c=:BuPu, frame=:none, cbar=false);
    dpi=600,
    size=(500, 400),
)
savefig("figures/imputed-subspaces.png")

# ## Generating a probabilistic network

# Because the Brownian motion model gives us a lower and upper bound, we will
# perform a series of random draws assuming that the values are uniformly
# distributed between these values. This step of the workflow can be adapted
# further. For example, one might want to account for the uncertainty in the
# phylogeny itself, or fit the distribution returned at each node rather than
# assuming a uniform distribution. We think that our approach introduces the
# least amount of guesses; it is likely to be over-estimating the chances of
# interactions a little, but this is the purpose of a metaweb: to give a list of
# possible interactions, to be later pared down.

𝓁ₗ = Array(canadian_rec[!, leftnames .* "_low"])
𝓇ₗ = transpose(Array(canadian_rec[!, rightnames .* "_low"]))

𝓁ᵤ = Array(canadian_rec[!, leftnames .* "_up"])
𝓇ᵤ = transpose(Array(canadian_rec[!, rightnames .* "_up"]))

# The distributions are expressed as actual Uniform distributions from the
# `Distributions` package.

ℒ = Matrix{Uniform}(undef, size(𝓁))
for i in eachindex(ℒ)
    ℒ[i] = Uniform(𝓁ₗ[i], 𝓁ᵤ[i])
end

ℛ = Matrix{Uniform}(undef, size(𝓇))
for i in eachindex(ℛ)
    ℛ[i] = Uniform(𝓇ₗ[i], 𝓇ᵤ[i])
end

# We will do a large enough number of draws:

draws = 20_000

𝐋 = [rand.(ℒ) for i in 1:draws]
𝐑 = [rand.(ℛ) for i in 1:draws]

# There are two pieces of information to keep in mind here. The first is that a
# RDPG is a matrix multiplication, so we simply need to multiply the 20000
# random subspaces, to get 20000 random matrices. The second is that these
# matrices give results not in {0,1} but in ℝ, but we have estimated an optimal
# threshold for this projection. Doing all this is a one-liner:

Ns = [(𝐋[i] * 𝐑[i]) .> threshold for i in 1:length(𝐋)]

# We can finally generate a probabilistic metaweb, in which the probability is
# defined as the proprtion of samples in which the interaction was inferred:

P = UnipartiteProbabilisticNetwork(
    reduce(.+, Ns) ./ draws, replace.(canadian_rec.tipNames, "_" => " ")
)

# We can have a little look at the interactions sorted by probabilities:

sort(interactions(P); by=(x) -> x.probability, rev=true)

# ## Visualising the results

# The next figures are very simple plots of the adjacency matrices:

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

# ## Writing output files for the raw predictions

# We will store the results in a data frame - the information we care about is
# the probability for the species pair, and whether the pair was also found in
# Europe, and if so, whether it interacted:

output = DataFrame(; from=String[], to=String[], score=Float64[], pair=Bool[], int=Bool[])
for int in interactions(P)
    pair = (int.from in species(M)) & (int.to in species(M))
    mint = pair ? M[int.from, int.to] : false
    push!(output, (int.from, int.to, int.probability, pair, mint))
end
sort!(output, [:score, :from, :to]; rev=[true, false, false])

# Save the basic network (no corrections)

CSV.write("artifacts/canadian_uncorrected.csv", output)

# ## Exploration of the relationship between subspaces and network properties

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

# ## Basic corrections

# We will directly bring knowledge from the European metaweb, meaning that if
# two species interact in Europe, we assume they also do in Canada (remember,
# these are metawebs, we only care about the biological feasibility of the
# interaction); if two species do not interact in Europe, we prevent them from
# interacting in Canada - this later point could be reversed, by inflating the
# European metaweb using the simulation results, and whether to apply this step
# at all can be considered on a case by case basis.

N = copy(P)
shared_species = filter(s -> s in species(M), species(P))
for s1 in shared_species
    for s2 in shared_species
        N[s1, s2] = M[s1, s2] ? 1.0 : 0.0
    end
end

# We may have introduced a number of 0s in the sparse matrix, and it is good
# hygiene to remove them.

SparseArrays.dropzeros!(N.edges)
simplify!(N)

# ## Writing the final metaweb

final = DataFrame(; from=String[], to=String[], score=Float64[])
for int in interactions(N)
    push!(final, (int.from, int.to, int.probability))
end
sort!(final, [:score, :from, :to]; rev=[true, false, false])

# Save the corrected network
CSV.write("artifacts/canadian_corrected.csv", final)

# ## Plots for the core results

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

# MaxEnt configuration model
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

