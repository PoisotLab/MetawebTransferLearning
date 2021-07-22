# Phylogenetic inference of ecological interactions through network embedding

**Executive summary**: 

1. network decomposition through a truncated-SVD capture the evolutionary signal of species interactions
2. the decomposition into left/right subspaces is unique, represents latent traits for resp. outgoing and incoming edges, and can be used to predict interactions
3. we can infer latent traits for unobserved species based on phylogenetic proximity with observed species
4. we use this information to transfer information from the trophic interaction of European mammals to Canadian mammals for which we have no *a priori* data


# Step 1 - name matching

The names in the parts of the datasets are using different versions of the
taxonomy, and so we will reconcile everything using GBIF. The process is
entirely automated, and required no human decision. We rely on a mix of strict
matching, and then search through synonyms.

## Dependencies

Nothing here is out of the extroardinary - we do some processing on multiple
threads, and use `ProgressMeter` to provide visual feedback on the time to
completion.

````julia
using Phylo
using GBIF
using DataFrames
using CSV: CSV
using ProgressMeter
using DelimitedFiles
using Base.Threads
````

## Cleaning the European metaweb

This step is extracting the species names from the species list in the
European metaweb - we will need to go through two steps: removing the
non-mammals by filtering on the species code, and then extracting the names
these species were given in the European metaweb.

````julia
metaweb_names = String.(readdlm(joinpath("data", "Spp_Id.txt")))
mammal_positions = findall(startswith("M"), metaweb_names[:, 1])
````

This line specifically is going to create the list of mammal names are
recorded in the European metaweb. As such, this is the first one we will
clean.

````julia
mammal_names = metaweb_names[mammal_positions, 2]
````

To ensure that we can join dataframes together, we will create a data frame
with the species code, the name and ID of the matched species in the GBIF
backbone, and finally, a flag to check that the two names are equal. This flag
is useful for manual inspection; most often, names are unequal because te
European metaweb uses deprecated taxonomic names.

Note that we are actually creating *one* dataframe per thread. This avoids
several threads pushing to the same position at the dataframe, which always
results in an error. There are, of course, other ways to avoid the issue, but
this one is simple and takes a single line to reconcile.

````julia
gbif_cleanup_components = [
    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for
    i in 1:nthreads()
]
````

With the dataframe stores in place, we can distribute the name reconciliation
on the different threads. The core of this loop is a search on the GBIF web
API; this step is therefore going to be limited by the responses of GBIF. Note
that although we use unstrict search (to allow for synonyms, etc), we restrict
the search to *Mammalia*.

````julia
p = Progress(length(mammal_names))
@threads for i in 1:length(mammal_names)
    cname = replace(mammal_names[i], '_' => ' ')
    try
        tax = GBIF.taxon(cname; strict=false, class="Mammalia")
        push!(
            gbif_cleanup_components[threadid()],
            (mammal_names[i], tax.species[1], tax.species[2], cname == tax.species[1]),
        )
    catch
        continue
    end
    next!(p)
end
````

Data frames are now filled, and we need to bring them together in a single one
-- this is a one-line call to `vcat`, as we simply need to stack them
vertically.

````julia
gbif_cleanup = vcat(gbif_cleanup_components...)
````

## Clean the GBIF matches

Names in the European metaweb that are not matched are removed as a
precaution:

````julia
select!(gbif_cleanup, Not(:equal))
````

We finally rename the columns, to facilitate the joins later in the pipeline:

````julia
rename!(gbif_cleanup, :code => :metaweb)
rename!(gbif_cleanup, :gbifname => :name)
rename!(gbif_cleanup, :gbifid => :id)
````

We save the file in the `artifacts` folders, as a csv:

````julia
CSV.write(joinpath("artifacts", "metaweb_gbif.csv"), gbif_cleanup)
````

## Cleaning the phylogeny

The phylogeny suffers from much the same problem as the European metaweb. We
will therefore use much of the same solution to fix the names.

The first step is to read the tree as a nexus file using `Phylo`, and to get
rid of all the internal nodes -- we are only interested in the leaves of the
tree.

````julia
tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]
treenodes = [n.name for n in tree.nodes if !startswith(n.name, "Node ")]
````

This step is a little dirty. Remember that data cleaning is a sin eater, and
we need to eat other people's data sins at some point. This is this point.
Essentially, we will iterate through the rows of the cleaned metaweb names,
and match on either `Gen sp` or `Gen_sp`.

There is an edge case that was not solved automatically, and so the *S.
musicus* species gets its correct tree name manually. This name is caught
correctly by *e.g.* `NCBITaxonomy.jl`, but this introduced another dependency
in the project, and requires another re-harmonization with the GBIF backbone.

````julia
n = []
for row in eachrow(gbif_cleanup)
    if replace(row.name, ' ' => '_') in treenodes
        push!(n, replace(row.name, ' ' => '_'))
    elseif row.metaweb in treenodes
        push!(n, row.metaweb)
    elseif row.name == "Spermophilus musicus"
        push!(n, "Spermophilus_pygmaeus")
    else
    end
end
````

## Merging the cleaned names together

Because the positions in `n` match with the row indices in the cleaned metaweb
names, we add this object as a column:

````julia
gbif_cleanup.upham = n
````

We finally save this new file as another artifact. Generally, we will err on
the side of caution and save multiple files with redudancy in them. This is
good practice, and it allows for easier data dumpster diving if something goes
wrong.

````julia
CSV.write(joinpath("artifacts", "names_metaweb_tree_gbif.csv"), gbif_cleanup)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*



````julia
using Phylo
using GBIF
using CSV, DataFrames
using ProgressMeter
using Base.Threads

tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]
treenodes = [n.name for n in tree.nodes if !startswith(n.name, "Node ")]

tree_cleanup_components = [
    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for
    i in 1:nthreads()
]

p = Progress(length(treenodes))
@threads for i in 1:length(treenodes)
    cname = replace(treenodes[i], '_' => ' ')
    try
        tax = GBIF.taxon(cname; strict=false, class="Mammalia")
        push!(
            tree_cleanup_components[threadid()],
            (treenodes[i], tax.species[1], tax.species[2], cname == tax.species[1]),
        )
    catch
        continue
    end
    next!(p)
end

tree_cleanup = vcat(tree_cleanup_components...)
CSV.write(joinpath("artifacts", "upham_gbif_names.csv"), tree_cleanup)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*



````julia
using GBIF
using CSV, DataFrames
using ProgressMeter
using Base.Threads

checklist = DataFrame(CSV.File(joinpath("data", "taxonomy.csv")))
valid_rows = map(
    fam ->
        !(
            fam ∈ [
                "BALAENIDAE",
                "PHYSETERIDAE",
                "DELPHINIDAE",
                "BALAENOPTERIDAE",
                "OTARIIDAE",
                "PHOCIDAE",
                "ODOBENIDAE",
                "ZIPHIIDAE",
                "MONODONTIDAE",
                "ESCHRICHTIIDAE",
                "KOGIIDAE",
                "PHOCOENIDAE"
            ]
        ),
    checklist.familyName,
)
checklist = checklist[findall(valid_rows), :]

extinct_sp = map(
    sp ->
        !(
            sp ∈ [
                "Neovison macrodon",
                "Enhydra lutris"
            ]
        ),
    checklist.scientificName,
)
checklist = checklist[findall(extinct_sp), :]

checklist_cleanup_components = [
    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for
    i in 1:nthreads()
]

checklist_cleanup_components = [
    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for
    i in 1:nthreads()
]

p = Progress(length(checklist.scientificName))
@threads for i in 1:length(checklist.scientificName)
    cname = replace(checklist.scientificName[i], '_' => ' ')
    try
        tax = GBIF.taxon(cname; strict=false, class="Mammalia")
        push!(
            checklist_cleanup_components[threadid()],
            (
                checklist.scientificName[i],
                tax.species[1],
                tax.species[2],
                cname == tax.species[1],
            ),
        )
    catch
        continue
    end
    next!(p)
end

checklist_cleanup = vcat(checklist_cleanup_components...)
CSV.write(joinpath("artifacts", "iucn_gbif_names.csv"), checklist_cleanup)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*



````julia
using Phylo
using EcologicalNetworks
using DelimitedFiles
using CSV
using DataFrames
using ProgressMeter
````

Read the tree

````julia
tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]
````

Read the names

````julia
namelist = DataFrame(CSV.File(joinpath("artifacts", "names_metaweb_tree_gbif.csv")))
````

Prepare the metaweb

````julia
mwspecies = unique(namelist.name)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
````

Read the metaweb

````julia
speciescodes = readdlm(joinpath("data", "Spp_Id.txt"))[2:end, :]
speciesdict = Dict([
    speciescodes[i, 1] => speciescodes[i, 2] for i in 1:size(speciescodes, 1)
])
mwlines = readlines(joinpath("data", "Metaweb_adults.csv"));
mwhead = [speciesdict[sp] for sp in replace.(split(mwlines[1], ","), '"' => "")[2:end]]
for row in mwlines[2:end]
    splitrow = replace.(split(row, ","), '"' => "")
    from = speciesdict[splitrow[1]]
````

Real name?

````julia
    realname = namelist[isequal(from).(namelist.metaweb), :name]
    if length(realname) == 0
        continue
    else
        sp_from = only(realname)
        int = findall(isequal("1"), splitrow[2:end])
        if !isempty(int)
            to = mwhead[int]
            to_names = namelist[map(n -> n in to, namelist.metaweb), :name]
            for t in to_names
                M[sp_from, t] = true
            end
        end
    end
end

simplify!(M)

open(joinpath("artifacts", "europeanmetaweb.csv"), "w") do euio
    for int in sort(interactions(M); by=x -> x.from)
        println(euio, "$(int.from),$(int.to)")
    end
end
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*



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

theme(:mute)
default(; frame=:box)
Random.seed!(01189998819991197253)
````

Get ancillary functions

````julia
include("lib/pwar.jl")
````

Read the tree

````julia
tree_net = readTopology(joinpath("data", "mammals.newick"));
````

Read the names

````julia
namelist = DataFrame(CSV.File(joinpath("artifacts", "names_metaweb_tree_gbif.csv")))
````

Read the metaweb

````julia
eurometa = readdlm(joinpath("artifacts", "europeanmetaweb.csv"), ',', String)
mwspecies = unique(eurometa)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(eurometa, 1)
    M[eurometa[i, :]...] = true
end
````

Get the correct rank to cut the matrix at

````julia
rnk = rank(Array(M.edges))
eig = svd(M).S[1:rnk]
neig = eig ./ sum(eig)
````

Plot the eigenvalues

````julia
scatter(eig; lab="", dpi=600, size=(500, 500))
xaxis!("Dimension", (1, 38))
yaxis!("Eigenvalue", (0, 40))
savefig("figures/screeplot.png")
````

Plot the variance explained

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

Partition the metaweb into two latent subspaces, using enough dimensions
to explain 60% of variance

````julia
L, R = rdpg(M, 12)
````

Get the correct threshold

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

Y = tp .+ tn .- 1
maxY, posY = findmax(Y)
threshold = thresholds[posY]

plot(thresholds, Y; lw=0.0, fill=(0, 0.2), lab="", dpi=600, size=(500, 500))
scatter!([threshold], [maxY]; lab="", c=:darkgrey, msw=0.0)
yaxis!("Youden's J", (0, 1))
xaxis!("Cutoff", extrema(L * R))
savefig("figures/optimalcutoff.png")
````

Plot the subspaces

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

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*



````julia
#%% Dependencies
using DelimitedFiles
using DataFrames
using CSV: CSV
using GBIF
using NCBITaxonomy: NCBITaxonomy
using EcologicalNetworks

#%% Read the Strong & Leroux data
sl_raw = readdlm("data/NLfoodweb.csv", ',')

sl_sp = replace.(sl_raw[1, 2:end], "." => " ")
sl_A = Bool.(sl_raw[2:end, 2:end])

#%% Correct the species names from Strong & Leroux to make them match the GBIF taxonomy
nf = NCBITaxonomy.mammalfilter(true);
scinames = Dict{String,String}()
for s in sl_sp
    t = NCBITaxonomy.taxon(s; strict=false)
    if !isnothing(t)
        scinames[s] = t.name
    end
end

#%% Convert the original names in GBIF names through the matched NCBI names
valnames = Dict{String,String}()
for (s, t) in scinames
    gbifmatch = GBIF.taxon(t; strict=false)
    if !isnothing(gbifmatch)
        if !ismissing(gbifmatch.species)
            if gbifmatch.class.first == "Mammalia"
                valnames[s] = gbifmatch.species.first
            end
        end
    end
end

#%% Get the correct keys in the original metaweb
idxmatch = findall(x -> x in keys(valnames), sl_sp)

#%% Assemble the network
spnames = [valnames[s] for s in sl_sp[idxmatch]]
A = sl_A[idxmatch, idxmatch]'
NL = UnipartiteNetwork(A, spnames)

#%% Save as a CSV
df = DataFrame(; from=String[], to=String[])
for i in interactions(NL)
    push!(df, (i.from, i.to))
end
CSV.write("artifacts/newfoundland.csv", df)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*



````julia
#%% Get the dependencies
using HTTP
using ProgressMeter
using JSON
using DataFrames
using CSV: CSV
using StatsPlots

theme(:mute)
default(; frame=:box)

#%% API info
_globi_api = "https://api.globalbioticinteractions.org/taxon"
relevant_types = ["eats", "preysOn"]

#%% Read the corrected Canadian metaweb
canmet = DataFrame(CSV.File("artifacts/canadian_corrected.csv"))

#%% Get the data from GLOBI
allsp = unique(vcat(canmet.from, canmet.to))

diet = DataFrame(; from=String[], to=String[])
@showprogress for sp in allsp
    url = "$(_globi_api)/$(sp)/eats/"
    r = HTTP.request("GET", url)
    globidiet = JSON.parse(String(r.body))
    if !isempty(globidiet["data"])
        for intlist in globidiet["data"]
            if intlist[2] in relevant_types
                for s in intlist[3]
                    if s in allsp
                        if s != sp
                            push!(diet, (sp, s))
                        end
                    end
                end
            end
        end
    end
end
diet = unique(diet)

sort!(diet, [:from, :to])
CSV.write("artifacts/globi_diet.csv", diet)

#%% Get the data we missed from GLOBI
intcode = canmet.from .* canmet.to
diet.intcode = diet.from .* diet.to

missedint = select(diet[findall([!(d in intcode) for d in diet.intcode]), :], [:from, :to])
sort!(missedint, [:from, :to])
CSV.write("artifacts/globi_newinteractions.csv", missedint)

matchedint = dropmissing(leftjoin(diet, canmet; on=[:from => :from, :to => :to]))

@info "GLOBI: found $(size(missedint, 1)) new interactions out of $(size(diet, 1))"
@info "GLOBI: $(length(unique(vcat(diet.from, diet.to)))) species"

#%% Get the data we missed from Strong & Leroux
canmetsp = unique(vcat(canmet.from, canmet.to))
sl = DataFrame(CSV.File("artifacts/newfoundland.csv"))
slshared = intersect(canmetsp, unique(vcat(sl.from, sl.to)))

slkeep = [(i.from in slshared) & (i.to in slshared) for i in eachrow(sl)]
sl = sl[findall(slkeep), :]

sl.intcode = sl.from .* sl.to

missedslint = select(sl[findall([!(d in intcode) for d in sl.intcode]), :], [:from, :to])

sort!(missedslint, [:from, :to])
CSV.write("artifacts/newfoundland_newinteractions.csv", missedslint)

matchedslint = dropmissing(leftjoin(sl, canmet; on=[:from => :from, :to => :to]))

density(matchedint.score; dpi=600, size=(500, 500), lab="GLOBI", lw=2.0)
density!(matchedslint.score; lab="Newfoundland", lw=2.0)
density!(canmet.score; lw=0.0, fill=(0.2, 0), c=:black, lab="All predictions")
xaxis!("Imputed probability", (0, 1))
yaxis!("Density", (0, 8))

savefig("figures/inflation-comparison.png")

@info "NFLD : found $(size(missedslint, 1)) new interactions out of $(size(sl, 1))"
@info "NFLD : $(length(unique(vcat(sl.from, sl.to)))) species"

#%% Merge all interactions
aug = leftjoin(unique(vcat(missedint, missedslint)), canmet; on=[:from, :to])
aug.score .= 1.0

inflated = vcat(canmet, aug)

sort!(inflated, [:score, :from, :to]; rev=[true, false, false])
````

Save the corrected network

````julia
CSV.write("artifacts/canadian_inflated.csv", inflated)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*



````julia
#%% Dependencies
using SparseArrays
using EcologicalNetworks
using DataFrames
using CSV: CSV
using StatsPlots

#%% Update the theme defaults
theme(:mute)
default(; frame=:box)

#%% Load the metaweb
Mij = DataFrame(CSV.File("artifacts/canadian_inflated.csv"))
Si = unique(vcat(Mij.from, Mij.to))

P = UnipartiteProbabilisticNetwork(zeros(Float64, (length(Si), length(Si))), Si)

for int in eachrow(Mij)
    P[int.from, int.to] = int.score
end

#%% Look at the cutoff for kept interactions
ρ = LinRange(extrema(P.edges.nzval)..., 500)
U = zeros(Float64, length(ρ))
L = zeros(Float64, length(ρ))
S = zeros(Float64, length(ρ))

for (i, cutoff) in enumerate(ρ)
    kept = P.edges.nzval .≥ cutoff
    U[i] = sum(kept) / length(kept)
    L[i] = sum(P.edges.nzval[kept]) / links(P)
    S[i] = richness(simplify(P ≥ cutoff)) / richness(P)
end

#%% Central difference as a proxy for second derivative
∂U = zeros(Float64, length(U))
∂L = zeros(Float64, length(U))
for i in 2:(length(U) - 1)
    ∂U[i] = U[i + 1] + U[i - 1] - 2U[i]
    ∂L[i] = L[i + 1] + L[i - 1] - 2L[i]
end

#%% Euclidean distance for U
dU = zeros(Float64, length(U))
for i in eachindex(U)
    dU[i] = sqrt(U[i] * U[i] + ρ[i] * ρ[i])
end

#%% Plot the results
plot(ρ, U; dpi=600, size=(500, 500), lab="Non-zero")
plot!(ρ, L; lab="Expected")
xaxis!("Cutoff", (0, 1))
yaxis!("Proportion of links left", (0, 1))
vline!([ρ[findlast(S .== 1)]]; c=:grey, ls=:dash, lab="")

savefig("figures/cutoff-interactions.png")

plot(ρ, 1.0 .- S; dpi=600, size=(500, 500), lab="Disconnected species", legend=:topleft)
l = L .* links(P)
s = S .* richness(P)
plot!(ρ, l ./ (s .* s); lab="Connectance")
xaxis!("Cutoff", (0, 1))
yaxis!("  ", (0, 0.5))
vline!([ρ[findlast(S .== 1)]]; c=:grey, ls=:dash, lab="")

savefig("figures/cutoff-connectance.png")

#%% Get the threshold
thrind = findlast(S .== 1)
@info "Optimal cutoff based on remaining species: $(ρ[thrind])"
@info "Optimal cutoff based on central differences: $(ρ[last(findmax(∂U))])"

#%% Cleaned-up network
K = copy(P)
K.edges[P.edges .< ρ[thrind]] .= 0.0
dropzeros!(K.edges)

int = DataFrame(; from=String[], to=String[], score=Float64[])
for i in interactions(K)
    push!(int, (i.from, i.to, i.probability))
end
sort!(int, [:score, :from, :to]; rev=[true, false, false])

CSV.write("artifacts/canadian_thresholded.csv", int)

#%% Write the functional classification of species
rls = DataFrame(;
    sp=String[],
    gen=Float64[],
    gen_var=Float64[],
    vul=Float64[],
    vul_var=Float64[],
    role=Symbol[],
)

din = degree(K; dims=2)
dinv = degree_var(K; dims=2)
dout = degree(K; dims=1)
doutv = degree_var(K; dims=1)

for s in species(K)
    rl = :int
    iszero(din[s]) && (rl = :top)
    iszero(dout[s]) && (rl = :prd)
    push!(rls, (s, dout[s], doutv[s], din[s], dinv[s], rl))
end

sort!(rls, :gen; rev=true)
CSV.write("artifacts/species_roles.csv", rls)

plot(log1p.(sort(rls.gen; rev=true)); lab="Out-degree", size=(500, 500), dpi=600)
plot!(log1p.(sort(rls.vul; rev=true)); lab="In-degree")
xaxis!("Rank", (0, 180))
yaxis!("log(degree + 1)", (0, 5))

savefig("figures/final-degree.png")
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

