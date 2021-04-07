using Phylo
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

# Read the tree
tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]

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

A = adjacency(M)
thresholds = LinRange(extrema(L*R)..., 100)
tp = zeros(Float64, length(thresholds))
tn = similar(tp)
fp = similar(tp)
fn = similar(tp)
for (i,t) in enumerate(thresholds)
    PN = (L*R).>= t
    tp[i] = sum((A).&(PN))/sum(A)
    tn[i] = sum((.!(A)).&(.!(PN)))/sum(.!(A))
    fp[i] = sum((.!(A)).&(PN))/sum(.!(A))
    fn[i] = sum((A).&(.!(PN)))/sum(A)
end

Y = tp .+ tn .- 1
threshold = thresholds[findmax(Y)[2]]

plot(fp, tp, lab="", aspectratio=1, xlim=(0,1), fill=(0, 0.2), ylim=(0,1), dpi=600)
xaxis!("False positives")
yaxis!("True positives")
savefig("figures/roc.png")

heatmap(L*R)
heatmap((L*R).>threshold)

# Get a random species pool from the tree
treenodes = [n.name for n in tree.nodes if !startswith(n.name, "Node ")]
canada = DataFrame(CSV.File(joinpath("artifacts", "iucn_gbif_names.csv")))
cancodes = replace.(unique(filter(!ismissing, canada.gbifname)), " " => "_")
tree_cleanup = DataFrame(CSV.File(joinpath("artifacts", "upham_gbif_names.csv")))

csp = dropmissing!(DataFrame(gbifname = canada.gbifname))
csp = dropmissing(leftjoin(csp, tree_cleanup, on=:gbifname))

canmammals = unique(csp.code)

k = 3
pool = canmammals
poolsize = length(pool)
l = zeros(Float64, poolsize, size(L, 2))
r = permutedims(l)
prog = Progress(length(pool))
Threads.@threads for i in 1:length(pool)
    target = pool[i]
    D = [distance(tree, target, first(namelist[isequal(n).(namelist.name),:upham])) for n in species(M)]
    p = sortperm(D)
    if iszero(minimum(D))
        r[:,i] = R[:,p][:,1]
        l[i,:] = L[p,:][1,:]
    else
        r[:,i] = mean(R[:,p][:,1:k]; dims=2)
        l[i,:] = mean(L[p,:][1:k,:]; dims=1)
    end
    next!(prog)
end

P = simplify(UnipartiteQuantitativeNetwork(l*r, pool))
adjacency(P) |> heatmap

s_int = sort(filter(x -> x.strength >= threshold, interactions(P)), by = (x) -> x.strength, rev=true)
predint = DataFrame(predator = String[], prey = String[], evidence = Float64[])
canpool = unique(csp.gbifname)
N = UnipartiteNetwork(zeros(Bool, length(canpool), length(canpool)), canpool)
for i in s_int
    fr = first(unique(csp[csp.code .==i.from, :gbifname]))
    to = first(unique(csp[csp.code .==i.to, :gbifname]))
    push!(predint, (fr, to, i.strength))
    N[fr,to] = true
end
simplify!(N)

CSV.write(joinpath("artifacts", "mammals_metaweb.csv"), predint)

adj = N |> adjacency
spy(adj, frame=:grid, xlab="Preys", ylab="Predators")
savefig(joinpath("figures", "canmammals.png"))

# Net stats

netstats = DataFrame(
    species = String[],
    generality = Int64[],
    vulnerability = Int64[],
    omnivory = Float64[],
    tl = Float64[]
)

din = degree(N; dims=2)
dot = degree(N; dims=1)
omn = omnivory(N)
trl = trophic_level(N)

for s in species(N)
    push!(
        netstats,
        (
            s,
            dot[s],
            din[s],
            omn[s],
            trl[s]
        )
    )
end

sort(netstats, :omnivory, rev=true)[1:10,:]

#

I = initial(UnravelledInitialLayout, N)

function random_omnivory(N::T) where {T <: UnipartiteNetwork}
  o = omnivory(N)
  for s in species(N)
    o[s] += (rand()-0.5)*0.1
  end
  return o
end

UL = UnravelledLayout(x=random_omnivory, y=trophic_level)
position!(UL, I, N)

plot(I, N, lab="", framestyle=:box)
scatter!(I, N, nodefill=degree(N), colorbar=true, framestyle=:box, dpi=600)
xaxis!("Omnivory")
yaxis!("Trophic level")
savefig("figures/omni-can.png")


J = initial(UnravelledInitialLayout, M)


UL = UnravelledLayout(x=random_omnivory, y=trophic_level)
position!(UL, J, M)

plot(J, M, lab="", framestyle=:box)
scatter!(J, M, nodefill=degree(M), colorbar=true, framestyle=:box, dpi=600)
xaxis!("Omnivory")
yaxis!("Trophic level")
savefig("figures/omni-eur.png")
