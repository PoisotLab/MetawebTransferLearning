# # Step 8 - finding the interaction cutoff for the final results

# In the final step, we will remove interactions that have a low probability, by
# examining different thresholds - these are refered to as "prediction
# thresholds" in the manuscript.

using SparseArrays
using EcologicalNetworks
using DataFrames
using CSV: CSV
using StatsPlots

theme(:mute)
default(; frame=:box)

# We can load the inflated metaweb, which has all the interactions we predicted
# plus the ones we collected from GLOBI and the Newfoundland dataset (a little
# less than 40 interactions).

Mij = DataFrame(CSV.File("artifacts/canadian_inflated.csv"))
Si = unique(vcat(Mij.from, Mij.to))

P = UnipartiteProbabilisticNetwork(zeros(Float64, (length(Si), length(Si))), Si)

for int in eachrow(Mij)
    P[int.from, int.to] = int.score
end

# ## Finding the prediction threshold

# Our technique for cutoff search will be to get 500 intermediate points between
# the smallest and largest probabilities (resp. 1/n and 1.0), and examine their
# effect on the network. We will specifically look for the number of non-zero
# probability interactions, the expected number of interactions, and the number
# of species with at least one non-zero probability interaction.

ρ = LinRange(extrema(P.edges.nzval)..., 500)
U = zeros(Float64, length(ρ))
L = zeros(Float64, length(ρ))
S = zeros(Float64, length(ρ))

# This is calculated in the following loop:

for (i, cutoff) in enumerate(ρ)
    kept = P.edges.nzval .≥ cutoff
    U[i] = sum(kept) / length(kept)
    L[i] = sum(P.edges.nzval[kept]) / links(P)
    S[i] = richness(simplify(P ≥ cutoff)) / richness(P)
end

# For reference, we have attempted to find a threshold with the central
# difference technique to identify the curvative in the link/species
# relationship, and this gives a threshold that is too low (basically only
# removing the singleton interactions):

∂U = zeros(Float64, length(U))
∂L = zeros(Float64, length(U))
for i in 2:(length(U) - 1)
    ∂U[i] = U[i + 1] + U[i - 1] - 2U[i]
    ∂L[i] = L[i + 1] + L[i - 1] - 2L[i]
end

# Instead, we rely on a visualisation of the relationships, and specifically of
# the point where we remove as many interactions as possible but still keep all
# species connected:

plot(ρ, U; dpi=600, size=(500, 500), lab="Non-zero")
plot!(ρ, L; lab="Expected")
xaxis!("Cutoff", (0, 1))
yaxis!("Proportion of links left", (0, 1))
vline!([ρ[findlast(S .== 1)]]; c=:grey, ls=:dash, lab="")

savefig("figures/cutoff-interactions.png")

# This next plot examines the (lack of an) effect on connectance, which is a
# fairly obvious result, but still interesting to confirm:

plot(ρ, 1.0 .- S; dpi=600, size=(500, 500), lab="Disconnected species", legend=:topleft)
l = L .* links(P)
s = S .* richness(P)
plot!(ρ, l ./ (s .* s); lab="Connectance")
xaxis!("Cutoff", (0, 1))
yaxis!("  ", (0, 0.5))
vline!([ρ[findlast(S .== 1)]]; c=:grey, ls=:dash, lab="")

savefig("figures/cutoff-connectance.png")

# Based on the above, we set the prediction threshold at the point where all
# species remain connected.

thrind = findlast(S .== 1)
@info "Optimal cutoff based on remaining species: $(ρ[thrind])"
@info "Optimal cutoff based on central differences: $(ρ[last(findmax(∂U))])"

# ## Finalizing the network

# We will finally remove all of the thresholded interactions, and this will be
# the final Canadian metaweb.

K = copy(P)
K.edges[P.edges .< ρ[thrind]] .= 0.0
dropzeros!(K.edges)

# We now convert the network into a data frame, which we sort by probability and
# then by species name:

int = DataFrame(; from=String[], to=String[], score=Float64[])
for i in interactions(K)
    push!(int, (i.from, i.to, i.probability))
end
sort!(int, [:score, :from, :to]; rev=[true, false, false])

# Finally, we write the really final Canadian metaweb to a file!

CSV.write("artifacts/canadian_thresholded.csv", int)

# ## Some final cleanup

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
