# Step 8 - finding the interaction cutoff for the final results

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

