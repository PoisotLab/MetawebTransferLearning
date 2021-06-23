#%% Dependencies
using EcologicalNetworks
using EcologicalNetworksPlots
using DataFrames
using CSV: CSV
using StatsPlots
using ProgressMeter

#%% Read the final metaweb into a quantitative network at the genus level
mw = DataFrame(CSV.File("artifacts/canadian_thresholded.csv"))

gen = first.(split.(unique(vcat(mw.from, mw.to)), " ")) |> unique
A = zeros(Float64, length(gen), length(gen))

G = UnipartiteQuantitativeNetwork(A, gen)
for i in eachrow(mw)
    gfrom = first(split(i.from, " "))
    gto = first(split(i.to, " "))
    G[gfrom, gto] += i.score
end

nodiagonal!(G)

#%% Plot the network
I = initial(RandomInitialLayout, G)

L = SpringElectric(1.2; gravity=0.2)
L.degree = true
L.δ = 0.2
@showprogress for step in 1:(100richness(G))
  position!(L, I, G)
end

P = UnipartiteProbabilisticNetwork(G.edges ./ maximum(G.edges), EcologicalNetworks._species_objects(G)...)

plot(I, P, size=(500,500), dpi=600)
scatter!(I, G, nodesize=degree(G), nodefill=fractional_trophic_level(convert(UnipartiteNetwork, G)))

