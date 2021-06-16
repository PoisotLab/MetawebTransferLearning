using DelimitedFiles
using Plots
using EcologicalNetworks
using EcologicalNetworksPlots
using DataFrames

eurometa = readdlm(joinpath("artifacts", "europeanmetaweb.csv"), ',', String)
mwspecies = unique(eurometa)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(eurometa, 1)
    M[eurometa[i, :]...] = true
end

canmeta = readdlm(joinpath("artifacts", "canadianmetaweb.csv"), ',', String)
mwspecies = unique(canmeta)
N = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(canmeta, 1)
    N[canmeta[i, :]...] = true
end

# Net stats
netstats = DataFrame(;
    species=String[],
    generality=Int64[],
    vulnerability=Int64[],
    omnivory=Float64[],
    tl=Float64[],
)

din = degree(N; dims=2)
dot = degree(N; dims=1)
omn = omnivory(N)
trl = trophic_level(N)

for s in species(N)
    push!(netstats, (s, dot[s], din[s], omn[s], trl[s]))
end

sort(netstats, :tl; rev=true)[1:10, :]

#

I = initial(UnravelledInitialLayout, N)

function random_omnivory(N::T) where {T<:UnipartiteNetwork}
    o = omnivory(N)
    for s in species(N)
        o[s] += (rand() - 0.5) * 0.1
    end
    return o
end

UL = UnravelledLayout(; x=random_omnivory, y=trophic_level)
position!(UL, I, N)

plot(I, N; lab="", framestyle=:box)
scatter!(I, N; nodefill=degree(N), colorbar=true, framestyle=:box, dpi=600)
xaxis!("Omnivory")
yaxis!("Trophic level")
savefig("figures/omni-can.png")

J = initial(UnravelledInitialLayout, M)

UL = UnravelledLayout(; x=random_omnivory, y=trophic_level)
position!(UL, J, M)

plot(J, M; lab="", framestyle=:box)
scatter!(J, M; nodefill=degree(M), colorbar=true, framestyle=:box, dpi=600)
xaxis!("Omnivory")
yaxis!("Trophic level")
savefig("figures/omni-eur.png")
