using DelimitedFiles
using ProgressMeter
using Statistics
using Plots
using EcologicalNetworks
using EcologicalNetworksPlots
using DataFrames
using CSV

eurometa = readdlm(joinpath("artifacts", "europeanmetaweb.csv"), ',', String)
mwspecies = unique(eurometa)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(eurometa, 1)
    M[eurometa[i, :]...] = true
end

canmeta = DataFrame(CSV.File("artifacts/canadian_inflated.csv"))
mwspecies = unique(vcat(canmeta.from, canmeta.to))
P = UnipartiteProbabilisticNetwork(
    zeros(Float64, length(mwspecies), length(mwspecies)), mwspecies
)
for int in eachrow(canmeta)
    P[int.from, int.to] = int.score
end

extrema(unique(Array(P.edges)))
thr = collect(-3:0.1:0.0)
lnk = zeros(Float64, length(thr))
ric = zeros(Float64, length(thr))
for (i,t) in enumerate(thr)
    this = simplify(P > Float64(10.0^t))
    lnk[i] = links(this)
    ric[i] = richness(this)
end

# Net stats
netstats = DataFrame(;
    species=String[],
    omnivory=Float64[],
    tl=Float64[],
)
@showprogress for rep in 1:100
    R = simplify(rand(P))
    om = omnivory(R)
    tl = trophic_level(R)
    for s in species(R)
        push!(netstats, (s, om[s], tl[s]))
    end
end

gf = groupby(netstats, :species)
avg = combine(gf, :omnivory => mean => :omnivory, :tl => mean => :trophiclevel)

scatter(avg.omnivory, avg.trophiclevel)

sort(netstats, :tl; rev=true)[1:10, :]

I = initial(UnravelledInitialLayout, P>0.0)

function random_omnivory(N::T) where {T<:UnipartiteNetwork}
    o = omnivory(N)
    for s in species(N)
        o[s] += (rand() - 0.5) * 0.1
    end
    return o
end

UL = UnravelledLayout(; x=random_omnivory, y=trophic_level)

position!(UL, I, P>0.0)

plot(I, P; lab="", framestyle=:box)

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
