using EcologicalNetworks
using Plots
using MultivariateStats
using DataFrames
using CSV: CSV
using DelimitedFiles
using Random
using StatsBase
using Statistics
using LinearAlgebra

eurometa = readdlm(joinpath("artifacts", "europeanmetaweb.csv"), ',', String)
mwspecies = unique(eurometa)
EUR = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(eurometa, 1)
    EUR[eurometa[i, :]...] = true
end

Mij = DataFrame(CSV.File("artifacts/canadian_thresholded.csv"))
Si = unique(vcat(Mij.from, Mij.to))
CAN = UnipartiteProbabilisticNetwork(zeros(Float64, (length(Si), length(Si))), Si)
for int in eachrow(Mij)
    CAN[int.from, int.to] = int.score
end

Mij = DataFrame(CSV.File("artifacts/all_mean_subset.csv"))
Si = unique(vcat(Mij.from, Mij.to))
ALL = UnipartiteNetwork(zeros(Bool, (length(Si), length(Si))), Si)
for int in eachrow(Mij)
    ALL[int.from, int.to] = true
end

# Sub-sample
CANs = [
    simplify(nodiagonal(CAN[sample(species(CAN), rand(30:70); replace=false)])) for
    i in 1:400
]
EURs = [
    simplify(nodiagonal(EUR[sample(species(EUR), rand(30:70); replace=false)])) for
    i in 1:400
]
ALLs = [
    simplify(nodiagonal(ALL[sample(species(ALL), rand(30:70); replace=false)])) for
    i in 1:400
]

filter!(n -> (richness(n)>3)&(links(n)>3), CANs)
filter!(n -> (richness(n)>3)&(links(n)>3), EURs)
filter!(n -> (richness(n)>3)&(links(n)>3), ALLs)

# Functions
function tib(n)
    dout = collect(values(degree(n; dims=1)))
    din = collect(values(degree(n; dims=2)))
    top = count(iszero, din)
    bottom = count(iszero, dout)
    inter = richness(n) - (top + bottom)
    return [top, inter, bottom] ./ richness(n)
end

# Get vec
function pvec(N::UnipartiteNetwork)
    try
        tl = collect(values(trophic_level(N)))
        return [
            mean(tl),
            std(tl),
            std(collect(values(degree(N)))),
            std(collect(values(degree(N; dims=1)))),
            std(collect(values(degree(N; dims=2)))),
            tib(N)...,
        ]
    catch
        return nothing
    end
end

pvec(N::UnipartiteProbabilisticNetwork) = pvec(rand(N))

# Get data

EURv = hcat(pvec.(EURs)...)
CANv = hcat(pvec.(CANs)...)
ALLv = hcat(pvec.(ALLs)...)

Y = hcat(CANv, EURv, ALLv)
ny = size(Y, 2)

# PCA

W = fit(
    PCA, Y[:, sample(1:ny, 300; replace=false)]
)

PCAe = MultivariateStats.transform(W, EURv)
PCAc = MultivariateStats.transform(W, CANv)
PCAa = MultivariateStats.transform(W, ALLv)

function ellipse(xy)
    μ = mean(xy; dims=2)
    C = cov(xy')
    λ = eigvals(C)
    Λ = sqrt.(5.991 .* λ)
    𝐯 = eigvecs(C)
    α = atan(𝐯[2, 1] / 𝐯[1, 1])
    t = LinRange(0.0, 2π, 200)

    x = @. μ[1] + Λ[1] * cos(t) * cos(α) - Λ[2] * sin(t) * sin(α)
    y = @. μ[2] + Λ[2] * sin(t) * cos(α) + Λ[1] * cos(t) * sin(α)
    return Shape(x, y)
end

plot(; frame=:origin, dpi=600, size=(500, 500), legend=:bottomleft)
plot!(ellipse(PCAe[1:2, :]); lw=2.5, alpha=0.05, c=:orange, lc=:orange, lab="")
scatter!(PCAe[1, :], PCAe[2, :]; c=:orange, msw=0.0, ms=2, lab="European sub-samples")
plot!(ellipse(PCAc[1:2, :]); lw=2.5, alpha=0.05, c=:teal, lc=:teal, lab="")
scatter!(PCAc[1, :], PCAc[2, :]; c=:teal, msw=0.0, ms=2, lab="Canadian sub-samples")
xaxis!("Princ. Comp. 1 ($(round(Int, (W.prinvars[1]/W.tvar)*100))%)")
yaxis!("Princ. Comp. 2 ($(round(Int, (W.prinvars[2]/W.tvar)*100))%)")
savefig(joinpath("figures", "supplementary", "variation_pca.png"))

res = DataFrame(;
    var=["meanTL", "stdTL", "stdK", "stdKout", "stdKin", "T", "I", "B"],
    pc1=W.proj[:, 1],
    pc2=W.proj[:, 2]
)
sort!(res, :pc1)
CSV.write("artifacts/supplementary_pca.csv")