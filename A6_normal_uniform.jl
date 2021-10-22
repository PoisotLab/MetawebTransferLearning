using EcologicalNetworks
using DataFrames
import CSV
using DelimitedFiles
using StatsPlots
using StatsBase
using ProgressMeter
using Random
using SparseArrays
using Statistics

theme(:mute)
default(; frame=:box, aspectratio=1, dpi=600, size=(500, 500))
Random.seed!(01189998819991197253)

Mij = DataFrame(CSV.File("artifacts/canadian_uncorrected.csv"))
Si = unique(vcat(Mij.from, Mij.to))
U = UnipartiteProbabilisticNetwork(zeros(Float64, (length(Si), length(Si))), Si)
N = UnipartiteProbabilisticNetwork(zeros(Float64, (length(Si), length(Si))), Si)
for int in eachrow(Mij)
    U[int.from, int.to] = int.score
    N[int.from, int.to] = int.score_normal
end

eurometa = readdlm(joinpath("artifacts", "europeanmetaweb.csv"), ',', String)
mwspecies = unique(eurometa)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(eurometa, 1)
    M[eurometa[i, :]...] = true
end


histogram2d(Mij.score, Mij.score_normal, bins=20, c=:batlow)
xaxis!("U(m,M)", (0,1))
yaxis!("N(μ,σ)", (0,1))
savefig(joinpath("figures", "supplementary", "comparison_models.png"))

density(Mij.score, aspectratio=1/130)
density!(Mij.score_normal)
xaxis!("Probability", (0, 1))
yaxis!("Density", (0, 130))

function tib(n)
    dout = collect(values(degree(n; dims=1)))
    din = collect(values(degree(n; dims=2)))
    top = count(iszero, din)
    bottom = count(iszero, dout)
    inter = richness(n) - (top + bottom)
    return [top, inter, bottom] ./ richness(n)
end

thr = LinRange(0.0, 1.0, 150)
lu = [connectance(U>=t) for t in thr]
ln = [connectance(N>=t) for t in thr]
plot(thr, lu, lab="Uniform model", ls=:dash, lw=2.0)
plot!(thr, ln, lab="Normal model", ls=:dot, lw=2.0)
yaxis!((0, 1), "Threshold")
xaxis!((0, 1), "Connectance")
savefig(joinpath("figures", "supplementary", "comparison_connectance.png"))

tibu = [tib(simplify(U>=t)) for t in thr]
tibn = [tib(simplify(N>=t)) for t in thr]

p1 = plot(thr, rotl90(hcat(tibn...))[:,1], lab="Apex")
plot!(p1, thr, rotl90(hcat(tibn...))[:,2], lab="Intermediate")
plot!(p1, thr, rotl90(hcat(tibn...))[:,3], lab="Basal")
xaxis!(p1, (0, 1), "Threshold")
yaxis!(p1, (0, 1), "Proportion")
title!(p1, "Normal model")

p2 = plot(thr, rotl90(hcat(tibu...))[:,1], lab="Apex")
plot!(p2, thr, rotl90(hcat(tibu...))[:,2], lab="Intermediate")
plot!(p2, thr, rotl90(hcat(tibu...))[:,3], lab="Basal")
xaxis!(p2, (0, 1), "Threshold")
yaxis!(p2, (0, 1), "Proportion")
title!(p2, "Uniform model")

plot(p1, p2, size=(800,550))
savefig(joinpath("figures", "supplementary", "comparison_tib.png"))

function TL(n)
    try
        tl = collect(values(trophic_level(n)))
        return [mean(tl), std(tl), maximum(tl)]
    catch
        return [NaN, NaN, NaN]
    end
end

tlu = [TL(simplify(U>=t)) for t in thr]
tln = [TL(simplify(N>=t)) for t in thr]

p1 = plot(thr, rotl90(hcat(tlu...))[:,1], ribbon=rotl90(hcat(tlu...))[:,2], lab="Average", aspectratio=1/5)
plot!(p1, thr, rotl90(hcat(tlu...))[:,3], lab="Maximum", c=:black, ls=:dash)
xaxis!(p1, (0,1), "Threshold")
yaxis!(p1, (0,5), "Trophic chain length")
title!(p1, "Uniform model")

p2 = plot(thr, rotl90(hcat(tln...))[:,1], ribbon=rotl90(hcat(tln...))[:,2], lab="Average", aspectratio=1/5)
plot!(p2, thr, rotl90(hcat(tln...))[:,3], lab="Maximum", c=:black, ls=:dash)
xaxis!(p2, (0,1), "Threshold")
yaxis!(p2, (0,5), "Trophic chain length")
title!(p2, "Normal model")

plot(p2, p1, size=(800,550))
savefig(joinpath("figures", "supplementary", "comparison_rophicchain.png"))
