using EcologicalNetworks
using DataFrames
import CSV
using DelimitedFiles
using StatsPlots
using StatsBase
using ProgressMeter
using Random
using SparseArrays

theme(:mute)
default(; frame=:box, aspectratio=1, dpi=600)
Random.seed!(01189998819991197253)

suppfig = joinpath("figures", "supplementary")
ispath(suppfig) || mkpath(suppfig)

namelist = DataFrame(CSV.File(joinpath("artifacts", "names_metaweb_tree_gbif.csv")))

# Finally, we grab the European metaweb *not* from the original file, but from
# our edgelist artifact.

eurometa = readdlm(joinpath("artifacts", "europeanmetaweb.csv"), ',', String)
mwspecies = unique(eurometa)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(eurometa, 1)
    M[eurometa[i, :]...] = true
end

ℒ, ℛ = rdpg(M, 12)

links_changed = 10:250:(links(M)-1)

function threshold_network(predicted::Matrix{Float64}, observed::Matrix{Bool})
    thresholds = LinRange(extrema(predicted)..., size(predicted, 1))
    tp = zeros(Float64, length(thresholds))
    tn = similar(tp)
    fp = similar(tp)
    fn = similar(tp)
    for (i, t) in enumerate(thresholds)
        prediction = predicted .>= t
        tp[i] = sum((observed) .& (prediction)) / sum(observed)
        tn[i] = sum((.!(observed)) .& (.!(prediction))) / sum(.!(observed))
        fp[i] = sum((.!(observed)) .& (prediction)) / sum(.!(observed))
        fn[i] = sum((observed) .& (.!(prediction))) / sum(observed)
    end
    n = tp .+ fp .+ tn .+ fn;
    tpr = tp ./ (tp .+ fn);
    fpr = fp ./ (fp .+ tn);
    acc = (tp .+ tn) ./ (n)
    racc = ((tn .+ fp) .* (tn .+ fn) .+ (fn .+ tp) .* (fp .+ tp)) ./ (n .* n);
    J = (tp ./ (tp .+ fn)) + (tn ./ (tn .+ fp)) .- 1.0;
    κ = (acc .- racc) ./ (1.0 .- racc);
    dx = [reverse(fpr)[i] - reverse(fpr)[i - 1] for i in 2:length(fpr)]
    dy = [reverse(tpr)[i] + reverse(tpr)[i - 1] for i in 2:length(tpr)]
    AUC = sum(dx .* (dy ./ 2.0))
    best_at = last(findmax(J))
    return (accuracy=acc[best_at], J=maximum(J), AUC=AUC, κ=κ[best_at], τ=thresholds[best_at])
end

function drop_interactions(network, n)
    m = copy(network)
    pool_to_remove = sample(findall(isequal(true), m.edges), n; replace=false)
    m.edges[pool_to_remove] .= false
    SparseArrays.dropzeros!(m.edges)
    return m
end

function add_interactions(network, n)
    m = copy(network)
    pool_to_add = sample(findall(isequal(false), m.edges), n; replace=false)
    m.edges[pool_to_add] .= true
    return m
end


add_results = []
rem_results = []
for n in links_changed
    m = add_interactions(M, n)
    t = threshold_network(prod(rdpg(m,12)), adjacency(M))
    rm = UnipartiteNetwork(prod(rdpg(m, 12)) .>= t.τ, species(m))
    detected = count(isnothing, indexin(setdiff(interactions(m), interactions(M)), interactions(rm)))
    push!(add_results, merge(t, (detected=detected, changed=n)))
end

for n in links_changed
    m = drop_interactions(M, n)
    t = threshold_network(prod(rdpg(m,12)), adjacency(M))
    rm = UnipartiteNetwork(prod(rdpg(m, 12)) .>= t.τ, species(m))
    detected = count(!isnothing, indexin(setdiff(interactions(M), interactions(m)), interactions(rm)))
    push!(rem_results, merge(t, (detected=detected, changed=n)))
end


plot([0, 1, 1, 0], [0, 0, 0.5, 0.5], st=:shape, c=:grey, alpha=0.1, lw=0.0, aspectratio=1, legend=:bottomright, lab="")
annotate!(0.05, 0.25, text("Worse than random", :black, :left, 8))
plot!([0, 1, 1, 0], [0.5, 0.5, 0.5+1/6, 0.5+1/6], st=:shape, c=:red, alpha=0.1, lw=0.0, lab="")
annotate!(0.05, 0.5+1/12, text("Close to random", :red, :left, 8))
plot!([0, 1, 1, 0], [0.5+1/6, 0.5+1/6, 0.5+2/6, 0.5+2/6], st=:shape, c=:orange, alpha=0.1, lw=0.0, lab="")
annotate!(0.05, 0.5+3/12, text("Fair classifier", :orange, :left, 8))
plot!([0, 1, 1, 0], [0.5+2/6, 0.5+2/6, 1.0, 1.0], st=:shape, c=:green, alpha=0.1, lw=0.0, lab="")
annotate!(0.95, 0.5+5/12, text("Excellent classifier", :green, :right, 8))
scatter!(links_changed./links(M), [r.AUC for r in add_results], lab="Added links", msw=1.0, c=:white)
scatter!(links_changed./links(M), [r.AUC for r in rem_results], lab="Removed links", msw=1.0, c=:white, m=:diamond)
xaxis!("Interactions changed", (0, 1))
yaxis!("ROC-AUC", (0,1))
savefig(joinpath(suppfig, "sensibility_rocauc.png"))

plot([0, 1, 1, 0], [0, 0, 0.5, 0.5], st=:shape, c=:grey, alpha=0.2, lw=0.0, aspectratio=1, lab="", legend=:bottomright)
annotate!(0.05, 0.25, text("Worse than random", :black, :left, 8))
scatter!(links_changed./links(M), [r.accuracy for r in add_results], lab="Added links", msw=1.0, c=:white)
scatter!(links_changed./links(M), [r.accuracy for r in rem_results], lab="Removed links", msw=1.0, c=:white, m=:diamond)
annotate!(0.05, 0.75, text("Better than random", :black, :left, 8))
xaxis!("Interactions changed", (0, 1))
yaxis!("Accuracy", (0,1))
savefig(joinpath(suppfig, "sensibility_accuracy.png"))

plot([0, 1, 1, 0], [0, 0, 0.2, 0.2], st=:shape, c=:grey, alpha=0.5, lw=0.0, aspectratio=1, legend=:bottomright, lab="")
annotate!(0.05, 0.1, text("Poor agreement", :black, :left, 8))
plot!([0, 1, 1, 0], [0.2, 0.2, 0.4, 0.4], st=:shape, c=:grey, alpha=0.4, lw=0.0, lab="")
annotate!(0.05, 0.3, text("Fair agreement", :black, :left, 8))
plot!([0, 1, 1, 0], [0.4, 0.4, 0.6, 0.6], st=:shape, c=:grey, alpha=0.3, lw=0.0, lab="")
annotate!(0.05, 0.5, text("Moderate agreement", :black, :left, 8))
plot!([0, 1, 1, 0], [0.6, 0.6, 0.8, 0.8], st=:shape, c=:grey, alpha=0.2, lw=0.0, lab="")
annotate!(0.05, 0.7, text("Good agreement", :black, :left, 8))
annotate!(0.95, 0.9, text("Very good agreement", :black, :right, 8))
scatter!(links_changed./links(M), [r.κ for r in add_results], lab="Added links", msw=1.0, c=:white)
scatter!(links_changed./links(M), [r.κ for r in rem_results], lab="Removed links", msw=1.0, c=:white, m=:diamond)
xaxis!("Interactions changed", (0, 1))
yaxis!("Cohen's κ", (0,1))
savefig(joinpath(suppfig, "sensibility_kappa.png"))

scatter(links_changed./links(M), [r.τ for r in add_results], lab="Added links", msw=1.0, c=:white)
scatter!(links_changed./links(M), [r.τ for r in rem_results], lab="Removed links", msw=1.0, c=:white, m=:diamond)
hline!([0.206], ls=:dash, c=:grey, lab="")
xaxis!("Interactions changed", (0, 1))
yaxis!("Best threshold", (0,1))
savefig(joinpath(suppfig, "sensibility_threshold.png"))

scatter(links_changed./links(M), [r.detected/r.changed for r in add_results], lab="Added links", msw=1.0, c=:white, legend=:bottomleft)
scatter!(links_changed./links(M), [r.detected/r.changed for r in rem_results], lab="Removed links", msw=1.0, c=:white, m=:diamond)
xaxis!("Interactions changed", (0, 1))
yaxis!("Interactions recovered", (0,1))
savefig(joinpath(suppfig, "sensibility_recovery.png"))

scatter(drop_at./links(M), drop_co, leg=false, msw=0.0, ms=3, c=:black, aspectratio=1)
hline!([connectance(M)], ls=:dash, c=:grey)
xaxis!("Interactions withheld", (0, 1))
yaxis!("Connectance", (0,1))
savefig(joinpath(suppfig, "sensibility_connectance.png"))

scatter(drop_at./links(M), drop_rho, leg=false, msw=0.0, ms=3, c=:black, aspectratio=1)
hline!([ρ(M)], ls=:dash, c=:grey)
xaxis!("Interactions withheld", (0, 1))
yaxis!("Spectral radius", (0,1))
savefig(joinpath(suppfig, "sensibility_spectralradius.png"))

scatter(drop_at./links(M), drop_s4./drop_s5, leg=false, msw=0.0, ms=3, c=:black, aspectratio=1)
s_ratio = length(find_motif(M, unipartitemotifs()[:S4])) / length(find_motif(M, unipartitemotifs()[:S5]))
hline!([s_ratio], ls=:dash, c=:grey)
xaxis!("Interactions withheld", (0, 1))
yaxis!("Competition type ratio", (0,1))
savefig(joinpath(suppfig, "sensibility_motifs.png"))