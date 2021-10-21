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


links_changed = 100:100:(round(Int, 2342/100)*100)

add_results = []
rem_results = []
mix_results = []

for n in links_changed
    m = add_interactions(M, n)
    t = threshold_network(prod(rdpg(m,12)), adjacency(M))
    rm = UnipartiteNetwork(prod(rdpg(m, 12)) .>= t.τ, species(m))
    detected = count(isnothing, indexin(setdiff(interactions(m), interactions(M)), interactions(rm)))
    push!(add_results, merge(t, (detected=detected, changed=n, network=rm)))
end

for n in links_changed
    m = drop_interactions(M, n)
    t = threshold_network(prod(rdpg(m,12)), adjacency(M))
    rm = UnipartiteNetwork(prod(rdpg(m, 12)) .>= t.τ, species(m))
    detected = count(!isnothing, indexin(setdiff(interactions(M), interactions(m)), interactions(rm)))
    push!(rem_results, merge(t, (detected=detected, changed=n, network=rm)))
end

for n in links_changed
    m1 = drop_interactions(M, round(Int, n/2))
    m = add_interactions(m1, round(Int, n/2))
    t = threshold_network(prod(rdpg(m,12)), adjacency(M))
    rm = UnipartiteNetwork(prod(rdpg(m, 12)) .>= t.τ, species(m))
    detected_add = count(isnothing, indexin(setdiff(interactions(m), interactions(M)), interactions(rm)))
    detected_rem = count(!isnothing, indexin(setdiff(interactions(M), interactions(m)), interactions(rm)))
    push!(mix_results, merge(t, (detected=detected_add+detected_rem, changed=n, network=rm)))
end

plot([0, 1, 1, 0], [0, 0, 0.5, 0.5], st=:shape, c=:grey, alpha=0.1, lw=0.0, aspectratio=1, legend=:bottomright, lab="")
annotate!(0.05, 0.25, text("Worse than random", :black, :left, 8))
plot!([0, 1, 1, 0], [0.5, 0.5, 0.5+1/6, 0.5+1/6], st=:shape, c=:red, alpha=0.1, lw=0.0, lab="")
annotate!(0.05, 0.5+1/12, text("Close to random", :red, :left, 8))
plot!([0, 1, 1, 0], [0.5+1/6, 0.5+1/6, 0.5+2/6, 0.5+2/6], st=:shape, c=:orange, alpha=0.1, lw=0.0, lab="")
annotate!(0.05, 0.5+3/12, text("Fair classifier", :orange, :left, 8))
plot!([0, 1, 1, 0], [0.5+2/6, 0.5+2/6, 1.0, 1.0], st=:shape, c=:green, alpha=0.1, lw=0.0, lab="")
annotate!(0.05, 0.5+5/12, text("Excellent class.", :green, :left, 8))
scatter!(links_changed./links(M), [r.AUC for r in add_results], lab="Added", msw=1.0, c=:darkgrey, m=:utriangle)
scatter!(links_changed./links(M), [r.AUC for r in mix_results], lab="Both", msw=1.0, c=:white, m=:circle)
scatter!(links_changed./links(M), [r.AUC for r in rem_results], lab="Removed", msw=1.0, c=:lightgrey, m=:dtriangle)
xaxis!("Interactions changed", (0, 1))
yaxis!("ROC-AUC", (0,1))
savefig(joinpath(suppfig, "sensibility_rocauc.png"))

plot([0, 1, 1, 0], [0, 0, 0.5, 0.5], st=:shape, c=:grey, alpha=0.2, lw=0.0, aspectratio=1, lab="", legend=:bottomright)
annotate!(0.05, 0.25, text("Worse than random", :black, :left, 8))
scatter!(links_changed./links(M), [r.accuracy for r in add_results], lab="Added", msw=1.0, c=:darkgrey, m=:utriangle)
scatter!(links_changed./links(M), [r.accuracy for r in mix_results], lab="Both", msw=1.0, c=:white, m=:circle)
scatter!(links_changed./links(M), [r.accuracy for r in rem_results], lab="Removed", msw=1.0, c=:lightgrey, m=:dtriangle)
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
scatter!(links_changed./links(M), [r.κ for r in add_results], lab="Added", msw=1.0, c=:darkgrey, m=:utriangle)
scatter!(links_changed./links(M), [r.κ for r in mix_results], lab="Both", msw=1.0, c=:white, m=:circle)
scatter!(links_changed./links(M), [r.κ for r in rem_results], lab="Removed", msw=1.0, c=:lightgrey, m=:dtriangle)
xaxis!("Interactions changed", (0, 1))
yaxis!("Cohen's κ", (0,1))
savefig(joinpath(suppfig, "sensibility_kappa.png"))

scatter(links_changed./links(M), [r.τ for r in add_results], lab="Added", msw=1.0, c=:darkgrey, m=:utriangle)
scatter!(links_changed./links(M), [r.τ for r in mix_results], lab="Both", msw=1.0, c=:white, m=:circle)
scatter!(links_changed./links(M), [r.τ for r in rem_results], lab="Removed", msw=1.0, c=:lightgrey, m=:dtriangle)
hline!([0.206], ls=:dash, c=:grey, lab="")
xaxis!("Interactions changed", (0, 1))
yaxis!("Best threshold", (0,1))
savefig(joinpath(suppfig, "sensibility_threshold.png"))

scatter(links_changed./links(M), [r.detected/r.changed for r in add_results], lab="Added", msw=1.0, c=:darkgrey, m=:utriangle, legend=:bottomleft)
scatter!(links_changed./links(M), [r.detected/r.changed for r in mix_results], lab="Both", msw=1.0, c=:white, m=:circle)
scatter!(links_changed./links(M), [r.detected/r.changed for r in rem_results], lab="Removed", msw=1.0, c=:lightgrey, m=:dtriangle)
xaxis!("Interactions changed", (0, 1))
yaxis!("Changes corrected", (0,1))
savefig(joinpath(suppfig, "sensibility_recovery.png"))

scatter(links_changed./links(M), [connectance(r.network) for r in add_results], lab="Added", msw=1.0, c=:darkgrey, m=:utriangle, legend=:topright)
scatter!(links_changed./links(M), [connectance(r.network) for r in mix_results], lab="Both", msw=1.0, c=:white, m=:circle)
scatter!(links_changed./links(M), [connectance(r.network) for r in rem_results], lab="Removed", msw=1.0, c=:lightgrey, m=:dtriangle)
hline!([connectance(M)], ls=:dash, c=:grey, lab="")
xaxis!("Interactions changed", (0, 1))
yaxis!("Connectance", (0,1))
savefig(joinpath(suppfig, "sensibility_connectance.png"))

scatter(links_changed./links(M), [ρ(r.network) for r in add_results], lab="Added", msw=1.0, c=:darkgrey, m=:utriangle, legend=:topright)
scatter!(links_changed./links(M), [ρ(r.network) for r in mix_results], lab="Both", msw=1.0, c=:white, m=:circle)
scatter!(links_changed./links(M), [ρ(r.network) for r in rem_results], lab="Removed", msw=1.0, c=:lightgrey, m=:dtriangle)
hline!([ρ(M)], ls=:dash, c=:grey, lab="")
xaxis!("Interactions changed", (0, 1))
yaxis!("Spectral radius", (0,1))
savefig(joinpath(suppfig, "sensibility_spectralradius.png"))

m1, m2 = unipartitemotifs()[:S4], unipartitemotifs()[:S5]

function motifratio(network, m1, m2)
    return length(find_motif(network, m1)) / length(find_motif(network, m2))
end

function motifratio(network)
    return motifratio(network, m1, m2)
end

scatter(links_changed./links(M), [motifratio(r.network) for r in add_results], lab="Added", msw=1.0, c=:darkgrey, m=:utriangle, legend=:topright)
scatter!(links_changed./links(M), [motifratio(r.network) for r in mix_results], lab="Both", msw=1.0, c=:white, m=:circle)
scatter!(links_changed./links(M), [motifratio(r.network) for r in rem_results], lab="Removed", msw=1.0, c=:lightgrey, m=:dtriangle)
hline!([motifratio(M)], ls=:dash, c=:grey, lab="")
xaxis!("Interactions changed", (0, 1))
yaxis!("Competition type ratio", (0,1))
savefig(joinpath(suppfig, "sensibility_motifs.png"))

# What if we use a subset of the network to estimate the threshold? - Part 2 of the supp mat

species_used = 20:10:richness(M)

ind_results = []

for n in species_used
    nodes_kept = sample(species(M), n; replace=false)
    m = M[nodes_kept]
    t = threshold_network(prod(rdpg(m,12)), adjacency(m))
    prediction = prod(rdpg(M, 12)) .>= t.τ
    observed = adjacency(M)
    tp = sum((observed) .& (prediction)) / sum(observed)
    tn = sum((.!(observed)) .& (.!(prediction))) / sum(.!(observed))
    fp = sum((.!(observed)) .& (prediction)) / sum(.!(observed))
    fn = sum((observed) .& (.!(prediction))) / sum(observed)
    J = (tp ./ (tp .+ fn)) + (tn ./ (tn .+ fp)) .- 1.0;
    push!(ind_results, (changed=n, thr=t.τ, tp=tp, fp=fp, tn=tn, fn=fn, J=J))
end

p1 = scatter(species_used./richness(M), [r.tp for r in ind_results], lab="", c=:darkgrey, msw=0.0, ylab="True positives")
p2 = scatter(species_used./richness(M), [r.tn for r in ind_results], lab="", c=:darkgrey, msw=0.0, ylab="True negatives")
p3 = scatter(species_used./richness(M), [r.fp for r in ind_results], lab="", c=:darkgrey, msw=0.0, ylab="False positives")
p4 = scatter(species_used./richness(M), [r.fn for r in ind_results], lab="", c=:darkgrey, msw=0.0, ylab="False negatives")

for p in [p1, p2, p3, p4]
    xaxis!(p, (0, 1), "Relative richness")
    yaxis!((0, 1))
end

plot(p1, p2, p3, p4)
savefig(joinpath(suppfig, "sensibility_species.png"))