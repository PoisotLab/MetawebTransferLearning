using EcologicalNetworks
using DataFrames
import CSV
using DelimitedFiles
using StatsPlots
using StatsBase
using ProgressMeter
using Random

theme(:mute)
default(; frame=:box)
Random.seed!(01189998819991197253)

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

drop_at = 10:150:(links(M)-1)
drop_auc = zeros(Float64, length(drop_at))
drop_acc = zeros(Float64, length(drop_at))
drop_J = zeros(Float64, length(drop_at))
drop_thr = zeros(Float64, length(drop_at))
drop_kappa = zeros(Float64, length(drop_at))
drop_recov = zeros(Float64, length(drop_at))
drop_L = zeros(Float64, length(drop_at))
drop_R = zeros(Float64, length(drop_at))


@showprogress for (dridx, dropint) in enumerate(drop_at)
    m = copy(M)
    pool_to_remove = sample(interactions(M), dropint; replace=false)
    for i in pool_to_remove
        m[i.from, i.to] = false
    end
    L, R = rdpg(m, 12)
    A = adjacency(M)
    thresholds = LinRange(extrema(L * R)..., 500)
    tp = zeros(Float64, length(thresholds))
    tn = similar(tp)
    fp = similar(tp)
    fn = similar(tp)
    for (i, t) in enumerate(thresholds)
        PN = (L * R) .>= t
        tp[i] = sum((A) .& (PN)) / sum(A)
        tn[i] = sum((.!(A)) .& (.!(PN))) / sum(.!(A))
        fp[i] = sum((.!(A)) .& (PN)) / sum(.!(A))
        fn[i] = sum((A) .& (.!(PN))) / sum(A)
    end
    n = tp .+ fp .+ tn .+ fn;
    tpr = tp ./ (tp .+ fn);
    fpr = fp ./ (fp .+ tn);
    acc = (tp .+ tn) ./ (n)
    racc = ((tn .+ fp) .* (tn .+ fn) .+ (fn .+ tp) .* (fp .+ tp)) ./ (n .* n);
    J = (tp ./ (tp .+ fn)) + (tn ./ (tn .+ fp)) .- 1.0;
    ppv = tp ./ (tp .+ fp);
    κ = (acc .- racc) ./ (1.0 .- racc);
    dx = [reverse(fpr)[i] - reverse(fpr)[i - 1] for i in 2:length(fpr)]
    dy = [reverse(tpr)[i] + reverse(tpr)[i - 1] for i in 2:length(tpr)]
    AUC = sum(dx .* (dy ./ 2.0))
    drop_auc[dridx] = AUC
    drop_J[dridx] = first(findmax(J))
    drop_kappa[dridx] = first(findmax(κ))
    drop_thr[dridx] = thresholds[last(findmax(J))]
    drop_acc[dridx] = acc[last(findmax(J))]
    drop_L[dridx] =  mean(sqrt.((L .- ℒ).^2.0))
    drop_R[dridx] =  mean(sqrt.((R .- ℛ).^2.0))
    drop_recov[dridx]= count(!isnothing, indexin(pool_to_remove, interactions(UnipartiteQuantitativeNetwork(L*R, species(m)) > t)))/dropint
end

plot([0, 1, 1, 0], [0, 0, 0.5, 0.5], st=:shape, c=:grey, alpha=0.1, lw=0.0, aspectratio=1)
annotate!(0.05, 0.25, text("Worse than random", :black, :left, 8))
plot!([0, 1, 1, 0], [0.5, 0.5, 0.5+1/6, 0.5+1/6], st=:shape, c=:red, alpha=0.1, lw=0.0)
annotate!(0.05, 0.5+1/12, text("Close to random", :red, :left, 8))
plot!([0, 1, 1, 0], [0.5+1/6, 0.5+1/6, 0.5+2/6, 0.5+2/6], st=:shape, c=:orange, alpha=0.1, lw=0.0)
annotate!(0.05, 0.5+3/12, text("Fair classifier", :orange, :left, 8))
plot!([0, 1, 1, 0], [0.5+2/6, 0.5+2/6, 1.0, 1.0], st=:shape, c=:green, alpha=0.1, lw=0.0)
annotate!(0.95, 0.5+5/12, text("Excellent classifier", :green, :right, 8))
scatter!(drop_at./links(M), drop_auc, lab="", leg=false, msw=2.0, msc=:black, c=:white)
xaxis!("Interactions withheld", (0, 1))
yaxis!("ROC-AUC", (0,1))
savefig("figures/sensibility_rocauc.png")

plot([0, 1, 1, 0], [0, 0, 0.5, 0.5], st=:shape, c=:grey, alpha=0.2, lw=0.0, aspectratio=1)
scatter!(drop_at./links(M), drop_acc, lab="", leg=false, msw=2.0, msc=:black, c=:white)
xaxis!("Interactions withheld", (0, 1))
yaxis!("Accuracy", (0,1))
savefig("figures/sensibility_accuracy.png")

plot([0, 1, 1, 0], [0, 0, 0.2, 0.2], st=:shape, c=:grey, alpha=0.5, lw=0.0, aspectratio=1)
annotate!(0.05, 0.1, text("Poor agreement", :black, :left, 8))
plot!([0, 1, 1, 0], [0.2, 0.2, 0.4, 0.4], st=:shape, c=:grey, alpha=0.4, lw=0.0)
annotate!(0.05, 0.3, text("Fair agreement", :black, :left, 8))
plot!([0, 1, 1, 0], [0.4, 0.4, 0.6, 0.6], st=:shape, c=:grey, alpha=0.3, lw=0.0)
annotate!(0.05, 0.5, text("Moderate agreement", :black, :left, 8))
plot!([0, 1, 1, 0], [0.6, 0.6, 0.8, 0.8], st=:shape, c=:grey, alpha=0.2, lw=0.0)
annotate!(0.05, 0.7, text("Good agreement", :black, :left, 8))
annotate!(0.95, 0.9, text("Very good agreement", :black, :right, 8))
scatter!(drop_at./links(M), drop_kappa, lab="", leg=false, msw=2.0, msc=:black, c=:white)
xaxis!("Interactions withheld", (0, 1))
yaxis!("Cohen's κ", (0,1))
savefig("figures/sensibility_kappa.png")

scatter(drop_at./links(M), drop_thr, leg=false, msw=2.0, msc=:black, c=:white)
hline!([0.206], ls=:dash, c=:grey)
xaxis!("Interactions withheld", (0, 1))
yaxis!("Best threshold", (0,1))
savefig("figures/sensibility_threshold.png")
