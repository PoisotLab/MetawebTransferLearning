using EcologicalNetworks
using DataFrames
import CSV
using DelimitedFiles
using StatsPlots
using StatsBase
using ProgressMeter
using Random

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

drop_at = 10:250:(links(M)-1)
drop_auc = zeros(Float64, length(drop_at))
drop_acc = zeros(Float64, length(drop_at))
drop_racc = zeros(Float64, length(drop_at))
drop_bacc = zeros(Float64, length(drop_at))
drop_J = zeros(Float64, length(drop_at))
drop_thr = zeros(Float64, length(drop_at))
drop_kappa = zeros(Float64, length(drop_at))
drop_recov = zeros(Float64, length(drop_at))
drop_L = zeros(Float64, length(drop_at))
drop_R = zeros(Float64, length(drop_at))
drop_co = zeros(Float64, length(drop_at))
drop_rho = zeros(Float64, length(drop_at))
drop_s4 = zeros(Float64, length(drop_at))
drop_s5 = zeros(Float64, length(drop_at))


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
    bacc = ((tp ./ (tp .+ fn)) .+ (tn ./ (fp .+ tn))) ./ 2.0
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
    drop_racc[dridx] = racc[last(findmax(J))]
    drop_bacc[dridx] = bacc[last(findmax(J))]
    drop_L[dridx] =  mean(sqrt.((L .- ℒ).^2.0))
    drop_R[dridx] =  mean(sqrt.((R .- ℛ).^2.0))
    N = UnipartiteQuantitativeNetwork(L*R, species(m)) > thresholds[last(findmax(J))]
    drop_recov[dridx]= count(!isnothing, indexin(pool_to_remove, interactions(N)))/dropint
    drop_co[dridx] = connectance(N)
    drop_rho[dridx] = ρ(N)
    drop_s4[dridx] = length(find_motif(N, unipartitemotifs()[:S4]))
    drop_s5[dridx] = length(find_motif(N, unipartitemotifs()[:S5]))
end

plot([0, 1, 1, 0], [0, 0, 0.5, 0.5], st=:shape, c=:grey, alpha=0.1, lw=0.0, aspectratio=1)
annotate!(0.05, 0.25, text("Worse than random", :black, :left, 8))
plot!([0, 1, 1, 0], [0.5, 0.5, 0.5+1/6, 0.5+1/6], st=:shape, c=:red, alpha=0.1, lw=0.0)
annotate!(0.05, 0.5+1/12, text("Close to random", :red, :left, 8))
plot!([0, 1, 1, 0], [0.5+1/6, 0.5+1/6, 0.5+2/6, 0.5+2/6], st=:shape, c=:orange, alpha=0.1, lw=0.0)
annotate!(0.05, 0.5+3/12, text("Fair classifier", :orange, :left, 8))
plot!([0, 1, 1, 0], [0.5+2/6, 0.5+2/6, 1.0, 1.0], st=:shape, c=:green, alpha=0.1, lw=0.0)
annotate!(0.95, 0.5+5/12, text("Excellent classifier", :green, :right, 8))
scatter!(drop_at./links(M), drop_auc, lab="", leg=false, msw=0.0, ms=3, c=:black)
xaxis!("Interactions withheld", (0, 1))
yaxis!("ROC-AUC", (0,1))
savefig(joinpath(suppfig, "sensibility_rocauc.png"))

plot([0, 1, 1, 0], [0, 0, 0.5, 0.5], st=:shape, c=:grey, alpha=0.2, lw=0.0, aspectratio=1)
annotate!(0.05, 0.25, text("Worse than random", :black, :left, 8))
scatter!(drop_at./links(M), drop_bacc, lab="", leg=false, msw=0.0, ms=3, c=:black)
annotate!(0.05, 0.75, text("Better than random", :black, :left, 8))
xaxis!("Interactions withheld", (0, 1))
yaxis!("Accuracy", (0,1))
savefig(joinpath(suppfig, "sensibility_accuracy.png"))

plot([0, 1, 1, 0], [0, 0, 0.2, 0.2], st=:shape, c=:grey, alpha=0.5, lw=0.0, aspectratio=1)
annotate!(0.05, 0.1, text("Poor agreement", :black, :left, 8))
plot!([0, 1, 1, 0], [0.2, 0.2, 0.4, 0.4], st=:shape, c=:grey, alpha=0.4, lw=0.0)
annotate!(0.05, 0.3, text("Fair agreement", :black, :left, 8))
plot!([0, 1, 1, 0], [0.4, 0.4, 0.6, 0.6], st=:shape, c=:grey, alpha=0.3, lw=0.0)
annotate!(0.05, 0.5, text("Moderate agreement", :black, :left, 8))
plot!([0, 1, 1, 0], [0.6, 0.6, 0.8, 0.8], st=:shape, c=:grey, alpha=0.2, lw=0.0)
annotate!(0.05, 0.7, text("Good agreement", :black, :left, 8))
annotate!(0.95, 0.9, text("Very good agreement", :black, :right, 8))
scatter!(drop_at./links(M), drop_kappa, lab="", leg=false, msw=0.0, ms=3, c=:black)
xaxis!("Interactions withheld", (0, 1))
yaxis!("Cohen's κ", (0,1))
savefig(joinpath(suppfig, "sensibility_kappa.png"))

scatter(drop_at./links(M), drop_thr, leg=false, msw=0.0, ms=3, c=:black, aspectratio=1)
hline!([0.206], ls=:dash, c=:grey)
xaxis!("Interactions withheld", (0, 1))
yaxis!("Best threshold", (0,1))
savefig(joinpath(suppfig, "sensibility_threshold.png"))

scatter(drop_at./links(M), drop_recov, leg=false, msw=0.0, ms=3, c=:black, aspectratio=1)
xaxis!("Interactions withheld", (0, 1))
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