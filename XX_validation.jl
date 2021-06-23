using DelimitedFiles
using Plots
using EcologicalNetworks
using EcologicalNetworksPlots
using DataFrames
import CSV
using GBIF

eurometa = readdlm(joinpath("artifacts", "europeanmetaweb.csv"), ',', String)
mwspecies = unique(eurometa)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(eurometa, 1)
    M[eurometa[i,:]...] = true
end

canmeta = readdlm(joinpath("artifacts", "canadianmetaweb.csv"), ',', String)
mwspecies = unique(canmeta)
N = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(canmeta, 1)
    N[canmeta[i,:]...] = true
end

eltontraits = DataFrame(CSV.File(joinpath("data", "EltonMammals.txt")))

masses = dropmissing(select(eltontraits, [:Scientific, Symbol("BodyMass-Value")]))
rename!(masses, :Scientific => :species)
rename!(masses, Symbol("BodyMass-Value") => :bodymass)
bm = masses[filter(!isnothing, indexin(species(M), masses.species)),:]
bmdict = Dict{String,Float64}([x.species => x.bodymass for x in eachrow(masses)])

pl = plot(; dpi=600)
for i in interactions(N)
    if haskey(bmdict, i.from)
        if haskey(bmdict, i.to)
            scatter!(pl, (bmdict[i.from], bmdict[i.to]), lab="", c=:lightgrey, ms=5)
        end
    end
end

bm = masses[filter(!isnothing, indexin(species(M), masses.species)),:]
bmdict = Dict{String,Float64}([x.species => x.bodymass for x in eachrow(masses)])

for i in interactions(M)
    if haskey(bmdict, i.from)
        if haskey(bmdict, i.to)
            scatter!(pl, (bmdict[i.from], bmdict[i.to]), lab="", c=:black, ms=3)
        end
    end
end


xaxis!(pl, "Predator bodymass", (1e2, 1e6), :log10)
yaxis!(pl, "Prey bodymass", (1e0, 1e7), :log10)
savefig(joinpath("figures", "bodymass.png"))
