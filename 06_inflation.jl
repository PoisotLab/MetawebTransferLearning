#%% Get the dependencies
using HTTP
using ProgressMeter
using JSON
using DataFrames
using CSV: CSV
using StatsPlots

theme(:mute)
default(; frame=:box)

#%% API info
_globi_api = "https://api.globalbioticinteractions.org/taxon"
relevant_types = ["eats", "preysOn"]

#%% Read the corrected Canadian metaweb
canmet = DataFrame(CSV.File("artifacts/canadian_corrected.csv"))

#%% Get the data from GLOBI
allsp = unique(vcat(canmet.from, canmet.to))

diet = DataFrame(; from=String[], to=String[])
@showprogress for sp in allsp
    url = "$(_globi_api)/$(sp)/eats/"
    r = HTTP.request("GET", url)
    globidiet = JSON.parse(String(r.body))
    if !isempty(globidiet["data"])
        for intlist in globidiet["data"]
            if intlist[2] in relevant_types
                for s in intlist[3]
                    if s in allsp
                        if s != sp
                            push!(diet, (sp, s))
                        end
                    end
                end
            end
        end
    end
end
diet = unique(diet)

sort!(diet, [:from, :to])
CSV.write("artifacts/globi_diet.csv", diet)

#%% Get the data we missed from GLOBI
intcode = canmet.from .* canmet.to
diet.intcode = diet.from .* diet.to

missedint = select(diet[findall([!(d in intcode) for d in diet.intcode]), :], [:from, :to])
sort!(missedint, [:from, :to])
CSV.write("artifacts/globi_newinteractions.csv", missedint)

matchedint = dropmissing(leftjoin(diet, canmet; on=[:from => :from, :to => :to]))

@info "GLOBI: found $(size(missedint, 1)) new interactions out of $(size(diet, 1))"
@info "GLOBI: $(length(unique(vcat(diet.from, diet.to)))) species"

#%% Get the data we missed from Strong & Leroux
canmetsp = unique(vcat(canmet.from, canmet.to))
sl = DataFrame(CSV.File("artifacts/newfoundland.csv"))
slshared = intersect(canmetsp, unique(vcat(sl.from, sl.to)))

slkeep = [(i.from in slshared)&(i.to in slshared) for i in eachrow(sl)]
sl = sl[findall(slkeep),:]

sl.intcode = sl.from .* sl.to

missedslint = select(sl[findall([!(d in intcode) for d in sl.intcode]), :], [:from, :to])

sort!(missedslint, [:from, :to])
CSV.write("artifacts/newfoundland_newinteractions.csv", missedslint)

matchedslint = dropmissing(leftjoin(sl, canmet; on=[:from => :from, :to => :to]))

density(matchedint.score; dpi=600, size=(500, 500), lab="GLOBI", lw=2.0)
density!(matchedslint.score; lab="Newfoundland", lw=2.0)
density!(canmet.score; lw=0.0, fill=(0.2, 0), c=:black, lab="All predictions")
xaxis!("Imputed probability", (0, 1))
yaxis!("Density", (0, 8))

savefig("figures/globi-nfl-comparison.png")

@info "NFLD : found $(size(missedslint, 1)) new interactions out of $(size(sl, 1))"
@info "NFLD : $(length(unique(vcat(sl.from, sl.to)))) species"

#%% Merge all interactions
aug = leftjoin(unique(vcat(missedint, missedslint)), canmet; on=[:from, :to])
aug.score .= 1.0

inflated = vcat(canmet, aug)

sort!(inflated, [:score, :from, :to]; rev=[true, false, false])

# Save the corrected network
CSV.write("artifacts/canadian_inflated.csv", inflated)
