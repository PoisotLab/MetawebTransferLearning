using HTTP
using ProgressMeter
using JSON
using DataFrames
using CSV: CSV
using StatsPlots

theme(:mute)
default(; frame=:box)

# API info
_globi_api = "https://api.globalbioticinteractions.org/taxon"
relevant_types = ["eats", "preysOn"]

canmet = DataFrame(CSV.File("artifacts/canadian_corrected.csv"))

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

intcode = canmet.from .* canmet.to
diet.intcode = diet.from .* diet.to

missedint = select(diet[findall([!(d in intcode) for d in diet.intcode]), :], [:from, :to])
sort!(missedint, [:from, :to])
CSV.write("artifacts/globi_newinteractions.csv")

matchedint = dropmissing(leftjoin(diet, canmet; on=[:from => :from, :to => :to]))

density(matchedint.score; dpi=600, size=(500, 500), lab="GLOBI", fill=(0.2, 0), lw=0.0)
density!(canmet.score; ls=:dot, lw=2, c=:grey, lab="All predictions")
xaxis!("Imputed probability", (0, 1))
yaxis!("Density", (0, 8))
savefig("figures/globi-comparison.png")

aug = leftjoin(missedint, canmet; on=[:from, :to])
aug.score .= 1.0

inflated = vcat(canmet, aug)

sort!(inflated, [:score, :from, :to]; rev=[true, false, false])

# Save the corrected network
CSV.write("artifacts/canadian_inflated.csv", inflated)
