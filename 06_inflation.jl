# # Step 7 - inflating the predictions

# In this step, we will add the data from Newfoundland, and get the data from
# the Global Biotic Interactions Database (GLOBI), in order to inflate the
# predictions made during the t-SVD/RDPG step.

using HTTP
using ProgressMeter
using JSON
using DataFrames
using CSV: CSV
using StatsPlots

# There will be some plots, with the same visual aspect as the main text ones.

theme(:mute)
default(; frame=:box)

# ## Writing a paper-thin GLOBI wrapper

# GLOBI uses their own types for interactions, which is at best very losely
# defined - to clarify, it follows an ontology, but this ontology happens to be
# largely disconnected from the concerns of ecologists, and the differences
# between `a preysOn b` and `b preysUpon a` are a great mystery. To solve this
# conundrum, we performed a series of manual data searches using some of the
# terms related to predation, and picked two that gave the least spurious
# results.

_globi_api = "https://api.globalbioticinteractions.org/taxon"
relevant_types = ["eats", "preysOn"]

# What we will do next is get the species from the Canadian metaweb (we do not
# really care for the European metaweb anymore, as it has been embedded and
# transfered), and see if interactions between any pairs of them exist in GLOBI.

canmet = DataFrame(CSV.File("artifacts/canadian_corrected.csv"))
allsp = unique(vcat(canmet.from, canmet.to))

# We will again aggregate everything into a data frame.

diet = DataFrame(; from=String[], to=String[])

# Our "paper-thin" GLOBI client is here: it constructs an API query URL based on
# the name of the taxa, and extract the "meat" from the JSON response.

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

# We now save the GLOBI diet as a CSV file

sort!(diet, [:from, :to])
CSV.write("artifacts/globi_diet.csv", diet)

# ## Inflation of the Canadian metaweb with the GLOBI data

# We now add all GLOBI interactions in the Canadian predicted metaweb, with a
# probability of 1.

intcode = canmet.from .* canmet.to
diet.intcode = diet.from .* diet.to

# We split these interactions for additional examinations as needed, and save as
# an artifact.

missedint = select(diet[findall([!(d in intcode) for d in diet.intcode]), :], [:from, :to])
sort!(missedint, [:from, :to])
CSV.write("artifacts/globi_newinteractions.csv", missedint)

# To get the proportion of validated interactions, we now do the opposite: keep
# the matched interactions:

matchedint = dropmissing(leftjoin(diet, canmet; on=[:from => :from, :to => :to]))

# And we print the proportions:

@info "GLOBI: found $(size(missedint, 1)) new interactions out of $(size(diet, 1))"
@info "GLOBI: $(length(unique(vcat(diet.from, diet.to)))) species"

# ## Inflation of the Canadian metaweb with the Newfoundland metaweb

# We now do the same things as above using the Newfoundland data.

canmetsp = unique(vcat(canmet.from, canmet.to))
sl = DataFrame(CSV.File("artifacts/newfoundland.csv"))
slshared = intersect(canmetsp, unique(vcat(sl.from, sl.to)))

slkeep = [(i.from in slshared) & (i.to in slshared) for i in eachrow(sl)]
sl = sl[findall(slkeep), :]

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

savefig("figures/inflation-comparison.png")

@info "NFLD : found $(size(missedslint, 1)) new interactions out of $(size(sl, 1))"
@info "NFLD : $(length(unique(vcat(sl.from, sl.to)))) species"

# ## Merging all the additional interaction sources

# We can finally join the different missed interactions from GLOBI and the
# Newfoundland datasets:

aug = leftjoin(unique(vcat(missedint, missedslint)), canmet; on=[:from, :to])
aug.score .= 1.0

# We simply add them to the bottom of the Canadian metaweb, to get the inflated
# version:

inflated = vcat(canmet, aug)

sort!(inflated, [:score, :from, :to]; rev=[true, false, false])

# We are now ready to save this metaweb, which is the penultimate data product
# of this pipeline:

CSV.write("artifacts/canadian_inflated.csv", inflated)
