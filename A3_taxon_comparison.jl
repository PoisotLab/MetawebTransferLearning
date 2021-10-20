using GBIF
using DataFrames
import CSV
using DelimitedFiles

# Get the list of the European mammals
function mammal_europe()
    return sort(unique(vec(readdlm(joinpath("artifacts", "europeanmetaweb.csv"), ','))))
end

function mammal_canada()
    return sort(unique(vec(readdlm(joinpath("artifacts/canadian_thresholded.csv"), ',')[2:end,1:2])))
end

list_can = mammal_canada()
list_eur = mammal_europe()

master_list = Dict{String, GBIFTaxon}()
for tax in string.(list_can ∪ list_eur)
    master_list[tax] = GBIF.taxon(tax)
end

genus_can = unique([master_list[s].genus.first for s in list_can])
genus_eur = unique([master_list[s].genus.first for s in list_eur])

family_can = unique([master_list[s].family.first for s in list_can])
family_eur = unique([master_list[s].family.first for s in list_eur])

order_can = unique([master_list[s].order.first for s in list_can])
order_eur = unique([master_list[s].order.first for s in list_eur])

congenerics = 0
for sp in list_can
    has_congeneric = false
    gen = master_list[sp].genus.first
    for seur in list_eur
        if sp != seur
            if master_list[seur].genus.first == gen
                has_congeneric = true
                continue
            end
        end
    end
    if has_congeneric
        congenerics += 1
    end
end

results = DataFrame(level=String[], can=Int[], common=Int[], eur=Int[], overlap=Float64[])
push!(results,
    ("species", length(setdiff(list_can, list_eur)), length(list_can ∩ list_eur), length(setdiff(list_can, list_eur)), length(list_can ∩ list_eur)/length(list_can ∪ list_eur))
)

push!(results,
    ("genus", length(setdiff(genus_can, genus_eur)), length(genus_can ∩ genus_eur), length(setdiff(genus_can, genus_eur)), length(genus_can ∩ genus_eur)/length(genus_can ∪ genus_eur))
)

push!(results,
    ("family", length(setdiff(family_can, family_eur)), length(family_can ∩ family_eur), length(setdiff(family_can, family_eur)), length(family_can ∩ family_eur)/length(family_can ∪ family_eur))
)

push!(results,
    ("order", length(setdiff(order_can, order_eur)), length(order_can ∩ order_eur), length(setdiff(order_can, order_eur)), length(order_can ∩ order_eur)/length(order_can ∪ order_eur))
)

CSV.write("artifacts/taxonomic_overlap.csv", results)