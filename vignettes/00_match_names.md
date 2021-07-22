# Step 1 - name matching

The names in the parts of the datasets are using different versions of the
taxonomy, and so we will reconcile everything using GBIF.

````julia
using Phylo
using GBIF
using DataFrames
using CSV: CSV
using ProgressMeter
using DelimitedFiles
using Base.Threads

metaweb_names = String.(readdlm(joinpath("data", "Spp_Id.txt")))
mammal_positions = findall(startswith("M"), metaweb_names[:, 1])
mammal_names = metaweb_names[mammal_positions, 2]

gbif_cleanup_components = [
    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for
    i in 1:nthreads()
]

p = Progress(length(mammal_names))
@threads for i in 1:length(mammal_names)
    cname = replace(mammal_names[i], '_' => ' ')
    try
        tax = GBIF.taxon(cname; strict=false, class="Mammalia")
        push!(
            gbif_cleanup_components[threadid()],
            (mammal_names[i], tax.species[1], tax.species[2], cname == tax.species[1]),
        )
    catch
        continue
    end
    next!(p)
end

gbif_cleanup = vcat(gbif_cleanup_components...)

select!(gbif_cleanup, Not(:equal))
rename!(gbif_cleanup, :code => :metaweb)
rename!(gbif_cleanup, :gbifname => :name)
rename!(gbif_cleanup, :gbifid => :id)
CSV.write(joinpath("artifacts", "metaweb_gbif.csv"), gbif_cleanup)

tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]
treenodes = [n.name for n in tree.nodes if !startswith(n.name, "Node ")]

n = []
for row in eachrow(gbif_cleanup)
    if replace(row.name, ' ' => '_') in treenodes
        push!(n, replace(row.name, ' ' => '_'))
    elseif row.metaweb in treenodes
        push!(n, row.metaweb)
    elseif row.name == "Spermophilus musicus"
        push!(n, "Spermophilus_pygmaeus")
    else
    end
end

gbif_cleanup.upham = n
CSV.write(joinpath("artifacts", "names_metaweb_tree_gbif.csv"), gbif_cleanup)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

