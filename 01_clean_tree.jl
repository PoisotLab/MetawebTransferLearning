using Phylo
using GBIF
using CSV, DataFrames
using ProgressMeter
using Base.Threads

tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]
treenodes = [n.name for n in tree.nodes if !startswith(n.name, "Node ")]

tree_cleanup_components = [
    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for
    i in 1:nthreads()
]

p = Progress(length(treenodes))
@threads for i in 1:length(treenodes)
    cname = replace(treenodes[i], '_' => ' ')
    try
        tax = GBIF.taxon(cname; strict=false, class="Mammalia")
        push!(
            tree_cleanup_components[threadid()],
            (treenodes[i], tax.species[1], tax.species[2], cname == tax.species[1]),
        )
    catch
        continue
    end
    next!(p)
end

tree_cleanup = vcat(tree_cleanup_components...)
CSV.write(joinpath("artifacts", "upham_gbif_names.csv"), tree_cleanup)
