# # Step 2 - tree cleaning

using Phylo
using GBIF
using CSV, DataFrames
using ProgressMeter
using Base.Threads

# While we have made a number of matches in the previous step, we want to make
# sure that the tree names are *entirely* reconciled to the GBIF version of the
# European metaweb names.

tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]
treenodes = [n.name for n in tree.nodes if !startswith(n.name, "Node ")]

# We use the same basic approach as for the metaweb name matching, *i.e.* a
# collection of data frames meant to make the name cleaning thread-safe.

tree_cleanup_components = [
    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for
    i in 1:nthreads()
]

# This code is, again, similar to the previous step - the only difference is
# that we need to get rid of the `_` that phylogeny files love so much.

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

# As previously, we create a new artifact (note that it is not merged to the
# reconciled metaweb names - we are mostly accumulating keys for joins at this
# point).

tree_cleanup = vcat(tree_cleanup_components...)
CSV.write(joinpath("artifacts", "upham_gbif_names.csv"), tree_cleanup)
