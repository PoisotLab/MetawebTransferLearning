# # Step 4 - bringing the names together

# At this point, we are ready to bring the different cleaned named together, in
# order to move on to the actual prediction. By the end of this step, we will
# have a cleaned name of mammals, corresponding to the tree, the Canadian
# species pool, and the European metaweb, all reconciled to the GBIF taxonomy
# backbone.

using Phylo
using EcologicalNetworks
using DelimitedFiles
using CSV
using DataFrames
using ProgressMeter

# We start by re-reading the tree from its nexus file - note that in the next
# steps, we will move away from `Phylo` to use `PhyloNetworks`, which handles
# the actual simulation of characters. But for name manipulation, `Phylo` is
# slightly easier to work with.

tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]

# We can not read the file with the names of the European metaweb, their tree
# equivalent, and the corect GBIF name.

namelist = DataFrame(CSV.File(joinpath("artifacts", "names_metaweb_tree_gbif.csv")))

# ## Prepare the metaweb from European data

# We will create a `UnipartiteNetwork` to store the European interactions. This
# is a much sparser version of the way the data are originally presented, and we
# will furthermore ensure that the species have the correct (*i.e.* GBIF) names.
# `EcologicalNetworks.jl` can handle having taxa objects from GBIF as nodes, but
# this is not something we will do here; strings are more than enough to do the
# matching, and we will not get back to the GBIF functions past this step.

mwspecies = unique(namelist.name)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)

# ## Filling the European metaweb 

# We need to start by preparing a dictionary of values, linking species codes to
# species names (metaweb names, that is, not the correct GBIF names):

speciescodes = readdlm(joinpath("data", "Spp_Id.txt"))[2:end, :]
speciesdict = Dict([
    speciescodes[i, 1] => speciescodes[i, 2] for i in 1:size(speciescodes, 1)
])

# The next step is to read the adjacency matrix for the European metaweb, which
# has species codes as identifiers, and `0` or `1` as values:

mwlines = readlines(joinpath("data", "Metaweb_adults.csv"))

# This next line will read the elements of the first row, from columns 2 to the
# end, and replace the codes by the names - this is our sorted list of species
# we can use:

mwhead = [speciesdict[sp] for sp in replace.(split(mwlines[1], ","), '"' => "")[2:end]]

# We then walk through the rows one by one, splitting them on the separator
# (`,`), and using the first element to identify the species. Everything that
# has a `1` is added as an interaction to the network object:

for row in mwlines[2:end]
    splitrow = replace.(split(row, ","), '"' => "")
    from = speciesdict[splitrow[1]]
    realname = namelist[isequal(from).(namelist.metaweb), :name]
    if length(realname) == 0
        continue
    else
        sp_from = only(realname)
        int = findall(isequal("1"), splitrow[2:end])
        if !isempty(int)
            to = mwhead[int]
            to_names = namelist[map(n -> n in to, namelist.metaweb), :name]
            for t in to_names
                M[sp_from, t] = true
            end
        end
    end
end

# Finally, out of precaution, we drop the species without interactions (there
# are none), and drop the zeros from the sparse matrix in which interactions are
# stored:

simplify!(M)

# We then save this artifact as a much more readable CSV edgelist, which we will
# use for the rest of the analysis:

open(joinpath("artifacts", "europeanmetaweb.csv"), "w") do euio
    for int in sort(interactions(M); by=x -> x.from)
        println(euio, "$(int.from),$(int.to)")
    end
end
