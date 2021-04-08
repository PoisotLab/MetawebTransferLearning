using Phylo
using EcologicalNetworks
using DelimitedFiles
using CSV
using DataFrames
using ProgressMeter

# Read the tree
tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]

# Read the names
namelist = DataFrame(CSV.File(joinpath("artifacts", "names_metaweb_tree_gbif.csv")))

# Prepare the metaweb
mwspecies = unique(namelist.name)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)

# Read the metaweb
speciescodes = readdlm(joinpath("data", "Spp_Id.txt"))[2:end,:]
speciesdict = Dict([speciescodes[i,1] => speciescodes[i,2] for i in 1:size(speciescodes,1)])
mwlines = readlines(joinpath("data", "Metaweb_adults.csv"));
mwhead = [speciesdict[sp] for sp in replace.(split(mwlines[1], ","), '"' => "")[2:end]]
for row in mwlines[2:end]
    splitrow = replace.(split(row, ","), '"' => "")
    from = speciesdict[splitrow[1]]
    # Real name?
    realname = namelist[isequal(from).(namelist.metaweb),:name]
    if length(realname) == 0
        continue
    else
        sp_from = only(realname)
        int = findall(isequal("1"), splitrow[2:end])
        if !isempty(int)
            to = mwhead[int]
            to_names = namelist[map(n -> n in to, namelist.metaweb),:name]
            for t in to_names
                M[sp_from,t] = true
            end
        end
    end
end

simplify!(M)

open(joinpath("artifacts", "europeanmetaweb.csv"), "w") do euio
    for int in sort(interactions(M), by=x->x.from)
        println(euio, "$(int.from),$(int.to)")
    end
end
