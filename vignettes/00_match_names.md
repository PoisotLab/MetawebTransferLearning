# Step 1 - name matching

The names in the parts of the datasets are using different versions of the
taxonomy, and so we will reconcile everything using GBIF. The process is
entirely automated, and required no human decision. We rely on a mix of strict
matching, and then search through synonyms.

The scripts are going to require a number of folders, so we will create them
here. Note that it is assume you will download the `data/` folder, and the
`Project.toml`. More details on Julia package management can be found on the
`Pkg.jl` website, [https://pkgdocs.julialang.org/v1/].

````julia
if !isdir("data")
    error("You need to download the data folder from https://github.com/PoisotLab/MetawebTransferLearning/tree/main/data")
end

if !isfile("Project.toml")
    error("You need to download the project file from https://github.com/PoisotLab/MetawebTransferLearning/blob/main/Project.toml")
end

isdir("artifacts") || mkdir("artifacts")
isdir("figures") || mkdir("figures")
````

## Dependencies

Nothing here is out of the extroardinary - we do some processing on multiple
threads, and use `ProgressMeter` to provide visual feedback on the time to
completion.

````julia
using Phylo
using GBIF
using DataFrames
using CSV: CSV
using ProgressMeter
using DelimitedFiles
using Base.Threads
````

## Cleaning the European metaweb

This step is extracting the species names from the species list in the
European metaweb - we will need to go through two steps: removing the
non-mammals by filtering on the species code, and then extracting the names
these species were given in the European metaweb.

````julia
metaweb_names = String.(readdlm(joinpath("data", "Spp_Id.txt")))
mammal_positions = findall(startswith("M"), metaweb_names[:, 1])
````

This line specifically is going to create the list of mammal names that are
recorded in the European metaweb. As such, this is the first one we will
clean.

````julia
mammal_names = metaweb_names[mammal_positions, 2]
````

To ensure that we can join dataframes together, we will create a data frame
with the species code, the name and ID of the matched species in the GBIF
backbone, and finally, a flag to check that the two names are equal. This flag
is useful for manual inspection; most often, names are unequal because the
European metaweb uses deprecated taxonomic names.

Note that we are actually creating *one* dataframe per thread. This avoids
several threads pushing to the same position at the dataframe, which always
results in an error. There are, of course, other ways to avoid the issue, but
this one is simple and takes a single line to reconcile.

````julia
gbif_cleanup_components = [
    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for
    i in 1:nthreads()
]
````

With the dataframe stores in place, we can distribute the name reconciliation
on the different threads. The core of this loop is a search on the GBIF web
API; this step is therefore going to be limited by the responses of GBIF. Note
that although we use unstrict search (to allow for synonyms, etc), we restrict
the search to *Mammalia*.

````julia
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
````

Data frames are now filled, and we need to bring them together in a single one
-- this is a one-line call to `vcat`, as we simply need to stack them
vertically.

````julia
gbif_cleanup = vcat(gbif_cleanup_components...)
````

## Clean the GBIF matches

Names in the European metaweb that are not matched are removed as a
precaution:

````julia
select!(gbif_cleanup, Not(:equal))
````

We finally rename the columns, to facilitate the joins later in the pipeline:

````julia
rename!(gbif_cleanup, :code => :metaweb)
rename!(gbif_cleanup, :gbifname => :name)
rename!(gbif_cleanup, :gbifid => :id)
````

We save the file in the `artifacts` folders, as a csv:

````julia
CSV.write(joinpath("artifacts", "metaweb_gbif.csv"), gbif_cleanup)
````

## Cleaning the phylogeny

The phylogeny suffers from much the same problem as the European metaweb. We
will therefore use much of the same solution to fix the names.

The first step is to read the tree as a nexus file using `Phylo`, and to get
rid of all the internal nodes -- we are only interested in the leaves of the
tree.

````julia
tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]
treenodes = [n.name for n in tree.nodes if !startswith(n.name, "Node ")]
````

This step is a little dirty. Remember that data cleaning is a sin eater, and
we need to eat other people's data sins at some point. This is this point.
Essentially, we will iterate through the rows of the cleaned metaweb names,
and match on either `Gen sp` or `Gen_sp`.

There is an edge case that was not solved automatically, and so the *S.
musicus* species gets its correct tree name manually. This name is caught
correctly by *e.g.* `NCBITaxonomy.jl`, but this introduced another dependency
in the project, and requires another re-harmonization with the GBIF backbone.

````julia
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
````

## Merging the cleaned names together

Because the positions in `n` match with the row indices in the cleaned metaweb
names, we add this object as a column:

````julia
gbif_cleanup.upham = n
````

We finally save this new file as another artifact. Generally, we will err on
the side of caution and save multiple files with redudancy in them. This is
good practice, and it allows for easier data dumpster diving if something goes
wrong.

````julia
CSV.write(joinpath("artifacts", "names_metaweb_tree_gbif.csv"), gbif_cleanup)
````

