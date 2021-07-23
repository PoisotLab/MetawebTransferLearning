# Phylogenetic inference of ecological interactions through network embedding

## Executive summary

1. network decomposition through a truncated-SVD capture the evolutionary signal of species interactions
2. the decomposition into left/right subspaces is unique, represents latent traits for resp. outgoing and incoming edges, and can be used to predict interactions
3. we can infer latent traits for unobserved species based on phylogenetic proximity with observed species
4. we use this information to transfer information from the trophic interaction of European mammals to Canadian mammals for which we have no *a priori* data

## About this README

This README is generated using `Literate.jl` and contains *every line* or *every
file* used to produce the entire results. This is, essentially, the director's
commentary version of the analysis - there are discussions of the purpose, but
also discussion of technical choices.

## Reproducing the code

| Step                           | Script | Annotated MD | Notebook |
| ------------------------------ | ------ | ------------ | -------- |
| Match names                    |        |              |          |
| Cleanup phylogeny              |
| Cleanup IUCN                   |
| Reconcile mammal names         |
| Make the predictions           |
| Compare with Newfoundland data |
| Compare with GLOBI data        |
| Produce the Canadian metaweb   |


{
 "cells": [
  {
   "cell_type": "markdown",
   "source": [
    "# Step 1 - name matching\n",
    "\n",
    "The names in the parts of the datasets are using different versions of the\n",
    "taxonomy, and so we will reconcile everything using GBIF. The process is\n",
    "entirely automated, and required no human decision. We rely on a mix of strict\n",
    "matching, and then search through synonyms."
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Dependencies\n",
    "\n",
    "Nothing here is out of the extroardinary - we do some processing on multiple\n",
    "threads, and use `ProgressMeter` to provide visual feedback on the time to\n",
    "completion."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "using Phylo\n",
    "using GBIF\n",
    "using DataFrames\n",
    "using CSV: CSV\n",
    "using ProgressMeter\n",
    "using DelimitedFiles\n",
    "using Base.Threads"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Cleaning the European metaweb"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "This step is extracting the species names from the species list in the\n",
    "European metaweb - we will need to go through two steps: removing the\n",
    "non-mammals by filtering on the species code, and then extracting the names\n",
    "these species were given in the European metaweb."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "metaweb_names = String.(readdlm(joinpath(\"data\", \"Spp_Id.txt\")))\n",
    "mammal_positions = findall(startswith(\"M\"), metaweb_names[:, 1])"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "This line specifically is going to create the list of mammal names that are\n",
    "recorded in the European metaweb. As such, this is the first one we will\n",
    "clean."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "mammal_names = metaweb_names[mammal_positions, 2]"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "To ensure that we can join dataframes together, we will create a data frame\n",
    "with the species code, the name and ID of the matched species in the GBIF\n",
    "backbone, and finally, a flag to check that the two names are equal. This flag\n",
    "is useful for manual inspection; most often, names are unequal because the\n",
    "European metaweb uses deprecated taxonomic names."
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "Note that we are actually creating *one* dataframe per thread. This avoids\n",
    "several threads pushing to the same position at the dataframe, which always\n",
    "results in an error. There are, of course, other ways to avoid the issue, but\n",
    "this one is simple and takes a single line to reconcile."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "gbif_cleanup_components = [\n",
    "    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for\n",
    "    i in 1:nthreads()\n",
    "]"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "With the dataframe stores in place, we can distribute the name reconciliation\n",
    "on the different threads. The core of this loop is a search on the GBIF web\n",
    "API; this step is therefore going to be limited by the responses of GBIF. Note\n",
    "that although we use unstrict search (to allow for synonyms, etc), we restrict\n",
    "the search to *Mammalia*."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "p = Progress(length(mammal_names))\n",
    "@threads for i in 1:length(mammal_names)\n",
    "    cname = replace(mammal_names[i], '_' => ' ')\n",
    "    try\n",
    "        tax = GBIF.taxon(cname; strict=false, class=\"Mammalia\")\n",
    "        push!(\n",
    "            gbif_cleanup_components[threadid()],\n",
    "            (mammal_names[i], tax.species[1], tax.species[2], cname == tax.species[1]),\n",
    "        )\n",
    "    catch\n",
    "        continue\n",
    "    end\n",
    "    next!(p)\n",
    "end"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "Data frames are now filled, and we need to bring them together in a single one\n",
    "-- this is a one-line call to `vcat`, as we simply need to stack them\n",
    "vertically."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "gbif_cleanup = vcat(gbif_cleanup_components...)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Clean the GBIF matches"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "Names in the European metaweb that are not matched are removed as a\n",
    "precaution:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "select!(gbif_cleanup, Not(:equal))"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We finally rename the columns, to facilitate the joins later in the pipeline:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "rename!(gbif_cleanup, :code => :metaweb)\n",
    "rename!(gbif_cleanup, :gbifname => :name)\n",
    "rename!(gbif_cleanup, :gbifid => :id)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We save the file in the `artifacts` folders, as a csv:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "CSV.write(joinpath(\"artifacts\", \"metaweb_gbif.csv\"), gbif_cleanup)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Cleaning the phylogeny"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "The phylogeny suffers from much the same problem as the European metaweb. We\n",
    "will therefore use much of the same solution to fix the names."
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "The first step is to read the tree as a nexus file using `Phylo`, and to get\n",
    "rid of all the internal nodes -- we are only interested in the leaves of the\n",
    "tree."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "tree = open(parsenexus, joinpath(\"data\", \"mammals.nex\"))[\"*UNTITLED\"]\n",
    "treenodes = [n.name for n in tree.nodes if !startswith(n.name, \"Node \")]"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "This step is a little dirty. Remember that data cleaning is a sin eater, and\n",
    "we need to eat other people's data sins at some point. This is this point.\n",
    "Essentially, we will iterate through the rows of the cleaned metaweb names,\n",
    "and match on either `Gen sp` or `Gen_sp`."
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "There is an edge case that was not solved automatically, and so the *S.\n",
    "musicus* species gets its correct tree name manually. This name is caught\n",
    "correctly by *e.g.* `NCBITaxonomy.jl`, but this introduced another dependency\n",
    "in the project, and requires another re-harmonization with the GBIF backbone."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "n = []\n",
    "for row in eachrow(gbif_cleanup)\n",
    "    if replace(row.name, ' ' => '_') in treenodes\n",
    "        push!(n, replace(row.name, ' ' => '_'))\n",
    "    elseif row.metaweb in treenodes\n",
    "        push!(n, row.metaweb)\n",
    "    elseif row.name == \"Spermophilus musicus\"\n",
    "        push!(n, \"Spermophilus_pygmaeus\")\n",
    "    else\n",
    "    end\n",
    "end"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Merging the cleaned names together"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "Because the positions in `n` match with the row indices in the cleaned metaweb\n",
    "names, we add this object as a column:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "gbif_cleanup.upham = n"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We finally save this new file as another artifact. Generally, we will err on\n",
    "the side of caution and save multiple files with redudancy in them. This is\n",
    "good practice, and it allows for easier data dumpster diving if something goes\n",
    "wrong."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "CSV.write(joinpath(\"artifacts\", \"names_metaweb_tree_gbif.csv\"), gbif_cleanup)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "---\n",
    "\n",
    "*This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*"
   ],
   "metadata": {}
  }
 ],
 "nbformat_minor": 3,
 "metadata": {
  "language_info": {
   "file_extension": ".jl",
   "mimetype": "application/julia",
   "name": "julia",
   "version": "1.6.2"
  },
  "kernelspec": {
   "name": "julia-1.6",
   "display_name": "Julia 1.6.2",
   "language": "julia"
  }
 },
 "nbformat": 4
}


# Step 1 - name matching

The names in the parts of the datasets are using different versions of the
taxonomy, and so we will reconcile everything using GBIF. The process is
entirely automated, and required no human decision. We rely on a mix of strict
matching, and then search through synonyms.

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



{
 "cells": [
  {
   "cell_type": "markdown",
   "source": [
    "# Step 2 - tree cleaning"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "using Phylo\n",
    "using GBIF\n",
    "using CSV, DataFrames\n",
    "using ProgressMeter\n",
    "using Base.Threads"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "While we have made a number of matches in the previous step, we want to make\n",
    "sure that the tree names are *entirely* reconciled to the GBIF version of the\n",
    "European metaweb names."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "tree = open(parsenexus, joinpath(\"data\", \"mammals.nex\"))[\"*UNTITLED\"]\n",
    "treenodes = [n.name for n in tree.nodes if !startswith(n.name, \"Node \")]"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We use the same basic approach as for the metaweb name matching, *i.e.* a\n",
    "collection of data frames meant to make the name cleaning thread-safe."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "tree_cleanup_components = [\n",
    "    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for\n",
    "    i in 1:nthreads()\n",
    "]"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "This code is, again, similar to the previous step - the only difference is\n",
    "that we need to get rid of the `_` that phylogeny files love so much."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "p = Progress(length(treenodes))\n",
    "@threads for i in 1:length(treenodes)\n",
    "    cname = replace(treenodes[i], '_' => ' ')\n",
    "    try\n",
    "        tax = GBIF.taxon(cname; strict=false, class=\"Mammalia\")\n",
    "        push!(\n",
    "            tree_cleanup_components[threadid()],\n",
    "            (treenodes[i], tax.species[1], tax.species[2], cname == tax.species[1]),\n",
    "        )\n",
    "    catch\n",
    "        continue\n",
    "    end\n",
    "    next!(p)\n",
    "end"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "As previously, we create a new artifact (note that it is not merged to the\n",
    "reconciled metaweb names - we are mostly accumulating keys for joins at this\n",
    "point)."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "tree_cleanup = vcat(tree_cleanup_components...)\n",
    "CSV.write(joinpath(\"artifacts\", \"upham_gbif_names.csv\"), tree_cleanup)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "---\n",
    "\n",
    "*This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*"
   ],
   "metadata": {}
  }
 ],
 "nbformat_minor": 3,
 "metadata": {
  "language_info": {
   "file_extension": ".jl",
   "mimetype": "application/julia",
   "name": "julia",
   "version": "1.6.2"
  },
  "kernelspec": {
   "name": "julia-1.6",
   "display_name": "Julia 1.6.2",
   "language": "julia"
  }
 },
 "nbformat": 4
}


# Step 2 - tree cleaning

````julia
using Phylo
using GBIF
using CSV, DataFrames
using ProgressMeter
using Base.Threads
````

While we have made a number of matches in the previous step, we want to make
sure that the tree names are *entirely* reconciled to the GBIF version of the
European metaweb names.

````julia
tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]
treenodes = [n.name for n in tree.nodes if !startswith(n.name, "Node ")]
````

We use the same basic approach as for the metaweb name matching, *i.e.* a
collection of data frames meant to make the name cleaning thread-safe.

````julia
tree_cleanup_components = [
    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for
    i in 1:nthreads()
]
````

This code is, again, similar to the previous step - the only difference is
that we need to get rid of the `_` that phylogeny files love so much.

````julia
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
````

As previously, we create a new artifact (note that it is not merged to the
reconciled metaweb names - we are mostly accumulating keys for joins at this
point).

````julia
tree_cleanup = vcat(tree_cleanup_components...)
CSV.write(joinpath("artifacts", "upham_gbif_names.csv"), tree_cleanup)
````



{
 "cells": [
  {
   "cell_type": "markdown",
   "source": [
    "# Step 3 - IUCN cleanup"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "using GBIF\n",
    "using CSV, DataFrames\n",
    "using ProgressMeter\n",
    "using Base.Threads"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We downloaded a checklist of mammals reported to be in Canada from the IUCN\n",
    "database. Before deciding on this solution, we examined a few alternatives,\n",
    "notably the use of GBIF occurrences. GBIF occurrences had a few issues,\n",
    "including spurious records, museum specimens incorrectly tagged, captive\n",
    "exotic species being reported as occurrences, etc."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "checklist = DataFrame(CSV.File(joinpath(\"data\", \"taxonomy.csv\")))"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Taxonomy filtering"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "The European metaweb is limited to \"terrestrial\" mammals. For this reason, we\n",
    "identified a number of taxonomic groups (mostly families) that are present in\n",
    "Canada but were excluded from the source dataset, and remove them."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "valid_rows = map(\n",
    "    fam ->\n",
    "        !(\n",
    "            fam ∈ [\n",
    "                \"BALAENIDAE\",\n",
    "                \"PHYSETERIDAE\",\n",
    "                \"DELPHINIDAE\",\n",
    "                \"BALAENOPTERIDAE\",\n",
    "                \"OTARIIDAE\",\n",
    "                \"PHOCIDAE\",\n",
    "                \"ODOBENIDAE\",\n",
    "                \"ZIPHIIDAE\",\n",
    "                \"MONODONTIDAE\",\n",
    "                \"ESCHRICHTIIDAE\",\n",
    "                \"KOGIIDAE\",\n",
    "                \"PHOCOENIDAE\"\n",
    "            ]\n",
    "        ),\n",
    "    checklist.familyName,\n",
    ")\n",
    "checklist = checklist[findall(valid_rows), :]"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Species-specific removal"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "*Neovison macrodon* (considered to be extinct) and *Enhydra lutris* (considered\n",
    "a marine mammal) are removed as well."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "extinct_sp = map(\n",
    "    sp ->\n",
    "        !(\n",
    "            sp ∈ [\n",
    "                \"Neovison macrodon\",\n",
    "                \"Enhydra lutris\"\n",
    "            ]\n",
    "        ),\n",
    "    checklist.scientificName,\n",
    ")\n",
    "checklist = checklist[findall(extinct_sp), :]"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Reconciliation on the GBIF names"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "By this point, the approach should be familiar: we will create a thread-safe\n",
    "structure for the name cleaning, and use the GBIF API to find the correct\n",
    "matches."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "checklist_cleanup_components = [\n",
    "    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for\n",
    "    i in 1:nthreads()\n",
    "]"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "Again, we get rid of `_` before doing the matching. This is actually *not*\n",
    "something we want built into the name cleaning function itself, because some\n",
    "taxa have underscores as valid identifiers. None of the taxa from this\n",
    "specific dataset do, but it is better to keep the low-level tools general, and\n",
    "make the specific changes in user-code."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "p = Progress(length(checklist.scientificName))\n",
    "@threads for i in 1:length(checklist.scientificName)\n",
    "    cname = replace(checklist.scientificName[i], '_' => ' ')\n",
    "    try\n",
    "        tax = GBIF.taxon(cname; strict=false, class=\"Mammalia\")\n",
    "        push!(\n",
    "            checklist_cleanup_components[threadid()],\n",
    "            (\n",
    "                checklist.scientificName[i],\n",
    "                tax.species[1],\n",
    "                tax.species[2],\n",
    "                cname == tax.species[1],\n",
    "            ),\n",
    "        )\n",
    "    catch\n",
    "        continue\n",
    "    end\n",
    "    next!(p)\n",
    "end"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We finally write the artifact:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "checklist_cleanup = vcat(checklist_cleanup_components...)\n",
    "CSV.write(joinpath(\"artifacts\", \"iucn_gbif_names.csv\"), checklist_cleanup)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "---\n",
    "\n",
    "*This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*"
   ],
   "metadata": {}
  }
 ],
 "nbformat_minor": 3,
 "metadata": {
  "language_info": {
   "file_extension": ".jl",
   "mimetype": "application/julia",
   "name": "julia",
   "version": "1.6.2"
  },
  "kernelspec": {
   "name": "julia-1.6",
   "display_name": "Julia 1.6.2",
   "language": "julia"
  }
 },
 "nbformat": 4
}


# Step 3 - IUCN cleanup

````julia
using GBIF
using CSV, DataFrames
using ProgressMeter
using Base.Threads
````

We downloaded a checklist of mammals reported to be in Canada from the IUCN
database. Before deciding on this solution, we examined a few alternatives,
notably the use of GBIF occurrences. GBIF occurrences had a few issues,
including spurious records, museum specimens incorrectly tagged, captive
exotic species being reported as occurrences, etc.

````julia
checklist = DataFrame(CSV.File(joinpath("data", "taxonomy.csv")))
````

## Taxonomy filtering

The European metaweb is limited to "terrestrial" mammals. For this reason, we
identified a number of taxonomic groups (mostly families) that are present in
Canada but were excluded from the source dataset, and remove them.

````julia
valid_rows = map(
    fam ->
        !(
            fam ∈ [
                "BALAENIDAE",
                "PHYSETERIDAE",
                "DELPHINIDAE",
                "BALAENOPTERIDAE",
                "OTARIIDAE",
                "PHOCIDAE",
                "ODOBENIDAE",
                "ZIPHIIDAE",
                "MONODONTIDAE",
                "ESCHRICHTIIDAE",
                "KOGIIDAE",
                "PHOCOENIDAE"
            ]
        ),
    checklist.familyName,
)
checklist = checklist[findall(valid_rows), :]
````

## Species-specific removal

*Neovison macrodon* (considered to be extinct) and *Enhydra lutris* (considered
a marine mammal) are removed as well.

````julia
extinct_sp = map(
    sp ->
        !(
            sp ∈ [
                "Neovison macrodon",
                "Enhydra lutris"
            ]
        ),
    checklist.scientificName,
)
checklist = checklist[findall(extinct_sp), :]
````

## Reconciliation on the GBIF names

By this point, the approach should be familiar: we will create a thread-safe
structure for the name cleaning, and use the GBIF API to find the correct
matches.

````julia
checklist_cleanup_components = [
    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for
    i in 1:nthreads()
]
````

Again, we get rid of `_` before doing the matching. This is actually *not*
something we want built into the name cleaning function itself, because some
taxa have underscores as valid identifiers. None of the taxa from this
specific dataset do, but it is better to keep the low-level tools general, and
make the specific changes in user-code.

````julia
p = Progress(length(checklist.scientificName))
@threads for i in 1:length(checklist.scientificName)
    cname = replace(checklist.scientificName[i], '_' => ' ')
    try
        tax = GBIF.taxon(cname; strict=false, class="Mammalia")
        push!(
            checklist_cleanup_components[threadid()],
            (
                checklist.scientificName[i],
                tax.species[1],
                tax.species[2],
                cname == tax.species[1],
            ),
        )
    catch
        continue
    end
    next!(p)
end
````

We finally write the artifact:

````julia
checklist_cleanup = vcat(checklist_cleanup_components...)
CSV.write(joinpath("artifacts", "iucn_gbif_names.csv"), checklist_cleanup)
````



{
 "cells": [
  {
   "cell_type": "markdown",
   "source": [
    "# Step 4 - bringing the names together"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "At this point, we are ready to bring the different cleaned named together, in\n",
    "order to move on to the actual prediction. By the end of this step, we will\n",
    "have a cleaned name of mammals, corresponding to the tree, the Canadian\n",
    "species pool, and the European metaweb, all reconciled to the GBIF taxonomy\n",
    "backbone."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "using Phylo\n",
    "using EcologicalNetworks\n",
    "using DelimitedFiles\n",
    "using CSV\n",
    "using DataFrames\n",
    "using ProgressMeter"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We start by re-reading the tree from its nexus file - note that in the next\n",
    "steps, we will move away from `Phylo` to use `PhyloNetworks`, which handles\n",
    "the actual simulation of characters. But for name manipulation, `Phylo` is\n",
    "slightly easier to work with."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "tree = open(parsenexus, joinpath(\"data\", \"mammals.nex\"))[\"*UNTITLED\"]"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We can not read the file with the names of the European metaweb, their tree\n",
    "equivalent, and the corect GBIF name."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "namelist = DataFrame(CSV.File(joinpath(\"artifacts\", \"names_metaweb_tree_gbif.csv\")))"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Prepare the metaweb from European data"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "We will create a `UnipartiteNetwork` to store the European interactions. This\n",
    "is a much sparser version of the way the data are originally presented, and we\n",
    "will furthermore ensure that the species have the correct (*i.e.* GBIF) names.\n",
    "`EcologicalNetworks.jl` can handle having taxa objects from GBIF as nodes, but\n",
    "this is not something we will do here; strings are more than enough to do the\n",
    "matching, and we will not get back to the GBIF functions past this step."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "mwspecies = unique(namelist.name)\n",
    "M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Filling the European metaweb"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "We need to start by preparing a dictionary of values, linking species codes to\n",
    "species names (metaweb names, that is, not the correct GBIF names):"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "speciescodes = readdlm(joinpath(\"data\", \"Spp_Id.txt\"))[2:end, :]\n",
    "speciesdict = Dict([\n",
    "    speciescodes[i, 1] => speciescodes[i, 2] for i in 1:size(speciescodes, 1)\n",
    "])"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "The next step is to read the adjacency matrix for the European metaweb, which\n",
    "has species codes as identifiers, and `0` or `1` as values:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "mwlines = readlines(joinpath(\"data\", \"Metaweb_adults.csv\"))"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "This next line will read the elements of the first row, from columns 2 to the\n",
    "end, and replace the codes by the names - this is our sorted list of species\n",
    "we can use:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "mwhead = [speciesdict[sp] for sp in replace.(split(mwlines[1], \",\"), '\"' => \"\")[2:end]]"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We then walk through the rows one by one, splitting them on the separator\n",
    "(`,`), and using the first element to identify the species. Everything that\n",
    "has a `1` is added as an interaction to the network object:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "for row in mwlines[2:end]\n",
    "    splitrow = replace.(split(row, \",\"), '\"' => \"\")\n",
    "    from = speciesdict[splitrow[1]]\n",
    "    realname = namelist[isequal(from).(namelist.metaweb), :name]\n",
    "    if length(realname) == 0\n",
    "        continue\n",
    "    else\n",
    "        sp_from = only(realname)\n",
    "        int = findall(isequal(\"1\"), splitrow[2:end])\n",
    "        if !isempty(int)\n",
    "            to = mwhead[int]\n",
    "            to_names = namelist[map(n -> n in to, namelist.metaweb), :name]\n",
    "            for t in to_names\n",
    "                M[sp_from, t] = true\n",
    "            end\n",
    "        end\n",
    "    end\n",
    "end"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "Finally, out of precaution, we drop the species without interactions (there\n",
    "are none), and drop the zeros from the sparse matrix in which interactions are\n",
    "stored:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "simplify!(M)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We then save this artifact as a much more readable CSV edgelist, which we will\n",
    "use for the rest of the analysis:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "open(joinpath(\"artifacts\", \"europeanmetaweb.csv\"), \"w\") do euio\n",
    "    for int in sort(interactions(M); by=x -> x.from)\n",
    "        println(euio, \"$(int.from),$(int.to)\")\n",
    "    end\n",
    "end"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "---\n",
    "\n",
    "*This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*"
   ],
   "metadata": {}
  }
 ],
 "nbformat_minor": 3,
 "metadata": {
  "language_info": {
   "file_extension": ".jl",
   "mimetype": "application/julia",
   "name": "julia",
   "version": "1.6.2"
  },
  "kernelspec": {
   "name": "julia-1.6",
   "display_name": "Julia 1.6.2",
   "language": "julia"
  }
 },
 "nbformat": 4
}


# Step 4 - bringing the names together

At this point, we are ready to bring the different cleaned named together, in
order to move on to the actual prediction. By the end of this step, we will
have a cleaned name of mammals, corresponding to the tree, the Canadian
species pool, and the European metaweb, all reconciled to the GBIF taxonomy
backbone.

````julia
using Phylo
using EcologicalNetworks
using DelimitedFiles
using CSV
using DataFrames
using ProgressMeter
````

We start by re-reading the tree from its nexus file - note that in the next
steps, we will move away from `Phylo` to use `PhyloNetworks`, which handles
the actual simulation of characters. But for name manipulation, `Phylo` is
slightly easier to work with.

````julia
tree = open(parsenexus, joinpath("data", "mammals.nex"))["*UNTITLED"]
````

We can not read the file with the names of the European metaweb, their tree
equivalent, and the corect GBIF name.

````julia
namelist = DataFrame(CSV.File(joinpath("artifacts", "names_metaweb_tree_gbif.csv")))
````

## Prepare the metaweb from European data

We will create a `UnipartiteNetwork` to store the European interactions. This
is a much sparser version of the way the data are originally presented, and we
will furthermore ensure that the species have the correct (*i.e.* GBIF) names.
`EcologicalNetworks.jl` can handle having taxa objects from GBIF as nodes, but
this is not something we will do here; strings are more than enough to do the
matching, and we will not get back to the GBIF functions past this step.

````julia
mwspecies = unique(namelist.name)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
````

## Filling the European metaweb

We need to start by preparing a dictionary of values, linking species codes to
species names (metaweb names, that is, not the correct GBIF names):

````julia
speciescodes = readdlm(joinpath("data", "Spp_Id.txt"))[2:end, :]
speciesdict = Dict([
    speciescodes[i, 1] => speciescodes[i, 2] for i in 1:size(speciescodes, 1)
])
````

The next step is to read the adjacency matrix for the European metaweb, which
has species codes as identifiers, and `0` or `1` as values:

````julia
mwlines = readlines(joinpath("data", "Metaweb_adults.csv"))
````

This next line will read the elements of the first row, from columns 2 to the
end, and replace the codes by the names - this is our sorted list of species
we can use:

````julia
mwhead = [speciesdict[sp] for sp in replace.(split(mwlines[1], ","), '"' => "")[2:end]]
````

We then walk through the rows one by one, splitting them on the separator
(`,`), and using the first element to identify the species. Everything that
has a `1` is added as an interaction to the network object:

````julia
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
````

Finally, out of precaution, we drop the species without interactions (there
are none), and drop the zeros from the sparse matrix in which interactions are
stored:

````julia
simplify!(M)
````

We then save this artifact as a much more readable CSV edgelist, which we will
use for the rest of the analysis:

````julia
open(joinpath("artifacts", "europeanmetaweb.csv"), "w") do euio
    for int in sort(interactions(M); by=x -> x.from)
        println(euio, "$(int.from),$(int.to)")
    end
end
````



{
 "cells": [
  {
   "cell_type": "markdown",
   "source": [
    "# Step 5 - the actual prediction"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "This is the largest step in the entire pipeline, but not necessarily a complex\n",
    "one. In short, this will generate the entire paper, minus some of the data\n",
    "inflation and post-processing steps. As such, there are a few dependencies in\n",
    "play."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "using PhyloNetworks\n",
    "using ProgressMeter\n",
    "using LinearAlgebra\n",
    "using EcologicalNetworks\n",
    "using DelimitedFiles\n",
    "using Random\n",
    "using StatsPlots\n",
    "using StatsBase\n",
    "using Statistics\n",
    "using SparseArrays\n",
    "using CSV\n",
    "using DataFrames\n",
    "using GBIF\n",
    "using Base.Threads\n",
    "using StatsBase\n",
    "using StatsModels\n",
    "using Distributions"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "These lines are ensuring that the figures all look uniform, and we also set a\n",
    "seed for reproducibility."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "theme(:mute)\n",
    "default(; frame=:box)\n",
    "Random.seed!(01189998819991197253)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "Some functions that are required for phylogenetic imputation are stored in\n",
    "another file, and we will load them here to be done with that."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "include(\"lib/pwar.jl\")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Reading the data pieces"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "Reading the tree is done with `PhyloNetworks` this time -- this is where we\n",
    "will get the Brownian motion code from:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "tree_net = readTopology(joinpath(\"data\", \"mammals.newick\"));"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We will grab back the corrected Europe/tree/GBIF names, as they will be used\n",
    "quite a lot to match simulations to species names:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "namelist = DataFrame(CSV.File(joinpath(\"artifacts\", \"names_metaweb_tree_gbif.csv\")))"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "Finally, we grab the European metaweb *not* from the original file, but from\n",
    "our edgelist artifact."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "eurometa = readdlm(joinpath(\"artifacts\", \"europeanmetaweb.csv\"), ',', String)\n",
    "mwspecies = unique(eurometa)\n",
    "M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)\n",
    "for i in 1:size(eurometa, 1)\n",
    "    M[eurometa[i, :]...] = true\n",
    "end"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Finding the rank for truncation"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "The ideal way to find a cutoff for the rank at which the matrix should be cut\n",
    "would involve the profile likelihood, or approximating the maximum curvature\n",
    "point on the screeplot using the central difference approximation for the\n",
    "second order partial derivative. Sadly all of these approaches will conclude\n",
    "that only the first dimension is required (there is a biological reason for\n",
    "this, to which we will return when we dig into the biological meaning of the\n",
    "latent variables)."
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "The only eigenvalues we care about are the ones that are up to the rank of the\n",
    "adjacency matrix - we will extract them, and range them so that they sum to\n",
    "unity. It makes no difference in the results, and allows us to read the\n",
    "proportion of variance explained directly."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "rnk = rank(adjacency(M))\n",
    "eig = svd(M).S[1:rnk]\n",
    "neig = eig ./ sum(eig)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "The first piece of information is a screeplot of the eigenvalues. As a rule,\n",
    "the document rendered from the literate files will not include the figures, as\n",
    "we do not want to go into the analysis (there is, after all, a whole paper for\n",
    "that)."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "scatter(eig; lab=\"\", dpi=600, size=(500, 500))\n",
    "xaxis!(\"Dimension\", (1, 38))\n",
    "yaxis!(\"Eigenvalue\", (0, 40))\n",
    "savefig(\"figures/screeplot.png\")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "The second diagnosis plot is the proportion of variance explained. We added\n",
    "lines marking the points that explain between 50% and 90% of the total\n",
    "variance, every 10%."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "scatter(cumsum(neig); lab=\"\", dpi=600, size=(500, 500))\n",
    "for ve in [0.5, 0.6, 0.7, 0.8, 0.9]\n",
    "    i = findfirst(cumsum(neig) .> ve)\n",
    "    x = cumsum(neig)[i]\n",
    "    plot!([0, i, NaN, i, i], [x, x, NaN, x, 0]; lab=\"\", c=:grey, ls=:dash)\n",
    "end\n",
    "xaxis!(\"Dimension\", (1, 38))\n",
    "yaxis!(\"Cumulative variance explained\", (0, 1.0))\n",
    "savefig(\"figures/varexplained.png\")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "At this point, we need to make an arbitrary call, which is to say to decide on\n",
    "a proportion of variance explained. We decided on 60%, which means that the\n",
    "first 12 ranks will be used."
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Embedding the European metaweb"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "L, R = rdpg(M, 12)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "The `L` and `R` variables are the left and right subspaces from the random dot\n",
    "product graph (which in turn are the left and right subspace from the t-SVD\n",
    "multiplied by the square root of the diagonal matrix containing the\n",
    "eigenvalues). Multiplying `L` and `R` (using a matrix multiplication!) will\n",
    "give a rank-12 approximation of the European metaweb."
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Thresholding the embedded network"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "The reconstucted network (`L*R`) will not have Boolean values. Essentially,\n",
    "the multiplication will approximate something with values in {0,1} by\n",
    "something with values in ℝ, and so we need to find a cutoff to separate\n",
    "interactions from non-interactions. This threshold is important because it\n",
    "holds for any reconstruction made using these latent variables - that is to\n",
    "say, when we infer the Canadian interactions using reconstrusted left and\n",
    "right subspaces, we will apply the same threshold."
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "Our basic approach to threshdolding is to examine 500 cutoffs points, get the\n",
    "confusion matrix, and retain as a threshold the value that maximizes\n",
    "Informedness:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "A = adjacency(M)\n",
    "thresholds = LinRange(extrema(L * R)..., 500)\n",
    "tp = zeros(Float64, length(thresholds))\n",
    "tn = similar(tp)\n",
    "fp = similar(tp)\n",
    "fn = similar(tp)\n",
    "for (i, t) in enumerate(thresholds)\n",
    "    PN = (L * R) .>= t\n",
    "    tp[i] = sum((A) .& (PN)) / sum(A)\n",
    "    tn[i] = sum((.!(A)) .& (.!(PN))) / sum(.!(A))\n",
    "    fp[i] = sum((.!(A)) .& (PN)) / sum(.!(A))\n",
    "    fn[i] = sum((A) .& (.!(PN))) / sum(A)\n",
    "end"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We can calculate Youden's J at every cutoff, and the threshold is going to be\n",
    "the threshold associated to the maximum Y:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "Y = tp .+ tn .- 1\n",
    "maxY, posY = findmax(Y)\n",
    "threshold = thresholds[posY]"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We can plot the optimal cutoff - this is not presented in the manuscript, but\n",
    "useful in case you need to satisfy your curiosity. The best cutoff gives an\n",
    "informedness essentially = 1, which is expected."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "plot(thresholds, Y; lw=0.0, fill=(0, 0.2), lab=\"\", dpi=600, size=(500, 500))\n",
    "scatter!([threshold], [maxY]; lab=\"\", c=:darkgrey, msw=0.0)\n",
    "yaxis!(\"Youden's J\", (0, 1))\n",
    "xaxis!(\"Cutoff\", extrema(L * R))\n",
    "savefig(\"figures/optimalcutoff.png\")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Visual examination of the subspaces"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "The left and right subspaces *do* hold ecological information, and so it is a\n",
    "good idea to check them out visually. One striking result is that species with\n",
    "a value of 0 in the left subspace also have no preys: this is a strong clue\n",
    "that the left subspace is associated to generality (in the sense of Schoener\n",
    "1989), a fact we will epxloit later on."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "plot(\n",
    "    heatmap(L; c=:GnBu, frame=:none, cbar=false),\n",
    "    heatmap(R'; c=:BuPu, frame=:none, cbar=false);\n",
    "    dpi=600,\n",
    "    size=(500, 400),\n",
    ")\n",
    "savefig(\"figures/subspaces.png\")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Preparing the tree for inference"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "We are now ready to move on to the next step: infering the values of the left\n",
    "and right subspaces for the species that are not in the European metaweb, but\n",
    "are in the tree. To do so, we will need the tip names:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "treeleaves = tipLabels(tree_net)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We are now ready to match the tree to the Canadian species pool names:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "canada = DataFrame(CSV.File(joinpath(\"artifacts\", \"iucn_gbif_names.csv\")))\n",
    "cancodes = replace.(unique(filter(!ismissing, canada.gbifname)), \" \" => \"_\")\n",
    "tree_cleanup = DataFrame(CSV.File(joinpath(\"artifacts\", \"upham_gbif_names.csv\")))\n",
    "\n",
    "csp = dropmissing!(DataFrame(; gbifname=canada.gbifname))\n",
    "csp = dropmissing(leftjoin(csp, tree_cleanup; on=:gbifname))"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We can store the cleaned Canadian species name into a data frame - this is\n",
    "mandated by the `PhyloNetworks` interface:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "canmammals = unique(csp.code)\n",
    "pool = DataFrame(; tipNames=canmammals)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We finally convert the network names to their underscored versions:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "metaweb_can_names = replace.(M.S, \" \" => \"_\")\n",
    "filter((x) -> (x) ∉ (metaweb_can_names ∩ treeleaves), metaweb_can_names)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "This prepare a data frame for the trait values:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "traitframe = DataFrame(; tipNames=treeleaves)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Infering the subspaces from the phylogeny"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "We will prepare the data frame required by `PhyloNetworks` to store the\n",
    "reconstructed traits:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "matching_tree_reconstruction = DataFrame(;\n",
    "    tipNames=treeleaves, nodeNumber=range(1, length(treeleaves); step=1)\n",
    ")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "To avoid writing long code, we prepare a series of columns with L and R as a\n",
    "prefix, and the number (from 1 to 12) of the corresponding dimension as a\n",
    "suffix:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "leftnames = \"L\" .* string.(1:size(L, 2))\n",
    "traits_L = DataFrame(L, Symbol.(leftnames))\n",
    "traits_L[!, \"tipNames\"] = metaweb_can_names\n",
    "\n",
    "rightnames = \"R\" .* string.(1:size(R, 1))\n",
    "traits_R = DataFrame(R', rightnames)\n",
    "traits_R[!, \"tipNames\"] = metaweb_can_names"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We can now merge these dataframes, and have something fully ready to be filled\n",
    "by the phylogenetic simulation:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "traits = leftjoin(traitframe, traits_L; on=:tipNames)\n",
    "traits = leftjoin(traits, traits_R; on=:tipNames)\n",
    "\n",
    "imputedtraits = DataFrame(;\n",
    "    tipNames=[treeleaves; fill(missing, tree_net.numNodes - tree_net.numTaxa)]\n",
    ")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "The last step is to reconstruct each trait (L1 to L12, R1 to R12) for the\n",
    "entire tree. This is by far the longest part of the script, but it is not\n",
    "terribly long. This could probably be made thread parallel fairly easily but,\n",
    "this would also take more time to do than it takes to run."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "for coord in 1:size(L, 2)\n",
    "    for prefix in [\"L\", \"R\"]\n",
    "        @info \"Reconstructing $(prefix)$(coord)\"\n",
    "        lower, upper, mean_trait = leaf_traits_reconstruction(\n",
    "            traits[!, [\"$(prefix)$(coord)\", \"tipNames\"]], tree_net\n",
    "        )\n",
    "        imputedtraits[!, \"$(prefix)$(coord)_low\"] = lower\n",
    "        imputedtraits[!, \"$(prefix)$(coord)_up\"] = upper\n",
    "        imputedtraits[!, \"$(prefix)$(coord)_mean\"] = mean_trait\n",
    "    end\n",
    "end"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "When the loop is done, we extract the values for the species in the Canadian\n",
    "species pool:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "canadian_rec = innerjoin(dropmissing(imputedtraits), pool; on=:tipNames)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Extracting the left and right subspaces from canada"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "This part of the script uses the same characters as the equations in the paper\n",
    "- if you do not have a typeface with good unicode mathematical support\n",
    "installed, you might not get the full effect."
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "Recall that a RDGP approximation is a matrix multiplication - we will\n",
    "therefore get the reconstructed average values after the Brownian motion\n",
    "model, and have a little look at them:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "𝓁 = Array(canadian_rec[!, leftnames .* \"_mean\"])\n",
    "𝓇 = transpose(Array(canadian_rec[!, rightnames .* \"_mean\"]))\n",
    "\n",
    "plot(\n",
    "    heatmap(𝓁; c=:GnBu, frame=:none, cbar=false),\n",
    "    heatmap(𝓇'; c=:BuPu, frame=:none, cbar=false);\n",
    "    dpi=600,\n",
    "    size=(500, 400),\n",
    ")\n",
    "savefig(\"figures/imputed-subspaces.png\")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Generating a probabilistic network"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "Because the Brownian motion model gives us a lower and upper bound, we will\n",
    "perform a series of random draws assuming that the values are uniformly\n",
    "distributed between these values. This step of the workflow can be adapted\n",
    "further. For example, one might want to account for the uncertainty in the\n",
    "phylogeny itself, or fit the distribution returned at each node rather than\n",
    "assuming a uniform distribution. We think that our approach introduces the\n",
    "least amount of guesses; it is likely to be over-estimating the chances of\n",
    "interactions a little, but this is the purpose of a metaweb: to give a list of\n",
    "possible interactions, to be later pared down."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "𝓁ₗ = Array(canadian_rec[!, leftnames .* \"_low\"])\n",
    "𝓇ₗ = transpose(Array(canadian_rec[!, rightnames .* \"_low\"]))\n",
    "\n",
    "𝓁ᵤ = Array(canadian_rec[!, leftnames .* \"_up\"])\n",
    "𝓇ᵤ = transpose(Array(canadian_rec[!, rightnames .* \"_up\"]))"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "The distributions are expressed as actual Uniform distributions from the\n",
    "`Distributions` package."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "ℒ = Matrix{Uniform}(undef, size(𝓁))\n",
    "for i in eachindex(ℒ)\n",
    "    ℒ[i] = Uniform(𝓁ₗ[i], 𝓁ᵤ[i])\n",
    "end\n",
    "\n",
    "ℛ = Matrix{Uniform}(undef, size(𝓇))\n",
    "for i in eachindex(ℛ)\n",
    "    ℛ[i] = Uniform(𝓇ₗ[i], 𝓇ᵤ[i])\n",
    "end"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We will do a large enough number of draws:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "draws = 20_000\n",
    "\n",
    "𝐋 = [rand.(ℒ) for i in 1:draws]\n",
    "𝐑 = [rand.(ℛ) for i in 1:draws]"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "There are two pieces of information to keep in mind here. The first is that a\n",
    "RDPG is a matrix multiplication, so we simply need to multiply the 20000\n",
    "random subspaces, to get 20000 random matrices. The second is that these\n",
    "matrices give results not in {0,1} but in ℝ, but we have estimated an optimal\n",
    "threshold for this projection. Doing all this is a one-liner:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "Ns = [(𝐋[i] * 𝐑[i]) .> threshold for i in 1:length(𝐋)]"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We can finally generate a probabilistic metaweb, in which the probability is\n",
    "defined as the proprtion of samples in which the interaction was inferred:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "P = UnipartiteProbabilisticNetwork(\n",
    "    reduce(.+, Ns) ./ draws, replace.(canadian_rec.tipNames, \"_\" => \" \")\n",
    ")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We can have a little look at the interactions sorted by probabilities:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "sort(interactions(P); by=(x) -> x.probability, rev=true)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Visualising the results"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "The next figures are very simple plots of the adjacency matrices:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "sporder = sortperm(vec(sum(adjacency(P); dims=2)))\n",
    "h1 = heatmap(\n",
    "    adjacency(P)[sporder, sporder];\n",
    "    c=:Greys,\n",
    "    frame=:none,\n",
    "    cbar=false,\n",
    "    dpi=600,\n",
    "    size=(500, 500),\n",
    "    aspectratio=1,\n",
    ")\n",
    "\n",
    "sporder = sortperm(vec(sum(adjacency(M); dims=2)))\n",
    "h2 = heatmap(\n",
    "    adjacency(M)[sporder, sporder];\n",
    "    c=:Greys,\n",
    "    frame=:none,\n",
    "    cbar=false,\n",
    "    dpi=600,\n",
    "    size=(500, 500),\n",
    "    aspectratio=1,\n",
    ")\n",
    "\n",
    "plot(h2, h1; size=(1000, 500))\n",
    "savefig(\"figures/adjacencymatrices.png\")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Writing output files for the raw predictions"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "We will store the results in a data frame - the information we care about is\n",
    "the probability for the species pair, and whether the pair was also found in\n",
    "Europe, and if so, whether it interacted:"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "output = DataFrame(; from=String[], to=String[], score=Float64[], pair=Bool[], int=Bool[])\n",
    "for int in interactions(P)\n",
    "    pair = (int.from in species(M)) & (int.to in species(M))\n",
    "    mint = pair ? M[int.from, int.to] : false\n",
    "    push!(output, (int.from, int.to, int.probability, pair, mint))\n",
    "end\n",
    "sort!(output, [:score, :from, :to]; rev=[true, false, false])"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "Save the basic network (no corrections)"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "CSV.write(\"artifacts/canadian_uncorrected.csv\", output)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Exploration of the relationship between subspaces and network properties"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "kout = degree(P; dims=1)\n",
    "kin = degree(P; dims=2)\n",
    "ordered_sp = replace.(canadian_rec.tipNames, \"_\" => \" \")\n",
    "\n",
    "scatter(\n",
    "    𝓁[:, 1], [kout[s] / richness(P) for s in ordered_sp]; dpi=600, size=(500, 500), lab=\"\"\n",
    ")\n",
    "xaxis!(\"Position in the left subspace\", extrema(vcat(𝓇', 𝓁)))\n",
    "yaxis!(\"Probabilistic generality\", (0, 1))\n",
    "savefig(\"figures/left-gen.png\")\n",
    "\n",
    "scatter(\n",
    "    𝓇'[:, 1],\n",
    "    [kin[s] / richness(P) for s in ordered_sp];\n",
    "    dpi=600,\n",
    "    size=(500, 500),\n",
    "    lab=\"\",\n",
    "    legend=:bottomright,\n",
    ")\n",
    "xaxis!(\"Position in the right subspace\", extrema(vcat(𝓇', 𝓁)))\n",
    "yaxis!(\"Probabilistic vulnerability\")\n",
    "savefig(\"figures/right-vuln.png\")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Basic corrections"
   ],
   "metadata": {}
  },
  {
   "cell_type": "markdown",
   "source": [
    "We will directly bring knowledge from the European metaweb, meaning that if\n",
    "two species interact in Europe, we assume they also do in Canada (remember,\n",
    "these are metawebs, we only care about the biological feasibility of the\n",
    "interaction); if two species do not interact in Europe, we prevent them from\n",
    "interacting in Canada - this later point could be reversed, by inflating the\n",
    "European metaweb using the simulation results, and whether to apply this step\n",
    "at all can be considered on a case by case basis."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "N = copy(P)\n",
    "shared_species = filter(s -> s in species(M), species(P))\n",
    "for s1 in shared_species\n",
    "    for s2 in shared_species\n",
    "        N[s1, s2] = M[s1, s2] ? 1.0 : 0.0\n",
    "    end\n",
    "end"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "We may have introduced a number of 0s in the sparse matrix, and it is good\n",
    "hygiene to remove them."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "SparseArrays.dropzeros!(N.edges)\n",
    "simplify!(N)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Writing the final metaweb"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "final = DataFrame(; from=String[], to=String[], score=Float64[])\n",
    "for int in interactions(N)\n",
    "    push!(final, (int.from, int.to, int.probability))\n",
    "end\n",
    "sort!(final, [:score, :from, :to]; rev=[true, false, false])"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "Save the corrected network"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "CSV.write(\"artifacts/canadian_corrected.csv\", final)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "## Plots for the core results"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "l = @layout [\n",
    "    a{0.3w} b\n",
    "    c{0.7h} d\n",
    "]\n",
    "\n",
    "sporder = sortperm(vec(sum(adjacency(N); dims=2)))\n",
    "\n",
    "plot(\n",
    "    plot(; legend=false, axes=false, frame=:none),\n",
    "    heatmap(𝓇[:, sporder]; frame=:none, legend=false, c=:BrBG, clim=(-1, 1)),\n",
    "    heatmap(𝓁[sporder, :]; frame=:none, legend=false, c=:PRGn, clim=(-1, 1)),\n",
    "    heatmap(adjacency(N)[sporder, sporder]; c=:Greys, frame=:none, legend=false);\n",
    "    layout=l,\n",
    "    size=(500, 500),\n",
    "    dpi=600,\n",
    ")\n",
    "\n",
    "savefig(\"figures/combined-prediction.png\")\n",
    "\n",
    "sporder = sortperm(vec(sum(adjacency(M); dims=2)))\n",
    "\n",
    "l = @layout [\n",
    "    a{0.3w} b\n",
    "    c{0.7h} d\n",
    "]\n",
    "\n",
    "plot(\n",
    "    plot(; legend=false, axes=false, frame=:none),\n",
    "    heatmap(R[:, sporder]; frame=:none, legend=false, c=:BrBG, clim=(-1, 1)),\n",
    "    heatmap(L[sporder, :]; frame=:none, legend=false, c=:PRGn, clim=(-1, 1)),\n",
    "    heatmap(adjacency(M)[sporder, sporder]; c=:Greys, frame=:none, legend=false);\n",
    "    layout=l,\n",
    "    size=(500, 500),\n",
    "    dpi=600,\n",
    ")\n",
    "\n",
    "savefig(\"figures/combined-empirical.png\")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "MaxEnt configuration model"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "a = zeros(Float64, size(adjacency(P)))\n",
    "C = UnipartiteProbabilisticNetwork(a, EcologicalNetworks._species_objects(P)...)\n",
    "\n",
    "for s1 in species(P; dims=1), s2 in species(P; dims=2)\n",
    "    C[s1, s2] = 0.5(kout[s1] / richness(P) + kin[s2] / richness(P))\n",
    "end\n",
    "\n",
    "rec = 𝓁[:, 1] * hcat(𝓇'[:, 1]...)\n",
    "\n",
    "Ceff = adjacency(C)\n",
    "Cinf = rec\n",
    "\n",
    "Zeff = (Ceff .- mean(Ceff)) ./ std(Ceff)\n",
    "Zinf = (Cinf .- mean(Cinf)) ./ std(Cinf)\n",
    "\n",
    "density(sqrt.(vec(Zeff .- Zinf) .^ 2.0); size=(500, 500), dpi=600, fill=(0, 0.2), lab=\"\")\n",
    "xaxis!(\"Mean squared error\", (0, 5))\n",
    "yaxis!(\"Density\", (0, 3))\n",
    "savefig(\"figures/distance-configuration.png\")\n",
    "\n",
    "sporder = sortperm(vec(sum(adjacency(C); dims=2)))\n",
    "\n",
    "ΔZ = Zeff .- Zinf\n",
    "\n",
    "heatmap(\n",
    "    ΔZ[sporder, sporder];\n",
    "    clim=(-3, 3),\n",
    "    c=:PuOr,\n",
    "    size=(500, 500),\n",
    "    dpi=600,\n",
    "    frame=:none,\n",
    "    aspectratio=1,\n",
    ")\n",
    "savefig(\"figures/heatmap-configuration.png\")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "---\n",
    "\n",
    "*This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*"
   ],
   "metadata": {}
  }
 ],
 "nbformat_minor": 3,
 "metadata": {
  "language_info": {
   "file_extension": ".jl",
   "mimetype": "application/julia",
   "name": "julia",
   "version": "1.6.2"
  },
  "kernelspec": {
   "name": "julia-1.6",
   "display_name": "Julia 1.6.2",
   "language": "julia"
  }
 },
 "nbformat": 4
}


# Step 5 - the actual prediction

This is the largest step in the entire pipeline, but not necessarily a complex
one. In short, this will generate the entire paper, minus some of the data
inflation and post-processing steps. As such, there are a few dependencies in
play.

````julia
using PhyloNetworks
using ProgressMeter
using LinearAlgebra
using EcologicalNetworks
using DelimitedFiles
using Random
using StatsPlots
using StatsBase
using Statistics
using SparseArrays
using CSV
using DataFrames
using GBIF
using Base.Threads
using StatsBase
using StatsModels
using Distributions
````

These lines are ensuring that the figures all look uniform, and we also set a
seed for reproducibility.

````julia
theme(:mute)
default(; frame=:box)
Random.seed!(01189998819991197253)
````

Some functions that are required for phylogenetic imputation are stored in
another file, and we will load them here to be done with that.

````julia
include("lib/pwar.jl")
````

## Reading the data pieces

Reading the tree is done with `PhyloNetworks` this time -- this is where we
will get the Brownian motion code from:

````julia
tree_net = readTopology(joinpath("data", "mammals.newick"));
````

We will grab back the corrected Europe/tree/GBIF names, as they will be used
quite a lot to match simulations to species names:

````julia
namelist = DataFrame(CSV.File(joinpath("artifacts", "names_metaweb_tree_gbif.csv")))
````

Finally, we grab the European metaweb *not* from the original file, but from
our edgelist artifact.

````julia
eurometa = readdlm(joinpath("artifacts", "europeanmetaweb.csv"), ',', String)
mwspecies = unique(eurometa)
M = UnipartiteNetwork(zeros(Bool, length(mwspecies), length(mwspecies)), mwspecies)
for i in 1:size(eurometa, 1)
    M[eurometa[i, :]...] = true
end
````

## Finding the rank for truncation

The ideal way to find a cutoff for the rank at which the matrix should be cut
would involve the profile likelihood, or approximating the maximum curvature
point on the screeplot using the central difference approximation for the
second order partial derivative. Sadly all of these approaches will conclude
that only the first dimension is required (there is a biological reason for
this, to which we will return when we dig into the biological meaning of the
latent variables).

The only eigenvalues we care about are the ones that are up to the rank of the
adjacency matrix - we will extract them, and range them so that they sum to
unity. It makes no difference in the results, and allows us to read the
proportion of variance explained directly.

````julia
rnk = rank(adjacency(M))
eig = svd(M).S[1:rnk]
neig = eig ./ sum(eig)
````

The first piece of information is a screeplot of the eigenvalues. As a rule,
the document rendered from the literate files will not include the figures, as
we do not want to go into the analysis (there is, after all, a whole paper for
that).

````julia
scatter(eig; lab="", dpi=600, size=(500, 500))
xaxis!("Dimension", (1, 38))
yaxis!("Eigenvalue", (0, 40))
savefig("figures/screeplot.png")
````

The second diagnosis plot is the proportion of variance explained. We added
lines marking the points that explain between 50% and 90% of the total
variance, every 10%.

````julia
scatter(cumsum(neig); lab="", dpi=600, size=(500, 500))
for ve in [0.5, 0.6, 0.7, 0.8, 0.9]
    i = findfirst(cumsum(neig) .> ve)
    x = cumsum(neig)[i]
    plot!([0, i, NaN, i, i], [x, x, NaN, x, 0]; lab="", c=:grey, ls=:dash)
end
xaxis!("Dimension", (1, 38))
yaxis!("Cumulative variance explained", (0, 1.0))
savefig("figures/varexplained.png")
````

At this point, we need to make an arbitrary call, which is to say to decide on
a proportion of variance explained. We decided on 60%, which means that the
first 12 ranks will be used.

## Embedding the European metaweb

````julia
L, R = rdpg(M, 12)
````

The `L` and `R` variables are the left and right subspaces from the random dot
product graph (which in turn are the left and right subspace from the t-SVD
multiplied by the square root of the diagonal matrix containing the
eigenvalues). Multiplying `L` and `R` (using a matrix multiplication!) will
give a rank-12 approximation of the European metaweb.

## Thresholding the embedded network

The reconstucted network (`L*R`) will not have Boolean values. Essentially,
the multiplication will approximate something with values in {0,1} by
something with values in ℝ, and so we need to find a cutoff to separate
interactions from non-interactions. This threshold is important because it
holds for any reconstruction made using these latent variables - that is to
say, when we infer the Canadian interactions using reconstrusted left and
right subspaces, we will apply the same threshold.

Our basic approach to threshdolding is to examine 500 cutoffs points, get the
confusion matrix, and retain as a threshold the value that maximizes
Informedness:

````julia
A = adjacency(M)
thresholds = LinRange(extrema(L * R)..., 500)
tp = zeros(Float64, length(thresholds))
tn = similar(tp)
fp = similar(tp)
fn = similar(tp)
for (i, t) in enumerate(thresholds)
    PN = (L * R) .>= t
    tp[i] = sum((A) .& (PN)) / sum(A)
    tn[i] = sum((.!(A)) .& (.!(PN))) / sum(.!(A))
    fp[i] = sum((.!(A)) .& (PN)) / sum(.!(A))
    fn[i] = sum((A) .& (.!(PN))) / sum(A)
end
````

We can calculate Youden's J at every cutoff, and the threshold is going to be
the threshold associated to the maximum Y:

````julia
Y = tp .+ tn .- 1
maxY, posY = findmax(Y)
threshold = thresholds[posY]
````

We can plot the optimal cutoff - this is not presented in the manuscript, but
useful in case you need to satisfy your curiosity. The best cutoff gives an
informedness essentially = 1, which is expected.

````julia
plot(thresholds, Y; lw=0.0, fill=(0, 0.2), lab="", dpi=600, size=(500, 500))
scatter!([threshold], [maxY]; lab="", c=:darkgrey, msw=0.0)
yaxis!("Youden's J", (0, 1))
xaxis!("Cutoff", extrema(L * R))
savefig("figures/optimalcutoff.png")
````

## Visual examination of the subspaces

The left and right subspaces *do* hold ecological information, and so it is a
good idea to check them out visually. One striking result is that species with
a value of 0 in the left subspace also have no preys: this is a strong clue
that the left subspace is associated to generality (in the sense of Schoener
1989), a fact we will epxloit later on.

````julia
plot(
    heatmap(L; c=:GnBu, frame=:none, cbar=false),
    heatmap(R'; c=:BuPu, frame=:none, cbar=false);
    dpi=600,
    size=(500, 400),
)
savefig("figures/subspaces.png")
````

## Preparing the tree for inference

We are now ready to move on to the next step: infering the values of the left
and right subspaces for the species that are not in the European metaweb, but
are in the tree. To do so, we will need the tip names:

````julia
treeleaves = tipLabels(tree_net)
````

We are now ready to match the tree to the Canadian species pool names:

````julia
canada = DataFrame(CSV.File(joinpath("artifacts", "iucn_gbif_names.csv")))
cancodes = replace.(unique(filter(!ismissing, canada.gbifname)), " " => "_")
tree_cleanup = DataFrame(CSV.File(joinpath("artifacts", "upham_gbif_names.csv")))

csp = dropmissing!(DataFrame(; gbifname=canada.gbifname))
csp = dropmissing(leftjoin(csp, tree_cleanup; on=:gbifname))
````

We can store the cleaned Canadian species name into a data frame - this is
mandated by the `PhyloNetworks` interface:

````julia
canmammals = unique(csp.code)
pool = DataFrame(; tipNames=canmammals)
````

We finally convert the network names to their underscored versions:

````julia
metaweb_can_names = replace.(M.S, " " => "_")
filter((x) -> (x) ∉ (metaweb_can_names ∩ treeleaves), metaweb_can_names)
````

This prepare a data frame for the trait values:

````julia
traitframe = DataFrame(; tipNames=treeleaves)
````

## Infering the subspaces from the phylogeny

We will prepare the data frame required by `PhyloNetworks` to store the
reconstructed traits:

````julia
matching_tree_reconstruction = DataFrame(;
    tipNames=treeleaves, nodeNumber=range(1, length(treeleaves); step=1)
)
````

To avoid writing long code, we prepare a series of columns with L and R as a
prefix, and the number (from 1 to 12) of the corresponding dimension as a
suffix:

````julia
leftnames = "L" .* string.(1:size(L, 2))
traits_L = DataFrame(L, Symbol.(leftnames))
traits_L[!, "tipNames"] = metaweb_can_names

rightnames = "R" .* string.(1:size(R, 1))
traits_R = DataFrame(R', rightnames)
traits_R[!, "tipNames"] = metaweb_can_names
````

We can now merge these dataframes, and have something fully ready to be filled
by the phylogenetic simulation:

````julia
traits = leftjoin(traitframe, traits_L; on=:tipNames)
traits = leftjoin(traits, traits_R; on=:tipNames)

imputedtraits = DataFrame(;
    tipNames=[treeleaves; fill(missing, tree_net.numNodes - tree_net.numTaxa)]
)
````

The last step is to reconstruct each trait (L1 to L12, R1 to R12) for the
entire tree. This is by far the longest part of the script, but it is not
terribly long. This could probably be made thread parallel fairly easily but,
this would also take more time to do than it takes to run.

````julia
for coord in 1:size(L, 2)
    for prefix in ["L", "R"]
        @info "Reconstructing $(prefix)$(coord)"
        lower, upper, mean_trait = leaf_traits_reconstruction(
            traits[!, ["$(prefix)$(coord)", "tipNames"]], tree_net
        )
        imputedtraits[!, "$(prefix)$(coord)_low"] = lower
        imputedtraits[!, "$(prefix)$(coord)_up"] = upper
        imputedtraits[!, "$(prefix)$(coord)_mean"] = mean_trait
    end
end
````

When the loop is done, we extract the values for the species in the Canadian
species pool:

````julia
canadian_rec = innerjoin(dropmissing(imputedtraits), pool; on=:tipNames)
````

## Extracting the left and right subspaces from canada

This part of the script uses the same characters as the equations in the paper
- if you do not have a typeface with good unicode mathematical support
installed, you might not get the full effect.

Recall that a RDGP approximation is a matrix multiplication - we will
therefore get the reconstructed average values after the Brownian motion
model, and have a little look at them:

````julia
𝓁 = Array(canadian_rec[!, leftnames .* "_mean"])
𝓇 = transpose(Array(canadian_rec[!, rightnames .* "_mean"]))

plot(
    heatmap(𝓁; c=:GnBu, frame=:none, cbar=false),
    heatmap(𝓇'; c=:BuPu, frame=:none, cbar=false);
    dpi=600,
    size=(500, 400),
)
savefig("figures/imputed-subspaces.png")
````

## Generating a probabilistic network

Because the Brownian motion model gives us a lower and upper bound, we will
perform a series of random draws assuming that the values are uniformly
distributed between these values. This step of the workflow can be adapted
further. For example, one might want to account for the uncertainty in the
phylogeny itself, or fit the distribution returned at each node rather than
assuming a uniform distribution. We think that our approach introduces the
least amount of guesses; it is likely to be over-estimating the chances of
interactions a little, but this is the purpose of a metaweb: to give a list of
possible interactions, to be later pared down.

````julia
𝓁ₗ = Array(canadian_rec[!, leftnames .* "_low"])
𝓇ₗ = transpose(Array(canadian_rec[!, rightnames .* "_low"]))

𝓁ᵤ = Array(canadian_rec[!, leftnames .* "_up"])
𝓇ᵤ = transpose(Array(canadian_rec[!, rightnames .* "_up"]))
````

The distributions are expressed as actual Uniform distributions from the
`Distributions` package.

````julia
ℒ = Matrix{Uniform}(undef, size(𝓁))
for i in eachindex(ℒ)
    ℒ[i] = Uniform(𝓁ₗ[i], 𝓁ᵤ[i])
end

ℛ = Matrix{Uniform}(undef, size(𝓇))
for i in eachindex(ℛ)
    ℛ[i] = Uniform(𝓇ₗ[i], 𝓇ᵤ[i])
end
````

We will do a large enough number of draws:

````julia
draws = 20_000

𝐋 = [rand.(ℒ) for i in 1:draws]
𝐑 = [rand.(ℛ) for i in 1:draws]
````

There are two pieces of information to keep in mind here. The first is that a
RDPG is a matrix multiplication, so we simply need to multiply the 20000
random subspaces, to get 20000 random matrices. The second is that these
matrices give results not in {0,1} but in ℝ, but we have estimated an optimal
threshold for this projection. Doing all this is a one-liner:

````julia
Ns = [(𝐋[i] * 𝐑[i]) .> threshold for i in 1:length(𝐋)]
````

We can finally generate a probabilistic metaweb, in which the probability is
defined as the proprtion of samples in which the interaction was inferred:

````julia
P = UnipartiteProbabilisticNetwork(
    reduce(.+, Ns) ./ draws, replace.(canadian_rec.tipNames, "_" => " ")
)
````

We can have a little look at the interactions sorted by probabilities:

````julia
sort(interactions(P); by=(x) -> x.probability, rev=true)
````

## Visualising the results

The next figures are very simple plots of the adjacency matrices:

````julia
sporder = sortperm(vec(sum(adjacency(P); dims=2)))
h1 = heatmap(
    adjacency(P)[sporder, sporder];
    c=:Greys,
    frame=:none,
    cbar=false,
    dpi=600,
    size=(500, 500),
    aspectratio=1,
)

sporder = sortperm(vec(sum(adjacency(M); dims=2)))
h2 = heatmap(
    adjacency(M)[sporder, sporder];
    c=:Greys,
    frame=:none,
    cbar=false,
    dpi=600,
    size=(500, 500),
    aspectratio=1,
)

plot(h2, h1; size=(1000, 500))
savefig("figures/adjacencymatrices.png")
````

## Writing output files for the raw predictions

We will store the results in a data frame - the information we care about is
the probability for the species pair, and whether the pair was also found in
Europe, and if so, whether it interacted:

````julia
output = DataFrame(; from=String[], to=String[], score=Float64[], pair=Bool[], int=Bool[])
for int in interactions(P)
    pair = (int.from in species(M)) & (int.to in species(M))
    mint = pair ? M[int.from, int.to] : false
    push!(output, (int.from, int.to, int.probability, pair, mint))
end
sort!(output, [:score, :from, :to]; rev=[true, false, false])
````

Save the basic network (no corrections)

````julia
CSV.write("artifacts/canadian_uncorrected.csv", output)
````

## Exploration of the relationship between subspaces and network properties

````julia
kout = degree(P; dims=1)
kin = degree(P; dims=2)
ordered_sp = replace.(canadian_rec.tipNames, "_" => " ")

scatter(
    𝓁[:, 1], [kout[s] / richness(P) for s in ordered_sp]; dpi=600, size=(500, 500), lab=""
)
xaxis!("Position in the left subspace", extrema(vcat(𝓇', 𝓁)))
yaxis!("Probabilistic generality", (0, 1))
savefig("figures/left-gen.png")

scatter(
    𝓇'[:, 1],
    [kin[s] / richness(P) for s in ordered_sp];
    dpi=600,
    size=(500, 500),
    lab="",
    legend=:bottomright,
)
xaxis!("Position in the right subspace", extrema(vcat(𝓇', 𝓁)))
yaxis!("Probabilistic vulnerability")
savefig("figures/right-vuln.png")
````

## Basic corrections

We will directly bring knowledge from the European metaweb, meaning that if
two species interact in Europe, we assume they also do in Canada (remember,
these are metawebs, we only care about the biological feasibility of the
interaction); if two species do not interact in Europe, we prevent them from
interacting in Canada - this later point could be reversed, by inflating the
European metaweb using the simulation results, and whether to apply this step
at all can be considered on a case by case basis.

````julia
N = copy(P)
shared_species = filter(s -> s in species(M), species(P))
for s1 in shared_species
    for s2 in shared_species
        N[s1, s2] = M[s1, s2] ? 1.0 : 0.0
    end
end
````

We may have introduced a number of 0s in the sparse matrix, and it is good
hygiene to remove them.

````julia
SparseArrays.dropzeros!(N.edges)
simplify!(N)
````

## Writing the final metaweb

````julia
final = DataFrame(; from=String[], to=String[], score=Float64[])
for int in interactions(N)
    push!(final, (int.from, int.to, int.probability))
end
sort!(final, [:score, :from, :to]; rev=[true, false, false])
````

Save the corrected network

````julia
CSV.write("artifacts/canadian_corrected.csv", final)
````

## Plots for the core results

````julia
l = @layout [
    a{0.3w} b
    c{0.7h} d
]

sporder = sortperm(vec(sum(adjacency(N); dims=2)))

plot(
    plot(; legend=false, axes=false, frame=:none),
    heatmap(𝓇[:, sporder]; frame=:none, legend=false, c=:BrBG, clim=(-1, 1)),
    heatmap(𝓁[sporder, :]; frame=:none, legend=false, c=:PRGn, clim=(-1, 1)),
    heatmap(adjacency(N)[sporder, sporder]; c=:Greys, frame=:none, legend=false);
    layout=l,
    size=(500, 500),
    dpi=600,
)

savefig("figures/combined-prediction.png")

sporder = sortperm(vec(sum(adjacency(M); dims=2)))

l = @layout [
    a{0.3w} b
    c{0.7h} d
]

plot(
    plot(; legend=false, axes=false, frame=:none),
    heatmap(R[:, sporder]; frame=:none, legend=false, c=:BrBG, clim=(-1, 1)),
    heatmap(L[sporder, :]; frame=:none, legend=false, c=:PRGn, clim=(-1, 1)),
    heatmap(adjacency(M)[sporder, sporder]; c=:Greys, frame=:none, legend=false);
    layout=l,
    size=(500, 500),
    dpi=600,
)

savefig("figures/combined-empirical.png")
````

MaxEnt configuration model

````julia
a = zeros(Float64, size(adjacency(P)))
C = UnipartiteProbabilisticNetwork(a, EcologicalNetworks._species_objects(P)...)

for s1 in species(P; dims=1), s2 in species(P; dims=2)
    C[s1, s2] = 0.5(kout[s1] / richness(P) + kin[s2] / richness(P))
end

rec = 𝓁[:, 1] * hcat(𝓇'[:, 1]...)

Ceff = adjacency(C)
Cinf = rec

Zeff = (Ceff .- mean(Ceff)) ./ std(Ceff)
Zinf = (Cinf .- mean(Cinf)) ./ std(Cinf)

density(sqrt.(vec(Zeff .- Zinf) .^ 2.0); size=(500, 500), dpi=600, fill=(0, 0.2), lab="")
xaxis!("Mean squared error", (0, 5))
yaxis!("Density", (0, 3))
savefig("figures/distance-configuration.png")

sporder = sortperm(vec(sum(adjacency(C); dims=2)))

ΔZ = Zeff .- Zinf

heatmap(
    ΔZ[sporder, sporder];
    clim=(-3, 3),
    c=:PuOr,
    size=(500, 500),
    dpi=600,
    frame=:none,
    aspectratio=1,
)
savefig("figures/heatmap-configuration.png")
````



{
 "cells": [
  {
   "cell_type": "markdown",
   "source": [
    "# Step 6 - inflating the predictions with the Newfoundland data"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "using DelimitedFiles\n",
    "using DataFrames\n",
    "using CSV: CSV\n",
    "using GBIF\n",
    "using NCBITaxonomy: NCBITaxonomy\n",
    "using EcologicalNetworks"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "..."
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "sl_raw = readdlm(\"data/NLfoodweb.csv\", ',')\n",
    "\n",
    "sl_sp = replace.(sl_raw[1, 2:end], \".\" => \" \")\n",
    "sl_A = Bool.(sl_raw[2:end, 2:end])\n",
    "\n",
    "#%% Correct the species names from Strong & Leroux to make them match the GBIF taxonomy\n",
    "nf = NCBITaxonomy.mammalfilter(true);\n",
    "scinames = Dict{String,String}()\n",
    "for s in sl_sp\n",
    "    t = NCBITaxonomy.taxon(s; strict=false)\n",
    "    if !isnothing(t)\n",
    "        scinames[s] = t.name\n",
    "    end\n",
    "end\n",
    "\n",
    "#%% Convert the original names in GBIF names through the matched NCBI names\n",
    "valnames = Dict{String,String}()\n",
    "for (s, t) in scinames\n",
    "    gbifmatch = GBIF.taxon(t; strict=false)\n",
    "    if !isnothing(gbifmatch)\n",
    "        if !ismissing(gbifmatch.species)\n",
    "            if gbifmatch.class.first == \"Mammalia\"\n",
    "                valnames[s] = gbifmatch.species.first\n",
    "            end\n",
    "        end\n",
    "    end\n",
    "end\n",
    "\n",
    "#%% Get the correct keys in the original metaweb\n",
    "idxmatch = findall(x -> x in keys(valnames), sl_sp)\n",
    "\n",
    "#%% Assemble the network\n",
    "spnames = [valnames[s] for s in sl_sp[idxmatch]]\n",
    "A = sl_A[idxmatch, idxmatch]'\n",
    "NL = UnipartiteNetwork(A, spnames)\n",
    "\n",
    "#%% Save as a CSV\n",
    "df = DataFrame(; from=String[], to=String[])\n",
    "for i in interactions(NL)\n",
    "    push!(df, (i.from, i.to))\n",
    "end\n",
    "CSV.write(\"artifacts/newfoundland.csv\", df)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "---\n",
    "\n",
    "*This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*"
   ],
   "metadata": {}
  }
 ],
 "nbformat_minor": 3,
 "metadata": {
  "language_info": {
   "file_extension": ".jl",
   "mimetype": "application/julia",
   "name": "julia",
   "version": "1.6.2"
  },
  "kernelspec": {
   "name": "julia-1.6",
   "display_name": "Julia 1.6.2",
   "language": "julia"
  }
 },
 "nbformat": 4
}


# Step 6 - inflating the predictions with the Newfoundland data

````julia
using DelimitedFiles
using DataFrames
using CSV: CSV
using GBIF
using NCBITaxonomy: NCBITaxonomy
using EcologicalNetworks
````

...

````julia
sl_raw = readdlm("data/NLfoodweb.csv", ',')

sl_sp = replace.(sl_raw[1, 2:end], "." => " ")
sl_A = Bool.(sl_raw[2:end, 2:end])

#%% Correct the species names from Strong & Leroux to make them match the GBIF taxonomy
nf = NCBITaxonomy.mammalfilter(true);
scinames = Dict{String,String}()
for s in sl_sp
    t = NCBITaxonomy.taxon(s; strict=false)
    if !isnothing(t)
        scinames[s] = t.name
    end
end

#%% Convert the original names in GBIF names through the matched NCBI names
valnames = Dict{String,String}()
for (s, t) in scinames
    gbifmatch = GBIF.taxon(t; strict=false)
    if !isnothing(gbifmatch)
        if !ismissing(gbifmatch.species)
            if gbifmatch.class.first == "Mammalia"
                valnames[s] = gbifmatch.species.first
            end
        end
    end
end

#%% Get the correct keys in the original metaweb
idxmatch = findall(x -> x in keys(valnames), sl_sp)

#%% Assemble the network
spnames = [valnames[s] for s in sl_sp[idxmatch]]
A = sl_A[idxmatch, idxmatch]'
NL = UnipartiteNetwork(A, spnames)

#%% Save as a CSV
df = DataFrame(; from=String[], to=String[])
for i in interactions(NL)
    push!(df, (i.from, i.to))
end
CSV.write("artifacts/newfoundland.csv", df)
````



{
 "cells": [
  {
   "cell_type": "markdown",
   "source": [
    "# Step 7 - inflating the predictions with the GLOBI data"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "#%% Get the dependencies\n",
    "using HTTP\n",
    "using ProgressMeter\n",
    "using JSON\n",
    "using DataFrames\n",
    "using CSV: CSV\n",
    "using StatsPlots\n",
    "\n",
    "theme(:mute)\n",
    "default(; frame=:box)\n",
    "\n",
    "#%% API info\n",
    "_globi_api = \"https://api.globalbioticinteractions.org/taxon\"\n",
    "relevant_types = [\"eats\", \"preysOn\"]\n",
    "\n",
    "#%% Read the corrected Canadian metaweb\n",
    "canmet = DataFrame(CSV.File(\"artifacts/canadian_corrected.csv\"))\n",
    "\n",
    "#%% Get the data from GLOBI\n",
    "allsp = unique(vcat(canmet.from, canmet.to))\n",
    "\n",
    "diet = DataFrame(; from=String[], to=String[])\n",
    "@showprogress for sp in allsp\n",
    "    url = \"$(_globi_api)/$(sp)/eats/\"\n",
    "    r = HTTP.request(\"GET\", url)\n",
    "    globidiet = JSON.parse(String(r.body))\n",
    "    if !isempty(globidiet[\"data\"])\n",
    "        for intlist in globidiet[\"data\"]\n",
    "            if intlist[2] in relevant_types\n",
    "                for s in intlist[3]\n",
    "                    if s in allsp\n",
    "                        if s != sp\n",
    "                            push!(diet, (sp, s))\n",
    "                        end\n",
    "                    end\n",
    "                end\n",
    "            end\n",
    "        end\n",
    "    end\n",
    "end\n",
    "diet = unique(diet)\n",
    "\n",
    "sort!(diet, [:from, :to])\n",
    "CSV.write(\"artifacts/globi_diet.csv\", diet)\n",
    "\n",
    "#%% Get the data we missed from GLOBI\n",
    "intcode = canmet.from .* canmet.to\n",
    "diet.intcode = diet.from .* diet.to\n",
    "\n",
    "missedint = select(diet[findall([!(d in intcode) for d in diet.intcode]), :], [:from, :to])\n",
    "sort!(missedint, [:from, :to])\n",
    "CSV.write(\"artifacts/globi_newinteractions.csv\", missedint)\n",
    "\n",
    "matchedint = dropmissing(leftjoin(diet, canmet; on=[:from => :from, :to => :to]))\n",
    "\n",
    "@info \"GLOBI: found $(size(missedint, 1)) new interactions out of $(size(diet, 1))\"\n",
    "@info \"GLOBI: $(length(unique(vcat(diet.from, diet.to)))) species\"\n",
    "\n",
    "#%% Get the data we missed from Strong & Leroux\n",
    "canmetsp = unique(vcat(canmet.from, canmet.to))\n",
    "sl = DataFrame(CSV.File(\"artifacts/newfoundland.csv\"))\n",
    "slshared = intersect(canmetsp, unique(vcat(sl.from, sl.to)))\n",
    "\n",
    "slkeep = [(i.from in slshared) & (i.to in slshared) for i in eachrow(sl)]\n",
    "sl = sl[findall(slkeep), :]\n",
    "\n",
    "sl.intcode = sl.from .* sl.to\n",
    "\n",
    "missedslint = select(sl[findall([!(d in intcode) for d in sl.intcode]), :], [:from, :to])\n",
    "\n",
    "sort!(missedslint, [:from, :to])\n",
    "CSV.write(\"artifacts/newfoundland_newinteractions.csv\", missedslint)\n",
    "\n",
    "matchedslint = dropmissing(leftjoin(sl, canmet; on=[:from => :from, :to => :to]))\n",
    "\n",
    "density(matchedint.score; dpi=600, size=(500, 500), lab=\"GLOBI\", lw=2.0)\n",
    "density!(matchedslint.score; lab=\"Newfoundland\", lw=2.0)\n",
    "density!(canmet.score; lw=0.0, fill=(0.2, 0), c=:black, lab=\"All predictions\")\n",
    "xaxis!(\"Imputed probability\", (0, 1))\n",
    "yaxis!(\"Density\", (0, 8))\n",
    "\n",
    "savefig(\"figures/inflation-comparison.png\")\n",
    "\n",
    "@info \"NFLD : found $(size(missedslint, 1)) new interactions out of $(size(sl, 1))\"\n",
    "@info \"NFLD : $(length(unique(vcat(sl.from, sl.to)))) species\"\n",
    "\n",
    "#%% Merge all interactions\n",
    "aug = leftjoin(unique(vcat(missedint, missedslint)), canmet; on=[:from, :to])\n",
    "aug.score .= 1.0\n",
    "\n",
    "inflated = vcat(canmet, aug)\n",
    "\n",
    "sort!(inflated, [:score, :from, :to]; rev=[true, false, false])"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "Save the corrected network"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "CSV.write(\"artifacts/canadian_inflated.csv\", inflated)"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "---\n",
    "\n",
    "*This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*"
   ],
   "metadata": {}
  }
 ],
 "nbformat_minor": 3,
 "metadata": {
  "language_info": {
   "file_extension": ".jl",
   "mimetype": "application/julia",
   "name": "julia",
   "version": "1.6.2"
  },
  "kernelspec": {
   "name": "julia-1.6",
   "display_name": "Julia 1.6.2",
   "language": "julia"
  }
 },
 "nbformat": 4
}


# Step 7 - inflating the predictions with the GLOBI data

````julia
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

#%% Merge all interactions
aug = leftjoin(unique(vcat(missedint, missedslint)), canmet; on=[:from, :to])
aug.score .= 1.0

inflated = vcat(canmet, aug)

sort!(inflated, [:score, :from, :to]; rev=[true, false, false])
````

Save the corrected network

````julia
CSV.write("artifacts/canadian_inflated.csv", inflated)
````



{
 "cells": [
  {
   "cell_type": "markdown",
   "source": [
    "# Step 8 - finding the interaction cutoff for the final results"
   ],
   "metadata": {}
  },
  {
   "outputs": [],
   "cell_type": "code",
   "source": [
    "#%% Dependencies\n",
    "using SparseArrays\n",
    "using EcologicalNetworks\n",
    "using DataFrames\n",
    "using CSV: CSV\n",
    "using StatsPlots\n",
    "\n",
    "#%% Update the theme defaults\n",
    "theme(:mute)\n",
    "default(; frame=:box)\n",
    "\n",
    "#%% Load the metaweb\n",
    "Mij = DataFrame(CSV.File(\"artifacts/canadian_inflated.csv\"))\n",
    "Si = unique(vcat(Mij.from, Mij.to))\n",
    "\n",
    "P = UnipartiteProbabilisticNetwork(zeros(Float64, (length(Si), length(Si))), Si)\n",
    "\n",
    "for int in eachrow(Mij)\n",
    "    P[int.from, int.to] = int.score\n",
    "end\n",
    "\n",
    "#%% Look at the cutoff for kept interactions\n",
    "ρ = LinRange(extrema(P.edges.nzval)..., 500)\n",
    "U = zeros(Float64, length(ρ))\n",
    "L = zeros(Float64, length(ρ))\n",
    "S = zeros(Float64, length(ρ))\n",
    "\n",
    "for (i, cutoff) in enumerate(ρ)\n",
    "    kept = P.edges.nzval .≥ cutoff\n",
    "    U[i] = sum(kept) / length(kept)\n",
    "    L[i] = sum(P.edges.nzval[kept]) / links(P)\n",
    "    S[i] = richness(simplify(P ≥ cutoff)) / richness(P)\n",
    "end\n",
    "\n",
    "#%% Central difference as a proxy for second derivative\n",
    "∂U = zeros(Float64, length(U))\n",
    "∂L = zeros(Float64, length(U))\n",
    "for i in 2:(length(U) - 1)\n",
    "    ∂U[i] = U[i + 1] + U[i - 1] - 2U[i]\n",
    "    ∂L[i] = L[i + 1] + L[i - 1] - 2L[i]\n",
    "end\n",
    "\n",
    "#%% Euclidean distance for U\n",
    "dU = zeros(Float64, length(U))\n",
    "for i in eachindex(U)\n",
    "    dU[i] = sqrt(U[i] * U[i] + ρ[i] * ρ[i])\n",
    "end\n",
    "\n",
    "#%% Plot the results\n",
    "plot(ρ, U; dpi=600, size=(500, 500), lab=\"Non-zero\")\n",
    "plot!(ρ, L; lab=\"Expected\")\n",
    "xaxis!(\"Cutoff\", (0, 1))\n",
    "yaxis!(\"Proportion of links left\", (0, 1))\n",
    "vline!([ρ[findlast(S .== 1)]]; c=:grey, ls=:dash, lab=\"\")\n",
    "\n",
    "savefig(\"figures/cutoff-interactions.png\")\n",
    "\n",
    "plot(ρ, 1.0 .- S; dpi=600, size=(500, 500), lab=\"Disconnected species\", legend=:topleft)\n",
    "l = L .* links(P)\n",
    "s = S .* richness(P)\n",
    "plot!(ρ, l ./ (s .* s); lab=\"Connectance\")\n",
    "xaxis!(\"Cutoff\", (0, 1))\n",
    "yaxis!(\"  \", (0, 0.5))\n",
    "vline!([ρ[findlast(S .== 1)]]; c=:grey, ls=:dash, lab=\"\")\n",
    "\n",
    "savefig(\"figures/cutoff-connectance.png\")\n",
    "\n",
    "#%% Get the threshold\n",
    "thrind = findlast(S .== 1)\n",
    "@info \"Optimal cutoff based on remaining species: $(ρ[thrind])\"\n",
    "@info \"Optimal cutoff based on central differences: $(ρ[last(findmax(∂U))])\"\n",
    "\n",
    "#%% Cleaned-up network\n",
    "K = copy(P)\n",
    "K.edges[P.edges .< ρ[thrind]] .= 0.0\n",
    "dropzeros!(K.edges)\n",
    "\n",
    "int = DataFrame(; from=String[], to=String[], score=Float64[])\n",
    "for i in interactions(K)\n",
    "    push!(int, (i.from, i.to, i.probability))\n",
    "end\n",
    "sort!(int, [:score, :from, :to]; rev=[true, false, false])\n",
    "\n",
    "CSV.write(\"artifacts/canadian_thresholded.csv\", int)\n",
    "\n",
    "#%% Write the functional classification of species\n",
    "rls = DataFrame(;\n",
    "    sp=String[],\n",
    "    gen=Float64[],\n",
    "    gen_var=Float64[],\n",
    "    vul=Float64[],\n",
    "    vul_var=Float64[],\n",
    "    role=Symbol[],\n",
    ")\n",
    "\n",
    "din = degree(K; dims=2)\n",
    "dinv = degree_var(K; dims=2)\n",
    "dout = degree(K; dims=1)\n",
    "doutv = degree_var(K; dims=1)\n",
    "\n",
    "for s in species(K)\n",
    "    rl = :int\n",
    "    iszero(din[s]) && (rl = :top)\n",
    "    iszero(dout[s]) && (rl = :prd)\n",
    "    push!(rls, (s, dout[s], doutv[s], din[s], dinv[s], rl))\n",
    "end\n",
    "\n",
    "sort!(rls, :gen; rev=true)\n",
    "CSV.write(\"artifacts/species_roles.csv\", rls)\n",
    "\n",
    "plot(log1p.(sort(rls.gen; rev=true)); lab=\"Out-degree\", size=(500, 500), dpi=600)\n",
    "plot!(log1p.(sort(rls.vul; rev=true)); lab=\"In-degree\")\n",
    "xaxis!(\"Rank\", (0, 180))\n",
    "yaxis!(\"log(degree + 1)\", (0, 5))\n",
    "\n",
    "savefig(\"figures/final-degree.png\")"
   ],
   "metadata": {},
   "execution_count": null
  },
  {
   "cell_type": "markdown",
   "source": [
    "---\n",
    "\n",
    "*This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*"
   ],
   "metadata": {}
  }
 ],
 "nbformat_minor": 3,
 "metadata": {
  "language_info": {
   "file_extension": ".jl",
   "mimetype": "application/julia",
   "name": "julia",
   "version": "1.6.2"
  },
  "kernelspec": {
   "name": "julia-1.6",
   "display_name": "Julia 1.6.2",
   "language": "julia"
  }
 },
 "nbformat": 4
}


# Step 8 - finding the interaction cutoff for the final results

````julia
#%% Dependencies
using SparseArrays
using EcologicalNetworks
using DataFrames
using CSV: CSV
using StatsPlots

#%% Update the theme defaults
theme(:mute)
default(; frame=:box)

#%% Load the metaweb
Mij = DataFrame(CSV.File("artifacts/canadian_inflated.csv"))
Si = unique(vcat(Mij.from, Mij.to))

P = UnipartiteProbabilisticNetwork(zeros(Float64, (length(Si), length(Si))), Si)

for int in eachrow(Mij)
    P[int.from, int.to] = int.score
end

#%% Look at the cutoff for kept interactions
ρ = LinRange(extrema(P.edges.nzval)..., 500)
U = zeros(Float64, length(ρ))
L = zeros(Float64, length(ρ))
S = zeros(Float64, length(ρ))

for (i, cutoff) in enumerate(ρ)
    kept = P.edges.nzval .≥ cutoff
    U[i] = sum(kept) / length(kept)
    L[i] = sum(P.edges.nzval[kept]) / links(P)
    S[i] = richness(simplify(P ≥ cutoff)) / richness(P)
end

#%% Central difference as a proxy for second derivative
∂U = zeros(Float64, length(U))
∂L = zeros(Float64, length(U))
for i in 2:(length(U) - 1)
    ∂U[i] = U[i + 1] + U[i - 1] - 2U[i]
    ∂L[i] = L[i + 1] + L[i - 1] - 2L[i]
end

#%% Euclidean distance for U
dU = zeros(Float64, length(U))
for i in eachindex(U)
    dU[i] = sqrt(U[i] * U[i] + ρ[i] * ρ[i])
end

#%% Plot the results
plot(ρ, U; dpi=600, size=(500, 500), lab="Non-zero")
plot!(ρ, L; lab="Expected")
xaxis!("Cutoff", (0, 1))
yaxis!("Proportion of links left", (0, 1))
vline!([ρ[findlast(S .== 1)]]; c=:grey, ls=:dash, lab="")

savefig("figures/cutoff-interactions.png")

plot(ρ, 1.0 .- S; dpi=600, size=(500, 500), lab="Disconnected species", legend=:topleft)
l = L .* links(P)
s = S .* richness(P)
plot!(ρ, l ./ (s .* s); lab="Connectance")
xaxis!("Cutoff", (0, 1))
yaxis!("  ", (0, 0.5))
vline!([ρ[findlast(S .== 1)]]; c=:grey, ls=:dash, lab="")

savefig("figures/cutoff-connectance.png")

#%% Get the threshold
thrind = findlast(S .== 1)
@info "Optimal cutoff based on remaining species: $(ρ[thrind])"
@info "Optimal cutoff based on central differences: $(ρ[last(findmax(∂U))])"

#%% Cleaned-up network
K = copy(P)
K.edges[P.edges .< ρ[thrind]] .= 0.0
dropzeros!(K.edges)

int = DataFrame(; from=String[], to=String[], score=Float64[])
for i in interactions(K)
    push!(int, (i.from, i.to, i.probability))
end
sort!(int, [:score, :from, :to]; rev=[true, false, false])

CSV.write("artifacts/canadian_thresholded.csv", int)

#%% Write the functional classification of species
rls = DataFrame(;
    sp=String[],
    gen=Float64[],
    gen_var=Float64[],
    vul=Float64[],
    vul_var=Float64[],
    role=Symbol[],
)

din = degree(K; dims=2)
dinv = degree_var(K; dims=2)
dout = degree(K; dims=1)
doutv = degree_var(K; dims=1)

for s in species(K)
    rl = :int
    iszero(din[s]) && (rl = :top)
    iszero(dout[s]) && (rl = :prd)
    push!(rls, (s, dout[s], doutv[s], din[s], dinv[s], rl))
end

sort!(rls, :gen; rev=true)
CSV.write("artifacts/species_roles.csv", rls)

plot(log1p.(sort(rls.gen; rev=true)); lab="Out-degree", size=(500, 500), dpi=600)
plot!(log1p.(sort(rls.vul; rev=true)); lab="In-degree")
xaxis!("Rank", (0, 180))
yaxis!("log(degree + 1)", (0, 5))

savefig("figures/final-degree.png")
````

