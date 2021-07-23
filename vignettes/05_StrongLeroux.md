# Step 6 - inflating the predictions with the Newfoundland data

````julia
using DelimitedFiles
using DataFrames
using CSV: CSV
using GBIF
using NCBITaxonomy: NCBITaxonomy
using EcologicalNetworks
````

We start by reading the Newfoundland food web, and check the names that are
mammals. The actual code to read the network looks exactly like the one to
read the European metaweb.

````julia
sl_raw = readdlm("data/NLfoodweb.csv", ',')

sl_sp = replace.(sl_raw[1, 2:end], "." => " ")
sl_A = Bool.(sl_raw[2:end, 2:end])
````

Because the original data use a mix of scientific and vernacular names, we are
going to rely on `NCBITaxonomy.jl. synonym matching abilities to first get the
taxonomic names, and then pass those to GBIF.  Please do keep in mind that
unless the `NCBITAXONOMY_PATH` environmental variable is set, the raw taxonomy
dump will be stored in the project folder (and this is a rather big file).

````julia
scinames = Dict{String,String}()
````

Note that we do *not* restrict the name matching to only mammals, as there are
non-mammal species in the Newfoundland metaweb.

````julia
for s in sl_sp
    try
        t = NCBITaxonomy.taxon(s; strict=false)
        scinames[s] = t.name
    catch
        @info "Newfoundland taxon $(s) unmatched on NCBI"
        continue
    end
end
````

The next step is to get the names from NCBI, and match them to the GBIF
backbone. We ended up relying on this two-step solution because using the GBIF
name matching directly missed a handful of species, and the Newfoundland
dataset is relatively small.

This loop will go through all nodes in the Newfoundland metaweb, match them at
the species level, and only return them if they are part of the *Mammalia*
class. There may be a few info messages about unmatched taxa, which are nodes
from the original data that are at a higher rank than species.

````julia
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
````

With the two dictionaries, we can get the positions of species from the
Newfoundland metaweb that are mammals:

````julia
idxmatch = findall(x -> x in keys(valnames), sl_sp)
````

And we can now assemble the network:

````julia
spnames = [valnames[s] for s in sl_sp[idxmatch]]
A = sl_A[idxmatch, idxmatch]'
NL = UnipartiteNetwork(A, spnames)
````

We finally save the network as a CSV - note that we do not add interactions
here, as this will be done as part of the thresholding step, which is the very
last in the pipeline.

````julia
df = DataFrame(; from=String[], to=String[])
for i in interactions(NL)
    push!(df, (i.from, i.to))
end
CSV.write("artifacts/newfoundland.csv", df)
````

