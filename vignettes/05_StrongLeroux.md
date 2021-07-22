````julia
#%% Dependencies
using DelimitedFiles
using DataFrames
using CSV: CSV
using GBIF
using NCBITaxonomy: NCBITaxonomy
using EcologicalNetworks

#%% Read the Strong & Leroux data
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

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

