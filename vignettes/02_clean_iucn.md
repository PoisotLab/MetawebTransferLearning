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

## Extinct species removal

Two species in the IUCN dataset are considered to be extinct, and we therefore
remove them as well.

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

