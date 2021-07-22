````julia
using GBIF
using CSV, DataFrames
using ProgressMeter
using Base.Threads

checklist = DataFrame(CSV.File(joinpath("data", "taxonomy.csv")))
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

checklist_cleanup_components = [
    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for
    i in 1:nthreads()
]

checklist_cleanup_components = [
    DataFrame(; code=String[], gbifname=String[], gbifid=Int64[], equal=Bool[]) for
    i in 1:nthreads()
]

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

checklist_cleanup = vcat(checklist_cleanup_components...)
CSV.write(joinpath("artifacts", "iucn_gbif_names.csv"), checklist_cleanup)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
