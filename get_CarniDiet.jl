using CSV
using DataFrames
using HTTP

# For now, let's keep the carniDiet db (>30 mb) out of the repo
carniDiet = DataFrame(CSV.File(HTTP.get("https://raw.githubusercontent.com/osmiddleton/CarniDIET-Database/master/Version%201.0/CarniDIET%201.0.csv").body))

# Replace "_" by " " in species names...
carniDiet.scientificNameCarni = replace.(carniDiet.scientificNameCarni, "_" => " ")
carniDiet.scientificNamePrey = replace.(carniDiet.scientificNamePrey, "_" => " ")

canada = DataFrame(CSV.File(joinpath("artifacts", "iucn_gbif_names.csv")))

# Option 1: Use interactions between two canadian mammals, observation made in Canada
carniDiet_canada = filter(row -> row[:country] == "Canada", carniDiet) # 1025 interactions
carniDiet_canada = filter(row -> row[:scientificNameCarni] in canada.gbifname, carniDiet_canada) # 956 interactions (name mismatch ?)
carniDiet_canada = filter(row -> row[:scientificNamePrey] in canada.gbifname, carniDiet_canada) # 370 interactions


# Option 2: Use interactions between two canadian mammals, observation made in anywhere
carniDiet_canada2 = filter(row -> row[:scientificNameCarni] in canada.gbifname, carniDiet) # 12 236 interactions (name mismatch ?)
carniDiet_canada2 = filter(row -> row[:scientificNamePrey] in canada.gbifname, carniDiet) # 2 443 interactions
