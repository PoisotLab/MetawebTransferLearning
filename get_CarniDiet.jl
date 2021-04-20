using CSV
using DataFrames
using HTTP

# For now, let's keep the carniDiet db (>30 mb) out of the repo
carniDiet = DataFrame(CSV.File(HTTP.get("https://raw.githubusercontent.com/osmiddleton/CarniDIET-Database/master/Version%201.0/CarniDIET%201.0.csv").body))

#TODO: Correct sp names in CarniDiet db to match gbif names.

# Replace "_" by " " in species names...
carniDiet.scientificNameCarni = replace.(carniDiet.scientificNameCarni, "_" => " ")
carniDiet.scientificNamePrey = replace.(carniDiet.scientificNamePrey, "_" => " ")

canada = DataFrame(CSV.File(joinpath("artifacts", "iucn_gbif_names.csv")))

# Option 1: Use interactions between two canadian mammals, observation made in Canada
carniDiet_canada = filter(row -> row[:country] == "Canada", carniDiet) # 1025 interactions
carniDiet_canada = filter(row -> row[:scientificNameCarni] in canada.gbifname, carniDiet_canada) # 956 interactions
carniDiet_canada = filter(row -> row[:scientificNamePrey] in canada.gbifname, carniDiet_canada) # 370 interactions
select!(carniDiet_canada, :scientificNameCarni, :scientificNamePrey)
unique!(carniDiet_canada) # 72 unique interactions

# Option 2: Use interactions between two canadian mammals, observation made in anywhere
carniDiet_canada2 = filter(row -> row[:scientificNameCarni] in canada.gbifname, carniDiet) # 12 236 interactions
carniDiet_canada2 = filter(row -> row[:scientificNamePrey] in canada.gbifname, carniDiet_canada2) # 1 750 interactions
select!(carniDiet_canada2, :scientificNameCarni, :scientificNamePrey)
unique!(carniDiet_canada2) # 221 unique interactions

# Let's look at Europe for fun!
europe = DataFrame(CSV.File(joinpath("artifacts", "names_metaweb_tree_gbif.csv")))
carniDiet_europe = filter(row -> row[:scientificNameCarni] in europe.name, carniDiet) # 15 326 interactions 
carniDiet_europe = filter(row -> row[:scientificNamePrey] in europe.name, carniDiet_europe) # 3 900 interactions
select!(carniDiet_europe, :scientificNameCarni, :scientificNamePrey)
unique!(carniDiet_europe) # 412 unique interactions
