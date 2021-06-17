using HTTP
using JSON
using DataFrames
using CSV: CSV

# API info
_globi_api = "https://api.globalbioticinteractions.org/taxon"
relevant_types = ["eats", "preysOn", "kills"]

canmet = DataFrame(CSV.File("artifacts/canadian_corrected.csv"))

allsp = unique(vcat(canmet.from, canmet.to))

url = "$(_globi_api)/$(rand(allsp))/eats/"

diet = DataFrame(; from=String[], to=String[], type=String[])

r = HTTP.request("GET", url; verbose=3)
globidiet = JSON.parse(String(r.body))

