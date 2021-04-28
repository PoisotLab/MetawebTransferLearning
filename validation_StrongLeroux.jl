using DataFrames
using CSV

sp_list = DataFrame(CSV.File(joinpath("data", "StrongLeroux","sp_transfer.csv")))
NLmetaweb = DataFrame(CSV.File(joinpath("data", "StrongLeroux","StrongLeroux_raw.csv")))

# Change predator names to scientific name
rename!(NLmetaweb, Dict("Focal species" => "FocalSpecies"))
NLmetaweb = leftjoin(NLmetaweb, sp_list, on = :FocalSpecies => :CommonName)
rename!(NLmetaweb, Dict("ScientificName" => "Predator"))

# Change prey names to scientific name
rename!(NLmetaweb, Dict("Food item" => "FoodItem"))
NLmetaweb = leftjoin(NLmetaweb, sp_list, on = :FoodItem => :CommonName)
rename!(NLmetaweb, Dict("ScientificName" => "Prey"))

# Clean the resulting metaweb
select!(NLmetaweb, :Predator, :Prey)
dropmissing!(NLmetaweb)
unique!(NLmetaweb) # 51 interactions

# Load Prediction
PredictedMetaweb = DataFrame(CSV.File(joinpath("artifacts/", "canadianmetaweb.csv"), header = 0))
rename!(PredictedMetaweb, ["Column1" => "Predator", "Column2" => "Prey"])
filter!(row -> row.Predator in sp_list.ScientificName, PredictedMetaweb)
filter!(row -> row.Prey in sp_list.ScientificName, PredictedMetaweb) #67 interactions

A = DataFrame()
A.Predator = repeat(sp_list.ScientificName, outer = 29)
A.Prey = repeat(sp_list.ScientificName, inner = 29)

NLmetaweb.Interaction = ones(size(NLmetaweb, 1))
NLmetaweb = leftjoin(A, NLmetaweb, on = [:Predator => :Predator, :Prey => :Prey])

PredictedMetaweb = leftjoin(A, PredictedMetaweb, on = [:Predator => :Predator, :Prey => :Prey])
# Confusion matrix
a = 
