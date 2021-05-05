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

NLmetaweb.Interaction = fill(true, size(NLmetaweb,1))
NLmetaweb = leftjoin(A, NLmetaweb, on = [:Predator => :Predator, :Prey => :Prey])
replace!(NLmetaweb.Interaction, missing => false)

PredictedMetaweb.Interaction = fill(true, size(PredictedMetaweb,1))
PredictedMetaweb = leftjoin(A, PredictedMetaweb, on = [:Predator => :Predator, :Prey => :Prey])
replace!(PredictedMetaweb.Interaction, missing => false)

# Confusion matrix
tp = sum(PredictedMetaweb.Interaction .& NLmetaweb.Interaction)
fp = sum(PredictedMetaweb.Interaction .& .!NLmetaweb.Interaction)
fn = sum(.!PredictedMetaweb.Interaction .& NLmetaweb.Interaction)
tn = sum(.!PredictedMetaweb.Interaction .& .!NLmetaweb.Interaction)

sensitivity = tp / (tp + fn) # true positive rate
specificity = tn / (tn + fp) # true negative rate
tss = specificity + sensitivity - 1 # True skill statistic
accuracy = (tp + tn) / size(NLmetaweb, 1) # accuracy
