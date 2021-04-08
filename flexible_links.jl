import CSV
using DataFrames
using Distributions
using StatsBase
using StatsPlots
using Turing

# Load data from Mangal (ecological interactions database)
# S: number of species (nodes)
# L: number of interactions (edges)
# P: number of interactions of predation
# H: number of interactions of herbivory
ls = CSV.read(joinpath("data","ls.csv"), DataFrame)

# Filter for food webs
ls = ls[ls.P .+ ls.H .> 0,:]

# Build model
# S: number of species
# R: number of flexible links (i.e. number of links above the minimum)
@model FL(S,R) = begin
  N = length(S)
  # Number of trials
  F = S.*S .- (S.-1)
  # Parameters (priors)
  ϕ ~ Normal(3.0, 0.5)
  μ ~ Beta(3.0, 7.0)
  for i in 1:N
    R[i] ~ BetaBinomial(F[i], μ*exp(ϕ), (1-μ)*exp(ϕ))
  end
  return μ, ϕ
end

# Train model
S = vec(ls[:,:S])
L = vec(ls[:,:L])
R = L .- (S.-1)
chain = sample(FL(S,R), HMC(0.01,10), 3000)

# Visualize posterior distribution
plot(chain[200:end,:,:])


# Number of species, links and flexible links in the Canadian metaweb
s = 224
l = 7098
r = l-s+1

# Simulate n values of L using the entire posterior distribution
n = 10000

function prediction(chain, S)
  p = get_params(chain[200:end,:,:])
  i = rand(1:length(p.μ))
  μ, ϕ = p.μ[i], p.ϕ[i]
  return rand(BetaBinomial(S^2-(S-1), μ*exp(ϕ), (1-μ)*exp(ϕ))) + (S-1)
end

L_hat = [prediction(chain, s) for i in 1:n]

# Proportion of simulated links smaller than the number of links in the metaweb
quantile_post = sum(l.>L_hat) / n

# Parameter values (point estimates)
μ = median(get_params(chain).μ)
ϕ = median(get_params(chain).ϕ)

# MAP distribution of flexible links
bb = BetaBinomial(s^2-(s-1), μ*exp(ϕ), (1-μ)*exp(ϕ))

# Quantile of l (using MAP distribution)
quantile_map = cdf(bb, r)
quantile_map = round(quantile_map; digits=4)

# z-score of l (using normal approximation of MAP distribution)
means = (s^2-s+1)*μ + s-1
sds = sqrt((s^2-s+1)*μ*(1-μ)*(1 + s*(s-1)/(1+ϕ)))

z_score = (l-means)/sds
z_score = round(z_score; digits=4)


# Does the number of links in the Canadian metaweb make sense?
println("Number of links in the Canadian metaweb in light of the flexible links model

    Number of species: $s
    Number of links: $l
    Quantile (posterior distribution): $quantile_post
    Quantile (MAP): $quantile_map
    z-score: $z_score")
