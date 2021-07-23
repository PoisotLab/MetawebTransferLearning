# Phylogenetic inference of ecological interactions through network embedding

## Executive summary

1. network decomposition through a truncated-SVD capture the evolutionary signal of species interactions
2. the decomposition into left/right subspaces is unique, represents latent traits for resp. outgoing and incoming edges, and can be used to predict interactions
3. we can infer latent traits for unobserved species based on phylogenetic proximity with observed species
4. we use this information to transfer information from the trophic interaction of European mammals to Canadian mammals for which we have no *a priori* data

## About this README

This README is generated using `Literate.jl` and contains *every line* or *every
file* used to produce the entire results. This is, essentially, the director's
commentary version of the analysis - there are discussions of the purpose, but
also discussion of technical choices.

## Reproducing the code

| Step                           | Script | Annotated MD | Notebook |
| ------------------------------ | ------ | ------------ | -------- |
| Match names                    |        |              |          |
| Cleanup phylogeny              |
| Cleanup IUCN                   |
| Reconcile mammal names         |
| Make the predictions           |
| Compare with Newfoundland data |
| Compare with GLOBI data        |
| Produce the Canadian metaweb   |