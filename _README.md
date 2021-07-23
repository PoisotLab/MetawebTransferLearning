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

| Step                                       |               Script                |                      Vignette                       |                                Notebook                                |
| ------------------------------------------ | :---------------------------------: | :-------------------------------------------------: | :--------------------------------------------------------------------: |
| Correct names from the European metaweb    |   [:computer:](00_match_names.jl)   |   [:page_facing_up:](vignettes/00_match_names.md)   |   [:notebook_with_decorative_cover:](vignettes/00_match_names.ipynb)   |
| Correct names from the reference phylogeny |   [:computer:](01_clean_tree.jl)    |   [:page_facing_up:](vignettes/01_clean_tree.md)    |   [:notebook_with_decorative_cover:](vignettes/01_clean_tree.ipynb)    |
| Harmonize IUCN names                       |   [:computer:](02_clean_iucn.jl)    |   [:page_facing_up:](vignettes/02_clean_iucn.md)    |   [:notebook_with_decorative_cover:](vignettes/02_clean_iucn.ipynb)    |
| Reconcile all mammal names                 | [:computer:](03_extract_mammals.jl) | [:page_facing_up:](vignettes/03_extract_mammals.md) | [:notebook_with_decorative_cover:](vignettes/03_extract_mammals.ipynb) |
| Make the predictions                       |   [:computer:](04_prediction.jl)    |   [:page_facing_up:](vignettes/04_prediction.md)    |   [:notebook_with_decorative_cover:](vignettes/04_prediction.ipynb)    |
| Compare with Newfoundland data             |  [:computer:](05_StrongLeroux.jl)   |  [:page_facing_up:](vignettes/05_StrongLeroux.md)   |  [:notebook_with_decorative_cover:](vignettes/05_StrongLeroux.ipynb)   |
| Compare with GLOBI data                    |    [:computer:](06_inflation.jl)    |    [:page_facing_up:](vignettes/06_inflation.md)    |    [:notebook_with_decorative_cover:](vignettes/06_inflation.ipynb)    |
| Produce the Canadian metaweb               |  [:computer:](07_thresholding.jl)   |  [:page_facing_up:](vignettes/07_thresholding.md)   |  [:notebook_with_decorative_cover:](vignettes/07_thresholding.ipynb)   |