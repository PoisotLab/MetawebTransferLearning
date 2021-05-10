# Phylogenetic inference of ecological interactions through network embedding

## Introduction

Why we should care about ecological networks yet information 
is scarce and collecting it is hard so we need reliable ways 
to predict them (Poisot 2021 100% fits in here somewhere) 
networks are difficult to adequately sample in nature 
Banašek‐Richter 2004 Chacoff 2012, Gibson 2011, Jordano, 2016 
--- both (both is good *i.e.* Jordano a and b) 
*why do we care about predicting networks*

Other methods such as the fourth corner model-based approach (Brown_Warton_Andrew_Binns_Cassis_Gibb_2014)
could be used to reproduce the network, by considering phylogeny. 
However, "Model selection methods tend to perform best when given a relatively 
small number of variables (and hence a small number of candidate models) 
to select from Burnham & Anderson (2002).It is especially important to try to keep the number 
of explanatory variables small in fourth‐corner models because of the use of interaction terms – 
a model with ten environmental terms and ten trait terms would then have 100 interaction terms, 
which is usually too many to reliably identify the important interactions when using just a single 
data set. Hence, one should carefully consider prior to analysis which variables to include and 
attempt to compile a ‘shortlist’. Particularly, useful for shortlisting variables is prior knowledge 
from the literature, and diagnostic tools for collinearity (Zuur, Ieno & Elphick 2010)." (from Brown&al.2014)

The tools we decided to use are more flexible and 'only need limited data'(-> define). 
The network is predicted using phylogenetic inference of ecological interactions through truncated-SVD, 
RDPG and network embedding (refer to Stock 2021 in here, Tanya).
SVD reduces information to important ones, while RDPG (uses latent space to..) and network-embedding (to learn...).

![Much art.](figures/conceptual.png){#fig:conceptual}

Goal is to use phylogeny to predict interactions

## Methods

### Dataset preparation 

[Data cleaning toss in supplementary; only mention sources]

Brief overview of the data *i.e* the European metaweb [@Maiorano2020DatTet] 
*this is the actual dataset* as well 
as the Upham tree [@Upham2019InfMam]. Using GBIF to get a 
mammals species list as 
well as for taxonomic harmonisation/resolution and how we cleaned 
it. We reconciled the names in the European metaweb to the GBIF 
taxonomic backbone, resulting in 260 species with an estimated 
2342 trophic interactions. We then reconciled the entire Upham 
et al. phylogenetic tree to the GBIF taxonomic backbone, ensuring 
that we could merge the data. 
We downloaded a list of mammals having been observed in Canada at 
least twice according to GBIF - this list includes fossils which 
are removed when matched against the tree. We filtered out all 
marine mammals.

### Prediction

#### Network Decomposition

[I think we can find a way to cut straight to t-SVD]

*What it is*

Singular Value Decomposition [SVD, @Forsythe1967ComSol; @Golub1971SinVal] 
is the factorisation of a matrix $\mathbf{A}$
(where $\mathbf{A}_{m,n} \in\mathbb{B}$, in the context of adjacency matrices) into the form:

$$ \mathbf{U}\cdot\mathbf{\Sigma}\cdot\mathbf{V}^T $$

Where $\sigma_{i} = \Sigma{ii}$, which contains the singular values of 
$\mathbf{A}$. When the values of $\mathbf{\sigma}$ are arranged in 
descending order, the singular values ($\mathbf{\Sigma}$) are
unique. 
It is also possible to do this in the truncated form where $\mathbf{A}$ becomes 
$\mathbf{U}\cdot\mathbf{\Sigma}\cdot\mathbf{V}^*$
where $\mathbf{\Sigma}$ is a square diagonal of size $k$ (*i.e* the rank of 
$\mathbf{A}$) and has only non-zero 
$\sigma$ values. $\mathbf{U}$ is now a $m \times k$ matrix and $\mathbf{V}$ 
is a $n \times k$ matrix.
The decomposition into left/right subspaces is unique, represents 
latent traits for resp. outgoing and incoming edges, and can be
 used to predict interactions. Predator traits on the one side 
 and prey (vulnerability traits) on the other

 *Selecting the rank:* mention the elbow method here ? Used to determine the number of clusters.

*Network embedding:* We extracted the left and right subspaces 
of a rank-5 *t-SVD* of the European mammals metaweb. When 
multiplied, these two matrices give an approximation of the 
network; we determined a threshold indicating that an interactions 
is likely by finding the value maximizing Youden's J statistic. 
At rank 5, this was approx. 0.11, giving J > 0.95.

*Latent traits inference:* For all species in the Canadian pool 
that had also been reported in the European metaweb, we kept their 
latent traits as is. For other species, we reconstructed their 
latent traits by averaging the values of those of their 3 closest 
neighbors (based on the cophenetic matrix of the entire tree). 
This ensure that shared species are at the same position in both 
latent subspaces.

*Prediction:* The prediction of interactions was done through a 
random dot product graph (RDPG), by multiplication of the inferred left 
and right subspaces and thresholding at the value estimated during 
step 2. This yields a network with 223 species of mammals and approx. 
6000 interactions.

### Validation

*Validation:* Network structure (omnivory), body mass ratio 
(biological relationships), and Newfoundland tree [@Strong2014ImpNon]

In the absence of empirical networks, we used ecological concepts as measures of trust to assess 
the efficacy/accuracy/performance of the predicted network.
The predator-prey bodymass/(size?) ratio influences the predator-prey interactions 
in food-web studies (Nakazawa_Ushio_Kondoh_2011). Predation matrices are used 
to verify the body mass ratio.
Visualizing the network helps to identify omnivory (trophic levels) in the network structure.
A sample of the model's prediction is compared to an empirical network 
(Newfoundland tree [@Strong2014ImpNon], beavers?) to evaluate the model's performance. 


## Results and discussion

> Figure 1 is a multi-paneled output of results: ROC, matrix??

Hey look lots of similarities between Europe and Canada and 
here are some metrics

> Figure 2 is a multi-paneled output of validation: body mass 
> ratio, omnivory/trophic comparison and Newfoundland metaweb

Challenges of validation - even empirical webs are imperfect 
Need for a suite of validation techniques?


## Concluding notes

This method is the bestest - use it

## References

## Supplementary stuff