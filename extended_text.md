



**Motivation** 
1. This is what we want. 

The extreme difficulty in empirically capturing species interactions 
and being able to build interaction networks poses a considerable 
challenge within the context of studying and understanding ecological 
networks. This means that real world datasets are sparse and geographically 
biased [@Poisot2021GloKno]. However, there is feasibility in taking 
the knowledge of species interactions for a given location and transferring 
this to a different species pool or location.




1. This is what is missing. 
2. But this is what we have. 
3. We can take what we have to a new place (we can bring it with 
us - "transfer knowledge across space") . This is what we do
5. This is what we got
6. How to make sure that what we got makes sense.

How to bring knowledge from one place to another

> Figure 1: A = conceptual figure; B = Euro web; C = CA web

![Much art.](figures/conceptual.png){#fig:conceptual}

*Next three points will be very short maybe even 1 paragraph for all*

**Phylogeny** We have species occurrence records *i.e.* so we know who is 
there. Even though we don't have the same species we have phylogenetic 
overlap *i.e.* so we know how they are connected/related.

**SVD** It reduces dimensionality (truncates) but more importantly 
it is very good for making predictions.

**RDPG** This is the transferring of knowledge part *i.e* predicting. This 
is how we make predictions in a new 'space'
"RDPG is a latent position generative model, in which the probability of an edge 
existing between pairs of vertices is determined by the dot product of the associated 
latent position vectors."

**Validation** Hey we found some really cool result but how de we validate them? 
We can look at structure (omnivory), biological relationships (body mass) and 
empirical data [@Strong2014ImpNon]

> Figure 2: A = comparison trophic levels, omnivory; B = body mass?

*Future considerations* Degree of phylogenetic overlap - how much do we need? 
Possibly woven into validation section


## 'Throw away' text
## Introduction

Why we should care about ecological networks yet information 
is scarce and collecting it is hard so we need reliable ways 
to predict them (Poisot 2021 100% fits in here somewhere) 
networks are difficult to adequately sample in nature 
Banašek‐Richter 2004 Chacoff 2012, Gibson 2011, Jordano, 2016 
--- both (both is good *i.e.* Jordano a and b) 
*why do we care about predicting networks*



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

## References



