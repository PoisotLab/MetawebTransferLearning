# Phylogenetic inference of ecological interactions through network embedding

The extreme difficulty in empirically documenting species interactions 
and being able to build interaction networks poses a considerable 
challenge within the context of studying and understanding ecological 
networks. This means that real world datasets are sparse and geographically 
biased [@Poisot2021GloKno]. However, there is feasibility in taking 
the knowledge of species interactions for a given location, transferring 
it to a different species pool or location and making a prediction of what 
the interaction network may look like. In other words can we 
transfer the knowledge contained in one network by using phylogenetic 
inference of the ecological interactions through using network embedding 
to make accurate predictions of interaction networks for a species pool 
or location for which we have no *a priori* knowledge of interactions?

The probability of an interaction occurring between two species depends 
on a set of conditions being met, including that they need to co-occur to 
some extent **and** that the two species are 'compatible' in the sense that 
an interaction is actually possible based on a physical or behavioural 
trait *e.g.* the relationship between predator and prey gape size 
(**better e.g here**) [@Jordano2016SamNet]. This trait-driven selection 
results in interactions being conserved at the evolutionary level and this 
signal will be captured within a phylogeny. Thus using the phylogenetic 
relatedness for a given community can inform as to how they may
interact with each other 
[@Davies2021EcoRed; @Elmasri2020HieBay; @Gomez2010EcoInt]. Given these 
relationships it is thus possible to predict what an interaction network 
may look like using two (relatively) easily accessible data sources - 
namely the community composition *i.e.* species occurrence records as 
well as some measure of how they are related *i.e.* a taxonomic 
backbone/measure of phylogenetic relatedness for the given species pool. 
In addition a known species interaction is needed to allow for the model 
to learn how communities and traits determine interactions - importantly 
the network need not be a perfect 1:1 match with regards to species overlap.
(**still some kinks in this last section**).

In this manuscript we present a proof-of-concept using the mammalian 
component of the European 
metaweb [@Maiorano2020TetEu] to predict the Canadian mammalian metaweb 
using only species occurence records for the region as well as a resolved 
mammalian phylogeny [@Upham2019InfMam] as additional data sources (Fig 1), 
for which both data types/sources are easily accessible.

> Figure 1: A: We have a workflow overview which we can potentially 
> subdivide/section to follow the narrative more closely. B: The Euro 
> network C: The CA netwrok???

*Embedding the knowledge:*
Using a truncated singular value 
decomposition (*t-SVD*) of the European metaweb we are able to extract the 
left and right subspaces, these subspaces are representative of the latent 
traits of prey and predator species respectively **CHECK**. This provides 
us with information as to how species' latent traits determine their 
interactions. We then map these latent traits to the mammalian phylogenetic 
tree *i.e* matching traits to the phylogeny. For species shared between 
Europe and Canada it is possible to directly infer their latent traits 
for novel species, we reconstruct their latent traits by averaging the 
values of those of their 3 closest neighbours (based on the cophenetic 
matrix of the entire tree). This ensures that shared species are at the 
same position in both latent subspaces.
**Potentially controversial - went present tense here**

*Transferring the knowledge:*
Through the process of inferring the latent traits of the Canadian mammals 
we are also recreating/inferring the left and right subspaces for the 
Canadian community. Leveraging the predictive potential from *t-SVD* 
subspaces as well as a random dot product graph (RDPG), we multiply the 
inferred left and right Canadian subspaces at a given threshold to predict 
the mammalian Canadian metaweb (see the extended methods for a more 
comprehensive breakdown of the methods). 

*The outcomes:*
**I do think this can be a one-liner we throw in with the previous/next section?**

> Figure 2: A: Possibly ROC, B: Body mass, C: omnivory (unless we have this
> as panels in fig 1)

*The validity of the knowledge:*
At face value this technique is shown to preform excellently (fig 2) 
it does raise the question of validation techniques with regards to 
biological plausibility. In the absence of empirical data to use as a 
test/validation dataset we need to turn to ecological concepts as a means 
of validating model outputs. For example using the estimated body size 
relationships between predator and prey [@Brose2006ConRes] or looking 
as network structure/motifs [**more on this**] as well as using other 
predictions of network structure *e.g.* number of links 
[@MacDonald2020RevLina] and compare those to the predicted network 
could be used as an 
indication of the ecological validity for a predicted network (ref fig2).
As a caveat even if we were to have an 
empirical dataset we need to query the completeness of that dataset 
given that observations may be overlooked in the field *e.g.*
using a similar pairwise learning approach @Stock2021PaiLea made 
predictions for interactions that had previously be unobserved in their 
network. This does also raise the possibility of using this transfer 
learning *within* a dataset as a means of evaluating its completeness 
as well as filling in potential gaps. (**I think this is valid???**). 

**This section can be set up can be better**

*Future directions:*
In our use-case we are predicting a new network using a network for which 
there is a fair amount of phylogenetic overlap between species to very 
successful results. However, the question remains as to the efficacy of 
this methodology when the overlap decreases - do we see a decline in 
predictive performance with a decline in overlap? Our species pool was 
constrained to only mammals - again what happens when we start to 
incorporate a greater diversity taxa?