# Food web reconstruction using phylogenetic inference of latent variables on ecological interactions through network embedding

The extreme difficulty in documenting species interactions 
and being able to build networks poses a considerable 
challenge to understanding the emergent properties of ecological communities
[@Jordano2016SamNet; @Jordano2016ChaEco].
Real world datasets are sparse and geographically 
biased [@Poisot2021GloKno]. Here we reconstruct the 
mammalian Canadian food web using interaction data from the 
European metaweb, network embedding, a species list 
and the mammalian phylogenetic tree. Thereby highlighting 
this as both a data and computationally inexpensive methodology 
for predicting interaction networks for which there is no data.
**Do we throw a sentence in here about validation?**

The probability of an interaction occurring between two species depends 
on a set of conditions being met, i) they need to co-occur
**and** ii) species are 'compatible' based on behavioural 
or physical traits [@Jordano2016SamNet]. This trait-driven selection 
results in interactions being conserved at the evolutionary level and, 
by extension, the phylogeny, which can be used to infer interactions 
between species [@Davies2021EcoRed; @Elmasri2020HieBay; @Gomez2010EcoInt]. 
Using a known interaction network we can *learn* which traits are determining 
what interactions between those species. This knowledge is then *transferred* 
to a different species pool using their phylogenetic relationship (relatedness) 
with the original species to infer traits, and by extension, the likelihood of 
an interaction occurring (Fig 1). 
**I think we need to mention that we are working with latent traits somehow...**

> Figure 1: A: We have a workflow overview which we can potentially 
> subdivide/section to follow the narrative more closely. B: The Euro 
> network C: The CA network???

*Retcon I think we should just have the conceptual figure here and have results relegated to figure 2?*

*Embedding the knowledge:*
Using a truncated singular value decomposition (*t-SVD*) of the European 
metaweb to reduce dimensionality, we are able to extract the left and right
subspaces which are the latent traits of predator and prey species respectively. 
We then map these latent traits to the mammalian phylogenetic 
tree *i.e* matching traits to the phylogeny. For species shared between 
Europe and Canada it is possible to directly infer their latent traits.
For novel species, we reconstruct their latent traits by averaging the 
values of those of their three closest neighbours (based on the cophenetic 
matrix of the entire tree). This ensures that shared species are at the 
same position in both latent subspaces.

*Transferring the knowledge:*
Through the process of inferring the latent traits of the Canadian mammals 
we are also reconstructing the left and right subspaces for the 
Canadian community. Leveraging the predictive potential from *t-SVD* 
subspaces as well as a random dot product graph (RDPG), we multiply the 
inferred left and right Canadian subspaces at a given threshold to predict 
the mammalian Canadian metaweb (see the extended methods for a more 
comprehensive breakdown of the methods). This yields the Canadian mammalian
food web, consisting of 223 species and approximately 6000 interactions
(refer either to fig 1 or 2?).

> Figure 2: A: Possibly ROC-AUC, B: Body mass, C: omnivory (unless we have this
> as panels in fig 1)

*The validity of the knowledge:*
At face value this technique is shown to perform excellently (fig 2) 
however it does raise the question as to how we can validate the predicted 
network if we were to extend this workflow to a region where we have no 
data. We can turn to ecological concepts as a means 
of validating model outputs *e.g.* using the estimated body size 
relationships between predator and prey [@Brose2006ConRes] or looking 
at network structure/motifs (**ref**), as well as using other 
predictions of network structure *e.g.* number of links 
[@MacDonald2020RevLina] and comparing those to the predicted network 
could be used as an indication of the ecological validity for a 
predicted network (ref fig2).
As a caveat even if we were to have an 
empirical dataset we need to query the completeness of that dataset 
given that observations may be overlooked in the field *e.g.*
using a similar pairwise learning approach @Stock2021PaiLea made 
predictions for interactions that had previously been unobserved in their 
network. This does also raise the possibility of using this transfer 
learning *within* a dataset as a means of evaluating its completeness 
as well as filling in potential gaps. 
**Depending on our overlap with the Newfoundland network we can weave those results in here as well**

*Future directions:*
In our use-case we are predicting a new network using a network for which 
there is a fair amount of phylogenetic overlap between species to very 
successful results. However, the question remains as to the efficacy of 
this methodology when the overlap decreases - do we see a decline in 
predictive performance with a decline in overlap? Our species pool was 
constrained to only mammals - again what happens when we start to 
incorporate a greater diversity taxa?

**The idea/feasibility of reconstruction i.e. historic networks**

**Lack of interaction strength - we only work with binary - should be fine since SVD takes $\mathbb{R}$ numbers along with $\mathbb{B}$**

**Differences in diet in different region - but so different that there is actually no overlap... different populations of the same species i.e. phylogeny might not catch this**
