# Food web reconstruction using phylogenetic inference of latent variables on ecological interactions through network embedding

The extreme difficulty in documenting species interactions 
and being able to build networks poses a considerable 
challenge to understanding the emergent properties of ecological communities
[@Jordano2016SamNet; @Jordano2016ChaEco].
With the added challenge that real world datasets are sparse 
and geographically biased [@Poisot2021GloKno]. Here, using 
interaction data from the European metaweb, a species list 
and the mammalian phylogenetic tree along with network embedding 
we are able reconstruct a mammalian Canadian food web. This 
knowledge transfer approach is computationally inexpensive 
and can help us make predictions when local data are missing 
by leveraging the fact that traits are conserved at the 
phylogenetic level and determine the nature of an interaction 
between species.
**Do we throw a sentence in here about validation?**
**I feel like we need to bring up the Stock 2021 paper somehow**

The probability of an interaction occurring between two species depends 
on a set of conditions being met, i) they need to co-occur
**and** ii) species are 'compatible' based on behavioural 
or physical traits [@Jordano2016SamNet]. This trait-driven selection 
results in interactions being conserved at the evolutionary level and, 
by extension, the phylogeny, which can be used to infer interactions 
between species [@Davies2021EcoRed; @Elmasri2020HieBay; @Gomez2010EcoInt]. 
Using a known interaction network we can *learn* which latent traits are determining 
what interactions between species. This knowledge is then *transferred* 
to a different species pool using phylogenetic relatedness 
infer their latent traits, this is then used to predict the 
interaction network (Fig 1).

> Figure 1: A: We have a workflow overview which we can potentially 
> subdivide/section to follow the narrative more closely. B: The Euro 
> network C: The CA network???

*Retcon I think we should just have the conceptual figure here and have results relegated to figure 2?*

*Embedding the knowledge:*
Using a truncated singular value decomposition (*t-SVD*) of the European 
metaweb, we are able to extract the left and right
subspaces - representing the latent traits of predator and prey species respectively. 
We then map these latent traits to the phylogenetic 
tree *i.e* matching latent traits to the phylogeny. For species shared between 
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
comprehensive breakdown of the methods). This yields a Canadian mammalian
food web, consisting of 223 species and approximately 6000 interactions
(refer either to fig 1 or 2?).

> Figure 2: A: Possibly ROC-AUC, B: Body mass, C: omnivory (unless we have this
> as panels in fig 1)

*The validity of the knowledge:*
At face value this technique is shown to perform excellently when 
compared to a known network for Newfoundland [@Strong2014ImpNon; fig 2]. 
However, if we were to extend this workflow to a region where we have no 
data it does raise the question as to how we can validate the predicted 
network. Turning to ecological concepts and relationships as a means 
of validating model outputs *e.g.* using the estimated body size 
relationships between predator and prey [@Brose2006ConRes] or looking 
at network structure/motifs (**ref**), as well as using other 
predictions of network structure *e.g.* number of links 
[@MacDonald2020RevLina] may be a promising alternative (ref fig2). 
A lack of an obvious approach to validating model outputs in the absence 
of empirical data suggests a need to place a greater focus on the 
mechanistic drivers of interactions between species. *I too am judging this statement*

*Future directions:*
Although we present a very compelling use-case for network embedding to 
make network predictions there is the caveat that our two species pools 
share a fair amount of phylogenetic overlap and are restricted to the same 
class and the question remains if this methodology can 'scale-up' to more 
disparate communities or more branches of the tree of life. There is also 
a need to validate alternative approaches to model validation if we are 
to maximise the potential of using network embedding to make predictions 
for data deficient localities. Caveats aside, network embedding has the 
potential extend its use-case beyond predicting across space and could 
be used to transfer knowledge across time *i.e.* reconstruct historic 
networks or as a means of validating empirical networks by detecting 
interactions that may have been unobserved in the field.
