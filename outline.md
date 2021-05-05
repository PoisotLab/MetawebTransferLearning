# Phylogenetic inference of ecological interactions through network embedding

## Introduction

Why we should care about ecological networks yet information 
is scarce and collecting it is hard so we need reliable ways 
to predict them (Poisot 2021 100% fits in here somewhere) 
networks are difficult to adequately sample in nature 
Banašek‐Richter 2004 Chacoff 2012, Gibson 2011, Jordano, 2016 
--- both (both is good *i.e.* Jordano a and b) 
*why do we care about predicting networks* -> **Tanya**

But hey we can predict them

Hey lets lean on phylogenies because they are useful, 
available, tends to be 'global', somewhat 'static'
Maybe something on the value of using phylogenies? They are 
more readily available (and cover a greater geographic 
extent) AND phylogenies tell us a lot about the species - 
trait conservatism. Linking phylogeny, traits and 
interactions. do we discuss the value of traits in determining interactions??
if yes turn to here: Bartomeus 2016, Brousseau 2018, 
Desjardins‐Proulx 2017, Gravel 2013 data‐poor environments (Beauchesne et al., 2016)
-> **Tanya**

Here some ideas using fancy maths:
Why other methods are potentially not as great 
4th corner (abundance, site characteristics, 
species traits) - need 3 different datasets, RLQ method, Need phylogeny and 
interactions matrix. 
*Discuss why they might not work* -> **Salomé**

mention tools we are using in light of their flexibility and 
how they 'only' need 'limited' data *could go above or below*

Phylogenetic inference of ecological interactions through 
network embedding
mention tools we are using in light of their flexibility and 
how they 'only' need 'limited' data. Segue this into needing 
only phylogeny (refer to Stock 2021 in here) Mention *RDPG* and 
*SVD*, *Network embedding* just name dropping 
Brief methods in the sense of we SVD (reduce
information to important info)
we RDPG (use latent space) and then we embed (to learn) -
refer to box one here *possibly* -> **Salomé**

How can we validate things if we don't have 
empirical networks. Either a ~~'throw away'~~ or a paragraph 
But there are some other general rules we could leverage. 
-> **Salomé**

Then some sort of 'mission statement' *i.e.* to test this 
method we wanted to see how well we could predict the Canadian 
mammal foodweb using the European foodweb [@Maiorano2020TetEu]
and maybe even a brief 
summary of result/value performance maybe also on validation 
and method?? **pending**

> Here we have Box 1: that provides a visual scaffold of the 
> methods workflow
## Methods

### Dataset preparation 
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

**SVD** what it is and t-SVD

**RDPG** what is it and how we apply/use it - 

Elbow method - check where this slots in

**Network embedding:** We extracted the left and right subspaces 
of a rank-5 *t-SVD* of the European mammals metaweb. When 
multiplied, these two matrices give an approximation of the 
network; we determined a threshold indicating that an interactions 
is likely by finding the value maximizing Youden's J statistic. 
At rank 5, this was approx. 0.11, giving J > 0.95.

**Latent traits inference** For all species in the Canadian pool 
that had also been reported in the European metaweb, we kept their 
latent traits as is. For other species, we reconstructed their 
latent traits by averaging the values of those of their 3 closest 
neighbors (based on the cophenetic matrix of the entire tree). 
This ensure that shared species are at the same position in both 
latent subspaces.

**Prediction** The prediction of interactions was done through a 
random dot product graph, by multiplication of the inferred left 
and right subspaces and thresholding at the value estimated during 
step 2. This yields a network with 223 species of mammals and approx. 
6000 interactions.

### Validation

**Validation** Omnivory v. trophic level diagram, body mass ratio,
Newfoundland tree [@Strong2014ImpNon]

## Results and discussion

> Figure 1 is a multi-paneled output of results: ROC, matrix??

Hey look lots of similarities between Europe and Canada and 
here are some metrics

> Figure 2 is a multi-paneled output of validation: body mass 
> ratio, omnivory/trophic comparison and Newfoundland metaweb

Newfoundland metaweb comparison as a matrix of sorts 
- colour interactions we missed 
and interactions we predicted that were not observed i.e. they 
could theoretically exist

Challenges of validation - even empirical webs are imperfect 
(FALSE negatives?? I always mix this up...) Need for a suite 
of validation techniques?

## Concluding notes

Something about how great this method is. Might be worth 
discussing ability to 'absorb' inclusion of new species. 
How can we make this prediction better? 

## References