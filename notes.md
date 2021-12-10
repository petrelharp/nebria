# Dec 6:

- Yosemite predictions with a GLM come out to 122 for the total range
- but this is on the low end, probably we want between 100 and 500;
    we could do a Bayesian thing to get a better posterior,
    but why not just uniform over that range

# Nov 22

- do ABC with just patch number, extrapolating from Yosemite to a number of patches
    for the whole range using the SDM
- change patch number to get number of pixels with more than 10 individuals
- ABC package abc

# Nov 15:

- count number of occupied patches
- use that for ABC along with genomic data
- want population (effective) size to decline with time
- main goal: recolonization path for where beetles came from:
    where were ancestors of today's beetles living at 21kya
- want to learn demographic history: how many beetles across time?

# Sep 27:

- add stochasticity for goodness of patch: perhaps draw independently by year 
    * on a coarse scale (like 50km?)? but no, because we don't want big patches to dissapear entirely
    * so, draw independently by pixel
- the glacier boundary is from about 21kya
- cut glacier out of only oldest map
- Ottoway has low heterozygosity, probably because it's very hard to find even the three beetles
    they sampled.

# Sep 13:

- note that when we compare to for instance pi, if it has a lot of local randomness
    we should not compare individual sites to each other, but rather compare the
    overall variance (or variance between neighbors)
- so: don't pair up populations!!
- estimate proportion of patches that are actually occupied in 60 lake basin
- upper limit on total size: maybe 1e5 or 2e5 today: maybe 5e5 over all time?
- stairwayplot got upper size limit of 120K at 100kya


# Aug 30

- lower diversity in some peripheral regions of yosemite
- sixty lakes and taboose pass lower dxy to neighbors - admixture/IBD
- LD decays from 0.25 to 0.125 over like 10kb, then flattens out to varying levels

# Aug 19:

- Fst between distant pops is like 0.5, between nearby is like 0.03.
- a few hundred mis-called alleles per individual
- many more mis-called homs that are really hets
- next meeting: discuss manuscript


# Jul 29:

- each sampling location is a collection of nearby locations,
    but a few of these (Lucy pass, Seldon pass, Paiute pass, ...) are a few km apart
- some other sites are pretty close but separated by a ridge, for instance
- Yi-Ming to send a sheet with the sub-sample locations (without duplicates)
- maybe want to compute e.g. AFS for bigger regions to better group samples?
- pay attention to occupancy of marginal pops
- take snapshots every 1000 years or so
- Yi-Ming: whole-genome including X: Fst is 0 to 0.36
- genome size:
    * VCF has 2M SNPs after removing alleles < 5%; with these it's 5M SNPs
    * don't know the recombination rate
    * has between 30 and 50 chromosomes
    * 147 Mb genome
- question: do we do LD thinning for comparison to sims?
- do "population-specific Fst": (one pop) vs (everyone else) Fst
- TODO: figure out stats like Fst that don't depend on singletons

# Jul 20:

Statistics:

- mean heterozygosity
- diversity by pop
- divergence between pops (or, Fst?) against distance
- AFS by pop without singletons
- PCA

Notes:

- use fact that the big area north of the San Joaquin was surveyed a bunch but none are found there
- dispersal over snow (any direction, as adults) and juveniles (in streams)
- maybe not long-distance dispersal because they are flightless (although: washed downriver over 1km)

To-do:

- make Yosemite subset
- Yi-Ming to send map & numbers of samples
- pre-compute spatial map of density?
- make each patch suitability fluctuate with time on a hundreds - 500 years time scale

Timing (rough!):

- 1 sec/generation with 5,000 adults and 35,000 juvys; scaling linearly; that's 5.5 hours for going back to 20Kya

From Yi-Ming:

- `1000M_plus.tif`: It's supposed to be a "geographic-only" raster with just slope but I have tried many different combinations of many geographic features and none of them improved the map (either similar or have less suitable habitat than the original one). Sean and I thought maybe we can try the blank raster with all pixels above elevation 1000 to be 1, and below 1000 to be NA. I don't know this is the best way to do this or not, so please let me know if you have different idea.
- `CCSM4_glm_boundary`: The LGM suitability map plus the glacier boundary. I merged the two rasters so now the values 0-1 are suitability and value 2 is the boundary. I am not sure this is the best format for you, if you need different format (like separated rasters for suitability and glacier boundary), please let me know, too. 

# Jul 15

- make dispersal depend on density (maybe not! shouldn't matter much)
- make some patches more than one pixel
- survivorship maybe 10-25%, so like half of adults dispersing
- note replacement rate is 2 offspring/female
- since 4% of habitat is good, no need to bias dispersal towards good places,
    since should have >= 1 surviving dispersers per patch
- effective migration rate can't be much lower than 1 per gen because they're not super inbred
- include glacier in habitat raster (and also dispersal rasters?)
- should differentiate suitability from dispersal rasters
- glacier was only longest-ago time?

# Jun 10 2021

- subdivide pixels into 5x5 and set 24/25th of the pixels to zero
- make probability of moving to a new spot equal to suitability map
- make them not stop at original habiat
- goal is to identify refugia and learn migration rates
- interpolate linearly between time points
- with each map at the midway of their time points

# From 29 April 2021

Peter, proposal:
- have one 'site' per grid cell by discretizing movement OR using NearestNeighborOfPoint
- give each site a composite 'suitability' by multiplying a "topography suitability"
    (obtained from slope & aspect or generated randomly, autocorrelated)
    by a "climate suitability" obtained from glaciers
- make dispersal stick to local features (water and/or suitability) by breaking it into
    several small steps
- do interpolation by computing distance to each boundary... ?
- Q: was current biodiversity shaped by refugia? or, just the landscape?
- Q: is it (really) secondary contact?
- Q: how patchy is the modern distribution? are current populations connected or not?
- Q: how does movement depend on continuous data,
     like how much snow there is present between two sites (measured by SWE)?

# From 22 Apr 2021

- fit maxent models: first one looks pretty good
- CCSM4: 21Kya - a bunch of ice-free patches within the glacier (nunataks) but very inhospitable
- nunataks are with ice on both sides; east of a block of ice
- except there's one larger ice-free patch
- higher-res version has 1km2 pixels
- two time snapshots of high-res; nine of low-res
- deglaciation ends and holocene starts at 10.5Kya (hundreds or 1000 years of deglaciation)
- how to interpolate between old map and new map?

- number of eggs per year: 15-20 (check it with yi-ming)
- females lay eggs once a year
- sex ratio: more females than males, 60/40 (common in beetles) - due to spiroplasma killing males?
- (but also males are bigger and take longer to develop so have higher risk (but not for this species))
- option1: egg laid in summer, overwinters as larvae, reproduces next year
- option2: or overwinters as an adult and reproduces at year 2
- probably reproduce *only once*, with some probability of doing it at year 1 versus 2
- dispersal: area of ideal habitat is like dozens of meters (smaller than a football field),
    not found outside of a few meters outside of the good habitat
- but maybe passive dispersal downstreams - like to walk in streams (and stick to streams)
- juveniles are not found in the field
- rare long distance: down/up streams; occasional long walkabouts of like 500m?
- pop size in a patch: cannibalism, food availability, predation (eg by other Nebria)
- multiple mating
- has rough census size estimates: tens to hundreds up to a thousand
- larvae are sit-and-wait predators in water and also in/under snow
- larvae can float on the water
- larvae live in tunnels


# From 12 April 2021

- live above 2900m in vertical seeps

- 100s to 1000s per site;
- 0.1 to 0.25 of good sites are occupied;
- roughly one site per sq mile in good habitat;
- roughly hundreds of square miles of good habitat

- one generation = one year

- see Mol Ecol paper
