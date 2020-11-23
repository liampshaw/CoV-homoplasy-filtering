[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# SARS-CoV-2 homoplasy filtering

This repository contains a set of R scripts which attempt to filter out spurious homoplasic sites from SARS-CoV-2 assemblies. 

For more information, please see the associated original publication (commit  `0fa7434`) which identified a list of possible recurrent mutations: 
 
**[Manuscript1]** *Emergence of genomic diversity and recurrent mutations in SARS-CoV-2*  
L. van Dorp, M. Acman\*, D. Richard\*, L. P. Shaw\*, C. E. Ford, L. Ormond, C. J. Owen, J. Pangae, C. C. S. Tan, F. A.T. Boshiere, A. Torres Ortiz, F. Balloux  

*Infection, Genetics and Evolution* (2020) doi: [10.1016/j.meegid.2020.104351](https://doi.org/10.1016/j.meegid.2020.104351)  
\* = equal contribution  

Please note that the current version of this repository (commit `e0f6dfe`) is associated with the below publication, which analyses the effect of recurrent mutations on transmissibility:

**[Manuscript2]** *No evidence for increased transmissibility from recurrent mutations in SARS-CoV-2*  
L. van Dorp\*, D. Richard\*, C. C. S. Tan, L. P. Shaw, M. Acman, F. Balloux  

*Nature Communications* (2020) doi: [10.1038/s41467-020-19818-2](https://doi.org/10.1038/s41467-020-19818-2)

Liam Shaw (liam.philip.shaw@gmail.com) & Damien Richard (richarddamienfr@gmail.com)

## Overview

Ongoing analyses of public SARS-CoV-2 genomes have identified a high number of recurrent mutations i.e. those occurring more than once on the phylogenetic tree ([homoplasies](https://en.wikipedia.org/wiki/Homoplasy)). These are potentially of interest for considering the ongoing evolution of SARS-CoV-2. However, apparent homoplasies can also be caused by sequencing errors. It is therefore important to take this into account in any analysis. 

 This repository contains a set of R scripts which attempt to filter out spurious homoplasic sites from SARS-CoV-2 assemblies. For example, in the analyses in Manuscript1 (see above) we went from an initial list of 1132 homoplasies identified using HomoplasyFinder to a reduced set of 198 (you can see this version at `git checkout 0fa7434`). 

**It is possible that even homoplasies which pass these filtering thresholds may be artefactual/spurious, so please use and interpret results with caution.** Additional sequencing artefacts are being identified from other independent analyses and we would expect the list of suspicious sites to grow. (**Edit 23/11/2020:** indeed, this happened.) For example, see the report of May 5 2020 from [de Maio et al.](http://virological.org/t/issues-with-sars-cov-2-sequencing-data/473) on Virological, subsequent report on 14 July 2020, and the maintained list of problematic sites maintained by this group [here](https://github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf).


## Inputs

The script takes input specified in `input-data/input-parameters.R`:

* A dataset of homoplasy counts - generated with HomoplasyFinder from a maximum parsimony tree (MPBoot), with some additional metadata added (e.g. which protein a site is in) [example included]
* A phylogenetic tree - generated with RaxML from a multiple sequence alignment [not included, see below]
* A multiple sequence alignment - generated from GISAID assemblies [not included, see below] 

The main script `filter-homoplasic-sites.R` is run from within the `scripts` directory. It assumes that you have a `input-data` directory with your input data and that the names of isolates in the tree and alignment are the same. 

As the script was developed with a specific analysis in mind (n=7666 assemblies downloaded from GISAID as of April 19 2020), it will **almost certainly break** on new data. We adapted it for our secondary publication. Please use with caution if this is your intention. 

**A note on filtering thresholds.** There are several filtering thresholds used. These values and notes below apply to the analyses in Manuscript1. The method was refined slightly for later analyses (see discussion in methods of Manuscript2)

* Number of isolates with homoplasy.  
**Default: >8**  
*Rationale*: small numbers of isolates more likely to be due to random sequencing error. 8 is an arbitrary number equivalent to 0.1% of dataset (n=7666) 
* Position *x* of homoplasy in genome.  
**Default: exclude masked sites: {{1, ..., 150}, 18529, 29849, 29851, 29853, {29853..29903}}** 
*Rationale*: sequencing error appears more common at the start and end of the SARS-CoV-2 genome, causing high density of apparent homoplasies, and at other known sites. Taken from [NextStrain masking](https://github.com/nextstrain/ncov/blob/20974e93a647cc17595718561a69eb02a42a4e5a/config/config.yaml). 
* Proportion of isolates with homoplasy which have a nearest neighbour in the tree with the homoplasy: ranges between 0 (singleton isolates with homoplasy throughout tree) and 1 (clusters of isolates with homoplasy).  
**Default: >0.1**  
*Rationale*: apparent homoplasies caused by random sequencing error are *a priori* unlikely to cluster with each other in the tree. 
* Proportion of isolates with the homoplasy which have at least one 'N' in the local region around the homoplasy.  
**Default: zero, with 'local' defined as +/- 5bp**  
*Rationale*: 'N' in local region could be suggestive of hard-to-sequence region. 
* Number of submitting and originating labs.  
**Default: >1**  
*Rationale*: homoplasies only present in isolates from a single originating or submitting lab could be due to batch effects.
 
These thresholds can be changed in the main script. See Manuscript1 for more discussion of the rationale. 

## Outputs

Outputs are stored in `figures` and `output-data`. They include:

* Histograms of homoplasy distribution
* Histograms of cophenetic distances between isolates sharing homoplasy for all filtered sites, compared to all cophenetic distances in tree, with useful summary statistics
* Dataset of all identified homoplasies, with statistics used for filtering
* Dataset of filtered putative homoplasic sites according to specified filtering thresholds

## Data availability

Underlying data (n=7666 assemblies) for the input files comes from the GISAID consortium. A full acknowledgements list for Manuscript1 of originating and submitting laboratories (as of April 23 2020) is available in `acknowledgements.tsv`.  

The tree and multiple sequence alignment are not included in this repository as per the terms of the GISAID consortium for sharing sequence data ('You agree not to distribute Data to any third party other than Authorized Users as contemplated by this Agreement.') In order to obtain access to the assemblies, you can register as a user of GISAID [here](https://www.gisaid.org/registration/register/). 

N.B. Protein annotations in the homoplasy count `input-data` file are based on the genbank file of Wuhan-Hu-1 (NCBI NC_045512.2) and might become out-of-date as future annotation changes. They are not used in any of the filtering. 
