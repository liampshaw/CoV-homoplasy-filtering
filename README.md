# Homoplasy filtering

This repository contains R scripts for filtering out potential sequencing artefacts from assemblies in order to identify putative homoplasic sites in SARS-CoV-2.

Associated publication: **to be included** 

Liam Shaw, liam.philip.shaw@gmail.com


## Inputs

The script takes as input:

* A phylogenetic tree - generated with RaxML from a multiple sequence alignment
* A dataset of homoplasy counts - generated with HomoplasyFinder from a maximum parsimony tree (MPBoot), with some additional metadata added (e.g. which protein a site is in)
* A multiple sequence alignment - generated from GISAID assemblies [not included, see below] 

The main script `filter-homoplasic-sites.R` is run from within the `scripts` directory. It assumes that you have a `input-data` directory with your input data and that the names of isolates in the tree and alignment are the same. As the script was developed with a specific analysis in mind, it will **almost certainly break** on new data. Please use with caution if this is your intention. 

**A note on filtering thresholds.** There are several filtering thresholds used (see manuscript for more discussion).

* Number of isolates with homoplasy. Default: >10  
Rationale: small numbers of isolates more likely to be due to random sequencing error. 
* Position **x** of homoplasy in genome. Default: 500<**x**<(29903-500).  
Rationale: sequencing error appears more common at the start and end of the SARS-CoV-2 genome, causing high density of apparent homoplasies.  
* Proportion of isolates with homoplasy which have a nearest neighbour in the tree with the homoplasy: ranges between 0 (singleton isolates with homoplasy throughout tree) and 1 (clusters of isolates with homoplasy). Default: >0.1
  Rationale: apparent homoplasies caused by random sequencing error are *a priori* unlikely to cluster with each other in the tree. 
* Proportion of isolates with the homoplasy which have at least one 'N' in the region +/- 2 bp around the homoplasy. Default: equal to zero
  Rationale: 'N' in local region could be suggestive of hard-to-sequence region. 
 
These thresholds can be changed in the main script. See manuscript for more discussion of the rationale.  

## Outputs

Outputs are stored in `figures` and `output-data`. They include:

* Histograms of homoplasy distribution
* Histograms of cophenetic distances between isolates sharing homoplasy for all filtered sites, compared to all cophenetic distances in tree, with useful summary statistics
* Dataset of all identified homoplasies, with statistics used for filtering
* Dataset of filtered putative homoplasic sites according to specified filtering thresholds

## Data availability

Underlying data (n=6971 assemblies) for the input files comes from the GISAID consortium. A full acknowledgements table of laboratories is available in `acknowledgements.tsv`.

The multiple sequence alignment is not included in this repository as per the terms of the GISAID consortium for sharing sequence data ('You agree not to distribute Data to any third party other than Authorized Users as contemplated by this Agreement.') In order to obtain access to the assemblies in the multiple sequence alignment, you can register as a user of GISAID [here](https://www.gisaid.org/registration/register/). 
