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

**A note on filtering thresholds.** These are available in the main script. See manuscript for more discussion of the rationale.  

## Outputs

Outputs are stored in `figures` and `output-data`. They include:

* Histograms of homoplasy distribution
* Histograms of cophenetic distances between isolates sharing homoplasy for all filtered sites, compared to all cophenetic distances in tree, with useful summary statistics
* Dataset of all identified homoplasies, with statistics used for filtering
* Dataset of filtered putative homoplasic sites according to specified filtering thresholds

## Data availability

Underlying data (assemblies) comes from the GISAID consortium. A full acknowledgements table of laboratories is available in `acknowledgements.tsv`.

The multiple sequence alignment is not included in this repository as per the terms of the GISAID consortium for sharing sequence data ('You agree not to distribute Data to any third party other than Authorized Users as contemplated by this Agreement.') In order to obtain access to the assemblies in the multiple sequence alignment, you can register as a user of GISAID [here](https://www.gisaid.org/registration/register/). 
