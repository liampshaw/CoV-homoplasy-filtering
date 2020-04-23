# Homoplasy filtering

This repository has a set of scripts for filtering out potential sequencing artefacts from assemblies in order to identify putative homoplasic sites in SARS-CoV-2.

It needs as input:

* A phylogenetic tree
* A multiple sequence alignment 
* A dataset of homoplasy counts (generated with MPBoot - HomoplasyFinder)

The main script `filter-homoplasic-sites.R` is run from the `scripts` directory. It assumes that you have a `data` directory with your input data and that the names of isolates in the tree and alignment are the same. As the script was developed for a specific analysis at speed, it will almost certainly break on new data. (Currently there is some hacky stuff to make sure names are the same, to be updated to work on a general dataset.) Please use with caution! 

Some of the outputs:

* Histograms of homoplasy distribution
* Histograms of cophenetic distances between isolates sharing homoplasy for all filtered sites, compared to all cophenetic distances in tree, with useful summary statistics
* Dataset of all identified homoplasies with some statistics used for filtering
* Dataset of filtered putative homoplasic sites

Liam Shaw, liam.philip.shaw [at] gmail
