# Homoplasy filtering

This repository has a set of scripts for filtering out potential sequencing artefacts from assemblies in order to identify putative homoplasic sites in SARS-CoV-2.

It needs as input:

* A phylogenetic tree
* A multiple sequence alignment 
* A dataset of homoplasy counts (generated with MPBoot - HomoplasyFinder)

Some of the outputs:

* Histograms of homoplasy distribution
* Histograms of cophenetic distances between isolates sharing homoplasy, compared to all cophenetic distances in tree
* Dataset of filtered putative homoplasic sites
