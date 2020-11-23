# Libraries
library(ape)
library(seqinr)
library(reshape2)
library(ggplot2)
library(broman)
library(cowplot)
library(dplyr)
library(chron)


alt_alleles_at_site_mod <- read.table("./dataset_23/gisaid_cov2020_sequences.30.07.2020.QC.human.nextstrain_filter.QC.NSmask.HPfinder.noambig.BASECOUNTS_edited.csv", row.names=NULL, header=T, sep=",",dec=".", stringsAsFactors=F,check.names=FALSE)

# Parameters for running the script. 
# Note that some listed input files in input-parameters.R are *NOT INCLUDED* in this repository
source('../input-data/input-parameters.R')

source('functions.R')
# Read in data
source('read-data.R')


# Make directories
dir.create(paste0(FIGURE.OUTPUT.FOLDER))
dir.create(paste0(DATA.OUTPUT.FOLDER))

# Select all homoplasies (should all have Min.No.ChangesonTree>0 anyway, but good to check)
all.homoplasic.sites <- snp.counts[which(snp.counts$Min.No.ChangesonTree>0),"bp"]

# 1. Histogram of all homoplasic sites by position in genome
pdf(paste0(FIGURE.OUTPUT.FOLDER, '/homoplasic-sites-histogram-position.pdf'))
hist(as.numeric(all.homoplasic.sites), col='black', xlab='Position in genome', main='',
      breaks=seq(0, 30000, 500))
dev.off()

# Create dataset, ignoring masked regions
homoplasic.counts <- snp.counts[which(snp.counts$Min.No.ChangesonTree>0 & 
                                                                !snp.counts$bp %in% MASKED.REGIONS),]

# Add nearest potential homoplasy 
homoplasic.counts$bp <- as.numeric(homoplasic.counts$bp)
homoplasic.counts$dist.nearest.homoplasy <- sapply(homoplasic.counts$bp, 
                                                   function(x) sort(abs(x-homoplasic.counts$bp), decreasing = FALSE)[2])

homoplasic.counts$N.nearest.neighbour.has.homoplasy <- NA
homoplasic.counts$N.nearest.neighbour.lacks.homoplasy <- NA
homoplasic.counts$N.nearest.neighbour.intermediate <- NA
homoplasic.counts$proportion.with.N.within.local.region <- NA
homoplasic.counts$isolates.in.dataset.with.N <- NA
homoplasic.counts$N.isolates.with.homoplasy <- NA
homoplasic.counts$originating_labs <- NA
homoplasic.counts$submitting_labs <- NA
homoplasic.counts$countries <- NA

# to gain time reduce the metadata table to the isolates actually in the dataset
metadata <- metadata[metadata$gisaid_epi_isl %in% tree$t,]

for (i in seq(1, nrow(homoplasic.counts))){
  h <- homoplasic.counts$bp[i]
  print(h)
  isolates <- getIsolatesWithMinorVariant(h)
  homoplasic.counts$N.isolates.with.homoplasy[i] <- length(isolates)
  # computing number of labs, countries etc takes time so i'll do it only for 

  if(length(isolates) > N.ISOLATES.WITH.HOMOPLASY){
  homoplasic.counts$originating_labs[i] <- length(unique(metadata[which(metadata$gisaid_epi_isl %in% isolates), "originating_lab"]))
homoplasic.counts$submitting_labs[i] <-  length(unique(metadata[which(metadata$gisaid_epi_isl %in% isolates), "submitting_lab"]))
homoplasic.counts$countries[i] <-  length(unique(metadata[which(metadata$gisaid_epi_isl %in% isolates), "country"]))
  }else{
    homoplasic.counts$originating_labs[i] <- 0
    homoplasic.counts$submitting_labs[i] <-  0
    homoplasic.counts$countries[i] <-  0
}
}

isolates <- tree$t[seq(from=1,to=10000)]

homoplasic.counts$proportion_alt2_vs_alt1 <- NA
homoplasic.counts$proportion_alt2_vs_alt1 <- alt_alleles_at_site_mod[homoplasic.counts$bp,"prop_alt2_relative_alt1"]

# Write to file
write.csv(homoplasic.counts, file=paste0(DATA.OUTPUT.FOLDER, '/all-homoplasic-sites-table.csv'))

# 4. Further homoplasy filtering thresholds
# Number of isolates with homoplasy 
homoplasic.counts.filt.HQ_0.05 <- homoplasic.counts[which(homoplasic.counts$proportion_alt2_vs_alt1<0.05 &
                               homoplasic.counts$submitting_labs>N.SUBMITTING.LABS & 
                               homoplasic.counts$originating_labs>N.SUBMITTING.LABS &
                               homoplasic.counts$N.isolates.with.homoplasy>N.ISOLATES.WITH.HOMOPLASY),]

write.csv(homoplasic.counts.filt.HQ_0.05, file=paste0(DATA.OUTPUT.FOLDER, '/filtered-homoplasic-sites-table_alt2_vs_alt1_0.05.csv'))

homoplasic.counts.filt.HQ_0.1 <- homoplasic.counts[which(homoplasic.counts$proportion_alt2_vs_alt1<0.1 &
                                                           homoplasic.counts$submitting_labs>N.SUBMITTING.LABS & 
                                                           homoplasic.counts$originating_labs>N.SUBMITTING.LABS &
                                                           homoplasic.counts$N.isolates.with.homoplasy>N.ISOLATES.WITH.HOMOPLASY),]

write.csv(homoplasic.counts.filt.HQ_0.1, file=paste0(DATA.OUTPUT.FOLDER, '/filtered-homoplasic-sites-table_alt2_vs_alt1_0.1.csv'))

homoplasic.counts.filt.HQ_0.2 <- homoplasic.counts[which(homoplasic.counts$proportion_alt2_vs_alt1<0.2 &
                                                           homoplasic.counts$submitting_labs>N.SUBMITTING.LABS & 
                                                           homoplasic.counts$originating_labs>N.SUBMITTING.LABS &
                                                           homoplasic.counts$N.isolates.with.homoplasy>N.ISOLATES.WITH.HOMOPLASY),]

write.csv(homoplasic.counts.filt.HQ_0.2, file=paste0(DATA.OUTPUT.FOLDER, '/filtered-homoplasic-sites-table_alt2_vs_alt1_0.2.csv'))



