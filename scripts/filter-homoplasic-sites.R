# Libraries
library(ape)
library(seqinr)
library(reshape2)
library(ggplot2)
library(broman)
library(cowplot)
library(dplyr)
library(chron)
source('functions.R')

# Parameters for running the script. 
# Note that some listed input files in input-parameters.R are *NOT INCLUDED* in this repository
source('../input-data/input-parameters.R')

# Make directories
dir.create(paste0(FIGURE.OUTPUT.FOLDER))
dir.create(paste0(FIGURE.OUTPUT.FOLDER, '/cophenetic-distributions/'))
dir.create(paste0(DATA.OUTPUT.FOLDER))

# Read in data
source('read-data.R')

# Select all homoplasies (should all have Min.No.ChangesonTree>0 anyway, but good to check)
all.homoplasic.sites <- snp.counts[which(snp.counts$Min.No.ChangesonTree>0),"bp"]

# 1. Histogram of all homoplasic sites by position in genome
pdf(paste0(FIGURE.OUTPUT.FOLDER, '/homoplasic-sites-histogram-position.pdf'))
hist(all.homoplasic.sites, col='black', xlab='Position in genome', main='',
      breaks=seq(0, 30000, 500))
dev.off()

# Create dataset, ignoring masked regions
homoplasic.counts <- snp.counts[which(snp.counts$Min.No.ChangesonTree>0 & 
                                                                !snp.counts$bp %in% MASKED.REGIONS),]

# Add nearest potential homoplasy 
homoplasic.counts$dist.nearest.homoplasy <- sapply(homoplasic.counts$bp, 
                                                   function(x) sort(abs(x-homoplasic.counts$bp), decreasing = FALSE)[2])


# 2. Homoplasic sites 
# excluding first and last N bp of genome


# 3. Go through and record whether an isolate's nearest neighbour has the homoplasy. 
homoplasic.counts$N.nearest.neighbour.has.homoplasy <- NA
homoplasic.counts$N.nearest.neighbour.lacks.homoplasy <- NA
homoplasic.counts$N.nearest.neighbour.intermediate <- NA
homoplasic.counts$proportion.with.N.within.local.region <- NA
homoplasic.counts$isolates.in.dataset.with.N <- NA

for (i in seq(1, nrow(homoplasic.counts))){
  h <- homoplasic.counts$bp[i]
  print(h)
  isolates <- getIsolatesWithMinorVariant(h)
  if (length(isolates)!=0){
    results <- sapply(isolates, 
                      function(x) nearestNeighbourHasVariant(x, h))
    homoplasic.counts$N.nearest.neighbour.lacks.homoplasy[i] <- length(results[which(results==0)])
    homoplasic.counts$N.nearest.neighbour.has.homoplasy[i] <- length(results[which(results==1)])
    homoplasic.counts$N.nearest.neighbour.intermediate[i] <- length(results[which(results!=1 & results!=0)])
    
    # Add adjacent N score
    n.scores <- sapply(isolates, 
                       function(x) getAdjacentNscore(isolate = x, site = h, region = LOCAL.REGION))
    homoplasic.counts$proportion.with.N.within.local.region[i] <- length(n.scores[which(n.scores>0)])/length(n.scores)
    
  }

  homoplasic.counts$isolates.in.dataset.with.N[i] <- length(grep("n", sapply(aln$nam, function(x) substr(aln$seq[[which(aln$nam ==x)]], h, h))))
}
# Calculate proportion of isolates with homoplasy with nearest neighbour with the homoplasy
homoplasic.counts$N.isolates.with.homoplasy <- homoplasic.counts$N.nearest.neighbour.has.homoplasy+homoplasic.counts$N.nearest.neighbour.lacks.homoplasy+homoplasic.counts$N.nearest.neighbour.intermediate
homoplasic.counts$proportion.nearest.neighbour.has.homoplasy <- homoplasic.counts$N.nearest.neighbour.has.homoplasy/homoplasic.counts$N.isolates.with.homoplasy

# Histogram of these distances
pdf(paste0(FIGURE.OUTPUT.FOLDER, '/histogram-homoplasic-sites-nearest-neighbour-proportion-all-homoplasies.pdf'))
hist(homoplasic.counts$proportion.nearest.neighbour.has.homoplasy, 
     breaks=100, 
     col='black', 
     xlab='Proportion of isolates where nearest neighbour in tree has homoplasy',
     main='')
dev.off()
# Show excluding start/end of genome
# N <- 500 
# homoplasic.counts.filt <- homoplasic.counts[which(homoplasic.counts$bp>N & homoplasic.counts$bp<GENOME.LENGTH-N),]
# pdf(paste0(FIGURE.OUTPUT.FOLDER, '/histogram-homoplasic-sites-nearest-neighbour-proportion-excluding-first-last-500-bp.pdf'))
# hist(homoplasic.counts$proportion.nearest.neighbour.has.homoplasy, 
#      breaks=100, 
#      col='black', 
#      xlab='Proportion of isolates where nearest neighbour in tree has homoplasy',
#      main='')
# dev.off()

# Write to file
write.csv(homoplasic.counts, file=paste0(DATA.OUTPUT.FOLDER, '/all-homoplasic-sites-table.csv'))


# 4. Further homoplasy filtering thresholds

# Proportion of nearest neighbours with homoplasy (p_nn in manuscript)
homoplasic.counts.filt.HQ <- homoplasic.counts[which(homoplasic.counts$proportion.nearest.neighbour.has.homoplasy>NEAREST.HOMOPLASY.PROP),]
# Number of isolates with homoplasy 
homoplasic.counts.filt.HQ <- homoplasic.counts.filt.HQ[which(  homoplasic.counts.filt.HQ$N.isolates.with.homoplasy>N.ISOLATES.WITH.HOMOPLASY),]
# No isolates with homoplasy with N in local region
homoplasic.counts.filt.HQ <- homoplasic.counts.filt.HQ[which(  homoplasic.counts.filt.HQ$proportion.with.N.within.local.region==0),]
homoplasic.counts.filt.HQ$proportion.with.N.within.local.region <- NULL # Don't need variable any more

# Add number of originating labs & countries
originating_labs <- c()
submitting_labs <- c()
countries <- c()
for (site in homoplasic.counts.filt.HQ$bp){
  isolates <- getIsolatesWithMinorVariant(site)
  originating_labs <- c(originating_labs, length(table(metadata[which(metadata$gisaid_epi_isl %in% isolates), "originating_lab"])))
  submitting_labs <- c(submitting_labs, length(table(metadata[which(metadata$gisaid_epi_isl %in% isolates), "submitting_lab"])))
  
  countries <- c(countries, length(table(metadata[which(metadata$gisaid_epi_isl %in% isolates), "country"])))
  print(site)
}
homoplasic.counts.filt.HQ$N.countries <- countries
homoplasic.counts.filt.HQ$N.originating.labs <- originating_labs
homoplasic.counts.filt.HQ$N.submitting.labs <- submitting_labs

# Filter to more than one submitting lab
homoplasic.counts.filt.HQ <- homoplasic.counts.filt.HQ[which(homoplasic.counts.filt.HQ$N.submitting.labs>N.SUBMITTING.LABS & 
                                                               homoplasic.counts.filt.HQ$N.originating.labs>N.SUBMITTING.LABS),]

# Write to file
write.csv(homoplasic.counts.filt.HQ, file=paste0(DATA.OUTPUT.FOLDER, '/filtered-homoplasic-sites-table.csv'))

# 5. Make plots 
# For these filtered homoplasies, make plots
for (h in homoplasic.counts.filt.HQ$bp){
  print(h)
  h.df <- getCopheneticDistributionForHomoplasy(site = h)
  h.plot <- plotHomoplasyCopheneticDistribution(h.df, title=getTitleString(h))
  ggsave(h.plot, file=paste0(FIGURE.OUTPUT.FOLDER, '/cophenetic-distributions/', h, '.pdf'), 
         width=9, height=6)
}

