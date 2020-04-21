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
source('read-data.R')

# Select all homoplasies (Min.No.ChangesonTree>2,
GENOME.LENGTH <- 29903
all.homoplasic.sites <- snp.counts[which(snp.counts$Min.No.ChangesonTree>0),"bp"]

# 1. Histogram of all homoplasic sites by position in genome
pdf('../figures/homoplasic-sites-histogram-position.pdf')
hist(all.homoplasic.sites, col='black', xlab='Position in genome', main='',
      breaks=seq(0, 30000, 500))
dev.off()

# Create dataset
homoplasic.counts <- snp.counts[all.homoplasic.sites,]

# Add nearest potential homoplasy 
homoplasic.counts$dist.nearest.homoplasy <- sapply(homoplasic.counts$bp, 
                                                   function(x) sort(abs(x-homoplasic.counts$bp), decreasing = FALSE)[2])


# 2. Homoplasic sites 
# excluding first and last N bp of genome
N <- 500 
homoplasic.counts.filt <- homoplasic.counts[which(homoplasic.counts$bp>N & homoplasic.counts$bp<GENOME.LENGTH-N),]



# 3. Go through and record whether an isolate's nearest eighbour has the homoplasy. 
homoplasic.counts.filt$N.nearest.neighbour.has.homoplasy <- NA
homoplasic.counts.filt$N.nearest.neighbour.lacks.homoplasy <- NA
homoplasic.counts.filt$N.nearest.neighbour.intermediate <- NA

for (i in seq(1, nrow(homoplasic.counts.filt))){
  h <- homoplasic.counts.filt$bp[i]
  print(h)
  isolates <- getIsolatesWithMinorVariant(h)
  results <- sapply(isolates, 
                    function(x) nearestNeighbourHasVariant(x, h))
  homoplasic.counts.filt$N.nearest.neighbour.lacks.homoplasy[i] <- length(results[which(results==0)])
  homoplasic.counts.filt$N.nearest.neighbour.has.homoplasy[i] <- length(results[which(results==1)])
  homoplasic.counts.filt$N.nearest.neighbour.intermediate <- length(results[which(results!=1 & results!=0)])
}
# Calculate proportion of isolates with homoplasy with nearest neighbour with the homoplasy
homoplasic.counts.filt$N.isolates.with.homoplasy <- homoplasic.counts.filt$N.nearest.neighbour.has.homoplasy+homoplasic.counts.filt$N.nearest.neighbour.lacks.homoplasy+homoplasic.counts.filt$N.nearest.neighbour.intermediate
homoplasic.counts.filt$proportion.nearest.neighbour.has.homoplasy <- homoplasic.counts.filt$N.nearest.neighbour.has.homoplasy/homoplasic.counts.filt$N.isolates.with.homoplasy

# Histogram of these distances
pdf('../figures/histogram-homoplasic-sites-nearest-neighbour-proportion.pdf')
hist(homoplasic.counts.filt$proportion.nearest.neighbour.has.homoplasy, 
     breaks=100, 
     col='black', 
     xlab='Proportion of isolates where nearest neighbour in tree has homoplasy',
     main='')
dev.off()

# Write to file
write.csv(homoplasic.counts.filt, file='../output-data/homoplasic-SNP-counts-filtered.csv')

# 4. Make plots 
# Exclude the homoplasies with < 0.1 proportion (i.e. the large peak at 0)
homoplasic.counts.filt <- homoplasic.counts.filt[which(homoplasic.counts.filt$proportion.nearest.neighbour.has.homoplasy>0.1),]

# Further high-quality homoplasy thresholds
#DISTANCE.TO.HOMOPLASY <- 10
# High-quality homoplasy thresholds
PROPORTION.NEAREST.NEIGHBOUR.HAS.HOMOPLASY <- 0.4
N.ISOLATES.WITH.HOMOPLASY <- 10
homoplasic.counts.filt.HQ <- homoplasic.counts.filt[which( homoplasic.counts.filt$proportion.nearest.neighbour.has.homoplasy>PROPORTION.NEAREST.NEIGHBOUR.HAS.HOMOPLASY &
                                                             homoplasic.counts.filt$N.isolates.with.homoplasy>N.ISOLATES.WITH.HOMOPLASY),]
write.csv(homoplasic.counts.filt.HQ, file='../output-data/homoplasic-SNP-counts-filtered-HQ.csv')



# For these, make plots
for (h in homoplasic.counts.filt.HQ$bp){
  print(h)
  h.df <- getCopheneticDistributionForHomoplasy(site = h)
  h.plot <- plotHomoplasyCopheneticDistribution(h.df, title=getTitleString(h))
  ggsave(h.plot, file=paste0('../figures/cophenetic-distributions/', h, '.pdf'))
}





