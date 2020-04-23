# Libraries
library(ape)
library(seqinr)
library(reshape2)
library(ggplot2)
library(broman)
library(cowplot)
library(dplyr)
library(chron)

GENOME.LENGTH <- 29903

source('functions.R')
source('read-data.R')

# Select all homoplasies (Min.No.ChangesonTree>2,
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


# 3. Go through and record whether an isolate's nearest neighbour has the homoplasy. 
homoplasic.counts$N.nearest.neighbour.has.homoplasy <- NA
homoplasic.counts$N.nearest.neighbour.lacks.homoplasy <- NA
homoplasic.counts$N.nearest.neighbour.intermediate <- NA

for (i in seq(1, nrow(homoplasic.counts))){
  h <- homoplasic.counts$bp[i]
  print(h)
  isolates <- getIsolatesWithMinorVariant(h)
  results <- sapply(isolates, 
                    function(x) nearestNeighbourHasVariant(x, h))
  homoplasic.counts$N.nearest.neighbour.lacks.homoplasy[i] <- length(results[which(results==0)])
  homoplasic.counts$N.nearest.neighbour.has.homoplasy[i] <- length(results[which(results==1)])
  homoplasic.counts$N.nearest.neighbour.intermediate[i] <- length(results[which(results!=1 & results!=0)])
  
  # Add adjacent N score
  n.scores <- sapply(isolates, 
                    function(x) getAdjacentNscore(isolate = x, site = h))
  homoplasic.counts$proportion.with.N.within.2.bp[i] <- length(n.scores[which(n.scores>0)])/length(n.scores)
  
  homoplasic.counts$isolates.in.dataset.with.N[i] <- length(grep("n", sapply(aln$nam, function(x) substr(aln$seq[[which(aln$nam ==x)]], h, h))))
}
# Calculate proportion of isolates with homoplasy with nearest neighbour with the homoplasy
homoplasic.counts$N.isolates.with.homoplasy <- homoplasic.counts$N.nearest.neighbour.has.homoplasy+homoplasic.counts$N.nearest.neighbour.lacks.homoplasy+homoplasic.counts$N.nearest.neighbour.intermediate
homoplasic.counts$proportion.nearest.neighbour.has.homoplasy <- homoplasic.counts$N.nearest.neighbour.has.homoplasy/homoplasic.counts$N.isolates.with.homoplasy

# Histogram of these distances
pdf('../figures/histogram-homoplasic-sites-nearest-neighbour-proportion-all-homoplasies.pdf')
hist(homoplasic.counts$proportion.nearest.neighbour.has.homoplasy, 
     breaks=100, 
     col='black', 
     xlab='Proportion of isolates where nearest neighbour in tree has homoplasy',
     main='')
dev.off()

# Add correct SNP counts (calculated in other script)
snp.counts.correct <- read.csv('data/SNP_counts_new.csv', header=T, row.names = 1)
homoplasic.counts$SNP_count <- snp.counts.correct$SNP.count[homoplasic.counts$bp]
# Write to file
write.csv(homoplasic.counts, file='../output-data/all-homoplasic-sites-table.csv')

# 4. Make plots 

# Further high-quality homoplasy thresholds
#DISTANCE.TO.HOMOPLASY <- 10
# High-quality homoplasy thresholds
#PROPORTION.NEAREST.NEIGHBOUR.HAS.HOMOPLASY <- 0.4
# Exclude the homoplasies with < 0.1 proportion (i.e. the large peak at 0)
N <- 500 
homoplasic.counts.filt <- homoplasic.counts[which(homoplasic.counts$bp>N & homoplasic.counts$bp<GENOME.LENGTH-N),]
pdf('../figures/histogram-homoplasic-sites-nearest-neighbour-proportion-excluding-first-last-500-bp.pdf')
hist(homoplasic.counts$proportion.nearest.neighbour.has.homoplasy, 
     breaks=100, 
     col='black', 
     xlab='Proportion of isolates where nearest neighbour in tree has homoplasy',
     main='')
dev.off()

NEAREST.HOMOPLASY.PROP <- 0.1
homoplasic.counts.filt.HQ <- homoplasic.counts.filt[which(homoplasic.counts.filt$proportion.nearest.neighbour.has.homoplasy>NEAREST.HOMOPLASY.PROP),]
# Number of isolates with homoplasy threshold
N.ISOLATES.WITH.HOMOPLASY <- 10
homoplasic.counts.filt.HQ <- homoplasic.counts.filt.HQ[which(  homoplasic.counts.filt.HQ$N.isolates.with.homoplasy>N.ISOLATES.WITH.HOMOPLASY),]
# No isolates with homoplasy with N in surround 4 bp
homoplasic.counts.filt.HQ <- homoplasic.counts.filt.HQ[which(  homoplasic.counts.filt.HQ$proportion.with.N.within.2.bp==0),]

write.csv(homoplasic.counts.filt.HQ, file='../output-data/filtered-homoplasic-sites-table.csv')


# For these, make plots
for (h in homoplasic.counts.filt.HQ$bp){
  print(h)
  h.df <- getCopheneticDistributionForHomoplasy(site = h)
  h.plot <- plotHomoplasyCopheneticDistribution(h.df, title=getTitleString(h))
  ggsave(h.plot, file=paste0('../figures/cophenetic-distributions/', h, '.pdf'))
}


