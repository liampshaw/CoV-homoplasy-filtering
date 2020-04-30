# Makes supplementary figure 7 in the associated manuscript
source('functions.R')
homoplasic.counts.filt.HQ <- read.csv('../output-data/filtered-homoplasic-sites-table.csv',
                                      row.names = 1, 
                                      stringsAsFactors = F,
                                      header=T)
h <- sapply(homoplasic.counts.filt.HQ$bp, function(x) {print(x); getNsAtSite(x)})


# Also do for random sites in genome (excluding masked sites)
MASKED.REGIONS <- c(seq(1,150), seq(29853, 29903), 18529, 29849, 29851, 29853) 
random.sites <- floor(runif(min=151, max=29852, n=nrow(homoplasic.counts.filt.HQ)))
# Check none in masked regions
random.sites[which(random.sites %in% MASKED.REGIONS)]

random.dist <- sapply(random.sites, function(x) {print(x); getNsAtSite(x)})

par(mfrow=c(2,1))
hist(h, breaks=seq(0, 500, 10), col='black', xlab='Number of isolates with N at site', ylab='Proportion of homoplasic sites', main='a. Filtered homoplasic sites (n=198)', probability = TRUE, xlim=c(0,300), ylim=c(0,0.1))
hist(random.dist, col='red',  breaks=seq(0, 500, 10),  xlab='Number of isolates with N at site', ylab='Proportion of homoplasic sites', main='b. Random sites (excluding masked regions)', probability = TRUE, xlim=c(0,300), ylim=c(0,0.11))
