# Function for a given SNP. 
# Identify all isolates with the SNP
# Get cophenetic distances
getMinorVariantSite <- function(site){
  # Get allele counts
  allele.counts <- as.numeric(strsplit(split = ":", fixed=TRUE, x = snp.counts[which(snp.counts$bp==site),"CountsACGT"])[[1]])
  names(allele.counts) <- c("A", "C", "G", "T")
  
  # Sort by counts
  allele.counts <- sort(allele.counts, decreasing = TRUE)
  # Return second-highest
  return(tolower(names(allele.counts)[2]))
}
getMinorVariantFraction <- function(site){
  if (snp.counts[site,"Min.No.ChangesonTree"]==0){
    return(NA)
  }
  # Get allele counts
  allele.counts <- as.numeric(strsplit(split = ":", fixed=TRUE, x = snp.counts[which(snp.counts$bp==site),"CountsACGT"])[[1]])
  names(allele.counts) <- c("A", "C", "G", "T")
  
  # Sort by counts
  allele.counts <- sort(allele.counts, decreasing = TRUE)
  # Return second-highest
  return(as.numeric(allele.counts[2]/sum(allele.counts)))
}

getCopheneticDistributionForHomoplasy <- function(site, snp.counts=snp.counts){
  # Find isolates with minor variant
  minor.variant.base <- getMinorVariantSite(site)
  minor.variant <- unlist(lapply(aln$seq, function(x) substr(x, site, site)==minor.variant.base))
  isolates.with.minor.variant <- aln$nam[minor.variant]
  
  # Look at cophenetic distances
  # Need to first remove those not present in tree
  # Grep for EPI ISL
  #isolates.with.minor.variant.epi.isl <- gsub("\\|.*|", "", gsub(".*EPI", "EPI", isolates.with.minor.variant))
  tree.isolates.select <- tree$tip.label[unlist(sapply(isolates.with.minor.variant,
                                                       function(x) grep(x,tree$tip.label )))]
  
  
  homoplasy.cophenetic.df <- melt(cophenetic.d[tree.isolates.select, tree.isolates.select])
  return(homoplasy.cophenetic.df)
}
plotHomoplasyCopheneticDistribution <- function(homoplasy.df, title=""){
  ggplot(homoplasy.df, aes(value))+
    geom_histogram(binwidth=5e-5, stat="density", data=cophenetic.df.random, aes(value), fill="black", alpha=0.5)+
    geom_histogram(binwidth=5e-5, stat="density", fill='red')+
    geom_vline(xintercept = median(homoplasy.df$value), colour="red")+
    geom_vline(xintercept = median(cophenetic.d), colour="black")+
    geom_vline(xintercept = as.numeric(quantile(cophenetic.d, probs = 0.025)), colour="black", linetype='dashed')+
    geom_vline(xintercept = as.numeric(quantile(cophenetic.d, probs = 0.975)), colour="black", linetype='dashed')+
    xlab("cophenetic distance")+
    ylab("density")+
    ggtitle(title)+
    theme_basic()
}
getTitleString <- function(site, snp.df=homoplasic.counts.filt.HQ){
  # Useful information for a plot title
  consistency <- myround(snp.df[which(snp.df$bp==site), "consistency_index"], 3)
  protein <- snp.df[which(snp.df$bp==site), "Protein"]
  changes <- snp.df[which(snp.df$bp==site), "Min.No.ChangesonTree"]
  proportion.nearest.neighbour <- myround(snp.df[which(snp.df$bp==site), "proportion.nearest.neighbour.has.homoplasy"], 3)
  neighbour.has.homoplasy <- snp.df[which(snp.df$bp==site), "N.nearest.neighbour.has.homoplasy"]
  total.isolates <- snp.df[which(snp.df$bp==site), "N.nearest.neighbour.has.homoplasy"]+snp.df[which(snp.df$bp==site), "N.nearest.neighbour.lacks.homoplasy"]+snp.df[which(snp.df$bp==site), "N.nearest.neighbour.intermediate"]
  number.labs <- snp.df[which(snp.df$bp==site), "N.submitting.labs"]
  number.countries <- snp.df[which(snp.df$bp==site), "N.countries"]
  
  
  nearest.homoplasic.site <- snp.df[which(snp.df$bp==site), "dist.nearest.homoplasy"]
  
  return(paste0("Site ", site,
                ", Region: ", protein, 
                "\nInferred minimum changes on tree: ", changes,
                "\nConsistency: ", consistency, ", Nearest homoplasic site is ", nearest.homoplasic.site, "bp away",
                "\nProp. nearest neighbour has homoplasy: ", proportion.nearest.neighbour, " (", neighbour.has.homoplasy, "/", total.isolates, ")", 
                "\nSeen from ", number.labs, " originating labs across ", number.countries, " countries."))
}

# Nearest neighbour in tree for an isolate
nearestNeighboursIsolate <- function(isolate){
  # Get its corresponding unique isolate
  identical.isolates <- names(which(cophenetic.d[isolate,]==0))
  unique.isolate <- identical.isolates[which(identical.isolates %in% rownames(cophenetic.d.unique))]
  
  distances.to.others <- cophenetic.d.unique[unique.isolate,]
  distances.to.others <- distances.to.others[which(distances.to.others!=0)]
  nearest.neighbours <- names(distances.to.others)[which(distances.to.others==min(distances.to.others))]
  return(nearest.neighbours)
}
nearestNeighbourHasVariant <- function(isolate, site){
  nearest.neighbours <-nearestNeighboursIsolate(isolate)
  isolate.site <- substr(aln$seq[[which(aln$nam  ==isolate)]], site, site)
  
  if (length(nearest.neighbours)==1){
    same.variant <- substr(aln$seq[[which(aln$nam ==nearest.neighbours)]], site, site)==isolate.site
    if (same.variant){
      return(1)
    }
    else{
      return(0)
    }
    
  }
  else{
    nearest.neighbour.sites <- 0
    for (n in nearest.neighbours){
      nearest.neighbour.sites <-  nearest.neighbour.sites+as.numeric(substr(aln$seq[[which(aln$nam  ==n)]], site, site)==isolate.site)
    }
    return(nearest.neighbour.sites/length(nearest.neighbours))
  }
}



getIsolatesWithMinorVariant <- function(site){
  minor.variant.base <- getMinorVariantSite(site)
  minor.variant <- unlist(lapply(aln$seq, function(x) substr(x, site, site)==minor.variant.base))
  isolates.with.minor.variant <- aln$nam[minor.variant]
  
  # Look at cophenetic distances
  # Need to first remove those not present in tree
  # Grep for EPI ISL
  #isolates.with.minor.variant.epi.isl <- gsub("\\|.*|", "", gsub(".*EPI", "EPI", isolates.with.minor.variant))
  tree.isolates.select <- tree$tip.label[unlist(sapply(isolates.with.minor.variant,
                                                       function(x) grep(x,tree$tip.label )))]
  return(tree.isolates.select)
}

# Function to plot the metadata of isolates with a homoplasy
plotMetadataIsolates <- function(set.of.EPI.ids, metadata=isolate.metadata){
  subset.metadata <- metadata[which(metadata$gisaid_epi_isl %in% set.of.EPI.ids),]
  # Add a 'month' field
  subset.metadata$month <- gsub("^[0-9][0-9]\\/", "", subset.metadata$date)
  subset.metadata.plot.df <- subset.metadata %>% dplyr::group_by(month, region_exposure) %>%
    summarise(n=length(month))
  # Remove those with unknown month
  subset.metadata.plot.df <- subset.metadata.plot.df[which(subset.metadata.plot.df$month!="2020"),]
  #subset.metadata.plot.df$month.chron <- chron(dates=subset.metadata.plot.df$month,
        #                                       format=c(dates="m/y"))
  ggplot(subset.metadata.plot.df, aes(month, n, fill=region_exposure))+
    geom_bar(position="stack", stat="identity", width=0.5)+
    theme_basic()+
    xlab("date")
}

# Basic plot theme
theme_basic <- function () { 
  theme_bw(base_size=14) %+replace% 
    theme(
      axis.text=element_text(colour="black")
    ) %+replace% 
    theme(
      panel.grid=element_blank()
    )
}

# Count number of Ns adjacent to a site
getAdjacentNscore <- function(site, isolate, region=2, fasta=aln){
  bp.start <- site-region
  if (bp.start<0){ bp.start <- 0}
  bp.end <- site+region
  if (bp.start>29903){bp.end <- 29903}
  character.region <- strsplit(substr(aln$seq[[which(aln$nam==isolate)]], start=bp.start, stop=bp.end), split='')[[1]]
  character.region <- character.region[-round(length(character.region)/2)]
  return(as.numeric(length(grep("n", character.region, ignore.case = TRUE))))
}

# Nearest N
getDistanceNearestN <- function(site, isolate, fasta=aln){
  genome <- strsplit(aln$seq[[which(aln$nam==isolate)]], split='')[[1]]
  if (length(N_positions)==0){
    return(29903)
  }
  else{
      N_positions <- grep("n", genome, ignore.case = TRUE)
      return(min(abs(site-N_positions)))
  }
}


# Get SNP count
getSNPcount <- function(fasta, site){
  table.alleles <- sort(table(unlist(lapply(fasta$seq, function(x) substr(x, start=site, stop=site)))), decreasing = TRUE)
  
  # Only keep a/t/c/g, but count total proportion
  total.alleles <- sum(table.alleles)
  table.alleles <- table.alleles[which(names(table.alleles) %in% c("a", "t", "c", "g"))]
  # If table not all one allele, then return sum of minor variants (this is the SNP count for that site)
  if (length(table.alleles)>1){
    return(sum(table.alleles[2:length(table.alleles)]))
  }
  else{
    return(0)
  }
  
}

getNumberOfNsInGenome <- function(isolate, fasta=aln, removeKnownDodgyRegions=TRUE){
  genome <- strsplit(aln$seq[[which(aln$nam==isolate)]], split='')[[1]]
  if (removeKnownDodgyRegions==TRUE){
    dodgy_regions <- c(seq(1,150), seq(29853, 29903), 18529, 29849, 29851, 29853)
    positions <- seq(1, 29903)
    positions <- positions[which(!positions %in% dodgy_regions)]
    genome <- genome[positions]
    return(length(grep("n", genome, ignore.case = TRUE)))
  }
  else{
    return(length(grep("n", genome, ignore.case = TRUE)))
  }
}
getNsAtSite <- function(site, fasta=aln){
  bases <- sapply(aln$nam, function(x)
    substr(aln$seq[[which(aln$nam==x)]], site, site))
  return(length(grep("n", bases, ignore.case = TRUE)))
}
getNonATCGCount <- function(fasta, site){
  table.alleles <- sort(table(unlist(lapply(fasta$seq, function(x) substr(x, start=site, stop=site)))), decreasing = TRUE)

  # Only keep not a/t/c/g
  total.alleles <- sum(table.alleles)
  table.alleles <- table.alleles[which(!names(table.alleles) %in% c("a", "t", "c", "g"))]
  # If table not all one allele, then return sum of these
  if (length(table.alleles)>0){
    return(as.numeric(sum(table.alleles)))
  }
  else{
    return(0)
  }

}

