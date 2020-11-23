
getMinorVariantFraction <- function(site){
  if (snp.counts[site,"Min.No.ChangesonTree"]==0){
    return(NA)
  }
  # Get allele counts
  allele.counts <- as.numeric(alt_alleles_at_site_mod[site,"N.isolates.with.homoplasy"])
  return(as.numeric(allele.counts/aln$nb))
}

getIsolatesWithMinorVariant <- function(site){
  # Get allele counts from the table providing the most represented alt allele
  minor.variant.base <- alt_alleles_at_site_mod[site,"ALT"]
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


# below = not used






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

