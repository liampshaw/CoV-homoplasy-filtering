# Read in tree
tree <- read.tree('../data/gisaid_cov2020_sequences.QC.human.filter.aln.raxml.ROOT.tree')
tree <- drop.tip(tree, tip='hCoV-19/Wuhan-Hu-1/2019|EPI_ISL_402125') # drop outgroup
# rename tree tips to use GISAID EPI ID
tree$tip.label <- gsub("\\|.*", "", gsub(".*\\|EPI", "EPI", tree$tip.label))

# Read in alignment
aln <- read.alignment('../data/gisaid_cov2020_sequences.QC.human.filter.rename.noambig.aln.fasta',
                      format = 'fasta')
# rename aln names to use GISAID EPI ID
aln$nam <- gsub("_2020-.*", "", gsub(".*_EPI", "EPI", aln$nam))
aln$nam <- gsub("_2019-.*", "", aln$nam)


# Only keep the overlap
overlap.isolates <- tree$tip.label[which(tree$tip.label %in% aln$nam )]
# Prune tree
tree <- drop.tip(tree, tip=tree$tip.label[which(!tree$tip.label %in% overlap.isolates)])
cophenetic.d <- cophenetic(tree)
cophenetic.df <- as.data.frame(melt(cophenetic.d))
# Prune alignment
aln.isolates.to.keep <- which(aln$nam %in% overlap.isolates)
prune.aln.seqs <- list()
for (i in seq(1, length(aln.isolates.to.keep))){
  isolate.index <- aln.isolates.to.keep[i]
  prune.aln.seqs[[i]] <- aln$seq[[isolate.index]]
}
write.fasta(sequences = prune.aln.seqs, names = aln$nam[aln.isolates.to.keep],
            file='../data/pruned-alignment.aln')
#aln <- read.alignment('../data/pruned-alignment.aln', format='fasta')

# Read in SNP counts (tsv despite extension)
snp.counts <- read.csv('../data/SNP_homoplasy_counts_table.filter.csv',
                       header = T, 
                       stringsAsFactors = F,
                       sep = '\t')

# Regenerate SNP counts from alignment
# Get all SNP counts
counts <- c()
for (i in 1:GENOME.LENGTH){
  print(i)
  counts <- c(counts, getSNPcount(aln, i))
}
snp.counts$SNP <- counts


# Read in isolate metadata
isolate.metadata <- read.csv('../data/metadata.tsv',
                             header=T, 
                             stringsAsFactors = F,
                             sep = '\t')

# Add minor variant proportion
snp.counts$minor.variant.fraction <- sapply(snp.counts$bp, 
                                            function(x) getMinorVariantFraction(x))

# Random subset of cophenetic distances (for plotting distribution)
cophenetic.df.random <- cophenetic.df[sample(nrow(cophenetic.df),size = 100000),]
