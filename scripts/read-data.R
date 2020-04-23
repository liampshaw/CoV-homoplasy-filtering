# Read in tree
tree <- read.tree(TREE.FILE)
#tree <- drop.tip(tree, tip='402125') # drop outgroup (sometimes necessary)
# Make cophenetic distances
cophenetic.d <- cophenetic(tree)
cophenetic.df <- as.data.frame(melt(cophenetic.d))
# Read in alignment
aln <- read.alignment(ALIGNMENT.FILE,
                      format = 'fasta')

# Read in SNP counts 
snp.counts <- read.csv(HOMOPLASY.COUNTS.FILE,
                       header = T, 
                       stringsAsFactors = F, 
                       row.names = 1)
rownames(snp.counts) <- as.character(snp.counts$bp)
#snp.counts$Min.No.ChangesonTree <- snp.counts$MinimumNumberChangesOnTree # for consistency with rest of script



# Add minor variant proportion
snp.counts$minor.variant.fraction <- sapply(as.character(snp.counts$bp), 
                                            function(x) getMinorVariantFraction(x))

# Random subset of cophenetic distances (for plotting distribution, otherwise
# takes a long time for each homoplasy plot
# NB if dataset is smaller this may throw an error
cophenetic.df.random <- cophenetic.df[sample(nrow(cophenetic.df),size = 100000),]
