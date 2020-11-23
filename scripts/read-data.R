# Read in tree
tree <- read.tree(TREE.FILE)

#tree <- drop.tip(tree, tip='402125') # drop outgroup (sometimes necessary)
# Read in alignment
aln <- read.alignment(ALIGNMENT.FILE,
                      format = 'fasta')

# Read in SNP counts 
snp.counts <- read.csv(HOMOPLASY.COUNTS.FILE,
                       header = T, 
                       stringsAsFactors = F, 
                       row.names = NULL, sep="\t")
snp.counts$Min.No.ChangesonTree <- snp.counts$MinimumNumberChangesOnTree # for consistency with rest of script
snp.counts$bp <- snp.counts$Position # for consistency with rest of script
snp.counts$bp <- as.character(snp.counts$bp)
rownames(snp.counts) <- as.character(snp.counts$bp)


# Add minor variant proportion
snp.counts$minor.variant.fraction <- sapply(as.character(rownames(snp.counts)), 
                                            function(x) getMinorVariantFraction(x))

# Read in metadata
metadata <- read.csv(METADATA.FILE, sep='\t', header=T, stringsAsFactors = F)



