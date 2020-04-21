# Read in tree
tree <- read.tree('../data/gisaid_cov2020_sequences.QC.human.filter.aln.treefile')
cophenetic.d <- cophenetic(tree)
cophenetic.df <- as.data.frame(melt(cophenetic.d))

# Read in alignment
aln <- read.alignment('../data/gisaid_cov2020_sequences.QC.human.filter.rename.noambig.aln.fasta',
                      format = 'fasta')

# Read in SNP counts (tsv despite extension)
snp.counts <- read.csv('../data/SNP_homoplasy_counts_table.filter.csv',
                       header = T, 
                       stringsAsFactors = F,
                       sep = '\t')

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