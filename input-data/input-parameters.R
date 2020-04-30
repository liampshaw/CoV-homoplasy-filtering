# Length of the SARS-CoV-2 genome
GENOME.LENGTH <- 29903
# RaxML tree file
TREE.FILE <- '../input-data/raxML_input_tree.tree'
# Multiple sequence alignment file, with ambiguous sites removed and replaced with N
ALIGNMENT.FILE <- '../input-data/input_alignment.aln'
# Homoplasy file (from HomoplasyFinder output)
HOMOPLASY.COUNTS.FILE <- '../input-data/input_homoplasy_finder.csv'
# Metadata file
METADATA.FILE <- '../input-data/metadata.tsv'


# Output folders
FIGURE.OUTPUT.FOLDER <- '../figures'
DATA.OUTPUT.FOLDER <- '../output-data'

# Masked regions (suspect regions of alignment)
MASKED.REGIONS <- c(seq(1,150), seq(29853, 29903), 18529, 29849, 29851, 29853) 

# Filtering parameters
NEAREST.HOMOPLASY.PROP <- 0.1 # p_nn
N.ISOLATES.WITH.HOMOPLASY <- 8 # >0.1% of 7666 = 7.7, rounded up to 8
LOCAL.REGION <- 5 # +/- number of bp to check around homoplasy for Ns