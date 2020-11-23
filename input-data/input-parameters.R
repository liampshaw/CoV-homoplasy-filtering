# Length of the SARS-CoV-2 genome
GENOME.LENGTH <- 29903


# # RaxML tree file
# TREE.FILE <- '../input-data/raxML_input_tree.tree'
# # TREE.FILE <- '/home/bacterio/Downloads/26_07_2020_bundle/26_07_2020_bundle/annotatedNewickTree_31-07-20.tree'
# # Multiple sequence alignment file, with ambiguous sites removed and replaced with N
# ALIGNMENT.FILE <- '../input-data/input_alignment.aln'

input_path <- "../input-data/"

# RaxML tree file (not provided)
TREE.FILE <- paste(input_path,"raxML_input_tree.tree",sep="")
# Multiple sequence alignment file, with ambiguous sites removed and replaced with N (not provided)
ALIGNMENT.FILE <- paste(input_path,"input_alignment.aln",sep="")
# Homoplasy file from HomoplasyFinder output (provided)
HOMOPLASY.COUNTS.FILE <- paste(input_path,"input_homoplasy_finder.txt",sep="")
# Metadata file (not provided)
METADATA.FILE <- paste(input_path,"metadata.tsv",sep="")

# Output folders
FIGURE.OUTPUT.FOLDER <- '../figures'
DATA.OUTPUT.FOLDER <- '../output-data'

# Masked regions (suspect regions of alignment)
MASKED.REGIONS <- c(seq(1,150), seq(29853, 29903), 18529, 29849, 29851, 29853) 

# Filtering parameters
NEAREST.HOMOPLASY.PROP <- 0.1 # p_nn
N.ISOLATES.WITH.HOMOPLASY <- 46 # >0.1% of ~46,000 = 46
LOCAL.REGION <- 5 # +/- number of bp to check around homoplasy for Ns
N.SUBMITTING.LABS <- 1 # more than this many submitting labs in GISAID
