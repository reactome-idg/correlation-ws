# Perform normalization (using normalize.quantiles) on the expression data, save to new HDF
# 
# Author: sshorser
###############################################################################

# Set up packages
packages <- c("rhdf5", "preprocessCore")
if (length(setdiff(packages, rownames(installed.packages()))) > 0)
{
	print("Install required packages")
	source("https://bioconductor.org/biocLite.R")
	biocLite("rhdf5")
	biocLite("preprocessCore")
}

library("rhdf5")
library("preprocessCore")

# args from CLI
args <- commandArgs(TRUE)

print(args)

if (length(args) < 2)
{
  print("Two arguments are required: 1) path to source file 2) path to output file")
  stopifnot(length(args) < 2)
}


source_file <- args[1]
# source_file <- '/tmp/human_matrix.h5'

output_file <- args[2]
#output_file <- '/tmp/normalized_human_matrix.h5'

# Get the input RDA file.
if(!file.exists(source_file))
{
	print(paste("Downloading compressed gene expression matrix to ", source_file))
	url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.rda"
	download.file(url, source_file, quiet = FALSE)
} else {
	print("Local file already exists.")
}
print(paste("Loading ", source_file))
# load the RDA file
load(source_file)
print("Loading is complete!")
class(meta)
class(expression)
# read in the list of sample IDs
#samples <- h5read(source_file, "meta/Sample_geo_accession")
samples <- meta["Sample_geo_accession"]
class(samples)
genes <- meta["genes"]
class(genes)
# read in the list of gene symbols
#genes <- h5read(source_file, "meta/genes")

# Need to loop through samples, normalizing them, and then write a new HDF file.
length(genes)
#for (i in 1:1)
#{
	
#	print(samples[i])
#	sample_location <- which(samples %in% sample(samples, 100, replace=TRUE))
	#sample_location <- which(samples %in% c(samples[i], samples[i+1]))
	#sample_location <- which(samples %in% samples[i])
#	index4h5 <- list(1:length(genes), sample_location)
#	expression <- h5read(file=source_file, name="data/expression", index=index4h5)
#	expression <- log2(expression+1)
	normalized <- normalize.quantiles(expression[,1:3])

#}
# cleanup before attempting to create new file.
if(file.exists(output_file))
{
  print("Removing pre-existing output file.")
  file.remove(output_file)
}

h5createFile(output_file)
h5createGroup(output_file,"data")
h5createDataset(output_file, "data/normalized_expression", dim(normalized), storage.mode=storage.mode(normalized))
h5write(normalized, file=output_file, "data/normalized_expression")
h5closeAll()