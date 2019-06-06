#! /usr/bin/R

# Use this to convert an .rda file to a .csv file.
# Output file will be in the same place as the input file but with ".rda" replaced with ".csv".
# 
# Author: sshorser
###############################################################################


argv <- commandArgs(TRUE)
inFile <- toString(argv[1])
print(inFile)

outFile <- gsub(".rda$", ".csv", inFile)
print(outFile)

inData <- load(inFile)
write.csv(get(inData), file=outFile)