# packages
library(CellNOptR)

# R functions
source("/Users/eduati/CRC-pathway-biomarkers/R/processMIDAS.R")

# data
pkn<-readSIF("/Users/eduati/CRC-pathway-biomarkers/data/base/PKN.sif")
rawMIDASpath<-"/Users/eduati/CRC-pathway-biomarkers/data/base/MIDAS/"
load("/Users/eduati/CRC-pathway-biomarkers/data/base/metadata.RData")
processedMIDASpath<-"/Users/eduati/CRC-pathway-biomarkers/data/processed/MIDAS/"

cellLines=mapping_cellLines$name
processMIDAS(pkn=pkn, rawMIDASpath, processedMIDASpath=processedMIDASpath, cellLines=cellLines)
