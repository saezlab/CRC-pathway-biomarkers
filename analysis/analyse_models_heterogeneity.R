# packages
library(NMF)

# R functions
source("/Users/eduati/CRC-pathway-biomarkers/R/getAllOptResults_bootstrap.R")


# data
load("/Users/eduati/CRC-pathway-biomarkers/data/base/metadata.RData")
load("/Users/eduati/CRC-pathway-biomarkers/data/base/namesSignals.RData")


# bootstrap results
resultsFolder = "/Users/eduati/CRC-pathway-biomarkers/data/processed/optimisation/bootstrap/"
allOptResults<-getAllOptResults_bootstrap(resultsFolder=resultsFolder, namesSignals = namesSignals, optimisedPars=c("k", "tau"), N=150)
cellLinesNames<-names(allOptResults)

# compute median values across bootstrap runs
allParameters<-do.call(rbind, lapply(cellLinesNames, function(x){
  apply(allOptResults[[x]]$parameters, 2, median)
}))
rownames(allParameters)<-cellLinesNames
allParameters<-allParameters[,-which(colnames(allParameters)=="PLX_k_BRAF")] # not a real parameter, fixed at 0.5 to model PLX stimulation of BRAF in BRAF wt cell lines

# remove the parameters that are =0 across all cell lines
ix_0<-which(apply(allParameters,2,function(x){sum(x==0)})!=nrow(allParameters))
allParameters<-allParameters[,ix_0]
  
# compute dendogram
cc.col <- hclust(dist(allParameters))
dd.col <- as.dendrogram(cc.col)
col.ord <- order.dendrogram(dd.col)
cellLinesNames_ordered<-cellLinesNames[col.ord]

allParameters<-allParameters[cellLinesNames_ordered,]

# separate parameters k and tau for plot
ix_k<-grep("_k_", colnames(allParameters))
ix_tau<-grep("tau_", colnames(allParameters))

parameters_k<-allParameters[,ix_k]
parameters_tau<-allParameters[,ix_tau]

# plot heatmap for parameters k
# myTEXT<-round(myParameters_tau,2)
# myTEXT[myParameters_tau!=0]<-""
colnames(parameters_k)<-gsub("_k_", " -> ", colnames(parameters_k))
aheatmap(parameters_k, color = "Blues:30", Rowv=NA,
         txt = round(parameters_k,2), cellwidth = 20, cellheight = 20, fontsize = 7.5, cexRow=1.2, cexCol=1.2,
         filename="/Users/eduati/CRC-pathway-biomarkers/figs/parHeatmap_k.pdf")

# plot heatmap for parameters k
colnames(parameters_tau)<-gsub("tau_", "", colnames(parameters_tau))
aheatmap(parameters_tau, color = "Blues:30", Rowv=NA, 
         txt = round(parameters_tau,2), cellwidth = 20, cellheight = 20, fontsize = 7.5, cexRow=1.2, cexCol=1.2,
         filename="/Users/eduati/CRC-pathway-biomarkers/figs/parHeatmap_tau.pdf")

# load mutation table file (generated in analysis/overview_allPhosphoData.R)
load("/Users/eduati/CRC-pathway-biomarkers/data/processed/mutations_table.RData")
mutations.t<-t(mutations.t)
mutations.t<-mutations.t[cellLinesNames_ordered,]
aheatmap(mutations.t, Colv=NA, Rowv=NA, cellwidth = 10, cellheight = 20, fontsize = 7.5, cexRow=1.7, cexCol=1.2, border_color="white",
         filename="/Users/eduati/CRC-pathway-biomarkers/figs/parHeatmap_mut.pdf")


# plot dendogram
pdf(paste("/Users/eduati/CRC-pathway-biomarkers/figs/", "parHeatmap_dendo.pdf", sep=""),width=2,height=3)
plot(dd.col, horiz = TRUE)
dev.off()
