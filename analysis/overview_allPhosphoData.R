# packages
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

# R functions
source("/Users/eduati/CRC-pathway-biomarkers/R/dataFromMIDAS.R")

# data
load("/Users/eduati/CRC-pathway-biomarkers/data/base/metadata.RData")


cellLines=as.character(mapping_cellLines$name)
setwd("/Users/eduati/CRC-pathway-biomarkers/figs/")


#####
# # for the normalised data
# MIDASpath<-"/Users/eduati/CRC-pathway-biomarkers/data/processed/MIDAS/"
# allResults<-dataFromMIDAS(MIDASpath=MIDASpath, cellLines=cellLines, is.normalised = T)
# 
# # plot the heatmap
# pdf(paste("/Users/eduati/CRC-pathway-biomarkers/figs/", "dataHeatmap_normalised.pdf", sep=""),width=42,height=8)
# Reduce("+", allResults$allHeatmaps)
# dev.off()
# 
# pdf(paste("/Users/eduati/CRC-pathway-biomarkers/figs/", "dataHeatmap_normalised_dendo.pdf", sep=""),width=42,height=8)
# plot(allResults$dendogram)
# dev.off()


#####
# for the raw data
MIDASpath<-"/Users/eduati/CRC-pathway-biomarkers/data/base/MIDAS/"
allResults<-dataFromMIDAS(MIDASpath=MIDASpath, cellLines=cellLines, is.normalised = F)

# plot the heatmap
# pdf(paste("/Users/eduati/CRC-pathway-biomarkers/figs/", "dataHeatmap_raw.pdf", sep=""),width=42,height=8)
pdf(paste("/Users/eduati/CRC-pathway-biomarkers/figs/", "dataHeatmap_raw.pdf", sep=""),width=34,height=8)
Reduce("+", allResults$allHeatmaps)
dev.off()

pdf(paste("/Users/eduati/CRC-pathway-biomarkers/figs/", "dataHeatmap_raw_dendo.pdf", sep=""),width=42,height=8)
plot(allResults$dendogram)
dev.off()

col = colorRamp2(c(min(allResults$allData.m, na.rm = T), 0, max(allResults$allData.m, na.rm = T)), c("blue", "white", "red"))
pdf(paste("/Users/eduati/CRC-pathway-biomarkers/figs/", "dataHeatmap_raw_allInOne.pdf", sep=""),width=8,height=15)
allData.m.H<-allResults$allData.m
n_na<-apply(allData.m.H, 1, function(x){sum(is.na(x))})
allData.m.H<-allData.m.H[which(n_na==0),]
Heatmap(allData.m.H, col = col, na_col = "grey90",
        cluster_rows = FALSE, cluster_columns = T, show_heatmap_legend = T) 
dev.off()



######
# look at functional mutations
all_mutations<-read.table(file = "/Users/eduati/CRC-pathway-biomarkers/data/base/mutations.txt",  sep = "\t",header = T, stringsAsFactors = F)
all_mutations<-subset(all_mutations, !is.na(MUTATION_TARGET)) # limit to mutations with target in PKN or first neighbour
mapping_cellLines$id<-as.character(mapping_cellLines$id)
mapping_cellLines$name<-as.character(mapping_cellLines$name)
all_mutations$CELLLINE_NAME<-sapply(all_mutations$COSMIC_ID, function(x){
  mapping_cellLines$name[which(mapping_cellLines$id==x)]
})

mutations.t<-matrix(0, ncol=length(cellLines), nrow=length(unique(all_mutations$GENE_NAME)))
colnames(mutations.t)<-cellLines
rownames(mutations.t)<-sort(unique(all_mutations$GENE_NAME))

for (i in 1:nrow(all_mutations)){
  mutations.t[all_mutations$GENE_NAME[i], all_mutations$CELLLINE_NAME[i]]<-1
}
col.ord <- order.dendrogram(allResults$dendogram)
mutations.t<-mutations.t[,col.ord]
mutations.t<-rbind(mutations.t,
                   MSI=rep(0,ncol(mutations.t)))
mutations.t["MSI", c("CCK81", "HCT116", "SNUC2B", "SNUC5")]<-1

# save(mutations.t, file = "/Users/eduati/CRC-pathway-biomarkers/data/processed/mutations_table.RData")

pdf(paste("/Users/eduati/CRC-pathway-biomarkers/figs/", "dataHeatmap_raw_mutations.pdf", sep=""),width=38,height=6)
Heatmap(mutations.t, col = colorRamp2(c(0, 1), c("grey99", "grey20"), 0.6),
        cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
        rect_gp = gpar(col = "grey50", lty = 1, lwd = 0.5), name = "  ", column_title = "mutations") 
dev.off()

