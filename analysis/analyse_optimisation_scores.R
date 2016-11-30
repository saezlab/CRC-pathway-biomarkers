# packages
library(ggplot2)

# R functions
source("/Users/eduati/CRC-pathway-biomarkers/R/getAllOptResults_bootstrap.R")
source("/Users/eduati/CRC-pathway-biomarkers/R/getAllROC.R")



# data
load("/Users/eduati/CRC-pathway-biomarkers/data/base/metadata.RData")
load("/Users/eduati/CRC-pathway-biomarkers/data/base/namesSignals.RData")



######
# L1 regularisation on tau

######
# L1 regularisation k

######
# bootstrap
resultsFolder = "/Users/eduati/CRC-pathway-biomarkers/data/processed/optimisation/bootstrap/"
allOptResults<-getAllOptResults_bootstrap(resultsFolder=resultsFolder, namesSignals = namesSignals, optimisedPars=c("k", "tau"), N=150)
cellLinesNames<-names(allOptResults)

# take a look at the scores
allScores<-do.call(rbind, lapply(cellLinesNames, function(x){
  tmp_data<-allOptResults[[x]]$scores
  stat_means<-apply(tmp_data,2,mean)
  stat_sd<-apply(tmp_data,2,sd)
  data.frame(cellName=x, value=stat_means, sd=stat_sd, label=colnames(tmp_data))
}))
allScores_plot<-subset(allScores, label %in% c("MSE", "COD", "corrCoeff"))

# plot all scores
ggplot(allScores_plot, aes(x=cellName, y=value)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd),
                width=.2,                    
                position=position_dodge(.9)) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_grid(label ~ ., scales="free_y")
ggsave("/Users/eduati/CRC-pathway-biomarkers/figs/allScoresBarplot.pdf", width=6, height=4)

# plot total sum of squares (SS_tot) vs residual sum of squares (SS_res)
SS_tot_all<-subset(allScores, label=="SS_tot")
SS_res_all<-subset(allScores, label=="SS_res")
pdf("/Users/eduati/CRC-pathway-biomarkers/figs/SStot_vs_SSres.pdf",width=6,height=6)
plot(SS_tot_all$value, SS_res_all$value, xlim=c(0,30), ylim=c(0,30), pch=16, xlab="total sum of squares (SS_tot)", ylab="residual sum of squares (SS_res)")
text(SS_tot_all$value, SS_res_all$value, labels = SS_tot_all$cellName, pos = 1)
abline(a=0, b=1)
dev.off()

# plot only MSE (for Figure 2)
MSE_barplot<-subset(allScores, label=="MSE")
ggplot(MSE_barplot, aes(x=cellName, y=value)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd),
                width=.2,
                position=position_dodge(.9)) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="cell line", y="MSE (mean squared error)")
ggsave("/Users/eduati/CRC-pathway-biomarkers/figs/MSE.pdf", width=8, height=4)

# compute ROC curves (can be slow depending on th_vec)
resROCallCellLines<-getAllROC(resultsFolder, namesSignals, posClass=c("up", "down", "noc"), th_vec=seq(0, 0.5, 0.001))
