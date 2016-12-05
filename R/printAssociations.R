#  Copyright (c) 2016 - EMBL-EBI
#
#  File author(s): Federica Eduati (federica.eduati@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://github.com/saezlab/CRC-pathway-biomarkers
# --------------------------------------------------------

#' @title Print associations between model parameters and drug response
#'
#' @description
#' \code{printAssociations} show table and heatmap for strong associations
#'
#' @details
#' This function
#'
#' @param associations table as results of function \code{\link{ToDo}}
#' @param topN Number of associations to show
#' @param processedMIDASpath The path where to save the processed MIDAS
#' @param cellLines Names of the cell lines to process
#' @return This fuction saves the processed MIDAS files and returns nothing

printAssociations <- function(associations, topN=50, sortBy="effectSize"){
  
  associations<-associations[order(abs(associations[,sortBy]), decreasing = T),]
  associations<-associations[1:topN,]
  
  # prepare links names
  newLinkNames<-as.character(associations$link)
  ix_tau<-grep("tau_", newLinkNames)
  ix_k<-grep("_k_", newLinkNames)
  newLinkNames[ix_tau]<-sapply(newLinkNames[ix_tau], function(x){
    x_new<-gsub("tau", "\u03C4", x)
    x_new<-gsub("_", " ", x_new)
    return(x_new)
  })
  newLinkNames[ix_k]<-sapply(newLinkNames[ix_k], function(x){
    tmp<-strsplit(x, split = "_")[[1]]
    x_new<-paste(tmp[2], " ", tmp[1], ",", tmp[3], sep="")
    return(x_new)
  })
  
  # prepare names for genomic alterations
  newGenomicAlterations<-associations$mutationsTargetName
  newGenomicAlterations<-gsub("\\[", "", newGenomicAlterations)
  newGenomicAlterations<-gsub("\\]", "", newGenomicAlterations)
  newGenomicAlterations<-gsub(">!", "|", newGenomicAlterations)
  newGenomicAlterations<-gsub(";", "; ", newGenomicAlterations)
  
  # prepare drug targets
  newDrugTargets<-associations$drugTargets
  newDrugTargets<-gsub("MEK1;MEK2->ERK", "MEK", newDrugTargets)
  newDrugTargets<-gsub(";", "; ", newDrugTargets)
  
  
  newDF<-data.frame(associations$drug,
                    newLinkNames,
                    newDrugTargets,
                    newGenomicAlterations)
  
  colnames(newDF)<-c("Drug", "Associated \n model parameter", "Targets", "Associated \n genomic alterations")
  rownames(newDF)<-NULL
  
  
  
  cairo_pdf("/Users/eduati/CRC-pathway-biomarkers/figs/associationsTable.pdf", height=15, width=13)
  grid.table(newDF)
  # ss <- tableGrob(newDF)
  dev.off()
  
  TMP_df<-data.frame(association="association", values=associations$effectSize, names=associations$names)
  
  gg <- ggplot(TMP_df, aes(x=association, y=reorder(names, abs(values))))
  gg <- gg + geom_tile(aes(fill = values), color="white", size=0.1)
  gg <- gg + scale_fill_gradient2(low="darkblue", high="darkgreen", guide="colorbar")
  gg <- gg + coord_equal()
  gg <- gg + labs(x=NULL, y=NULL, title=NULL)  
  # gg <- gg + theme_tufte(base_family="Helvetica")
  gg <- gg + theme(axis.ticks=element_blank())
  ggsave(filename="/Users/eduati/CRC-pathway-biomarkers/figs/associationsHeatmap.pdf", height=15, width=13, plot=gg)
  
  
}
  
