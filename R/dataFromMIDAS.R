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

#' @title Format MIDAS data as matrix 
#'
#' @description
#' \code{dataFromMIDAS} reads the data (treated, not control) from MIDAS and reformat them
#'
#' @details
#' This function extract all data from MIDAS files and format it as matrix with cell lines on the columns
#' and conditions on the rows
#'
#' @param MIDASpath The path to the MIDAS
#' @param cellLines Names of the cell lines to process
#' @param is.normalised True if processed data, False if raw data
#' @return This fuction return: a matrix with all data (at time 30) (rows) for each cell line (column),
#' the plot of the histograms by cell line, dendogram of clustered cells

dataFromMIDAS <- function(MIDASpath, cellLines, is.normalised = NA){
  
  # list the names of the midas files in the selected folder
  allMIDASFiles<-list.files(path = MIDASpath, pattern = NULL, all.files = FALSE,
                            full.names = FALSE, recursive = FALSE,
                            ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  # plot the perturbations
  tmp_data<-read.csv(paste(MIDASpath, allMIDASFiles[[7]], sep=""))
  ix_30<-which(tmp_data[,grep("DA", colnames(tmp_data))]==30)
  ix_TR<-grep("TR\\.", colnames(tmp_data))
  ix_TR<-ix_TR[-c(1)] # first one is just cell line name, secod is the one used to model BRAF activation
  
  tmp_data<-as.matrix(tmp_data[ix_30, ix_TR])
  colnames(tmp_data)<-gsub(colnames(tmp_data), pattern = "TR.", replacement = "")
  rownames(tmp_data)=rep(" ", nrow(tmp_data))
  cueHeatmap<-Heatmap(tmp_data, col = colorRamp2(c(0, 1), c("grey99", "grey20"), 0.6),
                      cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
                      rect_gp = gpar(col = "grey50", lty = 1, lwd = 0.5), name = "  ", column_title = "perturbations") 
  
  
  if (is.normalised == T){
    col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  }else if(is.normalised == F){
    allDataTmp<-do.call(cbind, lapply(cellLines, function(myCellLine){
      dataPath<-paste(MIDASpath, allMIDASFiles[grep(myCellLine, allMIDASFiles)], sep="")
      cellLineData<-read.csv(dataPath)
      ix_30<-which(cellLineData[,grep("DA", colnames(cellLineData))]==30)
      ix_DV<-grep("DV", colnames(cellLineData))
      theData.v_tmp<-as.vector(as.matrix(cellLineData[ix_30, ix_DV]))
      return(theData.v_tmp)
    }))
    
    col = colorRamp2(c(min(allDataTmp, na.rm = T), 0, max(allDataTmp, na.rm = T)), c("blue", "white", "red"))
  }else{
    cat("ERROR: is.normalised must be either true or false")
  }
  
  #extract, for each cell line all the data at time 30 and return them as a vector and as a plotted heatmap
  allResults<-lapply(cellLines, function(myCellLine){
    # cat(myCellLine, "\n")
    dataPath<-paste(MIDASpath, allMIDASFiles[grep(myCellLine, allMIDASFiles)], sep="")
    cellLineData<-read.csv(dataPath)
    ix_30<-which(cellLineData[,grep("DA", colnames(cellLineData))]==30)
    ix_DV<-grep("DV", colnames(cellLineData))
    theData.m<-as.matrix(cellLineData[ix_30, ix_DV])
    tmp_na<-apply(theData.m,2,function(x){sum(is.na(x))==(length(x)-1)}) # for the transoformed midas remove the artificially addecolorRamp2(c(-3, 0, 3), c("green", "white", "red")d species (which are NA)
    if (any(tmp_na)==T){
      # print(colnames(theData.m)[which(tmp_na==T)])
      theData.m<-theData.m[,which(tmp_na==F)]
    }
  
    colnames(theData.m)<-gsub(colnames(theData.m), pattern = "DV.", replacement = "")
    rownames(theData.m)=rep("", nrow(theData.m))
    theHeatmap<-Heatmap(theData.m, col = col, na_col = "grey90",
                        cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
                        rect_gp = gpar(col = "white", lty = 1, lwd = 0.01), name = myCellLine, column_title = myCellLine) 
    
    theData.v<-as.vector(theData.m)
    return(list(theHeatmap=theHeatmap, theData.v=theData.v))
  })
  
  # retrieve all the heatmaps 
  allHeatmaps<-lapply(allResults, function(x){x$theHeatmap})
  # reformat the data as a matrix (with one column for each cell line)
  allData.m<-do.call(cbind, lapply(allResults, function(x){x$theData.v}))
  colnames(allData.m)<-cellLines
  # remove NA (NA are the same across all cell lines)
  n_na<-apply(allData.m, 1, function(x){sum(is.na(x))})
  allData.m_noNA<-allData.m[which(n_na==0),]
  
  
  cc.col <- hclust(dist(t(allData.m_noNA)))
  dd.col <- as.dendrogram(cc.col)
  col.ord <- order.dendrogram(dd.col)
  
  # pdf(paste("/Users/eduati/CRC-pathway-biomarkers/figs/dendo_", fileName, sep=""),width=42,height=8)
  # plot(dd.col)
  # dev.off()
  
  allHeatmaps<-allHeatmaps[col.ord]
  
  allHeatmaps<-c(cueHeatmap, allHeatmaps)
  
  # plot all the data
  # if (!is.na(fileName)){
  #   pdf(paste("/Users/eduati/CRC-pathway-biomarkers/figs/", fileName, sep=""),width=42,height=8)
  #   Reduce("+", allHeatmaps)
  #   dev.off()
  # }
  
  return(list(allData.m=allData.m, allHeatmaps=allHeatmaps, dendogram=dd.col))
}
  
