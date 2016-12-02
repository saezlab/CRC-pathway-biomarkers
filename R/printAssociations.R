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

#' @title Preprocess MIDAS files 
#'
#' @description
#' \code{processMIDAS} returns the process version for the row midas file for use in logic ODEs
#'
#' @details
#' This function:
#' - deals with the interpretation of drug PLX as BRAF inhibitor or stimulus depending if cell line is mutated or wild type
#' - scales the data between 0 and 1 (positive data are >0.5, negative are <0.5) whith basal value (initial condition) at 0.5 for all species
#' - add columns for species that are in the compressed model but not in the network with 0.5 at t0 and NA at t1 (just to set initial conditions)
#' - remove data for EGFRi (experiment not working properly)
#'
#' @param pkn The prior knowledge network imported using CellNOptR \code{\link{readSIF}}
#' @param rawMIDASpath The path to the raw MIDAS
#' @param processedMIDASpath The path where to save the processed MIDAS
#' @param cellLines Names of the cell lines to process
#' @return This fuction saves the processed MIDAS files and returns nothing

processMIDAS <- function(pkn, rawMIDASpath, processedMIDASpath, cellLines){
  
  for (myCellLine in cellLines){
    dataPath<-paste(rawMIDASpath, paste("MD-", myCellLine, "_Ktuned_v1_n4_raw.csv", sep=""), sep="")
    
    # model BRAF inhibitor (PLX) as inhbitior for BRAF mutants and as stimuli for BRAF wild type
    Mydata<-readMIDAS(dataPath)
    
    # this for cell lines whre BRAF is wild type, it was modeled as "no effect" now we model PLX drug as ligand
    ix1<-which(colnames(Mydata$dataMatrix)=="TR:PLXi")
    # this for cell lines whre BRAF is mutates, we keep modeling it as inhibited when PLX drug applied
    ix2<-which(colnames(Mydata$dataMatrix)=="TR:PLXmuti")
    # I now alway consider BRAFi
    if (length(ix2)>0){
      # this is for mutated cell line:
      # drug act as an inhibitor directly on the BRAF
      # I add a column for PLX stimulus but I leave it to 0
      colnames(Mydata$dataMatrix)[ix2]<-"TR:BRAFi"
      Mydata$dataMatrix<-cbind(rep(0, nrow(Mydata$dataMatrix)),
                               Mydata$dataMatrix)
      colnames(Mydata$dataMatrix)[1]<-"TR:PLX"
    }else if (length(ix1)>0){
      # this is for wild-type cell line:
      # drug act as stimuly and regulates BRAF
      # I add a column for PLX stimulus and leave to 0 the "TR:BRAFi" column
      colnames(Mydata$dataMatrix)[ix1]<-"TR:BRAFi"
      tmpBRAF<-Mydata$dataMatrix[,ix1]
      Mydata$dataMatrix[,ix1]<-0
      
      Mydata$dataMatrix<-cbind(tmpBRAF,
                               Mydata$dataMatrix)
      colnames(Mydata$dataMatrix)[1]<-"TR:PLX"
    }
    # generate the corresponding CNOlist
    Mydata$TRcol<-grep("TR:", colnames(Mydata$dataMatrix))
    Mydata$DAcol<-grep("DA:", colnames(Mydata$dataMatrix))
    Mydata$DVcol<-grep("DV:", colnames(Mydata$dataMatrix))
    MyCNOlist<-makeCNOlist(Mydata, subfield=F)
    
    # find which species will be compressed compressed to add them to the MIDAS (needed to set initial condition to 0.5)
    model<-preprocessing(data=MyCNOlist, model=pkn, compression=TRUE, expansion=FALSE)
    
    # read the MIDAS file as table
    MyMIDAS<-read.table(dataPath, header=T, sep=',', check.names=F, stringsAsFactors=F)
    MyMIDAS.DV<-MyMIDAS[,grep("DV:", colnames(MyMIDAS))]
    
    # 1. scale the t1 data between 0 and 1 (positive data are >0.5, negative are <0.5) and set the t0 data to 0.5 (can be done with same scaling)
    MyMIDAS.DV.scaled<-apply(MyMIDAS.DV, 2, function(x){
      # x<-(x-min(x))/(max(x)-min(x)) ## this was wrong!!
      x<-0.5/max(abs(x), na.rm = T)*x+0.5
      # x<-(x+1)*0.5
      return(x)
    })
    MyMIDAS.DV.scaled[(which(MyMIDAS[,grep("DA:", colnames(MyMIDAS))]==0)),]<-0.5 # to set to 0.5 also the NAs
    MyMIDAS.new<-MyMIDAS
    MyMIDAS.new[,grep("DV:", colnames(MyMIDAS))]<-MyMIDAS.DV.scaled
    
    # 2. add columns for species that are in the compressed model but not measured or stimulated (i.e. white and inhibited nodes)
    # with 0.5 at t0 and NA at t1
    measSpecies<-sub("DV:", "", colnames(MyMIDAS.DV.scaled))
    toAddSpecies<-setdiff(model$namesSpecies, c(measSpecies, MyCNOlist$namesStimuli))
    MyMIDAS.DV.toAdd<-matrix(NA, nrow = dim(MyMIDAS.DV.scaled)[1], ncol = length(toAddSpecies))
    MyMIDAS.DV.toAdd[(which(MyMIDAS[,grep("DA:", colnames(MyMIDAS))]==0)),]<-0.5
    colnames(MyMIDAS.DV.toAdd)<-paste("DV:", toAddSpecies, sep="")
    MyMIDAS.new<-cbind(MyMIDAS.new, MyMIDAS.DV.toAdd)
    
    # 3. remove EGFRi data
    # consider the inhibitors treatment columns (i.e. those that start with TR and finish with i)
    ix_TRi<-grep("^TR.*i$", colnames(MyMIDAS.new))
    # remove data with stimuli but no inhibior
    # (they are the one with EGFRi, that they thought was not working but it's actually run in different conditions)
    ix_remove<-which(apply(MyMIDAS.new[,ix_TRi], 1, sum)==0)
    # remove those lines and the column for EGFRi
    ix_EGFRi<-grep("EGFRi", colnames(MyMIDAS.new))
    MyMIDAS.new<-MyMIDAS.new[-ix_remove, -ix_EGFRi]
    
    # 4. add row for control experiment
    t0control<-rep(0, ncol(MyMIDAS.new))
    t0control[1]<-1 # cell line column
    t0control[grep("DV:", colnames(MyMIDAS.new))]<-0.5 # readouts
    t1control<-t0control
    t1control[grep("DA:", colnames(MyMIDAS.new))]<-30 # time point
    
    MyMIDAS.new<-rbind(t0control,
                       MyMIDAS.new[which(MyMIDAS.new[,grep("DA:", colnames(MyMIDAS.new))]==0),],
                       t1control,
                       MyMIDAS.new[which(MyMIDAS.new[,grep("DA:", colnames(MyMIDAS.new))]==30),])
    
    
    # this for cell lines whre BRAF is wild type, it was modeled as "no effect" now we model PLX drug as ligand
    ix1<-which(colnames(MyMIDAS.new)=="TR:PLXi")
    # this for cell lines whre BRAF is mutates, we keep modeling it as inhibited when PLX drug applied
    ix2<-which(colnames(MyMIDAS.new)=="TR:PLXmuti")
    
    
    if (length(ix2)>0){
      # this is for mutated cell line:
      # drug act as an inhibitor directly on the BRAF
      # I add a column for PLX stimulus but I leave it to 0
      colnames(MyMIDAS.new)[ix2]<-"TR:BRAFi"
      MyMIDAS.new<-cbind(MyMIDAS.new[,1, drop=F], # this is the cell line column
                         rep(0, nrow(MyMIDAS.new)),
                         MyMIDAS.new[,2:ncol(MyMIDAS.new)])
      colnames(MyMIDAS.new)[2]<-"TR:PLX"
    }else if (length(ix1)>0){
      # this is for wild-type cell line:
      # drug act as stimuly and regulates BRAF
      # I add a column for PLX stimulus and leave to 0 the "TR:BRAFi" column
      colnames(MyMIDAS.new)[ix1]<-"TR:BRAFi"
      tmpBRAF<-MyMIDAS.new[,ix1]
      MyMIDAS.new[,ix1]<-0
      
      MyMIDAS.new<-cbind(MyMIDAS.new[,1, drop=F], # this is the cell line column
                         tmpBRAF,
                         MyMIDAS.new[,2:ncol(MyMIDAS.new)])
      colnames(MyMIDAS.new)[2]<-"TR:PLX"
    }
    
    
    
    setwd(processedMIDASpath)
    fileName<-paste("MD-", myCellLine, "_Ktuned_v1_n4_all_noEGFRi_CNORode.csv", sep="")
    write.csv(MyMIDAS.new, file=paste(processedMIDASpath, fileName, sep=""), row.names=F)
    
  }
  
}
  
