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

#' @title ROC analysis
#'
#' @description
#' \code{getAllROC} plot ROC curves for model predictions
#'
#' @details
#' For computing ROC curves, comparison of observed and predicted values is considered as a multi-class classification problem
#' which are up (increase with respect to basal), down (decrease with respect to basal), nc (no chance with respect to basal)
#'
#' @param resultsFolder Path of the folder with the results of the optimisation (form the cluster)
#' @param namesSignals Name of the readouts (to compute scores)
#' @param posClass Which class to consider as positive in the one-vs-all comparison: can be "up", "down", "nc" (no change) or a combination of those
#' @param th_vec Numeric vector of thresholds to consider. Must be between 0 and 0.5, default is seq(0, 0.5, 0.001)
#' 
#' @return This fuction plots the ROC curve for each class and returs TPR and FPR

getAllROC <- function(resultsFolder, namesSignals, posClass="up", th_vec=seq(0, 0.5, 0.001)){
  currentFolder<-getwd()
  
  setwd(resultsFolder)
  
  AllFiles<-list.files(path = ".", pattern = NULL, all.files = FALSE,
                       full.names = FALSE, recursive = FALSE,
                       ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  currentPath<-resultsFolder
  
  allByCellLine<-lapply(AllFiles, function(cellLineFolder){
    setwd(paste(currentPath, cellLineFolder, sep=""))
    AllRunsFile<-list.files(path = ".", pattern = NULL, all.files = FALSE,
                            full.names = FALSE, recursive = FALSE,
                            ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
    AllRunsFile<-AllRunsFile[grep(".RData", AllRunsFile)]
    AllRunsFile<-setdiff(AllRunsFile, "eSSR_report.RData")
    
    allRunsROC<-lapply(AllRunsFile, function(runFile){  
      
      res.df<-tryCatch({
        load(runFile)
        ix_signals<-which(cnolist$namesSignals %in% namesSignals)
        
        
        if (any(simulatedData[[1]][,ix_signals] != 0.5) | any(cnolist$valueSignals[[1]][,ix_signals] != 0.5)){cat("something wrong t0 should always be 0.5")}
        
        y_meas<-as.vector(cnolist$valueSignals[[2]][,ix_signals])
        y_pred<-as.vector(simulatedData[[2]][,ix_signals])
        # remove NAs
        ix_na<-(is.na(y_meas) | is.na(y_pred))
        y_meas<-y_meas[ix_na==F]
        y_pred<-y_pred[ix_na==F]
        
        resROC<-computeROC(y_meas, y_pred, posClass=posClass, th_vec=th_vec)
        
        return(resROC)
      },
      error=function(cond){
        message(cond)
        # cat("pippo")
        cat("\n not able to load", runFile)
        
        return(NULL)
      })
      
      return(res.df)
      
    })
    
    TPR<-apply(do.call(cbind, lapply(allRunsROC, function(myClassROC){myClassROC$TPR})),1,mean)
    FPR<-apply(do.call(cbind, lapply(allRunsROC, function(myClassROC){myClassROC$FPR})),1,mean)
    
    th.df<-allRunsROC[[1]]$th
    class.df<-allRunsROC[[1]]$class
    # # this was just to double check that they have all the same order, so that I can do the mean as I did
    # sapply(allRunsROC, function(x){all(x$th==th.df)})
    # sapply(allRunsROC, function(x){all(x$class==class.df)})
    
    
    myROCstatistics<-data.frame(TPR=TPR, FPR=FPR, th=th.df, class=class.df)
    
    return(myROCstatistics)
  })
  
  # retrieve cell line from fileName
  cellLinesNames<-sapply(AllFiles, function(x) {
    strsplit(x, split="_")[[1]][1]
  })
  cellLinesNames<-sapply(cellLinesNames, function(x) {
    strsplit(x, split="-")[[1]][2]
  })
  
  names(allByCellLine)<-cellLinesNames
  
  resROCallCellLines<-do.call(rbind, lapply(cellLinesNames, function(x){
    tmp<-allByCellLine[[x]]
    tmp$cellLine=x
    return(tmp)
  }))
  
  setwd(currentFolder)
  
  ggplot(resROCallCellLines, aes(x=FPR, y=TPR, colour = cellLine)) + geom_line() + geom_point(shape=".", size=0.0000001) + geom_abline(intercept = 0, slope = 1) + theme_bw() + facet_grid(class ~ .)
  ggsave("/Users/eduati/CRC-pathway-biomarkers/figs/ROCcurves.pdf", width=6, height=11)
  
  return(resROCallCellLines)
}


computeROC <-function(y_meas, y_pred, posClass="up", th_vec=seq(0, 0.5, 0.001)){

    ROCdata<-do.call(rbind, lapply(th_vec, function(th){
    conMat<-computeConfusionMatrix(y_meas, y_pred, th=th, posClass=posClass)
    
    res.df<-do.call(rbind, lapply(names(conMat), function(myConMatClass){
      TPR=conMat[[myConMatClass]][["TP"]]/(conMat[[myConMatClass]][["TP"]]+conMat[[myConMatClass]][["FN"]])
      FPR=conMat[[myConMatClass]][["FP"]]/(conMat[[myConMatClass]][["FP"]]+conMat[[myConMatClass]][["TN"]])
      return(data.frame(TPR=TPR, FPR=FPR, class=myConMatClass))
    }))
    
    res.df$th<-th
    return(res.df)
  }))
  
  return(ROCdata)
}


computeConfusionMatrix <- function(y_meas, y_pred, th=0.05, posClass="up"){
  # compute ROC analysis (as in Vicky's slides)
  # define which data go up (>0.55), down (<0.45) and no chage [0.45-0.55] in both measured and predicted values
  
  th_ROC_l<-0.5-th  #0.45
  th_ROC_u<-0.5+th  #0.55
  
  meas<-list(up=which(y_meas > th_ROC_u),
             down=which(y_meas < th_ROC_l),
             noc=which(y_meas >= th_ROC_l & y_meas <= th_ROC_u))
  
  pred<-list(up=which(y_pred > th_ROC_u),
             down=which(y_pred < th_ROC_l),
             noc=which(y_pred >= th_ROC_l & y_pred <= th_ROC_u))
  
  # calculate confusion matrix
  confMatrix<-matrix(NA, ncol=3, nrow=3)
  allClasses<-names(meas)
  pred<-pred[names(meas)] # just to make sure they have the same order when computing confusion matrix
  colnames(confMatrix)<-allClasses # columns are predicted (optimized) data
  rownames(confMatrix)<-allClasses # rows are measured (experimental) data
  
  
  for (i in 1:nrow(confMatrix)){
    for (j in 1:ncol(confMatrix)){
      confMatrix[i,j]<-length(intersect(meas[[i]], pred[[j]]))   #rows are measured (experimental) data, columns are predicted (optimized) data
    }
  }
  
  # compute TP, TN, FN, FP for each class
  res<-lapply(posClass, function(myClass){
    TP<-confMatrix[myClass, myClass]
    FN<-sum(confMatrix[myClass, setdiff(allClasses, myClass)]) # total FN for the class is the sum of values in corresponding row, excluding the TP
    FP<-sum(confMatrix[setdiff(allClasses, myClass), myClass]) # total FP for the class is the sum of values in corresponding column, excluding the TP
    TN<-sum(confMatrix[setdiff(allClasses, myClass), setdiff(allClasses, myClass)]) #total TN for the class is the sum of all columns and rows excluding the class's columns and rows
    return(c(TP=TP, TN=TN, FN=FN, FP=FP))
  })
  names(res)<-posClass
  
  return(res)
}