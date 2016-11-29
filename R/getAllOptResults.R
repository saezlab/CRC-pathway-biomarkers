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

#' @title retrieve bootstrap optimisation results for all cell lines
#'
#' @description
#' \code{getAllOptResults_bootstrap} return the optimised parameters and some performance scores
#'
#' @details
#' This function reads the file from the optimisation run on the cluster and returns optimised
#' parameters and scores for each run and each cell line (e.g. bootstrap was run ~150 times for each cell line)
#'
#' @param resultsFolder Path of the folder with the results of the optimisation (form the cluster)
#' @param namesSignals Name of the readouts (to compute scores)
#' @param optimisedPars Parameters estimated: can be "k", "tau", "n" or a combination of those
#' @param N Optional number of runs to consider (to have the same number of runs for all cell lines)
#' 
#' @return This fuction returns a list one element for each cell line with "scores" and "parameters".
#'    Scores are organised as data frame with one run for each row and one metric for each column. Metrics are:
#'    - MSE: mean squared error
#'    - SS_tot: total sum of squares
#'    - SS_res: residual sum of squares
#'    - COD: coefficient of determination (or R^2)
#'    - corrCoeff: Pearson correlation
#'    - corrPval: p-value for Pearson correlation
#'    Parameters are organised as matrix with parameters on the columns and runs on the rows

getAllOptResults_bootstrap <- function(resultsFolder, namesSignals, optimisedPars, N="all"){
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
    
    #####
    #### compute null hypothesis for COD for each cell line (randomising data)
    #   load(AllRunsFile[[1]])
    #   ix_signals<-which(cnolist$namesSignals %in% namesSignals)
    #   
    #   y_meas<-as.vector(cnolist$valueSignals[[2]][,ix_signals])
    #   y_meas<-y_meas[!is.na(y_meas)]
    #   
    #   # compute total sum of squares (SS_tot)
    #   y_meas_mean<-mean(y_meas)
    #   SS_tot<-sum((y_meas-y_meas_mean)^2)
    #   
    #   # compute SS_res and COD for randomised predictions
    #   COD_random<-sapply(seq(1:10000), function(x){
    #     
    #     # generate random predictions from permutated data
    #     y_pred<-sample(y_meas, replace = F)
    #     
    #     # compute residual sum of squares
    #     SS_res<-sum((y_meas-y_pred)^2)
    #     
    #     # compute the coefficient of determination (COD or R2)
    #     COD=1-SS_res/SS_tot
    #     
    #     # return(COD)
    #     return(COD)
    #   })
    #   #########
    # # })
    
    allRunsScores<-do.call(rbind, lapply(AllRunsFile, function(runFile){
      
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
        
        # compute MSE
        MSE<-mean((y_meas-y_pred)^2)
        
        # compute total sum of squares (SS_tot)
        y_meas_mean<-mean(y_meas)
        SS_tot<-sum((y_meas-y_meas_mean)^2)
        
        # compute residual sum of squares
        SS_res<-sum((y_meas-y_pred)^2)
        
        # compute the coefficient of determination (COD or R2)
        COD=1-SS_res/SS_tot
        # COD_pval<-length(which(COD_random>-COD))/length(COD_random) # uncomment this to have the p value
        
        # compute r (from correlation)
        corrCoeff<-cor.test(y_meas, y_pred)$estimate
        corrPval<-cor.test(y_meas, y_pred)$p.value
        
        # res.df<-data.frame(MSE=MSE, SS_tot=SS_tot, SS_res=SS_res, COD=COD, COD_pval=COD_pval, corrCoeff=corrCoeff, corrPval=corrPval)  # this also shows pvalues for COD
        res.df<-data.frame(MSE=MSE, SS_tot=SS_tot, SS_res=SS_res, COD=COD, corrCoeff=corrCoeff, corrPval=corrPval)
        
        if (sum(is.na(simulatedData[[2]][,ix_signals]))>0){
          cat("\n---------\nThere are NAs in the simulation, check\n---------\n")
        }
        return(res.df)
      },
      error=function(cond){
        message(cond)
        # cat("pippo")
        cat("\n not able to load", runFile)
        
        return(NULL)
      })
      
      return(res.df)
      
    }))
    rownames(allRunsScores)<-seq(1, nrow(allRunsScores))
    
    
    allRunsParameters<-do.call(rbind, lapply(AllRunsFile, function(runFile){
      
      parsValues<-tryCatch({
        load(runFile)
        if (length(optimisedPars)==3 & all(optimisedPars %in% c("tau", "k", "n"))){
          parsValues<-opt_pars$parValues
          names(parsValues)<-opt_pars$parNames
        }else if (length(optimisedPars)==2 & all(optimisedPars %in% c("tau", "k"))){
          parsValues<-opt_pars$parValues[c(opt_pars$index_k, opt_pars$index_tau)]  
          names(parsValues)<-opt_pars$parNames[c(opt_pars$index_k, opt_pars$index_tau)]
        }else if (length(optimisedPars)==1 & all(optimisedPars %in% c("k"))){
          parsValues<-opt_pars$parValues[opt_pars$index_k]
          names(parsValues)<-opt_pars$parNames[opt_pars$index_k]
        }else{
          cat("Error: please specify which parameters has been optimised")
        }
        parsNames<-names(parsValues)
        parsValues[grepl("_k_", parsNames) & parsValues<=0.001]<-0
        parsValues[grepl("tau_", parsNames) & parsValues<=0.001]<-0
        return(parsValues)
      },
      error=function(cond){
        message(cond)
        # cat("pippo")
        cat("\n not able to load", runFile)
        
        return(NULL)
      })
      #       
      return(parsValues)
    }))
    
    # if N is provided consider only first N runs of the bootstrap
    if (is.numeric(N) & N < nrow(allRunsScores)){
      allRunsScores<-allRunsScores[1:N,]
      allRunsParameters<-allRunsParameters[1:N,]
    }
    
    return(list(scores=allRunsScores, parameters=allRunsParameters))
    
  })
  
  # retrieve cell line from fileName
  cellLinesNames<-sapply(AllFiles, function(x) {
    strsplit(x, split="_")[[1]][1]
  })
  cellLinesNames<-sapply(cellLinesNames, function(x) {
    strsplit(x, split="-")[[1]][2]
  })
  
  names(allByCellLine)<-cellLinesNames
  
  setwd(currentFolder)
  
  return(allByCellLine)
}
