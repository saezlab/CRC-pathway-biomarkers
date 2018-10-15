# #########################################################################################################
# Example code for model optimisation using CNORode2017
# please see also desciption at: https://github.com/saezlab/CNORode2017
# #########################################################################################################

# NOTE: libraries should be loaded in this exact order:
# (this is because CNORode is a dependency of MEIGOR and has same function names as CNORode2017
# so it will be prioritizes when calling the optimisation function if loaded after CNORode2017)
library(CellNOptR)
library(MEIGOR)
library(CNORode2017)

# ****************
# working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load Prior Knowledge Network (PKN)
pknmodel<-readSIF("../data/base/PKN.sif")

# load normalised perturbation data 
MIDASfile <- "../data/processed/MIDAS/selectMIDASfile.csv" # select MIDAS file for the desired cell line
Mydata<-readMIDAS(MIDASfile=MIDASfile)
cnolist<-makeCNOlist(Mydata, subfield=F)
cnolist$valueStimuli[cnolist$valueStimuli==0]=0.5

# compress the network (no expansion, only OR gates are considered)
model<-preprocessing(data=cnolist, model=pknmodel, compression=TRUE, expansion=FALSE)

# set initial parameters (here parameters 'k' and 'tau' are optimised and 'n' fixed to 3)
ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0,
                                   LB_tau = 0, UB_n = 3, UB_k = 1, UB_tau = 1, default_n = 3,
                                   default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                                   opt_tau = TRUE, random = TRUE)

# PLX -> BRAF is an artificial regulation used to model paratoxical effect of PLX4720,
# which works as selective BRAF inhibitor in cell-lines where BRAF is mutated in V600E (i.e. HT29 and SNUC5 in our panel),
# but induces a paradoxical activation of wild type BRAF cells (modeled as stimulus on those cell lines)
ode_parameters$parValues[which(ode_parameters$parNames=="PLX_k_BRAF")]<-0.5
ode_parameters$index_opt_pars<-setdiff(ode_parameters$index_opt_pars, which(ode_parameters$parNames=="PLX_k_BRAF"))

## Parameter Optimization
# essm
paramsSSm=defaultParametersSSm()
paramsSSm$local_solver = "DHC"
paramsSSm$maxtime = 36000;
paramsSSm$maxeval = Inf;
paramsSSm$atol=1e-6;
paramsSSm$reltol=1e-6;
paramsSSm$nan_fac=1000;
paramsSSm$dim_refset=30;
paramsSSm$n_diverse=1000;
paramsSSm$maxStepSize=Inf;
paramsSSm$maxNumSteps=10000;
paramsSSm$transfer_function = 4;

paramsSSm$lambda_tau=0
paramsSSm$lambda_k=0
paramsSSm$bootstrap=F
paramsSSm$SSpenalty_fac=10
paramsSSm$SScontrolPenalty_fac=1000

opt_pars=parEstimationLBode(cnolist, model, method="essm", ode_parameters=ode_parameters, paramsSSm=paramsSSm, lambda=lambda)


