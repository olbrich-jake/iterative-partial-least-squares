###Function: iterative Partial Least Squares (x)
###Author: Aditya Singh and Jake Olbrich
###R. Version 3.4.3
###Purpose: iterate through a series of possible models running partial least squares to determine
###models with the lowest r squared and RMSE. Univariate.


###Packages required for xPLS

if(!require(foreach)){
  install.packages("foreach")
  library("foreach")
  # parallel for looping during component selection
}

if(!require(doParallel)){
  install.packages("doParallel")
  library(doParallel)
  #parallel backloops
}

if(!require(plsdepot)){
  install.packages("plsdepot")
  library("plsdepot")
}

if(!require(plspm)){
  install.packages("plspm")
  library("plspm")
  # plsreg1 to run univariate pls
}

if(!require(pls)){
  install.packages("pls")
  library(pls)
}

# if(!require(selectRcomp1)){
#   install.packages("selectRcomp1")
#   library(selectRcomp1)
#   #shoved selectNcomp in a random package
# }

#,"selectRcomp1" add back into line 74

##Data Requirements
##YRange: Range of Y variables
##XRange: Range of X variables
#mydata: dataframe being used
##BaseDir: remove later, for now file location to store outputs
xPLS <- function(YData, XRange, data = NULL){
  yData = mydata[,YData] #Set Y variables
  xData = mydata[,XRange[1]:XRange[2]] #Set X variables
  iData = xData #copy X data
  varNames = names(xData) #Get names of xData
  newLen = dim(iData)[2] #Set iteration limit
  nVars = dim(iData)[2] #Set number of variables
  if (print.progress) {
    cat("\nPLS progress\n")
    pb <- txtProgressBar(min = 0, max = nVars-2, style = 3)
    step <- 1
  }
  oVarList = list() # Output variable list
  oRMSList = c() # Output RMS vector
  oRSQList = c() # Output RSQ vector
  oQSQList = c() # Output QSQ vector
  i=0 #output variable list counter
  while (newLen>2){ # run until variables left > 2
    i=i+1
    print(paste("Iteration",i,"Num. Vars. left:",dim(iData)[2])) #Prints what we are running
    flush.console()
    varList=c() # ...intermediate variable name vector
    rmsList=c() # ...intermediate RMSE vector
    rsqList=c() # ...intermediate RSQ vector
    qsqList=c() # ...intermediate QSQ vector
    minNcomps = c()
    cores=detectCores()
    cl <- makeCluster(cores[1]-1) #not to overload your computer
    registerDoSNOW(cl) 
    newY = as.matrix(yData, ncol = 1, nrow = length(yData)) #Turns yData to a matrix
    foreach(j = 1:length(varNames), .packages=c("pls","plsdepot")) %dopar%  {# Do for each variable
      paste("hello world")
      v = varNames[j]
      intData = iData # Copy over int Dataframe
      intData[v] = NULL # Remove current variable
      newX = as.matrix(intData, ncol=length(varNames), nrow = length(yData)) #Turns intData to a matrix
      nPls = plsr(newY~newX, data=mydata, scale =TRUE,validation = "LOO") 
      newcomp = selectNcomp(nPls,  data=mydata, method="onesigma", plot = FALSE) #select the optimal number of components
      intPls = plsreg1(intData,yData,comps = newcomp, crosval=TRUE) # Run PLS regression
      intRsq = sum(intPls$R2) # Get r-squared
      intQsq = intPls$Q2[,5][length(intPls$Q2[,5])] # Get Q-squared
      intRMS = intPls$Q2[,1][length(intPls$Q2[,1])] # Get PRESS
      print(paste("Testing",v,"PRESS:",round(intRMS,3),"RSQ:",round(intRsq,3),"QSQ:",round(intQsq,3)))
      flush.console()
      varList = c(varList,v) # Store dropped variable
      rmsList = c(rmsList,intRMS) # Store PRESS
      rsqList = c(rsqList,intRsq) # Store RSQ
      qsqList = c(rsqList,intRsq) # Store QSQ
      minNcomps = c(minNcomps, newcomp)
    }
    stopCluster(cl)
    remVar=varList[qsqList==max(qsqList)][1] # Find variable with maximum QSQ
    # remVar=varList[rsqList==max(rsqList)][1] # Find variable with maximum RSQ
    # remVar=varList[rmsList==min(rmsList)][1] # Find variable with minimum RMSE
    iData[remVar]=NULL # Drop from dataframe
    varNames=names(iData) # Refresh varNames
    intcomps = min(minNcomps)
    intPls1 = plsreg1(iData,yData,comps = , crosval=TRUE) # Refit model with dropped variable
    intRsq1 = sum(intPls1$R2) # Get r-squared
    intRMS1 = intPls1$Q2[,1][length(intPls1$Q2[,1])] # Get PRESS
    intQsq1 = intPls1$Q2[,5][length(intPls1$Q2[,5])] # Get QSQ
    oVarList[i] = list(names(iData)) # Append current variables to list
    oRMSList = c(oRMSList,intRMS1) # Append current RMS to output vector
    oRSQList = c(oRSQList,intRsq1) # Append current RSQ to output vector
    oQSQList = c(oQSQList,intQsq1) # Append current QSQ to output vector
    print(paste("Iteration",i,"Variable removed",remVar,"current RMSE",round(intRMS1,3),"current RSQ",round(intRsq1,3),"current QSQ",round(intQsq1,3)))
    flush.console()
    newLen = dim(iData)[2] # Update number of variables left
  print(paste(YData,"PROCESS COMPLETE!!"))
  oFinList=list()
  oFinList[1]=paste("Minimum RMS",min(oRMSList))
  oFinList[2]=paste("Rsq of model with Min RMS:",oRSQList[oRMSList==min(oRMSList)])
  outVAR=oVarList[oRMSList==min(oRMSList)][[1]]
  oFinList[3]=paste("Variables in model with minimum RMS, (",length(outVAR),")")
  oFinList[4]=list(outVAR)
  oFinList[5]=paste("Maximum RSQ",max(oRSQList))
  oFinList[6]=paste("RMS of model with Max RSQ:",oRMSList[oRSQList==max(oRSQList)])
  outVRS=oVarList[oRSQList==max(oRSQList)][[1]]
  oFinList[7]=paste("Variables in model with maximum RSQ, (",length(outVRS),")")
  oFinList[8]=list(outVRS)
  }
  z = list(call=match.call(), mod=oFinList, qsqr=oQSQList, rmse=oRMSList, rsqr=oRSQList, vars=oVarList)
  return(z)
}