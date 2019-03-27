if(!require(foreach)){
  install.packages("foreach")
  library("foreach")
}

if(!require(plsdepot)){
  install.packages("plsdepot")
  library("plsdepot")
}

if(!require(plspm)){
  install.packages("plspm")
  library("plspm")
}

if(!require(pls)){
  install.packages("pls")
  library(pls)
}

xPLS <- function(YRange,XRange, in_Pre, mydata, BaseDir){
    for (se in c(1:length(in_Pre))) {
      yData = mydata[,YRange] # Specify Y variable
      xData = mydata[,XRange[1]:XRange[2]] # Specify X variables
      iPre=in_Pre[se] # Get output prefix
      iData=xData # Make a copy of xdata
      varNames=names(xData) # Get a list of variable names
      newLen=dim(iData)[2] # Set iteration limit (total number of variables)
      
      nVars=dim(iData)[2]
      
      oVarList=list() # Output variable list
      oRMSList=c() # Output RMS vector
      oRSQList=c() # Output RSQ vector
      oQSQList=c() # Output QSQ vector
      oAICList= c() # Output AIC vector
      i=0 # ...Output variable list counter
      while (newLen>2){ # Do till variables left > 2
        i=i+1 # ...increment list counter
        print(paste("Iteration",i,"Num. Vars. left:",dim(iData)[2]))
        flush.console()
        varList=c() # ...intermediate variable name vector
        rmsList=c() # ...intermediate RMSE vector
        rsqList=c() # ...intermediate RSQ vector
        qsqList=c() # ...intermediate QSQ vector
        AICList = c() #...intermediate AIC vector
        for (v in varNames) { # Do for each variable
          intData = iData # Copy over int Dataframe
          intData[v]=NULL # Remove current variable
          intPls = plsreg1(intData,yData,crosval=TRUE) # Run PLS regression
          intRsq = sum(intPls$R2) # Get r-squared
          intQsq = intPls$Q2[,5][length(intPls$Q2[,5])] # Get Q-squared
          intRMS = intPls$Q2[,1][length(intPls$Q2[,1])] # Get PRESS
          #intAIC = AIC(lm(yData~intData))
          print(paste("Testing",v,"PRESS:",round(intRMS,3),"RSQ:",round(intRsq,3),"QSQ:",round(intQsq,3)))
          flush.console()
          varList = c(varList,v) # Store dropped variable
          rmsList = c(rmsList,intRMS) # Store PRESS
          rsqList = c(rsqList,intRsq) # Store RSQ
          qsqList = c(rsqList,intRsq) # Store QSQ
          #AICList = c(AICList, intAIC) #Store AIC
        }
        remVar=varList[qsqList==max(qsqList)][1] # Find variable with maximum QSQ
        # remVar=varList[rsqList==max(rsqList)][1] # Find variable with maximum RSQ
        # remVar=varList[rmsList==min(rmsList)][1] # Find variable with minimum RMSE
        iData[remVar]=NULL # Drop from dataframe
        varNames=names(iData) # Refresh varNames
        intPls1 = plsreg1(iData,yData,comps = 8, crosval=TRUE) # Refit model with dropped variable
        intRsq1 = sum(intPls1$R2) # Get r-squared
        intRMS1 = intPls1$Q2[,1][length(intPls1$Q2[,1])] # Get PRESS
        intQsq1 = intPls1$Q2[,5][length(intPls1$Q2[,5])] # Get QSQ
        #intAIC1 = AIC(lm(yData~iData))
        oVarList[i] = list(names(iData)) # Append current variables to list
        oRMSList= c(oRMSList,intRMS1) # Append current RMS to output vector
        oRSQList= c(oRSQList,intRsq1) # Append current RSQ to output vector
        oQSQList= c(oQSQList,intQsq1) # Append current QSQ to output vector
        #oAICList= c(oAICList, intAIC1)
        print(paste("Iteration",i,"Variable removed",remVar,"current RMSE",round(intRMS1,3),"current RSQ",round(intRsq1,3),"current QSQ",round(intQsq1,3)))
        flush.console()
        newLen=dim(iData)[2] # Update number of variables left
      }
      
      print(paste(YRange,iPre,"PROCESS COMPLETE!!"))
      
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
      
      outFStub=paste(BaseDir,YRange,"_",iPre,"_",sep="")
      
      outPlotF=paste(outFStub,"_plot.png",sep="")
      outRMSEF=paste(outFStub,"_rmse.txt",sep="")
      outRSQRF=paste(outFStub,"_rsqr.txt",sep="")
      outQSQRF=paste(outFStub,"_qsqr.txt",sep="")
      outVARSF=paste(outFStub,"_vars.txt",sep="")
      outMODSF=paste(outFStub,"_mods.txt",sep="")
      
      png(outPlotF,pointsize=12,width=8,height=8,units="in",bg="white",res=120)
      par(mfcol=c(3,2),mar=c(5,4,1,1),oma=c(0,1,1,1))
      plot(c(1:length(oRMSList)),oRMSList,xlab="Iteration",ylab="PRESS",type="b")
      plot(c(1:length(oRSQList)),oRSQList,xlab="Iteration",ylab="RSQ",type="b")
      plot(oRMSList,oRSQList,,xlab="PRESS",ylab="RSQ",type="b")
      
      plot(c(1:length(oQSQList)),oQSQList,xlab="Iteration",ylab="QSQ",type="b")
      plot(oQSQList,oRSQList,xlab="Q-squared",ylab="RSQ",type="b")
      plot(oQSQList,oRMSList,xlab="Q-squared",ylab="PRESS",type="b")
      dev.off()
      
      dput(oVarList,file=outVARSF)
      dput(oRMSList,file=outRMSEF)
      dput(oRSQList,file=outRSQRF)
      dput(oQSQList,file=outQSQRF)
      dput(oFinList,file=outMODSF)
  }
}



