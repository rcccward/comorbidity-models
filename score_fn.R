score_fn<-function(pred_out,dat,binary, elixmat, iterations, size,covar_flag){
  ##libraries
  library(randomForest)
  library(pROC)
  library(Matrix)
  library(glmnet)
  library(lmtest)
  library(miscTools)
  library(BayesTree)
  library(unbalanced)
  library(mice)
  library(missForest)
  library(nnet)
  
  
  meanvarimpRF<-pred_out[[1]]
  meanvarimpREG<-pred_out[[2]]
  meanvarimpBART<-pred_out[[3]]
  test_id<-pred_out[[6]]
  icd<-colnames(binary)
  combined<-data.frame(icd,meanvarimpREG,meanvarimpRF, meanvarimpBART )
  colnames(combined)<-c("icd"," REG"," RF"," BART")
  sign_reg<-  sign(combined[,2]) 
  
  c2<-quantile(combined[,2])
  c3<-quantile(combined[,3])
  c4<-quantile(combined[,4])
  
  k=3 #this is the 50th quantile for c2 thorugh c4 index
  id2=which((combined[,2]>c2[4]|combined[,2]< c2[2])  & combined[,3] > c3[k]  &  combined[,4]  > c4[k] )
  resultsElix<-resultsIndex<-matrix(0,ncol=15,nrow=0)
  ###################################
  #           FUNCTIONS             #
  ###################################
  CONF<-function(vec){
    upper<-ceiling(.975*length(vec))
    lower<-ceiling(.025*length(vec))
    mean<-mean(vec)
    UCL<-sort(vec)[upper]
    LCL<-sort(vec)[lower]
    return(list(UCL,mean,LCL)) 
  }#ENDFN
  ###############################  
  #       GLM AUC FN        #####
  ###############################
  glmfun<-function(dat,mort,pred.elix){  
    mod<-glm(mort~., data=as.data.frame(dat),family=binomial(link="logit"))
    result1<-pred_fn(mod, mort, pred.elix,"unadj")
    return(list(result1))
  }#END FN
  ###############################
  #  PRED FN                    #
  ###############################
  pred_fn<-function(mod, mort, pred.elix,flag){
    pred<-predict(mod,type=c("response"))
    auc<-roc(mort~pred)     
    auc_ci<-ci.auc(auc)
    confusion<-table(round(pred,0),mort)
    if(length(confusion)==4){
      sens<-confusion[2,2]/(confusion[1,2]+confusion[2,2])
      spec<-confusion[1,1]/(confusion[1,1]+confusion[2,1])
      ppv<-confusion[2,2]/(confusion[2,1]+confusion[2,2])
      npv<-confusion[1,1]/(confusion[1,1]+confusion[1,2])
    }
    if(length(confusion)==2 & rownames(confusion)[1]=="0"){
      sens<-0
      spec<-1
      ppv<-0
      npv<-1
      confusion<-rbind(confusion, c(0,0)); rownames(confusion)[2]<-"1"
    }
    if(length(confusion)==2 ){
      sens<-1
      spec<-0
      ppv<-1
      npv<-0
      confusion<-rbind(c(0,0),confusion); rownames(confusion)[1]<-"0"
    }
    mort1<-as.numeric(levels(mort))[mort]
    brier<-sum((mort1-pred)**2)/length(mort1)
    compare<-cbind(round(pred.elix,0),round(pred,0),as.numeric(levels(mort)[mort]))
    
    
    event.up=event.down=nonevent.up=nonevent.down=0
    
    for(i in 1:nrow(compare)){
      event.up<-ifelse(compare[i,3]==1 & compare[i,2]==1 & compare[i,1]==0, event.up+1, event.up)
      event.down<-ifelse(compare[i,3]==1 & compare[i,2]==0 & compare[i,1]==1, event.down+1, event.down)
      nonevent.up<-ifelse(compare[i,3]==0 & compare[i,2]==1 & compare[i,1]==0, nonevent.up+1, nonevent.up)
      nonevent.down<-ifelse(compare[i,3]==0 & compare[i,2]==0 & compare[i,1]==1, nonevent.down+1, nonevent.down)
    }
    nri.event<-(event.up-event.down)/sum(compare[,3])
    nri.nonevent<-(nonevent.down-nonevent.up)/(nrow(compare)-sum(compare[,3]))
    return(list("auc"=auc[9],"auc_ci"=auc_ci,"confusion"=confusion,"sens"=sens,"spec"=spec,
                "ppv"=ppv,"npv"=npv,"nri.event"=nri.event,"nri.nonevent"=nri.nonevent,"brier"=brier,
                "pred"=pred,"specificities"=auc$specificities,"sensitivities"=auc$sensitivities))
  }
  ###############################  
  #      SPARSE CONVERT         #
  #      PROB FEATURES          #
  ###############################
  sparse.convert<-function(A,prob){
    A= Matrix(A, sparse=T)
    #rm(list=ls()[10])
    #multiply features by prob
    A=t(A)
    A=prob*A
    A=t(A)
    return(A)
  }#END FUNCTION
  
  ###############################  
  #      ITERATIONS             #
  ###############################
  for(b in 1:iterations){
    cat("Elix  b=",b,"\n")
    
    test<-sample(test_id,size,replace=T)
    
    
    binary_test<-binary[test,]
    cov_test<-dat[test,c(-2,-6,-7,-8)]
    mort_test<-dat[test,6]
    
    elix_test<-elixmat[test,-c(11,17)] #omits HIV and DM for very high and low cell counts: complete separation possible        
    
    ###ELIX####
    if(covar_flag==0) {modElix<-glmfun(elix_test,mort_test,pred.elix=rep(0,nrow(binary_test)))
    } else {
      
      modElix<-glmfun(cbind(elix_test,cov_test),mort_test,pred.elix=rep(0,nrow(binary_test)))
    }
    pred.elix<-unlist(modElix[[1]][11])
    
    resultsElix<-rbind(resultsElix,unlist(modElix[[1]][1:10]))
    
    ####SUMMARY SCORE ######### 
    idx<-combined[id2,1]
    
    loc<-rep(NA,length(idx));  param<-rep(0,ncol(binary_test));loc1<-rep(NA,length(idx))
    for(i in 1:length(idx)){
      ixx<-which(colnames(binary_test)%in%as.character(idx[i]))#loc is index of binary_test columns
      if(length(ixx)>0)loc[i]<-ixx#loc is index of binary_test columns
      loc1[i]<-which(combined[,1]%in%as.character(idx[i])) #loc1 is index of combined rows and sign_reg
      param[loc[i]]<-sign_reg[loc1[i]]
      
    }
    
    binary_num<-apply(binary_test,2,function(x) as.numeric(as.character(x)))
    binnew<-sparse.convert(binary_num,param)
    
    test_index<-rowSums(binnew)
    if(covar_flag==0){
      modIndex<- glmfun(test_index,mort_test,pred.elix )
    } else{
      modIndex<- glmfun(cbind(test_index,cov_test),mort_test,pred.elix )
    }
    resultsIndex<-rbind(resultsIndex,unlist(modIndex[[1]][1:10]))
    
    
  }#next b 
  ###Generate output table from results:
  output<-matrix(0,nrow=18,ncol=2)
  colnames(output)<-c("Elix ","Summary Score")
  rownames(output)<-c("AUC-UCL","AUC","AUC-LCL","sens-UCL","sens","sens-LCL","spec-UCL","spec","spec-LCL",
                      "Brier-UCL","Brier","Brier-LCL","NRIevent-UCL","NRIevent","NRIevent-LCL",
                      "NRInonevent-UCL","NRInonevent","NRInonevent-LCL")
  collect<-function(datt){
    new<-c()
    for(i in c(1,9,10,15,13,14)){
      out<-CONF(datt[,i])
      new<-c(new,unlist(out))
    }
    return(new)
  }
  output[,1 ]<-collect(resultsElix)
  output[,2]<-collect(resultsIndex)
  
  output<-round(output,3)
  
  
  #####print out section#####
  cat("\n\n\n\n\n") 
  
  if(covar_flag==1)cat("Score performance (covariates included in model)\n")
  if(covar_flag==0)cat("Score performance (covariates not included in model)\n")
  cat(iterations," iterations\n") 
  cat(noquote(sprintf("%-17s% -1s% -5s%  s\n","","Elix","","Summary Score")))
  for(i in 1:nrow(output)){
    if(i%in%c(4,7,10,13,16))cat("\n")
    cat(noquote(sprintf("%-16s% .3f%-3s% .3f\n",
                        rownames(output)[i],output[i,1]," ",output[i,2])))
  }
  cat("\n\n\n\n\n")
  ########################
  comorbidities=idx
  weights=sign_reg[id2] 
  return(list(comorbidities,weights,output))
  #idx are the ICD codes in the index; param is the weight, output is the performance matrix
  
}#end function
