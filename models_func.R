#USER INSTRUCTIONS:
#LOAD DATASETS DAT, BINARY, ELIXMAT
#SET OTHER PARAMETERS:
covar_flag=1  # value=1 indicates covariates(age, gender, marital status, race) should be used in addition to ICD codes to build and validate models
              # value=0 indicates covarates should not be used
iterations=1000 #number of boostrapped samples used to find the distribution of performance statistics
size=2000 # obs er bootstrapped sample

#######################
#      FUNCTION       #
#######################

models_func<-function(dat,binary, elixmat, iterations, size,covar_flag){
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
  library(icd9)
  
   #ESTABLISH PARTITION OF DATASET INTO TRAINING AND VALIDATION SUBSETS
  #PARTITION MAINTAINED THROUGHOUT REMAINING STEPS FOR SECTIONS 1 AND 2 (MODELS AND SUMMARY SCORE)
  train_id<-sample(1:nrow(dat),.5*nrow(dat))
  test_id<-c(1:nrow(dat))[-train_id]
  
  #ESTABLISH EMPTY MATRIX / VECTOR TO RECORD RESULTS
  resultsElix<-  resultsRF<-resultsREG<- resultsComb<-resultsBART<-matrix(0,ncol=15,nrow=0)
  meanvarimpRF<-meanvarimpREG<- meanvarimpBART<-c(rep(0,ncol(binary)))
  
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
  #     FACTORIZE FUNCTION      #
  ###############################
  factorize<-function(fulldat){
    fulldat<-as.data.frame(fulldat)
    for(q in 1:ncol(fulldat)){
      fulldat[,q]<-as.factor(fulldat[,q])
    }
    return(fulldat)
  }#END FN
  
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
  #      ITERATIONS         #####
  ###############################
  
  for(b in 1:iterations){
    cat("Elix  iter=",b,"\n")
    train<-sample(train_id,  size, replace=T)
    test<-sample(test_id,size,replace=T)
    binary_train<-binary[train,]
    binary_test<-binary[test,]
    mort_train<-dat[train,6]
    mort_test<-dat[test,6]
    elix_test<-elixmat[test,-c(11,17)]#omits HIV and DM for very high and low cell counts: complete separation possible
    testval=colSums(apply(elix_test,2,as.numeric))
    if(min(testval)==0|max(testval)==nrow(elix_test))cat("perfect separation error\n")
    cov_train<-dat[train,c(-2,-6,-7,-8)]
    cov_test<-dat[test,c(-2, -6,-7,-8)]
    # ##adjust for unbalanced data:
    # newdat<-ubUnder(train,mort,method="percPos",perc=50)
    # binary_train<-newdat[[1]]
    # mort_train<-newdat[[2]]
    ###ELIX####
    if(covar_flag==0) {modElix<-glmfun(elix_test,mort_test,pred.elix=rep(0,nrow(binary_test)))
    } else {
      
      modElix<-glmfun(cbind(elix_test,cov_test),mort_test,pred.elix=rep(0,nrow(binary_test)))
    }
    pred.elix<-unlist(modElix[[1]][11])
    
    resultsElix<-rbind(resultsElix,unlist(modElix[[1]][1:10]))
    
    ##RF##### 
    cat("RF  iter=",b,"\n")
    if(covar_flag==0){
      forest<-randomForest(x=binary_train,y=mort_train, ntree = 200)
      rf.pred<-predict(forest,newdata=binary_test,type="prob")
    } else{
      forest<-randomForest(x=cbind(binary_train,cov_train),y=mort_train, ntree = 200)
      rf.pred<-predict(forest,newdata=cbind(binary_test,cov_test),type="prob")
    }
    modRF<-glmfun(rf.pred[,2],mort_test,pred.elix )
    resultsRF<-rbind(resultsRF, unlist(modRF[[1]][1:10]))
    
    meanvarimpRF<-(meanvarimpRF*(b-1)+forest$importance[1:1000])/b
    ### BART   ####
    cat("bart iter=",b,"\n")
    if(covar_flag==0){
      mod<-bart(x.train=binary_train,y.train=mort_train,x.test=binary_test,verbose=F)
    } else {
      mod<-bart(x.train=cbind(binary_train,cov_train),y.train=mort_train,x.test=cbind(binary_test,cov_test),verbose=F)
    }
    
    med.out<-c()
    for(i in  1:ncol(mod$yhat.test)){
      output<-pnorm(mod$yhat.test[,i])
      med.out<-c(med.out,median(output))
    }
    
    modBART<-glmfun(med.out,mort_test,pred.elix)
    resultsBART<-rbind(resultsBART, unlist(modBART[[1]][1:10]))
    
    meanvarimpBART<-(meanvarimpBART*(b-1) + colSums(mod$varcount[,1:1000]))/b
    
    ###switch to numeric binary matrices from factor:
    binary_test<-Matrix(data.matrix(binary_test),sparse=T)
    binary_test<-as(binary_test,"matrix")-1
    binary_train<-Matrix(data.matrix(binary_train),sparse=T)
    binary_train<-as(binary_train,"matrix")-1
    cov_train<-model.matrix(~race+age+single+male,data=cov_train)[,-1]
    cov_train<-Matrix(data.matrix(cov_train),sparse=T)
    cov_test<-model.matrix(~race+age+single+male,data=cov_test)[,-1]
    cov_test<-Matrix(data.matrix(cov_test),sparse=T)
    ####REGULARIZED REGRESSION####
     
    cat("reg iter=",b,"\n")
    if(covar_flag==0){
      fit<-cv.glmnet(as.matrix(binary_train), mort_train, family="binomial",alpha=.5,intercept=1)
      fit2<-glmnet(data.matrix(binary_train), mort_train,family="binomial",alpha=.5, lambda=fit$lambda.min,intercept=1)
      beta_reg<-fit2$beta
      binary_test<- Matrix(data.matrix(binary_test), sparse=T)
    } else{
      fit<-cv.glmnet(as.matrix(cbind(binary_train,cov_train)), mort_train, family="binomial",alpha=.5,intercept=1)
      fit2<-glmnet(data.matrix(cbind(binary_train,cov_train)), mort_train,family="binomial",alpha=.5, lambda=fit$lambda.min,intercept=1)
      beta_reg<-fit2$beta
      binary_test<- Matrix(data.matrix(cbind(binary_test,cov_test)), sparse=T)
    }
    bx<-(binary_test)%*%(beta_reg)
    predreg<-(1/(1+exp(-bx)))
    predreg<-as(predreg,"matrix")
    modREG<-glmfun(predreg,mort_test,pred.elix)
    resultsREG<-rbind(resultsREG, unlist(modREG[[1]][1:10]))
    meanvarimpREG<-(meanvarimpREG*(b-1) + as(beta_reg,"matrix")[1:1000,1])/b
    ##################################
    cat("comb iter=",b,"\n")
    
    combined<-cbind(unlist(modElix[[1]][11]),rf.pred[,2],predreg,med.out)
    modComb<-glmfun(combined,mort_test,pred.elix)
    resultsComb<-rbind(resultsComb, unlist(modComb[[1]][1:10]))
     
    
    
  } #next b
  
  ##generate output table:
  output<-matrix(0,nrow=18,ncol=5)
  colnames(output)<-c("Elix ","RF","BART","REG","Combined")
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
  output[,2]<-collect(resultsRF)
  output[,3]<-collect(resultsBART)
  output[,4]<-collect(resultsREG)
  output[,5]<-collect(resultsComb)
  output<-round(output,3)
  
  
  #####print out section#####
  cat("\n\n\n\n\n")
  
  if(covar_flag==1)cat("Model Performance Comparison (covariates included in model)\n")
  if(covar_flag==0)cat("Model Performance Comparison (covariates not included in model)\n")
  cat(iterations," iterations\n")
  cat("size = ", size,"patients\n")
  cat(noquote(sprintf("%-17s% -1s% -5s% -3s%-6s% -1s% -4s% -2s% -6s% -1s  \n","","Elix","","RF","","BART","","REG","","Pool")))
  for(i in 1:nrow(output)){
    if(i%in%c(4,7,10,13,16))cat("\n")
    cat(noquote(sprintf("%-16s% .3f%-3s% .3f%-3s% .3f%-3s%.3f%-3s% .3f\n",
                        rownames(output)[i],output[i,1]," ",output[i,2]," ",output[i,3]," ",output[i,4],"",output[i,5])))
  }
  cat("\n\n\n\n\n")
  ######### 
  return(list(meanvarimpRF,meanvarimpREG,meanvarimpBART,output, train_id,test_id))
  
}#end models function
