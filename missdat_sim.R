#USER INSTRUCTIONS: Set missType, pctMiss, MItype, Nobs, iterations
#Possible values are missType, pctMiss, MIType are shown.
#missType<-c("MAR","MCAR","MNAR")
#pctMiss<-c(.1,.3,.5)
#MItype<-c("MICE-LR","MICE-RF","MICE-NNET","MissF","CC")    #MICE and Logistic Reg., MICE and RandomForest, MissForest, Complete Case
#Example:
MissType="MNAR"
pctMiss=.3
Nobs=1000
iterations=1000
missrun<-missdat_sim(dat,MissType,pctMiss,Nobs,iterations)




missdat_sim<-function(dat, MissType, pctMiss, Nobs,iterations){
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
  ########################################
  #   MISSING DATA ID GENERATOR          #
  ########################################
  miss_gen<-function(m,N,ncol,MissType,simdat){
    #MCAR GENERATOR:
    iter=0
    while(iter<1000){
      if(MissType=="MCAR"){
        miss_id<-matrix(1,N,ncol)
        for(i in 1:3){
          #id<-(unique(round(runif(n=2*m*N,min=1,max=nrow(simdat)))))[1:(m*N)]
          id<-sample(1:N,round(m*N/3,0), replace=F)
          miss_id[id,i]<-0
        }
      } #end MCAR
      
      #MAR GENERATOR:
      if(MissType=="MAR"){
        #this uses categorical variables to determine missingness
        #higher rank of determiner has greater  prob missingness
        #single determines meanmpr missingness
        #rural determines meana1c missingness
        #death determine race missingness
        #m*N missing values per column
        
        miss_id<-matrix(rbinom(N*ncol,1,1-m*.05),N,ncol)
        colid<-c(6,1,2)#death,single,rural>race,meanmpr,meana1c
        for(i in 1:3){
          missNeed<- m*N/3-(N-sum(miss_id[,i]))
          one_id<-which(simdat[,colid[i]]==1);zero_id<-which(simdat[,colid[i]]==0)
          max_fraction<-table(simdat[,colid[i]])[2]/N
          use_fraction<-min(missNeed/N,max_fraction)
          #missing race: use  death
          #1/3 missing from death=0,  2/3 missing  from  death=1
          id1<-zero_id[sample(1:length(zero_id),ceiling(N*(use_fraction/8)),replace=F)]
          id2<-one_id[sample(1:length(one_id),ceiling(N*(7*use_fraction/8)),replace=F)]
          miss_id[c(id1,id2),i]<-0	
          missNeed<- m*N/3-(N-sum(miss_id[,i]))
          if(missNeed>0){     zero_id<-which(miss_id[,i]==1)
          id3<-zero_id[sample(1:length(zero_id),missNeed,replace=F)]
          miss_id[id3,i]<-0
          }
          if(missNeed<0){     one_id<-which(miss_id[,i]==0)
          id3<-one_id[sample(1:length(one_id),abs(missNeed),replace=F)]
          miss_id[id3,i]<-1
          }
          
          
        }#nexti
        
      } #end MAR
      
      #MNAR generator: 
      #missingness depends on the variables themselves: race, a1c, mpr:
      if(MissType=="MNAR"){
        miss_id<-matrix(rbinom(N*ncol,1,1-m*.1),N,ncol)
        
        meda1c<-median(simdat$meana1c)
        medmpr<-median(simdat$meanmpr)
        for (i in  1:N){
          if((simdat$race[i]==2|simdat$race[i]==3) & simdat$death[i]==1)miss_id[i,1]<- rbinom(1,1,1-.9)
          if(simdat$meanmpr[i]>medmpr) miss_id[i,2]<-rbinom(1,1,1-.9)
          if(simdat$meana1c[i]>meda1c) miss_id[i,3]<-rbinom(1,1,1-.9)
          
        }
        for(j in 1:3){
          missNeed<- m*N/3-(N-sum(miss_id[,j]))
          if(missNeed>0){     zero_id<-which(miss_id[,j]==1)
          id3<-zero_id[sample(1:length(zero_id),missNeed,replace=F)]
          miss_id[id3,j]<-0
          }
          if(missNeed<0){     one_id<-which(miss_id[,j]==0)
          id3<-one_id[sample(1:length(one_id),abs(missNeed),replace=F)]
          miss_id[id3,j]<-1
          }
        }
        
      } #end MNAR
      colnames(miss_id)<-c("miss_race","miss_mpr","miss_a1c")
      simdat_miss<-cbind(simdat, miss_id)
      for(i in 1:nrow(simdat_miss)){
        if(simdat_miss$miss_race[i]==0)simdat_miss$race[i]<-NA
        if(simdat_miss$miss_a1c[i]==0)simdat_miss$meana1c[i]<-NA
        if(simdat_miss$miss_mpr[i]==0)simdat_miss$meanmpr[i]<-NA
      }
      #check complete cases:
      complete_id<-which(complete.cases(simdat_miss)==T)
      delta<-round(m*N-nrow(simdat_miss[!complete.cases(simdat_miss),]),0)
      if(delta>0){	  
        id3<-sample(complete_id,delta,replace=F)#this is added rows to make NA
        id4<-sample(c(5,7,8),delta,replace=T)#this is additional columns to make NA	
        for(i in 1:delta){		              
          simdat_miss[id3[i],id4[i]]<-NA
        }
      }
      #reset _miss ids to NA where needed:
      for(i in 1:nrow(simdat_miss)){
        if(is.na(simdat_miss$race[i])==T)simdat_miss$miss_race[i]<-0
        if(is.na(simdat_miss$meana1c[i])==T)simdat_miss$miss_a1c[i]<-0
        if(is.na(simdat_miss$meanmpr[i])==T)simdat_miss$miss_mpr[i]<-0
      }
      
      #test Odds Ratios
      flag=0
      #test missing  race
      mod_race<-glm(as.factor(abs(simdat_miss$miss_race))~simdat$male+
                      simdat$age+simdat$death+simdat$single+simdat$rural+simdat$meana1c+simdat$meanmpr+simdat$race, 
                    na.action=na.omit,family=binomial(link=logit))
      or_race<-exp(summary(mod_race)$coefficients[c(-1,-3),1])#exclude male because small number of females causes large OR
      
      #test missing  meanmpr
      mod_mpr<-glm(as.factor(simdat_miss$miss_mpr)~simdat$male+
                     simdat$age+simdat$death+simdat$single+simdat$rural+simdat$meana1c+simdat$race+simdat$meanmpr, 
                   na.action=na.omit,family=binomial(link=logit))
      
      or_mpr<-exp(summary(mod_mpr)$coefficients[c(-1,-3),1])
      
      #test missing  meana1c
      mod_a1c<-glm(as.factor(simdat_miss$miss_a1c)~simdat$male+
                     simdat$age+simdat$death+simdat$single+simdat$rural+simdat$race+simdat$meanmpr+simdat$meana1c, 
                   na.action=na.omit,family=binomial(link=logit))
      or_a1c<-exp(summary(mod_a1c)$coefficients[c(-1,-3),1])
      
      or_combined<-c(or_race,or_mpr,or_a1c)
      
      if(MissType=="MNAR"){if(or_race["simdat$death1"]>.667|or_race["simdat$race2"]>.667|or_race["simdat$race3"]>.667|or_mpr["simdat$meanmpr"]>.667|or_a1c["simdat$meana1c"]>.667)flag=1 }
      if(MissType=="MCAR"){ }#no conditions set
      if(MissType=="MAR"){if(or_race["simdat$death1"]>.667|or_mpr["simdat$single1"]>.667|or_a1c["simdat$rural1"]>.667)flag=1 }
      
      
      
      #check for zero cells in complete cases and simdat_miss:
      comp_case<-simdat_miss[complete_id,]
      if(0%in%table(comp_case$single,comp_case$death)| 0%in%table(comp_case$rural,comp_case$death)|
         0%in%table(comp_case$male,comp_case$death)| 0%in%table(comp_case$race,comp_case$death))flag=1
      if(0%in%table(simdat_miss$single,simdat_miss$death)| 0%in%table(simdat_miss$rural,simdat_miss$death)|
         0%in%table(simdat_miss$male,simdat_miss$death)|	0%in%table(simdat_miss$race,simdat_miss$death))flag=1
      
      if(flag==0)	{ return(simdat_miss)}
      
      iter=iter+1
    }#end iter
    return(flag)
    
  }# END FN
  #determine prob(missing race, mpr, a1c) for died/survived, rural/urban, single/married
  probs<-function(sim,dat ){#sim is the dataset with missingness imposed; dat is the original dataset
    
    aa=table(sim$miss_race,sim$death)
    p.mar.race.died= aa[3]/(aa[3]+aa[4])
    p.mar.race.surv=aa[1]/(aa[1]+aa[2])
    bb= table(sim$miss_mpr,sim$single)
    p.mar.mpr.single=bb[3]/(bb[3]+bb[4])
    p.mar.mpr.married=bb[1]/(bb[1]+bb[2])
    cc= table(sim$miss_a1c,sim$rural)
    p.mar.a1c.rural=cc[3]/(cc[3]+cc[4])
    p.mar.a1c.urban=cc[1]/(cc[1]+cc[2])
    
    a1c.ind<-which(dat$meana1c>=median(dat$meana1c))
    mpr.ind<-which(dat$meanmpr>=median(dat$meanmpr))
    a1c.binary<-mpr.binary<-rep(0,nrow(dat))
    a1c.binary[a1c.ind]<-1
    mpr.binary[mpr.ind]<-1
    dd=table(sim$miss_mpr,mpr.binary)
    p.mnar.mpr.below.median<-dd[1]/(dd[1]+dd[2])
    p.mnar.mpr.above.median<-dd[3]/(dd[3]+dd[4])
    ee=table(sim$miss_a1c,a1c.binary)
    p.mnar.a1c.below.median<-ee[1]/(ee[1]+ee[2])
    p.mnar.alc.above.median<-ee[3]/(ee[3]+ee[4])
    ff= table(sim$miss_race,dat$race,dat$death)
    p.mnar.nhb.died<-ff[1,2,2]/(ff[1,2,2]+ff[2,2,2])
    p.mnar.nhw.died<-ff[1,1,2]/(ff[2,1,2]+ff[1,1,2])
    p.mnar.hisp.died<-ff[1,3,2]/(ff[1,3,2]+ff[2,3,2])
    p.mnar.other.died<-ff[1,4,2]/(ff[1,4,2]+ff[2,4,2])
    p.mnar.nhb.surv<-ff[1,2,1]/(ff[1,2,1]+ff[2,2,1])
    p.mnar.nhw.surv<-ff[1,1,1]/(ff[2,1,1]+ff[1,1,1])
    p.mnar.hisp.surv<-ff[1,3,1]/(ff[1,3,1]+ff[2,3,1])
    p.mnar.other.surv<-ff[1,4,1]/(ff[1,4,1]+ff[2,4,1])
    labels<-c("p.mar.race.died","p.mar.race.surv","p.mar.mpr.single","p.mar.mpr.married","p.mar.a1c.rural","p.mar.a1c.urban",
              "p.mnar.mpr.below.median","p.mnar.mpr.above.median","p.mnar.a1c.below.median","p.mnar.alc.above.median","p.mnar.nhb.died","p.mnar.nhw.died",
              "p.mnar.hisp.died","p.mnar.other.died","p.mnar.nhb.surv","p.mnar.nhw.surv","p.mnar.hisp.surv","p.mnar.other.surv")
    values<-c(p.mar.race.died,p.mar.race.surv,p.mar.mpr.single,p.mar.mpr.married,p.mar.a1c.rural,p.mar.a1c.urban,
              p.mnar.mpr.below.median,p.mnar.mpr.above.median,p.mnar.a1c.below.median,p.mnar.alc.above.median,p.mnar.nhb.died,
              p.mnar.nhw.died,p.mnar.hisp.died,p.mnar.other.died,p.mnar.nhb.surv,p.mnar.nhw.surv,p.mnar.hisp.surv,p.mnar.other.surv)
    output<- data.frame(labels,values)
    return(output )
    
  }#end function
  ###############################  
  #     BART IMPUTE FUNCT       #
  ###############################
  
  mice.impute.bart<-function (y, ry, x, ...) 
  {
    if (!requireNamespace("BayesTree", quietly = TRUE)) 
      stop("Package 'BayesTree' needed for this function \n             to work. Please install it.", 
           call. = FALSE)
    if(is.factor(y)){y<-sapply(y, function(x) as.numeric(as.character(x)))}
    
    nmis <- sum(!ry)
    xobs <-  (x[ry, , drop = FALSE])
    xmis <- x[!ry, , drop = FALSE]
    yobs <- y[ry]
    mod<-bart(x.train=xobs,y.train=yobs,x.test=xmis,verbose=F)
    # if (nmis == 1) 
    #   forest <- array(forest, dim = c(1, ntree))
    impute<-apply(mod$yhat.test,2,function(x) median(pnorm(x)))
    
    return(impute)
  }
  ###############################  
  #    NNET IMPUTE FUNCT     ####
  ###############################
  mice.impute.Nnet<-function (y, ry, x, ...) {
    
    if (!requireNamespace("nnet", quietly = TRUE)) 
      stop("Package 'nnet' needed for this function \n             to work. Please install it.", 
           call. = FALSE)
    
    
    nmis <- sum(!ry)
    xobs <-  (x[ry, , drop = FALSE])
    xmis <- x[!ry, , drop = FALSE]
    yobs <- y[ry]
    if(is.numeric(y)){
      nn_mod<-nnet(xobs,yobs,   size=1, linout=T,rang=.1,decay=0, maxit=200)
      impute<-  (predict(nn_mod, xmis,type="raw"))
    }
    if(is.factor(y)){ 
      form<-as.formula(paste("yobs~", paste(names(xobs),collapse='+')))
      num<-length(yobs)
      wts<-num/table(y)
      data.weights<-rep(1,length(yobs))
      for(i in 1:length(yobs)){
        data.weights[i]<-wts[yobs[i]]
        
      }
      nn_mod<-nnet(form,data=xobs,size=1,   linout=T,rang=.1,decay=0, maxit=200)
      impute<- max.col(predict(nn_mod, xmis,type="raw"))
    }
    return(impute)
  }
  ###############################  
  #    MPBART IMPUTE FUNCT   ####
  ###############################
  #mpbart only considers multilevel y (multinomial)
  mice.impute.mpbart<-function (y, ry, x, ...) 
    #pp<-function (y, ry, x, ...) 
  {
    nmis <- sum(!ry)
    xobs <-  (x[ry, , drop = FALSE])
    xmis <- x[!ry, , drop = FALSE]
    yobs <- y[ry]
    
    #unbalanced data: randomly samples with replacement until there are equal numbers of NHW and other groups
    
    i.bin<-ifelse(yobs==1,0,1)
    race<-yobs
    xobs<-cbind(xobs,race)
    bal<-ubOver(xobs,i.bin,k=0,verbose=F)
    train<-bal[[1]] 
    
    p=4#number of classes in the outcome race
    #mpbart notes:
    cat(colnames(train),"/n")
    i.race<-which(colnames(train)=='race')
    colnames(xmis)<-colnames(train)[-i.race]
    prior1<-list(nu=p+2,V= 1*diag(p-1),ntrees=50,kfac=5, pbd=1, pb=.9, beta=1, alpha=.99,nc=400,priorindep=F,minobsnode=2)
    mcmc1<-list(sigma0=diag(p-1), keep=F, burn=100, ndraws=1000)
    out<-mpbart(race~1|.,
                train.data=train, test.data=xmis, base='1',Prior=prior1, Mcmc=mcmc1)
    
    #table(simTrain$race, as.character(out$predicted_class_train))
    
    impute<-out$predicted_class_test
    return(impute)
  }#end MPBART
  ###############################  
  #      COV PROB FUNCT         #
  ###############################
  
  #fit is the model fit from imputed data; mod_true is the model with true betas
  cov_prob<-function(fit_betas,fit_se, beta_true){
    UCL<-fit_betas+qt(.975,1000-9)*fit_se
    LCL<-fit_betas-qt(.975,1000-9)*fit_se
    cov_flag<-c(rep(0, length(fit_betas)))
    for(i in 1:length(fit_betas)){
      if(beta_true[i]<UCL[i] & beta_true[i]>LCL[i])cov_flag[i]<-1
    }
    return(cov_flag)
  }
  
  ####################### 
  #     ODDS RATIOS     #
  #######################
  
  odds_ratios<-function(simdat_miss,simdat){ #simdat_miss is the data with missingness imposed; simdat is the  original data
    #test missing  race
    mod_race<-glm(as.factor(abs(simdat_miss$miss_race))~simdat$male+
                    simdat$age+simdat$death+simdat$single+simdat$rural+simdat$meana1c+simdat$meanmpr+simdat$race, 
                  na.action=na.omit,family=binomial(link=logit))
    or_race<-exp(-summary(mod_race)$coefficients[-1,1]) 
    
    #test missing  meanmpr
    mod_mpr<-glm(as.factor(simdat_miss$miss_mpr)~simdat$male+
                   simdat$age+simdat$death+simdat$single+simdat$rural+simdat$meana1c+simdat$meanmpr+simdat$race,  
                 na.action=na.omit,family=binomial(link=logit))
    
    or_mpr<-exp(-summary(mod_mpr)$coefficients[-1,1])
    
    #test missing  meana1c
    mod_a1c<-glm(as.factor(simdat_miss$miss_a1c)~simdat$male+
                   simdat$age+simdat$death+simdat$single+simdat$rural+simdat$meana1c+simdat$meanmpr+simdat$race,  
                 na.action=na.omit,family=binomial(link=logit))
    or_a1c<-exp(-summary(mod_a1c)$coefficients[-1,1])
    output<-cbind(or_race,or_mpr,or_a1c)
    return(output)
  }#end function
  
  
  ##DETERMINE TRUE COEFFICIENTS FROM FULL DATA, ALL COMPLETE CASES:
  fulldat<- dat
  
  mod_true<-glm(death ~ single+rural+male+age+race+meana1c+meanmpr, data=fulldat , family="binomial"(link="logit"))
  beta_true<-summary(mod_true)$coefficients[,1]
  
  
  
  
  results<-list()
  odds.miss<-list()
  prob.miss<-list()
  
  #load("results_new.Rdata")
  #load("stats_boot.RData")
  #load("stats_50.RData")
  #load("stats_impute.RData")
  
  ############################
  #       ITERATIONS      ####
  ############################
  i<-which(c("MAR","MCAR","MNAR")%in%MissType)
  j<-which(c(.1,.3,.5) %in% pctMiss)
  if(length(j)==0|length(i)==0)stop("must enter correct options for MissType and/or\n pctMiss value of .1, .3, or .5 only\n")
  for(b in 1:iterations){
    
    
    rowid<-sample(1:nrow(dat), Nobs, replace=T)
    simdat<-dat[rowid,]
    
    #stats_boot[[b]]<-tables(simdat)
    results[[b]]<-matrix(0,nrow=1, ncol=13)
    colnames(results[[b]])<-c("missType","pctMiss","MIType","intercept","single1" ,"rural1","male1","age", "race2","race3",
                              "race4","meana1c", "meanmpr")
    prob.miss[[b]]<-list()
    odds.miss[[b]]<-list()
    count=1
    
    
    
    simdat_miss<-miss_gen(pctMiss, nrow(simdat), 3, MissType, simdat) 
    
    probdat<-probs(simdat_miss,simdat)
    odds<-odds_ratios(simdat_miss,simdat)
    prob.miss[[b]][[count]]<-list(i,j,probdat)
    odds.miss[[b]][[count]]<-list(i,j,odds)
    simdat_miss<-simdat_miss[,-c(9:11)]
    count=count+1
    
    #MICE WITH LOGISTIC IMPUTATION:
    k=1
    cat( "iteration",b,": ","missType = ",MissType, "pctMiss = ", pctMiss, "MItype= MICE-LR","\n")
    mice1<-mice(simdat_miss,
                m = 5,
                #method = vector("character", length = ncol(data)),
                #predictorMatrix = (1 - diag(1, ncol(data))),
                #visitSequence = (1:ncol(data))[apply(is.na(data), 2, any)],
                visitSequence='monotone',
                #form = vector("character", length = ncol(data)),
                #post = vector("character", length = ncol(data)),
                defaultMethod = c("pmm","logreg", "polyreg", "polr"),
                maxit = 10,
                diagnostics = TRUE,
                printFlag = TRUE,
                #seed = NA,
                #imputationMethod = NULL,
                #defaultImputationMethod = NULL,
                #data.init = NULL,
    )
    
    fit <- with(mice1, glm(death ~ single+rural+male+age+race+meana1c+meanmpr, family="binomial"(link="logit")))
    poolfit<-pool(fit)
    rel.bias<-(1-poolfit$qbar/beta_true)
    ratio_var<-diag(poolfit$t)/(summary(mod_true)$coefficients[,2])**2
    rmse<-sqrt((poolfit$qbar-beta_true)**2+diag(poolfit$t))
    fit_betas<-poolfit$qbar
    fit_se<-sqrt(diag(poolfit$t))
    covprob<-cov_prob(fit_betas,fit_se,beta_true)
    
    #update results:
    rel.bias<-c(i,j,k,rel.bias)
    ratio_var<-c(i,j,k,ratio_var)
    rmse<-c(i,j,k,rmse)
    covprob<-c(i,j,k,covprob)
    results[[b]]<-rbind(results[[b]],rel.bias,ratio_var,rmse,covprob)
    
    
    # #MICE WITH RF:
    k=2
    cat( "iteration",b,": ","missType = ",MissType, "pctMiss = ", pctMiss, "MItype= MICE-RF","\n")
    mice1<-mice(simdat_miss,
                m = 5,
                #method = vector("rf", length = ncol(simdat_miss)),
                #predictorMatrix = (1 - diag(1, ncol(data))),
                #visitSequence = (1:ncol(data))[apply(is.na(data), 2, any)],
                visitSequence='monotone',
                #form = vector("character", length = ncol(data)),
                #post = vector("character", length = ncol(data)),
                defaultMethod = c("rf","rf", "rf", "rf"),
                maxit = 5,
                diagnostics = TRUE,
                printFlag = TRUE,
                #seed = NA,
                #imputationMethod = NULL,
                #defaultImputationMethod = NULL,
                #data.init = NULL,
    )
    
    fit <- with(mice1, glm(death ~ single+rural+male+age+race+meana1c+meanmpr, family="binomial"(link="logit")))
    poolfit<-pool(fit)
    rel.bias<-(1-poolfit$qbar/beta_true)
    ratio_var<-diag(poolfit$t)/(summary(mod_true)$coefficients[,2])**2
    rmse<-sqrt((poolfit$qbar-beta_true)**2+diag(poolfit$t))
    fit_betas<-poolfit$qbar
    fit_se<-sqrt(diag(poolfit$t))
    covprob<-cov_prob(fit_betas,fit_se,beta_true)
    
    #update results:
    rel.bias<-c(i,j,k,rel.bias)
    ratio_var<-c(i,j,k,ratio_var)
    rmse<-c(i,j,k,rmse)
    covprob<-c(i,j,k,covprob)
    results[[b]]<-rbind(results[[b]],rel.bias,ratio_var,rmse,covprob)
    
    
    
    # MICE-BART+MPBART####
    #one iteration run: saved under "BART_NNET_Single_iteration_results.RData, opens as "results" object
    # k=1
    #cat( "iteration",b,": ",missType[i],pctMiss[j],MItype[k],"\n")
    # mice1<-mice(simdat_miss,
    #             m = 5,
    #             #method = vector("rf", length = ncol(simdat_miss)),
    #             #predictorMatrix = (1 - diag(1, ncol(data))),
    #             #visitSequence = (1:ncol(data))[apply(is.na(data), 2, any)],
    #             visitSequence='monotone',
    #             #form = vector("character", length = ncol(data)),
    #             #post = vector("character", length = ncol(data)),
    #             method = c(rep("bart",5),"mpbart",rep("bart",3)),
    #             maxit = 5,
    #             diagnostics = TRUE,
    #             printFlag = TRUE,
    #             #seed = NA,
    #             #imputationMethod = NULL,
    #             #defaultImputationMethod = NULL,
    #             #data.init = NULL,
    # )
    # 
    # fit <- with(mice1, glm(death ~ single+rural+male+age+race+meana1c+meanmpr, family="binomial"(link="logit")))
    # poolfit<-pool(fit)
    # tstat<-(poolfit$qbar-summary(mod_true)$coefficients[,1])/sqrt(diag(poolfit$t))
    # ratio_var<-diag(poolfit$t)/(summary(mod_true)$coefficients[,2])**2
    # rmse<-sqrt((poolfit$qbar-summary(mod_true)$coefficients[,1])**2+diag(poolfit$t))
    # fit_betas<-poolfit$qbar
    # fit_se<-sqrt(diag(poolfit$t))
    # 
    # covprob<-cov_prob(fit_betas,fit_se,beta_true)
    # 
    # #update results:
    # tstat<-c(i,j,k,tstat)
    # ratio_var<-c(i,j,k,ratio_var)
    # rmse<-c(i,j,k,rmse)
    # covprob<-c(i,j,k,covprob)
    # results[[b]]<-rbind(results[[b]],tstat,ratio_var,rmse,covprob)
    
    
    # MICE-NNET####
    
    k=3
    cat( "iteration",b,": ","missType = ",MissType, "pctMiss = ", pctMiss, "MItype= MICE-NNET","\n")
    mice1<-mice(simdat_miss,
                m = 5,
                #method = vector("rf", length = ncol(simdat_miss)),
                #predictorMatrix = (1 - diag(1, ncol(data))),
                #visitSequence = (1:ncol(data))[apply(is.na(data), 2, any)],
                visitSequence='monotone',
                #form = vector("character", length = ncol(data)),
                #post = vector("character", length = ncol(data)),
                method = "Nnet",
                maxit = 5,
                diagnostics = TRUE,
                printFlag = TRUE,
                #seed = NA,
                #imputationMethod = NULL,
                #defaultImputationMethod = NULL,
                #data.init = NULL,
    )
    
    fit <- with(mice1, glm(death ~ single+rural+male+age +race+meana1c+meanmpr, family="binomial"(link="logit")))
    poolfit<-pool(fit)
    rel.bias<-(1-poolfit$qbar/beta_true)
    ratio_var<-diag(poolfit$t)/(summary(mod_true)$coefficients[,2])**2
    rmse<-sqrt((poolfit$qbar-beta_true)**2+diag(poolfit$t))
    fit_betas<-poolfit$qbar
    fit_se<-sqrt(diag(poolfit$t))
    covprob<-cov_prob(fit_betas,fit_se,beta_true)
    
    #update results:
    rel.bias<-c(i,j,k,rel.bias)
    ratio_var<-c(i,j,k,ratio_var)
    rmse<-c(i,j,k,rmse)
    covprob<-c(i,j,k,covprob)
    results[[b]]<-rbind(results[[b]],rel.bias,ratio_var,rmse,covprob)
    
    #MISS FOREST:
    k=4
    cat( "iteration",b,": ","missType = ",MissType, "pctMiss = ", pctMiss, "MItype= MissForest","\n")
    # 
    MF<-missForest(simdat_miss , maxiter = 10, ntree = 100, variablewise = FALSE,
                   decreasing = FALSE, verbose = T,xtrue = simdat)
    mod_mf<-glm(as.factor(death) ~ single+rural+male+age+ race+meana1c+meanmpr, data=MF$ximp , family="binomial"(link="logit"))
    
    #MF$ximp is the  imputed set.
    # 
    # tstat<-(summary(mod_mf)$coefficients[,1]-summary(mod_true)$coefficients[,1])/summary(mod_mf)$coefficients[,2]
    rel.bias<-(1-summary(mod_mf)$coefficients[,1]/beta_true)
    ratio_var<-(summary(mod_mf)$coefficients[,2]**2)/(summary(mod_true)$coefficients[,2]**2)
    rmse<-sqrt((summary(mod_mf)$coefficients[,1]-summary(mod_true)$coefficients[,1])**2+(summary(mod_mf)$coefficients[,2])**2)
    fit_betas<-summary(mod_mf)$coefficients[,1]
    fit_se<-summary(mod_mf)$coefficients[,2]
    covprob<-cov_prob(fit_betas,fit_se,beta_true)
    
    #update results:
    rel.bias<-c(i,j,k,rel.bias)
    ratio_var<-c(i,j,k,ratio_var)
    rmse<-c(i,j,k,rmse)
    covprob<-c(i,j,k,covprob)
    results[[b]]<-rbind(results[[b]],rel.bias,ratio_var,rmse,covprob)
    # 
    # stats_impute[[count2]]<-list(b,i,j,k,tables(MF$ximp)); count2<-count2+1;
    # 
    # #COMPLETE CASE:
    k=5
    cat( "iteration",b,": ","missType = ",MissType, "pctMiss = ", pctMiss, "MItype= Complete Case","\n")
    cc_id<-which(complete.cases(simdat_miss)==T)
    
    mod_cc<-glm(as.factor(death) ~ single+rural+male+age +race+meana1c+meanmpr, data=simdat_miss[cc_id,] , family="binomial"(link="logit"))
    #tstat<-(summary(mod_cc)$coefficients[,1]-summary(mod_true)$coefficients[,1])/summary(mod_cc)$coefficients[,2]
    rel.bias<-(1-summary(mod_cc)$coefficients[,1]/beta_true)
    ratio_var<-(summary(mod_cc)$coefficients[,2]**2)/(summary(mod_true)$coefficients[,2]**2)
    rmse<-sqrt((summary(mod_cc)$coefficients[,1]-summary(mod_true)$coefficients[,1])**2+(summary(mod_cc)$coefficients[,2])**2)
    fit_betas<-summary(mod_cc)$coefficients[,1]
    fit_se<-summary(mod_cc)$coefficients[,2]
    covprob<-cov_prob(fit_betas,fit_se,beta_true)
    
    #update results:
    rel.bias<-c(i,j,k,rel.bias)
    ratio_var<-c(i,j,k,ratio_var)
    rmse<-c(i,j,k,rmse)
    covprob<-c(i,j,k,covprob)
    results[[b]]<-rbind(results[[b]],rel.bias,ratio_var,rmse,covprob)
    # stats_impute[[count2]]<-list(b,i,j,k,tables(simdat_miss[cc_id,])); count2<-count2+1;
    
    
    results[[b]]<-results[[b]][-1,]  #remove null first row
    #FOR NEW RUNS UNCOMMENT BELOW:
    
  }	#next b	
  ###################################
  #    SET UP OUTPUT MATRICES       #
  ###################################
  ## create 48 matrices as a list:
  MaxIter=length(results)
  output<-list()
  count<-1
  stats<-c("rel.bias","ratio_var","rmse","covprob")
  upper<-ceiling(.975*MaxIter)
  lower<-ceiling(.025*MaxIter)
  
  
  for(k in 1:5){
    #cycle through results to collect rows
    collect<-matrix(0,nrow=1,ncol=13); colnames(collect)<-colnames(results[[1]])
    # for(b in 1:MaxIter){
    for(b in 1:MaxIter){
      id<-which( results[[b]][,3]==k)
      collect<-rbind(collect,results[[b]][id,])
    }#nextb
    collect<-collect[-1,]
    
    
    output[[count]]<-matrix(0,nrow=1,ncol=13)
    
    
    for(e in 1:3){
      mean<-UCL<-LCL<<-c()
      id<-which(rownames(collect)==stats[e])
      for(d in 4:ncol(collect)){
        
        mean<-c(mean,median(collect[id,d]))
        UCL<-c(UCL,sort(collect[id,d])[upper])
        LCL<-c(LCL,sort(collect[id,d])[lower])
        
      }#next d
      mean<-c(i,j,k,mean)
      UCL<-c(i,j,k,UCL)
      LCL<-c(i,j,k,LCL)
      output[[count]]<-rbind(output[[count]],mean,UCL,LCL)
      idx<-nrow(output[[count]])
      rownames(output[[count]])[(idx-2):idx]<-c(paste(stats[e],"-median",sep=""),paste(stats[e],"-UCL",sep=""),
                                                paste(stats[e],"-LCL",sep=""))
    }#next e
    id<-which(rownames(collect)==stats[4])
    CP<-c(i,j,k,colSums(collect[id,4:13])/length(id))
    output[[count]]<-rbind(output[[count]],CP)
    idx<-nrow(output[[count]])
    rownames(output[[count]])[idx]<-"CovProb-pct"
    output[[count]]<-output[[count]][-1,]	
    #labels<-c(rep(c("median","UCL","LCL"),3),"pct")
    #output[[count]]<-as.data.frame(cbind(labels,output[[count]]))
    colnames(output[[count]])<-c( colnames(results[[b]]))
    count=count+1
    
  }#nextk
  return(list(output,prob.miss, odds.miss))
  
} #END FUNCTION