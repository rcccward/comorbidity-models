##################################################
# SIMULATION OF WEIGHTED SENSITIVITY ANALYSIS    #
# FOLLOWING MICE-LR ON DATA WITH                 #
# MNAR RACE MISSINGINESS IMPOSED                 #
#                                                #
#REF: Carpenter JR, et al.(2007)                 #
#Sensitivity analysis after multiple imputation  #
#under missing at random: a weighting approach.  #
#Stat Methods Med Res.16(3):259-75.              #
##################################################
 
library (mice)
library(randomForest)
library(pROC)
library(Matrix)
library(glmnet)
library(lmtest)
library(miscTools)
 
#demonstration requires patient_dat dataset from GitHub
mod_true<-glm(as.factor(death) ~ age+ male+single+ race+rural+ meana1c+meanmpr, data=patient_dat , family="binomial"(link="logit"))
 
# 'TRUE' PARAM ESTIMATES FROM COMPLETE CASE ANALYSIS
beta_true<-summary(mod_true)$coefficients[5:7,1]#nhw is reference
var_true<-diag(vcov(mod_true))[5:7]
or.true<-exp(beta_true)
ucl.true<-or.true+1.96*sqrt(or.true**2*var_true)
lcl.true<-or.true-1.96*sqrt(or.true**2*var_true)
resultsTrue<-cbind(or.true,lcl.true,ucl.true)

 


output<-matrix(0,ncol=21,nrow=0)
colnames(output)<-c("j","k","n","beta-nhb","beta-hisp","beta-other",
                    "relbias-nhb","relbias-hisp","relbias-other",
                    "eff-nhb","eff-hisp","eff-other",
                    "rmse-nhb","rmse-hisp","rmse-other",
                    "covp-nhb","covp-hisp","covp-other",
                    "ORnhb","ORhisp","ORother")
# #############################   
#      COV PROB FUNCT         #
###############################
#determine coverage probability
#fit is the model fit from imputed data; mod_true is the model with 'true' betas
cov_prob<-function(fit_betas,fit_se, beta_true){ 
  UCL<-fit_betas+qt(.975,1000-9)*fit_se
  LCL<-fit_betas-qt(.975,1000-9)*fit_se
  cov_flag<-c(rep(0, length(fit_betas)))
  for(i in 1:length(fit_betas)){
    if(beta_true[i]<UCL[i] & beta_true[i]>LCL[i])cov_flag[i]<-1
  }
  return(cov_flag)
}#END FUNCTION
 

simdat_miss<-patient_dat
m=.3
N=5000
ncol=3
MissType="MNAR"

 
simdat_miss<- miss_gen(m,N,ncol,MissType,simdat_miss)      

#MICE WITH LOGISTIC IMPUTATION:
 
mice1<-mice(simdat_miss[,1:9],
            m = 50 ,
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
 
fit <- with(mice1, glm(as.factor(death) ~ age+ male+single+ race+rural+ meana1c+meanmpr, family="binomial"(link="logit")))
poolfit<-pool(fit)
imputed = apply(mice1$imp$race,2,as.numeric)#one column per imputation, one row per imputed value
 
impute_table=(apply(imputed,2,function(x) table(x)))
 
betas<-poolfit$qhat[,c(5:7)]
 
var<-cbind(poolfit$u[,5,5],poolfit$u[,6,6],poolfit$u[,7,7])
 
       
 #########################################
 # ITERATIVE PROCESS TO  COMPARE DELTAS  #
 # FOR A SINGLE SET OF MICE-LR RESULTS:  #
 #########################################
 # REQUIRES POOLFIT, THE OBJECT OUTPUT FROM A MICE-LR IMPUTATION RUN
 # SEE ABOVE ITERATIONS CODE
 output<-matrix(0,ncol=21,nrow=0)
 output.wt<-list()
 values<-seq(from=-1,to=1 ,by=.1)
 
 betas<-poolfit$qhat[,c(5:7)]
 
 
 var<-cbind(poolfit$u[,5,5],poolfit$u[,6,6],poolfit$u[,7,7])
 
 count=1
 for(j in -2:2 ){#this is the non-Hispanic black adjustment
   for(k in -3.5:3.5){ #this is the Hispanic adjustment
     for(n in -1:0){ #this is the "other" racial group adjustment
       
       delta=c(j,k,n)
       # DELTA VALUES ARE THE ADJUSTMENTS DETERMINED BY HYPOTHESIS,
       # OR BY ITERATIVE PROCESS (BELOW)
       # SEE EQ. (1) IN THE REF. FOR DEFINTION OF DELTA
       delta<-matrix(rep(delta,mice1$m),byrow=F,ncol=mice1$m)
       delta_table= (delta)*impute_table[-1,]
       #center the values
       #change to columns:
       delta_table<-t(delta_table)
       delta_table<-t(scale(delta_table, center=T,scale=F))
       delta_exp=exp(-(delta_table))
       
       weights<-apply(delta_exp,1,function(x) x/sum(x))
       output.wt[[count]]<-weights
       
       count=count+1
       
       
       betas_mnar=colSums(weights*betas)
       
       w_mnar=colSums(weights*var)
       b_mnar=colSums(weights*(betas-matrix(rep(betas_mnar,mice1$m),byrow=T,nrow=mice1$m))**2)
       t_mnar=(1+1/nrow(betas))*b_mnar +w_mnar
       rel.bias<-1-(betas_mnar/beta_true)
       ratio_var<-t_mnar/(summary(mod_true)$coefficients[5:7,2])**2
       rmse<-sqrt((betas_mnar-beta_true)**2+t_mnar)
       covprob<-cov_prob(betas_mnar,sqrt(t_mnar),beta_true)
       OR_mnar=exp(betas_mnar)
       compare= c(delta[,1],betas_mnar,rel.bias,ratio_var,rmse,covprob,OR_mnar)
       output<-rbind(output,compare)
     }
   }
 }
 
 colnames(output)<-c("j","k","n","beta-nhb","beta-hisp","beta-other",
                     "relbias-nhb","relbias-hisp","relbias-other",
                     "eff-nhb","eff-hisp","eff-other",
                     "rmse-nhb","rmse-hisp","rmse-other",
                     "covp-nhb","covp-hisp","covp-other",
                     "ORnhb","ORhisp","ORother")
 rownames(output)<-c(rep("",nrow(output)))
 