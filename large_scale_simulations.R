#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE) 

#install.packages("Rcpp")
#install.packages("RcppArmadillo")
source("SLOBEC_BH.R")
library(glmnet)
library(mvtnorm)
library(dplyr)
library(SSLASSO)
library(missMDA)
library(SLOPE)

#### Parsing arguments

args <- as.numeric(args)

p <- args[1] #number of features

n <- args[2] #number of observations

k <- args[3] # number of significant features

pmis <- args[4] # probability of missingness

sig  <- args[5] # strength of the signal

rho  <- args[6] # correlation

#### Number of replications

m <- 5



### Creataing data frame with results

results <- data.frame(
			type=numeric(14*m),
			n=n,
			p=p,
			k=k,
			power=numeric(14*m),
			FDR=numeric(14*m),
			err_sigma=numeric(14*m),
			error_beta=numeric(14*m),
			MSE=numeric(14*m),
			pred_err=numeric(14*m)
)

#### Setting parameters of other methods

 lambda_alasso <- qnorm(1-0.1/2/p); #Penalty coefficient for adaptive lasso

### Main loop
for( i in 1:m){

  set.seed(i)

#### Generating data
 
  sigma <- matrix(rho,p,p) 
  diag(sigma) <- 1
  X <- (rmvnorm(n,numeric(p),sigma))/sqrt(n); # generating X and normalizing
  
  beta <-numeric(p)

  c <- sig*sqrt(2*log(p))
  beta[21:(k+20)] <- c
  # Centering design matrix 
   X <- t(t(X)- colMeans(X))

   Y <- X%*%beta+rnorm(n)
   # Centering observation vector
   Y <- Y - mean(Y)

   ### Generating missing data
   Xmis <- X
   if (pmis>0){

	misspatern <- matrix(rbinom(n*p,1,pmis),n,p)
	for(l in 1:n)
  		for(j in 1:p)
 			 {
    				if(misspatern[l,j]==1)
    					
    					  Xmis[l,j] <- NA
 			 }
### Imputation via PCA 			 
	X.simPCA <- imputePCA(Xmis)$completeObs

	} else {
  		
  	 	X.simPCA <- X
 }


### LASSO CV
	
	obj_LASSO_CV <- cv.glmnet(X.simPCA,Y,standardize=FALSE, intercept=FALSE)

	beta_lassoCV <- as.numeric(coef(obj_LASSO_CV,s="lambda.min"))[-1]

#### Estimating sigma by LASSO CV

    RSS <- sum((Y- X.simPCA%*%beta_lassoCV)^2)


    indsup <- which(abs(beta_lassoCV)>0)
    Xsup <- X.simPCA[,indsup];
    betasup <- beta_lassoCV[indsup];
    l2 <- length(indsup);
    sigma_lassoCV <- sqrt(RSS/(n-l2))


### SLOBE 
obj_slobe_sigma_unknown <- SLOBEC( beta_lassoCV,Xmis,X.simPCA,Y,a_prior = 0.01*n, b_prior = .1*n,FDR=0.1, max_iter=20,BH=TRUE,verbose = F,known_sigma = F,known_cov = F)  
beta_slobe_sigma_unknown <- obj_slobe_sigma_unknown$beta

obj_slobe_sigma_lassoCV <- SLOBEC( beta_lassoCV,Xmis,X.simPCA,Y,a_prior = 0.01*n, b_prior = .1*n,FDR=0.1, max_iter=20,BH=TRUE,verbose = F,sigma = sigma_lassoCV,known_sigma = T,known_cov = F)  
beta_slobe_sigma_lassoCV <- obj_slobe_sigma_lassoCV$beta

obj_slobe_sigma_known <- SLOBEC( beta_lassoCV,Xmis,X.simPCA,Y,a_prior = 0.01*n, b_prior = .1*n,FDR=0.1, max_iter=20,BH=TRUE,verbose = F,sigma=1,known_sigma = T,known_cov = T,Covmat = sigma/n)  
beta_slobe_sigma_known <- obj_slobe_sigma_known$beta


### ABLASSO
obj_ablasso_sigma_unknown <- SLOBEC( beta_lassoCV,Xmis,X.simPCA,Y,a_prior = 0.01*n, b_prior = .1*n,FDR=0.1, max_iter=20,BH=F,verbose = F,known_sigma = F,known_cov = F)  
beta_ablasso_sigma_unknown <- obj_ablasso_sigma_unknown$beta

obj_ablasso_sigma_lassoCV <- SLOBEC( beta_lassoCV,Xmis,X.simPCA,Y,a_prior = 0.01*n, b_prior = .1*n,FDR=0.1, max_iter=20,BH=F,verbose = F,sigma = sigma_lassoCV,known_sigma = T,known_cov = F)  
beta_ablasso_sigma_lassoCV <- obj_ablasso_sigma_lassoCV$beta

obj_ablasso_sigma_known <- SLOBEC( beta_lassoCV,Xmis,X.simPCA,Y,a_prior = 0.01*n, b_prior = .1*n,FDR=0.1, max_iter=20,BH=F,verbose = F,sigma=1,known_sigma = T,known_cov = T,Covmat = sigma/n)  
beta_ablasso_sigma_known <- obj_ablasso_sigma_known$beta


#### Spike and Slab Lasso

obj_SSLASSO_sigma_unknown <- SSLASSO(X.simPCA, Y, nlambda = 100, variance = "unknown",sigma = sigma_lassoCV)
beta_SSLASSO_sigma_unknown <- obj_SSLASSO_sigma_unknown$beta[,ncol(obj_SSLASSO_sigma_unknown$beta)]


obj_SSLASSO_sigma_known <- SSLASSO(X.simPCA, Y, nlambda = 100, variance = "fixed",sigma = 1)
beta_SSLASSO_sigma_known <- obj_SSLASSO_sigma_known$beta[,ncol(obj_SSLASSO_sigma_known$beta)]



obj_SSLASSO_sigma_lassoCV <- SSLASSO(X.simPCA, Y, nlambda = 100, variance = "fixed",sigma = sigma_lassoCV)
beta_SSLASSO_sigma_lassoCV <- obj_SSLASSO_sigma_lassoCV$beta[,ncol(obj_SSLASSO_sigma_lassoCV$beta)]



##### Adaptive Lasso



  W2<-abs(betasup)/sigma_lassoCV;
  Xtemp2<-sweep(Xsup,2,W2,'*');

 B=glmnet(Xtemp2,Y,intercept=FALSE,standardize=FALSE,lambda=sigma_lassoCV*lambda_alasso/n);
  betahatsup2=coefficients(B)[2:(l2+1)] * W2;
  beta_ALASSO_unknown<-rep(0,p);
  beta_ALASSO_unknown[indsup]<-betahatsup2;

  W2<-abs(betasup);
  Xtemp2<-sweep(Xsup,2,W2,'*');


B=glmnet(Xtemp2,Y,intercept=FALSE,standardize=FALSE,lambda=lambda_alasso/n);
  betahatsup2=coefficients(B)[2:(l2+1)] * W2;
  beta_ALASSO_known<-rep(0,p);
  beta_ALASSO_known[indsup]<-betahatsup2;


#### SLOPE ####
obj_SLOPE_unknown<-SLOPE(X.simPCA,Y,q=0.1, alpha=1/sqrt(n)*sigma_lassoCV, lambda='gaussian', solver='admm',max_passes=100)
beta_SLOPE_unknown<-rep(0,p);

if(length(obj_SLOPE_unknown$nonzeros)>0){
Xn<-X[,obj_SLOPE_unknown$nonzeros];
obj2<-lm(Y~Xn-1);
beta_SLOPE_unknown[obj_SLOPE_unknown$nonzeros]<-obj2$coefficients;
}


obj_SLOPE_known<-SLOPE(X.simPCA,Y,q=0.1, alpha=1/sqrt(n), lambda='gaussian', solver='admm',max_passes=100)
beta_SLOPE_known<-rep(0,p);

if(length(obj_SLOPE_known$nonzeros)>0){
Xn<-X[,obj_SLOPE_known$nonzeros];
obj2<-lm(Y~Xn-1);
beta_SLOPE_known[obj_SLOPE_known$nonzeros]<-obj2$coefficients;
}




#### writing results
results$type[i]="BH unknown"
results$type[i+m]="BH unknown LASSO_cv"
results$type[i+2*m]="BH known"
results$type[i+3*m]="Bonferoni unknown"
results$type[i+4*m]="Bonferoni unknown LASSO_cv"
results$type[i+5*m]="Bonferoni known"
results$type[i+6*m]="ALASSO known"
results$type[i+7*m]="ALASSO unknown"
results$type[i+8*m]="SLOPE unknown"
results$type[i+9*m]="SLOPE known"
results$type[i+10*m]="LASSO CV"
results$type[i+11*m]="SSL known lasso CV "
results$type[i+12*m]="SSL known"
results$type[i+13*m]="SSL unknown"

results$power[i] = sum(beta!=0 &beta_slobe_sigma_unknown!=0 )/k
results$power[i+m] = sum(beta!=0 &beta_slobe_sigma_lassoCV!=0 )/k
results$power[i+2*m] = sum(beta!=0 &beta_slobe_sigma_known!=0 )/k
results$power[i+3*m] = sum(beta!=0 &beta_ablasso_sigma_unknown!=0 )/k
results$power[i+4*m] = sum(beta!=0 &beta_ablasso_sigma_lassoCV!=0 )/k
results$power[i+5*m] = sum(beta!=0 &beta_ablasso_sigma_known!=0 )/k
results$power[i+6*m] = sum(beta!=0 &beta_ALASSO_known!=0 )/k
results$power[i+7*m] = sum(beta!=0 &beta_ALASSO_unknown!=0 )/k
results$power[i+8*m] = sum(beta!=0 &beta_SLOPE_unknown!=0 )/k
results$power[i+9*m] = sum(beta!=0 &beta_SLOPE_known!=0 )/k
results$power[i+10*m] = sum(beta!=0 &beta_lassoCV!=0 )/k
results$power[i+11*m] = sum(beta!=0 &beta_SSLASSO_sigma_lassoCV!=0 )/k
results$power[i+12*m] = sum(beta!=0 &beta_SSLASSO_sigma_known!=0 )/k
results$power[i+13*m] = sum(beta!=0 &beta_SSLASSO_sigma_unknown!=0 )/k

results$FDR[i] = sum(beta==0&beta_slobe_sigma_unknown!=0)/(max(1,sum(beta_slobe_sigma_unknown!=0)))
results$FDR[i+m] = sum(beta==0&beta_slobe_sigma_lassoCV!=0)/(max(1,sum(beta_slobe_sigma_lassoCV!=0)))
results$FDR[i+2*m] = sum(beta==0&beta_slobe_sigma_known!=0)/(max(1,sum(beta_slobe_sigma_known!=0)))
results$FDR[i+3*m] = sum(beta==0&beta_ablasso_sigma_unknown!=0)/(max(1,sum(beta_ablasso_sigma_unknown!=0)))
results$FDR[i+4*m] = sum(beta==0&beta_ablasso_sigma_lassoCV!=0)/(max(1,sum(beta_ablasso_sigma_lassoCV!=0)))
results$FDR[i+5*m] = sum(beta==0&beta_ablasso_sigma_known!=0)/(max(1,sum(beta_ablasso_sigma_known!=0)))
results$FDR[i+6*m] = sum(beta==0&beta_ALASSO_known!=0)/(max(1,sum(beta_ALASSO_known!=0)))
results$FDR[i+7*m] = sum(beta==0&beta_ALASSO_unknown!=0)/(max(1,sum(beta_ALASSO_unknown!=0)))
results$FDR[i+8*m] = sum(beta==0&beta_SLOPE_unknown!=0)/(max(1,sum(beta_SLOPE_unknown!=0)))
results$FDR[i+9*m] = sum(beta==0&beta_SLOPE_known!=0)/(max(1,sum(beta_SLOPE_known!=0)))
results$FDR[i+10*m] = sum(beta==0&beta_lassoCV!=0)/(max(1,sum(beta_lassoCV!=0)))
results$FDR[i+11*m] = sum(beta==0&beta_SSLASSO_sigma_lassoCV!=0)/(max(1,sum(beta_SSLASSO_sigma_lassoCV!=0)))
results$FDR[i+12*m] = sum(beta==0&beta_SSLASSO_sigma_known!=0)/(max(1,sum(beta_SSLASSO_sigma_known!=0)))
results$FDR[i+13*m] = sum(beta==0&beta_SSLASSO_sigma_unknown!=0)/(max(1,sum(beta_SSLASSO_sigma_unknown!=0)))


results$error_beta[i]= sum(abs(beta-beta_slobe_sigma_unknown))/sum(abs(beta))
results$error_beta[i+m]= sum(abs(beta-beta_slobe_sigma_lassoCV))/sum(abs(beta))
results$error_beta[i+2*m]= sum(abs(beta-beta_slobe_sigma_known))/sum(abs(beta))
results$error_beta[i+3*m]= sum(abs(beta-beta_ablasso_sigma_unknown))/sum(abs(beta))
results$error_beta[i+4*m]= sum(abs(beta-beta_ablasso_sigma_lassoCV))/sum(abs(beta))
results$error_beta[i+5*m]= sum(abs(beta-beta_ablasso_sigma_known))/sum(abs(beta))
results$error_beta[i+6*m]= sum(abs(beta- beta_ALASSO_known))/sum(abs(beta))
results$error_beta[i+7*m]= sum(abs(beta- beta_ALASSO_unknown))/sum(abs(beta))
results$error_beta[i+8*m]= sum(abs(beta- beta_SLOPE_unknown))/sum(abs(beta))
results$error_beta[i+9*m]= sum(abs(beta-beta_SLOPE_known))/sum(abs(beta))
results$error_beta[i+10*m]= sum(abs(beta- beta_lassoCV))/sum(abs(beta))
results$error_beta[i+11*m]= sum(abs(beta- beta_SSLASSO_sigma_lassoCV))/sum(abs(beta))
results$error_beta[i+12*m]= sum(abs(beta- beta_SSLASSO_sigma_known))/sum(abs(beta))
results$error_beta[i+13*m]= sum(abs(beta-beta_SSLASSO_sigma_unknown))/sum(abs(beta))



results$err_sigma[i] = 1- obj_slobe_sigma_unknown$sigma
results$err_sigma[i+m] = 1- obj_slobe_sigma_lassoCV$sigma
results$err_sigma[i+2*m] = 1-  obj_slobe_sigma_known$sigma
results$err_sigma[i+3*m] = 1- obj_ablasso_sigma_unknown$sigma
results$err_sigma[i+4*m] = 1- obj_ablasso_sigma_lassoCV$sigma
results$err_sigma[i+5*m] =1-  obj_ablasso_sigma_known$sigma
results$err_sigma[i+10*m] = 1- sigma_lassoCV
results$err_sigma[i+13*m] = 1- tail(obj_SSLASSO_sigma_unknown$sigma,1) 




results$MSE[i]=sum((beta-beta_slobe_sigma_unknown)^2)/sum(beta^2)
results$MSE[i+m]=sum((beta-beta_slobe_sigma_lassoCV)^2)/sum(beta^2)
results$MSE[i+2*m]=sum((beta-beta_slobe_sigma_known)^2)/sum(beta^2)
results$MSE[i+3*m]=sum((beta- beta_ablasso_sigma_unknown)^2)/sum(beta^2)
results$MSE[i+4*m]=sum((beta-beta_ablasso_sigma_lassoCV)^2)/sum(beta^2)
results$MSE[i+5*m]=sum((beta- beta_ablasso_sigma_known)^2)/sum(beta^2)
results$MSE[i+6*m]=sum((beta- beta_ALASSO_known)^2)/sum(beta^2)
results$MSE[i+7*m]=sum((beta- beta_ALASSO_unknown)^2)/sum(beta^2)
results$MSE[i+8*m]=sum((beta- beta_SLOPE_unknown)^2)/sum(beta^2)
results$MSE[i+9*m]=sum((beta-beta_SLOPE_known)^2)/sum(beta^2)
results$MSE[i+10*m]=sum((beta- beta_lassoCV)^2)/sum(beta^2)
results$MSE[i+11*m]=sum((beta-beta_SSLASSO_sigma_lassoCV)^2)/sum(beta^2)
results$MSE[i+12*m]=sum((beta- beta_SSLASSO_sigma_known)^2)/sum(beta^2)
results$MSE[i+13*m]=sum((beta-beta_SSLASSO_sigma_unknown)^2)/sum(beta^2)




results$pred_err[i] = sum((X%*%beta-X%*%beta_slobe_sigma_unknown)^2)/sum((X%*%beta)^2)
results$pred_err[i+m] = sum((X%*%beta-X%*%beta_slobe_sigma_lassoCV)^2)/sum((X%*%beta)^2)
results$pred_err[i+2*m] = sum((X%*%beta- X%*%beta_slobe_sigma_known)^2)/sum((X%*%beta)^2)
results$pred_err[i+3*m] = sum((X%*%beta- X%*%beta_ablasso_sigma_unknown)^2)/sum((X%*%beta)^2)
results$pred_err[i+4*m] = sum((X%*%beta- X%*%beta_ablasso_sigma_lassoCV)^2)/sum((X%*%beta)^2)
results$pred_err[i+5*m] = sum((X%*%beta- X%*%beta_ablasso_sigma_known)^2)/sum((X%*%beta)^2)
results$pred_err[i+6*m] = sum((X%*%beta- X%*%beta_ALASSO_known)^2)/sum((X%*%beta)^2)
results$pred_err[i+7*m] = sum((X%*%beta- X%*%beta_ALASSO_unknown)^2)/sum((X%*%beta)^2)
results$pred_err[i+8*m] = sum((X%*%beta-X%*%beta_SLOPE_unknown)^2)/sum((X%*%beta)^2)
results$pred_err[i+9*m] = sum((X%*%beta- X%*%beta_SLOPE_known)^2)/sum((X%*%beta)^2)
results$pred_err[i+10*m] = sum((X%*%beta- X%*%beta_lassoCV)^2)/sum((X%*%beta)^2)
results$pred_err[i+11*m] = sum((X%*%beta- X%*%beta_SSLASSO_sigma_lassoCV)^2)/sum((X%*%beta)^2)
results$pred_err[i+12*m] = sum((X%*%beta-X%*%beta_SSLASSO_sigma_known)^2)/sum((X%*%beta)^2)
results$pred_err[i+13*m] = sum((X%*%beta- X%*%beta_SSLASSO_sigma_unknown)^2)/sum((X%*%beta)^2)

print(i)

}

print(results%>%group_by(type)%>%summarise(FDR=mean(FDR),MSE=mean(MSE),power=mean(power),prediction =mean(pred_err)))
save(file=paste("s_", sig, "cor_",rho,"p_",p,"k_",k,"miss_",pmis,".Rdata",sep="",collapse = ""),results)

