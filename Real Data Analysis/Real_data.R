

library('mice')
load("trauma_new.Rdata")
library("dplyr")

Xtrauma = data.matrix(trauma)
meanplaq<-mean(trauma$Plaquettes.org);
normplaq<-sqrt(sum(trauma$Plaquettes^2));
normplaqorg<-sqrt(sum(trauma$Plaquettes.org^2));
n1<-nrow(Xtrauma);
Xtrauma.scale<-scale(Xtrauma)/sqrt(n1-1);
attr(Xtrauma.scale,'scaled:scale')
#imp = mice(Xtrauma,m=1, printFlag = FALSE)
#Xfull.sim  = complete(imp)
#save(Xfull.sim, file="imputed_mice.Rdata")
load('imputed_mice.Rdata')
Xtrauma.scale[,11]<-Xfull.simMice[,11]

*****all data are scaled********


##################repeat 10 times - function 
```{r}
library(SLOPE)
library(glmnet)
library(bigstep)
library(randomForest)
source("SLOBEC_BH.R")
library('SSLASSO')



compare_method = function(trauma,Xfull, trauma.org, nb_seed,percent_test,interaction_if=FALSE){
  n1<-nrow(trauma);
  set.seed(nb_seed) 
  sample <- sample.int(n = nrow(trauma), size = floor(percent_test*nrow(trauma)), replace = F)
  trauma.train  <- trauma[-sample, -1]
  trauma.test <-Xfull[sample, -1]
  Plaquettes.train <- trauma[-sample, 1]
  Plaquettes.test <- trauma[sample, 1]
  Plaq.org.test<-trauma.org[sample, 17]
  if(interaction_if==TRUE){
  trauma.train = interact(trauma.train)
  trauma.test = interact(trauma.test)
}
  
  


X = as.matrix(trauma.train[,1:15])
y = Plaquettes.train
p = ncol(X)
n = nrow(X)
lambda = qnorm(1-0.1/2/p);
normplaqorg.test<-sqrt(sum(Plaq.org.test^2))


X=X/sqrt(n-1)*sqrt(n1-1);
imp = mice(X,m=1, printFlag = FALSE)
X.sim  = complete(imp)
X.sim<-as.matrix(X.sim)

obj3<-cv.glmnet(X.sim,y,standardize=FALSE, intercept=FALSE)
betal<-coefficients(obj3, s='lambda.min');
betal<-betal[2:(p+1),1];

RSS<-sum((y- X.sim%*%betal)^2)


indsup<-which(abs(betal)>0);
Xsup<-X.sim[,indsup];
betasup<-betal[indsup];
l2<-length(indsup);
sigma2<-sqrt(RSS/(n-l2))


res<-SLOBEC(betal,X,X.sim,y,a_prior = 0.01*n, b_prior = .1*n,FDR=0.1, max_iter=100,BH=TRUE,verbose = F,sigma = sigma2,known_sigma = T,known_cov = F,Covmat = diag(p))  
beta.ABSLOPE<-res$beta;


pred.test<-as.matrix(trauma.test[,1:15]) %*% beta.ABSLOPE/sqrt(n-1)*sqrt(n1-1);
y.test<-Plaquettes.test

MSE.ABSLOPE<-sqrt(sum((pred.test-y.test)^2))*normplaq/normplaqorg.test
  
res<-SLOBEC(betal,X,X.sim,y,a_prior = 0.01*n, b_prior = .1*n,FDR=0.1, max_iter=100,BH=FALSE,verbose = F,sigma = sigma2,known_sigma = T,known_cov = F,Covmat = diag(p))  
beta.ABLAS<-res$beta;


pred.test<-as.matrix(trauma.test[,1:15]) %*% beta.ABLAS/sqrt(n-1)*sqrt(n1-1);


MSE.ABLAS<-sqrt(sum((pred.test-y.test)^2))*normplaq/normplaqorg.test
  

list.SLOPE = SLOPE(X.sim,y,q=0.1, alpha=1/sqrt(n)*sigma2, lambda='gaussian', solver='admm',max_passes=100, intercept=FALSE)

modelsel = list.SLOPE$nonzeros;

beta.SLOPE<-rep(0,15);
coef.SLOPE = summary(lm(y ~ X.sim[,modelsel]-1))$coef[,1]
beta.SLOPE[modelsel]= coef.SLOPE

pred.test<-as.matrix(trauma.test[,1:15]) %*% beta.SLOPE/sqrt(n-1)*sqrt(n1-1);


MSE.SLOPE<-sqrt(sum((pred.test-y.test)^2))*normplaq/normplaqorg.test

beta.LASSO<-betal
pred.test<-as.matrix(trauma.test[,1:15]) %*% betal/sqrt(n-1)*sqrt(n1-1)
MSE.LASSO<-sqrt(sum((pred.test-y.test)^2))*normplaq/normplaqorg.test


 W2<-abs(betasup)/sigma2;
  Xtemp2<-sweep(Xsup,2,W2,'*');

 B=glmnet(Xtemp2,y,intercept=FALSE,standardize=FALSE,lambda=sigma2*lambda/n);
  betahatsup2=coefficients(B)[2:(l2+1)] * W2;
  beta.ALAS<-rep(0,p);
  beta.ALAS[indsup]<-betahatsup2;

pred.test<-as.matrix(trauma.test[,1:15]) %*% beta.ALAS/sqrt(n-1)*sqrt(n1-1);

MSE.ALAS<-sqrt(sum((pred.test-y.test)^2))*normplaq/normplaqorg.test


list.SSLASSO=SSLASSO(X.sim, y, nlambda = 100, variance = "unknown",sigma=sigma2)
beta.SSL = list.SSLASSO$beta[,ncol(list.SSLASSO$beta)]

pred.test<-as.matrix(trauma.test[,1:15]) %*% beta.SSL/sqrt(n-1)*sqrt(n1-1);

MSE.SSL<-sqrt(sum((pred.test-y.test)^2))*normplaq/normplaqorg.test


rf = randomForest(X.sim, y)
XRF<-trauma.test[,1:15]/sqrt(n-1)*sqrt(n1-1);
y.test.rf = predict(rf, XRF)
MSE.RF<-sqrt(sum((y.test.rf-y.test)^2))*normplaq/normplaqorg.test


return(list(MSE.ABSLOPE=MSE.ABSLOPE, MSE.ABLAS=MSE.ABLAS,MSE.SLOPE=MSE.SLOPE, MSE.LASSO=MSE.LASSO, MSE.ALAS=MSE.ALAS,MSE.SSL=MSE.SSL, MSE.RF=MSE.RF, 
            beta.ABSLOPE=paste(round(beta.ABSLOPE,digit=6),collapse=" "),beta.ABLAS=paste(round(beta.ABLAS,digit=6),collapse=" "),
            beta.SLOPE=paste(round(beta.SLOPE,digit=6),collapse=" "), beta.LASSO=paste(round(beta.LASSO,digit=6),collapse=" "),
            beta.ALAS=paste(round(beta.ALAS,digit=6),collapse=" "),
            beta.SSL=paste(round(beta.SSL,digit=6),collapse=" ")))
}


##########repeat 10 times###################

res = NULL
simu.list = seq(100,190,10)
for (simu in simu.list){
 print(simu)
 res <- rbind(res,unlist(compare_method(trauma=Xtrauma.scale,Xfull=Xfull.simMice,trauma.org=Xtrauma, nb_seed=simu,percent_test=0.3)))
}
res.beta = res[,8:13]
res = res[,1:7]
save(res,res.beta, file='resPrehop_final.RData')

#######plot MSP #####

df = data.frame(res)
indx <- sapply(df, is.factor)
df[indx] <- lapply(df[indx], function(x) as.numeric(as.character(x)))
df1<-matrix(as.numeric(as.matrix(df)),nrow=10);
colnames(df1)<-c("SLOBE","ABLAS","SLOPE","LASSO","ALAS","SSL", "RF");
order.df1 = sort.int(apply(df1, 2, median, na.rm = TRUE),decreasing = FALSE, index.return	=TRUE)$ix

df1 = df1[,order.df1]
par(cex.axis=0.7)
boxplot(df1,ylab = 'relative l2 error of prediction')


#####numbers of variables statistics#####

numvar<-matrix(rep(0,10*6),nrow=10);
df2 = as.character(res.beta[,1])
df2 <- data.frame( do.call( rbind, strsplit( df2, ' ' ) ) )

print('ABSLOPE: number of selected times in 10 repetetions:')
colSums(df2!=0)

mean(rowSums(df2!=0))
numvar[,1]=rowSums(df2!=0)

df2 = as.character(res.beta[,2])
df2 <- data.frame( do.call( rbind, strsplit( df2, ' ' ) ) )
print('ABLAS: number of selected times in 10 repetetions:')
colSums(df2!=0)

mean(rowSums(df2!=0))
numvar[,2]=rowSums(df2!=0)
 
df2 = as.character(res.beta[,3])
df2 <- data.frame( do.call( rbind, strsplit( df2, ' ' ) ) )
print('SLOPE: number of selected times in 10 repetetions:')
colSums(df2!=0)

mean(rowSums(df2!=0))
numvar[,3]=rowSums(df2!=0)

df2 = as.character(res.beta[,4])
df2 <- data.frame( do.call( rbind, strsplit( df2, ' ' ) ) )
print('LASSO: number of selected times in 10 repetetions:')
colSums(df2!=0)

mean(rowSums(df2!=0))
numvar[,4]=rowSums(df2!=0)

df2 = as.character(res.beta[,5])
df2 <- data.frame( do.call( rbind, strsplit( df2, ' ' ) ) )
print('ALAS: number of selected times in 10 repetetions:')
colSums(df2!=0)
mean(rowSums(df2!=0))
numvar[,5]=rowSums(df2!=0)

df2 = as.character(res.beta[,6])
df2 <- data.frame( do.call( rbind, strsplit( df2, ' ' ) ) )
print('SSL: number of selected times in 10 repetetions:')
colSums(df2!=0)

mean(rowSums(df2!=0))
boxplot(rowSums(df2!=0))
numvar[,6]=rowSums(df2!=0)

colnames(numvar)<-c("SLOBE","ABLAS","SLOPE","LASSO","ALAS", "SSL");
order.num = sort.int(apply(numvar, 2, median, na.rm = TRUE),decreasing = FALSE, index.return	=TRUE)$ix

numvar = numvar[,order.num]
par(cex.axis=0.7)
boxplot(numvar,ylab = 'number of selected variables', ylim=c(0,15))



#####full data analysis

Xfull.sim<-as.matrix(Xfull.simMice);
X.sim<-Xfull.sim[,2:16];
p<-ncol(X.sim);
n<-nrow(X.sim);
y<-Xtrauma.scale[,1];
X<-Xtrauma.scale[,2:16];

obj3<-cv.glmnet(X.sim,y,standardize=FALSE, intercept=FALSE)
betal<-coefficients(obj3, s='lambda.min');
betal<-betal[2:(p+1),1];
pred.lasso<-X.sim%*%betal
plot(pred.lasso~y);

RSS<-sum((y- X.sim%*%betal)^2)


indsup<-which(abs(betal)>0);
Xsup<-X.sim[,indsup];
betasup<-betal[indsup];
l2<-length(indsup);
sigma2<-sqrt(RSS/(n-l2))
Xtrauma.scale<-as.matrix(Xtrauma.scale);


res<-SLOBEC(betal,X,X.sim,y,a_prior = 0.01*n, b_prior = .1*n,FDR=0.1, max_iter=100,BH=TRUE,verbose = F,sigma = sigma2,known_sigma = T,known_cov = F,Covmat = diag(p))  
beta.ABSLOPE<-res$beta;


