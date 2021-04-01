library('glmnet')
library('missMDA')
source('slope_admm.R')



p.miss<-0.1;
kor<-0;

n<-100;
p<-100;

sig<-1.3*sqrt(2*log(p));

q<-0.1;
beta <-numeric(p)
M<-200;
kvec<-c(3,6,9,12,15);
len<-length(kvec);


SSE <- matrix(rep(0,M*len), nrow=M);
FP <- matrix(rep(0,M*len), nrow=M);
TP <- matrix(rep(0,M*len), nrow=M);
SSP <- matrix(rep(0,M*len), nrow=M);



lambda = qnorm(1-q*seq(1:p)/2/p);


Covar<-matrix(kor,p,p)
diag(Covar)=1
invCov=solve(Covar);
Ch=t(chol(Covar))

for (j in 1:len)
{
beta<-rep(0,p);
k<-kvec[j];
beta[1:k+10]<-sig;



for (i in 1:M)
{
print(i);
X=matrix(rnorm(n*p),n,p)%*% t(Ch);
X<-X/sqrt(n);
Y<-X%*%beta+rnorm(n);

X.obs<-X;
X.sim<-X;

if (p.miss>0)
{
patterns <- runif(n*p)< p.miss # missing completely at random
      X.obs[patterns] <- NA

X.sim = imputePCA(X.obs)$completeObs
}


obj3<-cv.glmnet(X.sim,Y,standardize=FALSE, intercept=FALSE)
betal<-coefficients(obj3, s='lambda.min');
betal<-betal[2:(p+1),1];
RSS<-sum((Y-X.sim%*%betal)^2)
indsup<-which(abs(betal)>0);
  Xsup<-X.sim[,indsup];
  betasup<-betal[indsup];
  l2<-length(indsup);
  sigma2<-sqrt(RSS/(n-l2))


list.SLOB = ABSLOPE(X.obs, Y, lambda, a=0.01*n, b= 0.1*n, beta.start=betal,maxit = 100,print_iter =TRUE,tol_em=1e-3,impute='PCA',sigma.known=NA,sigma.init=sigma2,scale = FALSE,method_na= 'lineq')
ind2<-which(list.SLOB$gamma.avg>1/2)
hatbeta<-list.SLOB$beta.new
SSE[i,j] = sum((hatbeta-beta)^2)/sum(beta^2)
TP[i,j] = sum(hatbeta!=0&beta!=0)
FP[i,j] = sum(hatbeta!=0&beta==0)
SSP[i,j] = sum((X%*%hatbeta-X%*%beta)^2)/sum((X%*%beta)^2)


save(SSE, TP, FP, SSP, file='Wei_weak_ind_sigma_unknown.Rdata')


}}





