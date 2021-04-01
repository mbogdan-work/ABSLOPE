#source('fct/fct_others.R')
#source('fct/ABSLOPE.R')
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
theme_set(theme_bw())
library(gbm)
#source('slope_admm.R')
#source('initialization.R')
# Several curves of convergence
#----------with NA-------------------------

# Convergence of ABSLOPE

maxit=100;

# beta
set.seed(10)
n=100
p=100
nspr=10
mu = rep(0,p)
Sigma = diag(p)
signallevel=3
p.miss = 0.1
amplitude = signallevel*sqrt(2*log(p))  # signal amplitude (for noise level = 1)
nonzero = c(5,15,25,35,45,55,65,75,85,95)
beta = amplitude * (1:p %in% nonzero)
beta[45]=1.5*sqrt(2*log(p));
beta[55]=5*sqrt(2*log(p));
lambda = qnorm(1-0.1*seq(1:p)/2/p)
sigma = 1

# Design matrix 1
set.seed(100)
# normal distribution
X1 <- matrix(rnorm(n*p), nrow=n)%*%chol(Sigma) + matrix(rep(mu,n), nrow=n, byrow = TRUE)
X1 = scale(X1)/sqrt(n)
# response vectors
Y = X1 %*% beta +  sigma*rnorm(n)
# Missing values
X.obs <- X1
patterns <- runif(n*p)< p.miss # missing completely at random
X.obs[patterns] <- NA
X.sim = imputePCA(X.obs)
X.sim = X.sim$completeObs

obj3<-cv.glmnet(X.sim,Y,standardize=FALSE, intercept=FALSE)
betal<-coefficients(obj3, s='lambda.min');
betal<-betal[2:(p+1),1];
RSS<-sum((Y-X.sim%*%betal)^2)
indsup<-which(abs(betal)>0);
  Xsup<-X.sim[,indsup];
  betasup<-betal[indsup];
  l2<-length(indsup);
  sigma2<-sqrt(RSS/(n-l2))

list.ABSLOPE1 = ABSLOPE(X.obs, Y, lambda, a=0.01*n, b= 0.1*n, beta.start=betal,maxit = 100,print_iter =TRUE,tol_em=1e-3,impute='PCA',sigma.known=NA,sigma.init=sigma2,scale = FALSE,method_na= 'lineq')

res<-SLOBE_ADMM_approx_missing(betal,X.obs,X.sim, Y,a_prior = 0.01*n, b_prior = 0.1*n,FDR=0.1,max_iter=100, known_sigma=F,  Covmat=diag(p)/n, known_cov=F, BH=TRUE)



beta1<-res$beta[c(45,55,56)];




# Design matrix 2
set.seed(150)
# normal distribution
X2 <- matrix(rnorm(n*p), nrow=n)%*%chol(Sigma) + matrix(rep(mu,n), nrow=n, byrow = TRUE)
X2 = scale(X2)/sqrt(n)
# response vectors
Y2 = X2 %*% beta +  sigma*rnorm(n)
# Missing values
X.obs <- X2
patterns <- runif(n*p)< p.miss # missing completely at random
X.obs[patterns] <- NA
X.sim = imputePCA(X.obs)
X.sim = X.sim$completeObs

obj3<-cv.glmnet(X.sim,Y2,standardize=FALSE, intercept=FALSE)
betal<-coefficients(obj3, s='lambda.min');
betal<-betal[2:(p+1),1];
RSS<-sum((Y2-X.sim%*%betal)^2)
indsup<-which(abs(betal)>0);
  Xsup<-X.sim[,indsup];
  betasup<-betal[indsup];
  l2<-length(indsup);
  sigma2<-sqrt(RSS/(n-l2))


list.ABSLOPE2 = ABSLOPE(X.obs, Y2, lambda, a=0.01*n, b= 0.1*n, beta.start=betal,maxit = 100,print_iter =TRUE,tol_em=1e-3,impute='PCA',sigma.known=NA,sigma.init=sigma2,scale = FALSE,method_na= 'lineq')

res<-SLOBE_ADMM_approx_missing(betal,X.obs,X.sim, Y2,a_prior = 0.01*n, b_prior = 0.1*n,FDR=0.1,max_iter=100, known_sigma=F,  Covmat=diag(p)/n, known_cov=F, BH=TRUE)

beta2<-res$beta[c(45,55,56)];





# Design matrix 3
set.seed(400)
# normal distribution
X3 <- matrix(rnorm(n*p), nrow=n)%*%chol(Sigma) + matrix(rep(mu,n), nrow=n, byrow = TRUE)
X3 = scale(X3)/sqrt(n)
# response vectors
Y3 = X3 %*% beta +  sigma*rnorm(n)
# Missing values
X.obs <- X3
patterns <- runif(n*p)< p.miss # missing completely at random
X.obs[patterns] <- NA

X.sim = imputePCA(X.obs)
X.sim = X.sim$completeObs

obj3<-cv.glmnet(X.sim,Y3,standardize=FALSE, intercept=FALSE)
betal<-coefficients(obj3, s='lambda.min');
betal<-betal[2:(p+1),1];
RSS<-sum((Y3-X.sim%*%betal)^2)
indsup<-which(abs(betal)>0);
  Xsup<-X.sim[,indsup];
  betasup<-betal[indsup];
  l2<-length(indsup);
  sigma2<-sqrt(RSS/(n-l2))


list.ABSLOPE3 = ABSLOPE(X.obs, Y3, lambda, a=0.01*n, b= 0.1*n, beta.start=betal,maxit = 100,print_iter =TRUE,tol_em=1e-3,impute='PCA',sigma.known=NA,sigma.init=sigma2,scale = FALSE,method_na= 'lineq')

res<-SLOBE_ADMM_approx_missing(betal,X.obs,X.sim, Y3,a_prior = 0.01*n, b_prior = 0.1*n,FDR=0.1,max_iter=100, known_sigma=F,  Covmat=diag(p)/n, known_cov=F, BH=TRUE)

beta3<-res$beta[c(45,55,56)];





df1 = t(list.ABSLOPE1$seqbeta) %>% as.data.frame() %>% dplyr::select(V45,V55,V56) %>% mutate(iteration=0:maxit)  %>%  melt(id.vars = 'iteration', variable.name = 'beta')

df1$beta <-factor(df1$beta,labels=c("V45"=paste0('\u03b2',45),
                                    "V55"=paste0('\u03b2',55),
                                    "V56"=paste0('\u03b2',56)))
df1$Label = as.factor("1st ABSLOPE")
df2 = t(list.ABSLOPE2$seqbeta) %>% as.data.frame() %>% dplyr::select(V45,V55,V56) %>% mutate(iteration=0:maxit)  %>%  melt(id.vars = 'iteration', variable.name = 'beta')
df2$beta <-factor(df2$beta,labels=c("V45"=paste0('\u03b2',45),
                                    "V55"=paste0('\u03b2',55),
                                    "V56"=paste0('\u03b2',56)))
df2$Label = as.factor("2nd ABSLOPE")
df3 = t(list.ABSLOPE3$seqbeta) %>% as.data.frame() %>% dplyr::select(V45,V55,V56) %>% mutate(iteration=0:maxit)  %>%  melt(id.vars = 'iteration', variable.name = 'beta')
df3$beta <-factor(df3$beta,labels=c("V45"=paste0('\u03b2',45),
                                    "V55"=paste0('\u03b2',55),
                                    "V56"=paste0('\u03b2',56)))
df3$Label = as.factor("3rd ABSLOPE")
truebeta <- data.frame(iteration=rep(0:maxit,3), beta = c(rep(paste0('\u03b2',45),maxit+1), rep(paste0('\u03b2',55),maxit+1), rep(paste0('\u03b2',56),maxit+1)), value = c(rep(beta[45],maxit+1), rep(beta[55],maxit+1),rep(beta[56],maxit+1)))
truebeta$Label = as.factor("True value")
SLOB1 <- data.frame(iteration=rep(0:maxit,3), beta = c(rep(paste0('\u03b2',45),maxit+1), rep(paste0('\u03b2',55),maxit+1), rep(paste0('\u03b2',56),maxit+1)), value = c(rep(beta1[1],maxit+1), rep(beta1[2],maxit+1),rep(beta1[3],maxit+1)))
SLOB1$Label = as.factor("1st SLOBE")
SLOB2 <- data.frame(iteration=rep(0:maxit,3), beta = c(rep(paste0('\u03b2',45),maxit+1), rep(paste0('\u03b2',55),maxit+1), rep(paste0('\u03b2',56),maxit+1)), value = c(rep(beta2[1],maxit+1), rep(beta2[2],maxit+1),rep(beta2[3],maxit+1)))
SLOB2$Label = as.factor("2nd SLOBE")
SLOB3 <- data.frame(iteration=rep(0:maxit,3), beta = c(rep(paste0('\u03b2',45),maxit+1), rep(paste0('\u03b2',55),maxit+1), rep(paste0('\u03b2',56),maxit+1)), value = c(rep(beta3[1],maxit+1), rep(beta3[2],maxit+1),rep(beta3[3],maxit+1)))
SLOB3$Label = as.factor("3rd SLOBE")




df <- rbind(truebeta,df1, df2, df3, SLOB1, SLOB2, SLOB3)
pdf("convergence_beta_Gosia.pdf", width = 8, height = 4,  paper='special')
df %>%  ggplot() + aes(iteration,value,linetype=Label, shape=Label, color=Label, group=Label)   +
  geom_line() +
  #geom_point() +
  facet_grid(. ~ beta) +
  ylab(expression(paste("Estimate of ",beta))) +
  #theme(strip.text = element_text(size=12), axis.title=element_text(size=14))+
  scale_linetype_manual(values=c(2,1,1,1,2,2,2)) +
  scale_color_manual(values=c(1,2,3,4,2,3,4))
dev.off()


