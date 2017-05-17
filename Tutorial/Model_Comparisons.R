### Clear Console ###
cat("\014")

### Clear Environment ### 
rm(list = ls(all = TRUE))

### Load in the R libraries ###
library(coda)
library(MASS)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(BGLR)
library(monomvn)
library(kernlab)
library(adegenet)

### Load in the C++ functions from working directory ###
sourceCpp("BAKRGibbs.cpp")

######################################################################################
######################################################################################
######################################################################################

### Set the Seed for the analysis ###
set.seed(11151990)

ind = 500; nsnp = 2e3

# simulation parameters
pve=0.6; rho=0.2;
ncausal1= 10 #Set 1 of causal SNPs 
ncausal2 = 15 #Set 2 of Causal SNPs
ncausal3 = 50-(ncausal1+ncausal2) #Set 3 of Causal SNPs with only marginal effects

ncausal = ncausal1+ncausal2+ncausal3

X = glSim(ind,nsnp-ncausal,ncausal,parallel = TRUE, LD = TRUE)
X  = as.matrix(X)-1
s=(nsnp-ncausal+1):ncol(X)

#Select Causal SNPs
s1=sample(s, ncausal1, replace=F)
s2=sample(s[-s1], ncausal2, replace=F)
s3=sample(s[-c(s1,s2)], ncausal3, replace=F)

# Generate the ground-truth regression coefficients for the variables
# (X). Adjust the effects so that
# the variables (SNPs) explain x percent of the variance in the
# outcome.
Xcausal1=X[,s1]; Xcausal2=X[,s2]; Xcausal3=X[,s3]
Xepi=c()
for(i in 1:ncausal1){
  Xepi=cbind(Xepi,Xcausal1[,i]*Xcausal2)
}
dim(Xepi)

# Marginal Effects Only
Xmarginal=Xcausal3
beta=rnorm(dim(Xmarginal)[2])
y_marginal=Xmarginal%*%beta
beta=beta*sqrt(pve*rho/var(y_marginal))
y_marginal=Xmarginal%*%beta

#Pairwise Epistatic Effects
beta=rnorm(dim(Xepi)[2])
y_epi=Xepi%*%beta
beta=beta*sqrt(pve*(1-rho)/var(y_epi))
y_epi=Xepi%*%beta

# error
y_err=rnorm(ind)
y_err=y_err*sqrt((1-pve)/var(y_err))

y=y_marginal+y_epi+y_err
y=scale(y)

### Check dimensions and add SNP names ###
dim(X); dim(y)
colnames(X) = paste("SNP",1:ncol(X),sep="")
SNPs = colnames(X)[c(s1,s2)]

### Create the training and test sets ###
train.idc = sample(1:length(y), size=0.8*length(y), replace=FALSE)

X_train = X[train.idc,]; y_train = as.numeric(y[train.idc])
X_test = X[-train.idc,]; y_test = as.numeric(y[-train.idc])

######################################################################################
######################################################################################
######################################################################################

### Defined extra parameters needed to run the analysis ###
n = dim(X_train)[1] #Sample size
p = dim(X_train)[2] #Number of markers or genes

### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
Kn = ApproxGaussKernel(t(X_train),p,p)

### Center and Scale K_tilde ###
v=matrix(1, n, 1)
M=diag(n)-v%*%t(v)/n
Kn=M%*%Kn%*%M
Kn=Kn/mean(diag(Kn))

### Find the Eigenvalue Decomposition of K ###
evd = EigDecomp(Kn)

### Truncate the data based on the desired cumulative variance explained ###
explained_var = cumsum(evd$lambda/sum(evd$lambda))
q = 1:min(which(explained_var >= 0.99))
Lambda = diag(sort(evd$lambda,decreasing = TRUE)[q]^(-1)) # Matrix of Eigenvalues
U = evd$U[,q] # Unitary Matrix of Eigenvectors

### Define Inverse Mapping ###
B = InverseMap(t(X_train),U)

### Set up the number of MCMC samples and burn-ins ###
mcmc.iter = 2e3
mcmc.burn = 1e3

### Run BAKR ### 
Gibbs = BAKRGibbs(U,y_train,Lambda,mcmc.iter,mcmc.burn)

### Look at the Posterior Summaries ###
theta.out = PostMean(Gibbs$theta)
beta.out = PostBeta(B,theta.out); names(beta.out) = colnames(X)
BAKR_pred = X_test%*%beta.out

### Get Diagnostics ###
MSPE_BAKR = mean((y_test-BAKR_pred)^2)

######################################################################################
######################################################################################
######################################################################################

### Bayesian Ridge Regression ###
ETA = list(list(X = X_train, model="BRR"))

### Run the Gibbs Sampler ###
reg.BRR = BGLR(y=y_train, ETA=ETA, nIter=mcmc.iter, burnIn=mcmc.burn, verbose=FALSE)

### Get the posterior of the missing variables ###
BRR_pred = X_test%*%reg.BRR$ETA[[1]]$b

### Get Diagnostics ###
MSPE_BRR = mean((y_test-BRR_pred)^2)

######################################################################################
######################################################################################
######################################################################################

### Bayesian BLUP ###
K = GetLinearKernel(t(X_train))
ETA = list(list(K = K, model="RKHS"))

### Run the Gibbs Sampler ###
reg.BBLUP = BGLR(y=y_train, ETA=ETA, nIter=mcmc.iter, burnIn=mcmc.burn, verbose=FALSE)
reg.BBLUP_b = ginv(X_train)%*%reg.BBLUP$ETA[[1]]$u

### Get the posterior of the missing variables ###
BBLUP_pred = X_test%*%reg.BBLUP_b

### Get Diagnostics ###
MSPE_BBLUP = mean((y_test-BBLUP_pred)^2)

######################################################################################
######################################################################################
######################################################################################

### Set up List ###
ETA = list(list(X = X_train, model="BL"))

### Run the Gibbs Sampler ###
reg.BL = BGLR(y=y_train, ETA=ETA, nIter=mcmc.iter, burnIn=mcmc.burn, verbose=FALSE)

### Get the posterior of the missing variables ###
BL_pred = X_test%*%reg.BL$ETA[[1]]$b

### Find MSE and Correlations ###
MSPE_BL = mean((y_test-BL_pred)^2)

######################################################################################
######################################################################################
######################################################################################

### SVM Model ###
reg.svm  = ksvm(y=y_train,x=X_train,type="nu-svr",kernel="rbfdot")
SVM_pred = predict(reg.svm, X_test, type="response")

### Get Diagnostics ###
MSPE_SVM = mean((y_test-SVM_pred)^2)

######################################################################################
######################################################################################
######################################################################################

### Check the Results ###    
scores = c(MSPE_BRR,MSPE_BL,MSPE_BBLUP,MSPE_SVM,MSPE_BAKR)
names(scores) = c("BRR","BL","BBLUP","SVM","BAKR")

c(scores,names(scores)[which(scores == min(scores))])