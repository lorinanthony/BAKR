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

### Set the random seed to reproduce research ###
set.seed(111590)

### Call the available cores accross the computer for parallelization ###
registerDoParallel(cores = detectCores())

### Set up simulation parameters ###
nidv = 500; nsnp = 2e3; ncausal = 25; pve = 0.4; rho = 1; k = 500 # number of data points

### Generate the data ###
X = glSim(nidv,nsnp-ncausal,ncausal,parallel = TRUE,k = k, LD = FALSE, alpha = 0.3,pop.freq = rep(1/k,k))
X  = as.matrix(X)-1
colnames(X) = paste("SNP",1:ncol(X),sep = "")

### Select causal variables ###
s=(nsnp-ncausal+1):ncol(X)
Xcausal=X[,s]

### Generate the effects ###
beta1=rnorm(ncausal)
y_add=Xcausal%*%beta1
beta1 = beta1*sqrt(pve/var(y_add))
y_add=Xcausal%*%beta1

y_err=rnorm(nidv) #Random Error
y_err=y_err*sqrt((1-pve)/var(y_err))

y=y_add+y_err #Generate Responses
y=scale(y) #Standardize Responses

### Create the training and test sets ###
train.idc = sample(1:length(y), size=0.8*length(y), replace=FALSE)

X_train = X[train.idc,]; y_train = as.numeric(y[train.idc])
X_test = X[-train.idc,]; y_test = as.numeric(y[-train.idc])

######################################################################################
######################################################################################
######################################################################################

### Running BAKR ###

### Create the Approximate Kernel ###
n = dim(X_train)[1]; p = dim(X_train)[2]
G_tilde = GetApproxKernel(t(X_train),1e4)

v=matrix(1, n, 1)
M=diag(n)-v%*%t(v)/n
Gn=M%*%G_tilde%*%M
Gn=Gn/mean(diag(Gn))

evd = eigen(Gn)
### Choose the number of components that explain 99.995% of the cumulative variance in the data ###
explained_var = cumsum(evd$values/sum(evd$values))
q = 1:min(which(explained_var >= 0.75))
Lambda = diag(evd$values[q]^(-1)) # Matrix of Eigenvalues
U = evd$vectors[,q] # Unitary Matrix of Eigenvectors
### Define Inverse Mapping ###
B = InverseMap(t(X_train),U)

### Set up the number of MCMC samples and burn-ins ###
iter = 1e4
burn = 5e3

### Run BAKR ###
Gibbs = BAKRGibbs(U,y_train,Lambda,iter,burn)
beta.out = PostBeta(B,as.numeric(PostMean(Gibbs$alpha)))
BAKR_pred = BAKRPredict(t(X_test),beta.out)

### Get Diagnostics ###
MSPE_BAKR = mean((y_test-BAKR_pred)^2)

######################################################################################
######################################################################################
######################################################################################

### Bayesian Ridge Regression ###
ETA = list(list(X = X_train, model="BRR"))

### Run the Gibbs Sampler ###
reg.BRR = BGLR(y=y_train, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)

### Get the posterior of the missing variables ###
BRR_pred = X_test%*%reg.BRR$ETA[[1]]$b

### Get Diagnostics ###
MSPE_BRR = mean((y_test-BRR_pred)^2)

######################################################################################
######################################################################################
######################################################################################

### Bayesian LMM ###
K = GetLinearKernel(t(X_train))
ETA = list(list(K = K, model="RKHS"))

### Run the Gibbs Sampler ###
reg.BBLUP = BGLR(y=y_train, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)
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
reg.BL = BGLR(y=y_train, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)

### Get the posterior of the missing variables ###
BL_pred = X_test%*%reg.BL$ETA[[1]]$b

### Find MSE and Correlations ###
MSPE_BL = mean((y_test-BL_pred)^2)

######################################################################################
######################################################################################
######################################################################################

### SVM Model ###
reg.svm  = ksvm(y=y_train,x=X_train,type="eps-svr",kernel="rbfdot")
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
