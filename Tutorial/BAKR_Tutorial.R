### Clear Console ###
cat("\014")

### Clear Environment ### 
rm(list = ls(all = TRUE))

### Load in the necessary R libraries ###
library(coda)
library(MASS)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(BGLR)

### Load in the BAKR C++ functions from working directory ###
sourceCpp("BAKRGibbs.cpp")

######################################################################################
######################################################################################
######################################################################################

### Set the random seed to reproduce research ###
set.seed(111590)

### Set up simulation parameters ###
n = 500; p = 2e3; pve=0.5; rho=0.75;

### The Number of Causal Variables ###
ncausal = 30 
ncausal1= 10 #Set 1 of causal SNPs 
ncausal2 = 10 #Set 2 of Causal SNPs
ncausal3 = ncausal-(ncausal1+ncausal2) #Set 3 of Causal SNPs with only marginal effects

### Generate the data ###
maf <- 0.05 + 0.45*runif(p)
X   <- (runif(n*p) < maf) + (runif(n*p) < maf)
X   <- matrix(as.double(X),n,p,byrow = TRUE)
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); Geno=t((t(X)-Xmean)/Xsd)
colnames(X) = paste("SNP",1:ncol(X),sep="")

### Select Causal SNPs ###
s=sample(1:p,ncausal,replace = FALSE)
s1=sample(s, ncausal1, replace=F)
s2=sample(s[s%in%s1==FALSE], ncausal2, replace=F)
s3=sample(s[s%in%c(s1,s2)==FALSE], ncausal3, replace=F)

### Generate the Marginal Effects ###
Xmarginal=X[,s]
beta=runif(dim(Xmarginal)[2])
y_marginal=c(Xmarginal%*%beta)
beta=beta*sqrt(pve*rho/var(y_marginal))
y_marginal=Xmarginal%*%beta

### Generate the Pairwise Interaction Matrix W ###
Xcausal1=X[,s1]; Xcausal2=X[,s2];
W=c()
for(i in 1:ncausal1){
  W=cbind(W,Xcausal1[,i]*Xcausal2)
}
dim(W)

### Generate the Epistatic Effects ###
gamma=runif(dim(W)[2])
y_epi=c(W%*%gamma)
gamma=gamma*sqrt(pve*(1-rho)/var(y_epi))
y_epi=W%*%gamma

### Generate the Random Error Terms ###
y_err=rnorm(n)
y_err=y_err*sqrt((1-pve)/var(y_err))

### Generate the Phenotypes ###
y=c(y_marginal+y_epi+y_err)
y=(y-mean(y))/(sd(y))

### Create the training and test sets ###
train.idc = sample(1:length(y), size=0.8*length(y), replace=FALSE)

X_train = X[train.idc,]; y_train = as.numeric(y[train.idc])
X_test = X[-train.idc,]; y_test = as.numeric(y[-train.idc])

######################################################################################
######################################################################################
######################################################################################

### Running the Bayesian Approximate Kernel Regression Modeling Framework (BAKR) ###

### Create the Approximate Kernel Matrix ###
#This function takes on two arguments:
#(1) The Genotype matrix X. This matrix should be fed in as a pxn matrix. That is, SNPs are the rows and subjects/patients/cell lines are the columns.
#(2) The number of mcmc draws one wants to use in order to make K_tilde \approx K

n = dim(X_train)[1] #Sample size of the training set
p = dim(X_train)[2] #Number of markers or genes
K_tilde = ApproxGaussKernel(t(X_train),p,p)

### Center the Approximate Kernel Matrix for numerical stability ###
v=matrix(1, n, 1)
M=diag(n)-v%*%t(v)/n
Kn=M%*%K_tilde%*%M
Kn=Kn/mean(diag(Kn))

#Use the Eigenvalue Decomposition of the Kernel matrix in order to futher reduce dimensionality as described in the paper
#Choose the number of components that explain desired % of the cumulative variance in the data (i.e. q in the paper)
evd = EigDecomp(Kn)

explained_var = cumsum(evd$lambda/sum(evd$lambda))
q = 1:min(which(explained_var >= 0.99))
Lambda = diag(sort(evd$lambda,decreasing = TRUE)[q]^(-1)) # Matrix of Eigenvalues
U = evd$U[,q] # Unitary Matrix of Eigenvectors

#Compute the inverse mapping used to obtain the effect sizes of the original covariates 
B = InverseMap(t(X_train),U)

#Set up the desired number of MCMC samples and burn-ins
mcmc.iter = 2e3
mcmc.burn = 1e3

#Run the BAKR Gibbs Sampler. This function takes arguments: approximate kernel factor matrix (U), response (y); truncated covariance matrix (Lambda); and the number of iterations and burn-in for the MCMC. When y is binary, use the BAKRProbitGibbs() instead---it takes the same arguments. 
Gibbs = BAKRGibbs(U,y_train,Lambda,mcmc.iter,mcmc.burn)

#The function results in a list with MCMC samples for the kernel factor coefficients (Gibbs$theta), and the variance scale terms for theta and the random error terms (Gibbs$sigma_theta and Gibbs$sigma_e). Gibbs$sigma_e is not available when the BAKRProbitGibbs() since e ~ MVN(0,I)

#NOTE: Gibbs$theta is a (iter-burn.in) x q matrix; Gibbs$sigma_theta and Gibbs$sigma_e are (iter-burn.in) long vectors

######################################################################################
######################################################################################
######################################################################################

### Prediction in the BAKR framework ###

#Get the posterior mean of the original betas using the inverse mapping B and the posterior mean of theta (*)
beta.out = PostBeta(B,as.numeric(PostMean(Gibbs$theta)))

#Carry out prediction using the posterior mean of the beta. Again X is input as p x n matrix. When the outcome is binary, use the BAKRProbitPredict() insted---it takes the same arguments and produces classification probabilities.
BAKR_pred = BAKRPredict(t(X_test),beta.out)

#Check the predictive performance---in this case, mean square prediction error (MSPE)
MSPE_BAKR = mean((y_test-BAKR_pred)^2); MSPE_BAKR
