#NOTE: This script will walk through the BAKR functions. Specifically it shows:
#(1) How to compute the approximate Kernel matrix and its singular value decomposition form
#(2) Run the Gibbs Sampler for BAKR
#(3) Retrieve the beta estimates for the original variants/genes/obeserved variables
#(4) Conduct inference and/or out-of-sample prediction 

#NOTE: This script is based on a simple (and small) genetics example where we simulate genotype data for n = 500 individuals with p = 2000 measured SNPs. We will randomly select a small number (e.g. 25) of these SNPs to be causal and have true association with the generated (continuous) phenotype y

#NOTE: It is also noted how functions change when the variables are binary instead of continuous

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
nidv = 500; nvar = 2e3; ncausal = 25; pve = 0.4; rho = 1 # number of data points

### Generate the data ###
X = matrix(runif(nidv*nvar),nrow = nidv,ncol = nvar)
colnames(X) = paste("SNP",1:ncol(X),sep = "") #Give the SNPs names

### Randomly select causal variables ###
s=sample(dim(X)[2], ncausal, replace=F)
Xcausal=X[,s]

### Generate the genetic effects ###
beta1=rnorm(ncausal)
y_add=Xcausal%*%beta1

y_err=rnorm(nidv) #Random Error/Noise

y=y_add+y_err #Generate Phenotypes (i.e. Responses)
y=scale(y) #Standardize the Phenotypes

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

n = dim(X)[1] #Sample size
p = dim(X)[2] #Number of markers or genes
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

######################################################################################
######################################################################################
######################################################################################

### Inference in the BAKR framework ###

#Get all of the posterior draws for beta estimates. [We can then use the PostMean function to get the equivalent of beta.out in (*)] 
beta_draws = GetBeta(B,Gibbs$theta)

### Compute the Posterior Probability of Association Analog (PPAA) ###
sigval= 0.05 #Set the desired FDR
PPAA.out = PostMean(GetPPAAs(GetBeta(B,Gibbs$theta),sigval))
names(PPAA.out) = colnames(X)

#Check the PPAA of the causal variants 
summary(PPAA.out); PPAA.out[s]

#Compute the local false sign rates [Stephens, M. (2016). False discovery rates: A new deal. bioRxiv]. The local false sign rate is analogous to the local false discovery rate and provides a measure of confidence in the sign of an effect rather than confidence of the effect being non-zero. The lower the lfsr, the better. This function takes the beta draws from MCMC iteration. The desired significance threshold for these values maybe chosen subjectively
LFSR = as.numeric(lfsr(beta_draws)); names(LFSR) = colnames(X)

#Check the lfsr of the causal variants 
summary(LFSR); LFSR[s]