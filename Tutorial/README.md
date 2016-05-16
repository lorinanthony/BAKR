# Tutorial on Running BAKR
In the file BAKR_Tutorial.R, we walk through the BAKR functions. Specifically we show:
(1) How to compute the approximate kernel matrix and its singular value decomposition form;
(2) Run the Gibbs Sampler for BAKR;
(3) Retrieve the beta estimates for the original variants/genes/obeserved variables;
(4) Conduct inference and/or out-of-sample prediction.

This script is based on a simple (and small) genetics example where we simulate genotype data for n = 500 individuals with p = 2000 measured SNPs. Briefly, we randomly select a small number (e.g. 25) of these SNPs to be causal and have true association with the generated (continuous) phenotype y.

In this script, it is also noted how functions change when the response variables are binary (i.e. for classification problems) instead of continuous. 
