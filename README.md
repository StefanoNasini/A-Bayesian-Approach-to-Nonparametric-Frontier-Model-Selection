# A-Bayesian-Approach-to-Nonparametric-Frontier-Model-Selection

This repository contains the replication files associated to the paper “A Bayesian Approach to Nonparametric Frontier Model Selection: Application to U.S. Banks Cost Functions”, which proposes a novel Bayesian approach to nonparametric frontier model selection based on the probabilistic framework of Laplace-type estimators. 

To replicate the results related to the U.S. banks cost functions, one should first run the following:

Bayesian_DEA_MCMC_run_data_CRS_parrallel.R
Bayesian_DEA_MCMC_run_data_VRS_parrallel.R
Bayesian_DEA_MCMC_run_data_NDRS_parrallel.R
Bayesian_DEA_MCMC_run_data_NIRS_parrallel.R

These scripts needs to be run for each of years 2001 to 2010 by changing

myYear <- 2001

on the line 40.
