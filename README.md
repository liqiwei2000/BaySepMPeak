# BaySepMPeak
This repository is used for accessing the performance of Bayesian zero-inflated negative binomial linear regression model with spatial variable selection proposed in the submitted manuscript titled "..." to analyze MeRIP-Seq data (in multiple conditions).

When running "model_fitting.R" to fit the proposed model, please first load MeRIP-Seq data data. You can load your own data or one of the two examples in the repository. The necessary inputs should be a n-by-W count matrix denoted by Y, where n is the number of samples and W is the number of bins, and a design matrix to indicate the membership of each sample.

*Note 1, the notations in the code follow the notations in the manuscript.
