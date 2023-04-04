# MScThesis

This respository contains the code to reproduce the results and figures in my MSc thesis "Topics in Conformal Prediction"

Except for installing the necessary dependencies, all scripts should be runnable by simply including the utils.R file in the working directory. Note that some simulations have the option to parallelize over multiple cores. This feature only works on Linux machines and can cause trouble on Windows. 

The utils file contains all helper function that define all base learners, conformal prediction methods, and doubly robust prediction set approach. It also contains a method to construct prediction intervals for the rbart model (which only provides point predictions).

The other files can be run to reconstruct the results in the paper:

- visualize_conformal.R produces the plots in Chapter 2
- cov_hists_1D.R produces the histograms of empirical coerage of the simple 1D syntehtic examples in Chapter 4
- simulation_airfoil_results.R performs the simulations for the airfoil example in Chapter 5. Note that this simulation takes up much computational resources and can be run in parallel only on Linux.
- simulation_function_airfoil.R this script defines the functions needed to perform the simulations for the airfoil dataset in Chapter 5. It defines and sets up all the experiments
- boxplots_airfoil.R this script plots all the results for the airfoil simulations
- plot_weighted_quantiles.R plots the dummy example comparing the normal empirical quantile vs the weighted quantile in Chapter 5
- synthetic_increasing_N.R performs the synthetic simulations in Chapter 5 for growing sample sizes. It is defined for DCP, CQR and DRP with quantile regression forests as base learner. Different weight regimes are considered
- synthetic_1D_simulations.R constructs the results in Chapter 4 on the 1D syntethic example to visualize one run for each conformal predictor. It also considers the adaptivity and conditional coverage experiments.
- conditional_coverage_CQR_BIO.R performs CQR on the bio dataset with different nominal quantile levels for the base quantile regressor. And it plots the conditional coverage for these settings.
