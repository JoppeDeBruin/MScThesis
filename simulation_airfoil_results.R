#### This script performs the Simulations on the airfoil data in Chapter 5
rm(list = ls())

library(randomForest)
library(drf)
library(ggplot2)
library(grf)
library(quantreg)
library(rbart)
library(grf)
library(parallel)
library(pbmcapply)
library(purrr)
library(dplyr)

#etwd("~/ETH Zurich/MSc Thesis/thesis")
source("simulation_function_airfoil.R")
source("utils.R")

#### Load and prepare data
dat = read.table("airfoil.txt")
dim(dat)
colnames(dat) = c("Frequency",
                  "Angle",
                  "Chord",
                  "Velocity",
                  "Suction",
                  "Sound")

dat.x = as.matrix(dat[,1:5])
dat.y = as.numeric(dat[,6])
dat.x[,1] = log(dat.x[,1]) # Log transform
dat.x[,5] = log(dat.x[,5]) # Log transform
N = nrow(dat.x); p = ncol(dat.x)

alpha = 0.1

# Run simulations for the airfoil dataset. Note that parallelization is only possible on Linux and not on windows!
start = Sys.time()

res.mcp.rf<- run.airfoil.simulation(conformal.split, dat.x, dat.y, N, n_sim = 1000, all = T, 
                                    seeds = TRUE, method = "mean", y_model = "RF", parallel = FALSE)
res.cqr.rf <- run.airfoil.simulation(conformal.split, dat.x, dat.y, N, n_sim = 1000, all = T, 
                                     seeds = TRUE, method = "CQR", y_model = "RF", parallel = FALSE)
res.dcp.rf <- run.airfoil.simulation(conformal.split, dat.x, dat.y, N, n_sim = 1000, all = T, 
                                     seeds = TRUE, method = "DCP", y_model = "RF", parallel = FALSE)
res.naive.rf <- run.airfoil.simulation(naive.QR, dat.x, dat.y, N, n_sim = 1000, all = F, 
                                       seeds = TRUE, y_model = "RF", parallel = FALSE)

res.mcp.boost <- run.airfoil.simulation(conformal.split, dat.x, dat.y, N, n_sim = 1000, all = T, 
                                        seeds = TRUE, method = "mean", y_model = "boosting", parallel = FALSE)
res.cqr.boost <- run.airfoil.simulation(conformal.split, dat.x, dat.y, N, n_sim = 1000, all = T, 
                                        seeds = TRUE, method = "CQR", y_model = "boosting", parallel = FALSE)
res.dcp.boost <- run.airfoil.simulation(conformal.split, dat.x, dat.y, N, n_sim = 1000, all = T, 
                                        seeds = TRUE, method = "DCP", y_model = "boosting", parallel = FALSE)
res.naive.boost <- run.airfoil.simulation(naive.QR, dat.x, dat.y, N, n_sim = 1000, all = F, 
                                          seeds = TRUE, y_model = "boosting", parallel = FALSE)

# Note that due to parallelization of bart computation, it is better to run these sequentially!
res.mcp.bart <- run.airfoil.simulation(conformal.split, dat.x, dat.y, N, n_sim = 1000, all = T, 
                                       seeds = TRUE, method = "mean", y_model = "bart", parallel = FALSE)
res.cqr.bart <- run.airfoil.simulation(conformal.split, dat.x, dat.y, N, n_sim = 1000, all = T, 
                                       seeds = TRUE, method = "CQR", y_model = "bart", parallel = FALSE)


res.cb <- run.airfoil.simulation(conformal.split, dat.x, dat.y, N, n_sim = 1000, all = T, 
                                 seeds = TRUE, method = "CB", y_model = "bart_het", parallel = FALSE)
end = Sys.time()

print(end-start)
#save(res.mcp.rf, res.cqr.rf, res.dcp.rf, res.naive.rf, res.mcp.boost, res.cqr.boost, res.dcp.boost, res.naive.boost, res.mcp.bart,
#     res.cqr.bart, res.cb, file = paste0("results/res_airfoil_", Sys.Date(), "_1000.Rdata"))



start = Sys.time()

# Perform the simluations for DRP with different settings
# DO NOT PROVIDE THE WEIGHTS HERE! They are defined in the simulation function!

res.drp.mean.rf_r.rf_m <- run.airfoil.simulation.DRP(dat.x, dat.y, N, n_sim = 1000, seed = TRUE, score = "mean",
                                                     method = "RF", method_m = "binary", model_m = "RF", parallel = FALSE)
res.drp.mean.rf_r.boost_m <- run.airfoil.simulation.DRP(dat.x, dat.y, N, n_sim = 1000, seed = TRUE, score = "mean",
                                                        method = "RF", method_m = "binary", model_m = "boosting", parallel = FALSE)
res.drp.cqr.rf_r.rf_m <- run.airfoil.simulation.DRP(dat.x, dat.y, N, n_sim = 1000, seed = TRUE, score = "CQR",
                                                    method = "RF", method_m = "binary", model_m = "RF", parallel = FALSE)
res.drp.cqr.rf_r.boost_m <- run.airfoil.simulation.DRP(dat.x, dat.y, N, n_sim = 1000, seed = TRUE, score = "CQR",
                                                       method = "RF", method_m = "binary", model_m = "boosting", parallel = FALSE)

# Then with boosting as base learner
res.drp.mean.boost_r.rf_m <- run.airfoil.simulation.DRP(dat.x, dat.y, N, n_sim = 1000, seed = TRUE, score = "mean",
                                                        method = "boosting", method_m = "binary", model_m = "RF", parallel = FALSE)
res.drp.mean.boost_r.boost_m <- run.airfoil.simulation.DRP(dat.x, dat.y, N, n_sim = 1000, seed = TRUE, score = "mean",
                                                           method = "boosting", method_m = "binary", model_m = "boosting", parallel = FALSE)
res.drp.cqr.boost_r.rf_m <- run.airfoil.simulation.DRP(dat.x, dat.y, N, n_sim = 1000, seed = TRUE, score = "CQR",
                                                       method = "boosting", method_m = "binary", model_m = "RF", parallel = FALSE)
res.drp.cqr.boost_r.boost_m <- run.airfoil.simulation.DRP(dat.x, dat.y, N, n_sim = 1000, seed = TRUE, score = "CQR",
                                                          method = "boosting", method_m = "binary", model_m = "boosting", parallel = FALSE)

# Then with BART as base learner
start = Sys.time()
res.drp.mean.bart_r.rf_m <- run.airfoil.simulation.DRP(dat.x, dat.y, N, n_sim = 5, seed = TRUE, score = "mean",
                                                       method = "bart", method_m = "binary", model_m = "RF", parallel = FALSE)
res.drp.mean.bart_r.boost_m <- run.airfoil.simulation.DRP(dat.x, dat.y, N, n_sim = 5, seed = TRUE, score = "mean",
                                                          method = "bart", method_m = "binary", model_m = "boosting")
res.drp.cqr.bart_r.rf_m <- run.airfoil.simulation.DRP(dat.x, dat.y, N, n_sim = 5, seed = TRUE, score = "CQR",
                                                      method = "bart", method_m = "binary", model_m = "RF")
res.drp.cqr.bart_r.boost_m <- run.airfoil.simulation.DRP(dat.x, dat.y, N, n_sim = 5, seed = TRUE, score = "CQR",
                                                         method = "bart", method_m = "binary", model_m = "boosting")
end = Sys.time()
print(end-start)

#save(res.drp.mean.rf_r.rf_m, res.drp.mean.rf_r.boost_m, res.drp.cqr.rf_r.rf_m, res.drp.cqr.rf_r.boost_m,
#     res.drp.mean.boost_r.rf_m, res.drp.mean.boost_r.boost_m, res.drp.cqr.boost_r.rf_m, res.drp.cqr.boost_r.boost_m,
#     file = paste0("results/res_airfoil_DRP_", Sys.Date(), ".Rdata"))