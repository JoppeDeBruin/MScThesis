# This script provides insight in the experiments for Chapter 5. It compares DRP and WCP for increasing N

rm(list = ls())

library(randomForest)
library(gbm)
library(parallel)
library(pbmcapply)
library(purrr)
library(dplyr)

source("utils.R")

set.seed(1234)

alpha=0.1
p=4 # dimensionality of the regression problem

w = function(x) {
  #expit(x%*% c(-1,0.5,-0.25,-0.1))
  exp(x%*% c(-1,0.5,-0.25,-0.1))
}

# Create a covariate shifted sample using repeated weighted sampling
wsample = function(wts, frac=1) {
  n = length(wts)
  i = c()
  while(length(i) < n*frac) {
    j = which( runif(n) <= wts/sum(wts)  )[1] 
    if ( !is.na(j) ) {
      i = c(i, j )
    }
  }
  return(i)
}

beta = c(27.4,13.7,13.7,13.7)

# Define the variance structure
het = TRUE
sd.fun <- function(x, het){
  if (het){return(exp(0.3*x)+1)}
  return(1)
}

# Run the simulation for a specific setting of weights. It performs DCP, CQR and DRP with quantile regression forests
run_example<- function(R, N, n, weights = "RF", parallel = FALSE){
  
  cov.drp <- cov.cqr <- cov.dcp<- rep(NA, R)
  len.drp <- len.cqr <- len.dcp <- rep(NA, R)
  len.drp.med <- len.cqr.med <- len.dcp.med <- rep(NA, R)
  percentage.inf.cqr <- rep(NA, R)
  
  # One iteration of the simulation
  one_iteration <- function(r){
    
    if (r %% 2 == 0) cat(r,".. ")
    if (seeds) {set.seed(r)}
    
    # Generate the data
    dat.x=matrix(rnorm(N*p),N,p)
    dat.y = 210+dat.x%*%beta+sd.fun(dat.x[,1], het)*rnorm(N)
    
    i = sample(N,n)
    x = dat.x[i,]; y = dat.y[i]
    x0 = dat.x[-i,]; y0 = dat.y[-i] # This is going to be the test data (will still be tilted)

    
    # Tilting
    i0 = wsample(w(x0)) # We use a weighted sample of the validation data...
    
    n0 = floor(length(i0))
    x00 = x0[i0,]; y00 = y0[i0] # Now x00 is the weighted sample of the validation data 
    
    # Run the DRP method on the shifted data (weights are estimated internal in the DRP function)
    res1 <- DRP(rbind(x,x00), c(y, y00), c(rep(0, nrow(x)), rep(1, nrow(x00))), x00, method = "RF", score = "CQR", 
                method_pi = weights, model_m = "boosting")
    cov.drp <- mean(res1$intervals[,1] <= y00 & y00 <= res1$intervals[,2])
    len.drp <- mean(res1$intervals[,2] - res1$intervals[,1])
    len.drp.med <- median(res1$intervals[,2] - res1$intervals[,1])
    
    # Compute weights for the WCP method for different options
    if (weights == "RF"){
      model.weights <- randomForest(rbind(x, x00), as.factor(c(rep(0, nrow(x)), rep(1, nrow(x00)))), ntree = 500, nodesize = 10)
      preds.w.rf <- predict(model.weights, rbind(x, x00), type = "prob")[,2]
      preds.w.rf <- pmax(pmin(preds.w.rf, 0.99), 0.01)
      w <- preds.w.rf/(1-preds.w.rf)
    } else if (weights == "boosting"){
      model.weights.gbm <- gbm::gbm(Y ~ ., distribution = "bernoulli",
                                    data = data.frame(Y = c(rep(0, nrow(x)), rep(1, nrow(x00))), rbind(x, x00)), n.trees = 200)
      preds.w.gbm <- predict(model.weights.gbm, data.frame(rbind(x, x00)), type = "response", n.trees = 200)
      w <- preds.w.gbm/(1-preds.w.gbm)
    } else if (weights == "glm"){
      model.weights.gbm <- glm(Y ~ ., family = "binomial",
                               data = data.frame(Y = c(rep(0, nrow(x)), rep(1, nrow(x00))), rbind(x, x00)))
      preds.w.gbm <- predict(model.weights.gbm, data.frame(rbind(x, x00)), type = "response")
      w <- preds.w.gbm/(1-preds.w.gbm)
    }
    
    #Perform CQR 
    res2 <- conformal.split(x, y, x00, w = w, method = "CQR", y_model = "RF")
    cov.cqr <-mean(res2$intervals[,1] <= y00 & y00 <= res2$intervals[,2])
    len.cqr <- mean(res2$intervals[,2] - res2$intervals[,1])
    len.cqr.med <- median(res2$intervals[,2] - res2$intervals[,1])
    
    # Perform DCP
    res3 <- conformal.split(x, y, x00, w = w, method = "DCP", y_model = "RF")
    cov.dcp <-mean(res3$intervals[,1] <= y00 & y00 <= res3$intervals[,2])
    len.dcp <- mean(res3$intervals[,2] - res3$intervals[,1])
    len.dcp.med <- median(res3$intervals[,2] - res3$intervals[,1])
    
    # COmpute how many intervals have infiinite length for CQR
    percentage.inf.cqr <- mean(is.infinite(res2$intervals[,2] - res2$intervals[,1]))
    
    
    return(list(cov.drp = cov.drp, cov.cqr=cov.cqr, cov.dcp = cov.dcp, 
                len.drp = len.drp, len.cqr = len.cqr, len.dcp = len.dcp, 
                len.drp.med = len.drp.med, len.cqr.med = len.cqr.med, len.dcp.med = len.dcp.med, 
                percentage.inf.cqr = percentage.inf.cqr))
  }
  
  # Note: parallel can only be used on linux machines! 
  if (parallel){
    res <- pbmclapply(1:R, one_iteration, mc.cores = detectCores()[1]-2)  
  } else{
    res <- mclapply(1:R, one_iteration, mc.cores = 1)
  }
  
  # Construct the results matrices from the list return
  cov.drp <- res %>% map(1) %>% do.call(rbind, .)
  cov.cqr <- res %>% map(2) %>% do.call(rbind, .)
  cov.dcp <- res %>% map(3) %>% do.call(rbind, .)
  
  len.drp <- res %>% map(4) %>% do.call(rbind, .)
  len.cqr <- res %>% map(5) %>% do.call(rbind, .)
  len.dcp <- res %>% map(6) %>% do.call(rbind, .)
  
  len.drp.med <- res %>% map(7) %>% do.call(rbind, .)
  len.cqr.med <- res %>% map(8) %>% do.call(rbind, .)
  len.dcp.med <- res %>% map(9) %>% do.call(rbind, .)
  
  percentage.inf.cqr <- res %>% map(10) %>% do.call(rbind, .)
  
  return(list(cov.drp = cov.drp, cov.cqr=cov.cqr, cov.dcp = cov.dcp, 
              len.drp = len.drp, len.cqr = len.cqr, len.dcp = len.dcp, 
              len.drp.med = len.drp.med, len.cqr.med = len.cqr.med, len.dcp.med = len.dcp.med, 
              percentage.inf.cqr = percentage.inf.cqr))
}

set.seed(1234)
seeds = TRUE
par = TRUE

# Run the simulations for different settings (weights and sample sizes)
results_1500_400.rf <- run_example(1000, 2400, 400, weights = "RF", parallel = par)
results_1500_750.rf <- run_example(1000, 2750, 750, weights = "RF", parallel = par)
results_3000_1000.rf <- run_example(1000, 3000, 1000, weights = "RF", parallel = par)
results_5000_3000.rf <- run_example(1000, 5000, 3000, weights = "RF", parallel = par)
results_7000_5000.rf <- run_example(1000, 7000, 5000, weights = "RF", parallel = par)
results_9000_7000.rf <- run_example(1000, 9000, 7000, weights = "RF", parallel = par)
#save(results_1500_750.rf, results_3000_1000.rf,results_5000_3000.rf,results_7000_5000.rf,results_9000_7000.rf, file = paste0("results/DRPvsCQR_synthetic_increasing_n_QRF_RF_weights", Sys.Date(),"R=1000", ".Rdata"))

set.seed(1234)
results_1500_400.boost <- run_example(1000, 2400, 400, weights = "boosting", parallel = par)
results_1500_750.boost <- run_example(1000, 2750, 750, weights = "boosting", parallel = par)
results_3000_1000.boost <- run_example(1000, 3000, 1000, weights = "boosting", parallel = par)
results_5000_3000.boost <- run_example(1000, 5000, 3000, weights = "boosting", parallel = par)
results_7000_5000.boost <- run_example(1000, 7000, 5000, weights = "boosting", parallel = par)
results_9000_7000.boost <- run_example(1000, 9000, 7000, weights = "boosting", parallel = par)
#save(results_1500_750.boost, results_3000_1000.boost,results_5000_3000.boost,results_7000_5000.boost,results_9000_7000.boost, file = paste0("results/DRPvsCQR_synthetic_increasing_n_QRF_boosted_weights","R=1000", Sys.Date(), ".Rdata"))

set.seed(1234)
results_1500_400.glm <- run_example(1000, 2400, 400, weights = "glm", parallel = par)
results_1500_750.glm <- run_example(1000, 2750, 750, weights = "glm", parallel = par)
results_3000_1000.glm <- run_example(1000, 3000, 1000, weights = "glm", parallel = par)
results_5000_3000.glm <- run_example(1000, 5000, 3000, weights = "glm", parallel = par)
results_7000_5000.glm <- run_example(1000, 7000, 5000, weights = "glm", parallel = par)
results_9000_7000.glm <- run_example(1000, 9000, 7000, weights = "glm", parallel = par)
#save(results_1500_750.glm, results_3000_1000.glm,results_5000_3000.glm,results_7000_5000.glm,results_9000_7000.glm, file = paste0("results/DRPvsCQR_synthetic_increasing_n_QRF_glm_weights", Sys.Date(),"R=1000", ".Rdata"))


## ANALYSE THE RESULTS/Visualize

# Plots below are test plots to play around with
results_1500_750 <- results_1500_750.rf
results_3000_1000 <- results_3000_1000.boost
results_5000_3000 <- results_5000_3000.boost
results_7000_5000 <- results_7000_5000.boost
results_9000_7000 <- results_9000_7000.boost


mean(results_1500_750$cov.drp)
mean(results_1500_750$cov.cqr)

#pdf("figs/DRPvsCQR_synthetic_increasing_n.pdf", width = 12, height = 6)
par(mfrow = c(1,3))
names <- c("750", "1000", "3000", "5000", "7000")
boxplot(cbind(results_1500_750$cov.drp, results_3000_1000$cov.drp, results_5000_3000$cov.drp, results_7000_5000$cov.drp, results_9000_7000$cov.drp), 
        ylim = c(0.65, 1), names = names, main = "DRP", ylab = "Coverage")
abline(h = 0.9, lty = 2)
boxplot(cbind(results_1500_750$cov.cqr, results_3000_1000$cov.cqr, results_5000_3000$cov.cqr, results_7000_5000$cov.cqr, results_9000_7000$cov.cqr), 
        ylim = c(0.65, 1), names = names, main = "CQR", ylab = "Coverage")
abline(h = 0.9, lty = 2)
boxplot(cbind(results_1500_750$cov.dcp, results_3000_1000$cov.dcp, results_5000_3000$cov.dcp, results_7000_5000$cov.dcp, results_9000_7000$cov.dcp), 
        ylim = c(0.65, 1), names = names, main = "DCP", ylab = "Coverage")
abline(h = 0.9, lty = 2)


#graphics.off()

#pdf("figs/DRPvsCQR_synthetic_increasing_n_len.pdf", width = 12, height = 6)
#par(mfrow = c(1,2))
names <- c("750", "1000", "3000", "5000", "7000")
boxplot(cbind(results_1500_750$len.drp, results_3000_1000$len.drp.med, results_5000_3000$len.drp.med, results_7000_5000$len.drp.med, results_9000_7000$len.drp.med), 
        names = names, main = "DRP", ylim =c(30, 150), ylab = "Length")
abline(h = 0.9, lty = 2)
boxplot(cbind(results_1500_750$len.cqr, results_3000_1000$len.cqr.med, results_5000_3000$len.cqr.med,results_7000_5000$len.cqr.med, results_9000_7000$len.cqr.med), 
        names = names, main = "CQR", ylim = c(30, 150), ylab = "Length")
abline(h = 0.9, lty = 2)
boxplot(cbind(results_1500_750$len.dcp, results_3000_1000$len.dcp.med, results_5000_3000$len.dcp.med,results_7000_5000$len.dcp.med, results_9000_7000$len.dcp.med), 
        names = names, main = "DCP", ylim = c(30, 150), ylab = "Length")
abline(h = 0.9, lty = 2)
#graphics.off()


# These plots create the results from the paper
melt(results_1500_750.glm,)
res_1500_750.glm <- as.data.frame(results_1500_750.glm)

data.plot <- data.frame(rbind(cbind(melt(as.data.table(results_1500_400.glm), measure.vars = 1:10), rep("glm", 2500), rep("400", 2500)),
                              cbind(melt(as.data.table(results_1500_750.glm), measure.vars = 1:10), rep("glm", 2500), rep("750", 2500)),
                              cbind(melt(as.data.table(results_3000_1000.glm), measure.vars = 1:10), rep("glm", 2500), rep("1000", 2500)),
                              cbind(melt(as.data.table(results_5000_3000.glm), measure.vars = 1:10), rep("glm", 2500), rep("3000", 2500)),
                              cbind(melt(as.data.table(results_7000_5000.glm), measure.vars = 1:10), rep("glm", 2500), rep("5000", 2500)),
                              cbind(melt(as.data.table(results_9000_7000.glm), measure.vars = 1:10), rep("glm", 2500), rep("7000", 2500)),
                              cbind(melt(as.data.table(results_1500_400.boost), measure.vars = 1:10), rep("boost", 2500), rep("400", 2500)),
                              cbind(melt(as.data.table(results_1500_750.boost), measure.vars = 1:10), rep("boost", 2500), rep("750", 2500)),
                              cbind(melt(as.data.table(results_3000_1000.boost), measure.vars = 1:10), rep("boost", 2500), rep("1000", 2500)),
                              cbind(melt(as.data.table(results_5000_3000.boost), measure.vars = 1:10), rep("boost", 2500), rep("3000", 2500)),
                              cbind(melt(as.data.table(results_7000_5000.boost), measure.vars = 1:10), rep("boost", 2500), rep("5000", 2500)),
                              cbind(melt(as.data.table(results_9000_7000.boost), measure.vars = 1:10), rep("boost", 2500), rep("7000", 2500)),
                              cbind(melt(as.data.table(results_1500_400.rf), measure.vars = 1:10), rep("rf", 2500), rep("400", 2500)),
                              cbind(melt(as.data.table(results_1500_750.rf), measure.vars = 1:10), rep("rf", 2500), rep("750", 2500)),
                              cbind(melt(as.data.table(results_3000_1000.rf), measure.vars = 1:10), rep("rf", 2500), rep("1000", 2500)),
                              cbind(melt(as.data.table(results_5000_3000.rf), measure.vars = 1:10), rep("rf", 2500), rep("3000", 2500)),
                              cbind(melt(as.data.table(results_7000_5000.rf), measure.vars = 1:10), rep("rf", 2500), rep("5000", 2500)),
                              cbind(melt(as.data.table(results_9000_7000.rf), measure.vars = 1:10), rep("rf", 2500), rep("7000", 2500))))

data.plot.cov <- data.plot %>% filter(grepl("cov", variable)) %>%
  mutate(variable = recode_factor(variable, "cov.cqr" = "CQR", "cov.dcp" = "DCP", "cov.drp" = "DRP" )) %>%
  mutate(V2 = recode_factor(V2, "rf" = "RF","boost" = "Boosting", "glm" = "glm"))


plot1 <- ggplot(data = data.plot.cov, aes(x = V3, y = value)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.9, col = "blue") +
  facet_grid(V2 ~ variable) +
  scale_x_discrete(limits = (c("400", "750", "1000", "3000", "5000", "7000")))+
  #coord_flip() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  ylab("Empirical Coverage") + xlab("N") +
  coord_cartesian(ylim = c(0.77, 1))
plot1
#ggsave("figs/WCP_vs_DRP_increasing_n_coverage.pdf", width = 9, height = 9, plot1)

data.plot.len <- data.plot %>% filter(grepl("len", variable) & !grepl("med", variable)) %>%
  mutate(variable = recode_factor(variable, "len.cqr" = "CQR", "len.dcp" = "DCP", "len.drp" = "DRP" )) %>%
  mutate(V2 = recode_factor(V2, "rf" = "RF","boost" = "Boosting", "glm" = "glm"))


plot2 <- ggplot(data = data.plot.len, aes(x = V3, y = value)) +
  geom_boxplot() +
  #geom_hline(yintercept = 0.9, col = "blue") +
  facet_grid(V2 ~ variable) +
  scale_x_discrete(limits = (c("400", "750", "1000", "3000", "5000", "7000")))+
  #coord_flip() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  ylab("Average Length") + xlab("N") #+
#coord_cartesian(ylim = c(0.77, 1))
plot2

#ggsave("figs/WCP_vs_DRP_increasing_n_length.pdf", width = 9, height = 9, plot2)
