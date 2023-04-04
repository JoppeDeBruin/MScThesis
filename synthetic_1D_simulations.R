## This script simulates some 1D regression tasks with homoscedastic and heteroscedsatic cases and compares 
# the performance of the different conformal methods in terms of coverage, length and conditional coverage
# open question is how many simulations to use to smooth out some of the analyses

# Load requirements and help functions
rm(list = ls())
library(randomForest)
library(grf)
library(bartMachine)
library(rbart)
library(drf)

#setwd("~/ETH Zurich/MSc Thesis/thesis")
source("utils.R")

####

# Check options for noise structure
homoscedastic <- F
t.noise <- F
retrain <- T # When you want to average over some models to check conditional coverage, maybe should retrain
# otherwise it's just the same as taking more calibration data (not exactly, since the conformal scores change)

R <- 5 # Number of repetitions to smooth out results for conditional coverage graphs

n.train <- 2000
n.cal <- 300
n.test <- 2000
alpha <- 0.1
alpha.fit <- 0.1
distr.grid <- 200

set.seed(1)

# Generate training data
train.data <- generate.data.1D(n.train, homoscedastic = homoscedastic, tnoise = t.noise)
x.train <- train.data$x; y.train <- train.data$y
x.train.df <- as.data.frame(x.train)

plot(x.train, y.train)

if (!retrain){
  if (homoscedastic){
    CB.train <- rbart(x.train.df, y.train , ntree = 50, ntreeh = 1, pbd = c(0.7,0.0))
  } else{
    CB.train <- rbart(x.train.df, y.train, ntree = 50, ntreeh = 20)
  }
  
  CQR.train <- quantile_rf(x.train.df, y.train, quantiles = c(alpha.fit/2, 1-alpha.fit/2))$model
  DCP.train <- quantile_rf(x.train.df, y.train, quantiles = seq(0.001, 0.999, length = distr.grid))$model
  MCP.train <- rf(x.train.df, y.train, x.train.df)$model
}

intervals.mat <- NULL
x.test.vec <- c(); y.test.vec <- c()

for (r in 1:R){
  
  cat(r, "..")
  
  # Generate new training dataset if desired and retrain models
  if (retrain){
    train.data <- generate.data.1D(n.train, homoscedastic = homoscedastic, tnoise = t.noise)
    x.train <- train.data$x; y.train <- train.data$y
    x.train.df <- as.data.frame(x.train)
    
    if (homoscedastic){
      CB.train <- rbart(x.train.df, y.train , ntree = 50, ntreeh = 1, pbd = c(0.7,0.0))
    } else{
      CB.train <- rbart(x.train.df, y.train, ntree = 50, ntreeh = 20)
    }
    
    CQR.train <- quantile_rf(x.train.df, y.train, x.train.df, quantiles = c(0.05, 0.95))$model
    DCP.train <- quantile_rf(x.train.df, y.train, x.train.df, quantiles = seq(0.001, 0.999, length = 200))$model
    MCP.train <- rf(x.train.df, y.train, x.train.df)$model
    
    QR.train <- quantile_rf(x.train.df, y.train, x.train.df, quantiles = seq(0.01, 0.99, length = 200))$model
    
  }
  
  data.calibration <- generate.data.1D(n.cal, homoscedastic = homoscedastic, tnoise = t.noise)
  x.cal.df <- as.data.frame(data.calibration$x); y.cal <- data.calibration$y
  
  data.test <- generate.data.1D(n.test, homoscedastic = homoscedastic, tnoise = t.noise)
  x.test.df <- as.data.frame(data.test$x); y.test <- data.test$y
  
  # Perform the different predictio interval methods
  intervals.mcp <- conformal.split(rbind(x.train.df, x.cal.df), c(y.train, y.cal), x.test.df, split = 1:n.train, pretrained.model = MCP.train, 
                                   method = "mean", y_model = "RF")$intervals
  intervals.cqr <- conformal.split(rbind(x.train.df, x.cal.df), c(y.train, y.cal), x.test.df, split = 1:n.train, pretrained.model = CQR.train, 
                                   method = "CQR", y_model = "RF")$intervals
  intervals.distr <- conformal.split(rbind(x.train.df, x.cal.df), c(y.train, y.cal), x.test.df, split = 1:n.train, pretrained.model = DCP.train, 
                                     method = "DCP", y_model = "RF")$intervals
  intervals.cb <- conformal.split(rbind(x.train.df, x.cal.df), c(y.train, y.cal), x.test.df, split = 1:n.train, pretrained.model = CB.train, 
                                  method = "CB", y_model = "bart_het")$intervals
  
  # Also perform the naive calibration
  intervals.naive <- naive.QR(rbind(x.train.df, x.cal.df), c(y.train, y.cal), x.test.df, split = 1:n.train, pretrained.model = QR.train, 
                              y_model = "RF")#$intervals
  opt.alpha = intervals.naive$alpha.opt*2
  intervals.naive <- intervals.naive$intervals
  
  # And the non calibrated method
  intervals.qr <- predict(CQR.train, x.test.df)$predictions 
  
  intervals.mat <- rbind(intervals.mat, cbind(intervals.mcp, intervals.cqr, intervals.distr, intervals.cb, intervals.naive, intervals.qr))
  x.test.vec <- c(x.test.vec, data.test$x)
  y.test.vec <- c(y.test.vec, data.test$y)
  
}

## Visualize the last iteration of the intervals (note: only works for R>1)

par(mfrow = c(2,2))
x.test <- x.test.df$V1
o <- order(x.test)
first.ind <- 1:((R-1)*n.test)
titles = c("MCP", "CQR", "DCP", "CB")
plot.oracle = TRUE
for (j in 1:4){
  plot(x.test, y.test, xlab = "x", ylab = "y", main = titles[j])#, ylim = c(-5, 20))
  lines(x.test[o], intervals.mat[-first.ind,(1+2*(j-1))][o], col = "blue" )
  lines(x.test[o], intervals.mat[-first.ind,(2+2*(j-1))][o], col = "blue" )
  
  if (plot.oracle){
    lines(x.test[o], data.test$y_up[o], col = "purple")
    lines(x.test[o], data.test$y_lo[o], col = "purple")
  }
  
}

x.test <- x.test.df$V1
first.ind <- 1:((R-1)*n.test)
plot_list = list()
titles = c("MCP (RF)", "CQR (QRF)", "DCP (QRF)", "CB (HBART)")
subsample <- sample(1:length(x.test), 500)
for (j in 1:4){
  print(j)
  data.plot <- data.frame(
    Y = y.test,
    X = x.test,
    y.low = intervals.mat[-first.ind, (1+2*(j-1))],
    y.up = intervals.mat[-first.ind, (2+2*(j-1))],
    y.low.oracle = data.test$y_lo,
    y.up.oracle = data.test$y_up
  )
  plot_list[[j]] <- ggplot(data = data.plot, aes(x = X, y = Y)) + 
    geom_ribbon(aes(ymin = y.low, ymax = y.up), fill = "lightblue") +
    geom_point(data = data.frame(X = data.plot$X[subsample], Y = data.plot$Y[subsample]), shape = 1) + 
    geom_line(aes(y = y.up), col = "blue") + 
    geom_line(aes(y = y.low), col = "blue") +
    geom_line(aes(y = y.low.oracle), col = "red") +
    geom_line(aes(y = y.up.oracle), col = "red") +
    theme_light() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggtitle(titles[j]) + theme(plot.title = element_text(hjust = 0.5))
  
}
plot <- gridExtra::grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]])

# Compute the coverages. 
mean(intervals.mcp[,1] <= y.test & y.test <= intervals.mcp[,2])
mean(intervals.cqr[,1] <= y.test & y.test <= intervals.cqr[,2])
mean(intervals.distr[,1] <= y.test & y.test <= intervals.distr[,2])
mean(intervals.cb[,1] <= y.test & y.test <= intervals.cb[,2])

mean(intervals.naive[,1] <= y.test & y.test <= intervals.naive[,2])
mean(intervals.qr[,1] <= y.test & y.test <= intervals.qr[,2])

a = n.cal + 1 - floor((n.cal+1) * alpha)
b = floor((n.cal+1) * alpha)
a / (a+b)
# Change title potentially
#ggsave("figs/1D_conformal_example_heteroscedastic.pdf", width = 12, height = 8 , plot)
#ggsave("figs/1D_conformal_example_homoscedastic.pdf", width = 12, height = 12, plot)



## Now also visualize what would have happened when we use BART as base predictor for CQR (non adaptive!)

CQR.train <- quantile_bart(x.train.df, y.train, x.train.df, quantiles = c(0.05, 0.95))$model
intervals.cqr.bart <- conformal.split(rbind(x.train.df, x.cal.df), c(y.train, y.cal), x.test.df, split = 1:n.train, pretrained.model = CQR.train, 
                                      method = "CQR", y_model = "bart")$intervals

# Also using hbart! Since those will be adaptive, but not after conformalizing anymore!
intervals.hbart <- calc_hbart_intervals(CB.train, x.test.df, quantiles = c(0.05, 0.95))

intervals <- cbind(intervals.cqr.bart, intervals.hbart$quantiles)

titles = c("CQR (BART)", "HBART (non conformal)")
for (j in 1:2){
  print(j)
  data.plot <- data.frame(
    Y = y.test,
    X = x.test,
    y.low = intervals[,(1+2*(j-1))],
    y.up = intervals[, (2+2*(j-1))],
    y.low.oracle = data.test$y_lo,
    y.up.oracle = data.test$y_up
  )
  plot_list[[j]] <- ggplot(data = data.plot, aes(x = X, y = Y)) + 
    geom_ribbon(aes(ymin = y.low, ymax = y.up), fill = "lightblue") +
    geom_point(data = data.frame(X = data.plot$X[subsample], Y = data.plot$Y[subsample]), shape = 1) + 
    geom_line(aes(y = y.up), col = "blue") + 
    geom_line(aes(y = y.low), col = "blue") +
    geom_line(aes(y = y.low.oracle), col = "red") +
    geom_line(aes(y = y.up.oracle), col = "red") +
    theme_light() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggtitle(titles[j]) + theme(plot.title = element_text(hjust = 0.5))
  
}

plot <- gridExtra::grid.arrange(plot_list[[1]], plot_list[[2]], ncol = 2)
#ggsave("figs/1D_conformal_example_heteroscedastic_2.pdf", width = 12, height = 4 , plot)

### Visualize the results for conditional coverage
coverage.mat <- cbind(intervals.mat[,1] <= y.test.vec & y.test.vec <= intervals.mat[,2],
                      intervals.mat[,3] <= y.test.vec & y.test.vec <= intervals.mat[,4],
                      intervals.mat[,5] <= y.test.vec & y.test.vec <= intervals.mat[,6],
                      intervals.mat[,7] <= y.test.vec & y.test.vec <= intervals.mat[,8])
length.mat <- cbind(intervals.mat[,2] - intervals.mat[,1],
                    intervals.mat[,4] - intervals.mat[,3],
                    intervals.mat[,6] - intervals.mat[,5],
                    intervals.mat[,8] - intervals.mat[,7])


# Choose to plot conditional coverage based on X or Y
cond.cov.mat <- binning(x.test.vec, coverage.mat, 20)
cond.len.mat <- binning(x.test.vec, length.mat, 20)


#pdf(file = "figs/coverage_length_1D_homoscedastic.pdf", width = 12, height = 6)
#pdf(file = "figs/coverage_length_1D_heteroscedastic.pdf", width = 12, height = 6)
#graphics.off()

par(mfrow = c(1,2))
plot(cond.cov.mat$xs, cond.cov.mat$cond[,1], type = "b", #ylim = range(cond.cov.mat$cond),
     xlab = "x", ylab = "Coverage", ylim = c(0.5, 1))
lines(cond.cov.mat$xs, cond.cov.mat$cond[,2], type = "b", col = "red")
lines(cond.cov.mat$xs, cond.cov.mat$cond[,3], type = "b", col = "green")
lines(cond.cov.mat$xs, cond.cov.mat$cond[,4], type = "b", col = "blue")
abline(h = 0.9)
legend("bottomleft", legend = c("MCP", "CQR", "DCP", "CB"), col = c("black", "red", "green", "blue"), 
       lty = 1, cex = 0.5)

plot(cond.len.mat$xs, cond.len.mat$cond[,1], type = "b", ylim = range(cond.len.mat$cond),
     xlab = "x", ylab = "Length")
lines(cond.len.mat$xs, cond.len.mat$cond[,2], type = "b", col = "red")
lines(cond.len.mat$xs, cond.len.mat$cond[,3], type = "b", col = "green")
lines(cond.len.mat$xs, cond.len.mat$cond[,4], type = "b", col = "blue")

legend("bottomright", legend = c("MCP", "CQR", "DCP", "CB"), col = c("black", "red", "green", "blue"), 
       lty = 1, cex = 0.5)
#graphics.off()



# Also include what happens when you change the fitted quantiles in CQR! So we fit CQR at different values of alpha

plot.conditional = FALSE
alpha.fit <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
R = 5
par(mfrow = c(4,2))
coverage.mat <- length.mat<- NULL
x.test.mat <- y.test.mat <- intervals.mat.all <- y.low.oracle.mat <- y.up.oracle.mat <- NULL

for (i in 1:length(alpha.fit)){
  intervals.mat <- NULL
  x.test.vec <- c(); y.test.vec <- c()
  y.low.oracle.vec <- c(); y.up.oracle.vec <- c()
  
  # Note that a seed is set here
  set.seed(2023)
  
  for (r in 1:R){
    
    cat(r, "..")
    
    # Generate new training dataset if desired and retrain models
    if (retrain){
      train.data <- generate.data.1D(n.train, homoscedastic = homoscedastic, tnoise = t.noise)
      x.train <- train.data$x; y.train <- train.data$y
      x.train.df <- as.data.frame(x.train)
      
      CQR.train <- quantile_forest(x.train.df, y.train, num.trees = 500, quantiles = c(alpha.fit[i]/2, 1-alpha.fit[i]/2))
    }
    
    data.calibration <- generate.data.1D(n.cal, homoscedastic = homoscedastic, tnoise = t.noise)
    x.cal.df <- as.data.frame(data.calibration$x); y.cal <- data.calibration$y
    
    data.test <- generate.data.1D(n.test, homoscedastic = homoscedastic, tnoise = t.noise)
    x.test.df <- as.data.frame(data.test$x); y.test <- data.test$y
    
    intervals.cqr <- conformal.split(rbind(x.train.df, x.cal.df), c(y.train, y.cal), x.test.df, split = 1:n.train, pretrained.model = CQR.train, 
                                     method = "CQR", y_model = "RF")$intervals
    
    intervals.mat <- rbind(intervals.mat, intervals.cqr)
    x.test.vec <- c(x.test.vec, data.test$x)
    y.test.vec <- c(y.test.vec, data.test$y)
    y.low.oracle.vec <- c(y.low.oracle.vec, data.test$y_lo)
    y.up.oracle.vec <- c(y.up.oracle.vec, data.test$y_up)
    
  }
  
  x.test.mat <- cbind(x.test.mat, x.test.vec)
  y.test.mat <- cbind(y.test.mat, y.test.vec)
  y.low.oracle.mat <- cbind(y.low.oracle.mat, y.low.oracle.vec)
  y.up.oracle.mat <- cbind(y.up.oracle.mat, y.up.oracle.vec)
  intervals.mat.all <- cbind(intervals.mat.all, intervals.mat)
  
  
  if(plot.conditional){
    x.test <- x.test.df$V1
    o <- order(x.test)
    first.ind <- 1:((R-1)*n.test)
    #titles = c("MCP", "CQR", "DCP", "CB")
    plot.oracle = TRUE
    plot(x.test, y.test, xlab = "x", ylab = "y", main = paste0("alpha_fit = " , alpha.fit[i]))
    lines(x.test[o], intervals.mat[-first.ind,1][o], col = "blue" )
    lines(x.test[o], intervals.mat[-first.ind,2][o], col = "blue" )
    
    if (plot.oracle){
      lines(x.test[o], data.test$y_up[o], col = "purple")
      lines(x.test[o], data.test$y_lo[o], col = "purple")
    }
    
  } else{
    coverage.mat <- cbind(coverage.mat, intervals.mat[,1] <= y.test.vec & y.test.vec <= intervals.mat[,2])
    length.mat <- cbind(length.mat, intervals.mat[,2] - intervals.mat[,1])
  }
  
}


# Plot the different plots for some values of alpha
alphas.plot <- c(1,3,5,7)

plot.data <- data.frame(
  Y = rep(y.test.mat[1:1000,1], 2*length(alphas.plot)),
  X = rep(x.test.mat[1:1000,1], 2*length(alphas.plot)),
  y.low.oracle = rep(y.low.oracle.mat[1:1000, 1], 2*length(alphas.plot)),
  y.up.oracle = rep(y.up.oracle.mat[1:1000, 1], 2*length(alphas.plot)),
  y.low = c(intervals.mat.all[1:1000, 2*(alphas.plot - 1) + 1]),
  y.up = c(intervals.mat.all[1:1000, 2*(alphas.plot - 1) + 2]),
  class = as.factor(rep(paste0("Alpha = " ,alpha.fit[alphas.plot]), each = 1000))
)

ggplot(data = plot.data, aes(x = X, y = Y)) + 
  geom_ribbon(aes(ymin = y.low, ymax = y.up), fill = "lightblue") +
  geom_point(shape = 1) + 
  geom_line(aes(y=y.up), col = "blue") + 
  geom_line(aes(y = y.low), col = "blue") +
  geom_line(aes(y = y.low.oracle), col = "red") +
  geom_line(aes(y = y.up.oracle), col = "red") +
  theme_light() + 
  facet_wrap(~class) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ggsave("figs/CQR_1D_varying_alpha.pdf",width = 12,height = 8)


## Plot what happens to 
cond.cov.mat <- binning(x.test.vec, coverage.mat, 20)
cond.len.mat <- binning(x.test.vec, length.mat, 20)

data.plot <- matrix(NA, nrow = 20*8, ncol = 3)
data.plot[,1] <- rep(cond.cov.mat$xs, 8)
data.plot[,2] <- as.numeric(cond.cov.mat$cond)
data.plot <- as.data.frame(data.plot)
data.plot[,3] <- unlist(lapply(alpha.fit, function(x) rep(as.character(x), (20))))

names(data.plot) <- c("x", "y", "Alpha")

plot1 <- data.plot %>% ggplot( aes(x=x, y = y, group = Alpha, color = Alpha)) +
  geom_line() +
  geom_point(size=5, colour="white") + 
  geom_point(size=2, shape = 1) +
  theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("X") +
  ylab("Coverage") +
  geom_hline(yintercept = 0.9, linetype = "dashed", col = "black") + 
  theme(legend.position = "none")
plot1


data.plot <- matrix(NA, nrow = 20*8, ncol = 3)
data.plot[,1] <- rep(cond.len.mat$xs, 8)
data.plot[,2] <- as.numeric(cond.len.mat$cond)
data.plot <- as.data.frame(data.plot)
data.plot[,3] <- unlist(lapply(alpha.fit, function(x) rep(as.character(x), (20))))

names(data.plot) <- c("x", "y", "Alpha")

plot2 <- ggplot(data = data.plot) +
  geom_line(aes(x=x, y = y, group = Alpha, color = Alpha)) +
  geom_point(aes(x=x, y = y, group = Alpha, color = Alpha),size=5, colour="white") + 
  geom_point(aes(x=x, y = y, group = Alpha, color = Alpha),size=2, shape = 1) +
  theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("X") +
  ylab("Length") +geom_line(aes(x, y), data.frame(x = data.test$x, y = data.test$y_lo- data.test$y_up), color = alpha("black", 0.5), linetype = "dashed", size = 1)

plot2

plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 2)
#ggsave("figs/coverage_length_CQR_different_alphas.pdf", width = 12, height = 4, plot)
