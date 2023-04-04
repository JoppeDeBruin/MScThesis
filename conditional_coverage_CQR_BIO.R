# This script shows that conditional coverage can be lost for CQR on the CASP dataset

rm(list = ls())

library(grf)
library(randomForest)
library(gbm)
library(ggplot2)

source("utils.R")


data <- na.omit(read.csv(url("https://archive.ics.uci.edu/ml/machine-learning-databases/00265/CASP.csv")))

run.CASP.example <- function(alpha.fit){
  R = 5 # Repeat the experiment 5 times to smooth out some noise
  intervals.1 <- intervals.2 <- NULL
  x.vec <- NULL
  y.vec <- c()
  
  set.seed(1234)
  
  for (i in 1:R){
    cat(i, " .. ")
    # Sample 20000 data points randomly
    samples <- sample(nrow(data), 20000)
    
    x <- data[samples,2:10]
    y <- data[samples,1]
    
    # Split in train and test data
    train.ix <- sample(nrow(x), floor(0.6*nrow(x)))
    
    x.train <- x[train.ix, ]; x.test <- x[-train.ix,]
    y.train <- y[train.ix]; y.test <- y[-train.ix]
    
    ###
    #First: run MCP
    res1 <- conformal.split(x.train, y.train, x.test, y_model = "RF", method = "mean", 
                            regression.parameters = list(n.trees = 1000))
    #THen: run CQR with the specific value of alpha
    res2 <- conformal.split(x.train, y.train, x.test, y_model = "RF", method = "CQR", 
                            regression.parameters = list(n.trees = 1000), alpha.fit = alpha.fit)
    
    intervals.1 <- rbind(intervals.1, res1$intervals)
    intervals.2 <- rbind(intervals.2, res2$intervals)
    
    # We concatenate and add all the x and y values to one vector
    x.vec <- rbind(x.vec, x.test)
    y.vec <- c(y.vec, y.test)
    
    
  }  
  
  return(list("x.vec" = x.vec, "y.vec" = y.vec, "intervals.1" = intervals.1, "intervals.2" = intervals.2))
}

# Run the experiments for different values of alpha
res0.1 <- run.CASP.example(0.1)
res0.2 <- run.CASP.example(0.2)
res0.3 <- run.CASP.example(0.3)
res0.4 <- run.CASP.example(0.4)
res0.5 <-run.CASP.example(0.5)
res0.6 <-run.CASP.example(0.6)

#save(res0.1, res0.5, file = paste0("results/CASP_20000_", Sys.Date(), ".RData"))

# The minimum average length is achieved for alpha = 0.5. So compare 0.1 to 0.5

# plot all the graphs

res <- res0.5
intervals.1 <- res$intervals.1
intervals.2 <- res$intervals.2
y.vec <- res$y.vec
mean(intervals.1[,1] <= y.vec & y.vec <= intervals.1[,2])
mean(intervals.2[,1] <= y.vec & y.vec <= intervals.2[,2])

mean(intervals.1[,2]- intervals.1[,1])
mean(intervals.2[,2] - intervals.2[,1])

len1 <- intervals.1[,2] - intervals.1[,1]
len2 <- intervals.2[,2] - intervals.2[,1]

#pdf("figs/CASP_alpha_0dot5vs0dot1.pdf", height = 12, width = 10)
par(mfrow = c(4,2))
for (i in 1:8){
  # What is going on here? 
  res1.bins <- binning(res$x.vec[,i], cbind((intervals.1[,1] <= res$y.vec & res$y.vec <= intervals.1[,2]),(intervals.1[,1] <= res$y.vec & res$y.vec <= intervals.1[,2])), 20)
  res2.bins <- binning(res$x.vec[,i], cbind((intervals.2[,1] <= res$y.vec & res$y.vec <= intervals.2[,2]),(intervals.2[,1] <= res$y.vec & res$y.vec <= intervals.2[,2])), 20)
  res0.1.bins <- binning(res0.1$x.vec[,i], cbind((res0.1$intervals.2[,1] <= res0.1$y.vec & res0.1$y.vec <= res0.1$intervals.2[,2]),(res0.1$intervals.2[,1] <= res0.1$y.vec & res0.1$y.vec <= res0.1$intervals.2[,2])), 20)
  
  # The black line is the one using MCP vs CQR with certain alpha values
  plot(res1.bins$xs, res1.bins$cond[,1], type = "b", xlim = range(res1.bins$xs, res2.bins$xs), xlab = paste0("X", i), ylab = "Coverage", ylim = range(0.8, 1))
  lines(res2.bins$xs, res2.bins$cond[,1], type = "b", col = "blue", lwd = 2)
  lines(res0.1.bins$xs, res0.1.bins$cond[,1], type = "b", col = "#FFCC00", lwd = 2)
  abline(h = 0.9, col = "black", lty = 2)
  
}
#graphics.off()

##### Plot conditional coverage binned according to lenghts
bins = 10
res1.bins <- binning(res0.1$intervals.2[,2]-res0.1$intervals.2[,1], cbind((intervals.1[,1] <= y.vec & y.vec <= intervals.1[,2]),(intervals.2[,1] <= y.vec & y.vec <= intervals.2[,2]), (res0.1$intervals.2[,1] <= res0.1$y.vec & res0.1$y.vec <= res0.1$intervals.2[,2])), bins)


data.plot <- matrix(NA, nrow = bins*3, ncol = 3)
data.plot[,1] <- rep(res1.bins$xs, 3)
data.plot[,2] <- as.numeric(res1.bins$cond)
data.plot <- as.data.frame(data.plot)
data.plot[,3] <- unlist(lapply(c("MCP", "CQR 0.5", "CQR 0.1"), function(x) rep(as.character(x), (bins))))

names(data.plot) <- c("x", "y", "Method")

plot1 <- data.plot %>% ggplot( aes(x=x, y = y, group = Method, color = Method)) +
  geom_line() +
  geom_point(size=5, colour="white") + 
  geom_point(size=2, shape = 1) +
  theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Length (alpha = 0.1)") +
  ylab("Coverage") +
  geom_hline(yintercept = 0.9, linetype = "dashed", col = "black") + 
  #theme(legend.position = "none") +
  scale_color_manual(values = c("#FFCC00",
                                "blue",
                                "black"))
plot1


res2.bins <- binning(len2, cbind((intervals.1[,1] <= y.vec & y.vec <= intervals.1[,2]),(intervals.2[,1] <= y.vec & y.vec <= intervals.2[,2]), (res0.1$intervals.2[,1] <= res0.1$y.vec & res0.1$y.vec <= res0.1$intervals.2[,2])), bins)


data.plot <- matrix(NA, nrow = bins*3, ncol = 3)
data.plot[,1] <- rep(res2.bins$xs, 3)
data.plot[,2] <- as.numeric(res2.bins$cond)
data.plot <- as.data.frame(data.plot)
data.plot[,3] <- unlist(lapply(c("MCP", "CQR 0.5", "CQR 0.1"), function(x) rep(as.character(x), (bins))))

names(data.plot) <- c("x", "y", "Method")

plot2 <- data.plot %>% ggplot( aes(x=x, y = y, group = Method, color = Method)) +
  geom_line() +
  geom_point(size=5, colour="white") + 
  geom_point(size=2, shape = 1) +
  theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Length (alpha = 0.5)") +
  ylab("Coverage") +
  geom_hline(yintercept = 0.9, linetype = "dashed", col = "black") + 
  #theme(legend.position = "none") +
  scale_color_manual(values = c("#FFCC00",
                                "blue",
                                "black"))
plot2

plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 2)

#ggsave("figs/CASP_alpha_0dot5vs0dot1_length.pdf", width = 12, height = 4, plot)

