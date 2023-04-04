# This script reproduces the simple 1D results from Chapter 4

source("utils.R")

R = 1000
cov.cqr = cov.mcp = cov.dcp =cov.naive = cov.qr = cov.cb <- rep(NA, R)
opt.alpha = rep(NA, R)

set.seed(1)

retrain = FALSE
if (!retrain){
  train.data <- generate.data.1D(n.train, homoscedastic = homoscedastic, tnoise = t.noise)
  x.train <- train.data$x; y.train <- train.data$y
  x.train.df <- as.data.frame(x.train)  
  
  CB.train <- rbart(x.train.df, y.train, ntree = 50, ntreeh = 20) 
  
  CQR.train <- quantile_rf(x.train.df, y.train, x.train.df, quantiles = c(alpha.fit/2, 1-alpha.fit/2))$model
  DCP.train <- quantile_rf(x.train.df, y.train, x.train.df, quantiles = seq(0.001, 0.999, length = distr.grid))$model
  MCP.train <- rf(x.train.df, y.train, x.train.df)$model
}


one_iteration <- function(r, retrain = FALSE, CQR.train = NULL, MCP.train = NULL, DCP.train = NULL, CB.train = NULL, seed = FALSE){
  
  if(seed){set.seed(r)}
  cat(r, "..")
  
  data.calibration <- generate.data.1D(n.cal, homoscedastic = homoscedastic, tnoise = t.noise)
  x.cal.df <- as.data.frame(data.calibration$x); y.cal <- data.calibration$y
  
  data.test <- generate.data.1D(n.test, homoscedastic = homoscedastic, tnoise = t.noise)
  x.test.df <- as.data.frame(data.test$x); y.test <- data.test$y
  
  if (retrain){
    train.data <- generate.data.1D(n.train, homoscedastic = homoscedastic, tnoise = t.noise)
    x.train <- train.data$x; y.train <- train.data$y
    x.train.df <- as.data.frame(x.train)  
    
    CB.train <- rbart(x.train.df, y.train, ntree = 50, ntreeh = 20) 
    
    CQR.train <- quantile_rf(x.train.df, y.train, x.test.df, quantiles = c(alpha.fit/2, 1-alpha.fit/2))$model
    DCP.train <- quantile_rf(x.train.df, y.train, x.test.df, quantiles = seq(0.001, 0.999, length = distr.grid))$model
    MCP.train <- rf(x.train.df, y.train, x.test.df)$model
    
  }
  
  intervals.cqr <- conformal.split(rbind(x.train.df, x.cal.df), c(y.train, y.cal), x.test.df, split = 1:n.train, 
                                   method = "CQR", y_model = "RF", pretrained.model = CQR.train)$intervals
  intervals.mcp <- conformal.split(rbind(x.train.df, x.cal.df), c(y.train, y.cal), x.test.df, split = 1:n.train, 
                                   method = "mean", y_model = "RF", pretrained.model = MCP.train)$intervals
  intervals.dcp <- conformal.split(rbind(x.train.df, x.cal.df), c(y.train, y.cal), x.test.df, split = 1:n.train, 
                                   method = "DCP", y_model = "RF", pretrained.model = DCP.train)$intervals
  intervals.cb <- conformal.split(rbind(x.train.df, x.cal.df), c(y.train, y.cal), x.test.df, split = 1:n.train, 
                                  method = "CB", y_model = "bart_het")$intervals
  intervals.naive <- naive.QR(rbind(x.train.df, x.cal.df), c(y.train, y.cal), x.test.df, split = 1:n.train, 
                              y_model = "RF", pretrained.model = DCP.train)#$intervals
  
  opt.alpha[i] = intervals.naive$alpha.opt*2
  intervals.naive <- intervals.naive$intervals
  
  intervals.qr <- predict(CQR.train, x.test.df)$predictions
  
  return(list(cqr = mean(intervals.cqr[,1] <= y.test & y.test <= intervals.cqr[,2]),
              mcp = mean(intervals.mcp[,1] <= y.test & y.test <= intervals.mcp[,2]),
              dcp = mean(intervals.dcp[,1] <= y.test & y.test <= intervals.dcp[,2]),
              naive = mean(intervals.naive[,1] <= y.test & y.test <= intervals.naive[,2]),
              qr = mean(intervals.qr[,1] <= y.test & y.test <= intervals.qr[,2]),
              cb = mean(intervals.cb[,1] <= y.test & y.test <= intervals.cb[,2])))
}

parallel = TRUE
R=1000

set.seed(1)
start = Sys.time()

if (parallel){
  res <- pbmclapply(1:R, one_iteration, retrain = FALSE, CQR.train = CQR.train,
                    DCP.train = DCP.train, MCP.train = MCP.train, CB.train = CB.train, seed = TRUE,
                    mc.cores = detectCores()[1]-1)  
} else{
  res <- pbmclapply(1:R, one_iteration,retrain = FALSE, CQR.train = CQR.train,
                    DCP.train = DCP.train, MCP.train = MCP.train, CB.train = CB.train, seed = TRUE,
                    mc.cores = 1)
}
print(Sys.time()-start)

# Construct the results matrices from the list return
cov.cqr <- res %>% map(1) %>% do.call(rbind, .)
cov.mcp <- res %>% map(2) %>% do.call(rbind, .)
cov.dcp <- res %>% map(3) %>% do.call(rbind, .)
cov.naive <- res %>% map(4) %>% do.call(rbind, .)
cov.qr <- res %>% map(5) %>% do.call(rbind, .)
cov.cb <- res %>% map(6) %>% do.call(rbind, .)
#save(cov.cqr, cov.mcp, cov.dcp, cov.naive, cov.qr, cov.cb, file = paste0("results/1D_sim_all_methods_", Sys.Date(), ".Rdata"))


par(mfrow = c(1,1))
hist(cov.cqr, breaks = 50, probability = TRUE)
xs <- seq(0.82, 0.95, length = 2000)
a = n.cal + 1 - floor((n.cal+1) * alpha)
b = floor((n.cal+1) * alpha)
lines(xs, dbeta(xs, a, b), col = "red")
lines((1:2000)/2000, 2000*VGAM::dbetabinom.ab(1:2000, n.test, shape1 = a, shape2 = b, log = FALSE), col = "blue")

# Histogram 1
data.plot <- data.frame(matrix(cov.cqr.1, ncol = 1))
data.plot[,2] <- c(rep("CQR", R))
names(data.plot) <- c("x", "Method")
plot1 <- ggplot() +
  geom_histogram(data = data.plot, aes(x = x, y = ..density.., group = Method, colour = Method, fill = Method), position = "identity", alpha = 0.5)+
  scale_color_manual(values = c("#619CFF"))+
  scale_fill_manual(values = c("#619CFF")) +
  geom_line(data = data.frame(x = xs, y = dbeta(xs, a, b)), aes(x = x, y = y), col = "orange", size = 1) +
  theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  xlab("Coverage") +
  ylab("Density") + 
  theme(legend.position = c(0.2,0.9)) + 
  geom_vline(xintercept = 0.9, linetype = "dashed")+
  ylim(c(0, 45))
plot1
#
# Histogram 2
data.plot <- data.frame(rbind(matrix(cov.naive.1, ncol = 1), matrix(cov.qr.1, ncol = 1)))
data.plot[,2] <- c(rep("Naive cal. QR", R), rep("Naive QR", R))
names(data.plot) <- c("x", "Method")
plot2 <- ggplot() +
  geom_histogram(data = data.plot, aes(x = x, y = ..density.., group = Method, colour = Method, fill = Method),
                 position = "identity", alpha = 0.5, bins = 31) + 
  scale_color_manual(values = c(  "#619CFF", "#00BA38"))+
  scale_fill_manual(values = c(  "#619CFF", "#00BA38")) +
  geom_line(data = data.frame(x = xs, y = dbeta(xs, a, b)), aes(x = x, y = y), col = "orange", size = 1) +
  theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))+
  xlab("Coverage") +
  ylab("Density") + 
  theme(legend.position = c(0.2,0.9)) + 
  geom_vline(xintercept = 0.9, linetype = "dashed") + 
  ylim(c(0, 45))
plot2

plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 2)
# ggsave("figs/hist_cov_1D.pdf", width = 12, height = 4, plot)

# Histograms for all the methods

data.plot <- data.frame(matrix(c(cov.mcp, cov.cqr, cov.dcp, cov.cb), ncol = 1))
data.plot[,2] <- c(rep("MCP", R), rep("CQR", R), rep("DCP", R), rep("CB", R))
names(data.plot) <- c("x", "Method")
plot1 <- ggplot() +
  geom_histogram(data = data.plot, aes(x = x, y = ..density.., col = "#619CFF", fill = "#619CFF"), position = "identity", alpha = 0.5)+
  facet_wrap(. ~ Method, ncol = 2) +
  scale_color_manual(values = c("#619CFF"))+
  scale_fill_manual(values = c("#619CFF")) +
  geom_line(data = data.frame(x = xs, y = dbeta(xs, a, b)), aes(x = x, y = y), col = "orange", size = 1) +
  theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  xlab("Coverage") +
  ylab("Density") + 
  theme(legend.position = "none") + 
  geom_vline(xintercept = 0.9, linetype = "dashed")+
  ylim(c(0, 30))
plot1
#ggsave("figs/hist_cov_1D_all.pdf", width = 12, height = 8, plot1)


