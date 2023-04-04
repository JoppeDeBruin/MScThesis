# THis script reprocudes the graphs visualizing how conformal prediction works from Chapter 2
rm(list = ls())
library(ggplot2)
#setwd("~/ETH Zurich/MSc Thesis/thesis")

source("utils.R")

set.seed(1)
alpha = 0.1
data <- generate.data.1D(1000, homoscedastic = F, tnoise = FALSE,c_sd = 0.25 )

plot(data$x, data$y)
o <- order(data$x)
lines(data$x[o], data$y_noiseless[o], col = "green")

data.test <- generate.data.1D(1000, homoscedastic = F, c_sd = 0.25)

res <- conformal.split(data$x, data$y, data.frame(2.57), method = "mean", y_model = "boosting")

preds.test <- predict(res$model, data.frame(V1 = matrix(data$x, ncol = 1)), n.trees = 100)
preds.test.point <- predict(res$model, data.frame(V1 = matrix(2.57, ncol = 1)), n.trees = 100)
abs_resid <- abs(data$y[-res$split] - predict(res$model, data.frame( V1 = matrix((data$x[-res$split]), ncol = 1)), n.trees = 100))
all_resid <- c(abs_resid, -abs_resid)
hist(all_resid, breaks = 50)

quantile(abs_resid, 0.9)
mean( res$intervals[,1] <= data.test$y & data.test$y <= res$intervals[,2])

## Plot the results
data.plot <- data.frame(
  Y = data$y,
  X = data$x,
  y.pred = preds.test,
  y.up = res$intervals[,2]
)
data.plot$train = "cal"
data.plot$train[res$split] = "train"

x.test.point <- 2.57
x.test.point.index = 1

plot <- ggplot(data = data.plot, aes(x = X, y = Y, group = train)) + 
  geom_point(shape = 1, aes(col = train)) + 
  scale_color_manual(values = c(alpha("black", 0.5), alpha("red", 0.5)))+
  geom_line(aes(y = y.pred), col = "blue") + 
  theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "none")  + 
  scale_y_continuous(limits = c(-5, 15)) + 
  annotate("errorbar", x = x.test.point, y = preds.test[x.test.point.index], 
           ymin = res$intervals[x.test.point.index,1], 
           ymax = res$intervals[x.test.point.index,2],
           colour = "black", size = 1) +
  geom_point(aes(x = x.test.point, y = preds.test.point[x.test.point.index]), col ="black", size = 3)#, linewidth = 1)

plot.hist <- ggplot(data = data.frame(X=all_resid + preds.test.point[x.test.point.index]), aes(x = X)) +
  geom_histogram(col = alpha("red", 1), fill = alpha("red", 1), binwidth = 0.25) +
  theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(aes(xintercept = preds.test.point[x.test.point.index] + quantile(abs_resid, (1-alpha)*(1+1/length(abs_resid))))) + 
  geom_vline(aes(xintercept = preds.test.point[x.test.point.index] - quantile(abs_resid, (1-alpha)*(1+1/length(abs_resid))))) + 
  scale_x_continuous(limits = c(-5, 15)) +
  theme(axis.title.y = element_blank()) 
  
plot.hist <- plot.hist + coord_flip()

final_plot <- cowplot::plot_grid(plot,  plot.hist, rel_widths = c(1,1))
final_plot
#ggsave("figs/visualization_conformal_introduction.pdf", width = 8, height = 4, final_plot)


################# ALSO PLOT FOR FULL CONFORMAL #################


#set.seed(1234)
data.train <- data#generate.data.1D(200, homoscedastic = F, tnoise = F, c_sd = 0.25)

y_grid <- seq(-5, 15, length = 1000)
y_grid_include_q = y_grid_include_pi = rep(FALSE, length(y_grid))
pi_vec <- rep(NA, length(y_grid))
for (i in 1:length(y_grid)){
  if (i%%20 == 0){cat(i, "..")}
  model_y <- gbm::gbm(Y ~ . , data.frame(Y = c(data.train$y, y_grid[i]), c(data.train$x, 2.57)), 
                      n.trees = 100, distribution = "gaussian")
  
  resids <- abs(c(data.train$y, y_grid[i]) - predict(model_y, data.frame(c(data.train$x, 2.57)), n.trees = 100))
  
  q <- quantile(resids[-length(resids)], (1-alpha)*(1 + 1/length(data.train$y+1)))  
  
  if (resids[length(resids)] <= q){
    y_grid_include_q[i] = TRUE
  }
  
  pi <- 1/length(resids)*sum(resids <= resids[length(resids)])
  pi_vec[i] <- pi
  
  if (pi <= 1/length(resids) * ceiling((1-alpha)*length(resids))){
    y_grid_include_pi[i] = TRUE
  }
    
}


## Plot the results
test <- data.frame(c(data.test$x, 2.57))
names(test) <- names(data.frame(c(data.train$x, 2.57)))
preds.test <- predict(model_y, test, n.trees = 100)
data.plot <- data.frame(
  Y = c(data.test$y, 2.57^2 - 0.25*2.57^3),
  X = c(data.test$x, 2.57),
  y.pred =  preds.test#predict(model_y, data.frame(c(data.test$x, 2.57)), n.trees = 100)
  #y.up = res$intervals[,2]
)

x.test.point <- 2.57
#x.test.point.index = length(resids)

plot <- ggplot(data = data.plot, aes(x = X, y = Y)) + 
  geom_point(shape = 1, col = alpha("black", 0.5)) + 
  geom_line(aes(y = y.pred), col = "blue") + 
  theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "none")  + 
  scale_y_continuous(limits = c(-5, 15)) + 
  annotate("errorbar", x = x.test.point, y = preds.test[x.test.point.index], 
           ymin = min(y_grid[y_grid_include_pi]), 
           ymax = max(y_grid[y_grid_include_pi]),
           colour = "black", size = 1) 
 

plot.hist <- ggplot(data = data.frame(X = y_grid, Y = 1-pi_vec), aes(x = X, y = Y)) +
  geom_ribbon(aes(ymin = rep(0, length(y_grid)), ymax = Y), fill = "red") +
  geom_line(aes(y = Y), col = "red") + 
  geom_hline(aes(yintercept = 1- ceiling((1-alpha)*(length(data.train$y)+1))/(length(data.train$y+1))), linetype = "dotted") +   
  geom_vline(aes(xintercept = min(y_grid[y_grid_include_pi]))) +
  geom_vline(aes(xintercept = max(y_grid[y_grid_include_pi]))) + 
  geom_line(aes(y = rep(0, length(y_grid))), col = "red") +
  scale_x_continuous(limits = c(-5, 15)) +
  theme(axis.title.y = element_blank()) +
  theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_blank()) 

plot.hist <- plot.hist + coord_flip()

final_plot <- cowplot::plot_grid(plot,  plot.hist)
final_plot

#ggsave("figs/visualization_conformal_introduction_full.pdf", width = 8, height = 4, final_plot)
    
