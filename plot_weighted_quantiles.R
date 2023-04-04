# Generate the plot for dummy example of the weighted quantiles in Chapter 5

scores = runif(1000, 0, 1)

w = 10*abs(scores)^3 + 1
#plot(scores, w)

weighted.quantile(scores, 0.9, w)
weighted.quantile(scores, 0.9)


alphas = seq(0.01, 0.99, length = 200)
vals <- vals.w <- rep(NA, length(alphas))
for (i in 1:length(alphas)){
  
  vals[i] <- weighted.quantile(scores, alphas[i])
  vals.w[i] <- weighted.quantile(scores, alphas[i], w)
}

plot(vals, alphas)
lines(vals.w, alphas)


data.plot <- data.frame(X = c(vals, vals.w), Y = c(alphas, alphas))
data.plot[,3] <- c(rep("Unweighted", length(alphas)), rep("Weighted", length(alphas)))
names(data.plot) <- c("X", "Y", "Method")

plot <- ggplot(data = data.plot, aes(y = Y, x = X, group = Method))+
  geom_line(aes(color = Method), size = 1) + 
  theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Conformity Score") +
  ylab("(Weighted) CDF") + 
  #theme(legend.position = c(0.2,0.9)) + 
  geom_hline(yintercept = 0.75, linetype = "dashed")+
  geom_vline(xintercept = vals[which.min(abs(alphas - 0.75))], linetype = "dashed", col = "#F8766D" ) +
  geom_vline(xintercept = vals.w[which.min(abs(alphas - 0.75))], linetype = "dashed", col = "#00BFC4")
plot

#ggsave("figs/weighted_quantiles.pdf", width = 6, height = 3, plot)


