# This script produces the plots for the airfoil simulations in chapter 5

# First load the results! (created in simulation_airfoil_results.R)

### Plots for the WCP results

data.plot <- data.frame(rbind(cbind(as.numeric(res.mcp.rf$coverage), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.mcp.boost$coverage), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.mcp.bart$coverage), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.cqr.rf$coverage), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.cqr.boost$coverage), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.cqr.bart$coverage), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.dcp.rf$coverage), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.dcp.boost$coverage), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.cb$coverage), unlist(lapply(1:8, rep, 1000)))))
data.plot[,2] <- as.factor(data.plot[,2])
data.plot[,3] <- factor(c(rep("RF", 8000), rep("boost", 8000), rep("BART", 8000),
                          rep("RF", 8000), rep("boost", 8000), rep("BART", 8000),
                          rep("RF", 8000), rep("boost", 8000),
                          rep("BART", 8000)), levels = c("RF", "boost", "BART"))
data.plot[,4] <- factor(c(rep("MCP", 3*8000), rep("CQR", 3*8000), rep("DCP", 2*8000), rep("CB", 8000)),
                        levels = c("MCP", "CQR", "DCP", "CB"))
names(data.plot) <- c("cov", "x", "Regressor", "Conformal")

plot1 <- ggplot(data = data.plot, aes(x = x, y = cov)) +
  geom_boxplot() +
  geom_hline(data = data.frame(Conformal = c("MCP", "CQR", "DCP", "MCP", "CQR", "DCP", "MCP", "CQR", "CB"), 
                               Regressor = c("RF", "RF", "RF", "boost", "boost", "boost", "BART", "BART", "BART"), 
                               Z = rep(0.9, 9)), aes(yintercept = Z), col = "blue")+
  facet_grid(Regressor ~ Conformal) +
  scale_x_discrete(limits = rev(levels(data.plot$x)))+
  coord_flip() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  ylab("Empirical Coverage") + xlab("Experiment")
plot1

#ggsave("figs/boxplots_airfoil_coverage.pdf",width = 9,height = 9)

data.plot <- data.frame(rbind(cbind(as.numeric(res.mcp.rf$length.mean), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.mcp.boost$length.mean), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.mcp.bart$length.mean), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.cqr.rf$length.mean), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.cqr.boost$length.mean), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.cqr.bart$length.mean), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.dcp.rf$length.mean), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.dcp.boost$length.mean), unlist(lapply(1:8, rep, 1000))),
                              cbind(as.numeric(res.cb$length.mean), unlist(lapply(1:8, rep, 1000)))))
data.plot[,2] <- as.factor(data.plot[,2])
data.plot[,3] <- factor(c(rep("RF", 8000), rep("boost", 8000), rep("BART", 8000),
                          rep("RF", 8000), rep("boost", 8000), rep("BART", 8000),
                          rep("RF", 8000), rep("boost", 8000),
                          rep("BART", 8000)), levels = c("RF", "boost", "BART"))
data.plot[,4] <- factor(c(rep("MCP", 3*8000), rep("CQR", 3*8000), rep("DCP", 2*8000), rep("CB", 8000)),
                        levels = c("MCP", "CQR", "DCP", "CB"))
names(data.plot) <- c("cov", "x", "Regressor", "Conformal")

plot2 <- ggplot(data = data.plot, aes(x = x, y = cov)) +
  geom_boxplot() +
  facet_grid(Regressor ~ Conformal) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  ylab("Length") + xlab("Experiment") +
  coord_cartesian(ylim = c(5,30))
 
plot2
#ggsave("figs/boxplots_airfoil_length.pdf",width = 9,height = 6)


### Histograms of coverage for WCP (just one example)

# Histogram 1
data.plot <- data.frame(matrix(c(res.dcp.rf$coverage[,1], res.dcp.rf$coverage[,3]), ncol = 1))
data.plot[,2] <- c(rep("Unweighted", 1000), rep("Weighted", 1000))
names(data.plot) <- c("x", "Method")
plot1 <- ggplot() +
  geom_histogram(data = data.plot, aes(x = x, y = ..density.., group = Method, colour = Method, fill = Method), 
                 position = "identity", alpha = 0.6)+
  scale_color_manual(values = c("#619CFF", "#00BA38"))+
  scale_fill_manual(values = c("#619CFF", "#00BA38")) +
  theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  xlab("Empirical Coverage") +
  ylab("Density") + 
  theme(legend.position = c(0.2,0.9)) + 
  geom_vline(xintercept = 0.9, linetype = "dashed") +
  ylim(c(0, 25))
  
plot1

data.plot <- data.frame(matrix(c(res.dcp.rf$coverage[,3], res.dcp.rf$coverage[,4]), ncol = 1))
data.plot[,2] <- c(rep("Weighted", 1000), rep("Unweighted, smaller sample", 1000))
names(data.plot) <- c("x", "Method")
plot2 <- ggplot() +
  geom_histogram(data = data.plot, aes(x = x, y = ..density.., group = Method, colour = Method, fill = Method), 
                 position = "identity", alpha = 0.6)+
  scale_color_manual(values = c("#619CFF", "#00BA38"))+
  scale_fill_manual(values = c("#619CFF", "#00BA38")) +
  theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  xlab("Empirical Coverage") +
  ylab("Density") + 
  theme(legend.position = c(0.2,0.9)) + 
  geom_vline(xintercept = 0.9, linetype = "dashed") +
  ylim(c(0, 25))
plot2

plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 2)
#ggsave("figs/histograms_coverage_dcp_airfoil.pdf",width = 12,height = 4, plot)

### Boxplots for the DRP results
data.plot <- data.frame(rbind(cbind(as.numeric(res.drp.cqr.boost_r.rf_m$coverage[,-4]), unlist(lapply((1:8)[-4], rep, 1000))),
                        cbind(as.numeric(res.drp.cqr.boost_r.boost_m$coverage[,-4]), unlist(lapply((1:8)[-4], rep, 1000))),
                        cbind(as.numeric(res.drp.cqr.rf_r.rf_m$coverage[,-4]), unlist(lapply((1:8)[-4], rep, 1000))),
                        cbind(as.numeric(res.drp.cqr.rf_r.boost_m$coverage[,-4]), unlist(lapply((1:8)[-4], rep, 1000)))))

data.plot[,2] <- as.factor(data.plot[,2])
data.plot[,3] <- factor(c(rep("R_boost", 7000), rep("R_boost", 7000), rep("R_RF", 7000), rep("R_RF", 7000)))
data.plot[,4] <- as.factor(c(rep("M_RF", 7000), rep("M_boost", 7000), rep("M_RF", 7000), rep("M_boost", 7000)))
names(data.plot) <- c("cov", "x", "R", "M")

plot1 <- ggplot(data = data.plot, aes(x = x, y = cov)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.9, col = "blue")+
  facet_grid(R~M) +
  scale_x_discrete(limits = rev(levels(data.plot$x)))+
  coord_flip(ylim = c(0.6, 1)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  ylab("Empirical Coverage") + xlab("Experiment") +
  ggtitle("CQR")
  
plot1


data.plot <- data.frame(rbind(cbind(as.numeric(res.drp.mean.boost_r.rf_m$coverage[,-4]), unlist(lapply((1:8)[-4], rep, 1000))),
                              cbind(as.numeric(res.drp.mean.boost_r.boost_m$coverage[,-4]), unlist(lapply((1:8)[-4], rep, 1000))),
                              cbind(as.numeric(res.drp.mean.rf_r.rf_m$coverage[,-4]), unlist(lapply((1:8)[-4], rep, 1000))),
                              cbind(as.numeric(res.drp.mean.rf_r.boost_m$coverage[,-4]), unlist(lapply((1:8)[-4], rep, 1000)))))

data.plot[,2] <- as.factor(data.plot[,2])
data.plot[,3] <- factor(c(rep("R_boost", 7000), rep("R_boost", 7000), rep("R_RF", 7000), rep("R_RF", 7000)))
data.plot[,4] <- as.factor(c(rep("M_RF", 7000), rep("M_boost", 7000), rep("M_RF", 7000), rep("M_boost", 7000)))
names(data.plot) <- c("cov", "x", "R", "M")

plot2 <- ggplot(data = data.plot, aes(x = x, y = cov)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.9, col = "blue")+
  facet_grid(R~M) +
  scale_x_discrete(limits = rev(levels(data.plot$x)))+
  coord_flip(ylim = c(0.6, 1)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  ylab("Empirical Coverage") + xlab("Experiment") +
  ggtitle("MCP")
  
  
plot2

plot <- gridExtra::grid.arrange(plot2, plot1, ncol = 2)

#ggsave("figs/boxplots_airfoil_coverage_DRP.pdf",width = 9,height = 6, plot)



### Now plot the lengths
data.plot <- data.frame(rbind(cbind(as.numeric(res.drp.cqr.boost_r.rf_m$length.mean[,-4]), unlist(lapply((1:8)[-4], rep, 1000))),
                              cbind(as.numeric(res.drp.cqr.boost_r.boost_m$length.mean[,-4]), unlist(lapply((1:8)[-4], rep, 1000))),
                              cbind(as.numeric(res.drp.cqr.rf_r.rf_m$length.mean[,-4]), unlist(lapply((1:8)[-4], rep, 1000))),
                              cbind(as.numeric(res.drp.cqr.rf_r.boost_m$length.mean[,-4]), unlist(lapply((1:8)[-4], rep, 1000)))))

data.plot[,2] <- as.factor(data.plot[,2])
data.plot[,3] <- factor(c(rep("R_boost", 7000), rep("R_boost", 7000), rep("R_RF", 7000), rep("R_RF", 7000)))
data.plot[,4] <- as.factor(c(rep("M_RF", 7000), rep("M_boost", 7000), rep("M_RF", 7000), rep("M_boost", 7000)))
names(data.plot) <- c("cov", "x", "R", "M")

plot1 <- ggplot(data = data.plot, aes(x = x, y = cov)) +
  geom_boxplot() +
  facet_grid(R~M) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  ylab("Length") + xlab("Experiment") +
  ggtitle("CQR")

plot1


data.plot <- data.frame(rbind(cbind(as.numeric(res.drp.mean.boost_r.rf_m$length.mean[,-4]), unlist(lapply((1:8)[-4], rep, 1000))),
                              cbind(as.numeric(res.drp.mean.boost_r.boost_m$length.mean[,-4]), unlist(lapply((1:8)[-4], rep, 1000))),
                              cbind(as.numeric(res.drp.mean.rf_r.rf_m$length.mean[,-4]), unlist(lapply((1:8)[-4], rep, 1000))),
                              cbind(as.numeric(res.drp.mean.rf_r.boost_m$length.mean[,-4]), unlist(lapply((1:8)[-4], rep, 1000)))))

data.plot[,2] <- as.factor(data.plot[,2])
data.plot[,3] <- factor(c(rep("R_boost", 7000), rep("R_boost", 7000), rep("R_RF", 7000), rep("R_RF", 7000)))
data.plot[,4] <- as.factor(c(rep("M_RF", 7000), rep("M_boost", 7000), rep("M_RF", 7000), rep("M_boost", 7000)))
names(data.plot) <- c("cov", "x", "R", "M")

plot2 <- ggplot(data = data.plot, aes(x = x, y = cov)) +
  geom_boxplot() +
  facet_grid(R~M) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  ylab("Length") + xlab("Experiment") +
  ggtitle("MCP")




plot2

plot <- gridExtra::grid.arrange(plot2, plot1, ncol = 2)

#ggsave("figs/boxplots_airfoil_length_DRP.pdf",width = 9,height = 6, plot)