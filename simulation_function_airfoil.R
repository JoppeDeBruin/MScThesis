# This script defines the simulation functions for the airfoil simulations
run.airfoil.simulation <- function(conformal.function, X, Y, N, n_sim = 100, all = TRUE, seeds = FALSE, parallel = FALSE, ...){
  
  R = n_sim
  n = round(N/2)
  
  one_iteration <- function(r){
    cov.vec <- len.vec <- len.vec.mean <- rep(NA, 8)
    
    if (r %% 1 == 0) cat(r, "..")
    if (seeds) {set.seed(r)}
    i = sample(N,n)
    x = X[i,]; y = Y[i]
    x0 = X[-i,]; y0 = Y[-i]
    
    # No tilting
    out1 = conformal.function(x, y, x0, ...)
    
    cov.vec[1] = mean(out1$intervals[,1] <= y0 & y0 <= out1$intervals[,2])
    len.vec[1] = median(out1$intervals[,2] - out1$intervals[,1])
    len.vec.mean[1] = mean(out1$intervals[,2] - out1$intervals[,1])
    
    # Tilting
    i0 = weighted.sample(tilting.function(x0))
    x00 = x0[i0,]; y00 = y0[i0]
    i1 = out1$split; n1 = length(out1$split)
    ssize = sum(tilting.function(x[-i1,]))^2/sum(tilting.function(x[-i1,])^2)
    
    # Unweighted
    out2 = conformal.function(x, y, x00, ...)
    
    cov.vec[2] = mean(out2$intervals[,1] <= y00 & y00 <= out2$intervals[,2])
    len.vec[2] = median(out2$intervals[,2] - out2$intervals[,1])
    len.vec.mean[2] = mean(out2$intervals[,2] - out2$intervals[,1])
    
    # Weighted, oracle  
    out3 = conformal.function(x, y, x00, w = tilting.function(rbind(x, x00)), ...)
    
    cov.vec[3] = mean(out3$intervals[,1] <= y00 & y00 <= out3$intervals[,2])
    len.vec[3] = median(out3$intervals[,2] - out3$intervals[,1])
    len.vec.mean[3] = mean(out3$intervals[,2] - out3$intervals[,1])
    
    if (all){
      # No tilting, smaller sample size
      i2 = (1:n)[-i1]
      i3 = c(i1, sample(i2, floor(ssize)))
      
      out4 = conformal.function(x[i3, ], y[i3], x0, split = 1:n1, ...)
      
      cov.vec[4] = mean(out4$intervals[,1] <= y0 & y0 <= out4$intervals[,2])
      len.vec[4] = median(out4$intervals[,2] - out4$intervals[,1])
      len.vec.mean[4] = mean(out4$intervals[,2] - out4$intervals[,1])
      
      # Weighted, estimated with logistic regression
      zy = as.factor(c(rep(0,n),rep(1,nrow(x00))))
      zx = rbind(x,x00)
      obj.glm = glm(zy ~ zx, family="binomial")
      prob.glm = predict(obj.glm, type="response")
      wts.glm = prob.glm / (1-prob.glm)
      
      out5 = conformal.function(x, y, x00, w = wts.glm, ...)
      
      cov.vec[5] = mean(out5$intervals[,1] <= y00 & y00 <= out5$intervals[,2])
      len.vec[5] = median(out5$intervals[,2] - out5$intervals[,1])
      len.vec.mean[5] = mean(out5$intervals[,2] - out5$intervals[,1])
      
      # Weighted, estimated with random forest
      obj.rf = randomForest(zx, zy)
      prob.rf = predict(obj.rf, type="prob")[,2]
      prob.rf = pmax(pmin(prob.rf, 0.99), 0.01)
      wts.rf = prob.rf / (1-prob.rf)
      
      out6 = conformal.function(x, y, x00, w = wts.rf, ...)
      
      cov.vec[6] = mean(out6$intervals[,1] <= y00 & y00 <= out6$intervals[,2])
      len.vec[6] = median(out6$intervals[,2] - out6$intervals[,1])
      len.vec.mean[6] = mean(out6$intervals[,2] - out6$intervals[,1])
      ratio.inf <- mean((out6$intervals[,2] - out6$intervals[,1])==Inf)
      #print(mean(is.infinite(out6$intervals[,2] - out6$intervals[,1])))
      
      # No tilting, weighted, estimated with logistic
      zy = as.factor(c(rep(0,n),rep(1,nrow(x0))))
      zx = rbind(x,x0)
      obj.glm = glm(zy ~ zx, family="binomial")
      prob.glm = predict(obj.glm, type="response")
      wts.glm = prob.glm / (1-prob.glm)
      
      out7 = conformal.function(x, y, x0, w = wts.glm, ...)
      
      cov.vec[7] = mean(out7$intervals[,1] <= y0 & y0 <= out7$intervals[,2])
      len.vec[7] = median(out7$intervals[,2] - out7$intervals[,1])
      len.vec.mean[7] = mean(out7$intervals[,2] - out7$intervals[,1])
      
      # Weighted, estimated with random forest
      obj.rf = randomForest(zx, zy)
      prob.rf = predict(obj.rf, type="prob")[,2]
      prob.rf = pmax(pmin(prob.rf, 0.99), 0.01)
      wts.rf = prob.rf / (1-prob.rf)
      
      out8 = conformal.function(x, y, x0, w = wts.rf, ...)
      
      cov.vec[8] = mean(out8$intervals[,1] <= y0 & y0 <= out8$intervals[,2])
      len.vec[8] = median(out8$intervals[,2] - out8$intervals[,1]) # Length is taken as the mean. Is that a sensible thing to do? Can we do stratified medians?
      len.vec.mean[8] = mean(out8$intervals[,2] - out8$intervals[,1])
    }
    if (all){
      return(list(cov = cov.vec, len = len.vec, len.mean = len.vec.mean, ssize = ssize, ratio.inf = ratio.inf))
    }
    return(list(cov = cov.vec, len = len.vec, len.mean = len.vec.mean, ssize = ssize))
  }
  
  if (parallel){
    res <- mclapply(1:R, one_iteration, mc.cores = detectCores()[1]-1)  
  } else{
    res <- mclapply(1:R, one_iteration, mc.cores = 1)
  }

  # Construct the results matrices from the list return
  cov.mat <- res %>% map(1) %>% do.call(rbind, .)
  len.mat <- res %>% map(2) %>% do.call(rbind, .)
  len.mat.mean <- res %>% map(3) %>% do.call(rbind, .)
  ssize <- res %>% map(4) %>% do.call(rbind, .)
  
  if (all==T) {
    ratio.inf <- res %>% map(5) %>% do.call(rbind, .)
    return(list(coverage = cov.mat, length = len.mat, length.mean = len.mat.mean, effective.ss = ssize, ratio.inf = ratio.inf))
  }
  
  return(list(coverage = cov.mat, length = len.mat, length.mean = len.mat.mean, effective.ss = ssize))
}

# Quick simulations to show what happens only for DRP and CQR with QRF
run.airfoil.simulation.RF.only <- function(conformal.function, X, Y, N, n_sim = 100, seeds = FALSE, weights = "RF", ...){
  
  params <- list(...)
  
  R = n_sim
  cov.mat = len.mat = len.mat.mean =  matrix(0,R,8)
  ssize = numeric(R)
  n = round(N/2)
  intervals.drp <- intervals.cqr <- NULL
  y.vec <- c()
  
  for (r in 1:R) {
    if (r %% 1 == 0) cat(r,".. ")
    if (seeds) {set.seed(r)}
    i = sample(N,n)
    x = X[i,]; y = Y[i]
    x0 = X[-i,]; y0 = Y[-i]
    
    # Tilting
    i0 = weighted.sample(tilting.function(x0))
    x00 = x0[i0,]; y00 = y0[i0]
    
    # Weighted, estimated with random forest
    zy = as.factor(c(rep(0,n),rep(1,nrow(x00))))
    zx = rbind(x,x00)
    if (weights == "RF"){
      obj.rf = randomForest(zx, zy)
      prob.rf = predict(obj.rf, type="prob")[,2]
      prob.rf = pmax(pmin(prob.rf, 0.99), 0.01)
      wts.rf = prob.rf / (1-prob.rf)
      w = wts.rf
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
    
    out6 = conformal.function(x, y, x00, w = w, ...)
    out8 <- DRP(rbind(x, x00), c(y, y00), c(rep(0, nrow(x)), rep(1, nrow(x00))), x00, method_pi = weights, 
                method = "RF", score = "CQR", method_m = "binary", model_m = "RF", regression.parameters = list(n.trees = 200))
    
    intervals.cqr <- rbind(intervals.cqr, out6$intervals)
    intervals.drp <- rbind(intervals.drp, out8$intervals)
    y.vec <- c(y.vec, y00)
    
  }
  
  return(list(intervals.cqr = intervals.cqr, intervals.drp = intervals.drp, y.vec = y.vec))
}

run.airfoil.simulation.DRP <- function(X, Y, N, n_sim = 100, seed = FALSE, oracle_weights = FALSE, parallel = FALSE, ...){
  
  R = n_sim
  n = round(N/2)
  
  one_iteration <- function(r){
    
    cov.vec <- len.vec <- len.vec.mean <- rep(NA, 8)
    cat(r, "..")
    if (seed) {set.seed(r)}
    i = sample(N,n)
    x = X[i,]; y = Y[i]
    x0 = X[-i,]; y0 = Y[-i]
    
    # First without shift (and we know about it)
    out1 <- DRP(rbind(x, x0), c(y, y0), c(rep(0, nrow(x)), rep(1, nrow(x0))), x0, w = function(x){rep(0.5, length(x))}, ...)
    cov.vec[1] = mean(out1$interval[,1] <= y0 & y0 <= out1$interval[,2])
    len.vec[1] <- median(out1$interval[,2]-out1$interval[,1])
    len.vec.mean[1] = mean(out1$interval[,2] - out1$interval[,1])
    
    # Tilting
    i0 = weighted.sample(tilting.function(x0))
    x00 = x0[i0,]; y00 = y0[i0]
    
    # Then with shift without weights
    out2 <- DRP(rbind(x, x00), c(y, y00), c(rep(0, nrow(x)), rep(1, nrow(x00))), x00, w = function(x){rep(0.5, length(x))} ,...)
    cov.vec[2] = mean(out2$interval[,1] <= y00 & y00 <= out2$interval[,2])
    len.vec[2] <- median(out2$interval[,2]-out2$interval[,1])
    len.vec.mean[2] = mean(out2$interval[,2] - out2$interval[,1])
    
    # Oracle weights
    out3 <- DRP(rbind(x, x00), c(y, y00), c(rep(0, nrow(x)), rep(1, nrow(x00))), x00, w = tilting.function ,...)
    cov.vec[3] = mean(out3$interval[,1] <= y00 & y00 <= out3$interval[,2])
    len.vec[3] <- median(out3$interval[,2]-out3$interval[,1])
    len.vec.mean[3] = mean(out3$interval[,2] - out3$interval[,1])
    
    # DOn't do 4 (with less samples)
    
    out5 <- DRP(rbind(x, x00), c(y, y00), c(rep(0, nrow(x)), rep(1, nrow(x00))), x00, method_pi = "glm", ...)
    cov.vec[5] = mean(out5$interval[,1] <= y00 & y00 <= out5$interval[,2])
    len.vec[5] <- median(out5$interval[,2]-out5$interval[,1])
    len.vec.mean[5] = mean(out5$interval[,2] - out5$interval[,1])
    
    out6 <- DRP(rbind(x, x00), c(y, y00), c(rep(0, nrow(x)), rep(1, nrow(x00))), x00, method_pi = "RF", ...)
    cov.vec[6] = mean(out6$interval[,1] <= y00 & y00 <= out6$interval[,2])
    len.vec[6] <- median(out6$interval[,2]-out6$interval[,1])
    len.vec.mean[6] = mean(out6$interval[,2] - out6$interval[,1])
    
    out7 <- DRP(rbind(x, x0), c(y, y0), c(rep(0, nrow(x)), rep(1, nrow(x0))), x0, method_pi = "glm" , ...)
    cov.vec[7] = mean(out7$interval[,1] <= y0 & y0 <= out7$interval[,2])
    len.vec[7] <- median(out7$interval[,2]-out7$interval[,1])
    len.vec.mean[7] = mean(out7$interval[,2] - out7$interval[,1])
    
    out8 <- DRP(rbind(x, x0), c(y, y0), c(rep(0, nrow(x)), rep(1, nrow(x0))), x0, method_pi = "RF" , ...)
    cov.vec[8] = mean(out8$interval[,1] <= y0 & y0 <= out8$interval[,2])
    len.vec[8] <- median(out8$interval[,2]-out8$interval[,1])
    len.vec.mean[8] = mean(out8$interval[,2] - out8$interval[,1])
    
    return(list(cov = cov.vec, len = len.vec, len.mean = len.vec.mean))
  }
  
  if (parallel){
    res <- pbmclapply(1:R, one_iteration, mc.cores = detectCores()[1] - 1)  
  } else{
    res <- mclapply(1:R, one_iteration, mc.cores = 1)
  }
  
  
  cov.mat <- res %>% map(1) %>% do.call(rbind, .)
  len.mat <- res %>% map(2) %>% do.call(rbind, .)
  len.mat.mean <- res %>% map(3) %>% do.call(rbind, .)
  
  return(list(coverage = cov.mat, length = len.mat, length.mean = len.mat.mean))
}


## Exponential tilting functions
tilting.function = function(x, coefs = c(-1,1)) {
  exp(x[,c(1,5)] %*% c(coefs[1], coefs[2])) # can potentially add other weighting schemes too (see the other paper)
}

weighted.sample = function(wts, frac=0.25) {
  n = length(wts)
  i = c()
  # only sample a part of the data set according to the tilting!
  while(length(i) <= n*frac) {
    i = c(i, which(runif(n) <= wts/max(wts)))
  }
  return(i)
}




