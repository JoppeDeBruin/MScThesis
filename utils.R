# First define all base learners used in the conformal and double robust prediction approaches.

linear <- function(X.fit, Y.fit, X.test, X.cal = NULL, Y.cal = NULL, family = "regression", pretrained.model = NULL, ...){
  
  # Training
  if (is.null(pretrained.model)){
    if (family == "regression"){
      model <- lm(Y~., data.frame(Y = Y.fit, X.fit))
    } else if (family == "binary"){
      if (length(unique(Y.fit)) > 2){warning("Classification is attempted with more than 2 values for target! lin")}
      if (!is.factor(Y.fit)){
        Y.fit <- as.factor(Y.fit)
      }
      model <- glm(Y~., data.frame(Y = Y.fit, X.fit), family = "binomial") #randomForest::randomForest(X.fit, Y.fit, ntree = n.trees, ...)
    }
  } else{
    model <- pretrained.model
  }
  
  # Prediction
  if (family == "regression"){
    preds.test <- predict(model, data.frame(X.test))
    
    if (!is.null(X.cal) & !is.null(Y.cal)){
      preds.cal <- predict(model, data.frame(X.cal))
      return(list(cal = preds.cal, test = preds.test, model = model))
    }
  } else if (family == "binary"){
    preds.test <- predict(model, data.frame(X.test), type = "response") # only for class = 1!
    
    if (!is.null(X.cal) & !is.null(Y.cal)){
      preds.cal <- predict(model, data.frame(X.cal), type = "response")
      return(list(cal = preds.cal, test = preds.test, model = model))
    }
  }
  
  # these results are only returned if we do not have a calibration set passed
  res <- list(test = preds.test, model= model)
  return(res)
}

boosting <- function(X.fit, Y.fit, X.test, X.cal = NULL, Y.cal = NULL, n.trees = 100, family = "regression", pretrained.model = NULL, ...){
  
  if (class(X.fit)[1] != "data.frame"){
    X.fit <- as.data.frame(X.fit)
    X.test <- as.data.frame(X.test)
    names(X.test) <- names(X.fit)
    if (!is.null(X.cal)){
      X.cal <- as.data.frame(X.cal)
      names(X.cal) <- names(X.fit)
    }
    
  }
  
  # Training
  if (is.null(pretrained.model)){
    
    if (family == "binary"){
      if (length(unique(Y.fit))>2){warning("Classification is attempted with more than 2 values for target! boost")}
      if (is.factor(Y.fit)){
        Y.fit <- as.numeric(Y.fit)-1 # factor to numeric is cast as 1s and 2s
      }
      distribution <- "bernoulli"
    } else if (family == "regression"){
      distribution <- "gaussian"
    }
    
    model <- gbm::gbm(Y ~ . , data.frame(Y = Y.fit, X.fit), distribution = distribution, n.trees = n.trees)
    
  } else{
    model <- pretrained.model
  }
  preds.test <- predict(model, data.frame(X.test), type = "response", n.trees = n.trees)
  
  if (!is.null(X.cal) & !is.null(Y.cal)){
    preds.cal <- predict(model, data.frame(X.cal), type = "response", n.trees = n.trees)
    
    return(list(cal = preds.cal, test = preds.test, model = model))
  }
  
  # these results are only returned if we do not have a calibration set passed
  res <- list(test = preds.test, model= model)
  return(res)
  
}

rf <- function(X.fit, Y.fit, X.test, X.cal = NULL, Y.cal = NULL, n.trees = 500, family = "regression", pretrained.model = NULL, ...){
  # this is a function for mean rf
  if (class(X.fit)[1]!=class(X.test)[1] | !identical(names(X.fit), names(X.test))){
    names(X.test) <- names(X.fit)
  }
  
  # Training
  if (is.null(pretrained.model)){
    if (family == "regression"){
      model <- randomForest::randomForest(X.fit, Y.fit, ntree = n.trees, ...)    
    } else if (family == "binary"){
      if (length(unique(Y.fit))> 2){warning("Classification is attempted with more than 2 values for target! rf")}
      if (!is.factor(Y.fit)){
        Y.fit <- as.factor(Y.fit)
      }
      model <- randomForest::randomForest(X.fit, Y.fit, ntree = n.trees, ...)
    }
  } else{
    model <- pretrained.model
  }
  
  # Prediction
  if (family == "regression"){
    preds.test <- predict(model, X.test)
    
    if (!is.null(X.cal) & !is.null(Y.cal)){
      preds.cal <- predict(model, X.cal)
      return(list(cal = preds.cal, test = preds.test, model = model))
    }
  } else if (family == "binary"){
    preds.test <- predict(model, X.test, type = "prob")[,2] # only for class = 1!
    
    if (!is.null(X.cal) & !is.null(Y.cal)){
      preds.cal <- predict(model, X.cal, type = "prob")[,2]
      return(list(cal = preds.cal, test = preds.test, model = model))
    }
  }
  
  # these results are only returned if we do not have a calibration set passed
  res <- list(test = preds.test, model= model)
  return(res)
}

bart <- function(X.fit, Y.fit, X.test, X.cal = NULL, Y.cal = NULL, n.trees = 500, family = "regression", pretrained.model = NULL, ...){
  
  if (class(X.fit)[1] != "data.frame"){
    X.fit <- as.data.frame(X.fit)
    X.test <- as.data.frame(X.test)
    names(X.test) <- names(X.fit)
    if (!is.null(X.cal)){
      X.cal <- as.data.frame(X.cal)
      names(X.cal) <- names(X.fit)
    }
  }
  
  if (is.null(pretrained.model)){
    model <- bartMachine::bartMachine(data.frame(X.fit), Y.fit, verbose = FALSE)
  } else{
    model <- pretrained.model
  }
  
  preds.test <- predict(model, data.frame(X.test))
  
  if (!is.null(X.cal) & !is.null(Y.cal)){
    preds.cal <- predict(model, data.frame(X.cal))
    
    return(list(cal = preds.cal, test = preds.test, model = model))
  }
  return(list(test = preds.test, model = model))
}

bart_het <- function(X.fit, Y.fit, X.test, X.cal = NULL, Y.cal = NULL, n.trees = 500, family = "regression", pretrained.model = NULL, ...){
  
  if (is.null(pretrained.model)){
    sink("rbart_output.txt")
    model <- rbart::rbart(X.fit, Y.fit, ntree = 50, ntreeh = 20)
    sink()
  } else {
    model <- pretrained.model
  }
  
  preds.test <- predict(model, X.test)
  
  if (!is.null(X.cal) & !is.null(Y.cal)){
    preds.cal <- predict(model, X.cal)
    
    return(list(cal = preds.cal, test = preds.test, model = model))
  }
  
  return(list(test = preds.test, model = model))
}

quantile_boosting <- function(X.fit, Y.fit, X.test, quantiles, X.cal = NULL, Y.cal = NULL, n.trees = 100, pretrained.model = NULL, ...){
 
  if (class(X.fit)[1] != "data.frame"){
    X.fit <- as.data.frame(X.fit)
    X.test <- as.data.frame(X.test)
    names(X.test) <- names(X.fit)
    if (!is.null(X.cal)){
      X.cal <- as.data.frame(X.cal)
      names(X.cal) <- names(X.fit)
    }
  }
  
  if (is.null(pretrained.model)){
    model <- lapply(quantiles, function(x) {
      return(gbm::gbm(Y ~ ., data.frame(Y = Y.fit, X.fit), distribution = list(name = "quantile", alpha = x), n.trees = n.trees))
    })
  } else{
    model <- pretrained.model
  }
  
  preds.test <- do.call(cbind, lapply(model, function(x) predict(x, data.frame(X.test), type = "response", n.trees = n.trees)))
  
  # When a calibration set is provided also provide info for those!
  if (!is.null(X.cal) & !is.null(Y.cal)){
    preds.cal <- do.call(cbind, lapply(model, function(x) predict(x, data.frame(X.cal), type = "response", n.trees = n.trees)))
    res <- list(cal = preds.cal, test = preds.test, model = model)
    return(res) 
  }
  
  res <- list(test = preds.test, model= model)
  return(res)
}

quantile_rf <- function(X.fit, Y.fit, X.test, quantiles, X.cal = NULL, Y.cal = NULL, n.trees = 500, pretrained.model = NULL, ...){

  if (is.null(pretrained.model)){
    model <- grf::quantile_forest(X.fit, Y.fit, num.trees = n.trees, quantiles = quantiles)  
  } else{
    model <- pretrained.model
  }

  preds.test <- predict(model, data.frame(X.test))$predictions
  
  # When a calibration set is provided also provide info for those!
  if (!is.null(X.cal) & !is.null(Y.cal)){
    preds.cal <- predict(model, data.frame(X.cal))$predictions 
    res <- list(cal = preds.cal, test = preds.test, model = model)
    return(res) 
  }
  
  res <- list(test = preds.test, model= model)
  return(res)
}

quantile_bart <- function(X.fit, Y.fit, X.test, quantiles, X.cal = NULL, Y.cal = NULL, n.trees = 50, pretrained.model = NULL, ..){
  
  if (class(X.fit)[1] != "data.frame"){
    X.fit <- as.data.frame(X.fit)
    X.test <- as.data.frame(X.test)
    names(X.test) <- names(X.fit)
    if (!is.null(X.cal)){
      X.cal <- as.data.frame(X.cal)
      names(X.cal) <- names(X.fit)
    }
  }
  
  if (is.null(pretrained.model)){
    model <- bartMachine::bartMachine(data.frame(X.fit), Y.fit, num_trees = n.trees, verbose = FALSE)
  } else{
    model <- pretrained.model
  }
  
  # Note: this implementation only works for CQR until now. When for DCP, we would need all quantiles on a grid
  preds.test <- bartMachine::calc_prediction_intervals(model, data.frame(X.test), pi_conf = 1-2*min(quantiles))$interval
  
  # When a calibration set is provided also provide info for those!
  if (!is.null(X.cal) & !is.null(Y.cal)){
    preds.cal <- bartMachine::calc_prediction_intervals(model, data.frame(X.cal), pi_conf = 1-2*min(quantiles))$interval
    res <- list(cal = preds.cal, test = preds.test, model = model)
    return(res) 
  }
  
  res <- list(test = preds.test, model= model)
  return(res)
}

### CONFORMAL/DRP FUNCTIONS ###

conformal.split <- function(X.train, Y.train, X.test, w = NULL, split.prop = 0.5, split = NULL, alpha = 0.1, alpha.fit = 0.1, method = "CQR", y_model = "RF", pretrained.model = NULL, regression.parameters = NULL){

  # First split the data in training and calibration folds. If split is provided, use that one.
  if (!is.null(pretrained.model) & is.null(split)){
    warning("When using pretrained models, provide the split corresponding to the train split (part of X.train used for training the model)!")
  }
  
  splits <- create.split(X.train, Y.train, split.prop = split.prop, split = split)
  
  X.fit <- splits$X.fit; X.cal <- splits$X.cal
  Y.fit <- splits$Y.fit; Y.cal <- splits$Y.cal
  
  train.ix <- splits$train.ix
  
  # For some conformity scores, define some extra parameters needed for the base learners.
  interval.params = NULL
  if (method == "CQR"){
    # provide the quantiles to fit the base learner. Default is coverage target alpha/2. Note that deviating from default might deteriorate adaptivity!
    regression.parameters = c(regression.parameters, list(quantiles = c(alpha.fit/2, 1-alpha.fit/2)))
  } else if (method == "DCP"){
    # Fit quantile regressions on fine grid. If more precision is needed, change length of this grid.
    regression.parameters = c(regression.parameters, list(quantiles = seq(0.001, 0.999, length = 200)))
    interval.params = c(interval.params, list(dcp_grid = regression.parameters$quantiles))
  }
  
  # Obtain the correct base learner, based on the conformity score and the type of desired base learner
  regression.function <- get_regression_results(method, y_model) 
  
  # Fit the base learner to the proper training data and make predictions for the calibration and test data.
  results.regression <- obtain_model_results(regression.function, X.fit, Y.fit, X.test, c(regression.parameters, 
                                                               list(Y.cal = Y.cal, X.cal = X.cal, 
                                                                    pretrained.model = pretrained.model)))
  # Compute conformity scores for calibration data
  conformal.scores.cal <- compute_conformal_scores(results.regression$cal, Y.cal, method)
  
  # Then compute the intervals. If weights are provided, perform weighted conformal prediction
  if (is.null(w)){
    w <- rep(1, nrow(X.train) + nrow(X.test))
  }
  
  # Each test point has a custom weighted quantile when weighted conformal prediction is used!
  weighted.quantiles <- sapply(w[( nrow(X.train) + 1 ):length(w)], function (x) {
    weighted.quantile(c(conformal.scores.cal, Inf),
                      1-alpha,
                      w = c(w[1:nrow(X.train)][-train.ix], x),#c(w[-train.ix], w[nrow(X.train)+i]),
                      sorted = FALSE)
  })
  
  if (method == "DCP"){
    interval.params <- c(interval.params, list(quantiles.test = results.regression$test))
  } else if (method == "CB"){
    grid.length = 5000 # Grid in accordance with Tibshirani et al. (2019)
    grid = seq(-1.25*max(abs(Y.train)), 1.25*max(abs(Y.train)), length = grid.length)
    interval.params <- c(interval.params, list(samples.test = results.regression$test,
                                               grid = grid))
  }
  
  intervals <- compute.intervals(weighted.quantiles, results.regression$test, method, interval.params)
  
  return(list(intervals = intervals, split = train.ix, model = results.regression$model))
}

naive.QR <- function(X.train, Y.train, X.test, split.prop = 0.5, w = NULL, split = NULL, alpha = 0.1, alpha.fit = 0.1, y_model = "RF", pretrained.model = NULL, regression.parameters = NULL){
  # This method performs naive calibration of quantile regressions by choosing the nominal value of 
  # quantiles fitted to achieve marginal coverage on the calibration data
  
  # First split the data
  if (!is.null(pretrained.model) & is.null(split)){
    warning("When using pretrained models, provide the split corresponding to the train split (part of X.train used for training the model)!")
  }
  
  splits <- create.split(X.train, Y.train, split.prop = split.prop, split = split)
  
  X.fit <- splits$X.fit; X.cal <- splits$X.cal
  Y.fit <- splits$Y.fit; Y.cal <- splits$Y.cal
  
  train.ix <- splits$train.ix 
  
  # Provide the grid of quantiles that we want to fit
  alpha.grid <- seq(0.01, 0.99, length = 200)
  
  regression.parameters = c(regression.parameters, list(quantiles = alpha.grid))
  
  # Obtain the base learner (we want to do quantile regressions, so we pick CQR)
  regression.function <- get_regression_results("CQR", y_model) # This only returns the correct regression function...
  
  # Fit the base regressor and make predictions for the calibration data (for all fitted conditional quantiles)
  results.regression <- obtain_model_results(regression.function, X.fit, Y.fit, X.test, c(regression.parameters, 
                                                                                          list(Y.cal = Y.cal, X.cal = X.cal, 
                                                                                               pretrained.model = pretrained.model)))
  # We split the predictions in two parts. All alphas up to 0.5 and all alphas from 0.5-1
  preds.cal.1 <- results.regression$cal[,1:(length(alpha.grid)/2)]
  preds.cal.2 <- results.regression$cal[,(length(alpha.grid)/2+1):length(alpha.grid)]
  
  # Define function to compute the empirical coverage of a number of points
  compute_coverage <- function(preds.lower, preds.upper, y){
    coverages <- rep(NA, ncol(preds.lower))
    for (i in 1:ncol(preds.lower)){
      coverages[i] <- mean(preds.lower[, i ] <= y & y <= preds.upper[, ncol(preds.upper) - i + 1])
    }
    return(coverages)
  }
  
  # Determine the empirical coverage of the calibration points for each interval [q_alpha/2, q_1-alpha/2]
  covs <- compute_coverage(preds.cal.1, preds.cal.2, Y.cal)
  
  # Determine the value of alpha that provides intervals with empirical coverage closest to the target
  ind.opt <- which.min(abs(covs-(1-alpha)))
  alpha_opt <- 2 * alpha.grid[ ind.opt ]
  
  # Use that optimal alpha to create prediction intervals
  preds.alpha.opt <- results.regression$test[, c( ind.opt, length(alpha.grid) - ind.opt + 1)]
  
  return(list(intervals = preds.alpha.opt, model = results.regression$model, split = train.ix, alpha.opt = alpha_opt))
}

DRP <- function(X.train, Y.train, T.train, X.test, w = NULL, split.prop = 0.5, alpha = 0.1, alpha.fit = 0.1, split = NULL, pretrained.model = NULL, method = "boosting", score = "CQR", method_m = "binary", model_m = "boosting", method_pi = "RF", regression.parameters = NULL) {
  # This implements the 3-split approach of doubly robust prediction based on Yang et al (2022)
  
  # First split the data in train/calibration. We define P as the training distribution and Q as the shifted distribution (with missing Y values!)
  # We create a proper training and calibration set for each category (shifted and non-shifted data)
  
  
  X.train.P <- X.train[T.train == 0, ]; Y.train.P <- Y.train[T.train==0]
  X.train.Q <- X.train[T.train == 1, ]; Y.train.Q <- Y.train[T.train==1]
  
  splits.P <- create.split(X.train.P, Y.train.P, split.prop = split.prop, split = split)
  splits.Q <- create.split(X.train.Q, Y.train.Q, split.prop = split.prop, split = split)
  
  train.ix.P <- splits.P$train.ix
  train.ix.Q <- splits.Q$train.ix
  
  drp.train.prop = 0.5
  
  X.P.fit <- splits.P$X.fit; X.P.cal <- splits.P$X.cal
  Y.P.fit <- splits.P$Y.fit; Y.P.cal <- splits.P$Y.cal
  
  X.Q.fit <- splits.Q$X.fit; X.Q.cal <- splits.Q$X.cal
  Y.Q.fit <- splits.Q$Y.fit; Y.Q.cal <- splits.Q$Y.cal
  
  # Then split the training data also in two parts (this is the three split DRP. First split is used for computing the regression model, the 
  # second split is used to train the pi and m functions)
  ind.fit.P.1 <- sample(1:nrow(X.P.fit), floor(drp.train.prop * nrow(X.P.fit)))
  ind.fit.Q.1 <- sample(1:nrow(X.Q.fit), floor(drp.train.prop * nrow(X.Q.fit)))
  
  X.P.fit.1 <- X.P.fit[ind.fit.P.1,]; X.P.fit.2 <- X.P.fit[-ind.fit.P.1,]
  Y.P.fit.1 <- Y.P.fit[ind.fit.P.1]; Y.P.fit.2 <- Y.P.fit[-ind.fit.P.1]
  
  X.Q.fit.1 <- X.Q.fit[ind.fit.Q.1,]; X.Q.fit.2 <- X.Q.fit[-ind.fit.Q.1,]
  Y.Q.fit.1 <- Y.Q.fit[ind.fit.Q.1]; Y.Q.fit.2 <- Y.Q.fit[-ind.fit.Q.1]

  
  # Set parameters of different base learners (if necessary)
  interval.params = NULL
  if (score == "CQR"){
    regression.parameters = c(regression.parameters, list(quantiles = c(alpha.fit/2, 1-alpha.fit/2)))
  }
  
  # Create the appropiate regression model/base learner
  regression.function <- get_regression_results(score, method)
  
  # Run the regression function on the correct training sample and provide predictions for the calibration set(s)
  results.regression <- obtain_model_results(regression.function, X.P.fit.1, Y.P.fit.1, X.test, c(regression.parameters, 
                                                                                          list(Y.cal = c(Y.P.fit.2 , Y.P.cal), 
                                                                                               X.cal = rbind(X.P.fit.2, X.P.cal), 
                                                                                               pretrained.model = pretrained.model
                                                                                              )))
  
  # Compute the conformity scores. Note that .cal refers to the true calibration data and R.2 refers to the data used to train the m and pi functions
  conformal.scores.cal <- compute_conformal_scores(as.matrix(results.regression$cal)[-(1:nrow(X.P.fit.2)),], Y.P.cal, score)
  conformal.scores.R.2 <- compute_conformal_scores(as.matrix(results.regression$cal)[(1:nrow(X.P.fit.2)), ], Y.P.fit.2, score)

  # If no weights are provided, estimate the pi function using a classification model.
  if (is.null(w)){
    model.pi <- get_classification_results(method_pi) # note since this is a classification method, we provide "mean" as parameter
    classification.results <- obtain_model_results(model.pi, rbind(X.P.fit.2, X.Q.fit.2), c(rep(0, nrow(X.P.fit.2)), rep(1, nrow(X.Q.fit.2))), 
                                                   rbind(X.P.cal, X.Q.cal), list(
                                                                pretrained.model = pretrained.model,
                                                                family = "binary"))
    # Cut the probabilities to prevent weird behaviour 
    pi.hat <- pmax(pmin(classification.results$test, 0.99), 0.01)
    pi.hat <- pi.hat/(1-pi.hat) 
    
  } else {
    # Note: then the w should be a function and not vector with values!!!
    if (!is.function(w)){stop("Note that if w is provided for DRP, it should be a function!")}
    pi.hat <- w(rbind(X.P.cal, X.Q.cal))
    pi.hat <- pi.hat/(1-pi.hat)
  }
  
  # Two different ways of estimating the conditional distribution of the conformity scores. Binary is used for the results in the thesis/paper
  if (method_m == "binary"){
    model.m <- function(theta, Xnew, Xtrain, Ytrain){ # partially adopted from Yang et al. (2022)
      # First get the necessary model
      model <- get_classification_results(model_m)
      
      # Run the regression model and get predictions for Xnew
      model.results <- obtain_model_results(model, Xtrain, Ytrain, Xnew, list(family = "binary"))
      
      preds <- pmax(pmin(model.results$test, 0.99), 0.01)
      return(preds)
    }
  } else if (method_m == "quantile"){
    # In this approach, we directly estimate quantile functions and (discretely) approximate the value of the CDF with it. We need a grid of quantiles.
    model <- get_regression_results( method = "CQR", y_model = model_m)
    
    model.results <- obtain_model_results(model, X.P.fit.2, conformal.scores.R.2, 
                                          rbind(X.P.cal, X.Q.cal), 
                                          list(quantiles = seq(0.01, 0.99, length = 200)))
    
    preds.m <- model.results$test
  }
  
  # Define the efficient influence function (see paper for definition)
  EIF <- function(theta){
    if (method_m == "binary"){
      m.hat <- model.m(theta, rbind(X.P.cal, X.Q.cal), X.P.fit.2, (as.numeric(conformal.scores.R.2 <= theta)))
      1/(nrow(X.P.cal) + nrow(X.Q.cal)) * (sum( pi.hat[1:nrow(X.P.cal)] * ( I(conformal.scores.cal <= theta) - m.hat[1:nrow(X.P.cal)] ) ) + 
                                             sum( m.hat[ (1+nrow(X.P.cal)) : length(m.hat) ]  - (1-alpha) ))
    } else if (method_m == "quantile"){
      1/(nrow(X.P.cal) + nrow(X.Q.cal)) * (sum( pi.hat[1:nrow(X.P.cal)] * ( I(conformal.scores.cal <= theta) - rowMeans(preds.m[1:nrow(X.P.cal),] <= theta) ) ) + 
                                             sum( rowMeans(preds.m[ (1+ nrow(X.P.cal)) : (nrow(X.P.cal) + nrow(X.Q.cal)),] <= theta)  - (1-alpha) ))
    }
  }
  
  # Find the smallest value such that the empirical EIF = 0. Split for different conformity scores, to specify different search ranges.
  if (score == "mean"){ 
    r.hat <- try(uniroot(EIF, c(1, max(conformal.scores.cal, conformal.scores.R.2)/2), extendInt= "yes", tol = 0.1)$root)
    # If root cannot be found, change the search parameters (range)
    if(class(r.hat)!= "try-error"){
      r.hat=r.hat
    } else {
      r.bound.1 <- try(EIF(max(conformal.scores.cal, conformal.scores.R.2)/2))
      r.bound.2 <- try(EIF(max(conformal.scores.cal, conformal.scores.R.2)/1.5))
      if (class(r.bound.1)!= "try-error" & r.bound.1 > 0 ){
        print("Enter else if 1: ")
        r.hat = try(uniroot(EIF,c(min(conformal.scores.cal, conformal.scores.R.2)/2,max(conformal.scores.cal, conformal.scores.R.2)/2),tol=0.1)$root)
      } else if (class(r.bound.2)!= "try-error" & r.bound.2 > 0 ){
        print("Enter else if 2: ")
        r.hat = try(uniroot(EIF,c(max(conformal.scores.cal, conformal.scores.R.2)/2,max(conformal.scores.cal, conformal.scores.R.2)+0.1),tol=0.1, extendInt = "yes")$root)
      } else {
        print("enter here")
        r.hat = try(uniroot(EIF,c(max(conformal.scores.cal, conformal.scores.R.2)/1.5,max(conformal.scores.cal, conformal.scores.R.2)/1.25),tol=0.1, extendInt = "yes")$root)
        
        
      }
      if (class(r.hat) == "try-error"){
        r.hat <- Inf
      }
    }
  } else if (score == "CQR"){
    if(max(conformal.scores.cal) <= 0){
      print("smaller than 0")
    }
    r.hat <- try(uniroot(EIF, c(min(conformal.scores.cal, conformal.scores.R.2)/2, max(conformal.scores.cal, conformal.scores.R.2)/4), extendInt= "upX", tol = 0.01)$root)
    if(class(r.hat)!= "try-error"){
      r.hat=r.hat
    }
    else {
      r.bound.1 <- try(EIF(1))
      r.bound.2 <- try(EIF(max(conformal.scores.cal, conformal.scores.R.2)/2))
      if (class(r.bound.1)!= "try-error" & r.bound.1 > 0 ){ # There can be an error here!
        print("Enter else if 1: ")
        r.hat = try(uniroot(EIF,c(min(conformal.scores.cal)/2, 1),tol=0.01, extendInt = "yes")$root)
      } else{
        print("Enter else if 2: ")
        r.hat = try(uniroot(EIF,c(1,max(conformal.scores.cal, conformal.scores.R.2)/2),tol=0.01, extendInt = "upX")$root)
      }
    }
    if (class(r.hat) == "try-error"){
      r.hat <- Inf
    }
  }
  
  # Compute intervals
  intervals.drp <- compute.intervals(r.hat, results.regression$test, score, interval.params = NULL)
  
  return(list(intervals = intervals.drp, predictions = results.regression$test, model = results.regression$model, split = list(P = train.ix.P, Q = train.ix.Q)))
  
}

### Help functions for conformal prediction ###

# Return the correct base learner for each conformal predictor
get_regression_results <- function(method = "CQR", y_model = "RF"){
  if (method=="CQR"){
    if (y_model == "RF"){
      return(quantile_rf)
    } else if (y_model == "boosting"){
      return(quantile_boosting)
    } else if (y_model == "bart"){
      return(quantile_bart)
    }
  }
  if (method == "mean"){
    if (y_model == "RF"){
      return(rf)
    } else if (y_model == "bart"){
      return(bart)
    } else if (y_model == "boosting"){
      return(boosting)
    } 
  }
  if (method == "DCP"){
    if (y_model == "RF"){
      return(quantile_rf)
    } else if (y_model == "boosting"){
      return(quantile_boosting)
    } else if (y_model == "bart"){
      warning("BART cannot be used for DCP yet, since the quantiles need to be on a grid, and only the pi_conf can be provided...")
    }
  }
  if (method == "CB"){
    if (y_model == "bart_het"){
      return(bart_het)   
    } else if (y_model == "bart"){
      warning("CB cannot be used for BART yet, since the output format of the things is different! See if another level of abstraction is possible...")
    }
    
  }
  
  stop(paste0("Currently no other methods are supported for ", method))
}

# Return classifier for weight estimation
get_classification_results <- function(method){
  if (method == "RF"){
    return(rf)
  } else if (method == "boosting"){
    return(boosting)
  } else if (method == "glm"){
    return(linear)
  }
}

compute_conformal_scores <- function(Y.pred, Y, method){
  
  # See paper for the definiiton of the different conformity scores
  if (method == "CQR"){
    
    conformal.scores <- pmax(Y.pred[,1] - Y,  Y - Y.pred[,2] )
    
  } else if (method == "mean"){
    
    conformal.scores <- abs(Y - Y.pred)
    
  } else if (method == "DCP"){
    
    # Do DCP conformal score here
    predicted.quantiles <- Y.pred
    
    quantiles.sorted <- t(apply(predicted.quantiles, 1, FUN = sort)) # this is still before going through the logit function!
    U.hat <- rowMeans((quantiles.sorted <= matrix(Y, length(Y), ncol(quantiles.sorted))))
    
    conformal.scores <- abs(U.hat-0.5)
    
    
    
  } else if (method == "naive"){
    
    
  } else if (method == "CB"){
    
    conformal.scores <- rep(NA, length(Y))
    for ( j in 1:length(Y)){
      conformal.scores[j] <- -dnorm(Y[j], 
                                    mean = Y.pred$mmean[j], sd = Y.pred$smean[j]) # cannot be used for bart yet...
    }
    
  }
  
  return(conformal.scores)
}

# Compute the intervals for the conformal approach, using the weighted quantile as input
compute.intervals <- function(weighted.quantiles, Y.pred, method, interval.params){
  # See paper for the shape of the different prediction regions.
  if (method == "CQR"){
    res <- cbind(Y.pred[,1] - weighted.quantiles, Y.pred[,2] + weighted.quantiles)
  } else if (method == "mean"){
    res <- cbind(Y.pred - weighted.quantiles, Y.pred + weighted.quantiles)
  } else if (method == "DCP"){
    
    ci.grid <- abs(interval.params$dcp_grid-0.5)
    
    quantiles.test.sorted <- t(apply(interval.params$quantiles.test, 1, FUN = sort))
    
    intervals.dcp.qr <- matrix(NA, nrow = nrow(quantiles.test.sorted), ncol = 2)
    
    for (i in 1:nrow(quantiles.test.sorted)){
      ci <- quantiles.test.sorted[i,(ci.grid <= weighted.quantiles[i])] # See notes with paper!
      intervals.dcp.qr[i,1] <- min(ci)
      intervals.dcp.qr[i,2] <- max(ci)
    }
    res <- intervals.dcp.qr
  } else if (method == "naive"){
    res <- compute.intervals.naive.QR()
  } else if (method == "CB" ){
    intervals.bayes <- matrix(NA, nrow = length(weighted.quantiles), ncol = 2)
    conformal.scores.test <- matrix(NA, length(weighted.quantiles), ncol = length(interval.params$grid))

    for (j in 1:length(weighted.quantiles)){
      conformal.scores.test[j,] <- dnorm(interval.params$grid,  
                                         mean = interval.params$samples.test$mmean[j], # For the fast version, this is only one number!
                                         sd = interval.params$samples.test$smean[j])
    
      intervals.bayes[j,1] <- min(interval.params$grid[which(conformal.scores.test[j,] > -weighted.quantiles[j])])
      intervals.bayes[j,2] <- max(interval.params$grid[which(conformal.scores.test[j,] > -weighted.quantiles[j])])
    }
    res <- intervals.bayes
  }
  return(res)
}

# Run the regression model with the correct training and testing data
obtain_model_results <- function(reg.func, X.fit, Y.fit, X.test, reg.parameters){
  do.call(reg.func, c(reg.parameters, list(Y.fit = Y.fit, X.fit = X.fit, 
                                           X.test = X.test)
  ))
}

# This function is adopted from Tibshirani (2019)
weighted.quantile <- function(v, prob, w=NULL, sorted=FALSE) {
  if (is.null(w)) w = rep(1,length(v))
  if (!sorted) { o = order(v); v = v[o]; w = w[o] }
  i = which(cumsum(w/sum(w)) >= prob) # take cumsum, because it represents the CDF ("we sum all pdf values")
  if (length(i)==0) return(Inf) # Can happen with infinite weights
  else return(v[min(i)])
}

# Split the data in calibration and proper training
create.split <- function(X.train, Y.train, split.prop, split = NULL){
  
  # If a split is provided, use those indices for training
  if (!is.null(split)){
    train.ix <- split
  }
  else{
    # First split the data in train/calibration
    train.ix <- sample(1:nrow(X.train), floor(split.prop * nrow(X.train)))
  }
  
  X.fit <- as.matrix(X.train[train.ix, ], nrow = length(train.ix))
  X.cal <- as.matrix(X.train[-train.ix, ], nrow = nrow(X.train) - length(train.ix))
  
  Y.fit <- Y.train[train.ix]
  Y.cal <- Y.train[-train.ix]
  
  return(list(X.fit = X.fit, X.cal = X.cal, Y.fit = Y.fit, Y.cal = Y.cal, train.ix = train.ix))
  
}

# Function to create bins and compute the coverage (or length) metric across bins. This function uses equal numbers of obs in each bin
binning <- function(X,res.mat,num.seg){
  cond      <- matrix(NA,num.seg,dim(res.mat)[2])
  xs        <- matrix(NA,num.seg,1)
  quantiles <- seq(0,1,length=num.seg+1)
  for (i in 2:(num.seg+1)){
    q1                <- quantile(X,quantiles[i])
    q0                <- quantile(X,quantiles[i-1])
    ind               <- which((X<=q1)*(X>q0)==1)
    cond[(i-1),]      <- colMeans(res.mat[ind,],na.rm=TRUE)
    xs[(i-1),]        <- (q1+q0)/2 # midpoint
  }
  return(list(xs=xs,cond=cond))
}

# Note that this function only works for 1D X!
binning.equispaced <- function(X, res.mat, num.seg){
  cond <- matrix(NA, num.seg, dim(res.mat)[2])
  xs <- matrix(NA, num.seg, 1)
  boundaries <- seq(min(X), max(X), length = num.seg+1)
  
  for (i in 2:(num.seg+1)){
    ind <- (boundaries[i-1] <= X & X <= boundaries[i])
    
    cond[(i-1), ] <- colMeans(as.matrix(res.mat[ind,]), na.rm = T)
    xs[(i-1), ] <- (boundaries[i]+boundaries[i-1])/2 # midpoint
  }
  
  return(list(xs = xs, cond = cond))
}

# Generate 1D Data for the results in Chapter 4
generate.data.1D <- function(n, homoscedastic = T, tnoise = F, shift = FALSE, c_sd = 0.5){
  
  if (shift){
    xs <- matrix(runif(n, -3, 5), n)
    ps <- 0.25+exp(xs-2)/(1+exp(xs-2))*1/2
    x <- c()
    while( length(x) <= n){
      j <- which(runif(n) <= ps/sum(ps))[1]
      if (!is.na(j)){
        x <- c(x, xs[j])
      }
    }
    x <- matrix(x[1:n], nrow = n)
  }
  else {
    x <- matrix(runif(n, -3, 5), nrow = n)
  }
  
  y.noiseless <- as.numeric(x^2 - 0.25 * x^3)
  if (homoscedastic){
    if (!tnoise){
      y <- as.numeric(y.noiseless + rnorm(n))
      y.up <- y.noiseless + qnorm(alpha/2)
      y.low <- y.noiseless - qnorm(alpha/2)
    } else{
      y <- as.numeric(y.noiseless + rt(n, 2))
      y.up <- y.noiseless + qt(alpha/2, 2)
      y.low <- y.noiseless - qt(alpha/2, 2)
    }
  } else {
    y <- as.numeric(y.noiseless + rnorm(n, sd = exp(c_sd * x)))
    y.up <- as.numeric(y.noiseless + qnorm(alpha/2, sd = exp(c_sd * x)))
    y.low <- as.numeric(y.noiseless - qnorm(alpha/2, sd = exp(c_sd * x)))
  }
  
  
  return(list("x" = x, "y" = y, "y_noiseless" = y.noiseless, "y_up" = y.up, "y_lo" = y.low))
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

# Implement prediction intervals for the HBART algorithm. The ideas are adopted form the bartMachine implementation
calc_hbart_intervals <- function(hbart, x.test, quantiles){
  
  # Use the fitted hbart model to make predictions
  preds.hbart <- predict(hbart, x.test)
  
  # How many points to sample from the posterior draws
  num_samples_per_data_point = 1000
  n_test = nrow(x.test)
  all_prediction_samples = matrix(NA, nrow = n_test, ncol = num_samples_per_data_point)
  q_mat = matrix(NA, nrow = n_test, ncol = length(quantiles))
  for (i in 1:n_test){
    
    # Retrieve the posterior draws for the mean and variance model
    y_hats = preds.hbart$mdraws[, i]
    sigma_hats = preds.hbart$sdraws[, i]
    
    # From the posterio samples, sample 1000 poitns
    n_gs = sample(1:hbart$ndpost, num_samples_per_data_point, replace = TRUE)
    
    for (j in 1:num_samples_per_data_point){
      
      y_hat_draw = y_hats[n_gs[j]]
      sigma_hat_draw = sigma_hats[n_gs[j]]
      
      # Sample from a normal distribution with the draws from the posterior
      all_prediction_samples[i, j] = rnorm(1, mean = y_hat_draw, sd = (sigma_hat_draw))
    }
    
    # Take the relevant quantiles
    q_mat[i, ] = quantile(all_prediction_samples[i, ], quantiles)
  }
  
  return(list(quantiles = q_mat, preds = preds.hbart))
  
}
