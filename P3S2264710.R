##Aditya Prabaswara Mardjikoen (S2264710)
linmod <- function(formula,dat){
  ##initialize the QR decomposition calculation
  X <- model.matrix(formula,dat) ##model matrix
  qrx <- qr(X) ##get QR decomposition of model matrix
  
  ##find response variable and the number of parameters
  y <- model.frame(formula,dat)[[1]] ##response
  yname <- all.vars(formula)[1] ##response variable name
  p <- ncol(X) ##number of parameters
  
  ##calculate beta and sigma
  beta <- backsolve(qr.R(qrx),qr.qty(qrx,y)[1:p]) ##get beta
  names(beta) <- colnames(X) ##label beta column based on the column name
  mu <- X %*% beta ##fitted values
  mu <- drop(mu)
  sigma_resp <-
  sigma_resid <- sqrt(sum(mu-y)^2/(nrow(dat)-p))##estimated standard deviation
  sigma <- c(sigma_resp, sigma_resid)
  names(sigma) <- c('response', 'residual')
  
  ##calculate covariance matrix
  inv_Rt <- forwardsolve(t(qr.R(qrx)),diag(ncol(t(qr.R(qrx)))))
  inv_RtR <- backsolve(qr.R(qrx),inv_Rt)
  V <- sigma^2*inv_RtR ##covariance matrix
  rownames(V) <- colnames(X)
  colnames(V) <- colnames(X)
  
  ##finding factor and its levels
  factors <- names(Filter(is.factor,dat))
  flev <- list()##initialize factor level list
  flev_var <- c()
  factor_index <- which(all.vars(formula)==factors)
  for (i in all.vars(formula)[factor_index]){
    if (i %in% factors){
      flev <- c(flev,list(levels(dat[,i])))
      flev_var <- c(flev_var, i)
    }
  }
  names(flev) <- flev_var
  
  model_summary <- list(beta=beta, V=V, mu=mu, y=y,yname=yname,formula=formula, flev=flev, sigma=sigma)
  return(model_summary)
}

print.linmod <- function(x){
  
}

plot.linmod <- function(x){
  
}

predict.linmod <- function(x,newdata){
  
}