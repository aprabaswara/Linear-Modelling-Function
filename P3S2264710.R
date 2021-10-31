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
  sigma <- sqrt(sum((mu-y)^2)/(nrow(dat)-p))##estimated standard deviation
  
  ##calculate covariance matrix
  inv_Rt <- forwardsolve(t(qr.R(qrx)),diag(p))
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
  beta_estimate <- x$beta
  cov_matrix <- x$V
  reg_formula <- x$formula
  se_square <- c()
  for (i in 1:nrow(cov_matrix)){
    for (j in 1:nrow(cov_matrix)){
      if (i == j){
        se_square <- c(se_square, cov_matrix[i,i])
      }
    }
  }
  se <- sapply(se_square,sqrt)
  model_report <- cbind(beta_estimate,se)
  colnames(model_report) <- c('Estimate','s.e.')
  print(reg_formula)
  cat('\n')
  print(model_report)
}

plot.linmod <- function(x){
  resid <- x$y-x$mu
  fitted <- x$mu
  plot(x=fitted,y=resid,xlab='Fitted values',ylab='Residuals')
  lines(lowess(x=fitted,y=resid),col="red")
  abline(h=0,lty=3)
  mtext(text="Residuals vs Fitted",side=3)
}

predict.linmod <- function(x,newdata){
  x <- linmod(formula=len~supp+dose,dat=ToothGrowth)
  names(x$flev)
  x$yname
  x$formula
  newdata <- data.frame(supp=c('VC','VC','OJ','OJ'),dose=c(2.8,1,0.3,0.4))
  for (i in names(newdata)){
    if (i %in% names(x$flev) & is.factor(newdata[,i])==FALSE){
      newdata[,i] <- factor(newdata[,i])
    }
  }
}