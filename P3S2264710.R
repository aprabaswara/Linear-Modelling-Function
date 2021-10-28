##Aditya Prabaswara Mardjikoen (S2264710)
library(ggplot2)
dat <- mpg
formula <- cty~fl+displ+cyl
linmod <- function(formula,dat){
  X <- model.matrix(formula,dat) ##model matrix
  qrx <- qr(X) ##get QR decomposition of model matrix
  y <- model.frame(formula,dat)[[1]] ##response
  yname <- all.vars(formula)[1] ##response variable name
  p <- ncol(X) ##number of parameters
  V <- cov(X) ##covariance matrix
  
  beta <- backsolve(qr.R(qrx),qr.qty(qrx,y)[1:p]) ##get beta
  names(beta) <- colnames(X) ##label beta column based on the column name
  mu <- X %*% beta ##fitted values
  sigma <- sqrt(sum(mu-y)^2/(nrow(dat)-p-1))##estimated standard deviation
  
  flev <- vector(mode = "list")##initialize factor list
  flev_var <- c()
  for (i in colnames(dat)){
    if (is.factor(dat[,i])==TRUE){
      flev.append(i)
      flev_var <- c(flev_var, i)
      }
    }
  names(flev)<-flev_var
  
  model_summary <- list(beta=beta, mu=mu, y=y,yname=yname,formula=formula,flev=flev, sigma=sigma)
  return(model_summary)
}

print.linmod <- function(x){
  
}

plot.linmod <- function(x){
  
}

predict.linmod <- function(x,newdata){
  
}