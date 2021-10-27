##Aditya Prabaswara Mardjikoen (S2264710)
linmod <- function(formula,dat){
  X <- model.matrix(~group,dat) ##model matrix
  qrx <- qr(X) ##get QR decomposition of model matrix
  y <- model.frame(formula,dat)[[1]] ##response
  yname <- all.vars(formula)[1] ##response variable name
  p <- ncol(X) ##number of parameters
  V <- cov(dat) ##covariance matrix
  beta <- backsolve(qr.R(qrx),qr.qty(qrx,y)[1:p]) ##get beta
}

print.linmod <- function(x){
  
}

plot.linmod <- function(x){
  
}

predict.linmod <- function(x,newdata){
  
}