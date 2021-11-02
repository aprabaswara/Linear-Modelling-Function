##Aditya Prabaswara Mardjikoen (S2264710)


##Overview:
##------------------------------------------------------------------------------------------------------------------------
##The aim of this code is to create a function that can be use for fitting linear model by using the 
##QR decomposition. Overall,there are four task that this code will performed:
##
##1. Calculate and display the model summary after fitting the linear model by using linmod function
##
##2. Calculate and display the regression model parameter and standard error of its parameter using print.linmod function
##
##3. Display the residual vs predicted value plot using plot.linmod function
##
##4. Predict the value of response variable if a newdata is given by using predict.linmod function
##
##------------------------------------------------------------------------------------------------------------------------


linmod <- function(formula,dat){
  
  ##function to estimates the specified linear model using the QR decomposition of the model matrix approach.
  
  ##input: formula = regression model formula, dat = the data used to build regression model
  
  ##output: 
  ##list containing this following element:
  
  ##beta = vector of estimated regression parameter; V = covariance matrix;
  ##mu = vector of fitted values or predicted values of the response variable; y = vector containing the response variable data;
  ##yname = response variable name; formula = model formula;
  ##flev = a named list containing factor and level factor; sigma = the estimated standard deviation of the regression model
  
  
  ##initialize the QR decomposition calculation
  X <- model.matrix(formula,dat) ##model matrix
  qrx <- qr(X) ##get QR decomposition of model matrix
  
  
  
  ##find response variable and the number of parameters
  y <- model.frame(formula,dat)[[1]] ##response variable data
  yname <- all.vars(formula)[1] ##response variable name
  p <- ncol(X) ##number of parameters
  n <- nrow(X) ##number of observations
  
  
  
  ##calculate the estimated regression parameter
  beta <- backsolve(qr.R(qrx),qr.qty(qrx,y)[1:p]) ##regression parameter
  names(beta) <- colnames(X) ##label beta column based on the column name
  
  
  
  ##calculate fitted values and the estimated standard deviation
  mu <- X %*% beta ##fitted values
  mu <- drop(mu) ##return fitted value result as a vector
  sse <- sum((mu-y)^2)##calculate sum squared error
  sigma <- sqrt(sse/(n-p))##estimated standard deviation
  
  
  
  ##Recall that in QR decomposition R is an upper triangular matrix. Therefore, R^{T} (transpose of R) is an lower triangular matrix. Thus, we can calculate
  ##R^{-T} by solving R^{T}x=I using forwardsolve and calculate R^{-1}R^{-T} by solving Rx=R^{-T} because R is an upper triangular matrix. 
  
  
  ##calculate covariance matrix
  inv_Rt <- forwardsolve(t(qr.R(qrx)),diag(p))##calculate R^{-T}
  inv_RtR <- backsolve(qr.R(qrx),inv_Rt)##calculate R^{-1}R^{-T}
  V <- sigma^2*inv_RtR ##covariance matrix
  
  
  
  ##rename the column and row of covariance matrix by the regression parameter name
  rownames(V) <- colnames(X)
  colnames(V) <- colnames(X)
  
  
  
  ##finding factor variable in the data
  factors <- names(which(sapply(dat,is.factor)==TRUE))
  
  
  
  ##initialize factor level list and vector of factor variable name
  flev <- list()
  flev_var <- c()
  
  
  
  ##find level factor and factor variable name
  factor_index <- which(all.vars(formula)==factors)
  
  for (i in all.vars(formula)[factor_index]){
    if (i %in% factors){
      flev <- c(flev,list(levels(dat[,i])))
      flev_var <- c(flev_var, i)
    }
  }
  
  names(flev) <- flev_var ##label level factor vector with factor variable name
  
  
  
  ##store all required calculation result in a list
  model_summary <- list(beta=beta, V=V, mu=mu, y=y,yname=yname,formula=formula, flev=flev, sigma=sigma)
  class(model_summary) <- 'linmod'
  return(model_summary)
}



print.linmod <- function(x){
  
  ##function to display the model formula along with the estimated parameter and its standard error
  
  ##input: x = object of class linmod
  
  ##output: model summary containing regression formula, estimated value of each parameter with its standard error
  
  
  ##To calculate the estimated standard error of regression coefficient we can do it by calculating the square root
  ##of each entries in the diagonal of covariance matrix.
  
  
  se <- sapply(diag(x$V),sqrt) ##standard error of regression parameter
  
  
  ##displaying the regression model formula along with the parameter and standard error of the parameter
  model_report <- cbind(x$beta,se)##combine beta and standard error into matrix by column
  colnames(model_report) <- c('Estimate','s.e.')##give name to the matrix column
  print(x$formula)
  cat('\n')
  print(model_report)
}



plot.linmod <- function(x){
  
  ##function to plot residual vs fitted value
  
  ##input: x = object of class linmod
  
  ##output: residual vs fitted value plot
  
  
  resid <- x$y-x$mu ##residual
  fitted <- x$mu ##fitted value
  
  plot(x=fitted,y=resid,xlab='Fitted values',ylab='Residuals')##create plot
  lines(lowess(x=fitted,y=resid),col="red")##trend line of residual pattern
  abline(h=0,lty=3)##horizontal line indicating residual equal zero
  mtext(text="Residuals vs Fitted",side=3)##give title to the plot
}


predict.linmod <- function(x,newdata){
  
  ##function to estimate the response variable value given newdata
  
  ##input: x = object of class linmod , newdata = data without response variable
  
  ##ouput:vector containing the predicted value or response variable by using newdata
  
  x <- linmod(formula=len~supp+dose,dat=ToothGrowth)
  newdata1 <- data.frame(supp=c('VC','VC'),dose=c(2.8,1))
  if (is.null(newdata1[[x$yname]])==TRUE){
    newdata1[[x$yname]]<-sample(0:1,nrow(newdata1),replace=TRUE)
  }
  
  for (i in names(x$flev)){
    if (is.factor(newdata1[[i]])==FALSE){
      newdata1[[i]] <- factor(newdata1[[i]])
    }
  }
  
  
  if (colnames(X)==colnames(x$V)){
    X = model.matrix(x$formula,newdata1)
    mu_new = X %*% x$beta
  }
  else{
    if (length(colnames(X))>length(names(x$beta))){
      unmatch_index <- which(is.na(match(colnames(X),names(x$beta)))==TRUE)
      X = model.matrix(x$formula,newdata1)
      mu_new = X[,-unmatch_index] %*% x$beta
    }
    else{
      
    }
  }
  
}