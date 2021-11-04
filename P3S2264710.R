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
  
  
  X <- model.matrix(formula,dat) ##model matrix
  qrx <- qr(X) ##get QR decomposition of model matrix
  
  
  y <- model.frame(formula,dat)[[1]] ##response variable data
  yname <- all.vars(formula)[1] ##response variable name
  p <- ncol(X) ##number of parameters
  n <- nrow(X) ##number of observations
  
  
  beta <- backsolve(qr.R(qrx),qr.qty(qrx,y)[1:p]) ##regression parameter
  names(beta) <- colnames(X) ##label beta column based on the column name of X
  
  
  mu <- X %*% beta ##fitted values or predicted values
  mu <- drop(mu) ##return fitted value result as a vector
  sse <- sum((mu-y)^2) ##calculate sum squared error
  sigma <- sqrt(sse/(n-p)) ##estimated standard deviation
  
  
  ##Recall that in QR decomposition R is an upper triangular matrix. Therefore, R^{T} (transpose of R) is an lower triangular matrix. Thus, we can calculate
  ##R^{-T} by solving R^{T}x=I using forwardsolve and calculate R^{-1}R^{-T} by solving Rx=R^{-T} because R is an upper triangular matrix. 
  
  
  inv_Rt <- forwardsolve(t(qr.R(qrx)),diag(p)) ##calculate R^{-T}
  inv_RtR <- backsolve(qr.R(qrx),inv_Rt) ##calculate R^{-1}R^{-T}
  V <- sigma^2 * inv_RtR ##covariance matrix
  
  
  rownames(V) <- colnames(X) ##name the row of covariance matrix by the parameter name
  colnames(V) <- colnames(X) ##name the column of covariance matrix by the parameter name
  
  
  flev <- list() ##initialize factor level list
  flev_var <- c() ##initialize vector of factor variable name
  
  
  factors <- names(which(sapply(dat,is.factor))) ##finding factor variable in the data
  

  factor_index <- which(all.vars(formula) %in% factors) ##finding factor variable in the regression formula
  
  ##find level factor and the corresponding factor variable name that was used in the formula
  for (i in all.vars(formula)[factor_index]){
      flev <- c(flev,list(levels(dat[,i])))
      flev_var <- c(flev_var, i)
  } 
  
  names(flev) <- flev_var ##label level factor vector with factor variable name
  
  ##store all required calculation result in a list and return object of class linmod
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
  
  
  model_report <- cbind(x$beta,se)##combine beta and standard error into matrix by column
  colnames(model_report) <- c('Estimate','s.e.')##give name to the matrix column
  
  ##displaying the regression model formula along with the parameter and standard error of the parameter
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
  mtext(text="Residuals vs Fitted",side=3)##plot title
}


predict.linmod <- function(x,newdata){
  
  ##function to estimate the response variable value given newdata
  
  ##input: x = object of class linmod , newdata = data without response variable
  
  ##ouput:vector containing the predicted value or response variable by using newdata
  
  
  ##recall that newdata is not containing response variable. So, we have to add the dummy
  ##response variable(0 and 1) so when we specified the formula from x class model matrix
  ##can detect it in newdata and perform the correct calculation that we desired. 
  
  
  ##give dummy response variable to newdata if this data doesn't contain response variable
  if (is.null(newdata[[x$yname]])==TRUE){
    newdata[[x$yname]]<-sample(0:1,nrow(newdata),replace=TRUE)
  }
  
  ##convert all variable that supposed to be factor in the newdata if it has factor variable
  ii <- which(!all.vars(x$formula)[-1] %in% colnames(newdata))
  if (length(ii)!=0){
    output <- 'Error, Independent Variable In Formula Not Found In Input Data!'
  }
  else if(length(x$flev)==0){
    output <- 'No Error'
  }
  else{
    for (i in names(x$flev)){
      if (is.factor(newdata[[i]])==FALSE){
        newdata[[i]] <- factor(newdata[[i]])
      }
      
      jj <- which(!levels(newdata[[i]]) %in% x$flev[[i]])
      
      if (length(jj)!=0){
        output <- 'Error, Level Factor Is Difference!'
      }else{
        output <-'No Error'
      }
    }
  }
  
  ##display prediction result
  if (output != 'No Error'){
    result <- output
  }
  else{
    X_new <- model.matrix(x$formula, data=newdata, xlev=x$flev) ##model matrix for newdata
    mu_new <- X_new %*% x$beta
    mu_new <- drop(mu_new)
    result <- mu_new
  }
  return(result)
}