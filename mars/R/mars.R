#' Multivariate Adaptive Regression Splines (MARS)
#'
#' Fit Friedman's Multivariate Adaptive Regression Splines (MARS) model.
#'
#' @param formula an regression formula
#' @param data a data frame containing the data
#' @param control a mars control object (see mars.control to build one)
#' @usage first create a mars control object using mars.control then
#' mars(y~.,data, marscontrol)
#' @details This function creates a Multivariate Adaptive Regression Splines
#' model for set of data with a forward step and backward step based using Friedman's
#' method
#'
#'
#' @return a mars class containing the call, input formula, y values, coefficients,
#' and data frame containing the split points, and direction for a given coefficient
#'
#' @author Jaspal Raman, Ben Shires Nakamura, Jessica Kim
#' @references Multivariate Adaptive Regression Splines
#' @seealso anova.mars, print.mars, plot.mars, predict.mars, print.mars, summary.mars
#' @example mars(y~., data, marscontrol)
mars <- function(formula,data,control=mars.control()) {
  cc <- match.call() # save the call
  mf <- model.frame(formula,data)
  y <- model.response(mf)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)[,-1,drop=FALSE]
  x_names <- colnames(x)
  control <- validate_mars.control(control)
  fwd <- fwd_stepwise(y,x,control)
  bwd <- bwd_stepwise(fwd,control)
  fit <- lm(y~.-1,data=data.frame(y=y,bwd$B)) # notice -1 added
  out <- c(list(call=cc,formula=formula,y=y,B=bwd$B,Bfuncs=bwd$Bfuncs,
                x_names=x_names),fit)
  class(out) <- c("mars",class(fit))
  out
}

#'Helper function for mars function
#'@param y a vector of y values
#' @param x data frame with the x values
#'
#'@return returns a list with y values, B, and Bfuncs
#'@example fwd_stepwise(y,x,control)
fwd_stepwise <- function(y,x,control=mars.control()){
  Mmax = control$Mmax
  N <- length(y)
  n <- ncol(x)
  B <- init_B(N,Mmax)
  Bfuncs <- vector("list", length = Mmax+1)
  #---------------------------------------------------
  # Error checking for Mmax:(write your code below)
  if(Mmax<2){
    warning("Mmax must be >= 2; setting to 2")
    Mmax=2
  }
  #---------------------------------------------------
  # Initialize:
  N <- length(y) # sample size
  n <- ncol(x) # number of predictors = number of X
  # B: a data frame with optimal basis function as columns
  B <- init_B(N,Mmax) # Exercise: write init_B()
  # splits: a data frame records all the optimal (m, v, t):
  # m: parent basis func to split, v: splitting var, t: splitting point
  # splits <- data.frame(m=rep(NA,Mmax),v=rep(NA,Mmax),t=rep(NA,Mmax)) # this can be replaced by Bfuncs
  # Bfuncs(): a list records the optimal (m, v, t) for basis functions (regions)
  Bfuncs <- vector("list", length = Mmax+1)
  #
  #---------------------------------------------------
  # Looping for forward selection:
  for(i in 1:(Mmax/2)) { # contrast to indexing 2...Mmax in Friedman

    lof_best <- Inf
    M <- (i*2) -1

    for(m in 1:M) { # choose a basis function to split
      svars <- setdiff(1:n,Bfuncs[[m]][,"v"])
      for(v in svars){ # select a variable to split on
        tt <- split_points(x[,v],B[,m]) # Exercise: write split_points()
        for(t in tt) {

          Bnew <- data.frame(B[,(1:M)], # drop m-th col: B[,-m]
                             # replace parent B[,m] with Btem1,Btem2
                             Btem1=B[,m]*h( x[,v],+1,t),
                             Btem2=B[,m]*h(x[,v],-1,t))
          gdat <- data.frame(y=y,Bnew)
          lof <- LOF(y~.,gdat,control) #  Use your LOF() from week 4
          if(lof < lof_best) {
            lof_best <- lof
            best_split <- c(m=m,v=v,t=t)
            # splits[M,] <- c(m,v,t)
          } # end if
        }# end loop over splits

      } # end loop over variables
    } # end loop over basis functions to split
    # save optimal (m, v, t) and update basis functions
    mstar <- best_split["m"]; vstar <- best_split["v"]; tstar <- best_split["t"]
    B[,M+1] <- (B[,mstar]*h(x[,vstar],-1,tstar))
    B[,M+2] <- (B[,mstar]*h(x[,vstar],+1,tstar))
    Bfuncs[[M +1]] <- rbind(Bfuncs[[mstar]], c(s= -1, vstar, tstar)) # Add left child basis function
    Bfuncs[[M+2]] <- rbind(Bfuncs[[mstar]], c(s=+1, vstar, tstar)) # Update parent basis with the right cihld basis

  } # end loop over M
  return(list(y=y,B=B, Bfuncs=Bfuncs))
}

#'Helper function for fwd stepwise
#'initializes the b matrix size N,Mmax
#'@param N number of rows
#' @param Mmax number of columns
#' @return a NxMmax matrix
#' @example  init_B(1,1)
init_B <- function(N,Mmax) {
  B <- data.frame(matrix(NA,nrow=N,ncol=(Mmax+1)))
  B[,1] <- 1
  names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)

}

#'helper for mars function
#'does the backward step for mars
#'@param fwd the results for the fwd stepwise
#'@param  controls a mars control object
#'@return list with y values, B values and Bfuncs
#'@example bwd_stepwise(fwd,control)
bwd_stepwise <- function(fwd,control) {
  Mmax = ncol(fwd$B)-1
  Jstar = 2:(Mmax+1)
  Kstar = Jstar
  data = data.frame(y=fwd$y,fwd$B)
  lofstar = LOF(y~.-1,data,control)
  for(M in (Mmax+1):2){
    L<-Kstar
    b <- Inf
    for(m in L){
      K <- setdiff(L,m)
      dat <- data.frame(y=fwd$y,fwd$B[,c(K)])
      lof <- LOF(y~.,dat,control)
      if(lof < b){
        b <-lof
        Kstar <- K
      }
      if(lof < lofstar){
        Jstar <- K
        lofstar <- lof
      }
    }
  }
  Jstar <- c(1,Jstar)
  return(list(y=fwd$y,B=fwd$B[,Jstar], Bfuncs=fwd$Bfuncs[Jstar]))
}

#' helper for fwd_stepwise and bwd_stepwise
#' comptues the LOF for a form and conrol
#' @param form formula used
#' @param data a data frame
#' @param  control a mars control object
#'
#' @return number for error
#' @example LOF(form,data,control)
LOF <- function(form,data,control) {
  # update this LOF to GCV
  d <-3
  mod<-lm(form,data)
  RSS<- sum((mod$res)^2)
  Ctilde<- sum(diag(hatvalues(mod))) + d*(length(coefficients(mod))-1 )
  sample_size = nrow(data)
  return( RSS*(sample_size/(sample_size-Ctilde)^2))
}


h <- function(x,s,t) {
  return(pmax(0,s*(x-t)))
}



#-------------------------------------

#if(Mmax<2) {
#  warning("Input Mmax must be >= 2; setting to 2")
#  Mmax <- 2
# }


split_points <- function(xv,Bm) {
  out <- sort(unique(xv[Bm>0]))
  return(out[-length(out)])
}

#------------------------------------------------------------------------
# constructor, validator and helper for class mars.control
#------------------------------------------------------------------------
#
new_mars.control <- function(control) {
  structure(control, class="mars.control")
}

validate_mars.control <- function(control) {
  stopifnot(is.integer(control$Mmax),is.numeric(control$d),is.logical(control$trace))
  if(control$Mmax < 2){
    warning("Mmax should be greater than 2")
    control$Mmax = 2
  }
  if(control$Mmax %% 2 != 0){
    control$Mmax = (ceiling(control$Mmax) + 1)
    warning("Mmax should be even, setting to ", control$Mmax)
  }
  control
}


#' Constructor for `mars.control` objects
#'
#' This function constructs a `mars.control` object that specifies
#' parameters used in the model fitting procedure.
#'
#' @param Mmax Maximum number of basis functions. Should be an even integer. Default value is 2.
# .....
# ...
mars.control <- function(Mmax=2,d=3,trace=FALSE) {
  Mmax <- as.integer(Mmax)
  control <- list(Mmax=Mmax,d=d,trace=trace)
  control <- validate_mars.control(control)
  new_mars.control(control)
}

#------------------------------------------------------------------------
# Methods
#------------------------------------------------------------------------
## Predict mars
#'Conducts analysis of variance on a mars object
#'
#' @param object an mars objects
#'
#' @return A data frame containing the sum of squares, F-value, p-value, etc.
#' @examples
#' mars.anova(mars)
predict.mars <- function(object,newdata) {
  if(missing(newdata) || is.null(newdata)) {
    B <- as.matrix(object$B)
  }
  else {
    tt <- terms(object$formula,data=newdata)
    tt <- delete.response(tt)
    mf <- model.frame(tt,newdata)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)[,-1] # remove intercept
    B <- make_B(X,object$Bfuncs)
  }
  beta <- object$coefficients
  drop(B %*% beta)
}

h <- function(x,s,t) {
  return(pmax(0,s*(x-t)))
}

make_B <- function(X, Bfuncs){
  #initialize matrix
  B<-matrix(1, nrow = nrow(X), ncol = length(Bfuncs))
  #loop over each row of X
  for(i in seq_len(nrow(X))){
    #loop over each Bfunc
    for(j in seq_len(length(Bfuncs))){
      temp<-1
      hnum<-length(Bfuncs[[j]][,1])
      #loop over each value in an individual Bfunc
      for(k in seq_len(hnum)){
        temp<-temp*h(Bfuncs[[j]][k,1],X[i,Bfuncs[[j]][k,2]],Bfuncs[[j]][k,3])
      }
      B[i,j]<-B[i,j]*temp
    }
  }
  return(B)
}

## ANOVA mars

#recall f value is calculated with MSR and MSE
#'Conducts analysis of variance on a mars object
#'
#' @param object an mars objects
#'
#' @return A data frame containing the sum of squares, F-value, p-value, etc.
#' @examples
#' mars.anova(mars)
anova.mars <- function(object) {
  ncoeffs <- object$coefficients
  y <- object$y
  y_hat <- decomp(object)
  y_bar <- mean(y)

  SSR <- sum((y_hat-y_bar)^2)
  SSE <- sum((y_hat-y)^2)
  SST <- sum((y-y_bar)^2)

  Rdf <- 1
  Edf <- length(y)-2
  Tdf <- length(y)-1

  MSR <- SSR/Rdf
  MSE <- SSE/Edf

  Fstat <- MSR/MSE

  pval <- pf(Fstat, Rdf, Edf, lower.tail = FALSE)

  table <- data.frame(c(Rdf,Edf,Tdf), c(SSR,SSE,SST),
                      c(MSR,MSE,""), c(Fstat,"",""), c(pval,"",""))
  names(table) <- c("df","Sum Sq","Mean Sq","F value","Pr(>F)")
  row.names(table) <- c("Regression","Residual","Total")

  return(table)
}

#anova decomposition
decomp <- function(object){
  ncoeffs <- length(object$coefficients)
  coeffs <- object$coefficients
  X <- object$model[,-1] # remove intercept

  f <- seq_len(nrow(X))*0
  for(i in seq_len(nrow(X))){
    for(j in seq_len(ncoeffs)){
      #loop over each coefficent
      f[i] <- X[i,j]*coeffs[j]
    }
  }
  return(f)
}


## Plot mars

plot.mars <- function(object){
  y <- object$y
  fit <- fitted(object)
  res <- residuals(object)
  strdres <- rstandard(object)
  sqrtstrdres <- strdres/((1/(length(y)-1))*sum(strdres^2))

  #residuals vs fitted
  plot(fit, res, ylab="Residuals", xlab="Fitted Values")
  title("Residuals vs Fitted Values")
  abline(h=0,lty=2)

  #qqplot
  readline(prompt="Press enter to continue to next plot")
  qqnorm(y)
  qqline(y)

  #Scale location
  readline(prompt="Press enter to continue to next plot")
  plot(fit,sqrtstrdres, ylab="Square Root of Standaradized Residuals", xlab="Fitted Values")
  title("Scale-Location")

  #Residual leverage
  readline(prompt="Press enter to continue to next plot")
  plot(hatvalues(object),strdres, ylab="Standaradized Residuals", xlab="Leverage")
  title("Residual vs Leverage")
}
#'Prints the mars object
#'
#' @param mars an mars objects
#' @param data a data frame containing the data
#'
#' @return void
#' @examples
#' print(mars)
print.mars<-function(mars)  {
  coeff <- mars$coefficients
  bfuncs <- mars$Bfuncs
  name <- names(coeff)
  index = 1
  varnames <- mars$x_names
  cat("Call:\n")
  print(mars$call)
  cat("\n")
  cat("Coefficients:\n")
  for(i in bfuncs){
    if(index == 1){ #handles the intercept
      cat(paste(name[1],": ",round(coeff[1],digits=7), "\n",sep=""))
    }
    else{
      equation <-c()
      for(ii in 1:nrow(i)){
        bfuncsrow <- i[ii,]
        if(bfuncsrow['s'] == 1){ #if the sign is negative
          if(bfuncsrow["t"] >=0){ #if the split value is positive
              str <- paste("max(0,", varnames[bfuncsrow["v"]]," - " ,round(bfuncsrow["t"],digits=7), ")",
                         sep="")
            equation <-c(equation,str)
          }
          else{#if the split value is negative
            str <- paste("max(0,", varnames[bfuncsrow["v"]]," + " ,abs(round(bfuncsrow["t"],digits=7)), ")",
                         sep="")
            equation <-c(equation,str)
          }

        }
        else{ #if the sign is positive
          str <- paste("max(0,", round(bfuncsrow["t"],digits=7)," - ",varnames[bfuncsrow["v"]], ")",
                       sep="")
          equation <-c(equation,str)

        }

      }
      cat(paste(name[index],": ",round(coeff[index],digits=7),sep=""))
      for(i in equation){ #prints the bfuncs
        cat(paste(" * ", i,sep=""))
      }
      cat("\n")
    }
    index = index +1

  }
}
