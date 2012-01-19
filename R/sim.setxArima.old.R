sim.setxArima <- function(object, x, x1 = NULL, num = 1000, prev = NULL,
                     bootstrap = FALSE, bootfn = NULL, cond.data = NULL,
                     max.iter=10, ...) {
  require(mvtnorm)
#  library.dynam("stats")
  if (bootstrap | !is.null(bootfn) | !is.null(cond.data)){
    warning("boostrap, bootfn, and cond.data are ignored in ARIMA models")
  }
  t.effect <- x$t.effect
  ##extracting out the series
  if (is.data.frame(object$zelig.data))
    dat <- object$zelig.data
  else{
    envir<- attr(object$terms, ".Environment")
    dat <- eval(object$call$data, envir)
	}
  series <- eval(eval(object$call$formula[[2]])$name, envir=dat)
  pred.ahead <- x$pred.ahead
  if (!is.null(prev)){
    draw.parm <- prev
  }
  if (is.null(prev)){
    draw.parm <- rmvnorm(num, mean=object$coef, sigma=object$var.coef)
###inserting a function that will be used to ensure that the 
###ma portion is invertible.  This code is taken from 
##the implementation of ARIMA in the stats package
    maInvert <- function(ma) {
        q <- length(ma)
        q0 <- max(which(c(1, ma) != 0)) - 1
        if (!q0) 
            return(ma)
        roots <- polyroot(c(1, ma[1:q0]))
        ind <- Mod(roots) < 1
        if (all(!ind)) 
            return(ma)
        if (q0 == 1) 
            return(c(1/ma[1], rep(0, q - q0)))
        roots[ind] <- 1/roots[ind]
        x <- 1
        for (r in roots) x <- c(x, 0) - c(0, x)/r
        c(Re(x[-1]), rep(0, q - q0))
    }
##
temp<- which(object$arma[1:4]>0)
temp<- sort(temp)
for(j in 1:nrow(draw.parm)){
if(length(temp)==1){
	if(temp==2 | temp==4){
		draw.parm[j, 1:object$arma[temp]]<- maInvert(draw.parm[i, 1:object$arma[temp]])
		}
	if(temp==1){
		draw.parm[j, 1:object$arma[temp]]<- .Call("ARIMA_transPars", draw.parm[j, 1:sum(object$arma[temp[1:length(temp)]])], 
			as.integer(object$arma[1:5]), TRUE, PACKAGE = "stats")[[1]][1:object$arma[temp]]
		}
	if(temp==3){
		draw.parm[j, 1:object$arma[temp]]<- .Call("ARIMA_transPars", draw.parm[j, 1:sum(object$arma[temp[1:length(temp)]])], 
			as.integer(object$arma[1:5]), TRUE, PACKAGE = "stats")[[1]][(1:object$arma[3])*object$arma[5]]
		}
	}
if(length(temp)==2){
	if(temp[1]==2 | temp[1]==4){
		draw.parm[i, 1:object$arma[temp[1]]]<- maInvert(draw.parm[j, 1:object$arma[temp[1]]])
		}
	if(temp[1]==1){
		draw.parm[j, 1:object$arma[temp[1]]]<- .Call("ARIMA_transPars", draw.parm[j, 1:sum(object$arma[temp[1:length(temp)]])], as.integer(object$arma[1:5]), TRUE, PACKAGE = "stats")[[1]][1:object$arma[temp[1]]]
		}
	if(temp[1]==3){
		draw.parm[j, 1:object$arma[temp[1]]]<- .Call("ARIMA_transPars", draw.parm[j, 1:sum(object$arma[temp[1:length(temp)]])], 
			as.integer(object$arma[1:5]), TRUE, PACKAGE = "stats")[[1]][(1:object$arma[3])*object$arma[5]]
		}
	if(temp[2]==2 | temp[2]==4){
		draw.parm[j, (object$arma[temp[1]] + 1):(sum(object$arma[temp[1]:temp[2]]))]<- maInvert(draw.parm[j, (object$arma[temp[1]] + 1):(sum(object$arma[temp[1]:temp[2]]))])
		}
	if(temp[2]==1){
		draw.parm[j, (sum(object$arma[temp[1]]) + 1):(sum(object$arma[temp[1]:temp[2]]))]<- .Call("ARIMA_transPars",draw.parm[j, 1:sum(object$arma[temp[1:length(temp)]])],
													 as.integer(object$arma[1:5]), TRUE, PACKAGE = "stats")[[1]][(sum(object$arma[1]) + 1):(sum(object$arma[1:2]))]
}
	if(temp[2]==3){
		draw.parm[j, (sum(object$arma[temp[1]]) + 1):object$arma[temp[2]]]<- .Call("ARIMA_transPars", draw.parm[j, 1:sum(object$arma[temp[1:length(temp)]])], 
			as.integer(object$arma[1:5]), TRUE, PACKAGE = "stats")[[1]][(1:object$arma[3])*object$arma[5]]

  }
}
if(length(temp)>2){
if(temp[1]==2 | temp[1]==4){
		draw.parm[i, 1:object$arma[temp[1]]]<- maInvert(draw.parm[j, 1:object$arma[temp[1]]])
		}
	if(temp[1]==1){
		draw.parm[j, 1:object$arma[temp[1]]]<- .Call("ARIMA_transPars", draw.parm[j, 1:sum(object$arma[temp[1:length(temp)]])], 
				as.integer(object$arma[1:5]), TRUE, PACKAGE = "stats")[[1]][1:object$arma[temp[1]]]
		}
	if(temp[1]==3){
		draw.parm[j, 1:object$arma[temp[3]]]<- .Call("ARIMA_transPars", draw.parm[j, 1:sum(object$arma[temp[1:length(temp)]])], 
			as.integer(object$arma[1:5]), TRUE, PACKAGE = "stats")[[1]][(1:object$arma[3])*object$arma[5]]
		}
for(i in 2:length(temp)){
if(temp[i]==2 | temp[i]==4){
		draw.parm[j, (sum(object$arma[temp[1]:temp[(i-1)]]) + 1):(sum(object$arma[temp[1]:temp[(i)]]))]<- maInvert(draw.parm[j, 
						(sum(object$arma[temp[1]:temp[(i-1)]]) + 1):(sum(object$arma[temp[1]:temp[(i)]]))])
		}
	if(temp[i]==1){
		draw.parm[j, (sum(object$arma[temp[1]:temp[(i-1)]]) + 1):(sum(object$arma[temp[1]:temp[(i)]]))]<- .Call("ARIMA_transPars",draw.parm[j, 1:sum(object$arma[temp[1:length(temp)]])],
													 as.integer(object$arma[1:5]), TRUE, PACKAGE = "stats")[[1]][(sum(object$arma[temp[1]:temp[(i-1)]]) + 1):(sum(object$arma[temp[1]:temp[(i)]]))]
}
	if(temp[i]==3){
		draw.parm[j, (sum(object$arma[temp[1]:temp[(i-1)]]) + 1):(sum(object$arma[temp[1]:temp[(i)]]))]<- .Call("ARIMA_transPars", draw.parm[j, 1:sum(object$arma[temp[1:length(temp)]])], 
			as.integer(object$arma[1:5]), TRUE, PACKAGE = "stats")[[1]][(1:object$arma[3])*object$arma[5]]
		}

  }
}
}
}
  if (x$min.time==1 | x$min.time==2)
    stop("Counterfactuals can only be specified from the third observation and later \n")
  if (ncol(x$dta)==0){
    ev <- matrix(NA, nrow=nrow(draw.parm), ncol=(pred.ahead))
    se <- matrix(NA, nrow=nrow(draw.parm), ncol=(pred.ahead))
    for (i in 1:nrow(draw.parm)){
      temp <- arima(series, xreg=NULL,
                    order=c(object$arma[1], object$arma[6], object$arma[2]), 
                    seasonal=list(order=c(object$arma[3], object$arma[7], object$arma[4]),
                      period=object$arma[5]), fixed=draw.parm[i,], transform.pars=TRUE)
      temp2 <-predict(temp, newxreg=NULL, n.ahead=(pred.ahead))
      ev[i,] <- temp2$pred[1:(pred.ahead)]
      se[i,] <- temp2$se[1:(pred.ahead)]	
    }
    ev <- as.matrix(ev)
    se <- as.matrix(se)
    if (!is.null(x1)){
      warning("First differences are only calculated when external regressors are used \n")
    }
  }
  if (ncol(x$dta) > 0){
    x.obs <- as.matrix(x$dta[1:(x$min.time-1),])
    x.cf <- as.matrix(x$dta[x$min.time: (x$max.time), ])
    ev <- matrix(NA, nrow=nrow(draw.parm), ncol=(length(x$min.time:x$max.time)))
    se <- matrix(NA, nrow=nrow(draw.parm), ncol=(length(x$min.time:x$max.time)))
    if (x$min.time == nrow(x$dta)){
      x.cf<- t(x.cf)
    }
    for (i in 1:nrow(draw.parm)){
      temp <- arima(series[1:(x$min.time-1)], xreg=x.obs,
                    order=c(object$arma[1], object$arma[6], object$arma[2]), 
                    seasonal=list(order=c(object$arma[3], object$arma[7], object$arma[4]),
                      period=object$arma[5]), fixed=draw.parm[i,])
      temp2 <-predict(temp, newxreg=x.cf, n.ahead=nrow(x.cf))
      ev[i,] <- temp2$pred
      se[i,] <- temp2$se
    }
    ev <- as.matrix(ev)
  }
  if (!is.null(x1)){
    if (ncol(x1$dta) > 0){
      x1.obs <- as.matrix(x1$dta[1:(x1$min.time-1), ])
      x1.cf <- as.matrix(x1$dta[x1$min.time:(x1$max.time),])
      ev.1 <- matrix(NA, nrow=nrow(draw.parm), ncol=length(x$min.time:x$max.time))
      se.1 <- matrix(NA, nrow=nrow(draw.parm), ncol=length(x$min.time:x$max.time))
      pred.ahead <- x$pred.ahead
      pred.arima.x1 <- list()
      false.arima.x1 <- list() 
      if (x1$min.time == nrow(x1$dta)){
	x1.cf <- t(x1.cf)
      }
      for (i in 1:nrow(draw.parm)){
	temp3<- arima(series[1:(x1$min.time-1)], xreg=x1.obs,
                                     order=c(object$arma[1], object$arma[6], object$arma[2]),
                                    seasonal=list(order=c(object$arma[3], object$arma[7],
                                                    object$arma[4]), period=object$arma[5]),
                                     fixed=draw.parm[i,])
	temp4 <- predict(temp3, newxreg=x1.cf, n.ahead=nrow(x.cf))
	ev.1[i,] <- temp4$pred	
	se.1[i,] <- temp4$se
      }
      ev.1 <- as.matrix(ev.1)
    }
  }
  if (!is.null(x1) & t.effect){
    warning("First differences and treatment effects are not both calculated.Calculating first differences only.\n")
  }
  if (is.null(x1) & (t.effect)){
    t.eff <- matrix(NA, nrow=nrow(ev), ncol=ncol(ev))
    for (i in 1:nrow(ev)){
      t.eff[i,] <- series[(x$min.time):x$max.time] - ev[i,]
    }
  }
  if (!is.null(x1) | !(t.effect)){
    t.eff <- NULL
  }
  if (!is.null(x1)){
    se <- as.matrix(se)
    se.1 <- as.matrix(se.1)
    fd <- ev.1 -  ev
    qi <- list(ev=ev, se=se, fd=fd, t.eff=t.eff)
    qi.name <- list(ev="Expected Values, E(Y|X)", se="Prediction Standard Error",
                    fd="First Difference, E(Y|X1) - E(Y|X)", t.eff="Treatment Effect")
    res <- list(min.time=x$min.time, qi=qi, qi.name=qi.name,
                zelig.call = object$call, 
                t.series = eval(eval(object$call$formula[[2]],
                  envir=dat)$name, envir=dat))
    class(res) <- c("zelig.arima", "zelig")
    return(res)
  } 
  if (is.null(x1)){
    se <- as.matrix(se)
    qi <- list(ev=ev, se=se, t.eff=t.eff)
    qi.name <- list(ev="Expected Values, E(Y|X)", se="Prediction Standard Error",
                    t.eff="Treatment Effect") 
    res <- list(min.time=x$min.time, qi=qi, qi.name=qi.name,
                t.series=eval(eval(object$call$formula[[2]],
                  envir=dat)$name, envir=dat),
                zelig.call = object$call)
    class(res) <- c("zelig.arima", "zelig")
    return(res) 
  } 
} 

