zelig2arima <- function(formula, model, data, M, ...){
  mf <- match.call()
  ##arima.wrap allows arima to run without using a formula
  mf[[1]] <- arima.wrap
  ##assembling the order of the ARIMA model
  dep.var <- eval(mf[[2]][[2]])$name
  d <- eval(mf[[2]][[2]])$d
  d.s <- eval(mf[[2]][[2]])$ds
  per <- eval(mf[[2]][[2]])$per
  ##the don't use vector will be used below to assemble the actual formula employed
  dont.use <- vector()
  ##this for loop is used to find the values of p, p.s., q, and q.s.
  n.par <- length(unlist(strsplit(deparse(mf[[2]][[3]], width.cutoff=500), "\\+")))
  for(i in 1:n.par){
    if(class(eval(parse(text=unlist(strsplit(deparse(mf[[2]][[3]], width.cutoff=500), "\\+"))[i]),
                  envir=data))=="list"){
      if(eval(parse(text=unlist(strsplit(deparse(mf[[2]][[3]], width.cutoff=500), "\\+"))[i]),
              envir=data)$ep==TRUE){
        dont.use[1] <- i
        q <- eval(parse(text=unlist(strsplit(deparse(mf[[2]][[3]], width.cutoff=500), "\\+"))[i]),
                  envir=data)$q
        q.s <- eval(parse(text=unlist(strsplit(deparse(mf[[2]][[3]], width.cutoff=500), "\\+"))[i]),
                    envir=data)$qs
      }
      if(eval(parse(text=unlist(strsplit(deparse(mf[[2]][[3]], width.cutoff=500), "\\+"))[i]),
              envir=data)$y==TRUE){
        dont.use[2] <- i
        p <- eval(parse(text=unlist(strsplit(deparse(mf[[2]][[3]], width.cutoff=500), "\\+"))[i]),
                  envir=data)$p
        p.s <- eval(parse(text=unlist(strsplit(deparse(mf[[2]][[3]], width.cutoff=500), "\\+"))[i]),
                    envir=data)$ps
      }
    }
  }
  
  if (length(dont.use) < n.par){
    ##this vector then denotes what we are going to use on the right hand side
    use.vec <- c(1:length(unlist(strsplit(deparse(mf[[2]][[3]], width.cutoff=500), "\\+"))))
    use.vec <- use.vec[-dont.use]
    ##now, we are creating the right hand side portion of the formula
    rhs <- vector()
    for (i in 1:length(use.vec)){
      rhs[2*i-1] <- paste(unlist(strsplit(deparse(mf[[2]][[3]], width.cutoff=500), "\\+"))[use.vec[i]])
      if (i+1 <= length(use.vec))
        rhs[2*i] <- paste("+")
      if (i + 1 > length(use.vec))
       	break  
    }
    mf$formula <-formula <- as.formula(paste(paste(dep.var,"~", collapse=""),
                                             paste(rhs, collapse=""), collapse=""))
    D <- model.frame(mf$formula, data)
    cols.D <- colnames(D)
    mf$x <- D[,1]
    X <- as.matrix(cbind(D[,2:ncol(D)]))
    colnames(X) <- cols.D[2:length(cols.D)]
    mf$xreg <- X
  }
  if(length(dont.use)==n.par){
    mf$x <- as.matrix(eval(dep.var, envir=data))
    mf$xreg <- NULL
    mf$formula <- eval(dep.var, envir=data)~1 
  }	
  mf$order<- c(p,d,q) 
  if(!is.null(d.s) & !is.null(p.s) & !is.null(q.s) & !is.null(per)){	
    mf$seasonal$order <- c(p.s, d.s, q.s)
    mf$seasonal$period <- per
  }
  mf$model <- mf$M <- mf$data <- NULL
  as.call(mf) 
}
