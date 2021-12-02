##using new version of MOEstart ----- working good but error with some data
#using the second transformation calculating the rth moment
#all function to estimate parameters of MOE "Method of Moments"
MOEmmedist<-function (data, distr, memp, start = NULL, fix.arg = NULL, optim.method = "default", lower = -Inf, upper = Inf, custom.optim = NULL, weights = NULL, silent = TRUE, gradient = NULL, ...) 
{#MOEmmedist
	 if (!is.character(distr)) 
 	stop("distr must be a character string naming a distribution")
 	else distname <- distr

	if (!(is.numeric(data) & length(data) > 1)) 
	stop("data must be a numeric vector of length greater than 1")

	if (missing(memp)) 
	stop("the empirical moment function must be defined")

	if (is.character(memp)) 
	memp <- get0(memp, envir = pos.to.env(1))
        	
	if (!is.function(memp)) 
	stop("the empirical moment must be defined as a function")
meth <- optim.method
ddistname <- paste("d", distname, sep = "")
argddistname <- names(formals(ddistname))
	if (is.null(custom.optim)) 
	optim.method <- match.arg(optim.method, c("default", "Nelder-Mead", "BFGS", "CG", 	"L-BFGS-B", "SANN", "Brent"))

	if (!is.null(weights)) 
	{#if (!is.null(weights))
		if (any(weights < 0)) 
		stop("weights should be a vector of integers greater than 0")

		if (!is.allint.w(weights)) 
            		stop("weights should be a vector of (strictly) positive integers")

		if (length(weights) != NROW(data)) 
        		stop("weights should be a vector with a length equal to the observation 			number")
	warning("weights are not taken into account in the default initial values")
		# ----- weighted log liklehood function of MOE distribution ----   
		loglik <- function(par, fix.arg, obs, distr,weights) 
		{#loglik 
			if (!is.list(par))
			par<-as.list(par)
		#exclude tilt from par or fixed.arg
		tiltpar<-match("tilt",names(par))
			if(!is.na(tiltpar))
			{
			alpha<-par$tilt
			par<-par[-tiltpar] #exclude tilt parameter from par
			}
		tiltfix<-match("tilt",names(fix.arg))
			if(!is.na(tiltfix))
			{
			alpha<-fix.arg$tilt
			fix.arg<-fix.arg[-tiltfix] #exclude tilt parameter from par
			}
		ddistname <- paste("d", distr, sep = "")
		pdistname <- paste("p", distr, sep = "")
		fbase<-do.call(ddistname, c(list(obs), as.list(par),as.list(fix.arg)))
		Fbase<-do.call(pdistname, c(list(obs), as.list(par),as.list(fix.arg)))
		n<-length(obs)
		alpha_bar<-1-alpha
		LL <- log(alpha)*sum(weights)+ sum(weights*log(fbase))-2*sum(weights*log(1-		alpha_bar*(1-Fbase)))
		return(LL)
		} #loglik 
		# ------------------------------------------------------------------------------------- 
	txt <- "the empirical moment function must be a three-argument function of 'x', 	'order',weights'"
	if (length(formals(memp)) != 3) 
	stop(txt)

	if (any(names(formals(memp)) != c("x", "order", "weights"))) 
	stop(txt)

		DIFF2 <- function(par, fix.arg, order, obs, distr, memp, weights) 
		{
		momtheo <-mMOE(order,par,fix.arg,distr)
		momemp <- as.numeric(memp(obs, order, weights)) 
		(momemp - momtheo)^2
            		}

            		fnobj <- function(par, fix.arg, obs, distr, memp,  weights) 					sum(sapply(order, function(o) DIFF2(par,fix.arg, o, obs, distr, memp,weights)))  

	} #if (!is.null(weights))
	else
	{#else (!is.null(weights))
		# ----- log liklehood function of MOE distribution ----   
		loglik <- function(par, fix.arg, obs, distr) 
		{#loglik 
			if (!is.list(par))
			par<-as.list(par)
		#exclude tilt from par or fixed.arg
		tiltpar<-match("tilt",names(par))
			if(!is.na(tiltpar))
			{
			alpha<-par$tilt
			par<-par[-tiltpar] #exclude tilt parameter from par
			}
		tiltfix<-match("tilt",names(fix.arg))
			if(!is.na(tiltfix))
			{
			alpha<-fix.arg$tilt
			fix.arg<-fix.arg[-tiltfix] #exclude tilt parameter from par
			}
		ddistname <- paste("d", distr, sep = "")
		pdistname <- paste("p", distr, sep = "")
		fbase<-do.call(ddistname, c(list(obs), as.list(par),as.list(fix.arg)))
		Fbase<-do.call(pdistname, c(list(obs), as.list(par),as.list(fix.arg)))
		n<-length(obs)
		alpha_bar<-1-alpha
		LL <- n*log(alpha)+ sum(log(fbase))-2*sum(log(1-alpha_bar*(1-Fbase)))
		return(LL)
		}#loglik 
		# ------------------------------------------------------------------------------------- 
	txt <- "the empirical moment function must be a two-argument function of 'x', 'order'"
	if (length(formals(memp)) != 2) 
	stop(txt)

	if (any(names(formals(memp)) != c("x", "order"))) 
	stop(txt)

		DIFF2 <- function(par, fix.arg, order, obs, distr, memp) 
		{
               	momtheo <-mMOE(order,par,fix.arg,distr)
             		momemp <- as.numeric(memp(obs, order)) 
		(momemp - momtheo)^2
            		}

	fnobj <- function(par, fix.arg, obs, distr, memp) 
	sum(sapply(order, function(o) DIFF2(par, fix.arg, o, obs, distr, memp)))
	} #else (!is.null(weights))

	if (is.vector(start)) 
	start <- as.list(start)
arg_startfix <- MOEmanageparam(start.arg = start, fix.arg = fix.arg,obs = data, 	distname= distr)
hasnodefaultval <- sapply(formals(ddistname), is.name) # True for arguments does 	not have default value and False for other arguments (ex. 	sapply(formals(dgamma),is.name))
arg_startfix <- MOEcheckparamlist(arg_startfix$start.arg,arg_startfix$fix.arg, 	argddistname, hasnodefaultval) 
	if (is.function(fix.arg)) 
	fix.arg.fun <- fix.arg
	else fix.arg.fun <- NULL
vstart <- unlist(arg_startfix$start.arg)
fix.arg <- arg_startfix$fix.arg
order<-length(vstart)      	
owarn <- getOption("warn")
	if (is.null(custom.optim)) 
	{#if (is.null(custom.optim)) 
	hasbound <- any(is.finite(lower) | is.finite(upper))
		if (optim.method == "default") 
		{
               	opt.meth <- ifelse(length(vstart) > 1, "Nelder-Mead", "BFGS")
            		}
           		else opt.meth <- optim.method

            		if (opt.meth == "BFGS" && hasbound && is.null(gradient)) 
		{
	               opt.meth <- "L-BFGS-B"
               	txt1 <- "The BFGS method cannot be used with bounds without provided the 			gradient."
               	txt2 <- "The method is changed to L-BFGS-B."
               	warning(paste(txt1, txt2))
            		}
	options(warn = ifelse(silent, -1, 0))
		if (hasbound) 
		{#if (hasbound) 
               		if (!is.null(gradient)) 
			{
               		opt.fun <- "constrOptim"
                		}
                		else 
			{
                  			if (opt.meth == "Nelder-Mead") 
                    			opt.fun <- "constrOptim"
                  			else if (opt.meth %in% c("L-BFGS-B", "Brent")) 
                    			opt.fun <- "optim"
                  			else 
				{
                  			txt1 <- paste("The method", opt.meth,"cannot be used by 					constrOptim() nor optim() without gradient and bounds.")
                    			txt2 <- "Only optimization methods L-BFGS-B, Brent and 					Nelder-Mead can be used in such case."
                    			stop(paste(txt1, txt2))
                  			}
                		}
                		if (opt.fun == "constrOptim") 
			{
                 		npar <- length(vstart)
                  		lower <- as.double(rep_len(lower, npar))
                  		upper <- as.double(rep_len(upper, npar))
                  		haslow <- is.finite(lower)
                  		Mat <- diag(npar)[haslow, ]
                  		hasupp <- is.finite(upper)
                  		Mat <- rbind(Mat, -diag(npar)[hasupp, ])
                  		colnames(Mat) <- names(vstart)
                  		rownames(Mat) <- paste0("constr", 1:NROW(Mat))
                  		Bnd <- c(lower[is.finite(lower)], -upper[is.finite(upper)])
                  		names(Bnd) <- paste0("constr", 1:length(Bnd))
                  		initconstr <- Mat %*% vstart - Bnd
                  			if (any(initconstr < 0)) 
                    			stop("Starting values must be in the feasible region.")
                  		opttryerror <- try(opt <- constrOptim(theta = vstart,f = fnobj, ui = Mat, 				ci =Bnd, grad = gradient,fix.arg = fix.arg, obs = data, mdistnam = 				mdistname,memp = memp, hessian = !is.null(gradient),method = 				opt.meth, weights = weights, ...), silent = TRUE)
                  			if (!inherits(opttryerror, "try-error")) 
                    				if (length(opt$counts) == 1) 
                      				opt$counts <- c(opt$counts, NA)
			}
			else
			{
			opttryerror <- try(opt <- optim(par = vstart,fn = fnobj, fix.arg = fix.arg, 				obs =data, gr = gradient, distr = distr, memp = memp, 						hessian = TRUE, method = opt.meth, lower = lower,upper = upper, 				weights = weights, ...), silent = TRUE)
			}
 		}#if (hasbound) 
		else 
		{#else(hasbound) 
               	opt.fun <- "optim"
		opttryerror <- try(opt <- optim(par = vstart,fn = fnobj, fix.arg = fix.arg, obs = 			data,gr = gradient, distr = distr, memp = memp, hessian = TRUE, 				method =opt.meth, lower = lower,upper = upper, ...), silent = TRUE)
		}
	options(warn = owarn)
		if (inherits(opttryerror, "try-error")) 
		{
              		warnings("The function optim encountered an error and stopped.")
              			if (getOption("show.error.messages")) 
                		print(attr(opttryerror, "condition"))
              		return(list(estimate = rep(NA, length(vstart)), convergence = 100, value = NA, 			hessian = NA))
		}
		if (opt$convergence > 0) 
		{
		warnings("The function optim failed to converge, with the error code ", opt			$convergence)
            		}
		if (is.null(names(opt$par))) 
              	 	names(opt$par) <- names(vstart)
	res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 	hessian = opt$hessian, optim.function = opt.fun, optim.method = opt.meth, fix.arg = 	fix.arg, fix.arg.fun = fix.arg.fun, weights = weights, counts = opt$counts, 	optim.message = opt$message,loglik = ifelse(exists(ddistname), loglik(opt$par, fix.arg, 	data, distr), NULL), method = meth, order = order, memp = memp)
	} #if (is.null(custom.optim)) 
	else 
	{ #else is.null(custom.optim)
	opt.meth <- NULL
	options(warn = ifelse(silent, -1, 0))
	opttryerror <- try(opt <- custom.optim(fn = fnobj,fix.arg = fix.arg, obs = data, distr = 	distr, memp = memp, par = vstart, weights = weights, ...), silent =TRUE)
	options(warn = owarn)
	          	if (inherits(opttryerror, "try-error")) 
		{
		warnings("The customized optimization function encountered an error and 			stopped.")
			if (getOption("show.error.messages")) 
			print(attr(opttryerror, "condition"))
		return(list(estimate = rep(NA, length(vstart)), convergence = 100, value = NA, 			hessian = NA))
		}
		if (opt$convergence > 0) 
		{
		warnings("The customized optimization function failed to converge, with the			error code ", opt$convergence)
		}
		if (is.null(names(opt$par))) 
		names(opt$par) <- names(vstart)
	argdot <- list(...)
	method.cust <- argdot$method
	res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 	hessian = opt$hessian, optim.function = custom.optim, optim.method = method.cust, 	fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, weights = weights, counts = opt$counts, 	optim.message = opt$message, loglik = ifelse(exists(ddistname), loglik(opt$par, fix.arg, 	data, distr), NULL), method = meth,  order = order, memp = memp)
	}#else is.null(custom.optim)
return(res)
}#MOEmmedist
#---------------1-MOEstart------------
#------- new version of MOEstart ---
MOEstart<- function(x, distr)
{
pdistname <- paste("p", distr, sep = "")
if (!exists(pdistname, mode = "function")) 
  stop(paste("The ", ddistname, " function must be defined"))
n<-length(x)
alpha0<-2
#objfn<-function(obs,alpha)
#{
#n<-length(obs)
#minusloglik <- -n*log(alpha)+sum(obs)+2*sum(log(1-(1-alpha)*exp(-obs)))
#return(minusloglik)
#}
if (distr == "norm") 
{
sd0 <- sqrt((n - 1)/n) * sd(x)
mx <- mean(x)
par1<-list(mean=mx, sd=sd0)
}
else if (distr == "exp") 
{
rate0=1/mean(x)
par1<-list(rate=rate0)
}
else if (distr == "lnorm") 
{
    if (any(x <= 0)) 
   stop("values must be positive to fit a lognormal distribution")
    lx <- log(x)
    sd0 <- sqrt((n - 1)/n) * sd(lx)
    ml <- mean(lx)
    par1<- list(meanlog=ml, sdlog=sd0)
  }
else if (distr == "gamma") 
{
    if (any(x < 0)) 
      stop("values must be positive to fit an gamma  distribution")
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    par1<- list(shape=m^2/v, rate=m/v)
  }
else if (distr == "beta") 
{
    if (any(x < 0) | any(x > 1)) 
      stop("values must be in [0-1] to fit a beta distribution")
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    aux <- m*(1-m)/v - 1
    par1<- list(shape1=m*aux, shape2=(1-m)*aux)
  }
else if (distr == "weibull") 
{
    if (any(x < 0)) 
      stop("values must be positive to fit an Weibull  distribution")
    m <- mean(log(x))
    v <- var(log(x))
    shape <- 1.2/sqrt(v)
    scale <- exp(m + 0.572/shape)
    par1<- list(shape = shape, scale = scale)
  }
else if (distr == "logis") 
{
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    par1<- list(location=m, scale=sqrt(3*v)/pi)
}
else if (distr == "cauchy") 
{
    par1<- list(location=median(x), scale=IQR(x)/2)
}
else if (distr == "unif")
{
par1<- list(min=0, max=1)
}
else if (distr == "invgamma")
{
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse gamma  distribution")
    #http://en.wikipedia.org/wiki/Inverse-gamma_distribution
    m1 <- mean(x)
    m2 <- mean(x^2)
    shape <- (2*m2-m1^2)/(m2-m1^2)
    scale <- m1*m2/(m2-m1^2)
    par1<- list(shape=shape, scale=scale)
  }
else if (distr == "llogis")
 {
    if (any(x < 0)) 
      stop("values must be positive to fit a log-logistic  distribution")
    p25 <- as.numeric(quantile(x, 0.25))
    p75 <- as.numeric(quantile(x, 0.75))
    shape <- 2*log(3)/(log(p75)-log(p25))
    scale <- exp(log(p75)+log(p25))/2
    par1<- list(shape=shape, scale=scale)
}
else if (distr == "invweibull")
{
    if (any(x < 0)) 
    stop("values must be positive to fit an inverse Weibull distribution")
    g <- log(log(4))/(log(log(4/3)))
    p25 <- as.numeric(quantile(x, 0.25))
    p75 <- as.numeric(quantile(x, 0.75))
    shape <- exp((g*log(p75)-log(p25))/(g-1))
    scale <-log(log(4))/(log(shape)-log(p25))
    par1<- list(shape=shape, scale=max(scale, 1e-9))
}
else if (distr == "pareto1")
 {
    if (any(x < 0)) 
    stop("values must be positive to fit a Pareto distribution")
    #http://www.math.umt.edu/gideon/pareto.pdf
    x1 <- min(x)
    m1 <- mean(x)
    n <- length(x)
    shape <- (n*m1-x1)/(n*(m1-x1))
    min <- x1*(n*shape - 1)/(n*shape)
    par1<- list(shape=shape, min=min)
  }
else if (distr == "pareto")
  {
    if (any(x < 0)) 
    stop("values must be positive to fit a Pareto distribution")
    m1 <- mean(x)
    m2 <- mean(x^2)
    scale <- (m1*m2)/(m2-2*m1^2)
    shape <- 2*(m2-m1^2)/(m2-2*m1^2)
    par1<- list(shape=shape, scale=scale)
  }
else if (distr == "lgamma")
  {
    if (any(x < 0)) 
    stop("values must be positive to fit a log-gamma distribution")
    #p228 of Klugmann and Hogg (1984)
    m1 <- mean(log(x))
    m2 <- mean(log(x)^2)
    alpha <- m1^2/(m2-m1^2)
    lambda <- m1/(m2-m1^2)
    par1<- list(shapelog=alpha, ratelog=lambda)
  }
else if (distr == "trgamma") 
{
    if (any(x < 0)) 
    stop("values must be positive to fit an trans-gamma  distribution")
    #same as gamma with shape2=tau=1
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    par1<- list(shape1=m^2/v, shape2=1, rate=m/v)
  }
else if (distr == "invtrgamma") 
{
    if (any(x < 0)) 
    stop("values must be positive to fit an inverse trans-gamma  distribution")
    #same as gamma with shape2=tau=1
    n <- length(1/x)
    m <- mean(1/x)
    v <- (n - 1)/n*var(1/x)
    par1<- list(shape1=m^2/v, shape2=1, rate=m/v)
  }
else
stop(paste0("Unknown starting values for distribution ", distr, "."))
CDF<-do.call(pdistname, c(list(median(data)), as.list(par1)))
alpha<-CDF/(1-CDF)
start <- c(tilt=alpha,par1)
return(start)
}
# ------------------2- MOEmanageparam---------------------
#typical to manageparam
# checkparam function checks start.arg and fix.arg that parameters are named correctly
# INPUTS 
# start.arg : starting values for optimization or the function to compute them from data or NULL
# fix.arg : fixed values of paramaters or the function to compute them from data or NULL
# obs : the full dataset
# distname : name of the distribution
# OUTPUTS 
# two named list with untested components
MOEmanageparam <- function(start.arg, fix.arg, obs, distname)
{
# if clause with 3 different cases:
# start.arg : NULL | named list | a function
if (is.null(start.arg))
{      
trystart <- try(MOEstart(obs, distname), silent = TRUE)
	if (class(trystart) == "try-error")     
	stop("Error in computing default starting values:\n" ,trystart)
lstart <- trystart
# lstart should be a named list but check it
hasnoname <- is.null(names(lstart)) || !is.list(lstart)
	if(hasnoname)                              
	stop("Starting values must be a named list, error in default starting value.")       
}	                                
else if (is.list(start.arg))      
{ 
# start.arg should be a named list but check it
hasnoname <- is.null(names(start.arg))
    	if(hasnoname)
      	stop("Starting values must be a named list (or a function returning a 			named list).")
lstart <- start.arg
}
else if(is.function(start.arg))
{
trystart <- try(start.arg(obs), silent = TRUE)
	if(class(trystart) == "try-error")       
	stop("Error in computing starting values with your function.			\n",trystart)          
lstart <- trystart              
hasnoname <- is.null(names(lstart)) || !is.list(lstart)
    	if(hasnoname)                 
      	stop("Starting values must be a named list, your function does not 			return that.")         
}                                       
else                                                                  
stop("Wrong type of argument for start")   
# if clause with 3 different cases:
# fix.arg : NULL | named list | a function
if(is.null(fix.arg))
{                                        
lfix <- NULL                          
}                                   
else if(is.list(fix.arg))
{                                
hasnoname <- is.null(names(fix.arg))          
	if(hasnoname)          
	stop("Fixed parameter values must be a named list (or a function returning a 	named list).")
lfix <- fix.arg
}
else if(is.function(fix.arg))
{
tryfix <- try(fix.arg(obs), silent = TRUE)
	if(class(tryfix) == "try-error")
      	stop("Error in computing fixed parameter values with your function.\n",tryfix)
lfix <- tryfix
hasnoname <- is.null(names(lfix)) || !is.list(lfix)
	if(hasnoname)
	stop("Fixed parameter values must be a named list, your function does not 	return that.")
 }                                     
else                                      
stop("Wrong type of argument for fix.arg")
#eliminate arguments both in lstart and lfix (when start.arg was NULL)
if(is.null(start.arg) && !is.null(lfix))
 {
    lstart <- lstart[!names(lstart) %in% names(lfix)]
    if(length(lstart) == 0)
      stop("You don't need to fit the distiburtion if all parameters have fixed values")
 }
list("start.arg"=lstart, "fix.arg"=lfix)
}
#--------------3- MOEcheckparamlist --------------------------------
# checkparam function checks start.arg and fix.arg that parameters are named correctly
# INPUTS 
# start.arg : a named list
# fix.arg : NULL or a named list
# argdistname : argument names of the distribution
# hasnodefaultval : vector of logical indicating no default value of argument
# OUTPUTS 
# a named list with components: ok (TRUE or FALSE), txt (NULL or the error message), 
# start.arg : a named list of starting values for optimization 
# or a function to compute them from data
MOEcheckparamlist <- function(start.arg, fix.arg, argdistname, hasnodefaultval)
{
errtxt <- list(t3="'start' must specify names which are arguments to 'distr'.", t4="'fix.arg' must specify names which are arguments to 'distr'.", t5="A distribution parameter cannot be specified both in 'start' and 'fix.arg'.", t6="'start' should not have NA or NaN values.", t7="'fix.arg' should not have NA or NaN values.", t8="Some parameter names have no starting/fixed value: ", t9="Some parameter names have no starting/fixed value but have a default value: ")
#detect tilt parameter in start.arg
tiltstart<-match("tilt",names(start.arg))
if(!is.na(tiltstart))
nottiltstart.arg<-start.arg[-tiltstart] #exclude tilt parameter from start.arg
else
nottiltstart.arg<-start.arg
vstart <- unlist(start.arg) # including tilt parameter
vnottiltstart <- unlist(nottiltstart.arg) #not including tilt parameter
# check unexpected names
m <- match(names(vnottiltstart ), argdistname)
if (any(is.na(m))) 
stop(errtxt$t3)
#check NA/NaN values
if(any(is.na(vstart)) || any(is.nan(vstart)))
stop(errtxt$t6)

if(!is.null(fix.arg))
{
#detect tilt parameter in fix.arg
tiltfix<-match("tilt",names(fix.arg))
	if (is.na(tiltstart) & is.na(tiltfix))
	stop(errtxt$t8, " tilt parameter.")
                                                                                                                                                           
	if(!is.na(tiltfix))
	{ 
		if (fix.arg[tiltfix]==1)
		warning("MOE distribution with tilt parameter=1 is identical to original 				distribution")
	nottiltfix.arg<-fix.arg[-tiltfix] #exclude tilt parameter from fix.arg
	}
	else
	nottiltfix.arg<-fix.arg
vfix <- unlist(fix.arg) # including tilt parameter
vnottiltfix <- unlist(nottiltfix.arg)# not including tilt parameter
# check unexpected names
mfix <- match(names(vnottiltfix ), argdistname)
	if (any(is.na(mfix))) 
      	stop(errtxt$t4)
# check that some parameters are not both in fix.arg and start.arg
minter <- match(names(vstart), names(vfix))
	if (any(!is.na(minter)))
	stop(errtxt$t5)
#check NA/NaN values
	if(any(is.na(vfix)) || any(is.nan(vfix)))
 	stop(errtxt$t7)
allparname <- names(c(vstart, vfix))# including tilt parameter
nottiltallparam<-names(c(vnottiltstart , vnottiltfix)) # not including tilt parameter
}
else
{
allparname <- names(vstart) # including tilt parameter
nottiltallparam<-names(vnottiltstart) # not including tilt parameter
}
theoparam <- MOEcomputegetparam(argdistname)
#special case where both scale and rate are allowed, see ?dgamma
if("scale" %in% theoparam && "rate" %in% theoparam)
errt8 <- any(!nottiltallparam %in% theoparam) || length(nottiltallparam) != length(theoparam)-1
else
errt8 <- any(!theoparam %in% nottiltallparam)
#only make a warning if unset arguments have a default value
if(errt8)
{
unsetarg <- theoparam[!theoparam %in% nottiltallparam] 
	if(any(hasnodefaultval[unsetarg]))
	stop(errtxt$t8, " ",unsetarg, ".")
	else
	warning(errtxt$t9, " ",unsetarg, ".")
}
list("start.arg"=start.arg, "fix.arg"=fix.arg)
}
#--------------4- MOEcomputegetparam --------------------------------
# typical to computegetparam
# INPUTS 
# argdistname : argument names of the distribution from names(formals())
# OUTPUTS 
# parameter names (as a vector) of the distribution (excluding non parameter argument)
MOEcomputegetparam <- function(argdistname)
{
  #remove first argument, that should be "x", "p", "q", or "n", see ?dgamma, pgamma, qgamma
  argdistname <- argdistname[-1]
  nonparaminR <- c("x", "p", "q", "n") #defensive programming
  #remove other arguments, see ?dgamma, pgamma, qgamma, dbeta
  nonparaminR <- c(nonparaminR, "log", "log.p", "lower.tail", "ncp")
  nonparaminActuar <- c("limit", "order", "t")
  nonparaminGamlssdist <- "fast"
  nonparamspecial <- c("...", "..1", "..2")
  #see ?dnig, dhyperb, dskewlap, dgig,...
  nonparaminGenHyperbolic <- c("param", "KOmega", "ibfTol", "nmax", "method", "intTol","valueOnly", "nInterpol", "uniTol", "subdivisions", "logPars")
  #see ?dsn
  nonparamsn <- "dp"
  plist <- setdiff(argdistname, nonparaminR)
  plist <- setdiff(plist, nonparaminActuar)
  plist <- setdiff(plist, nonparaminGamlssdist)
  plist <- setdiff(plist, nonparamspecial)
  plist <- setdiff(plist, nonparaminGenHyperbolic)
  plist <- setdiff(plist, nonparamsn)
  plist
}
#--------------5- memp --------------------------------
memp<-function(x,order)
{
mean(x^order)
}
#--------------6- memp2 --------------------------------
memp2<-function(x,order,weights)
{
sum(x^order * weights)/sum(weights)
}
#--------------7- mMOE --------------------------------
mMOE<-function(order,par,fix.arg,distr)
{#mMOE
	if (!is.list(par))
	par<-as.list(par)
	if (is.null(names(par)))
	stop("'par' must be a named list")
	if (!is.null(fix.arg))
	{#if (!is.null(fix.arg))
		if (!is.list(fix.arg))
		fix.arg<-as.list(fix.arg)
		if (is.null(names(fix.arg)))
		stop("'fix.arg', if exist, it must be a named list")
	}#if (!is.null(fix.arg))
#exclude tilt from par or fixed.arg
tiltpar<-match("tilt",names(par))
	if(!is.na(tiltpar))
	{#if(!is.na(tiltpar))
	alpha<-par$tilt
	par<-par[-tiltpar] #exclude tilt parameter from par
	}#if(!is.na(tiltpar))
tiltfix<-match("tilt",names(fix.arg))
	if(!is.na(tiltfix))
	{#if(!is.na(tiltfix))
	alpha<-fix.arg$tilt
	fix.arg<-fix.arg[-tiltfix] #exclude tilt parameter from par
	}#if(!is.na(tiltfix))
ddistname <- paste("d", distr, sep = "")
pdistname <- paste("p", distr, sep = "")
qdistname <- paste("q", distr, sep = "")
	if (!exists(ddistname, mode = "function")) 
	paste("The " , distr, " disribution", "is not defined")
alpha_bar<-1-alpha
parnm <- names(par)
args <- names(formals(ddistname))
m <- match(parnm,args)
	if (any(is.na(m))) 
	stop("'par' specifies names which are not arguments to ",distr)
	integrand <- function(z,alpha,qdistname,par,order)
	{#integrand 
	alpha_bar<-1-alpha
	x<-alpha*(1-z)/(alpha+alpha_bar*z)
	L1<-do.call(qdistname, c(list(x), as.list(par)))
	L1<-L1^order
	}#integrand 
a<-integrate(integrand, lower = 0, upper =1,alpha=alpha,qdistname=qdistname ,par=par,order=order)
return(a$value)
}#mMOE
