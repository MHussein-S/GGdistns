#all function to estimate parameters of MOE "Maximum Liklehood"
#-------input starting value >> under construction ---

MOEmledist<-function (data, distr, start.arg = NULL, fix.arg = NULL, optim.method = "default",
    lower = -Inf, upper = Inf, silent = TRUE, gradient = NULL,  ...)
{
ddistname <- paste("d", distr, sep = "")
argddistname <- names(formals(ddistname))
if (!exists(ddistname, mode = "function"))
stop("The ", ddistname, " function must be defined")

if (is.vector(start.arg))
start.arg <- as.list(start.arg)

if (!(is.vector(data)&is.numeric(data) & length(data) > 1))
stop("data must be a numeric vector of length greater than 1")
arg_startfix <- MOEmanageparam(start.arg = start.arg, fix.arg = fix.arg,obs = data, distname= distr)
hasnodefaultval <- sapply(formals(ddistname), is.name) # True for arguments does not have default value and False for other arguments (ex. sapply(formals(dgamma),is.name))
arg_startfix <- MOEcheckparamlist(arg_startfix$start.arg,arg_startfix$fix.arg, argddistname, hasnodefaultval)
if (is.function(fix.arg))
fix.arg.fun <- fix.arg
else fix.arg.fun <- NULL
vstart <- unlist(arg_startfix$start.arg)
fix.arg <- arg_startfix$fix.arg
# ----- Minus log liklehood function of MOE distribution ----
objfn <- function(par, fix.arg, obs, distr)
{
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
LL <- -n*log(alpha)- sum(log(fbase))+2*sum(log(1-alpha_bar*(1-Fbase)))
return(LL)
}
# -------------------------------------------------------------------------------------
hasbound <- any(is.finite(lower) | is.finite(upper))
if (optim.method == "default")
meth <- ifelse(length(vstart) > 1, "Nelder-Mead", "BFGS")
else
meth <- optim.method

if (meth == "BFGS" && hasbound && is.null(gradient))
{
meth <- "L-BFGS-B"
txt1 <- "The BFGS method cannot be used with bounds without provided the gradient."
txt2 <- "The method is changed to L-BFGS-B."
warning(txt1, " ",txt2)
}

if (hasbound)
{
	if (!is.null(gradient))
	opt.fun <- "constrOptim"
	else
	{
                   	if (meth == "Nelder-Mead")
		opt.fun <- "constrOptim"
             		else if (meth %in% c("L-BFGS-B", "Brent"))
             		opt.fun <- "optim"
              		else
		{
		txt1 <- paste("The method", meth, "cannot be used by constrOptim() nor 			optim() without gradient and bounds.")
          		txt2 <- "Only optimization methods L-BFGS-B, Brent and Nelder-Mead can be 			used in 	such case."
		stop(paste(txt1, txt2))
                 	}
	}

	if (opt.fun == "constrOptim")
	{
	par <- length(vstart)
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

	#	if (!cens)
      opttryerror <- try(opt <- constrOptim(theta = vstart,f = objfn, ui = Mat, ci = Bnd, grad = gradient, fix.arg = fix.arg, obs = data, distr = distr,hessian = !is.null(gradient), method = meth, ...), silent = TRUE)
      # 	else
	#	opttryerror <- try(opt <- constrOptim(theta = vstart,  f = fnobjcens, ui = Mat, 		#	ci = Bnd, grad = gradient, ddistnam = ddistname, rcens = rcens, lcens = lcens,                   	#	icens = icens, ncens = ncens, pdistnam = pdistname, fix.arg = fix.arg, hessian = !		#	is.null(gradient),  method = meth, ...), silent = TRUE)

		if (!inherits(opttryerror, "try-error"))
			if (length(opt$counts) == 1)
                    		opt$counts <- c(opt$counts, NA)

	}
	else
	{
	#	if (!cens)
opttryerror <- try(opt <- optim(par = vstart, fn = objfn, fix.arg = fix.arg, obs = data, gr = gradient, distr = distr, hessian = TRUE, method = meth,lower = lower, 	upper = upper, ...), silent = TRUE)
              	#	else
	#	opttryerror <- try(opt <- optim(par = vstart, fn = fnobjcens, fix.arg = fix.arg, gr =gradient, cens = rcens, lcens = lcens, icens = icens, ncens = ncens, ddistnam = ddistname,  distnam = pdistname, hessian = TRUE, method = meth, 	#	lower = lower, upper = 	upper, ...), silent = TRUE)
	}
}
else
{
opt.fun <- "optim"
	#if (!cens)
opttryerror <- try(opt <- optim(par = vstart, fn = objfn, fix.arg = fix.arg, obs = data, gr = 	gradient, distr = distr, hessian = TRUE, method = meth, lower = lower, 	upper = upper,  ...), silent = TRUE)

           	#else
	#opttryerror <- try(opt <- optim(par = vstart, fn = fnobjcens, fix.arg = fix.arg, gr = gradient,rcens = rcens, lcens = lcens, icens = icens, ncens = ncens, ddistnam = 	ddistname, distnam= pdistname, hessian = TRUE, method = meth, lower = lower, 	#upper = upper, ...), 	silent = TRUE)
}

if (inherits(opttryerror, "try-error"))
{
warnings("The function optim encountered an error and stopped.")
	if (getOption("show.error.messages"))
	print(attr(opttryerror, "condition"))
return(list(estimate = rep(NA, length(vstart)), convergence = 100, loglik = NA, hessian 	=NA, optim.function = opt.fun,fix.arg = fix.arg, optim.method = meth, fix.arg.fun =fix.arg.fun,  counts = c(NA, NA)))
}

if (opt$convergence > 0)
warnings("The function optim failed to converge, with the error code ", opt$convergence)

if (is.null(names(opt$par)))
names(opt$par) <- names(vstart)
loglik <- -opt$value
k<-length(opt$par)
n<-length(data)
AIC<--2*loglik+2*k
BIC<- -2*loglik+ k*log(n)
HQIC<- -2*loglik+2*k*log(log(n))

res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, hessian = opt$hessian, optim.function = opt.fun, optim.method = meth, fix.arg = fix.arg, fix.arg.fun =fix.arg.fun, counts = opt$counts, optim.message = opt$message, loglik = -opt$value,
 AIC=AIC, BIC=BIC, HQIC=HQIC)
return(res)
}
#---------------1-MOEstart------------
#------- new version of MOEstart ---
MOEstart<- function(x, distr)
{
pdistname <- paste("p", distr, sep = "")
if (!exists(pdistname, mode = "function"))
  stop(paste("The ", ddistname, " function must be defined"))
n<-length(x)
#alpha0<-2    #made as comment recently
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

