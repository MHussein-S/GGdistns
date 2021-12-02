# ------ plots SF, CDF, PDF, HRF and PRHRF form x = xmin to x = xmax} --- working
plotMOE <- function (distr,para,xmin,xmax)
{
if (missing(distr) || missing(para))
      stop("distr and parameters must defined")

if(is.null(names(para)) || !is.list(para))
	stop("Parameters must be a named list.")



tiltpar<-match("tilt",names(para))

	if (is.na(tiltpar))
	stop("Tilt parameter is not defined")

apdistr <- switch(distr, 'uniform' = "unif", 'exponential' = "exp", 'normal' = "norm", 'chi-squared' = "chisq", 'chisquared' = "chisq",'log-normal' = "lnorm", 'lognormal' = "lnorm", 'logistic' = "logis", 't' = "mydt")


	if(!is.null(apdistr))
	distr<-apdistr


	if (distr %in% c("exp", "chisq","weibull","gamma")) #need to check all distributions
	{
	if (xmin < 0)
	stop("negative domain not supported for ",distr)
	}

	if (distr %in% c("beta"))
	{
	if (xmin < 0 | xmax>1)
	stop("(",xmin,",",xmax,")", "is not supported domain for ",distr)
	}

alpha<-para$tilt
para<-para[-tiltpar]

pdistname <- paste("p", distr, sep = "")
ddistname <- paste("d", distr, sep = "")
if (!exists(pdistname, mode = "function"))
         stop(paste("The " , distr, " disribution", "is not defined"))

parnm <- names(para)
args <- names(formals(ddistname))
m <- match(parnm,args)
if (any(is.na(m)))
	  stop("'para' specifies names which are not arguments to ",distr)
#theoparam<-MOEcomputegetparam(args)
#print(theoparam)
#hasnodefaultval <- sapply(formals(ddistname), is.name)# True for arguments does not have default value and False for other arguments (ex. sapply(formals(dgamma),is.name))
#special case where both scale and rate are allowed, see ?dgamma
#	if("scale" %in% theoparam && "rate" %in% theoparam)
#	err <- any(!parnm %in% theoparam) || length(parnm) != length(theoparam)-1
#	else
#	err <- any(!theoparam %in% parnm)
#print(err)
#only make a warning if unset arguments have a default value
#	if(err)
#	{
#	unsetarg <- theoparam[!theoparam %in% parnm]
#print("unsetarg: ",unsetarg)

#		if(any(hasnodefaultval[unsetarg]))
#		stop("Some parameter names have no value: ",unsetarg, ".")
#		else
#		warning("Some parameter names have no value but have a default value:" ,unsetarg, ".")
#	}
x <- seq(xmin, xmax, by = (xmax - xmin)/1000)
ycdf <- do.call(pdistname, c(list(x), as.list(para)))
ypdf <- do.call(ddistname, c(list(x), as.list(para)))
alpha_hat<-1-alpha
MOCDF<-ycdf/(1-alpha_hat*(1-ycdf))
MOSF<-1-MOCDF
MOPDF<-(alpha*ypdf)/(1-alpha_hat*(1-ycdf))^2
if(distr=="exp" && alpha==1)
{
MOHRF<-rep(para,length(x))
}
else
{
MOHRF<-ypdf/(1-ycdf)
MOHRF<-MOHRF/(1-alpha_hat*(1-ycdf))
}
MOPRHRF<-alpha*ypdf/(ycdf*(1-alpha_hat*(1-ycdf)))
par(mfrow = c(3, 2))
plot(x, MOSF,ylab="SF",col = "red",main = paste("Theo. MO SF for ",distr ,"basline distr"),type='l')
plot(x, MOCDF,ylab="CDF",col = "red",main = paste("Theo. MO CDF for ",distr ,"basline distr"),type='l')
plot(x, MOPDF,ylab="PDF",col = "red",main = paste("Theo. MO PDF for ",distr ,"basline distr"),type='l')
plot(x, MOHRF,ylab="HRF",col = "red",main = paste("Theo. MO HRF for ",distr ,"basline distr"),type='l')
plot(x, MOPRHRF,ylab="PRHRF",col = "red",main = paste("Theo. MO PRHRF for ",distr ,"basline distr"),type='l')

}
#--------------4- MOEcheckparamlist --------------------------------
# typical to computegetparam
# INPUTS
# argdistname : argument names of the distribution from names(formals())

# OUTPUTS
# parameter names (as a vector) of the distribution (excluding non parameter argument)
#MOEcomputegetparam <- function(argdistname)
#{
  #remove first argument, that should be "x", "p", "q", or "n", see ?dgamma, pgamma, qgamma
  #argdistname <- argdistname[-1]
  #nonparaminR <- c("x", "p", "q", "n") #defensive programming
  #remove other arguments, see ?dgamma, pgamma, qgamma, dbeta
  #nonparaminR <- c(nonparaminR, "log", "log.p", "lower.tail", "ncp")
  #nonparaminActuar <- c("limit", "order", "t")
  #nonparaminGamlssdist <- "fast"
  #nonparamspecial <- c("...", "..1", "..2")
  #see ?dnig, dhyperb, dskewlap, dgig,...
  #nonparaminGenHyperbolic <- c("param", "KOmega", "ibfTol", "nmax", "method", "intTol","valueOnly", "nInterpol", "uniTol", "subdivisions", "logPars")
  #see ?dsn
  #nonparamsn <- "dp"
  #plist <- setdiff(argdistname, nonparaminR)
  #plist <- setdiff(plist, nonparaminActuar)
  #plist <- setdiff(plist, nonparaminGamlssdist)
  #plist <- setdiff(plist, nonparamspecial)
  #plist <- setdiff(plist, nonparaminGenHyperbolic)
  #plist <- setdiff(plist, nonparamsn)
  #plist
#}


