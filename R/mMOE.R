mMOE<-function(order=1, distr, par)
{

if (missing(distr) | missing(par))   
	stop("distr and par  must defined")
if (!is.list(par)||is.null(names(par))) 
	stop("'par' must be a named list")
ddistname <- paste("d", distr, sep = "")
pdistname <- paste("p", distr, sep = "")
qdistname <- paste("q", distr, sep = "")
if (!exists(ddistname, mode = "function")) 
paste("The " , distr, " disribution", "is not defined")

tiltpar<-match("tilt",names(par))
if (is.na(tiltpar))
stop("tilt parameter not defined")
else
alpha<-par$tilt
par<-par[-tiltpar]
alpha_bar=1-alpha
parnm <- names(par)
args <- names(formals(ddistname))
m <- match(parnm,args)
if (any(is.na(m))) 
	  stop("'par' specifies names which are not arguments to ",distr)

integrand <- function(z,alpha,qdistname,par,order)
{
alpha_bar<-1-alpha
x<-alpha*(1-z)/(alpha+alpha_bar*z)
do.call(qdistname, c(list(x), as.list(par)))^order
}
a<-integrate(integrand, lower = 0, upper = 1,alpha=alpha,qdistname=qdistname ,par=par,order=order)
return(a$value)
}
