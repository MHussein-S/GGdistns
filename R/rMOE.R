#random number generator for MOE for specific distribution (done)
rMOE<-function (n, distr,par)  
{
	if(is.null(names(par)) || !is.list(par))                              
	stop("Parameters must be a named list.")       
tiltpar<-match("tilt",names(par))
if (is.na(tiltpar))
stop("Tilt parameter not defined")

y<-runif(n)

alpha<-par$tilt
qdistname <- paste("q", distr, sep = "")

if (!exists(qdistname, mode = "function")) 
stop(paste("The " , distr, " disribution", "is not defined"))

qf <- get(qdistname, mode = "function")
par<-par[-tiltpar]
nm <- names(par)
args <- names(formals(qf))
m <- match(nm,args)
if (any(is.na(m))) 
	  stop(paste("'par' specifies names which are not arguments to ",distr))
alpha_bar=1-alpha
y<-alpha*y/(1-alpha_bar*y)
x<-do.call(qf, c(list(y), as.list(par)))
}
