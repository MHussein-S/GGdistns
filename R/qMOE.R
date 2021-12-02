#q-th quantile of MOE for specific distribution (done)
qMOE<-function (q, distr,par)  
{
	if(is.null(names(par)) || !is.list(par))                              
	stop("Parameters must be a named list.")       
tiltpar<-match("tilt",names(par))
if (is.na(tiltpar))
stop("Tilt parameter not defined")

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
q1<-alpha*q/(1-alpha_bar*q)
x<-do.call(qf, c(list(q1), as.list(par)))
}
