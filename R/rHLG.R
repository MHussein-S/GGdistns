#n variate from Half logistic G for any base distribution
rHLG<-function (n, distr,par)
{
  if(!is.list(par))
    par<- as.list(par)
  if (!is.character(distr))
    stop("distr must be a character string naming a distribution")
  if(is.null(names(par)))
  	stop("Parameters must be a named list or vector.")
  apar<-match("a",names(par))
  if (is.na(apar))
    stop("parameter \"a\" is not defined")
  a<-par$a
  par1<-par[-apar]
#--------
  pdistname <- paste("p",distr,sep="")
  if (!exists(pdistname, mode="function"))
    stop(paste("The ", distr, " distribution is not defined"))
  parnm <- names(par1)
  args <- names(formals(pdistname))
  m <- match(parnm,args)
  if (any(is.na(m)))
    stop("'par' specifies names which are not arguments to ",distr)
#--------
  u<-runif(n)
  r<- qHLG(u, distr,par)
}
