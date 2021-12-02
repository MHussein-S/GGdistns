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
		war("MOE distribution with tilt parameter=1 is identical to original 				distribution")
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