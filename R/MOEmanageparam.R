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
