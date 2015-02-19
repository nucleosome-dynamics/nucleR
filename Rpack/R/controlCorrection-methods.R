#############
setMethod("controlCorrection", signature(exp="SimpleRleList"), function(exp, ctr, mc.cores=1)  {

	if(class(exp) != class(ctr)) stop("'exp' and 'ctr' classes must be equal") 
	if(length(which(!(names(exp) %in% names(ctr)))) != 0 | 
		 length(which(!(names(ctr) %in% names(exp)))) != 0 ) stop("names should be equal in both datasets")

	#Check cores
	mc.cores = .check.mc(mc.cores)

  if(mc.cores > 1)
  {
		expt = lapply(exp, identity)
		ctrt = lapply(ctr, identity)
		res = mclapply(names(exp), function(chr) controlCorrection(expt[[chr]], ctrt[[chr]]), mc.cores=mc.cores)
  }
  else
  {
  	res = lapply(names(exp), function(chr) controlCorrection(exp[[chr]], ctr[[chr]]))
	}

	names(res) = names(exp)
	return(RleList(res, compress=FALSE))
})

################
setMethod("controlCorrection", signature(exp="Rle"), function(exp, ctr)  {

  if(class(exp) != class(ctr)) stop("'exp' and 'ctr' classes must be equal")

	return(Rle(controlCorrection(as.vector(exp), as.vector(ctr))))
})

################
setMethod("controlCorrection", signature(exp="list"), function(exp, ctr, mc.cores=1)  {

  if(class(exp) != class(ctr)) stop("'exp' and 'ctr' classes must be equal")

  if(length(which(!(names(exp) %in% names(ctr)))) != 0 |  
     length(which(!(names(ctr) %in% names(exp)))) != 0 ) stop("names should be equal in both datasets")

  #Check cores
  mc.cores = .check.mc(mc.cores)

  if(mc.cores > 1)
  {
    res = mclapply(names(exp), function(chr) controlCorrection(exp[[chr]], ctr[[chr]]), mc.cores=mc.cores)
  }
  else
  {
    res = lapply(names(exp), function(chr) controlCorrection(exp[[chr]], ctr[[chr]]))
  }
    
  names(res) = names(exp)
  return(res)
})

################
setMethod("controlCorrection", signature(exp="numeric"), function(exp, ctr)  {

  if(class(exp) != class(ctr)) stop("'exp' and 'ctr' classes must be equal")

	res = exp - (ctr - mean(ctr[ctr != 0]))

	#If value is under 0, set to 0
	res[res  < 0] = 0
	return(res)
})
