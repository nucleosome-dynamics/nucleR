#One shot version
pcKeepCompDetect<-function(data, pc.min=0.01, pc.max=0.1, max.iter=20, verbose=FALSE,
               cor.target=0.999, cor.tol=1e-4, smpl.num=25, smpl.min.size=2^10, smpl.max.size=2^14)
{
	#Sample data
	samp = .sampleData(data, smpl.num, smpl.min.size, smpl.max.size, verbose=verbose)

  #Calculates the correlation between fft and normal for a given pcKeepComp
  .meancor<-function(pc)
  {
    if(pc==1) return(1)
		if(pc==0) return(0)
    cors = sapply(samp, function(x) cor(filterFFT(x, pcKeepComp=pc), x))
    return(mean(cors, na.rm=TRUE))
  }

	#Initial evaluation	
	eval.min = .meancor(pc.min) 
	eval.max = .meancor(pc.max)
	desv.min = eval.min - cor.target
	desv.max = eval.max - cor.target

	if(abs(desv.min) < abs(desv.max))
	{
		pc.best = pc.min
		dv.best = desv.min
	}else{
		pc.best = pc.max
		dv.best = desv.max
	}

	#Iterate
	i = 0
	while((abs(dv.best) > cor.tol) & (i < max.iter))
	{
		i = i + 1
		if(verbose) cat("Iteration", i, ".  pcKeepComp=", pc.best, ".  Correlation=", dv.best+cor.target, "\n")

		#Lineal interpolation of optimal value
		pc.new = pc.min-(desv.min*((pc.max-pc.min)/(desv.max-desv.min)))
		if(pc.new > 1) pc.new = 1 else if(pc.new < 0) pc.new = 0
		
		#New evaluation
		eval.new = .meancor(pc.new)
		desv.new = eval.new - cor.target

		#Update limits	
		if(desv.min<0 & desv.max>0)
		{
			if(desv.new > 0){ pc.max = pc.new ; desv.max = desv.new }
			else{ pc.min = pc.new ; desv.min = desv.new }

		}else if(desv.min>0 & desv.max>0)
		{
			pc.max = pc.min; desv.max = desv.min
			pc.min = pc.new; desv.min = desv.new

		}else if(desv.min<0 & desv.max<0)
		{
			pc.min = pc.max; desv.min = desv.new
			pc.max = pc.new; desv.max = desv.new		
		}

		#New best value update
		if(abs(desv.new) < abs(dv.best)){ pc.best = pc.new; dv.best = desv.new}
	}

	return(pc.best)

}

.sampleData <- function(data, smpl.num, smpl.min.size, smpl.max.size, verbose=FALSE)
{
	#Check coherency
	if(smpl.min.size > smpl.max.size)
	{
		warning("smpl.min.size > smpl.max.size, using only smpl.min.size")
		smpl.max.size = smpl.min.size
	}

	res = list()

	#For short sequences, use all the data
  if(length(data) < smpl.min.size*(smpl.num/2))
  {
		if(verbose) message("Short fragment. Using all data")
		res[[1]] = data

	#For long sequence, use sampling
	}else{

		if(verbose) message("Long sequence. Trying sampling...")

		#Select ranges <> from 0 or NA and longer than minimum size
		rang = IRanges(data != 0 & !is.na(data))
		rang = rang[width(rang) > smpl.min.size]
		
		tota = sum(width(rang))
		#If the overall useful bases don't satisfy the criteria, return all
		if(tota < smpl.min.size*(smpl.num/2))
		{
			if(verbose) message(" No enough covered bases to apply sampling. Using all data") 
			res[[1]] = data
		}
		else
		{
			if(verbose) message(" Selecting regions for sampling")
			#Sampling with weighted probabilities according range's width
			reps = sample(1:length(rang), size=smpl.num, replace=TRUE, prob=width(rang) / tota)
			
			#Random selection of start points
			# In short, if a range has a width higher than smpl.max.size we will pickup a random
			# start position inside it that still it's in the limits
			# The width of the final ranges will be the minimum of the range's size and smpl.max.size
			wids = width(rang[reps])
			marg = unlist(sapply(wids - smpl.max.size, max, 0))
			rnof = unlist(sapply(marg[marg > 0], function(x) floor(runif(n=1, max=x, min=1))))
			marg[marg > 0] = rnof
			
			fina = IRanges(start=start(rang[reps])+marg, width=sapply(wids, min, smpl.max.size))

			res = lapply(fina, function(x) data[x])
		}
	}
	if(verbose) message("Returning ", length(res), " regions")
	return(res)
}

