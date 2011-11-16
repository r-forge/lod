# TODO: Add comment
# 
# Author: adehuson
###############################################################################
pruneFunc<-function(strainX){
	## "strainX":data
	
	strainX2<-strainX
	cnt.between <- sum(strainX2$proportion > 0 & strainX2$proportion < 1)

	if (cnt.between < 2) {
			suff <- "Inadequate"
		}else if (cnt.between < 3) {
			suff <- "Weak"
		}else suff <- "Adequate"


	## Data Editing to Remove Singularity when Data is Inadequate
	if (suff == "Inadequate") {
		## Prune redundant 0% or 100% proportions and retain at most one extreme case.
		## Assume the data is sorted by DNA in ascending order!
		prop.start.idx <- max(1, match(TRUE, strainX2$proportion > 0) - 1)
                prop.end.idx <- min(nrow(strainX2),
                                    match(TRUE, strainX2$proportion == 1),
                                    na.rm=T)

		strainX2<- strainX2[prop.start.idx:prop.end.idx,] # pruning extreme cases
	}

	return(strainX2=strainX2)
}
