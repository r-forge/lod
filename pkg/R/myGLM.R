

## A helper function for calling different models for binary model fitting
myGLM<- function(formula, link.fun, strainX)  {
	## "link.fun": Link Function for GLM modeling. 
	## "strainX: data

	## RESPONSE: data$prop. 

	r <-try(glm(formula, data=strainX,weights=strainX$trials,family=binomial(link=link.fun)), silent=T)
        error<-inherits(r, "try-error")
	
	##fixes possible error by running r again with pruned data
	if (error){
		##Prune redundant points at the upper and lower bound
		strainX2<-pruneFunc(strainX)
		
		r <-try(glm(formula, data=strainX2,weights=strainX$trials,family=binomial(link=link.fun)), silent=T)
			error<-inherits(r, "try-error")
	}

	##return r for criterion calculations
        if (error) ans <- list(error=TRUE)
        else ans <- list(par.est=coef(r), var.cov=vcov(r), r=r, error=FALSE)

        ans
}




