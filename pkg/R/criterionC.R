# TODO: Add comment
# 
# Author: adehuson

##find AIC, BIC, AICc of Model C
criterionC<-function(strainX, type=c("AIC", "BIC", "AICc")){
	## "strainX":data
	## "type": select criterion

	##Prune redundant points at the upper and lower bound
	strainX2<-pruneFunc(strainX)
	
	##run optim to find estimates of parameters
	optimC <- try(optim(c(0.693/median(strainX2$DNA), 1),
                            modelCLik, gr=NULL, strainX2,
                            control=list(fnscale=-1), # maximization, NOT minimization! 
                            method = "Nelder-Mead", hessian=TRUE),
                      silent=T)
        ## BFGS requires gradient calculation, which cannot be done on our function.

        if (inherits(optimC, "try-error"))
          return(list(crit=NA))
          
	par.est <- optimC$par #parameter estimates found by optim function
	beta<-par.est[1]		
	theta<-par.est[2]
	
	fisher <- - optimC$hessian              	# Fisher Information Matrix
	var.cov <- try(solve(fisher + diag(c(1e-5, 1e-5))),
                             silent=T)  # guard against singularity
                                        # Variance-Covariance Matrix
	
	##log likelihood equation using new parameter estimates
	L<-modelCLik(mu=c(beta, theta), strainX2)	
	
	##select between AIC, BIC, or AICc
	type <- match.arg(type)						
	switch(type,			
			AIC = {  
				crit<-4-2*L
			},
			BIC = { 
				crit<--2*L+2*log(nrow(strainX))							
			},
			AICc = { 
				crit<-(4-2*L)+12/(nrow(strainX)-3) 
			}
	)
	## placeholder for LODFunc
	r<-NULL
	
	#return criterion, parameter estimates, and vcov matrix
	return(list(crit=crit, par.est=par.est, var.cov=var.cov, r=r))		
}


