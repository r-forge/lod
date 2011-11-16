# TODO: Add comment
# 
# Author: adehuson
###############################################################################
##Model C log likelihood equation
modelCLik<-function(mu=c(0.1, 1), strainX) { 	
	## "mu":starting values for the parameters, ALWAYS positive! 
	## "strainX": data
	
	delta <- 1e-5		#small number to help avoid log(0) error
	beta <- mu[1]			
	theta <- mu[2]
        if (beta <=0 || theta <=0)
          return(-1e10)                 # restrict beta and theta > 0
        else {
          strainX$DNA <- strainX$DNA + delta

          ## bernoulli prob.
          pp <- (1 - exp(- beta * strainX$DNA)) ^ theta + delta
          
          ##log liklihood equation
          log.lik<- sum(strainX$positive * log(pp) + 
                        (strainX$trials - strainX$positive) *
                        log(1 + 2*delta - pp))
          
          return(log.lik)
        }
}


  
