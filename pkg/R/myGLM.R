## A helper function for calling different models for binary model fitting
myGLM<- function(formula, link.fun, dat)  {
  ## "link.fun": Link Function for GLM modeling. 
  ## "dat": data, assume 1st column is
  ##           covariate (analyte concentration), 2nd col is proportions
  ##           (positives / trials), 3rd col is number of trials, (optional) 4th
  ##           col is group ID's (multiple sets of experiments). 
  ##           The names of these columns are:
  ##           ("DNA", "proportion", "trials", "groups")

  require(MASS)
  
  ## Use GLM model. If "groups" exist, ignore it. 
  r <-try(glm(formula, data=dat, weights=trials,
              family=binomial(link=link.fun)), silent=T)
  error <-inherits(r, "try-error")

  ## Fixes possible error by running r again with pruned data
  if (error){
    ##Prune redundant points at the upper and lower bound
    dat2<-pruneFunc(dat)
    
    r <-try(glm(formula, data=dat2, weights=trials,
                family=binomial(link=link.fun)), silent=T)
    error <- inherits(r, "try-error")
  }

  rg <- NULL
  if (ncol(dat)==4)  {
    ## Has Group ID's, use GLMM model via "glmmPQL" from MASS library.
    require(MASS) 
    rg <-try(glmmPQL(formula, random = ~ 1 | groups, 
                     data=dat, weights=trials,
                     family=binomial(link=link.fun), verbose=F), silent=T)
    ## fail if any model (GLM or GLMM) fails
    error <- error && inherits(rg, "try-error")
  }

  ##return r for criterion calculations
  if (error) ans <- list(error=TRUE)
  else ans <- list(r=r, rg=rg, error=FALSE)

  ans
}

if (F) {
  strainX <- structure(list(DNA = c(0, 0, 0, 2.5, 2.5, 2.5, 5, 5, 5, 10, 10, 
10, 15, 15, 15, 20, 20, 20), proportion = c(0, 0, 0.1, 0.7, 0.2, 
0.3, 0.9, 0.5, 0.4, 0.9, 0.7, 0.6, 1, 0.9, 0.8, 1, 0.9, 1), trials = c(10, 
10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 
10), groups = c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 
1L, 2L, 3L, 1L, 2L, 3L)), .Names = c("DNA", "proportion", "trials", 
"groups"), row.names = c(1L, 7L, 13L, 2L, 8L, 14L, 3L, 9L, 15L, 
4L, 10L, 16L, 5L, 11L, 17L, 6L, 12L, 18L), class = "data.frame")

  r <- myGLM(proportion ~ DNA, "logit", strainX)
 
}




