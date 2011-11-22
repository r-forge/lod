                                        # TODO: Add comment
                                        # 
                                        # Author: adehuson
###############################################################################

LODFunc<- function(model, obj, lod.fun, lod.level, strainX, max.range) {
  ## "model": name of model used
  ## "obj": a GLM model object or Model C object. 
  ## "lod.fun" is the FUNCTION (not the evaluation of the function)
  ##           for 95% LOD calculation.
  ##           If it is not missing, do the Confidence Interval calculation. 
  ## "lod.level":level for LOD(95%LOD, 5%LOD, etc.)should be input as a list
  ## "strainX": data
  ## "max.range": acceptable range for the 95% CI
  ## Return message about adequacy of the data

  ## Extract "var.cov" (variance covariance matrix) and
  ##         "par.est" (parameter etimates)

  ## Class Description as a single string
  class.desc <- paste(class(obj), collapse=" ")
  require(MASS)
  
  glmm <- FALSE
  if (grepl("glmmPQL", class.desc)) {   # GLMM 
    par.est <- fixef(obj)
    var.cov <- vcov(obj)
    glmm <- TRUE
  } else if (grepl("glm", class.desc)) { # GLM
    par.est <- coef(obj)
    var.cov <- vcov(obj)
  } else {                              # for Model C 
    par.est <- obj$par.est
    var.cov <- obj$var.cov
  }    
  
  cnt.between <- sum(strainX$proportion > 0 & strainX$proportion < 1)
  
  if (cnt.between < 2) {
    suff <- "Inadequate"
  } else if (cnt.between < 3) {
    suff <- "Weak"
  } else suff <- "Adequate"
  
  ## output plot
  plot(proportion ~ DNA, data = strainX, 
       main=paste(model, ": Data is", suff, "!"))

  grid <- seq(0,  max(strainX$DNA), length.out=100)

  if (glmm) {
    ## Plot individual line for each random intercept level
    for (g in unique(strainX$groups)) {
      nd <- data.frame(DNA=1e-5 + grid, groups=g)
      lines(grid, predict(obj, newdata=nd, type="response"), lty=2)
    }
  }    
    
  ##prediction lines for glm equations
  if(class.desc != "list") {
    lines(grid, predict(obj, newdata=data.frame(DNA=1e-5 + grid),
                        type="response", level=0))
  }  else {   ## prediction lines for Model C 
    beta<-par.est[1]		
    theta<-par.est[2]
    lines(grid, (1-exp(-beta*grid))^theta)
  }

  ##loop through all LOD levels
  for(i in seq(along=lod.level)){
    ce.arg <- as.list(par.est) 								# used for assigning function default parameter values
    names(ce.arg) <- par.names <- names(formals(lod.fun)) 	# get parameter names
    
    ## Calc Symbolic Derivatives for the function "lod.deriv".
    lod.deriv <-
      deriv(parse(text=deparse(lod.fun)[3]), 			# get LOD expression.
            par.names,       
            function.arg=T) # return a function with derivative attribute

    pp <- lod.level[[i]]
    ## for Probit model, use "Phi^(-1)(pp). 
    pp <- ifelse(model=="Probit",  qnorm(pp), pp)
    formals(lod.deriv) <- c(ce.arg, list(p=pp)) # assign parameter estimates as default values
    
    ## 95% LOD estimate from model parameters, its attr has derivatives.
    lod.est <- lod.deriv()									 # using values that are assigned from estimates!
    
    ## evaluate 1st derivative of lod.deriv with respect to ce (parameter)
    lod.der.ce <- attr(lod.est, "gradient")
    
    ## calc SD of 95% LOD using Delta Method
    if (inherits(var.cov, "try-error")) {
      ci.lower <- -Inf
      ci.upper <- Inf
    } else {
      lod.sd <- sqrt(lod.der.ce %*% var.cov %*% t(lod.der.ce))
      
      ## 95% CI for 95% LOD is ...
      ci.lower <- round(lod.est - 1.96 * lod.sd, 4)
      ci.upper <- round(lod.est + 1.96 * lod.sd, 4)
    }
    
    ##Check for an extremely large CI
    ci.range <- ci.upper - ci.lower
    
    range.error <- (ci.range < max.range)

    ## output LOD and 95% CI into legend based on validity of CI
    if(range.error) {
      leg.text1 <- c(paste(100*lod.level[[i]],"% LOD =", round(lod.est, 1),", ",
                           100*lod.level[[i]],"% LOD CI = (", round(ci.lower, 1), ",", round(ci.upper, 1), ")"))
      segments(-1, lod.level[[i]], lod.est, lod.level[[i]], col="blue")
      segments(lod.est, lod.level[[i]], lod.est, -1, col="blue")

      output.table1<-matrix(c(round(lod.est, 4), ci.lower, ci.upper), ncol=3)
      dimnames(output.table1)<-list(Selected.Model=model, c("LOD", "Lower", "Upper"))
      
      
      ##put all text output values into the same variable
      if (i==1){
        leg.text<-leg.text1
        output.table<-output.table1
      }else{
        leg.text<-rbind(leg.text, leg.text1)
        output.table<-rbind(output.table, output.table1)						
      }
    }else{
      leg.text<-c(paste("Cannot accurately determine the LOD based on the 95% CI."),
                  paste("Collect more data or select a different model."))
      output.table<-paste(model, "is not valid")
    }
  }
  legend("bottomright", leg.text, bg="gray", cex=0.8)
  cat("\n")
  print(output.table)
  return(range.error)		#return T/F to loop in selectModel		
}


if (F) {  
  ## Following is Unit Test
  LODFunc(model, r, var.cov, par.est, lod.fun, lod.level, strainX, max.range, graph.glm=T)
  ## Folowing is the expected result. 
  ## F (and a plot)
}




