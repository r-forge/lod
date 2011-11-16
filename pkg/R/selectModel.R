                                        # TODO: Add comment
                                        # 
                                        # Author: adehuson
##############################################################################

##Select the best model based on AIC, BIC, or AICc
selectModel  <-  function(positive, lod.level, trials, max.range=Inf,
                        type=c("AIC", "BIC", "AICc"),
                        A=T, B=T, C=T, L=T, CL=T, G=T, P=T, verbose=F) {
  ## March 11, 2011: Add: P (Probit Model) as the 7th model. 
  ## "positive": data
  ## "lod.level":level for LOD(95%LOD, 5%LOD, etc.)
  ## "trials":number of trials in data
  ## "max.range":the acceptable size of the 95% CI
  ## "type": desired criterion
  ## "verbose": print all criteria
  
  ##create data for use by functions
  proportion <- positive[1]/trials
  strainX <- cbind(positive, proportion, trials)
  colnames(strainX) <- list("positive", "DNA", "proportion", "trials")
  strainX  <-  strainX[order(strainX$DNA),] # sort by DNA increasing. 

  ##user input error messages
  if (max(strainX$proportion)>1)
    stop("Check number of trials entered.  Proportion greater than 1!")
  ## else if (max(strainX$proportion)<1)
  ##   stop("Check number of trials entered. Proportion less than 1!")
  
  if (A==F && B==F && C==F && L==F && CL==F && G==F && P==F)
    stop("No model selected!")
  
  if (max(lod.level)>=1||min(lod.level<0))
    stop("LOD must be between 0 and 1!")
  
  ##glm equations calling myGLM function
  glm.a  <-  myGLM(proportion ~ 0 + DNA, model.A(), strainX)
  glm.log  <-  myGLM(proportion ~ DNA, "logit", strainX)

  glm.b <-  myGLM(proportion ~ log(DNA + 1e-5), model.B(delta=1e-5), strainX)
  ## add delta to aviod log(0) error if first Model B returns error
  if (glm.b$error) {	
    glm.b  <-  glm.log  # copy model "glm.log" to model B to avoid complication
    B  <-  FALSE        # Drop Model B. 
  }

  glm.clog  <-  myGLM(proportion ~ DNA, "cloglog", strainX)
  if (glm.clog$error) {	
    glm.clog  <-  glm.log  # copy model "glm.log" over to avoid complication
  }

  glm.gom  <-  myGLM(proportion ~ DNA, Gompertz(), strainX)
  if (glm.clog$error) {	
    glm.gom  <-  glm.log  # copy model "glm.log" over to avoid complication
  }

  glm.probit  <-  myGLM(proportion ~ DNA, "probit", strainX)
  if (glm.probit$error) {	
    glm.probit  <-  glm.log  # copy model "glm.log" over to avoid complication
  }
  
  ##list glm equations for model selection later on
  all.models <- list(glm.a$r, glm.b$r, glm.log$r, glm.clog$r, glm.gom$r, glm.probit$r)
  
  ##function for calculating AIC 
  calcAic <- function(z){ 
    ifelse(is.null(z), Inf, AIC(z))
  } 
  ##function for calculating BIC 
  calcBic <- function(z){ 
    ifelse(is.null(z), Inf, AIC(z, k=log(nrow(strainX)))) 
  }  
  ##function for calculating AICc 
  calcAicc <- function(z){
    if (is.null(z))
      ans  <-  Inf
    else {
      k <- NROW(z$coefficients)
      ans  <-  AIC(z)+ 2*k*(k+1)/(nrow(strainX)-k-1)
    }
    ans
  }

  ##switch between user selected criterion
  type  <-  match.arg(type)
  switch(type,
         AIC = { 
           crit.all  <-  sapply(all.models,  calcAic)			#aic for glm models
           crit.c <- criterionC(strainX, type="AIC")				#aic for Model C calling criterionC
           if (is.na(crit.c))  C  <-  FALSE
           model.list <- order(c(crit.all, crit.c$crit))			#ordered criterion		
         },
         BIC = { 
           crit.all <- sapply(all.models, calcBic)
           crit.c <- criterionC(strainX, type="BIC")
           if (is.na(crit.c))  C  <-  FALSE
           model.list <- order(c(crit.all, crit.c$crit))		
         },
         AICc = { 
           crit.all <- sapply(all.models, calcAicc)
           crit.c <- criterionC(strainX, type="AICc")
           if (is.na(crit.c))  C  <-  FALSE
           model.list <- order(c(crit.all, crit.c$crit))
         }
         )

  ##print all criteria
  if (verbose==T){

    crit.list <- matrix(signif (c(crit.all[1], crit.all[2], crit.c$crit,
                               crit.all[3], crit.all[4], crit.all[5],
                                  crit.all[6]), 3), 7, 1,
                      dimnames=list(
                        Model= c("A", "B", "C", "Log", "CLog", "Gom", "Probit"),
                        Criterion= type))

    crit.list[, 1][!c(A, B, C, L, CL, G, P)]  <-  NA
    
    print(crit.list)
  }

  ##loop to select the best model and do plotting.
  for(i in seq(along=model.list)){

    ## Run best model through LODFunc
    if (model.list[i]==1 && A){
      LODFunc.a <- LODFunc("Model A", glm.a$r, glm.a$var.cov, glm.a$par.est, lod.A, lod.level, 
                         strainX, max.range)
      break			
    } else if (model.list[i]==2 && B){
      LODFunc.b <- LODFunc("Model B", glm.b$r, glm.b$var.cov, glm.b$par.est, lod.B, lod.level, 
                         strainX, max.range)
      break
    } else if (model.list[i]==3 && L){
      LODFunc.l <- LODFunc("Logistic", glm.log$r, glm.log$var.cov, glm.log$par.est, lod.D, lod.level, 
                         strainX, max.range)
      break
    } else if (model.list[i]==4 && CL){
      LODFunc.cl <- LODFunc("Comp. log-log", glm.clog$r, glm.clog$var.cov,
                            glm.clog$par.est, lod.E, lod.level,  
                            strainX, max.range)
      break
    } else if (model.list[i]==5 && G){
      LODFunc.g <- LODFunc("Gompertz", glm.gom$r, glm.gom$var.cov,
                           glm.gom$par.est, lod.F, lod.level, strainX, max.range)
      break
    } else if (model.list[i]==6 && P){
      LODFunc.p <- LODFunc("Probit", glm.probit$r, glm.probit$var.cov,
                         glm.probit$par.est, lod.P, lod.level, strainX, max.range)
      break			
    } else if (model.list[i]==7 && C){
      LODFunc.c <- LODFunc("Model C", crit.c$r, crit.c$var.cov,crit.c$par.est,
                           lod.C, lod.level, strainX, max.range, graph.glm=F)
      break			
    }
  }
}






