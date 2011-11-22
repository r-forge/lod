###  Select the best model based on AIC, BIC, or AICc
selectModel  <-  function(data, trials, group=NULL,
                          lod.level=0.95, 
                          max.range=Inf,
                          type=c("AIC", "BIC", "AICc"),
                          A=T, B=T, C=T, L=T, CL=T, G=T, P=T, verbose=F) {
  ## March 11, 2011: Add: P (Probit Model) as the 7th model. 
  ## "data": data, expects 1st column is covariate (analyte concentration),
  ##                       2nd column is number of positives. 
  ## "group": if non-NULL, a vector of group ID's of length nrow(data).
  ##           This is where we want to fit multiple sets of LOD curves together,
  ##           treating the groups/sets as random effects (as random intercept
  ##           on Systematic Components of GLMM model. 
  ## "lod.level":level for LOD(95%LOD, 5%LOD, etc.)
  ## "trials":number of trials in data, can be a single number or vector. 
  ## "max.range":the acceptable size of the 95% CI
  ## "type": desired criterion
  ## "verbose": print all criteria

  ## Create data for use by functions
  proportion <- data[2]/trials  # (num of positives / total)
  strainX <- cbind(data[1], proportion, trials) # DNA, prop, trials
  colnames(strainX) <- list("DNA", "proportion", "trials")
  has.group <- !is.null(group)
  if (has.group) {
    strainX <- cbind(strainX, group)
    colnames(strainX) <- list("DNA", "proportion", "trials", "groups")
  }

  strainX  <-  strainX[order(strainX$DNA),] # sort by DNA increasing. 

  ## user input error messages
  if (max(strainX$proportion)>1)
    stop("Check number of trials entered.  Proportion greater than 1!")
  ## else if (max(strainX$proportion)<1)
  ##   stop("Check number of trials entered. Proportion less than 1!")
  
  if (A==F && B==F && C==F && L==F && CL==F && G==F && P==F)
    stop("No model selected!")
  
  if (max(lod.level)>=1||min(lod.level<0))
    stop("LOD must be between 0 and 1!")
  
  ## glm equations calling myGLM function
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
           crit.all  <-  sapply(all.models,  calcAic) #AIC for glm models
           crit.c <- criterionC(strainX, type="AIC") #AIC for Model C criterionC
         },
         BIC = { 
           crit.all <- sapply(all.models, calcBic)
           crit.c <- criterionC(strainX, type="BIC")
         },
         AICc = { 
           crit.all <- sapply(all.models, calcAicc)
           crit.c <- criterionC(strainX, type="AICc")
         })
  if (is.na(crit.c$crit))  C  <-  FALSE
  model.list <- order(c(crit.all, crit.c$crit)) # ordered criterion


  model.names <-c("Model A", "Model B", "Logistic", "Comp-Log",
                  "Gompertz", "Probit", "Model C")
  model.range <- c(A, B, L, CL, G, P, C) # which models are considered (from User)
    
  ## Select the best model, estimate LOD (with CI), and do plotting.
  ## If "group" is specified, then re-estimate models with GLMM where "group"
  ##  serves as random intercept effect in the systematic component. Exception
  ##  for Model C: no random effects is assumed, use the global model (ignorming
  ##  grouping effects). 
  best.id <- model.list[model.range][1] # according to smallest AIC, BIC, or AICc
  best.model <- model.names[best.id]

  wLODFunc <- function(lod.fun) {
    ## Wraper to the "LODFunc"
    ## "lod.fun": is the FUNCTION for 95% LOD calculation.
    LODFunc(best.model, r, lod.fun, lod.level, strainX, max.range)
  }
  
  if (best.model== "Model A") {
    r <- ifelse(has.group, list(glm.a$rg), list(glm.a$r))[[1]]
    wLODFunc(lod.A)
  } else if (best.model== "Model B") {
    r <- ifelse(has.group, list(glm.b$rg), list(glm.b$r))[[1]]
    wLODFunc(lod.B)
  } else if (best.model== "Logistic") {
    r <- ifelse(has.group, list(glm.log$rg), list(glm.log$r))[[1]]
    wLODFunc(lod.D)
  } else if (best.model== "Comp-Log"){
    r <- ifelse(has.group, list(glm.clog$rg), list(glm.clog$r))[[1]]
    wLODFunc(lod.E)
  } else if (best.model== "Gompertz"){
    r <- ifelse(has.group, list(glm.gom$rg), list(glm.gom$r))[[1]]
    wLODFunc(lod.F)
  } else if (best.model== "Probit") {
    r <- ifelse(has.group, list(glm.probit$rg), list(glm.probit$r))[[1]]
    wLODFunc(lod.P)
  } else if (best.model== "Model C") {
    r <- crit.c
    wLODFunc(lod.C)
  }

  ## print all criteria
  if (verbose==T){
    crit.list <- matrix(signif (c(crit.all, crit.c$crit), 3), 7, 1,
                        dimnames=list(Model= model.names, Criterion= type))
    crit.list[, 1][!model.range]  <-  NA
    cat("\n")
    print(crit.list)
    cat("\n")
    print(r)
  }

  return(r)
}

if (F) {
  strain1 <- data.frame(DNA = c(0, 2.5, 5, 10, 15, 20),
                        positive = c(0, 7, 9, 10, 10, 10))
  selectModel(strain1, trials = 10, lod.level = c(.05, 0.5, .95),  verbose=T)

  (strain5 <- data.frame(DNA = rep(c(0, 2.5, 5, 10, 15, 20), 3),
                         positive = c(0, 7, 9, 9, 10, 10,
                           0, 2, 5, 7, 9, 9,
                           1, 3, 4, 6, 8, 10),
                         grp=rep(1:3, each=6)))
  selectModel(strain5[c(1,2)], lod.level = .95, trials = 10,
              group=strain5[3], type = "AIC", verbose=T)

  selectModel(strain5[c(1,2)], lod.level = .95, trials = 10,
              group=strain5[3], type = "AIC", verbose=T,
              A=T, B=T, C=T, L=T, CL=T, G=T, P=T)

}




