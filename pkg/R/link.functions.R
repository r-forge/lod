# TODO: Add comment
# 
# Author: adehuson
###############################################################################

p <- 0.95                               # default global variable. 

###### MODEL A #############
## Custom Link Function
model.A <- function(delta = 1e-5) { 
	linkfun <- function(mu) -log(1 + delta - mu) 
	linkinv <- function(eta) 1 + delta - exp(-eta) 
	mu.eta <- function(eta) exp(-eta) 
	valideta <- function(eta) TRUE 
	link <- "Model-A-Link-Function" 
	structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta =
							mu.eta, valideta = valideta, name = link), class = "link-glm")
}

lod.A <- function(x) {
	## The 95% LOD function required by the Delta Method.
	## "x"(beta) is the coefficients estimated from the GLM model.
	## NOTE: D(parse(text=deparse(lod.A)[3]), "x") will return the deriavive
	##       of this function with respect to "x". 
	- log(1 - p)/x
}



###### MODEL B #############
model.B <- function(delta = 0) { 
	linkfun <- function(mu) log(delta - log(1 + delta - mu))         
	linkinv <- function(eta) 1 + delta - exp(delta - exp(eta))
	mu.eta <- function(eta) exp(eta - exp(eta) + delta)
	valideta <- function(eta) TRUE         
	link <- "Model B Link Function"         
	structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = 
							mu.eta, valideta = valideta, name = link), class = "link-glm") 
} 

lod.B <- function(x, y) {
	## The 95% LOD function required by the Delta Method.
	## "x"(beta), "y"(theta) is the coefficients estimated from the GLM model.
	exp((log(-log(1 - p)) - x) / y)
}



##### Model C #####	
lod.C<-function(x, y){
	## The 95% LOD function required by the Delta Method.
	## "x"(beta), "y"(theta) are the coefficients.
	-log(1-p^y)/x
}



##### Model D (logit) #####

lod.D <- function(x, y) {
	## The 95% LOD function required by the Delta Method.
	## "x"(beta), "y"(theta) is the coefficients estimated from the GLM model.
	(log(p/(1-p)) - x) / y
}



##### Model E (Complimentary log-log) #####
## NOTE: It is different from Model B since the predictor is DNA, not log(DNA)!

lod.E <- function(x, y) {
	## The 95% LOD function required by the Delta Method.
	## "x"(beta), "y"(theta) is the coefficients estimated from the GLM model.
	(log(-log(1-p)) - x) / y
}



##### Model F (Gompertz Growth Curve as link function) ####

lod.F <- function(x, y) {
	## The 95% LOD function required by the Delta Method.
	## "x"(beta), "y"(theta) is the coefficients estimated from the GLM model.
	(log(-log(p)) - x) / y
}

Gompertz <- function (delta=0) {                                                 
	linkfun <- function(mu) log(delta - log(mu + delta))       
	linkinv <- function(eta) exp(delta - exp(eta)) - delta     
	mu.eta <- function(eta) -exp(eta + delta - exp(eta)) 
	valideta <- function(eta) TRUE
	link <- "Gompertz Link Function"                                         
	structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
					valideta = valideta, name = link), class = "link-glm") 
} 

## ADD: March  2, 2011 
##### Model P (Probit Model) ####

lod.P <- function(x, y) {
  ## The 95% LOD function required by the Delta Method.
  ## "x"(beta), "y"(theta) is the coefficients estimated from the GLM model.
  ## PHI^{-1}(p) = beta.x + theta.y * LOD
  ## ==> LOD = (PHI^{-1}(p) - beta.x) / theta.y
  ## Let p = qnorm(q), where q=0.95 (for 95% LOD).
  (p - x) / y
}



