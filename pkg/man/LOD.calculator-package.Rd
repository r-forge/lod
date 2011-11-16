\name{LOD.calculator-package}
\alias{LOD.calculator-package}
\alias{LOD.calculator}
\docType{package}
\title{
Limit of Detection (LOD) Calculator with Binary Data
}
\description{
  Limit of Detection(LOD) is a means of measuring an
instrument's level of sensitivity. The smaller the LOD, the higher the
sensitivity of the instrument to small quantities of an analyte (in
this case copies of bacteria DNA). An accurate estimate of LOD will
give a level of cut-off that separates positive and negative results
(in our case the detection or lack of detection for bacteria DNA). Our
goal is to consolidate various LOD models into a more powerful and
user-friendly web application featuring easier access, less manual
intervention, and more information in the form of warnings regarding
inadequacies of input data.
}
\details{
\tabular{ll}{
Package: \tab LOD.calculator\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-09-10\cr
License: \tab BD Proprietary\cr
LazyLoad: \tab yes\cr
}

  The R Code is broken up into a number of functions.  The user will
 only interact with the function \code{selectModel}, which will then
 call the other functions needed to perform the calculations.  The
 calculations are performed primarily by fitting the models to a
 Generalized Linear Model (glm) from which an information criteria (AIC,
 BIC, or AICc) is found.  The model is selected based on the value of
 the criteria and the LOD and 95\% CI (or other probability level) is
 found using the Delta Method.

}

\author{
Nick Huso

Maintainer: Feiming Chen <feiming\_chen@bd.com>
}

\references{
\url{http://bancs.bdx.com/wiki/LOD}
}

\keyword{ package }
\examples{
strain1 <- data.frame(positive = c(0, 7, 9, 10, 10, 10),
                      DNA = c(0, 2.5, 5, 10, 15, 20))
                      
selectModel(strain1, lod.level = c(.05, 0.5, .95), trials = 10, type = "AIC")
}
