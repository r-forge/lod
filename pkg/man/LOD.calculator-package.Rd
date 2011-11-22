\name{LOD.calculator-package}
\alias{LOD.calculator-package}
\alias{LOD.calculator}
\docType{package}
\title{
Limit of Detection (LOD) Calculator with Binary Data
}
\description{
  Limit of Detection (LOD) reflects an (medical diagnostic)
  instrument's level of sensitivity.  The smaller is the LOD, the
  higher is the sensitivity of the instrument to small quantities of an analyte
  (e.g. copies of DNA per reaction).  Seven binomial regression models
  are fitted to (binary) LOD measurement data and the best model is selected
  from AIC, BIC, or AICc criteria.  Multiple groups of LOD data can be
  incorporated into the model as Random Effect.  The 95% LOD, which is
  the level of analyte concentration that produces 95% positive
  results, and its Confidence Interval, is estimated from the
  selected model via the Delta Method.
}
\details{
\tabular{ll}{
Package: \tab LOD.calculator\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-11-22\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}

  The R Code is broken up into a number of functions.  The user will
 only interact with the function \code{selectModel}, which will then
 call the other functions needed to perform the calculations.  The
 calculations are performed primarily by fitting the models to a
 Generalized Linear Model (glm) from which an information criteria (AIC,
 BIC, or AICc) is found.  The model is selected based on the value of
 the criteria and the LOD and its 95\% CI (or other probability level) is
 found using the Delta Method.

}

\author{
Feiming Chen, Nick Huso
Maintainer: Feiming Chen <feiming\_chen@bd.com>
}

\references{
\url{http://en.wikipedia.org/wiki/Detection_limit}
}

\keyword{ package }
\examples{
strain1 <- data.frame(DNA = c(0, 2.5, 5, 10, 15, 20),
                      positive = c(0, 7, 9, 10, 10, 10))
selectModel(strain1, lod.level = c(.05, 0.5, .95), trials = 10, type = "AIC")
}
