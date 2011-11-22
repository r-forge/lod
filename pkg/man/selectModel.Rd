\name{selectModel}
\alias{selectModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate LOD level.}
\description{
  It calculate LOD (Limit of Detection) level and its Confidence
  Interval for a given probabilty level from binary data.
}
\usage{
selectModel(data, trials, group=NULL, lod.level=0.95, 
            max.range=Inf, type=c("AIC", "BIC", "AICc"),
            A=T, B=T, C=T, L=T, CL=T, G=T, P=T, verbose=F)
%- maybe also 'usage' for other objects documented here.
}
\arguments{
  \item{data}{Data frame with 2 columns. 1st column is covariate
(e.g. analyte concentration, copies per reaction), 2nd column is number of positive results.}
  \item{trials}{Number of trials used at each concentration level. If it
  is a single number, it is replicated for each concentration level.}
\item{group}{If non-NULL, a vector of length nrow(data) for group ID's.
     This is useful when we want to fit multiple sets of LOD curves together,
     treating the groups/sets as random effects (specifically, as random intercept
     on Systematic Component of the GLMM model).} 
  \item{lod.level}{Consists of a list or vector of the levels the user wishes to see.  }
  \item{max.range}{Meta parameter which sets the maximum acceptable range for the 95\% CI. }
  \item{type}{Meta parameter that controls which criterion will be used to select the model based on a switch function.  }
  \item{A}{Optional user input which controls which models will be
    omitted. In this case, model A. }
  \item{B}{Ditto for model B}
  \item{C}{Ditto for model C}
  \item{L}{Ditto for Logistic Regression Model}
  \item{CL}{Ditto for Complimentary Log Model}
  \item{G}{Ditto for Gompertz Model}
  \item{P}{Ditto for Probit Model}
  \item{verbose}{Optional user input which allows the criteria for all the models to be displayed.}
}
\details{
  
  The R Code is broken up into a number of functions.  The user will
 only interact with this function, which will then call the other
 functions needed to perform the calculations.  The calculations are
 performed primarily by fitting the models to a Generalized Linear Model
 (glm) from which an information criteria (AIC, BIC, or AICc) is found.
 The model is selected based on the value of the criteria and the LOD
 and 95\% CI (or other probability level) is found using the Delta
 Method.
 
}
\value{
  A plot showing the fitted curve and estimated LOD. 
}
\author{Feiming Chen, Nick Huso}
\note{  
  
  There are several ways the user can force an error by inputting items
  incorrectly.

  \describe{
    
    \item{trials}{If the number of trials is less than the number of
      positive results, the application is stopped and we return the error
      message: "Check the number of trials entered".  If the number of trials
      is greater than the largest number of positive results, we return the
      same error message, but the application is not stopped as this may be
      due to poor data structure.}

    \item{models}{If the user attempts to omit all 6 models, the
      application is stopped and we return the error message: "No model
      selected!".}
      
      \item{lod.level}{If the user attempts to find an LOD greater than
	or equal to 1, the application is stopped and we return the error
	message: "LOD must be less than 1!".}
    }
      
 In some cases, the 95\% CI can become extremely wide or make little
sense. This makes estimation of the LOD inaccurate. In these cases, we
typically loop to the next "best" model with an acceptable 95\% CI. If
the user forces an incorrect output by omitting all valid models or if
there is no valid model we display a graph similar to the one seen
below. A list of invalid models will also be shown as text.

}
      
\examples{
strain1 <- data.frame(DNA = c(0, 2.5, 5, 10, 15, 20),
                      positive = c(0, 7, 9, 10, 10, 10))
selectModel(strain1, trials = 10, lod.level = c(.05, 0.5, .95), type = "AIC")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models} % __ONLY ONE__ keyword per line
