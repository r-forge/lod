## Test Probit Model.

strain1 <- data.frame(positive = c(0, 7, 9, 10, 10, 10),
                      DNA = c(0, 2.5, 5, 10, 15, 20))

library(LOD.calculator)
## Run
selectModel(strain1, c(0.05, 0.5, 0.95), trials=10, type="AIC",
            verbose=T)

## Run "Probit" only. 
selectModel(strain1, c(0.05, 0.5, 0.95), trials=10, type="AIC",
            A=F, B=F, C=F, L=F, CL=F, G=F, verbose=T)


## 
