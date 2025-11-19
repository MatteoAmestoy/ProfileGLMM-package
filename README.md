<h2>ProfileGLMM</h2>

This package implements a Bayesian profile regression using a generalized linear mixed model as output model. The package allows for binary (probit mixed model) and continuous (linear mixed model) outcomes and both continuous and categorical clustering variables. The package utilizes 'RcppArmadillo' and 'RcppDist' for high-performance statistical computing in C++.

<h4>Installation</h4>
  
Aspirational CRAN Installation
Once accepted to CRAN, you will be able to install the stable version directly using the following command:
```r
install.packages("ProfileGLMM")
```

<h4>Basic Usage</h4>


The primary function for sampling form the posterior of the parameters is profileGLMM_Gibbs()
```r
data("exposure_data")
exp_data = exposure_data$df

covList = {}
covList$FE = c('X')
covList$RE = c('t')
covList$REunit = c('indiv')
covList$Lat = c('X')
covList$Assign$Cont = c('Exp1','Exp2')
covList$Assign$Cat = NULL
covList$Y = c('Y')
dataProfile = profileGLMM_preprocess(regtype='linear',
                                     covList = covList,
                                     dataframe = exp_data,
                                     nC = 30,
                                     intercept = list(FE = T, RE = F, Lat = T))


MCMC_Obj = profileGLMM_Gibbs(model = dataProfile,
                             nIt = 5000,
                             nBurnIn = 2000)
```
