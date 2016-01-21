#Functions
library("foreach")
library("MASS")
library("quantreg")
library("KernSmooth")
library("doParallel")


kernelq = function(u){
  dnorm(u, mean = 0, sd = 1)
}

# M-type smoother
lnrob = function(x, y, h, maxiter, x0 = seq(0,1, length.out = 100)){
  
    xx = sort(x0)
    xx = (xx - min(xx))/(max(xx) - min(xx)) 
    fv = xx
    dv = xx
    x  = (x - min(x))/(max(x)-min(x))
  
    for (i in 1:length(xx)) {
        z     = x - xx[i]
        wx    = dnorm(z/h)
        r     =  rlm(y ~ z, weights = wx, method = "M", maxit = maxiter)
        u     =  r$wresid
        fv[i] =  r$coefficients[[1]]
        dv[i] =  r$coefficients[[2]]
    }
  
    "psi1" = r$psi

    return( list(xx = xx, fv = fv, dv = dv, "psi1" = psi1) )
}

# Quantile regression with specific tau
lprq2 = function(x, y, h, tau, x0) {
    xx = sort(x0) 
    xx = (xx - min(xx))/(max(xx) - min(xx)) 
    fv = xx
    dv = xx
    x  = (x - min(x)) / (max(x) - min(x)) 
    
    for(i in 1:length(xx)){
        z     = x - xx[i]
        wx    = dnorm(z/h)
        r     = rq(y ~ z, tau = tau, weights = wx, method = "br")
        fv[i] = r$coef[1.]
        dv[i] = r$coef[2.]
    }
    list(xx = xx, fv = fv, dv = dv)
}

# Quantile regression with random tau
lprq3 = function(x, y, h, x0){

    xx  = sort(x0) 
    xx  = (xx - min(xx)) / (max(xx) - min(xx))
    fv  = xx
    dv  = xx
    x   = (x - min(x)) / (max(x) - min(x)) 
    tau = runif(1) 
    
    for(i in 1:length(xx)) {
        z     = x - xx[i]
        wx    = dnorm(z/h)
        r     = rq(y ~ z, weights = wx, tau = runif(1), ci = FALSE) 
        fv[i] = r$coef[1.]
        dv[i] = r$coef[2.]
    }
    list(xx = xx, fv = fv, dv = dv)
}
