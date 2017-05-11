## from popABC, Aug 05, 2009 version
## Naoki fixed that $x2 prints correct categories instead of mu1, mu2,...

calmod <- function(target,x,sumstat,tol,gwt,rejmethod=T)
{
# this function uses categorical regression to estimate the posterior probability of a particular model 
#      under the ABC framework P(Y=y | S = s)
#
# target is the set of target summary stats - i.e. s, what you measured from the data.
# x is the parameter vector, Y (long vector of numbers from the simulations) and is the dependent variable for the regression
#          This is a categorical variable, taking values 1, 2, .. n where there are n categories (models)
# sumstat is an array of simulated summary stats, S (i.e. independent variables).
# tol is the required proportion of points nearest the target values
# gwt is a vector with T/F weights, weighting out any 'bad' values (determined by the simulation program - i.e. nan's etc)
# if rejmethod=T it doesn't bother with the regression, and just does rejection.


# If rejmethod=F it returns a list with the following components:-

# $x1 expected value on a logit scale, with standard errors of the estimate
# $x2 expected value on a natural scale - i.e. p(Y=y | S = s) This is what we would normally report.
# $vals - the Ys in the rejection region. The proportion of the different Ys (different models) is a Pritchard etal-style, rejection-based
#                                         estimate of the posterior probability of the different models. You might get some improvement by 
#                                          weighting the frequencies with $wt.
# $wt - the Epanechnikov weight. 
# $ss - the sumstats corresponding to the points in $vals. 

  if(missing(gwt))gwt <- rep(T,length(sumstat[,1]))
  
  nss <- length(sumstat[1,])
  # scale everything 
  scaled.sumstat <- sumstat
    
  for(j in 1:nss){
    
    scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
  }
  target.s <- target
  
  for(j in 1:nss){
    target.s[j] <- normalise(target[j],sumstat[,j][gwt])
  }
    
  # calc euclidean distance
  sum1 <- 0
  for(j in 1:nss){
    sum1 <- sum1 + (scaled.sumstat[,j]-target.s[j])^2
  }
  dst <- sqrt(sum1)
  
  # includes the effect of gwt in the tolerance
  dst[!gwt] <- floor(max(dst[gwt])+10)
  
  # wt1 defines the region we're interested in 
  abstol <- quantile(dst,tol)
  # making sure all simulated results are included when tol = 1.
  if(tol == 1) {abstol <- abstol * 1.1}
  wt1 <- dst < abstol
  nit <- sum(wt1)
  
  if(rejmethod){
    l1 <- list(x=x[wt1],wt=0)
  }
  else{
    regwt <- 1-dst[wt1]^2/abstol^2
    
    catx <- as.numeric(names(table(x[wt1])))
    ncat <- length(catx)
    yy <- matrix(nrow=nit,ncol=ncat)
    for(j in 1:ncat){
      yy[,j] <- as.numeric(x[wt1] == catx[j])
    }
    # Used to print generic column names: mu1, mu2, ... for $x2
    # now it prints the "cagegory names" from x, Naoki
    colnames(yy)  <- names(table(x[wt1]))
    
    tr <- list()
    
    for(j in 1:nss){
      tr[[j]] <- scaled.sumstat[wt1,j]
    }
    
    xvar.names <- paste("v",as.character(c(1:nss)),sep="")
    
    names(tr) <- xvar.names
    
    fmla <- as.formula(paste("yy ~ ", paste(xvar.names, collapse= "+")))
    
#        fit1 <- vglm(fmla,data=tr,multinomial) this is the version described in the 
#               manuscript,which did not use the Epanechnikov weights. 
        

    fit1 <- vglm(fmla,data=tr,multinomial,weights=regwt)
        
    target.df <- list()
    for(j in 1:nss){
      target.df[[j]] <- target.s[j]
    }
    names(target.df) <- xvar.names
    
    prediction1 <- predict(fit1,target.df,se.fit=T)
    prediction2 <- predict(fit1,target.df,type="response")

    colnames(prediction2) <- colnames(yy) # added by Naoki
    
    l1 <- list(x1=prediction1,x2=prediction2,vals=x[wt1],wt=regwt,ss=sumstat[wt1,])
    
  }
  l1
}


normalise <- function(x,y){
  if(var(y) == 0)return (x - mean(y))
  (x-(mean(y)))/sqrt(var(y))
}
