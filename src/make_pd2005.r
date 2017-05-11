makepdANY <- function(target,x,sumstat,tol,gwt,rejmethod=T)
{
# target is the set of target summary stats
# x is the parameter vector (long vector of numbers from the simulations) and is the dependent variable for the regression
# sumstat is an array of simulated summary stats (i.e. independent variables).
# NBB this function assumes 4 summary stats - I edited by hand for other numbers
# tol is the required proportion of points nearest the target values
# gwt is a vector with T/F weights, weighting out any 'bad' values (determined by the simulation program - i.e. nan's etc)
# if rejmethod=T it doesn't bother with the regression, and just does rejection.


# If rejmethod=F it returns a list with the following components:-

# $x regression adjusted values
# $vals - unadjusted values in rejection region (i.e. normal rejection)
# $wt - the regression weight (i.e. the Epanechnikov weight)
# $ss - the sumstats corresponding to these points
# $predmean - estimate of the posterior mean
# $fv - the fitted value from the regression

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

    if(rejmethod){
        l1 <- list(x=x[wt1],wt=0)
    }
    else{
        # this is proportional to k_delta(t) expression (5) of
        # Beaumont et al (2002)
        regwt <- 1-dst[wt1]^2/abstol^2

        # Note that scaled.sumstat is used instead of difference between
        # scaled.sumstat and target.s.  Since target.s (constant) is
        # subtracted from each simulated sumstat in Beaumont et al (2002),
        # we do not need to include it in the lsfit below.  The estimated
        # slopes are the same.  However intercept needs to be adjusted
        # below (see predmean).  predmean = hat(alpha) of the paper.
        fit1 <- lsfit(scaled.sumstat[wt1,],x[wt1],wt=regwt)
        predmean <- fit1$coeff %*% c(1,target.s)

        # x, correspond to phi*_i (between expressions (3) and (4) of the paper)
        l1 <- list(x=fit1$residuals+predmean,vals=x[wt1],wt=regwt,ss=sumstat[wt1,],predmean=predmean,fv = x[wt1]-fit1$residuals)

    }
    l1
}


normalise <- function(x,y){

if(var(y) == 0)return (x - mean(y))
(x-(mean(y)))/sqrt(var(y))
}
