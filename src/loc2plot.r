loc2plot <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,...)
{
        sc1 <- sqrt(var(x))
        sc2 <- sqrt(var(y))
        if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        maxk=100,mint=100,cut=0.8,maxit=100)
        else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        xlim=xlim,maxk=100,mint=100,cut=0.8,maxit=100)
        lev <- sort(fitted(fit))[floor(cprob*length(x))]
        plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
        xlim=gxlim,...)
	
}

loc2plotw <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,wt,...)
{
        sc1 <- sqrt(var(x))
        sc2 <- sqrt(var(y))
        if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        maxk=100,mint=100,cut=0.8,maxit=100,weight=wt)
        else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        xlim=xlim,maxk=100,mint=100,cut=0.8,maxit=100,weight=wt)
        lev <- sort(fitted(fit))[floor(cprob*length(x))]
        plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
        xlim=gxlim,...)
}

gethpdprob2 <- function(x,y,px,py,alpha=0.5,xlim,gxlim,...)
{
        sc1 <- sqrt(var(x))
        sc2 <- sqrt(var(y))
        if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        maxk=100,mint=100,cut=0.8,maxit=100)
        else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        xlim=xlim,maxk=100,mint=100,cut=0.8,maxit=100)
#       d1 <- (x-px)^2+(y-py)^2
#       best <- d1 == min(d1)
#       lev <- mean(fitted(fit)[best])
        lev <- predict.locfit(fit,list(px,py))
        slev <- sort(fitted(fit))
        indic <- slev <= lev
        sum(indic)/length(x)
}
        
        


loc2mode <- function(x,y,alpha=0.5,xlim,...)
{
        sc1 <- sqrt(var(x))
        sc2 <- sqrt(var(y))
        if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        maxk=100,mint=100,cut=0.8,maxit=100)
        else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        xlim=xlim,maxk=100,mint=100,cut=0.8,maxit=100)
        tt <- max(fitted(fit))
        wt <- fitted(fit) == tt
        c(x[wt][1],y[wt][1])
}


loc1stats <- function(x,prob=0.95,alpha=0.5,xlim,...)
{
        if(missing(xlim)) fit <- locfit(~x,alpha=alpha)
        else fit <- locfit(~x,alpha=alpha,xlim=xlim)
        fx <- fitted(fit)
        x.modef <- max(fx)
        x.mode <- x[fx == x.modef]
        if(length(x.mode)>1)x.mode <- x.mode[1]
        lev <- sort(fx)[floor(prob*length(x))]
       
        cat("INFO: loc1stats: difference in log(prob. density) between the ",
            "max (mode) and the ends of confidence interval (at p=", prob,
            ") is ", log(x.modef)-log(lev), ".\n", sep="")
		
        l1 <- list()
        l1[[1]] <- x.mode
        indx <- order(x)
        ii <- 2
        flip <- TRUE
        for(j in 2:length(x)){
                if(flip && fx[indx[j]] > lev){
                        l1[[ii]] <- x[indx[j-1]]
                        flip <- FALSE
                        ii <- ii+1
                }
                else if(!flip && fx[indx[j]] < lev){
                        l1[[ii]] <- x[indx[j]]
                        flip <- TRUE
                        ii <- ii+1
                }
        }
        as.numeric(l1)
}
                        
                

tloc2plot <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,...)
{
        sc1 <- sqrt(var(x))
        sc2 <- sqrt(var(y))
        if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha)
        else fit <- locfit(~x+y,alpha=alpha,
        xlim=xlim)
        lev <- sort(fitted(fit))[floor(cprob*length(x))]
        plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
        xlim=gxlim,...)
}
