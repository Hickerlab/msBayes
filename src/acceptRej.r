#library and R code  and package stuff
source("make_pd2005.r")
source("loc2plot.r")
source("calmod.r")
library(VGAM, warn.conflicts=F)
library(KernSmooth)
library(locfit)  # this need to be installed from CRAN 
# 1. To install locfit, download a package locfit_*.tar.gz from CRAN.
#    (replace * with the most recent version number of locfit, currently 1.5-3)
#    <http://cran.r-project.org/>
# 2. As a super user
#      R CMD INSTALL locfit_*.tar.gz
#    Or if you are using Mac OS-X, you can use sudo
#      sudo R CMD INSTALL locfit_*.tar.gz
# For more details, see <http://cran.r-project.org/doc/manuals/R-admin.html>
# There is a setion "Add-on packages" describing how to install packages.

# The first several columns of simDat contain the parameter values
# used for the simulations.  These values are sampled from the prior
# distributions.
# Example:
# When Psi constrainted to be 3
# PRI.numTauClass PRI.Tau.1       PRI.Tau.2       PRI.Tau.3       PRI.Psi.1       
# PRI.Psi.2       PRI.Psi.3       PRI.Psi PRI.var.t       PRI.E.t PRI.omega

# If there are 9 taxon pairs (or 3 taxon pairs with 3 genes per taxon
# pair), there are 9 columns of "pi.b" stats for each taxon pairs,
# then 9 columns of "pi.w", ... etc.
#
# Examples of summary stats (may not be up to date)
# c("pi.b", "pi.w", "pi", "wattTheta", "pi.net", "tajD", "tajD.denom", "pi.wPop2", "pi.wPop1", "wattTheta.Pop2", "wattTheta.Pop1", "tajD.denomPop2", "tajD.denomPop1", "ShannonsIndex.Between", "ShannonsIndex.Net", "ShannonsIndex.Pop1", "ShannonsIndex.Pop2"  )


##### The main function to do the standard analysis
# sim.infile and obs.infile are the file names of the input files
# usedStatNames: a list of stats used to summarize
#                These names should be from summary.stat.names
# return.res: If set to TRUE, the function returns a huge results in a list,
#             so you can analyze the results further.
# rejmethod: if True, it doesn't boether with the regression, and uses
#            simple rejection
# If tol=NA, tolearance is set to select 1000 closest matches.
stdAnalysis <- function(obs.infile, sim.infile, prior.infile, 
	                pdf.outfile="figs.pdf",
                        posterior.tbl.file="posterior_table",
                        tol=NA,
                        used.stats=c("pi","wattTheta","pi.net","tajD.denom"),
                        rejmethod=T, pre.rejected=F,
                        return.res=F
                        ) {
  simDat <- getData(sim.infile)
  if(pre.rejected) {
    priorHeader <- scan(prior.infile,what="character", nlines=1, quiet=T)
    if (length(priorHeader) != length(simDat$prior.names)) {
      cat("ERROR: prior file doesn't have the correct number of columns\n")
      return(NA)
    } else if (any(priorHeader != simDat$prior.names)) {
      cat("ERROR: prior file doesn't match simDat\n")
      return(NA)
    }
    prior.dat <- scan(prior.infile, skip=1, quiet=T)
    prior.dat <- data.frame(matrix(prior.dat, ncol=length(priorHeader), byrow=T))
    names(prior.dat) <- priorHeader
  }
  
  nPairs <- simDat[["numTaxonPairs"]]
  params.from.priorDistn <- simDat[["prior.names"]]
  summary.stat.names <- simDat[["summary.stats"]]
  simDat <- simDat[["dat"]]

  # if tol is NA, set default tol to get 1000 best matches.
  if (is.na(tol)) {
    tol <- 1000/nrow(simDat)
  }

  # construct the column names
  usedColNames <- as.vector(sapply(used.stats, paste,
                                   1:nPairs, sep="."))

  #load OBSERVED summary stat vector
  obsDat <-getData(obs.infile)
  if (obsDat[["numTaxonPairs"]] != nPairs) {
    cat("ERROR: The number of taxon pairs are not same between the\n      ",
        "observed data,", obsDat$numTaxonPairs, "pairs, and simulated data,",
        nPairs, "pairs.\n")
    return(NA)
  }
  obsDat <- obsDat[["dat"]]

  # acceptance/regression, .... ie  the meat
  # The 1st column is PRI.numTauClass, which should be removed from analysis
  result <- list(prior.names=
                 params.from.priorDistn[params.from.priorDistn != "PRI.numTauClass"])

  noPsiAnalysis <- F
  # The column tells the model of Psi
  if(simDat[1,"PRI.numTauClass"] != 0) { # constrained
    constrained <- T
    # get rid of PRI.Psi, which is always constant in constrained Psi
    result$prior.names <- result$prior.names[result$prior.names != "PRI.Psi"]

    numTaxonPairs <- sum(simDat[1,paste("PRI.Psi.",
                                        1:simDat[1,"PRI.numTauClass"],
                                        sep="")])
    # This means there are 3 taxon pairs, and constrained to have 3 tau's.
    # No need to analyze Psi.1 etc because they are always set to 1.
    # If we don't remove these, the analysis chokes
    if(numTaxonPairs == simDat[1,"PRI.numTauClass"]) {
      noPsiAnalysis <- T
      result$prior.names <-
        result$prior.names[- grep("^PRI[.]Psi[.]", result$prior.names)]
    } else {
      noPsiAnalysis <- F
    }
  } else {
    constrained <- F
  }
                           
  # divide prior.names into discrete and continuous priors
  if (noPsiAnalysis) {
    prior.names.cont <- result$prior.names
    prior.names.discrete <- character(0)
  } else {
    prior.names.cont <- result$prior.names[- grep("^PRI[.](Psi|model)", result$prior.names)]
    prior.names.discrete <- result$prior.names[  grep("^PRI[.](Psi|model)", result$prior.names)]
  }
  
  # set min, max boundary of prior parameters, and verbose print message
  # PRI.omega >= 0, PRI.E.t >= 0, PRI.Psi > 0
  min.val <- list(PRI.Psi = 1, PRI.var.t = 0, PRI.E.t = 0, PRI.omega = 0 )
  max.val <- list(PRI.Psi = nPairs)
  verbose.print <- list(PRI.omega = "(=Var(t)/E(t))",
                        PRI.Psi="(= number of possible divtimes)",
                        PRI.E.t="(= E(t))")
  # PRI.Tau* >= 0
  tauNames <- params.from.priorDistn[grep("^PRI[.]Tau", params.from.priorDistn)]
  if(length(tauNames) > 0) {
    temp.val <- sapply(rep(0, length(tauNames)), list)
    names(temp.val) <- tauNames
    min.val <- c(min.val, temp.val)
  }
  # 1 <= PRI.Psi.* <= nPairs - (numTauClasss - 1)
  if(constrained) {
    psiNames <- params.from.priorDistn[grep("^PRI[.]Psi[.]", params.from.priorDistn)]
    temp.val <- sapply(rep(1, length(psiNames)), list)
    names(temp.val) <- psiNames
    min.val <- c(min.val, temp.val)
    # I could use length of psiNames below, instead of simDat[1,1] to get numTauClass
    temp.val <- sapply(rep(nPairs - (simDat[1,"PRI.numTauClass"]-1), length(psiNames)), list)
    names(temp.val) <- psiNames
    max.val <- c(max.val, temp.val)
    
    temp.val <- sapply(rep("(= number of taxon pairs that divergence at corresponding tau)", length(psiNames)), list)
    names(temp.val) <- psiNames
    verbose.print <- c(verbose.print, temp.val)
  }
  
  # run makepdANY for each para for continuous
  for (i in 1:length(prior.names.cont)) {
    thisPriorName <- prior.names.cont[i]
    # might need to work on constrained vs unconstrained here
    temp <- list(makepdANY(as.vector(obsDat[1,usedColNames],mode="numeric"),
                           simDat[,thisPriorName], simDat[,usedColNames], tol,
                           rep(T,len=nrow(simDat)),rejmethod=rejmethod))
    names(temp) <- thisPriorName
    
    # absorbing boundary
    if ( !  is.null(min.val[[thisPriorName]])) {
      temp[[thisPriorName]]$x[which(temp[[thisPriorName]]$x<min.val[[thisPriorName]])] <- min.val[[thisPriorName]]
    }
    if ( !  is.null(max.val[[thisPriorName]])) {
      temp[[thisPriorName]]$x[which(temp[[thisPriorName]]$x>max.val[[thisPriorName]])] <- max.val[[thisPriorName]]
    }

    result <- c(result, temp)
  }

  # deal with the discrete priors
  calmod.fail <- c()
  if(length(prior.names.discrete) > 0) {
    for (i in 1:length(prior.names.discrete)) {
      thisPriorName <- prior.names.discrete[i]
      calmod.res <- try(calmod(as.vector(obsDat[1,usedColNames],mode="numeric"),
                               simDat[,thisPriorName], simDat[,usedColNames], tol,
                               rep(T,len=nrow(simDat)),rejmethod=rejmethod))
      if(class(calmod.res) == "try-error") {
        calmod.res <- calmod(as.vector(obsDat[1,usedColNames],mode="numeric"),
                             simDat[,thisPriorName], simDat[,usedColNames], tol,
                             rep(T,len=nrow(simDat)),rejmethod=T)
        calmod.fail <- c(calmod.fail, thisPriorName)
        this.failed <- T
      } else {
        this.failed <- F
      }
      temp <- list(calmod.res)
      names(temp) <- thisPriorName
      
      if (rejmethod || this.failed) {
        # with simple rejection, $x contains the accepted values
        # Need to copy to $vals to make the later handling easier.
        temp[[thisPriorName]]$vals <- temp[[thisPriorName]]$x
      }
      result <- c(result, temp)
    }
  }

  # create tables for posterior_table file
  big.table<-c()
  for (i in 1:length(result$prior.names)){
    thisPriorName <- result$prior.names[i]
    if(any(thisPriorName==prior.names.cont)) {
      # transformed values (or untransfomred with simple rejection)
      big.table<-cbind(big.table,result[[thisPriorName]]$x)
      if (! rejmethod)  # untransformed values
        big.table<-cbind(big.table,result[[thisPriorName]]$vals)
      # following is needed because colname vector is NULL at first
      # and assignment of names won't work.
      if (i==1) {colnames(big.table) <- colnames(big.table,do.NULL=FALSE)}
      if (! rejmethod)
        colnames(big.table)[ncol(big.table)-1] <- sub("PRI[.]", "Pos.LLR.", thisPriorName)        
      colnames(big.table)[ncol(big.table)] <- sub("PRI[.]", "Pos.wo.LLR.", thisPriorName)
    } else {  # discrete doesn't have transformed values
      big.table<-cbind(big.table,result[[thisPriorName]]$vals)
      if (i==1) {
        colnames(big.table) <- sub("PRI[.]", "Pos.wo.LCR.", thisPriorName)
      } else {
        colnames(big.table)[ncol(big.table)] <- sub("PRI[.]", "Pos.wo.LCR.", thisPriorName)
      }
    }
  }

  # combining the selected priors and summary stats and print
  write.table(cbind(big.table,result$PRI.Psi$ss),file=posterior.tbl.file,row.names = FALSE)
  
  real.mode.mean.median <- NULL
  modeToPrint <- NULL
  fileToPrint <- paste("acceptedPriorSummary_")
  for(i in 1:length(used.stats))
    {
      fileToPrint <- paste(fileToPrint, used.stats[i], sep = "")
    }

  tempMatrix <- scan(file = obs.infile, what = double(0), skip = 1, nmax = 5)
  TempMatrix <- matrix(tempMatrix, ncol = 1, byrow = T)
  truePRI <- c(TempMatrix[2,1], TempMatrix[3,1], TempMatrix[4,1], TempMatrix[5,1])
  counter <- 0

  cat("######### results #######\n")
  # make print loop for Mode and quantile
  for (i in 1:length(result$prior.names)) {
    counter <- counter + 1
    thisPriorName <- result$prior.names[i]
    name.rm.PRI <- sub("PRI[.]", "", thisPriorName)
    if(! is.null(verbose.print[[thisPriorName]])) {
      additional.print <- verbose.print[[thisPriorName]]
    } else {
      additional.print <- ""
    }

    cat ("##### Summary of posterior distribution #####\n")
    cat ("#####",name.rm.PRI, additional.print, "#####\n")

    # When constarined with large numTauClasses, PRI.Psi.1 become
    # mostly 1, and locfit has following problem, see help(locfit.raw)
    # Warning: procv: density estimate, empty integration region
    # Error: newsplit: out of vertex space Error: Descend tree proble
    # So, simply printing the posterior mean, and mode from
    # accepted values
    if (length(grep("^PRI\\.Psi(|\\.[0-9]+)$", thisPriorName)) == 1 ||
	thisPriorName %in% prior.names.discrete) {
      if (thisPriorName %in% calmod.fail) {
        this.calmod.failed <- T
      } else {
        this.calmod.failed <- F
      }
      
      # prior distribution
      if(pre.rejected) {
        this.prior.p <- table(prior.dat[,thisPriorName])
      } else {
        this.prior.p <- table(simDat[,thisPriorName])
      }
      this.prior.p <- this.prior.p / sum(this.prior.p)

      # transformed posterior prob.
      if ((! rejmethod) && (! this.calmod.failed))  {
        cat ("\n### Posterior probability table with local multinomial logit regression.\n")
        transformed.posterior.p.tbl <- (result[[thisPriorName]])$x2
        # removing "mu" from mu1, mu2, mu3 ...
        # colnames(transformed.posterior.p.tbl) <- sub("mu", "", colnames(transformed.posterior.p.tbl))

        p.tbl <- merge.2tbl.byName(transformed.posterior.p.tbl[1,], this.prior.p)
        rownames(p.tbl) <- c("posterior.p", "prior.p")
        print(p.tbl)
        
        cat ("\n## Mode (from local multinomial logit regression):\n")
        print(colnames(transformed.posterior.p.tbl)[which.max(transformed.posterior.p.tbl)])

        # posterior mean and median, a little weird with categorical var
        this.posterior.mean <- sum(as.numeric(colnames(transformed.posterior.p.tbl)) * transformed.posterior.p.tbl)
        this.posterior.median <- as.numeric(colnames(transformed.posterior.p.tbl)[ncol(transformed.posterior.p.tbl)])
        this.cumu <- 0
        for(j in 1:ncol(transformed.posterior.p.tbl)) {
          this.cumu <- this.cumu + transformed.posterior.p.tbl[1,j]
          if(this.cumu >= 0.5) {
            this.posterior.median <- as.numeric(colnames(transformed.posterior.p.tbl)[j])
            break
          }
        }
        mean.median.vect <- c(this.posterior.mean, this.posterior.median)
        names(mean.median.vect) <- c("mean", "median")
        cat ("## Mean/Median  (from local multinomial logit regression):\n")
        print(mean.median.vect)
        cat ("\n### The above results should be better than simple rejection method below.\n")
      }

      if(this.calmod.failed) {
        cat ("\n### WARNING: local mulitnomial logit regression FAILED\n")
        cat ("### WARNING: reporting results of SIMPLE rejection method\n")
      }
      cat ("\n### Posterior probability table with SIMPLE rejection method, NO REGRESSION.\n")

      raw.accepted.tbl <- table((result[[thisPriorName]])$vals)
      raw.posterior.p <- raw.accepted.tbl / sum(raw.accepted.tbl)
      
      p.tbl <- merge.2tbl.byName(raw.posterior.p, this.prior.p)
      rownames(p.tbl) <- c("posterior.p", "prior.p")

      print (p.tbl)
            
      cat ("\n## Mode (from simple rejection):\n")      
      print(names(raw.posterior.p)[which.max(raw.posterior.p)])

      cat ("## Mean/Median (from simple rejection)\n")
      this.mean.median.vect <- c(mean((result[[thisPriorName]])$vals),
                                 median((result[[thisPriorName]])$vals))
      names(this.mean.median.vect) <- c("mean", "median")
      print(this.mean.median.vect)

      post.distn.accRej <- c()
      if (rejmethod || this.calmod.failed) {
        post.distn.accRej <- raw.posterior.p
        mean.median.vect <- this.mean.median.vect
      } else {
        post.distn.accRej <- transformed.posterior.p.tbl
      }
      
      real.mode.mean.median <- append(real.mode.mean.median, c(truePRI[counter],1,mean.median.vect), after = length(real.mode.mean.median))
      #temp.pd.ar <- c("frequencies", post.distn.accRej)
      #names(temp.pd.ar)[1] <- name.rm.PRI      
    } else {  # Do regular summary, Not PRI.Psi.*
      cat ("### Mode\n")
      # With locfit 1.5_4, newsplit error of locfit() will stop the
      # analysis.  So I'm doing the error handling by myself with try().
      res.mode <- try(loc1stats((result[[thisPriorName]])$x,prob=0.95),silent=T)
      if(class(res.mode) == "try-error") {
        cat("NA\n")
        cat("** locfit failed to find mode, this occurs with an " ,
            "extreme\n** L-shaped posterior distribution, and the mode is ",
            "at the boundary (e.g. 0)\n", sep="")

        if(name.rm.PRI == "Psi")
          { modeToPrint <- 1.0 }
        else if((name.rm.PRI == "var.t")||(name.rm.PRI == "omega"))
          { modeToPrint <- 0.0 }
	else
	  { modeToPrint <- NA }
      } else {
        cat ("MODE:\n" )
        res.mode <- res.mode[1]
        print(res.mode)
        modeToPrint <- res.mode
      }
      cat ("### Mean/Median\n")
      mean.median.vect <- c(mean((result[[thisPriorName]])$x),
                                 median((result[[thisPriorName]])$x))
      names(mean.median.vect) <- c("mean", "median")
      print(mean.median.vect)
      
      cat ("### 95 % quantile\n")
      print(quantile((result[[thisPriorName]])$x,prob=c(0.025,0.975)))

      cat ("### 95% Highest Postererior Density (HPD) interval\n")
      print(emp.hpd((result[[thisPriorName]])$x,conf=0.95))
      
      real.mode.mean.median <- append(real.mode.mean.median, c(truePRI[counter], modeToPrint, mean.median.vect),after = length(real.mode.mean.median))
    }
    cat("\n")
  }

  # Wen was using the following, Mike said that we can disable this
  # write(real.mode.mean.median, file = fileToPrint, ncol = 20, append = T)
                         
  # Print out figures for continuous variables
  pdf(pdf.outfile, width=7.5, height=10, paper="letter")
  layout(mat=matrix(1:2,2,1))
  for (i in 1:length(prior.names.cont)) {
    thisPriorName <- prior.names.cont[i]
    name.rm.PRI <- sub("PRI[.]", "", thisPriorName)
    if(! is.null(verbose.print[[thisPriorName]])) {
      additional.print <- verbose.print[[thisPriorName]]
    } else {
      additional.print <- ""
    }
    
    this.title <- paste(name.rm.PRI, additional.print, sep=" ")
    if (pre.rejected) {
      make.hist(prior.dat[,thisPriorName],result[[thisPriorName]], title=this.title, breaks=20)
      
    } else {
      make.hist(simDat[,thisPriorName],result[[thisPriorName]],title=this.title,breaks=20)
#      plot.bf(simDat[,thisPriorName],result[[thisPriorName]]$x,main="Bayes Support for true Hyper-parameter value < threshold")
    }
  }
  
  # figures for discrete
  if (length(prior.names.discrete) > 0) {
    for (i in 1:length(prior.names.discrete)) {
      thisPriorName <- prior.names.discrete[i]
      name.rm.PRI <- sub("PRI[.]", "", thisPriorName)
      if(! is.null(verbose.print[[thisPriorName]])) {
        additional.print <- verbose.print[[thisPriorName]]
      } else {
        additional.print <- ""
      }
      
      this.title <- paste(name.rm.PRI, additional.print, sep=" ")
      this.legend <- c("posterior probability", "prior probability")
      
      if(thisPriorName %in% calmod.fail) {
        this.calmod.failed <- T
      } else {
        this.calmod.failed <- F
      }
      
      # prior distribution
      if(pre.rejected) {
        this.prior.p <- table(prior.dat[,thisPriorName])
      } else {
        this.prior.p <- table(simDat[,thisPriorName])
      }
      this.prior.p <- this.prior.p / sum(this.prior.p)
      
      # transformed posterior prob.
      if ((!rejmethod) && (!this.calmod.failed)) {
        this.pp.tbl <- result[[thisPriorName]]$x2
        # note $x2 is 1 x N matrix, and converting it to a vector (similar to table() output)
        barplot(merge.2tbl.byName(this.pp.tbl[1,], this.prior.p),beside=T,ylab="Posterior probability",
                legend=this.legend, main=paste(this.title,"With categorical regression", sep="\n"), space=c(0,0.05))
      }
      
      # posterior probability from simple rejection    
      this.pp.tbl <- 
        table(result[[thisPriorName]]$vals)/sum(table(result[[thisPriorName]]$vals))
      barplot(merge.2tbl.byName(this.pp.tbl, this.prior.p),beside=T,ylab="Posterior probability",
              legend=this.legend, main=paste(this.title,"With Simple Rejection",sep="\n"), space=c(0,0.05))
    }
  }
  

  if(pre.rejected) {
    plot((result[["PRI.omega"]])$x,(result[["PRI.E.t"]])$x,lty=2,lwd=0.5,ylim=c(0,max(prior.dat[["PRI.E.t"]])),xlim=c(0,max(prior.dat[["PRI.omega"]])))
  } else {
    plot((result[["PRI.omega"]])$x,(result[["PRI.E.t"]])$x,lty=2,lwd=0.5,ylim=c(0,max(simDat[["PRI.E.t"]])),xlim=c(0,max(simDat[["PRI.omega"]])))

  }
  rc <- try(plotKernDensity(result[["PRI.omega"]],result[["PRI.E.t"]],
                            xlab="Omega", ylab="E(t)", title="Omega and E(t)"))
  

  if(class(rc) == "try-error") {
    cat("WARN: plotKernDensity failed for some reason, so the kernel density ",
        "plot was not created\n", file=stderr())
  }
  # this plot doesn't seem to work.
  ## pdf("Skink0.5Milltol0.002Na0.5.pdf") 
  #  loc2plot(result.omega$x,result.Psi$x,cprob=0.6,alpha=0.4,xlab="omega", ylab="Psi")
  #  points(result.omega$x,result.Psi$x,pch=1,cex=0.7,col="darkgray")
  ## dev.off()

  cat("######### end of results #######\n")
  dev.off()

#  aaa <- list(nPairs=nPairs, simDat=simDat, obsDat=obsDat, result=result)
#  aaa
	
  if (return.res)
    return (list(nPairs=nPairs, simDat=simDat, obsDat=obsDat, result=result))
	
  
}


# This function takes an file name as an argument, read in the data
# from the file, and assign appropriate names.
# Returns a list(dat, numTaxonPairs)
# dat is the data.frame, numTaxonPairs is an integer indicating the number
# of taxon pairs.

getData <- function (infile) {
  first.line <- scan(infile, what="character", nlines=1, quiet=T) #header
  dat <- scan(infile, skip=1, quiet=T)
  dat <- data.frame(matrix(dat, ncol=length(first.line), byrow=T))
  names(dat) <- first.line  # assign the column names to the data.frame

  prior.names <- first.line[grep("^PRI[.]", first.line)]
  num.prior <- length(prior.names)
  # sum stats column-names (header) have the following form
  #   c("pi.b.1", "pi.b.2", "pi.b.3", "pi.w.1", "pi.w.2", "pi.w.3", ...)
  # Here, I'm getting rid of .digits part and taking unique names.
  sum.stat.names <- unique(sub("[.][0-9]+$", "", first.line[(num.prior+1):length(first.line)],fixed=F))  

  # number of taxon pairs can be calculated from
  nTaxPairs <-
    (ncol(dat) - num.prior) / (length(sum.stat.names))

  return (list(dat=dat, numTaxonPairs=nTaxPairs, prior.names=prior.names, summary.stats=sum.stat.names))
}


# -----------------------------------------------------------------------
# 2D Kernel density estimates: q1 X q2
# -----------------------------------------------------------------------
# horizontal view
plotKernDensity <- function (res1, res2, title="q1 and q2", xlab="q1", ylab="q2") {
  bwq1 <- try(dpik(res1$x),silent=T)
  bwq2 <- try(dpik(res2$x),silent=T)
  # I think bandwidth choice by dpik may fail if there aren't enough unique
  # values (e.g. mostly 0), I'm not sure following is ok or not, but
  # bw.nrd0 seems to be more robust
  if(class(bwq1) == "try-error") {
    cat("INFO: In plotKernDensity(), simpler bandwidth selection used for ",
            xlab, "\n", file=stderr())
    bwq1 <- bw.nrd0(res1$x)
  }
  if(class(bwq2) == "try-error") {
    cat("INFO: In plotKernDensity(), simpler bandwidth selection used for ",
            ylab, "\n", file=stderr())
    bwq2 <- bw.nrd0(res1$x)
  }  
  x <- cbind(res1$x, res2$x)
  est <- bkde2D(x, bandwidth=c(bwq1,bwq2))
  par(ask = FALSE)
  persp(est$x1, est$x2, est$fhat, theta = 145, phi = 25, col = "lightblue1",
        xlab=xlab, ylab =ylab,
        zlab = paste("Pr(", xlab, ", ", ylab, "| X)",sep=""),
        main = paste("Joint Density of", title), axes = TRUE, nticks = 5,
        ticktype = "detailed", ltheta = -135, lphi = 145, shade = TRUE)
}


make.hist <-function(vect, res.makepd, title="", xlim, ...) {
  #old.mfcol <- par()$mfcol
#  par(mfcol=c(3,1))
bw_vect<-max(vect)/100
bw_res.makepd<-max(res.makepd$x)/50
  hist.col = "white"
  if(missing(xlim)) {
    hist(vect,col=hist.col,border ="white",xlim=c(0,max(vect)),ylim=c(0,max(density(vect,bw=bw_vect)$y,density(res.makepd$x,bw=bw_res.makepd)$y)),prob=TRUE,main=paste(title),xlab=title, ...)
  } else {
    hist(vect,col=hist.col,freq=F,xlim=xlim,prob=TRUE,main=paste(title),
       xlab=title, ...)
  }
  lines(density(vect,bw=bw_vect),lty=2,col="red",pch=3)
#  if(missing(xlim)) {
#    hist(res.makepd$x,col="blue",prob=TRUE,freq=F,xlim=c(0,max(vect)),
#         main=paste(title, ": Posterior Dist'n "), xlab=title, ...)
#  } else {
#    hist(res.makepd$x,col=hist.col,xlim=xlim,freq=F,prob=TRUE,
#         main=paste(title, ": Posterior Dist'n "), xlab=title, ...)
#  
  lines(density(res.makepd$x,bw=bw_res.makepd),lty=1,col="blue",pch=3)
  #par(mfcol=old.mfcol)
legend("topright",c("Prior","Posterior"),text.col=c("red","blue"),lty=c(2,1),col=c("red","blue"))

}

plot.bf <- function(prior,posterior,...) {
  bf.out <- make.bf.vect(prior,posterior)
  if (is.null(bf.out)) {
    return(NULL);
  }
  # finding some reasonable y.lim for the plot, all values here are arbitrary
  y.max <- 40
  if (max(bf.out$bf) * 1.1 < 40) {
    y.max <- max(bf.out$bf)* 1.1
  } else if (mean(bf.out$bf > 40) > 0.6) {
    y.max <- sort(bf.out$bf)[ceiling(length(bf.out$bf) * 0.4)] * 1.1
  }
  plot(bf.out$crit.vals, bf.out$bf, type="l",
       ylim=c(0,y.max), 
       xlab="Hyper Parameter Thresholds",ylab="Bayes Factor",...)
  line.type = 3
  abline (1, 0, lty=2)
  abline(3, 0, lty=line.type)
  abline(10, 0, lty=line.type)
  abline(1/3, 0, lty=line.type)
  abline(1/10, 0, lty=line.type)
}

# Compares two Model of x < crit.val and x >= crit.val, and
# calculate the bayes factor.
# It returns a data.frame with two columns
make.bf.vect <- function(prior, posterior, crit.vals = NULL, num.points=100) {

  prior.len <- length(prior)
  posterior.len <- length(posterior)
  if (prior.len * posterior.len == 0) {  # some error
    return(NULL);
  }
  
  if (is.null(crit.vals)) {

    post.range <- range(posterior)

    # A single value of posterior is observed (e.g. at the boundary)
    if(post.range[1] == post.range[2]) {
      if(post.range[1] == 0) {
        post.range <- c(-0.01, 0.01)
      } else  {
        post.range <- c(0.99, 1.01) * post.range[1]
      }
    }

    # following is probably not required
    pri.min <- min(prior)
    if (post.range[1] <= pri.min) {
      if(pri.min >= post.range[2]) {  # impossible, but better check
        cat("WARN: In make.bf.vect(), weird prior/posterior range encountered\n",
            file=stderr())
        return(NULL)
      }
      post.range[1] <- pri.min
    }
    crit.vals <- seq(post.range[1], post.range[2], length.out=num.points+2)
    # get rid of the ends to avoid bf = Inf
    crit.vals <- crit.vals[2:(length(crit.vals) - 1)]
  }
  
  if (is.null(crit.vals)) {
    temp.pos <- posterior[posterior != max(posterior) &
                          posterior != min(posterior)]
     post.range <- range(temp.pos)
     crit.vals <- seq(post.range[1], post.range[2], length.out=num.points)
  }
  
  bf.vect <- numeric (0)
  for(i in 1:length(crit.vals)) {
    prior.below <- length(which(prior < crit.vals[i]))
    posterior.below <- length(which(posterior < crit.vals[i]))
    if (prior.below == 0 || posterior.below == posterior.len) {
      this.bf <- Inf
    } else {
      this.bf <- posterior.below * (prior.len-prior.below)/
        ((posterior.len - posterior.below) * prior.below)
    }
    bf.vect <- c(bf.vect, this.bf)
  }

  return (data.frame(crit.vals=crit.vals, bf=bf.vect))
}

# Marging two tables created by table()
# this is not generic, only for 1 dimensional table
merge.2tbl.byName <- function(arr1, arr2) {
  #this.name <- dimnames(arr1)[[length(dimnames(arr1))]]
  this.name <- names(arr1)
  d1 <- data.frame(cbind(this.name, arr1))
  names(d1)[1] <- 'name'
  
  #this.name <- dimnames(arr2)[[length(dimnames(arr2))]]
  this.name <- names(arr2)
  d2 <- data.frame(cbind(this.name, arr2))
  names(d2)[1] <- 'name'
  
  m1 <- merge(d1,d2,by='name',all=T)
  m1 <- m1[order(as.numeric(levels(m1$name)[m1$name])),]  #sorting rows
  m1 <- as.matrix(m1)
  m1[is.na(m1)] <- 0  # replacing NA with 0
  result <- as.table(t(m1[,2:3]))
  # previous t() doesn't work if it has only 1 row, so do it here
  if(nrow(m1) == 1) {
    result <- t(result)
  }
  colnames(result) <- m1[,1]
  # convert character table to numeric
  return(type.convert(result))
}

## This is copied from emp.hpd() of TeachingDemos by
## Greg Snow <greg.snow@intermountainmail.org>

# This function computes the highest posterior density intervals
# (sometimes called minimum length confidence intervals) for a
# Bayesian posterior distribution. The emp.hpd function is used when
# you have realizations of the posterior (when you have results from
# an MCMC run).
emp.hpd <- function (x, conf = 0.95) {
  conf <- min(conf, 1 - conf)
  n <- length(x)
  nn <- round(n * conf)
  x <- sort(x)
  xx <- x[(n - nn + 1):n] - x[1:nn]
  m <- min(xx)
  nnn <- which(xx == m)[1]
  return(c(x[nnn], x[n - nn + nnn]))
}
