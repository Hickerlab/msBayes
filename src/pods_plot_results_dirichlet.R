source("//Volumes/WorkDrive/msbayes-buffering/hickerlab-repository/msbayes-buffering/src/loc2plot.r")
library(locfit)

################################################################
# Housekeeping variables that describe your environment and data
################################################################

# The working directory. This should be whatever directory contains
# directories of your results files.
DIRICHLET_DIR <- ("//Volumes/WorkDrive/msbayes-buffering/hickerlab-repository/msbayes-buffering/data/dirichlet")
setwd(DIRICHLET_DIR)
BUFFER_LEVELS <- c(0, 0.01, 0.05, 0.1)
# 0.01 results not done yet
#BUFFER_LEVELS <- c(0, 0.05, 0.1)

# Set xlim and ylim so all plots are on the same scale
# These values are specific to the wasp dataset, the number
# number of species pairs defines the max number of codivergence events
xlm <- c(1, 18)
ylm <- c(1, 18)
ticks <- seq(1, 18)

# In the psi_omega results files, these are the columns:
# 1 - True Psi
# 2 - Estimated Psi
# 3 - RMSPE Psi
# 4 - True Omega
# 5 - Estimated Omega
# 6 - RMSPE Omega
# 7 - Estimated Omega (unadjusted)
# 8 - RMSPE Omega (unadjusted)


#####################################
# Plot observed vs estimated Psi
#####################################
plot.new()
par(mfrow=c(2, 2), oma=c(0,0,2,0))
for (buff in BUFFER_LEVELS) {

  PRIOR_FILE <- paste("./priors/buffer", buff, sep="")
  RESULTS_DIR <- paste("./results/buffer", buff, sep="")

  psi_omega <- read.table(paste(RESULTS_DIR, "results_psi_omega_et_tol.out", sep="/"))

  # In the sorted-diff-buffer case trip the dataset so no values are above
  # the max value the buffered prior is capable of inferring
  # effectively these values for soft masking psi estimates in the sorted-diff-buffer case
  #psi_omega <- subset(psi_omega, V1 <= max(psi_omega[,2]))

#  plot( psi_omega[,1], psi_omega[,2], xlim=xlm, ylim=ylm, at=ticks, xlab="True Psi", ylab="Estimated Psi", main="Dirichlet")
  plot( psi_omega[,1], psi_omega[,2], xlim=xlm, ylim=ylm, xlab="True Psi", ylab="Estimated Psi", main=buff)
  print( par('usr'))
  fit <- lm(psi_omega[,2]~ psi_omega[,1])
  summary(fit)
  abline(fit, col="red")
  abline( 0,1 )

  pts <- par('usr')
  text( pts[2]*.9, pts[4]*.1, paste("R^2 =", summary(fit)$r.squared ), col="red")
}
  title(main=paste(results_dir, "Observed vs Estimated Psi"),font.main= 4, outer=T)

  #####################################
  # Plot observed vs estimated Omega
  #####################################
  plot.new()
  par(mfrow=c(2, 2), oma=c(0,0,2,0))
  for (buff in BUFFER_LEVELS) {
    
    PRIOR_FILE <- paste("./priors/buffer", buff, sep="")
    RESULTS_DIR <- paste("./results/buffer", buff, sep="")

    print (paste( "Estimated vs Observed Omega: Buffer value =", buff))
    psi_omega <- read.table(paste(RESULTS_DIR, "results_psi_omega_et_tol.out", sep="/"))
  
    plot( psi_omega[,4], psi_omega[,5], xlab="True Omega", ylab="Estimated Omega", xlim=c(0,0.5), ylim=c(0,0.5), main=buff)
    fit <- lm(psi_omega[,5]~ psi_omega[,4])

# Uncomment these lines to plot the unadjusted values (e.g. in the case of sorted-diff-buffer)
# plot( psi_omega[,4], psi_omega[,7], xlab="True Omega", ylab="Estimated Omega", xlim=c(0,0.5), ylim=c(0,0.5), main=paste("Buffer = ", i ))
# fit <- lm(psi_omega[,7]~ psi_omega[,4])
  
    print( par('usr'))
    summary(fit)
    abline(fit, col="red")  
    abline( 0,1 )
  
    pts <- par('usr')
    text( pts[2]*.9, pts[4]*.1, paste("R^2 =", summary(fit)$r.squared ), col="red")
  }
  title(main=paste(results_dir, "Observed vs Estimated Omega"),font.main= 4, outer=T)

#####################################
# Plot histograms of RMSPE for Psi
#####################################
plot.new()
par(mfrow=c(2, 2), oma=c(0,0,2,0))
for ( buff in BUFFER_LEVELS ) {
  
  RESULTS_DIR <- paste("./results/buffer", buff, sep="")
  print (paste( "histograms of RMSPE for Psi: Buffer value =", buff))
  psi_omega <- read.table( paste(RESULTS_DIR, "results_psi_omega_et_tol.out", sep="/"))
  
  hist( psi_omega[,3], xlab="RMSPE Psi", ylab="Frequency (%)", ylim=c(0,35), xlim=c(0,18), breaks=10, main=paste("Buffer = ", buff ))
}
title(main=paste(results_dir, "RMSPE for Psi"),font.main= 4, outer=T)


#####################################
# Plot histograms of RMSPE for Omega
#####################################
plot.new()
par(mfrow=c(2, 2), oma=c(0,0,2,0))
for ( buff in BUFFER_LEVELS ) {
  RESULTS_DIR <- paste("./results/buffer", buff, sep="")
  print (paste( "histograms of RMSPE for Omega: Buffer value =", buff))
  psi_omega <- read.table( paste(RESULTS_DIR, "results_psi_omega_et_tol.out", sep="/"))
  
  hist( psi_omega[,6], xlab="RMSPE Omega", ylab="Frequency (%)", ylim=c(0,50), xlim=c(0,0.5), breaks=6, main=paste("Buffer = ", buff ))
}
title(main=paste(results_dir, "RMSPE for Omega"),font.main= 4, outer=T)


#####################################
# Plot histograms of prior PSI values
#####################################
plot.new()
par(mfrow=c(2, 2))
for ( buff in BUFFER_LEVELS ) {
  PRIORS_DIR <- paste("./priors/buffer", buff, sep="")
  ## Fuckin gotta love paste()
  PRIORS_FILE <- paste(paste(paste(PRIORS_DIR, "buffer", sep="/"), buff, sep=""), ".prior", sep="")
  RESULTS_DIR <- paste("./results/buffer", buff, sep="")
  print(paste("histograms of priors: Buffer value =", buff))

  priors <- read.table(PRIORS_FILE, nrows=10000, header=TRUE)
  
  hist( priors[,2], xlab="Prior Psi", ylab="# of draws", breaks=c(0:18), main=paste("Buffer = ", buff))
}
title(main=paste(results_dir, "Priors for PSI"),font.main= 4, outer=T)

#######################################
# Plot histograms of prior Omega values
#######################################
plot.new()
par(mfrow=c(2, 2))
for (buff in BUFFER_LEVELS) {
  print (paste( "histograms of priors: Buffer value =", buff))
  PRIORS_DIR <- paste("./priors/buffer", buff, sep="")
  ## Fuckin gotta love paste()
  PRIORS_FILE <- paste(paste(paste(PRIORS_DIR, "buffer", sep="/"), buff, sep=""), ".prior", sep="")
  priors <- read.table( PRIORS_FILE, nrows=10000, header=TRUE)
  
  hist( priors[,5], xlab="Prior Omega", ylab="# of draws", breaks=c(0:18), main=paste("Buffer = ", buff))
}

#####################################
# Plot kernel densities of estimated psi values
#####################################

install.packages( "vioplot" )
library(vioplot)

plot.new()
par(mfrow=c(2, 2))
for ( buff in BUFFER_LEVELS ) {
  
  RESULTS_DIR <- paste("./results/buffer", buff, sep="")
  print (paste( "Estimated vs Observed psi: buffer value =", buff))
  psi_omega <- read.table( paste(RESULTS_DIR, "results_psi_omega_et_tol.out", sep="/"))

  # Get the range of values for Psi and iterate over it
  # We want to build a list of vectors of values for each observed Psi bin
    mypts <- c()
    for( j in seq( 1, max(psi_omega[,1]))) {
      this_psi <- subset(psi_omega, psi_omega[,1]==as.character(j))
      if( length(this_psi[,2]) > 1  ) {
        mypts <- c( mypts, list(this_psi[,2]))
      } else {
      # pass. less than one point.
      }
    }

    # This is a hack to set the axes to all be equal. Plot
    # the points, then plot them again after the violinplot
    # for some reason you can't pass in xlim to vioplot()
    plot(psi_omega[,1], psi_omega[,2], col="gray", xlab="True Psi", ylab="Estimated Psi", xlim=xlm,ylim=ylm)#, at=ticks)

      # Draw the densities
    do.call(vioplot, c(c( mypts, wex=0.5), col="white", add=T))

    # This is the stupid way
    # vioplot(mypts[[1]], mypts[[2]], mypts[[3]], mypts[[4]], mypts[[5]], mypts[[6]], mypts[[7]], wex=.5, col="yellow")
    # Could use boxplot and a formula alternatively
    # boxplot(psi_omega[,2]~psi_omega[,1])
    
    # Add the actual points
    points( psi_omega[,1], psi_omega[,2], col="gray")

    # Draw the "expected line"
    abline(0,1)

    print( par('usr'))
    # Draw the regression line
    fit <- lm(psi_omega[,2]~ psi_omega[,1])
    summary(fit)
    abline(fit, col="red")
    
    pts <- par('usr')
    text( pts[2]*.8, pts[4]*.2, paste("R^2 =", summary(fit)$r.squared ), col="red")
}

#####################################
# Plot euclidean distance between delta
# parameter values and delta sumstats
# X-axis is difference between pairs of
# summary stats from pods vs priors
# Y-axis is difference between pairs of
# summary stats.
#
# This should only really be done for
# sorted and unsorted runs where buffer
# size is congruent between priors and
# pods, so these two:
# unsorted_samebufferPrior_PODS_results
# samebufferPrior_PODS_results
#####################################
#results_dir <- results_dirs[1]
PODS_FILE <- "all_pods.obs"

plot.new()
par(mfrow=c(2, 2), oma=c(0,0,2,0))
for ( buff in BUFFER_LEVELS ) {
  RESULTS_DIR <- paste("./results/buffer", buff, sep="")
  print (paste( "Plot euclidean distance betweeen sumstats and parameters: Buffer value =", buff))
  pods <- read.table( paste(RESULTS_DIR, PODS_FILE, sep="/"))
  
  hist( psi_omega[,3], xlab="RMSPE Psi", ylab="Frequency (%)", ylim=c(0,35), xlim=c(0,18), breaks=10, main=paste("Buffer = ", buff ))
}
title(main=paste(results_dir, "RMSPE for Psi"),font.main= 4, outer=T)
