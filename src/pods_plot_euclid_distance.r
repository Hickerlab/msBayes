install.packages( "ks" )
install.packages( "tree" )
library(ks)
library(tree)

SUMSTATS <- c(42:95,114:131)
PSI_INDEX <- 2
OMEGA_INDEX <- 5

SORTED_PRIORS <- "new_priors/sorted_buffer0.prior"
UNSORTED_PRIORS <- "new_priors/unsorted_buffer0.prior"
TEST_PRIORS <- "new_priors/buffer_0.01_10.prior"
OUTFILE <- "results_euclidistance.txt"
NUM_SIMS <- 10000
LENGTH_OF_PRIORS_FILES = 3000000

# Output matrices. The second dimension is 3
# because we're keeping track of diff psi
# diff omega, and euclidean distance of sumstats
OUT_SORTED <- matrix( NUM_SIMS, 3 )
OUT_UNSORTED <- matrix( NUM_SIMS, 3 )
# Output matrices for the simulated sumstats
SUMSTATS_SORTED <- matrix( NUM_SIMS, length(SUMSTATS))
SUMSTATS_UNSORTED <- matrix( NUM_SIMS, length(SUMSTATS))

# Plotting variables
YLIM_psi <- c(1, 18)
YLIM_omega <- c(0, 0.5)
XLIM <- c( 0, 0.4 )
TICKS <- seq( 0, 0.4, 0.05)

# I tried to do this the clever R way, but dealing with the large file
# was abysmally slow, so we do it the hackish awk way and it's quicker
# Still takes several minutes to gather 
AWK_CMD <- paste("awk 'BEGIN {srand()} !/^$/ { if (rand() <=", (NUM_SIMS+100)/LENGTH_OF_PRIORS_FILES, ") print $0}' <", sep=" " )

################################################################
# The euclidean distance function
################################################################
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
sum.sq <- function( x1 ) sqrt(sum(x1^2))

################################################################
# Euclidistance vs difference in parameters (omega and psi)
# for sorted vs unsorted simulations
################################################################
plot.new()
par(mfrow=c(2, 2))
for ( PRIOR_FILE in c( SORTED_PRIORS, UNSORTED_PRIORS)) {
  ## Generate pairs of matrices simulated data of equal length  
#  SIM1 <- matrix(nrow=NUM_SIMS, ncol=383)
#  SIM2 <- matrix(nrow=NUM_SIMS, ncol=383)
  
  cmd <- paste(AWK_CMD, PRIOR_FILE, sep=" " )

  SIM1 <- read.table( pipe(cmd), nrows=NUM_SIMS )
  SIM2 <- read.table( pipe(cmd), nrows=NUM_SIMS )

  ## Get absolute values of differences between psi and omega
  ## for each simulation pair
  diff_psi <- abs( SIM1[,PSI_INDEX] - SIM2[,PSI_INDEX] )
  diff_omega <- abs( SIM1[,OMEGA_INDEX] - SIM2[,OMEGA_INDEX] )
  diff_sumstats <- SIM1[,SUMSTATS] - SIM2[,SUMSTATS]
  sumstats_euclidist <- numeric( NUM_SIMS )
  for (i in 1: NUM_SIMS)
  {
    ## Get euclidean distance for sumstats
    ss_sum <- sum(diff_sumstats[i,]^2)
    sumstats_euclidist[i]<- sqrt(ss_sum)
    #write(c(eD),file=paste("eDfrogs"),ncol=1,append=T)
  }

  # Save the outputs for downstream processing  
  if ( PRIOR_FILE == SORTED_PRIORS ){
    OUT_SORTED <- cbind( diff_psi, diff_omega, sumstats_euclidist )
    SUMSTATS_SORTED <- diff_sumstats
  } else {
    OUT_UNSORTED <- cbind( diff_psi, diff_omega, sumstats_euclidist )
    SUMSTATS_UNSORTED <- diff_sumstats
  }
    
  plot( sumstats_euclidist, diff_psi, xlim=XLIM, ylim=YLIM_psi, main=paste(PRIOR_FILE, "psi vs euclidean distance of sumstats") )
  plot( sumstats_euclidist, diff_omega, xlim=XLIM, ylim=YLIM_omega, main=paste(PRIOR_FILE, "omega vs euclidean distance of sumstats"))

  fit <- lm(diff_omega~sumstats_euclidist)
  summary(fit)
  abline(fit, col="red")
  pts <- par('usr')
  text( pts[2]*.9, pts[4]*.1, paste("R^2 =", summary(fit)$r.squared ), col="red")
}

################################################################
# Euclidistance vs difference in parameters (omega and psi)
# for sorted vs unsorted simulations plotted by sumstat
# category, which are indexed like this in the original data, but
# we have to kludge it here cuz we chopped these vals out of the 
# retained stats.
#   pi: 42-59
#   Watterson: 60-77
#   pi.net: 78-95
#   Tajima's D: 114-131
################################################################
PI_INDEX <- c(1:18)
WATTERSON_INDEX <- c( tail(PI_INDEX, n=1)+1:+18 )
PI_NET_INDEX <- c(tail(WATTERSON_INDEX, n=1)+1:+18)
TAJIMA_INDEX <- c(tail(PI_NET_INDEX, n=1)+1:+18)
SUMSTAT_CATEGORIES <- list( PI_INDEX, WATTERSON_INDEX, PI_NET_INDEX, TAJIMA_INDEX )
SUMSTAT_NAMES <- list( "Pi", "Watterson", "Pi_Net", "V(pi-watterson)")

## Just recycle the stats from the last simulation
for ( PRIOR_FILE in c( SORTED_PRIORS, UNSORTED_PRIORS)){
  
  if ( PRIOR_FILE == SORTED_PRIORS ){
    STATS <- SUMSTATS_SORTED
    PARAMS <- OUT_SORTED
  } else {
    STATS <- SUMSTATS_UNSORTED
    PARAMS <- OUT_UNSORTED
  }
  
  #########################
  # Don't forget, the values in STATS
  # are _differences_ between the sumstats
  # which is why you have to apply the sum.sq
  # function across each row of the matrix
  #########################
  plot.new()
  par(mfrow=c(2, 2), oma=c(0,0,2,0))
  plot( apply(STATS[PI_INDEX], 1, sum.sq ), PARAMS[,1], ylab="Difference in psi", xlab="Difference in Pi", xlim=XLIM, ylim=YLIM_psi, main=paste(PRIOR_FILE, "psi estimated from only - Pi") )  
  plot( apply(STATS[WATTERSON_INDEX], 1, sum.sq ), PARAMS[,1], ylab="Difference in psi", xlab="Difference in Watterson",xlim=XLIM, ylim=YLIM_psi, main=paste(PRIOR_FILE, "psi estimated from only - Watterson") )  
  plot( apply(STATS[PI_NET_INDEX], 1, sum.sq ), PARAMS[,1], ylab="Difference in psi", xlab="Difference in Pi_net",xlim=XLIM, ylim=YLIM_psi, main=paste(PRIOR_FILE, "psi estimated from only - Pi_net") )  
  plot( apply(STATS[TAJIMA_INDEX], 1, sum.sq ), PARAMS[,1], ylab="Difference in psi", xlab="Difference in Tajima",xlim=XLIM, ylim=YLIM_psi, main=paste(PRIOR_FILE, "psi estimated from only - Tajima") )  
  title(main=paste(PRIOR_FILE, "sumstats by category"),font.main= 4, outer=T)
  
  sumstats_euclidist <- apply(STATS[PI_INDEX], 1, sum.sq )
  plot( sumstats_euclidist, PARAMS[,2], ylab="Difference in omega", xlab="Difference in Pi", xlim=XLIM, ylim=YLIM_omega, main=paste(PRIOR_FILE, "omega estimated from only - Pi") )  
  fit <- lm(PARAMS[,2]~sumstats_euclidist)
  summary(fit)
  abline(fit, col="red")
  pts <- par('usr')
  text( pts[2]*.9, pts[4]*.1, paste("R^2 =", summary(fit)$r.squared ), col="red")

  sumstats_euclidist <- apply(STATS[WATTERSON_INDEX], 1, sum.sq )
  plot( sumstats_euclidist, PARAMS[,2], ylab="Difference in omega", xlab="Difference in Watterson", xlim=XLIM, ylim=YLIM_omega, main=paste(PRIOR_FILE, "omega estimated from only - Watterson") )  
  fit <- lm(PARAMS[,2]~sumstats_euclidist)
  summary(fit)
  abline(fit, col="red")
  pts <- par('usr')
  text( pts[2]*.9, pts[4]*.1, paste("R^2 =", summary(fit)$r.squared ), col="red")

  sumstats_euclidist <- apply(STATS[PI_NET_INDEX], 1, sum.sq )
  plot( sumstats_euclidist, PARAMS[,2], ylab="Difference in omega", xlab="Difference in Pi_net", xlim=XLIM, ylim=YLIM_omega, main=paste(PRIOR_FILE, "omega estimated from only - Pi_net") )  
  fit <- lm(PARAMS[,2]~sumstats_euclidist)
  summary(fit)
  abline(fit, col="red")
  pts <- par('usr')
  text( pts[2]*.9, pts[4]*.1, paste("R^2 =", summary(fit)$r.squared ), col="red")

  sumstats_euclidist <- apply(STATS[TAJIMA_INDEX], 1, sum.sq )
  plot( sumstats_euclidist, PARAMS[,2], ylab="Difference in omega", xlab="Difference in Tajima", xlim=XLIM, ylim=YLIM_omega, main=paste(PRIOR_FILE, "omega estimated from only - Tajima") )  
  fit <- lm(PARAMS[,2]~sumstats_euclidist)
  summary(fit)
  abline(fit, col="red")
  pts <- par('usr')
  text( pts[2]*.9, pts[4]*.1, paste("R^2 =", summary(fit)$r.squared ), col="red")  
  title(main=paste(PRIOR_FILE, "sumstats by category"),font.main= 4, outer=T)
}
################################################################
# 2-D Kernel Density Estimates
################################################################

## 2-d kernel density estimate to evaluate the probability that
## the samples are drawn from the same distribution
## Good example from here: http://stats.stackexchange.com/questions/25946/goodness-of-fit-for-2d-histograms
kde.test( OUT_SORTED[,2:3], OUT_UNSORTED[,2:3] )
kde.test( OUT_SORTED[,1:3], OUT_UNSORTED[,1:3] )

## Sanity check, these should be the same
kde.test( OUT_SORTED[,1:2], OUT_UNSORTED[,1:2] )

################################################################
# T-tests for each Psi bin
################################################################

## T-test for each psi bin?
## Not very convincing results
# "1 0.207124813190253"
# "2 0.66048996864914"
# "3 0.73279841553834"
# "4 0.801724839909159"
# "5 0.872168591053327"
# "6 0.89107411266157"
# "7 0.909479220073261"
# "8 0.905622506064999"
# "9 0.937853672712299"
# "10 0.960143105347159"
# "11 0.959443360858703"
# "12 0.964108449377407"
# "13 0.965985570431417"
# "14 0.971040473150165"
# "15 0.988887975436563"
# "16 0.98824529041035"
# "17 0.991093276122544"

N_PSI <- max(OUT_SORTED[,1])
for ( i in 1:N_PSI) {
  print(paste( i, t.test(OUT_SORTED[OUT_SORTED[,1]==i], OUT_UNSORTED[OUT_UNSORTED[,1]==i], paired = F )$p.value))
}

################################################################
# Tajima's D vs V(pi-watterson)
# This is garbage.
################################################################
plot.new()
par(mfrow=c(2, 2), oma=c(0,0,2,0))

for ( PRIOR_FILE in c( SORTED_PRIORS, UNSORTED_PRIORS)){
  
  if ( PRIOR_FILE == SORTED_PRIORS ){
    STATS <- SUMSTATS_SORTED
    PARAMS <- OUT_SORTED
  } else {
    STATS <- SUMSTATS_UNSORTED
    PARAMS <- OUT_UNSORTED
  }
  
  sumstats_euclidist <- apply(STATS[TAJIMA_INDEX], 1, sum.sq )
  plot( sumstats_euclidist, PARAMS[,2], ylab="Difference in omega", xlab="Difference in Tajima", xlim=XLIM, ylim=YLIM_omega, main=paste(PRIOR_FILE, "omega estimated from only - Tajima") )  
  fit <- lm(PARAMS[,2]~sumstats_euclidist)
  summary(fit)
  abline(fit, col="red")
  pts <- par('usr')
  text( pts[2]*.9, pts[4]*.1, paste("R^2 =", summary(fit)$r.squared ), col="red")  
  title(main=paste(PRIOR_FILE, "sumstats by category"),font.main= 4, outer=T)

  sumstats_euclidist <- apply(STATS[PI_INDEX]-STATS[WATTERSON_INDEX], 1, var)
  plot( sumstats_euclidist, PARAMS[,2], ylab="Difference in omega", xlab="Difference in Tajima", ylim=YLIM_omega, main=paste(PRIOR_FILE, "omega estimated from only - Tajima") )  
  fit <- lm(PARAMS[,2]~sumstats_euclidist)
  summary(fit)
  abline(fit, col="red")
  pts <- par('usr')
  text( pts[2]*.9, pts[4]*.1, paste("R^2 =", summary(fit)$r.squared ), col="red")  
  title(main=paste(PRIOR_FILE, "sumstats by category"),font.main= 4, outer=T)
}