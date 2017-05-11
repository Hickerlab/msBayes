## This script has been superceded by an ipython notebook
PODS_RESULTS <- ("//Volumes/WorkDrive/msbayes-buffering/hickerlab-repository/msbayes-buffering/data/empirical/results")
setwd(PODS_DIR)

orgs = c("frog", "fish")
#orgs = "frog"
for(organism in orgs){
  print(organism)

  DPP_SORTED <- paste(PODS_RESULTS, paste(organism, "-dirichlet-sort", sep=""), sep="/")
  DPP_UNSORTED <- paste(PODS_RESULTS, paste(organism, "-dirichlet-usort", sep=""), sep="/")
  UNIFORM_SORTED <- paste(PODS_RESULTS, paste(organism, "-uniform-sort", sep=""), sep="/")
  UNIFORM_UNSORTED <- paste(PODS_RESULTS, paste(organism, "-uniform-usort", sep=""), sep="/")

  OUTDIRS <- c(DPP_SORTED, DPP_UNSORTED, UNIFORM_SORTED, UNIFORM_UNSORTED)
  print(OUTDIRS)

  plot.new()
  par(mfrow=c(2, 2), oma=c(0,0,2,0))
  ## Plot Omega
  for (RESULTS_DIR in OUTDIRS) {
    RUN_NAME = tail((strsplit(RESULTS_DIR, "/"))[[1]], n=1)
    print (paste( "Estimated vs Observed Omega: Buffer value =", RUN_NAME))
    psi_omega <- read.table(paste(RESULTS_DIR, "results_psi_omega_et_tol.out", sep="/"))
  
    #plot( psi_omega[,4], psi_omega[,5], xlab="True Omega", ylab="Estimated Omega", xlim=c(0,1), ylim=c(0,1), main=RUN_NAME)
    #fit <- lm(psi_omega[,5]~ psi_omega[,4])
  
    # Uncomment these lines to plot the unadjusted values (e.g. in the case of sorted-diff-buffer)
    plot( psi_omega[,4], psi_omega[,7], xlab="True Omega", ylab="Estimated Omega", xlim=c(0,1), ylim=c(0,1), main=RUN_NAME)
    fit <- lm(psi_omega[,7]~ psi_omega[,4])
  
    print( par('usr'))
    summary(fit)
    abline(fit, col="red")  
    abline( 0,1 )
  
    pts <- par('usr')
    text( pts[2]*.9, pts[4]*.1, paste("R^2 =", summary(fit)$r.squared ), col="red")
  }
  title(main=paste("DPP vs Uniform / Sorted vs Unsorted", "Observed vs Estimated Omega"),font.main= 4, outer=T)

  ## Plot PSI
  for (RESULTS_DIR in OUTDIRS) {
    RUN_NAME = tail(strsplit(RESULTS_DIR, "/")[[1]], n=1)
    print (paste( "Estimated vs Observed Psi: Buffer value =", RUN_NAME))
    psi_omega <- read.table(paste(RESULTS_DIR, "results_psi_omega_et_tol.out", sep="/"))
  
    plot( psi_omega[,1], psi_omega[,2], xlab="True Psi", ylab="Estimated Psi", xlim=c(1,4), ylim=c(1,4), main=RUN_NAME)
    fit <- lm(psi_omega[,2]~ psi_omega[,1])
  
    print( par('usr'))
    summary(fit)
    abline(fit, col="red")  
    abline( 0,1 )
  
    pts <- par('usr')
    text( pts[2]*.9, pts[4]*.3, paste("R^2 =", summary(fit)$r.squared ), col="red")
  }
  title(main=paste("DPP vs Uniform / Sorted vs Unsorted", "Observed vs Estimated Psi"),font.main= 4, outer=T)
}
