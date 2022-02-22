# This used to be called: HUIB_PODS_PARSnewPsi18_3M_HUIBPODStest_neuralnet_psi_omega_et_tol0.1_061512_ne0.5_taumax0.75.r
#rm(list=ls(all=TRUE))
#print(getwd())
source("loc2plot.r")
source("make_pd2005.r")
source("calmod.r")


library(VGAM)
library(locfit)
library(abc)

Prior <- scan(file=paste("posterior_table"),skip=0)
POSTVEC<-matrix(Prior,ncol=383,byrow=T)
MZ_OBS<- scan(("pods.obs"),skip=0,nlines=1)
OBS<-matrix(MZ_OBS,ncol=383,byrow=T)

########function of displaying rmspe############s

### for abc function#####
rmspevalues<-function(true, values){
	SumSqure<-0
	nretained<-length(values)
	for (i in 1:nretained){
		SumSqure<-SumSqure+(values[i]-true)^2
	}
	return (sqrt(SumSqure/nretained))
}


##########start estimating Psi#################

Z1true1<- OBS[1,2]



####method 5########
# Original R script used 'rejection' but in the paper it says to use 'neuralnet'
#Z1noLL6 <- try(abc(OBS[1,c(42:95,114:131)], POSTVEC[,2], POSTVEC[,c(42:95,114:131)], tol=0.1,method="rejection"),T)
Z1noLL6 <- try(abc(OBS[1,c(42:95,114:131)], POSTVEC[,2], POSTVEC[,c(42:95,114:131)], tol=0.1,method="neuralnet"),T)
if(inherits(Z1noLL6, "try-error")){
	TDnoLLRmode6 <- NA
	rmspe6<-NA
} else { 
	TDnoLLRmode6 <- loc1stats(Z1noLL6$unadj.values, prob=0.95,alpha=10)[1]
#	rmspe6<-1111
	rmspe6<-rmspevalues(Z1true1, as.numeric(as.vector(Z1noLL6$unadj.values)))
}
TDnoLLRmode6
rmspe6

##########start estimating Omega#################

Z1true2<- OBS[1,5]

####method 4########
# Adjust so omega posterior is never negative, which happens in .1% of cases, rarely.
# OmegaNoNeg<-NA
# Also adjust so it is never bigger than the tau prior, in this case 0.75. In a better
# world you'd pass in this value, here it's just hard coded. Womp womp.
over<-0
under<-0
for(i in 1:dim(POSTVEC)[1]){
	if (POSTVEC[i,5]<0){ 
		POSTVEC[i,5]<-0
		under <- under+1
	}
	if ( POSTVEC[i,5] > 0.75 ){
		POSTVEC[i,5] <- 0.75
		over <- over+1
	}
}
over
under

# Original method was 'loclinear', but it fsck on small data?
#Z1noLL8 <- try(abc(OBS[1,c(42:95,114:131)], POSTVEC[,5], POSTVEC[,c(42:95,114:131)], tol=0.1,method="neuralnet"),T)
Z1noLL8 <- try(abc(OBS[1,c(42:95,114:131)], POSTVEC[,5], POSTVEC[,c(42:95,114:131)], tol=0.1,method="loclinear"),T)
if(inherits(Z1noLL8, "try-error")){
	TDnoLLRmode8 <- NA
	rmspe8<-NA

} else { 
	# get mode estimate of omega and rmspe for adjusted values
	TDnoLLRmode8 <- loc1stats(Z1noLL8$adj.values, prob=0.95,alpha=10)[1]
	rmspe8<-rmspevalues(Z1true2,as.numeric(as.vector(Z1noLL8$adj.values)))

	# We want to keep track of adjusted and unadjusted values for omega because
	# loclinear fsck on the unsorted sumstats for omega.

	# get mode estimate of omega and rmspe for unadjusted values
	TDnoLLRmode8_unadjusted <- loc1stats(Z1noLL8$unadj.values, prob=0.95,alpha=10)[1]
	rmspe8_unadjusted <-rmspevalues(Z1true2,as.numeric(as.vector(Z1noLL8$unadj.values)))
}
TDnoLLRmode8
rmspe8
TDnoLLRmode8_unadjusted
rmspe8_unadjusted


#Write out results files
write(c(Z1true1, TDnoLLRmode6, rmspe6, Z1true2, TDnoLLRmode8, rmspe8, TDnoLLRmode8_unadjusted, rmspe8_unadjusted), file=paste("results_psi_omega_et_tol.out"),ncol=8,append=T)
write(append(Z1true1,sort(as.numeric(as.vector(Z1noLL6$unadj.values)))), file=paste("results_psi.out"),ncol=901,append=T)
write(append(Z1true2,sort(as.numeric(as.vector(Z1noLL8$adj.values)))), file=paste("results_omega.out"),ncol=901,append=T)
write(append(Z1true2,sort(as.numeric(as.vector(Z1noLL8$unadj.values)))), file=paste("results_omega_unadjusted.out"),ncol=901,append=T)
write(OBS, file=paste("all_pods.obs"),ncol=383,append=T)
