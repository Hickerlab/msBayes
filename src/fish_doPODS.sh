#!/bin/bash
#this script will make 100 estimates from 100 PODS simulated from batch file HU_IB_batchNEWpars_psi18_ne0.5_taumax0.75.txt, which is of the data configuration of Stone et al. data from 18 species pairs of parisitoid wasps codistributed across Iberian peninsula and Hungary
# see 2012 Stone, G.N., K. Lohse, J. A. Nicholls, P. Fuentes-Utrilla, F. Sinclair, K. SchÃ¶nrogge, G. CsÃ³ka, G. Melika, J-L. Nieves-Aldrey, J. Pujade-Villar, M. Tavakoli, R. R. Askew and M. J. Hickerson. Reconstructing Community Assembly in Time and Space Reveals Enemy Escape in a Western Palearctic Insect Community. Current Biology; 22: 532-537 
#in this example, the tau-max value is 0.75, but it varies from 1.0 to 0.03 in Table 2

#the prior file (prior_IB_HUPARSnewPsi00_3Mill) contains 3 million draws from the prior given identical data configuration and priors except for tau-max which is 1.0. 

usage(){
	echo "-p : Priors file 				(required)"
	echo "-c : msbayes conf file 			(required)"
	echo "-n : number of PODS			<default 1>"
	echo "-e : Do not generate new PODS \
		   instead use existing pods file	<filname>"
	echo "-s : Set sorting type         <default 7>"
	echo "-t : Do test run. Uses subset of the priors file"
}
# Must pass these in now, prevents errors
#MSBAYES_CONF="conf_buffer0.1_psi18_ne0.5_taumax0.75.txt"
#PRIORS_FILE="new_priors/buffer0.05.prior"

PODS_OBS="pods.obs"
PODS_POSTERIOR="posterior_table"
REJECT_THRESH=0.003
ITERATIONS=1
SORTING=7
EXISTING_PODS_FILE=""

while getopts p:n:s:c:te: flag; do
  case $flag in
	p)
		echo "Using Priors: $OPTARG";
		PRIORS_FILE=$OPTARG
		;;
        n)
      		echo "Doing Iterations: $OPTARG";
      		ITERATIONS=$OPTARG;
      		;;
        s)
      		echo "Doing sorting: $OPTARG";
      		SORTING=$OPTARG;
      		;;
        c)
          	echo "Using config: $OPTARG";
          	MSBAYES_CONF=$OPTARG
          	;;
	e)
		echo "Using existing pods file: $OPTARG";
		EXISTING_PODS_FILE=$OPTARG;
		# Do some sanity checks, e.g. the file should exist
		if [ ! -f "$EXISTING_PODS_FILE" ]; then
			echo "Existing pods file doesn't exist: $EXISTING_PODS_FILE"
			exit -1
		fi
	  	;;
        t)
          	echo "Doing TEST run. Set defaults so it'll run quick."
          	ITERATIONS=1
	  	# If you just want to run a test on 1/10 of the priors
	  	# uncomment these two lines, uses a smaller dataset and
	  	# a less stringent rejection filter
	  	TMP_PRIOR="/tmp/tmp_msbuff_prior"
	  	head -n 100000 $PRIORS_FILE > $TMP_PRIOR
	  	PRIORS_FILE=$TMP_PRIOR
	  	REJECT_THRESH=0.03
          	;;
    	?) 	usage; exit;
      	  	;;
	:) 	echo "Missing option argument for -$OPTARG" >&2; exit 1;;
  esac
done

# Sanity check on conf and priors files
if [ ! -f "$MSBAYES_CONF" ] || [ ! -f "$PRIORS_FILE" ]; then
	echo "Problem with input files: $MSBAYES_CONF - $PRIORS_FILE"
	exit -1
fi

# If using existing pods file make sure the file has enough pods 
# for the number of iterations requested
if [ "$EXISTING_PODS_FILE" != "" ]; then

	NPODS=$(awk 'END {print NR}' $EXISTING_PODS_FILE)
	echo "Number of existing pods: $NPODS"
	if [ $NPODS -lt $ITERATIONS ]; then
		echo "Existing pods file has $NPODS pods, you asked" \
			"for $ITERATIONS iterations"
		exit -1
	fi 
fi

echo "Doing PODS - $ITERATIONS"
echo "MSBAYES_CONF - $MSBAYES_CONF"
echo "PRIORS_FILE - $PRIORS_FILE"
sleep 5
for ((i=0; i<$ITERATIONS; i++))
do
	#clean up products of each loop
	rm -f $PODS_OBS
	rm -f $PODS_POSTERIOR

	#produces a single PODS via a draw from the prior (as opposed to the "leave one out" procedure)
	./msbayes.pl -r 1  -o $PODS_OBS -c $MSBAYES_CONF -s $SORTING

	#this above gets rid of the header file
	#stupid way, but it works
	tail -r $PODS_OBS | head -n 1 > tmp.txt
	mv tmp.txt $PODS_OBS

	# calls up the R script that takes out the first column so that the PODS 
	# file has dimensions identical to the big prior file
	# I don't think this is useful for us now, plus it trims out way more than the first column
	#R --quiet --vanilla < pods_trimcol.r

	# This is somewhat hackish but it keeps the code cleaner
	# If we are using existing pods just overwrite the $PODS_OBS file
	# with the nth line from the existing pods file
	#
	# The weird $(($i+1)) thing is to add 1 to $i since i is zero indexed
	# but file line numbers are 1 indexed
	if [ "$EXISTING_PODS_FILE" != "" ]; then
		head -n $(($i+1)) $EXISTING_PODS_FILE | tail -n 1 > $PODS_OBS
	fi

	#this does the simple ABC rejection filter on each PODS
	#make sure the correct columns are used (the numbers). These correspond the summary statistics that re used for the ABC procedure
	# The sumstats are
	# pi: 12-14
	# Watterson: 15-17
	# pi.net: 18-20
	# Tajima's D denominator: 24-26
	./msReject $PODS_OBS $PRIORS_FILE $REJECT_THRESH 12 13 14 15 16 17 18 19 20 24 25 26 > $PODS_POSTERIOR

	# Do rejection sampling
	# Set sampling threshold to 0.003 for real runs
	# being more lenient for testing runs so it'll run quicker on small datasets
	# This outputs the accepted parameter values and summary stats by default to the
	# file called 'posterior_table'
	#./acceptRej.pl -t $1 $PODS_OBS $PRIORS_FILE

	# this must be a script that does local linear regression etc
	# reads in by default from the 'posterior_table' file.
	R --quiet --vanilla < fish_pods_post.r

	echo "Done with iteration - " $i

done

#list of files needed for PODS horse race

#HU_IB_batchNEWpars_psi18_ne0.5_taumax0.75.txt
#batch file to configure the PODS simulation

#HUIBObsTAIL_columnReduction_ne0.5_taumax0.75.r
#calls up the R script that takes out the first column so that the PODS file has dimensions identical to the big prior file named prior_IB_HUPARSnewPsi00_3Mill

#prior_IB_HUPARSnewPsi00_3Mill
#large prior simulation file (aka reference table you will have to make via simulations first)

#HUIB_PODS_PARSnewPsi18_3M_HUIBPODStest_neuralnet_psi_omega_et_tol0.1_061512_ne0.5_taumax0.75.r
# this must be a script that does local linear regression etc

#and of course this script above
#Script_Table2_HUIB_PODS_PARSnewPsi18_3M_HUIBPODStest_neuralnet_psi_omega_et_tol0.1_110812_ne0.5_taumax0.75

#note, a lot of these files are called within downstream scripts from within the scripts

