#!/bin/bash

usage(){
	echo "It's your responsibility to be sure the buffer size specified in the"
	echo "passed in conf file is the same as the value passed in for -b here (which"
	echo "is only used for naming the output priors file but could still be confusing."
	echo ""
        echo "-b : buffer size            		(required)"
        echo "-c : msbayes conf file      		(required)"
	echo "-o : Output dir				<default is ./priors>"
        echo "-n : number of 250k priors files      	<default 1>"
	echo "-s : msbayes sort sumstats		<default 7>"
        echo "-t : Do test run. Run fast and generate 1 small priors file"
}

PRIORS_FILES=1
BUFFER_SIZE=-1
PRIORS_SIZE=250000
MSBAYES_SORT_STATS=7
OUTDIR=./priors

while getopts b:c:n:o:s:t flag; do
  case $flag in
        b)
        	echo "Doing Buffer Size: $OPTARG";
        	BUFFER_SIZE=$OPTARG
        	;;
        c)
          	echo "Using config: $OPTARG";
          	MSBAYES_CONF=$OPTARG
          	;;
        n)
      		echo "Generating this many 250k priors files: $OPTARG";
      		PRIORS_FILES=$OPTARG;
      		;;
        o)
      		echo "Output directory: $OPTARG";
      		OUTDIR=$OPTARG;
      		;;
	s)
		echo "Using msbayes sort version: $OPTARG";
		MSBAYES_SORT_STATS=$OPTARG
		;;
        t)
          	echo "Doing TEST run. Set defaults so it'll run quick."
		echo "Generate 1 prior file with 1000 draws, automatically"
		echo "cram it into /tmp"
          	PRIORS_SIZE=1000
          	OUTDIR=/tmp
		;;
        ?) usage; exit;
          ;;
       	:) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
  esac
done

if [ ! -f "$MSBAYES_CONF" ] || [ "$BUFFER_SIZE" == "-1"  ]; then
        echo "Problem with input either conf file is bad or no buffer size provided: $MSBAYES_CONF - $BUFFER_SIZE"
        exit -1
fi

for i in `seq 1 $PRIORS_FILES`;
do
	# Seed with random value because the perl seed will generate
	# identical priors if you run it too close together.
	./msbayes.pl -s $MSBAYES_SORT_STATS -r $PRIORS_SIZE -c $MSBAYES_CONF -S $RANDOM -o $OUTDIR/buffer_"$BUFFER_SIZE"_$i.prior &
done    
