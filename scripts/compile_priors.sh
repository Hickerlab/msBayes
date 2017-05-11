#!/bin/bash

if [ ! -d $2 ]
then
	echo "Second arg must be the directory with the priors"
	exit -1
fi

BUFF_FILE=`ls -1 $2/buffer_$1_* 2> /dev/null | head -n 1`
echo $BUFF_FILE
if [ -z $BUFF_FILE ]
then
	echo "No buffer by that name in $2"
	exit -1
fi

# Get the header
ls -1 $2/buffer_$1_* | head -n 1 | xargs head -n 1 >> $2/buffer$1.prior

for i in `ls $2/buffer_$1_*`
do
        echo $i
        tail -n +2 $i >> $2/buffer$1.prior
done
