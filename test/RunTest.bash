#!/usr/bin/env bash

../absolve.bash --R1 ./simulatedseqs.fa --VH ../dbs/human/VHv.fa --outdir ./testrun

if !  diff -r -q --exclude=".svn" testrun testrun.ref ; then 
	echo TestFailed
	exit -1
else
	echo Success
	exit 0
fi



