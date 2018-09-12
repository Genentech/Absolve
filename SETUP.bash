#!/bin/bash
########################################################################


########################################################################

if [ ! -d "$1" ]
then
	>&2 echo "A valid location absolve application must be provided."
	>&2 echo "ERROR folder $1 not found"
fi

export ABSOLVE_HOME=$1




export ABSOLVE_HMMFILE_SCFV=$ABSOLVE_HOME/hmm/absolveSCFV.hmm
export ABSOLVE_HMMFILE_HEAVY=$ABSOLVE_HOME/hmm/absolveVH.hmm
export ABSOLVE_HMMFILE_LIGHT=$ABSOLVE_HOME/hmm/absolveVL.hmm

########################################################################


