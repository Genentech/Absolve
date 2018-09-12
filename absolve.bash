#!/usr/bin/env bash


SCRIPTDIR=`dirname $0`

source $SCRIPTDIR/SETUP.bash "$SCRIPTDIR"

$SCRIPTDIR/absolve_rel "$@"

