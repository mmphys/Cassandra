#!/usr/bin/env bash
OutDir=Pass1
mkdir -p "$OutDir"
echo Performing BCDMS run
Cmd="./Cassandra data/B* -wu12.5 -qu3.5 -sc -r101 -l'$OutDir'/B_High"
echo "$Cmd"
eval $Cmd &
Cmd="./Cassandra data/B* -wc -qc -sc -r101 -l'$OutDir'/B_Low"
echo "$Cmd"
eval $Cmd &
