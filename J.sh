#!/bin/bash
OutDir=Pass1
mkdir -p "$OutDir"
echo Performing JLAB run
Cmd="./Cassandra data/j* -r101 -l'$OutDir'/J"
echo "$Cmd"
eval $Cmd &
