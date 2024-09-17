#!/bin/bash
OutDir=Pass1
mkdir -p "$OutDir"
echo Performing SLAC run
./Cassandra data/NNPDF/*      -wu12.5 -qu3.5 -r101 -l"$OutDir"/S_NNPDF_High      &
./Cassandra data/reanalyzed/* -wu12.5 -qu3.5 -r101 -l"$OutDir"/S_reanalyzed_High &

./Cassandra data/NNPDF/*      -wc -qc -r101 -l"$OutDir"/S_NNPDF_Low      &
./Cassandra data/reanalyzed/* -wc -qc -r101 -l"$OutDir"/S_reanalyzed_Low &
