#!/bin/bash

for i in *.csv; do
 Name=`basename "$i" .csv`
 echo $Name
 ./sect <$i | grep -v "NaN" | sed -e "s/Normalized\ heat\ flow/Normalized_heat_flow/"> $Name.dat
 ./start-R.sh DeconvolutionKinetic.R $Name
done
