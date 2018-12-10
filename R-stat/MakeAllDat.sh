#!/bin/bash

for i in *.csv; do
 Name=`basename "$i" .csv`
 echo $Name
 ./sect <$i | grep -v "NaN" > $Name.dat
 ./start-R.sh DeconvolutionKinetic.R $Name
done
