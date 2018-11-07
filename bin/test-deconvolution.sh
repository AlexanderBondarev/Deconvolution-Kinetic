#!/bin/bash

Step=12000
DeltaTime=1500000

k1=0.65
A01=0.0035185
C01=0.0009
H1=-48000

k2=0.070
A02=0.0035185
C02=0.0027
H2=430000

./test-deconvolution -n $Step -dt $DeltaTime -k1 $k1 -A01 $A01 -C01 $C01 -H1 $H1 -k2 $k2 -A02 $A02 -C02 $C02 -H2 $H2 > test-deconvolution.csv
gnuplot test-deconvolution-svg.gnu
eog test-deconvolution.svg
