#!/bin/bash

../bin/test-deconvolution -n 12000 -dt 1500000 -k1 0.65 -A01 0.0035185 -C01 0.0009 -H1 -48000 -k2 0.070 -A02 0.0035185 -C02 0.0027 -H2 430000 > test-deconvolution-85.csv
#./test-deconvolution -n 12000 -dt 1500000 -k1 0.42 k2 0.14 -C01 0.0002 -C02 0.0036 -H1 -45000 -H2 540000 -A01 0.0035185 -A02 0.0035185 > test-deconvolution-85.csv
gnuplot test-deconvolution-svg-85.gnu
eog test-deconvolution-85.svg
