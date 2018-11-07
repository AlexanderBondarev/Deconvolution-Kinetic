#!/bin/bash

../bin/test-deconvolution -n 12000 -dt 4206000 -k1 0.30 -A01 0.0035185 -C01 0.0018 -H1 -47000 -k2 0.067 -A02 0.0035185 -C02 0.0018 -H2 430000 > test-deconvolution-75.csv
gnuplot test-deconvolution-svg-75.gnu
eog test-deconvolution-75.svg
