#!/bin/bash

../bin/test-deconvolution -n 12000 -dt 1500000 -k1 0.30 -A01 0.0035185 -C01 0.0020 -H1 -45000 -k2 0.055 -A02 0.0035185 -C02 0.0020 -H2 430000 > test-deconvolution-80.csv
gnuplot test-deconvolution-svg-80.gnu
eog test-deconvolution-80.svg
