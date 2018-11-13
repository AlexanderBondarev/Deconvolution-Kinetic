
cc deconvolution-BFGS.c -o deconvolution-BFGS -I /usr/local/include -L /usr/local/lib -lgsl -lgslcblas -lm

./deconvolution-BFGS -k1 0.65 -k2 0.070 -C01 0.0009 -C02 0.0027 -H1 40 -H2 -430 -f experimental.dat | tee deconvolution-BFGS.log
