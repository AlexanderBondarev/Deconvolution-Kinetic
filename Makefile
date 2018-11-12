
CC=gcc
CFLAGS=-c
BFGS_FLAGS=-I /usr/local/include -L /usr/local/lib -lgsl -lgslcblas -lm

all: deconvolution deconvolution-BFGS sect test-deconvolution

deconvolution: src/deconvolution.c
	$(CC) $(CFLAGS) src/deconvolution.c -o bin/deconvolution
	chmod a+x bin/deconvolution

deconvolution-BFGS: src/deconvolution-BFGS.c
	cc src/deconvolution-BFGS.c -o bin/deconvolution-BFGS $(BFGS_FLAGS)
	chmod a+x bin/deconvolution-BFGS

sect: src/sect.c
	$(CC) $(CFLAGS) src/sect.c -o bin/sect
	chmod a+x bin/sect

test-deconvolution: src/test-deconvolution.c
	$(CC) $(CFLAGS) src/test-deconvolution.c -o bin/test-deconvolution
	chmod a+x bin/test-deconvolution

clean:
	rm bin/deconvolution bin/deconvolution-BFGS bin/sect bin/test-deconvolution

