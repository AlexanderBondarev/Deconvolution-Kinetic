
CC=gcc
CFLAGS=-c

all: deconvolution sect test-deconvolution

deconvolution: src/deconvolution.c
	$(CC) $(CFLAGS) src/deconvolution.c -o bin/deconvolution
	chmod a+x bin/deconvolution

sect: src/sect.c
	$(CC) $(CFLAGS) src/sect.c -o bin/sect
	chmod a+x bin/sect

test-deconvolution: src/test-deconvolution.c
	$(CC) $(CFLAGS) src/test-deconvolution.c -o bin/test-deconvolution
	chmod a+x bin/test-deconvolution

clean:
	rm bin/deconvolution bin/sect bin/test-deconvolution
