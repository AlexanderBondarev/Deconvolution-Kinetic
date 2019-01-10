#!/bin/bash

name=`basename $1 .R`

Rscript --vanilla $name.R $2 $3 $4 $5 $6 $7 $8 $9 $10 | tee $2\_Concurent.log
##R -q --vanilla <$name.R $2 $3  | tee $name.log
