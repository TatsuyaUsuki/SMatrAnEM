#!/bin/sh
#gcc -Wall -O3 *.c -std=c11 -llapack -lblas -lgfortran -lm -o ../Scalc 2> error_gcc.txt && echo 'compiled Scalc'
gcc -Wall -O3 *.c -std=c11 -llapack -lblas -lm -o ../Scalc 2> error_gcc.txt && echo 'compiled Scalc'

