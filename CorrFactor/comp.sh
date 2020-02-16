#!/bin/sh
#gcc -Wall -O3 *.c -std=c11 -llapack -lblas -lgfortran -lm -o ../CorrF 2> error_gcc.txt && echo 'compiled CorrF by gcc'
gcc -Wall -O2 *.c -std=c11 -llapack -lblas -lm -o ../CorrF 2> error_gcc.txt && echo 'compiled CorrF by gcc'

