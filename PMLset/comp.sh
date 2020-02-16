#!/bin/sh
#gcc -Wall -O3 *.c -std=c11 -llapack -lblas -lgfortran -lm -o ../PMLgen 2> error_gcc.txt && echo 'compiled PMLgen by gcc'
gcc -Wall -O2 *.c -std=c11 -llapack -lblas -lm -o ../PMLgen 2> error_gcc.txt && echo 'compiled PMLgen by gcc by gcc'

