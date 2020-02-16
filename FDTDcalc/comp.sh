#!/bin/sh
#gcc -Wall -O3 *.c -std=c11 -llapack -lblas -lgfortran -lm -o ../FDTD 2> error_gcc.txt && echo 'compiled FDTD by gcc'
gcc -Wall -O2 *.c -std=c11 -llapack -lblas -lm -o ../FDTD 2> error_gcc.txt && echo 'compiled FDTD by gcc'

