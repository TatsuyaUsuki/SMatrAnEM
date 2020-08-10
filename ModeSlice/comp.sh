#!/bin/sh
#gcc -Wall -O3 *.c -std=c11 -llapack -lblas -lgfortran -lm -o ../Mcalc2 2> error_gcc.txt && echo 'compiled Mcalc2 by gcc'
#gcc -Wall -O2 *.c -std=c11 -llapack -lblas -lgfortran -lm -o ../Mcalc2 2> error_gcc.txt && echo 'compiled Mcalc2 by gcc'
#gcc -Wall -O2 *.c -std=gnu11 -llapack -lblas -lm -o ../Mcalc2 2> error_gcc.txt && echo 'compiled Mcalc2 by gcc'
gcc -Wall -O2 *.c -std=c11 -llapack -lblas -lm -o ../Mslice 2> error_gcc.txt && echo 'compiled Mslice by gcc'

