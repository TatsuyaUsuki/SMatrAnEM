#!/bin/sh
#clang -Wall -O3 *.c -std=c11 -llapack -lblas -lgfortran -lm -o ../Mcalc2 2> error_clang.txt && echo 'compiled Mcalc2 by clang'
#clang -Wall -O3 *.c -std=c11 -llapack -lblas -lm -o ../Mcalc2 2> error_clang.txt && echo 'compiled Mcalc2 by clang'
clang -Wall -O2 *.c -std=c99 -llapack -lblas -lm -o ../Mslice 2> error_clang.txt && echo 'compiled Mslice by clang'

