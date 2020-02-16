#!/bin/sh
clang -Wall -O2 *.c -std=c11 -llapack -lblas -lm -o ../Scalc 2> error_clang.txt && echo 'compiled Scalc by clang'

