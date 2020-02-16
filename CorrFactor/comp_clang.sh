#!/bin/sh
clang -Wall -O2 *.c -std=c11 -llapack -lblas -lm -o ../CorrF 2> error_clang.txt && echo 'compiled CorrF by clang'

