#!/bin/sh
clang -Wall -O2 *.c -std=c11 -llapack -lblas -lm -o ../PMLgen 2> error_clang.txt && echo 'compiled PMLgen by clang'

