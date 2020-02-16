#!/bin/sh
clang -Wall -O3 *.c -std=c11 -lm -o ../YeeSlice 2> error_clang.txt && echo 'compiled YeeSlice by clang'

