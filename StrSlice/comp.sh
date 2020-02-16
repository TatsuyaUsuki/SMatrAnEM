#!/bin/sh
gcc -Wall -O2 *.c -std=c11 -lm -o ../YeeSlice 2> error_gcc.txt && echo 'compiled YeeSlice by gcc'

