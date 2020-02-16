#!/bin/sh
gcc -Wall -O2 *.c -std=c11 -lm -o ../outer 2> error_gcc.txt && echo 'compiled outer by gcc'

