#!/bin/sh
gcc -Wall -O2 *.c -std=c11 -lm -o ../StrMap 2> error_gcc.txt && echo 'compiled StrMap by gcc'

