#!/bin/sh
gcc -Wall -O2 *.c -std=c11 -lm -o ../mkCurv 2> error_gcc.txt && echo 'compiled mkCurv by gcc'

