#!/bin/sh
gcc -Wall -O2 *.c -std=c11 -lm -o ../BindData 2> error_gcc.txt && echo 'compiled BindData by gcc'

