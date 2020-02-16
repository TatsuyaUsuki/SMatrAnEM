#!/bin/sh
gcc -Wall -O2 *.c -std=c11 -lm -o ../xi2u 2> error_gcc.txt && echo 'compiled xi2u by gcc'

