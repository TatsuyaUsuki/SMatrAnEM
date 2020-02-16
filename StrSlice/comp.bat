@echo off
gcc -Wall -O3 *.c -std=c11 -lm -o ../YeeSlice 2> error_wingcc.txt & echo 'compiled YeeSlice by gcc'

