@echo off
gcc -Wall -O3 *.c -std=c11 -lm -o ../CorrF 2> error_wingcc.txt & echo 'compiled CorrF by gcc'

