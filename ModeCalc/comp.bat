@echo off
gcc -Wall -O3 *.c -std=c11 -lm -o ../Mcalc 2> error_wingcc.txt & echo 'compiled Mcalc by gcc'

