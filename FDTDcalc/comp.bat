@echo off
gcc -Wall -O3 *.c -std=c11 -lm -o ../FDTD 2> error_wingcc.txt & echo 'compiled FDTD by gcc'

