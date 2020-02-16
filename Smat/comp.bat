@echo off
gcc -Wall -O3 *.c -std=c11 -lm -o ../Scalc 2> error_wingcc.txt & echo 'compiled Scalc by gcc'

