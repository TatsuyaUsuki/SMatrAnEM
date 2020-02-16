@echo off
gcc -Wall -O3 *.c -std=c11 -lm -o ../PMLgen 2> error_wingcc.txt & echo 'compiled PMLgen by gcc'

