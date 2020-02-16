@echo off
gcc -Wall -O3 *.c -std=c11 -lm -o ../StrMap 2> error_wingcc.txt & echo 'compiled StrMap by gcc'

