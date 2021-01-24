@echo off
gcc -Wall -O3 *.c -std=c11 -lm -o ../Med3DSi 2> error_wingcc.txt & echo 'compiled Med3DSi by gcc'

