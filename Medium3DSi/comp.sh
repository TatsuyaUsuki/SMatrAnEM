#!/bin/sh
gcc -Wall -O2 *.c -std=c11 -lm -o ../Med3DSi 2> error_gcc.txt && echo 'compiled Med3DSi by gcc'

