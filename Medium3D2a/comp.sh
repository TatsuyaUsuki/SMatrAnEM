#!/bin/sh
gcc -Wall -O2 *.c -std=c11 -lm -o ../Med3D2a 2> error_gcc.txt && echo 'compiled Med3D2a by gcc'

