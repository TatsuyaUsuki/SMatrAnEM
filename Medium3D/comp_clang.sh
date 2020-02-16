#!/bin/sh
clang -Wall -O3 *.c -std=c11 -lm -o ../Med3D 2> error_clang.txt && echo 'compiled Med3D by clang'

