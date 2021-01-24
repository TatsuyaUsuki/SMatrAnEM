#!/bin/sh
clang -Wall -O3 *.c -std=c11 -lm -o ../Med3DSi 2> error_clang.txt && echo 'compiled Med3DSi by clang'

