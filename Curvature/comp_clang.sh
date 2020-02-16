#!/bin/sh
clang -Wall -O3 *.c -std=c11 -lm -o ../mkCurv 2> error_clang.txt && echo 'compiled mkCurv by clang'

