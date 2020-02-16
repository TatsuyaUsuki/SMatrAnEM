#!/bin/sh
clang -Wall -O3 *.c -std=c11 -lm -o ../StrMap 2> error_clang.txt && echo 'compiled StrMap by clang'

