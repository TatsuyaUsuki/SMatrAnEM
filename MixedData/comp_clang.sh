#!/bin/sh
clang -Wall -O3 *.c -std=c11 -lm -o ../BindData 2> error_clang.txt && echo 'compiled BindData by clang'

