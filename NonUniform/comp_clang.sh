#!/bin/sh
clang -Wall -O2 *.c -std=c11 -lm -o ../xi2u 2> error_clang.txt && echo 'compiled xi2u by clang'

