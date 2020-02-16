@echo off
clang -Wall -O3 *.c -std=c11 -lm -o ..\YeeSlice 2> error_winclang.txt & echo 'compiled YeeSlice by clang'
