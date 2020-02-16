@echo off
clang -Wall -O3 *.c -std=c11 -lm -o ..\FDTD 2> error_winclang.txt & echo 'compiled FDTD by clang'
