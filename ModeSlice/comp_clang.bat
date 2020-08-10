@echo off
clang -Wall -O3 *.c -std=c11 -lm -o ..\Mcalc2 2> error_winclang.txt & echo 'compiled Mcalc2 by clang'
