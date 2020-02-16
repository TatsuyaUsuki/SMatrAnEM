@echo off
clang -Wall -O3 *.c -std=c11 -lm -o ..\Mcalc 2> error_winclang.txt & echo 'compiled Mcalc by clang'
