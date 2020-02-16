@echo off
clang -Wall -O3 *.c -std=c11 -lm -o ..\Scalc 2> error_winclang.txt & echo 'compiled Scalc by clang'
