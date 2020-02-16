@echo off
clang -Wall -O3 *.c -std=c11 -lm -o ..\CorrF 2> error_winclang.txt & echo 'compiled CorrF by clang'
