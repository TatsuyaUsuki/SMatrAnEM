@echo off
clang -Wall -O3 *.c -std=c11 -lm -o ..\StrMap 2> error_winclang.txt & echo 'compiled StrMap by clang'
