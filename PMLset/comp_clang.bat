@echo off
clang -Wall -O3 *.c -std=c11 -lm -o ..\PMLgen 2> error_winclang.txt & echo 'compiled PMLgen by clang'
