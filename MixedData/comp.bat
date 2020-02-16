@echo off
gcc -Wall -O3 *.c -std=c11 -lm -o ../BindData 2> error_wingcc.txt & echo 'compiled BindData by gcc'

