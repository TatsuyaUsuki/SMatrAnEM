@echo off
gcc -Wall -O3 *.c -std=c11 -lm -o ../u2r2x 2> error_wingcc.txt & echo 'compiled u2r2x by gcc'
#gcc -DNDEBUG -Wall -O3 *.c -std=c11 -lm -o ../u2r2x 2> error_wingcc.txt & echo 'compiled u2r2x by gcc'

