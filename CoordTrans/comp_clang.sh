#!/bin/sh
clang -Wall -O3 *.c -std=c11 -lm -o ../u2r2x 2> error_clang.txt && echo 'compiled u2r2x by clang'
#clang -DNDEBUG -Wall -O3 *.c -std=c11 -lm -o ../u2r2x 2> error_clang.txt && echo 'compiled u2r2x by clang'
