//This header file includes general header files, and it defines minimum required macro. 
//                                           Last updated on 2018/08/03.
#ifndef HEADERS_MACROS_H
#define HEADERS_MACROS_H

#include <float.h>
#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>
//#include <sys/time.h>// This files is only used in POSIX.

// To avoid "memory leak", "dangling pointer" and "double-free"
// https://www.ipa.go.jp/security/awareness/vendor/programmingv2/clanguage.html
#define SAFEFREE(A) {free(A);	A = NULL;}
//#define SAFEFREE(A) {if(A != NULL){free(A);	A = NULL;}}
//void SAFEFREE(void *A){if(A != NULL){free(A);	A = NULL;}}

#endif /* HEADERS_MACROS_H */

