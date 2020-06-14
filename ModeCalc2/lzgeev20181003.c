//------------------------------------------------------------------------------
//
//  This program is a high-precision solver of eigenproblem:
//                          A*V = V*diag(E),
//                             where A is an asymmetric complex matrix.
//------------------------------------------------------------------------------
//  Programming language is C99.
//
//  Tatsuya Usuki wrote this to compute precise results by using
//  the 80-bit extended precision on the x86 architecture. 
//  On i7-4770 and gcc 5.4.0 of x86_64-linux-gnu,
//  this can also output the 80-bit extended precision.
//  This is not optimized for speed.
//  When you need speed and efficiency rather than precision,
//  it is strongly recommended to use the "zgeev" in LAPACK.
//------------------------------------------------------------------------------
//  Program for extended precision:   long err = lzgeev_(N,A,E,V)
//  
//  Input matrix `A' changes in this program.
//
//  const long    N   : dimension of matrix, that is grater than 2,
//  long double complex  A  : N by N matrix,
//  long double complex  V  : N by N matrix, 
//  long double complex  E  : N by 1 matrix.
//
//  err = 0 : normal end,
//  err = -1:   N < 3,
//  err > 0 : it did not converge within `err' iterations.
//------------------------------------------------------------------------------
//  Please let me know any bug that you find.
//                                             e-mail: usuki-tty@smatran.org
//                                             URL   : http://www.smatran.org/
//
//                                                Last updated on Oct 03, 2018
//------------------------------------------------------------------------------
#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#define HDQR_EXECTIME 0 // When HDQR_EXECTIME = 1,
                        //   this outputs exec time at each step as stderr.
//------------------------------------------------------------------------------
//
//                      lzgeev_ -+-- Hessenberg_z
//                               |
//                               +-- QR_z <--> *
//                               |
//                               +-- Norm_eigenvector
//
//------------------------------------------------------------------------------
void Hessenberg_z(const long n, long double complex A[n*n],
                     long double complex V[n*n]);
long QR_z(const long n, long double complex A[n*n], long double complex E[n],
            long double complex V[n*n]);
void Norm_comp_eigenvector(const long N, long double complex V[N*N]);
//------------------------------------------------------------------------------
//  Main routine of the eigenproblem solver.
//                                                 Last updated on Aug 23, 2016
//------------------------------------------------------------------------------
long lzgeev_(const long N, long double complex A[N*N], long double complex E[N],
             long double complex V[N*N])
{
	if(N <= 2 ) return(-1);
  for(long j = 0 ; j < N ; j++){//initialize eigenvector matrix V
    for(long i = 0 ; i < N ; i++){V[i+j*N] = 0.L;}
    V[j+j*N] = 1.L;
  }
	time_t timer_s, timer_f;
	double t_span;
	
	timer_s = time(0);
  Hessenberg_z(N,A,V);
	timer_f = time(0);
	t_span = difftime(timer_f,timer_s);
	if(HDQR_EXECTIME == 1) {
    fprintf(stderr,"lzgeev: Hessenberg form. exec time= %1.0f sec\n", t_span);
  }
	timer_s = time(0);
  long ierr = QR_z(N,A,E,V);
	timer_f = time(0);
	t_span = difftime(timer_f,timer_s);
	if(HDQR_EXECTIME == 1) {
    fprintf(stderr,"lzgeev: QR method exec time= %1.0f sec\n", t_span);
  }
  timer_s = time(0);
  Norm_comp_eigenvector(N, V);
  timer_f = time(0);
  t_span = difftime(timer_f,timer_s);
  if(HDQR_EXECTIME == 1) fprintf(stderr,"Normalization of eigenvector exec time= %1.0f sec\n", t_span);
	return(ierr);
}
//------------------------------------------------------------------------------
//                        Normalization of eigenvectors.
//                         Last updated on Oct 04, 2014
//------------------------------------------------------------------------------
void Norm_comp_eigenvector(const long N, long double complex V[N*N])//revised on 20181004
{
	for(long h0 = 0 ; h0 < N ; h0++) {//	fprintf(stderr,"Order= %ld \n", h0);
		long double sum0 = 0.L;
		for(long j0 = 0 ; j0 < N ; j0++) {
			sum0 += creal(V[j0+N*h0])*creal(V[j0+N*h0]);
			sum0 += cimag(V[j0+N*h0])*cimag(V[j0+N*h0]);
		}
		sum0 = 1.L/sqrt(sum0);
		for(long j0 = 0 ; j0 < N ; j0++) {
			V[j0+N*h0] *= sum0;
		}
	}
}
//------------------------------------------------------------------------------
//                        "Hessenberg_z" and "QR_z"
//
//                                               Last updated on Aug 23, 2016
//
//     These functions are translated and simplied from the algorithm of 
//     "corth.f", "comqr2.f" and "cbabk.f" in EISPACK. If you need details,
//     see these FORTRAN files @ http://www.netlib.org/eispack/
//------------------------------------------------------------------------------
void Hessenberg_z(const long n, long double complex A[n*n],
                     long double complex V[n*n])
{// N > 2
  long double complex *cv;
	cv = malloc(sizeof(long double complex)*n);
  for(long m = 1 ; m <= n - 2 ; m++){
    long double h = 0.L;
    cv[m] = 0.L;
    long double norm = 0.L;
    for(long i = m ; i < n ; i++){
      norm += fabs(creal(A[i+(m-1)*n])) + fabs(cimag(A[i+(m-1)*n]));
    }
    if (fabs(norm) > DBL_MIN) {
      for(long i = n - 1 ; i >= m ; i--){
        cv[i] = (A[i+(m-1)*n]) / norm;
        h += creal(cv[i]) * creal(cv[i]) + cimag(cv[i]) * cimag(cv[i]);
      }
      long double g = sqrt(h);
      long double f = fabs(cv[m]);
      if(fabs(f) > DBL_MIN){
        h += f * g;
        g /= f;
        cv[m] *= (1.L + g);
      }else{
        cv[m] = g+I*cimag(cv[m]);
        A[m+(m-1)*n] = norm + I*cimag(A[m+(m-1)*n]);
      }
      for(long j = m ; j < n ; j++){
        long double complex fz = 0.L;
        for(long i = n - 1 ; i >= m ; i--){fz += (conj(cv[i])) * (A[i+j*n]);}
        fz /= h;
        for(long i = m ; i < n ; i++){A[i+j*n] -= fz * cv[i];}
      }
      for(long i = 0 ; i < n ; i++){
        long double complex fz = 0.L;
        for(long j = n - 1 ; j >= m ; j--){fz += (cv[j])*(A[i+j*n]);}
        fz /= h;
        for(long j = m ; j < n ; j++){A[i+j*n] -=  fz * conj(cv[j]);}
      }
      cv[m] *= norm;
      A[m+(m-1)*n] *= -g;
    }
  }
  for(long i = n - 2 ; i >= 1 ; i--){
    if((fabs(creal(cv[i])) > DBL_MIN || fabs(cimag(cv[i])) > DBL_MIN)&&
       (fabs(creal(A[i+(i-1)*n])) > DBL_MIN
        || fabs(cimag(A[i+(i-1)*n])) > DBL_MIN)){
      long double norm = creal(A[i+(i-1)*n]) * creal(cv[i])
                         + cimag(A[i+(i-1)*n]) * cimag(cv[i]);
      for(long j = i + 1 ; j < n ; j++){cv[j] = A[j+(i-1)*n];}
      for(long j = i ; j < n ; j++){
        long double complex sz = 0.L;
        for(long k = i ; k < n ; k++){sz += conj(cv[k]) * V[k+j*n];}
        sz /= norm;
        for(long k = i ; k < n ; k++){V[k+j*n] += sz * cv[k];}
      }
    }
  }
  free(cv);//CAUTION! Never remove it.
  for(long i = 1 ; i < n ; i++){
    long ll = i+1 ; if(ll >= n) ll = n-1;
    if(fabs(cimag(A[i+(i-1)*n])) > DBL_MIN) {
      long double norm = fabs(A[i+(i-1)*n]);
      long double complex yz = (A[i+(i-1)*n]) / norm;
      A[i+(i-1)*n] = norm;
      for(long j = i ; j < n ; j++){A[i+j*n] = conj(yz) * A[i+j*n];}
      for(long j = 0 ; j <= ll ; j++){A[j+i*n] =  yz*A[j+i*n];}
      for(long j = 0 ; j < n ; j++){V[j+i*n] = yz * V[j+i*n];}
    }
  }
}
//------------------------------------------------------------------------------
//
//           QR_z --+-- zero_check
//                  |
//                  +-- r_cal --+-- re_tr
//                  |           |
//                  |           +-- inv_ope
//                  |
//                  +-- gen_vec
//
//------------------------------------------------------------------------------
long zero_check(const long n, const long k, const long double complex A[n*n],
                long *ll);
long double complex r_cal(const long n, const long k, const long l, 
                          const long ip, long double complex A[n*n], 
                          long double complex E[n], long double complex V[n*n]);
long double complex re_tr(const long n, const long l, const long k,
                          long double complex A[n*n], long double complex E[n]);
void inv_ope(const long n, const long k, const long l, 
             const long double complex sz, long double complex A[n*n], 
             const long double complex E[n], long double complex V[n*n]);
void gen_vec(const long n, long double complex A[n*n], long double complex E[n], 
             long double complex V[n*n]);
//------------------------------------------------------------------------------
long QR_z(const long n, long double complex A[n*n], long double complex E[n],
            long double complex V[n*n])
{
  long double complex tz = 0.L;
  long im = 30*n;
  for(long k = n-1 ; k >= 0 ; k--) {
    long ip = 0;
    long l;
    while(zero_check(n, k, A, &l) != 0) {
      if(im == 0){return(k);}//All eigenvalues have not converged. 
      tz += r_cal(n, k, l, ip, A, E, V);
      ip++;
      im--;
    }
    A[k+k*n] += tz;
    E[k] = A[k+k*n];// A root found.
  }
  gen_vec(n,A,E,V);
  return(0);
}
long zero_check(const long n, const long k, const long double complex A[n*n],
                long *ll)
{
  long l;
  for(l = k ; l > 0 ; l--){
    long double a = fabs(creal(A[l-1+(l-1)*n])) + fabs(cimag(A[l-1+(l-1)*n]))
                    + fabs(creal(A[l+l*n])) + fabs(cimag(A[l+l*n]));
    long double b = a + fabs(creal(A[l+(l-1)*n]));
    if(fabs(b - a) <= DBL_MIN) goto l_break;
  }
  l_break:;
  *ll = l;
  if(l == k){return(0);} else {return(1);}
}
long double complex r_cal(const long n, const long k, const long l, 
                          const long ip, long double complex A[n*n], 
                          long double complex E[n], long double complex V[n*n])
{//Root calculation.
  long double complex sz;
  if (ip == 10 || ip == 20) {
    sz = fabs(creal(A[k+(k-1)*n])) + fabs(creal(A[(k-1)+(k-2)*n]));
  }else{
    sz = A[k+k*n];
    long double complex xz = A[(k-1)+k*n] * creal(A[k+(k-1)*n]);
    if(fabs(creal(xz)) > DBL_MIN || fabs(cimag(xz)) > DBL_MIN) {
      long double complex yz = (A[(k-1)+(k-1)*n] - sz) / 2.L;
      long double complex zz = sqrt(yz*yz+xz);
      if(creal(yz) * creal(zz) + cimag(yz) * cimag(zz) < 0.L){zz = -zz;}
      xz /= yz+zz;
      sz -= xz;
    }
  }  
  for(long i = 0 ; i <= k ; i++){A[i+i*n] -= sz;}
  long double complex sz1 = re_tr(n, l, k, A, E);
  inv_ope(n, k, l, sz1, A, E, V);
  return(sz);
}
long double complex re_tr(const long n, const long l, const long k,
                          long double complex A[n*n], long double complex E[n])
{//Reducing to triangle form.
  long double sr = 0.L;//For warning of gcc, it was revised on 20181003
  for(long i = l + 1 ; i <= k ; i++){
    sr = creal(A[i+(i-1)*n]);
    A[i+(i-1)*n] = I*cimag(A[i+(i-1)*n]);
    long double norm = fabs(fabs(A[i-1+(i-1)*n])+I*sr);
    long double complex xz = (A[i-1+(i-1)*n]) / norm;
    E[i-1] = xz;
    A[i-1+(i-1)*n] = norm;
    A[i+(i-1)*n] = I*sr / norm;
    for(long j = i ; j < n ; j++){
      long double complex yz = A[i-1+j*n];
      long double complex zz = A[i+j*n];
      A[i-1+j*n] = conj(xz) * yz + cimag(A[i+(i-1)*n]) * zz;
      A[i+j*n] = xz * zz - cimag(A[i+(i-1)*n]) * yz;
    }
  }
  long double si = cimag(A[k+k*n]);
  long double complex sz = sr+I*si;
  if (fabs(si) > DBL_MIN) {
    long double norm = fabs(creal(A[k+k*n])+I*si);
    sz = (creal(A[k+k*n])+I*si) / norm;
    A[k+k*n] = norm;
    if (k != n-1) {
      for(long j = k+1 ; j < n ; j++){A[k+j*n] *= conj(sz);}
    }
  }
  return(sz);
}
void inv_ope(const long n, const long k, const long l, 
             const long double complex sz, long double complex A[n*n], 
             const long double complex E[n], long double complex V[n*n])
{// Inverse operation.
  for(long j = l + 1 ; j <= k ; j++){
    long double complex xz = E[j-1];
    for(long i = 0 ; i <= j ; i++){
      long double complex yz = creal(A[i+(j-1)*n]);
      long double complex zz = A[i+j*n];
      if(i != j){
        yz = creal(yz) + I*cimag(A[i+(j-1)*n]);
        A[i+(j-1)*n] = creal(A[i+(j-1)*n])
                       +I*(creal(xz) * cimag(yz) + cimag(xz) * creal(yz) 
                                      + cimag(A[j+(j-1)*n]) * cimag(zz));
      }
      A[i+(j-1)*n] = creal(xz) * creal(yz) - cimag(xz) * cimag(yz) 
                      + cimag(A[j+(j-1)*n]) * creal(zz) 
                      + I*cimag(A[i+(j-1)*n]);
      A[i+j*n] = conj(xz) * zz - cimag(A[j+(j-1)*n]) * yz;
    }
    for(long i = 0 ; i < n ; i++){
      long double complex yz = V[i+(j-1)*n];
      long double complex zz = V[i+j*n];
      V[i+(j-1)*n] = xz*yz + cimag(A[j+(j-1)*n]) * zz;
      V[i+j*n] = conj(xz)*zz - cimag(A[j+(j-1)*n]) * yz;
    }
  }
  if(fabs(cimag(sz)) > DBL_MIN){
    for(long i = 0 ; i <= k ; i++){A[i+k*n] *= sz;}
    for(long i = 0 ; i < n ; i++){V[i+k*n] *= sz;}
  }  
}
void gen_vec(const long n, long double complex A[n*n], long double complex E[n], 
             long double complex V[n*n])
{
  long double norm = 0.L;
  for(long i = 0 ; i < n ; i++){
    for(long j = i ; j < n ; j++){
      long double tr = fabs(creal(A[i+j*n])) + fabs(cimag(A[i+j*n]));
      if (tr > norm) norm = tr;
    }
  }
  if (fabs(norm) <= DBL_MIN) return;
  for(long k = n - 1 ; k >= 1 ; k--) {
    long double complex xz = E[k];
    A[k+k*n] = 1.L;
    for(long i = k - 1 ; i >= 0 ; i--){
      long double t1,t2;
      long double complex zz = 0.L;
      for(long j = i + 1 ; j <= k ; j++){zz += (A[i+j*n]) * (A[j+k*n]);}
      long double complex yz = xz - E[i];
      if (fabs(creal(yz)) <= DBL_MIN && fabs(cimag(yz)) <= DBL_MIN) {
        t1 = norm;
        yz = t1 + I*cimag(yz);
        do{
          yz = 0.01L * creal(yz) + I*cimag(yz);
          t2 = norm + creal(yz);
        }while(t2 > t1);
      }
      A[i+k*n] = zz/yz;

      long double tr = fabs(creal(A[i+k*n])) + fabs(cimag(A[i+k*n]));
      if (fabs(tr) > DBL_MIN) {// Overflow control
        t1 = tr;
        t2 = t1 + 1.L/t1;
        if (t2 <= t1) {
          for(long j = i ; j <= k ; j++){A[j+k*n] /= tr;}
        }
      }
    }
  }
  for(long j = n - 1 ; j >= 0 ; j--){
    for(long i = 0 ; i < n ; i++){
      long double complex zz = 0.L;
      for(long k = 0 ; k <= j ; k++){zz += V[i+k*n] * A[k+j*n];}
      V[i+j*n] = zz;
    }// Multiplication by transformation matrix to give all vectors.
  }
}
