/*------------------------------------------------------------------------------
//
//  This program is a high-precision solver of eigenproblem:
//                          A*(VR+I*VI) = (VR+I*VI)*diag((ER+I*EI)),
//                             where A is an asymmetric real matrix.
//------------------------------------------------------------------------------
//  Programming language is C99.
//
//  Tatsuya Usuki wrote this to compute precise results by using
//  the 80-bit extended precision on the x86 architecture. 
//  On i7-4770 and gcc 4.8.2 of x86_64-linux-gnu,
//  this can also output the 80-bit extended precision.
//  This is not optimized for speed.
//  When you need speed and efficiency rather than precision,
//  it is strongly recommended to use the "dgeev" in LAPACK.
//------------------------------------------------------------------------------
//  Program for extended precision:   long err = ldgeev_(N,A,ER,EI,VR)
//  
//  Input matrix `A' changes in this program.
//
//  const long    N   : dimension of matrix, that is grater than 2,
//  long double   A   : N by N matrix,
//  long double   VR  : N by N matrix. VI for complex eigenvalue is
//   									   the next column for real part of VR.
//  long double ER,EI : N by 1 matrix. The complex conjugate eigenvalue
//                           with negative imaginary part sets in the next row.
//
//  err = 0 : normal end,
//  err = -1:   N < 3,
//  err > 0 : it did not converge within `err' iterations.
//------------------------------------------------------------------------------
//  Please let me know any bug that you find.
//                                             e-mail: usuki-tty@smatran.org
//                                             URL   : http://www.smatran.org/
//
//                                                Last updated on Aug 25, 2016
//------------------------------------------------------------------------------
//  References:
//  1) Section 11.10 "Double QR method", pp.219-221,
//    ELMHES(p.345),ELTRAN(p.346),HQR2(p.348),
//    Tsutomu Oguni(Ed.) "Matrix Computing Software (in Japanese)",1991,
//    ISBN 4-621-0654-8
//  2) Section 4.8.6 "Double QR method", pp.151-161,
//    I. Kawakami "Fundmentals of Numerical Computation (in Japanese)",2009,
//    http://www7.ocn.ne.jp/~kawa1/
//  3) "HEQRVS.f", Section 4, pp.1036-1037,
//    Ichizo Ninomiya, Yasuyo Hatano "Mathematical Library NUMPAC",
//    (in Japanese) Information Processing Society of Japan, vol.26,no.9,1985;
//    http://hdl.handle.net/10098/3316  (in English).
//  4) Sections 11.5 "Reduction to Hessenberg Form", pp.478-480,
//    NUMERICAL RECIPES in FORTRAN
//    The Art of Scientific Computing
//    Second Edition
//    W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. Flannery
//    CAMBRIDGE UNIVERSITY PRESS, 1992, ISBN 0 521 43064 X.
//----------------------------------------------------------------------------*/
#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#define HDQR_EXECTIME 0 // When HDQR_EXECTIME = 1, this outputs exec time at each step as stderr.
/*------------------------------------------------------------------------------
//
//                      lgeev_ -+-- Hessenberg_form
//                              |
//                              +-- DoubleQR_LD <--> *
//                              |
//                              +-- Norm_eigenvector
//
//----------------------------------------------------------------------------*/
void Hessenberg_form(const long N, long double A[N*N], long double Q[N*N]);
long DoubleQR_LD(const long N, long double A[N*N], long double E_real[N],
                 long double E_imag[N], long double V_real[N*N]);
void Norm_eigenvector(const long N, const long double E_imag[N], 
                      long double V_real[N*N]);
//------------------------------------------------------------------------------
//  Main routine of the eigenproblem solver.
//                                                 Last updated on Aug 25, 2016
//------------------------------------------------------------------------------
long ldgeev_(const long N, long double A[N*N], long double E_real[N], 
            long double E_imag[N], long double V_real[N*N])
{
	if(N < 3) return(-1);

	for(long i0 = 0 ; i0 < N ; i0++){
		for(long j0 = 0 ; j0 < N ; j0++) {V_real[j0+N*i0] = 0.L;}
		V_real[i0+N*i0] = 1.L;
	}

	time_t timer_s, timer_f;
	double t_span;
	
	timer_s = time(0);
	Hessenberg_form(N, A, V_real);//	Transformation to the Hessenberg form
	timer_f = time(0);
	t_span = difftime(timer_f,timer_s);
	if(HDQR_EXECTIME == 1) fprintf(stderr,"Hessenberg formulation exec time= %1.0f sec\n", t_span);
	
	timer_s = time(0);
	long ierr = DoubleQR_LD(N,A,E_real,E_imag,V_real);// Double QR method
	timer_f = time(0);
	t_span = difftime(timer_f,timer_s);
	if(HDQR_EXECTIME == 1) fprintf(stderr,"Double QR method exec time= %1.0f sec\n", t_span);
	if(ierr != 0) return(ierr);

  timer_s = time(0);
  Norm_eigenvector(N, E_imag, V_real);
  timer_f = time(0);
  t_span = difftime(timer_f,timer_s);
  if(HDQR_EXECTIME == 1) fprintf(stderr,"Normalization of eigenvector exec time= %1.0f sec\n", t_span);
	return(0);
}
//------------------------------------------------------------------------------
//                      Transformation to Hessenberg form.
//                     See section 4.5.2 in p.125, problem 4-6 in p.165,
//  I. Kawakami "Fundmentals of Numerical Computation (in Japanese)",2009,
//                                         http://www7.ocn.ne.jp/~kawa1/
//                                          Last updated on Jul 30, 2014
//------------------------------------------------------------------------------
void Hessenberg_form(const long N, long double A[N*N], long double V[N*N])
{
	long double *Qj0;
	Qj0 = malloc(sizeof(long double)*N);
	long double *Qj1;
	Qj1 = malloc(sizeof(long double)*N);
	long double *W;
	W = malloc(sizeof(long double)*N);

	for(long j0 = 0 ; j0 < N-2 ; j0++){
		long double tt = 0.L;
		for(long i0 = N-1 ; i0 >= j0+2 ; i0--) {tt += A[i0+N*j0]*A[i0+N*j0];}
		if(tt != 0.L){
			tt += A[j0+1+N*j0]*A[j0+1+N*j0];
			long double s0  = copysign(sqrt(tt),A[j0+1+N*j0]);
			long double c2  = 1.L/(tt + A[j0+1+N*j0]*s0);
			for(long i0 = j0+1 ; i0 < N ; i0++)	W[i0] = A[i0+N*j0];
			W[j0+1] += s0;
			for(long i0 = 0 ; i0 < N ; i0++){
				Qj0[i0] = 0.L;
				for(long k0 = j0+1 ; k0 < N ; k0++)	{Qj0[i0] += A[i0+N*k0] * W[k0];}
				Qj0[i0] *= c2;
			}
			for(long i0 = j0 ; i0 < N ; i0++){
				Qj1[i0] = 0.L;
				for(long k0 = j0+1 ; k0 < N ; k0++)	{Qj1[i0] += A[k0+N*i0] * W[k0];}
				Qj1[i0] *= c2;
			}	
			long double cwAw = 0.L;
			for(long i0 = j0+1 ; i0 < N ; i0++)	{cwAw += Qj0[i0] * W[i0];}
			cwAw *= ldexp(c2, -1);
			
			for(long i0 = j0+1 ; i0 < N ; i0++){
				Qj0[i0] -= cwAw*W[i0];
				Qj1[i0] -= cwAw*W[i0];
			}
			for(long i0 = 0 ; i0 < N ; i0++){
				for(long k0 = j0+1 ; k0 < N ; k0++)	{A[i0+N*k0] -= Qj0[i0] * W[k0];}
			}
			for(long k0 = j0 ; k0 < N ; k0++) {A[j0+1+N*k0] -= W[j0+1] * Qj1[k0];}
			for(long i0 = j0+2 ; i0 < N ; i0++){
				A[i0+N*j0] = 0.L;
				for(long k0 = j0+1 ; k0 < N ; k0++) {A[i0+N*k0] -= W[i0] * Qj1[k0];}
			}
// Matrix V accumulation			
			for(long i0 = 0 ; i0 < N ; i0++){
				Qj1[i0] = 0.L;// Qj1 is used as temporary vector.
				for(long k0 = j0+1 ; k0 < N ; k0++)	{Qj1[i0] += V[i0+N*k0] * W[k0];}
				Qj1[i0] *= c2;
			}
			for(long i0 = 0 ; i0 < N ; i0++){
				for(long k0 = j0+1 ; k0 < N ; k0++)	{V[i0+N*k0] -= Qj1[i0] * W[k0];}
			}
		}
	}
		free(Qj0);
		free(Qj1);
		free(W);     
}
//------------------------------------------------------------------------------
//                        Normalization of eigenvectors.
//                         Last updated on Aug 20, 2014
//------------------------------------------------------------------------------
void Norm_eigenvector(const long N, const long double E_imag[N], long double V_real[N*N])
{
	for(long h0 = 0 ; h0 < N ; h0++) {//	fprintf(stderr,"Order= %ld \n", h0);
		if(fabs(E_imag[h0]) < LDBL_MIN) {// Normalization of eigenvector
			long double sum0 = 0.L;
			for(long j0 = 0 ; j0 < N ; j0++) {
				sum0 += V_real[j0+N*h0]*V_real[j0+N*h0];
			}
			sum0 = 1.L/sqrt(sum0);
			for(long j0 = 0 ; j0 < N ; j0++) {
				V_real[j0+N*h0] *= sum0;
			}
		}else{// Normalization of eigenvector
			long double sum0 = 0.L;
			for(long j0 = 0 ; j0 < N ; j0++) {
				sum0 += V_real[j0+N*h0]*V_real[j0+N*h0] + V_real[j0+N*(h0+1)]*V_real[j0+N*(h0+1)];
			}
			sum0 = 1.L/sqrt(sum0);
			for(long j0 = 0 ; j0 < N ; j0++) {
				V_real[j0+N*h0] *= sum0;
				V_real[j0+N*(h0+1)] *= sum0;
			}
			h0++;
		}
	}
}
/*------------------------------------------------------------------------------
//                    Double QR method for long double precision
//
//                                               Last updated on Aug 25, 2016
//
//  This function is translated and modified from FORTRAN procedure "hqr2.f" 
//  which was translated by
//  Burton S. Garbow in Aug. 1983 at Math. Comp. Sci. div, Argonne National Lab.
//  The original hqr2 (ALGOL procedure) was from
//  G. Peters and J. H. Wilkinson, Numr. Math. vol.16 pp.181-204(1970);   
//  Volume. II: Linear Algebra, Handbook for Automatic Computation,
//                           pp.372-395(Springer-Verlag, New York 1971).
//  See details of "hqr.f" @ http://www.netlib.no/netlib/seispack/hqr2.f
------------------------------------------------------------------------------*/
/*
//  Main routine  <--> DoubleQR_LD --+-------------------- lmin                                          
//                                   |                                                                   
//                                   +-- real_vector                                                     
//                                   |                                                                   
//                                   +-- complex_vector -- cdiv                                          
//                                   |                                                                   
//                                   +-- vectors_fullmatrix                                              
*/
void real_vector(const long N, const long en, const long double norm, 
                 const long double E_real[N], const long double E_imag[N], 
                 long double A[N*N]);
void complex_vector(const long N, const long en, const long double norm, 
                    const long double E_real[N], const long double E_imag[N], 
                    long double A[N*N]);
void vectors_fullmatrix(const long N, const long double A[N*N], 
                        long double V_real[N*N]);
/*inline*/ long lmin(const long x, const long y);//modified on 20160825
#define ITER_MAX 30
long DoubleQR_LD(const long N, long double A[N*N], long double E_real[N], 
                 long double E_imag[N], long double V_real[N*N])
{
	long double norm = 0.L;
	{// store roots isolated by balancing routine and compute matrix norm
		long k0 = 0;
		for(long i0 = 0 ; i0 < N ; i0++){
			for(long j0 = k0 ; j0 < N ; j0++){norm += fabs(A[i0+N*j0]);}
			k0 = i0;
		}
	}
   
	long en = N-1;
	long double t = 0.L;
	while(en >= 0){//en loop: search for next eigenvalues
		long iter = 0;
		long na = en - 1;
		long enm2 = na - 1;
		while(iter <= ITER_MAX){
			long l0;
			for(l0 = en  ; l0 >= 0 ; l0--){// look for single small sub-diagonal element for l0=en step -1 until 0
				if (l0 == 0) break;
				long double s = fabs(A[l0-1 + N*(l0-1)]) + fabs(A[l0+N*l0]);
				if (s == 0.L) s = norm;
//				if(fabs(A[l0+N*(l0-1)]) < LDBL_EPSILON * s) break;
				if(fabs(s + A[l0+N*(l0-1)]) == s) break;
			}
//  	   .......... form shift ..........
			long double x = A[en+N*en];
			if (l0 == en){// one root found
				A[en + N*en] = x + t;
				E_real[en] = A[en + N*en];
				E_imag[en] = 0.L;
				en = na;
				goto en_loop_end;
			}
			long double y = A[na+N*na];
			long double w = A[en+N*na] * A[na+N*en];
			long double p = 0.L;
			long double q = 0.L;
			long double r = 0.L;
			if (l0 == na){// two roots found
				p = ldexp(y - x,-1);// / 2.0L;
				q = p * p + w;
				long double zz = sqrt(fabs(q));
				A[en + N*en] = x + t;
				x = A[en + N*en];
				A[na + N*na] = y + t;
				if (q >= 0.L){// real pair
					zz = p + copysign(zz, p);
					E_real[na] = x + zz;
					E_real[en] = E_real[na];
					if (zz != 0.L) E_real[en] = x - w / zz;
					E_imag[na] = 0.L;
					E_imag[en] = 0.L;
					x = A[en + N*na];
					long double s = fabs(x) + fabs(zz);
					p = x / s;
					q = zz / s;
					r = sqrt(p*p+q*q);
					p /= r;
					q /= r;
					for(long j0 = na ; j0 < N ; j0++){// row modification
						zz = A[na + N*j0];
						A[na + N*j0] = q * zz + p * A[en + N*j0];
						A[en + N*j0] = q * A[en + N*j0] - p * zz;
					}
					for(long i0 = 0 ; i0 <= en ; i0++){// column modification
						zz = A[i0 + N*na];
						A[i0 + N*na] = q * zz + p * A[i0 + N*en];
						A[i0 + N*en] = q * A[i0 + N*en] - p * zz;
					}
					for(long i0 = 0 ; i0 < N ; i0++){// accumulate transformations
						zz = V_real[i0 + N*na];
						V_real[i0 + N*na] = q * zz + p * V_real[i0 + N*en];
						V_real[i0 + N*en] = q * V_real[i0 + N*en] - p * zz;
					}
				}else{// complex pair
					E_real[na] = x + p;
					E_real[en] = x + p;
					E_imag[na] = zz;
					E_imag[en] = -zz;
				}
				en = enm2;
				goto en_loop_end;
			}
			if (iter == ITER_MAX) return(ITER_MAX);//set error: An eigenvalue has not converged after ITER_MAX iterations
			if (iter > 0 && iter % 10 == 0){// form exceptional shift
				t += x;
				for(long i0 = 0 ; i0 <= en ; i0++) A[i0 + N*i0] -= x;
	
				long double s = fabs(A[en + N*na]) + fabs(A[na + N*enm2]);
				x = 0.75L * s;
				y = x;
				w = -0.4375L * s * s;
			}

			iter++;
			long m0;
			for(m0 = enm2 ; m0 >= l0 ; m0--){// look for two consecutive small sub-diagonal elements.
				long double zz = A[m0 + N*m0];
				r = x - zz;
				long double s = y - zz;
				p = (r * s - w) / A[m0+1 + N*m0] + A[m0 + N*(m0+1)];
				q = A[m0+1 + N*(m0+1)] - zz - r - s;
				r = A[m0+2 + N*(m0+1)];
				s = fabs(p) + fabs(q) + fabs(r);
				p /= s;
				q /= s;
				r /= s;
				if (m0 == l0) break;
				long double tst1 = fabs(p)*(fabs(A[m0-1 + N*(m0-1)]) + fabs(zz) + fabs(A[m0+1 + N*(m0+1)]));
//				if(fabs(A[m0 + N*(m0-1)])*(fabs(q) + fabs(r)) < LDBL_EPSILON * tst1) break;
				if(fabs(A[m0 + N*(m0-1)])*(fabs(q) + fabs(r)) + tst1 == tst1) break;
			}
	
			for(long i0 = m0+2 ; i0 <= en ; i0++){
		         A[i0 + N*(i0-2)] = 0.L;
		         if (i0 != m0+2) A[i0 + N*(i0-3)] = 0.L;
			}
			for(long k0 = m0 ; k0 <= na ; k0++){
				if (k0 != m0){
					p = A[k0 + N*(k0-1)];
					q = A[k0+1 + N*(k0-1)];
					r = 0.L;
					if (k0 != na) r = A[k0+2 + N*(k0-1)];
					x = fabs(p) + fabs(q) + fabs(r);
					if (x == 0.L) goto pqr_norm1_eq_0;
					p /= x;
					q /= x;
					r /= x;
				}
				long double s = copysign(sqrt(p*p+q*q+r*r),p);
				if (k0 != m0){
					A[k0 + N*(k0-1)] = -s * x;
				}else{
					if (l0 != m0) A[k0 + N*(k0-1)] = -A[k0 + N*(k0-1)];
				}
				p += s;
				x = p / s;
				y = q / s;
				long double zz = r / s;
				q /= p;
				r /= p;
				if (k0 == na){
					for(long j0 = k0 ; j0 < N ; j0++){// row modification @ k0 == na
						p = A[k0 + N*j0] + q * A[k0+1 + N*j0];
						A[k0 + N*j0] -= p * x;
						A[k0+1 + N*j0] -= p * y;
					}
					for(long i0 = 0 ; i0 <= lmin(en,k0+3) ; i0++){// column modification @ k0 == na
						p = x * A[i0 + N*k0] + y * A[i0 + N*(k0+1)];
						A[i0 + N*k0] -= p;
						A[i0 + N*(k0+1)] -= p * q;
			        }
			        for(long i0 = 0 ; i0 < N ; i0++){// accumulate transformations @ k0 == na
						p = x * V_real[i0 + N*k0] + y * V_real[i0 + N*(k0+1)];
						V_real[i0 + N*k0] -= p;
						V_real[i0 + N*(k0+1)] -= p * q;
					}
				}else{
					for(long j0 = k0 ; j0 < N ; j0++){// row modification @ k0 != na
						p = A[k0 + N*j0] + q * A[k0+1 + N*j0] + r * A[k0+2 + N*j0];
						A[k0 + N*j0] -= p * x;
						A[k0+1 + N*j0] -= p * y;
						A[k0+2 + N*j0] -= p * zz;
					}
					for(long i0 = 0 ; i0 <= lmin(en,k0+3) ; i0++){// column modification @ k0 != na
						p = x * A[i0 + N*k0] + y * A[i0 + N*(k0+1)] + zz * A[i0 + N*(k0+2)];
						A[i0 + N*k0] -= p;
						A[i0 + N*(k0+1)] -= p * q;
						A[i0 + N*(k0+2)] -= p * r;
					}
					for(long i0 = 0 ; i0 < N ; i0++){// accumulate transformations @ k0 != na
						p = x * V_real[i0 + N*k0] + y * V_real[i0 + N*(k0+1)] + zz * V_real[i0 + N*(k0+2)];
						V_real[i0 + N*k0] -= p;
						V_real[i0 + N*(k0+1)] -= p * q;
						V_real[i0 + N*(k0+2)] -= p * r;
					}
				}
pqr_norm1_eq_0:;
			}
		}
en_loop_end:;
	}// .. all roots found.
	if (norm != 0.L){
		for(long en = N-1 ; en >= 0 ; en--){// backsubstitute to find vectors of upper triangular form
			if(E_imag[en] == 0.L) real_vector(N, en, norm, E_real, E_imag, A);
			if(E_imag[en] < 0.L)  complex_vector(N, en, norm, E_real, E_imag, A);
		}//	 end of back substitution.
		vectors_fullmatrix(N, A, V_real);
	}
	return(0);
}
/*inline*/ void cdiv(const long double AR, const long double a_imag, 
                 const long double b_real, const long double b_imag, 
                 long double *c_real, long double *c_imag);//modified on 20160825
void complex_vector(const long N, const long en, const long double norm, 
                    const long double E_real[N], const long double E_imag[N], 
                    long double A[N*N])
{
	long double p = E_real[en];
	long double q = E_imag[en];
	long na = en - 1;
	long m = na;
	if (fabs(A[en + N*na]) > fabs(A[na + N*en])){
		A[na + N*na] = q / A[en + N*na];
		A[na + N*en] = -(A[en + N*en] - p) / A[en + N*na];
	}else{
		cdiv(0.L,-A[na + N*en],A[na + N*na]-p,q,&A[na + N*na],&A[na + N*en]);
	}
	A[en + N*na] = 0.L;
	A[en + N*en] = 1.L;
	long double r = 0.L;
	long double s = 0.L;
	long double zz = 0.L;//revised on Jan 06, 2015
	for(long i0 = na-1 ; i0 >= 0 ; i0--){
		long double w = A[i0 + N*i0] - p;
		long double ra = 0.L;
		long double sa = 0.L;

		for(long j0 = m ; j0 <= en ; j0++){
			ra += A[i0 + N*j0] * A[j0 + N*na];
			sa += A[i0 + N*j0] * A[j0 + N*en];
		}
		if(E_imag[i0] < 0.L){
			zz = w;
			r = ra;
			s = sa;
		}else{
			m = i0;
			if (E_imag[i0] == 0.L) {
				cdiv(-ra,-sa,w,q,&A[i0 + N*na],&A[i0 + N*en]);
			}else{// solve complex equations
				long double x = A[i0 + N*(i0+1)];
				long double y = A[i0+1 + N*i0];
				long double vr = (E_real[i0] - p) * (E_real[i0] - p) 
                          + E_imag[i0] * E_imag[i0] - q * q;
				long double vi = ldexp((E_real[i0] - p)*q, 1);
				if (vr == 0.L && vi == 0.L){
					vr = LDBL_EPSILON * norm * 
                (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));
				}
				cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,&A[i0 + N*na],&A[i0 + N*en]);
				if (fabs(x) > fabs(zz) + fabs(q)){
					A[i0+1 + N*na] = (-ra - w * A[i0 + N*na] + q * A[i0 + N*en]) / x;
					A[i0+1 + N*en] = (-sa - w * A[i0 + N*en] - q * A[i0 + N*na]) / x;
				}else{
					cdiv(-r-y*A[i0 + N*na],-s-y*A[i0 + N*en],zz,q,&A[i0+1 + N*na],
                &A[i0+1 + N*en]);
				}
			}// end of solving complex equations

			long double t = fmax(fabs(A[i0 + N*na]), fabs(A[i0 + N*en]));// overflow control
			if (t != 0.L){
				if (1.L/t < t*LDBL_EPSILON){
					for(long j0 = i0 ; j0 <= en ; j0++){
						A[j0 + N*na] /= t;
						A[j0 + N*en] /= t;
					}
				}
			}// end of overflow control
		}
	}
}
void real_vector(const long N, const long en, const long double norm, 
                 const long double E_real[N], const long double E_imag[N], 
                 long double A[N*N])
{
	long double p = E_real[en];
	long m = en;
	A[en + N*en] = 1.L;
	long double s = 0.L;
	long double zz = 0.L;//revised on Jan 06, 2015
	for(long i0 = en-1 ; i0 >= 0 ; i0--){
		long double w = A[i0 + N*i0] - p;
		long double r = 0.L;
		for(long j0 = m ; j0 <= en ; j0++) r += A[i0 + N*j0] * A[j0 + N*en];
		if (E_imag[i0] < 0.L){
			zz = w;
			s = r;
		}else{
			m = i0;
			if(E_imag[i0] == 0.L){
				long double t = w;
				if(t == 0.L) {
					t = LDBL_EPSILON * norm;
				}
				A[i0 + N*en] = -r / t;
			}else{// solve real equations
				long double x = A[i0 + N*(i0+1)];
				long double y = A[i0+1 + N*i0];
				long double q = (E_real[i0] - p) * (E_real[i0] - p)
                          + E_imag[i0] * E_imag[i0];
				long double t = (x * s - zz * r) / q;
				A[i0 + N*en] = t;
				if (fabs(x) <= fabs(zz)) {
					A[i0+1 + N*en] = (-s - y * t) / zz;
				}else{
					A[i0+1 + N*en] = (-r - w * t) / x;
				}
			}// end of solve real equations

			long double t = fabs(A[i0 + N*en]);// overflow control
			if(t != 0.L){
				if (1.L/t < t*LDBL_EPSILON){
//				if(t + 1.L/t <= t) {
					for(long j0 = i0 ; j0 <= en ; j0++) {A[j0 + N*en] /= t;}
				}
			}// end of overflow control
		}
	}
}
void vectors_fullmatrix(const long N, const long double A[N*N], 
                        long double V_real[N*N])
{
	for(long j0 = N-1 ; j0 >= 0 ; j0--){
		for(long i0 = 0 ; i0 < N ; i0++){
			long double sum = 0.L;
			for(long k0 = 0 ; k0 <= j0 ; k0++) sum += V_real[i0 + N*k0] * A[k0 + N*j0];
			V_real[i0 + N*j0] = sum;// multiply by transformation matrix to give vectors of original full matrix.
		}
	}
}
/*inline*/ void cdiv(const long double a_real, const long double a_imag, 
                 const long double b_real, const long double b_imag, 
                 long double *c_real, long double *c_imag)//modified on 20160825
{
      long double norm1 = fabs(b_real) + fabs(b_imag);
      long double a_real_n = a_real/norm1;
      long double a_imag_n = a_imag/norm1;
      long double b_real_n = b_real/norm1;
      long double b_imag_n = b_imag/norm1;
      norm1 = b_real_n*b_real_n + b_imag_n*b_imag_n;
      *c_real = (a_real_n*b_real_n + a_imag_n*b_imag_n)/norm1;
      *c_imag = (a_imag_n*b_real_n - a_real_n*b_imag_n)/norm1;
}
/*inline*/ long lmin(const long x, const long y)//modified on 20160825
{
	if(x < y){
		return(x);
	}else{
		return(y);
	}
}

