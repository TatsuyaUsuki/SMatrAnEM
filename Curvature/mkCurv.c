//  Start of mkCurv.h
#include "header_macro.h"
#include "constant.h"
struct FnChar{
	char name[BUFSIZE];
	long double sum;
};
//======================================================================
//   main -+-- input_file1, input_file2 -- rm_space                                 
//         +-- input_file3                                  
//         +-- Calc_main -- FnName                                      
//======================================================================
int input_file1(FILE *fp_i1, long L[3], char xi_file[]);
int input_file2(FILE *fp_i1, struct FnChar *F_b, struct FnChar *F_w, char kappa_file[]);
int input_file3(const char xi_file[], const long L[3], long double u[], long double dudxi[]);
void Calc_main(char **argv, const char kappa_file[], const long L_2,
								const long double u[], const long double dudxi[], 
								const struct FnChar *F_b, const struct FnChar *F_w);
//  End of mkCurv.h
//======================================================================
//  Main routine of mkCurv.                                             
//======================================================================
int main(int argc, char **argv)
{
	if(argc != 2) {
		fprintf(stderr,"error: number of files \n");
		exit(1);
	}else if(strncmp(argv[1], "-v", 2) == 0 ) {// This scope was added on Jul 02, 2018
		fprintf(stderr,"The '%s' creates functions for two curvatures.\n", argv[0]);
		fprintf(stderr,"Version 18.08.03 is compiled at %s on %s\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(stderr," Source code : '%s'\n Author      : Tatsuya Usuki\n URL         : http://www.smatran.org\n", __FILE__);
		fprintf(stderr," Reference   : 'Wave scattering in frequency domain' as 'Formulation.pdf' on May 20, 2018\n");
		exit(0);//normal end
	}
	
	FILE *fp_i1;
	fp_i1 = fopen(argv[1],"r");
	if (fp_i1 == NULL){
		fprintf(stderr,"open error!: open input-file!\n");
		exit(1);
	}

	long L[3]; //division number: L0 = L[0], L1 = L[1], L2 = L[2]
	char xi_file[BUFSIZE];
	if(input_file1(fp_i1, L, xi_file) != 0) { 
		fprintf(stderr,"input_file1 error!\n");
		exit(1);
	}
	rewind(fp_i1);
	
	struct FnChar F_b, F_w;
	char kappa_file[BUFSIZE];
	if(input_file2(fp_i1, &F_b, &F_w, kappa_file) != 0) { 
		fprintf(stderr,"input_file2 error!\n");
		exit(1);
	}
	fclose(fp_i1);

	long double *u = calloc((2*L[2]+1),sizeof(long double));//revised on 20180725
	long double *dudxi = calloc((2*L[2]+1),sizeof(long double));//revised on 20180725
	if(u == NULL || dudxi == NULL){fprintf(stderr,"u and/or dudxi can not be secured!\n");	exit(EXIT_FAILURE);}//revised on 20180725
	if(input_file3(xi_file, L, u, dudxi) != 0) { // inputting parameters from u-mesh data 
		fprintf(stderr,"input_file3 error!\n");
		exit(1);
	}
	
	Calc_main(argv, kappa_file, L[2], u, dudxi, &F_b, &F_w); // calc. of EM and outputting results to an outputfile
	SAFEFREE(u);	SAFEFREE(dudxi);
}
//======================================================================
//  input_file3                                                         
//                                       Last updated on May 23, 2018.  
//======================================================================
int input_file3(const char xi_file[], const long L[3], long double u[], long double dudxi[])
{
	FILE *fp_i1;
	fp_i1 = fopen(xi_file,"r");
//	fprintf(stderr,"file name = %s@input_file3\n", xi_file);
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = -(2*L[0] + 1 + 2*L[1] + 1);
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 2*L[2] + 1) { 
		long double xi_dummy, u_dummy, dudxi_dummy;
		if(strncmp(buf, "#", 1) != 0 && sscanf(buf,"%Lf, %Lf, %Lf", &xi_dummy, &u_dummy, &dudxi_dummy) == 3){
//			fprintf(stderr,"buf = %s@input_file3\n", buf);
			if( j_count >= 0  && j_count < 2*L[2] + 1 ) {
				u[j_count] = u_dummy;
				dudxi[j_count] = dudxi_dummy;
			}
//			fprintf(stderr,"j_count = %d@input_file3\n", j_count);
			j_count++;
		}
	}
	fclose(fp_i1);
	return(j_count - (2*L[2] + 1));
}
//======================================================================
//  input_file1                                                         
//                           Last updated on Jun 16, 2018.  
//======================================================================
void rm_space( char *A );
int input_file1(FILE *fp_i1, long L[3], char xi_file[])
{
	char buf[BUFSIZE];	// buffer for fgets
	for(long j = 0 ; j < 3 ; j++) {L[j] = -1;}
	int j_count = 0;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 1) { 
		rm_space(buf);
		if(strncmp(buf, "xi2u", 4) == 0 && sscanf(buf,"%*[^=] %*[=] %s", xi_file) == 1){
			j_count++;//j_count == 1
		}
	}
	if(j_count != 1) {fprintf(stderr,"xi_file can not read @ input_file1!\n");	exit(1);}
	FILE *fp_ix;
	fp_ix = fopen(xi_file,"r");
	if (fp_ix == NULL){
		fprintf(stderr,"%s can not open!\n", xi_file);
		exit(1);
	}else{
		while(fgets(buf, sizeof( buf ), fp_ix) != NULL && j_count < 4){
			if(sscanf(buf,"%*[^=] %*[=] %ld", &L[j_count - 1]) == 1){	j_count++;}//j_count == 1+3
		}
	}
	fclose(fp_ix);
	return(j_count - 4);
}
//======================================================================
//  input_file2                                                         
//                          Last updated on Jun 25, 2018.  
//======================================================================
int input_file2(FILE *fp_i1, struct FnChar *F_b, struct FnChar *F_w, char kappa_file[])
{
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = -3;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 0) { 
		rm_space(buf);
		if(strncmp(buf, "F_b", 3) == 0 && sscanf(buf,"%*[^=] %*[=] %s %*[^=] %*[=] %Lf", F_b->name, &(F_b->sum)) == 2) {
			j_count++;
		}else if(strncmp(buf, "F_w", 3) == 0 && sscanf(buf,"%*[^=] %*[=] %s %*[^=] %*[=] %Lf", F_w->name, &(F_w->sum)) == 2) {
			j_count++;
		}else if(strncmp(buf, "kappa", 5) == 0 && sscanf(buf,"%*[^=] %*[=] %s", kappa_file) == 1) {
			j_count++;
		}
	}
	return(j_count);
}
//======================================================================
//  This program removes spaces of head in characters,
//                            it needs #include <ctype.h> 
//                          Last updated on Jun 05, 2018. 
//======================================================================
void rm_space(char *A)
{
//	fprintf(stderr,"Read data before removing spaces =%s\n", A);
//	A[BUFSIZE - 1] = '\0';// A[] has to include NULL character.
	while(isspace( A[0] ) != 0){
		int i = 1;
		while(A[i] != '\0') {
			A[i-1] = A[i];
			i++;
		}
		A[i-1] = '\0';
	}
//	fprintf(stderr,"Read data after removing spaces =%s\n", A);
}
//======================================================================
//  This program is Calc_main Program of mkCurv.                        
//                                       Last updated on Jun 22, 2018.  
//======================================================================
long double FnName(const char *S, const long double x);
void Calc_main(char **argv, const char kappa_file[], const long L_2,
								const long double u[], const long double dudxi[], 
								const struct FnChar *F_b, const struct FnChar *F_w)
{
	FILE *fp_o;
	fp_o = fopen(kappa_file,"w");
	if (fp_o == NULL){
		fprintf(stderr,"open error!:  output-test\n");
		exit(1);
	}
	
	time_t timer_s = time(0);
	clock_t startClock = clock();
 /* "sys/time.h" is not ANSI but POSIX.
	struct timeval startTime, endTime;// see http://8ttyan.hatenablog.com/entry/2015/02/03/003428
	gettimeofday(&startTime, NULL);
// end: sys/time.h */
	
	long double *kappa_b = calloc((2*L_2+1),sizeof(long double));//revised on 20180725
	long double *kappa_w = calloc((2*L_2+1),sizeof(long double));//revised on 20180725
	if(kappa_b == NULL || kappa_w == NULL){fprintf(stderr,"kappa_b and/or kappa_w can not be secured!\n");	exit(EXIT_FAILURE);}//revised on 20180725
	
	for(long J0 = 0 ; J0 < 2*L_2+1 ; J0++) { 
		long double x = u[J0]/u[2*L_2];
		kappa_b[J0] = FnName((*F_b).name, x);
		kappa_w[J0] = FnName((*F_w).name, x);
	}

	long double halfsum_b = 0.5L*kappa_b[L_2]*dudxi[L_2];
	long double halfsum_w = 0.5L*kappa_w[L_2]*dudxi[L_2];
	for(long J0 = L_2 + 1 ; J0 < 2*L_2 ; J0++) { 
		halfsum_b += kappa_b[J0]*dudxi[J0];
		halfsum_w += kappa_w[J0]*dudxi[J0];
	}
	halfsum_b += 0.5L*kappa_b[2*L_2]*dudxi[2*L_2];
	halfsum_w += 0.5L*kappa_w[2*L_2]*dudxi[2*L_2];
	long double c_b = Pi*(F_b->sum)/(0.5L*halfsum_b);//Unit of half-sum was set to 'pi' on 2018.06.14
	long double c_w = log(F_w->sum)/(0.5L*halfsum_w);//Note that xi = 0.5*J0 on 2018.06.14
	for(long J0 = 0 ; J0 < 2*L_2+1 ; J0++) { 
		kappa_b[J0] *= c_b;
		kappa_w[J0] *= c_w;
	}
	fprintf(fp_o,"# xi2,             k0*u2,             k0*du2/dxi2,          kappa_b/k0,          kappa_w/k0\n");
	for(long J0 = 0 ; J0 < 2*L_2+1 ; J0++) { 
		fprintf(fp_o,"%.1f, %.20LE, %.20LE, %.20LE, %.20LE\n", 0.5*((double) J0), u[J0], dudxi[J0], kappa_b[J0], kappa_w[J0]);
	}
	SAFEFREE(kappa_b);	SAFEFREE(kappa_w);
	
	time_t timer_f = time(0);
	double t_span = difftime(timer_f,timer_s);
	clock_t endClock = clock();
 /* "sys/time.h" is not ANSI but POSIX.
//	struct timeval endTime;// see http://8ttyan.hatenablog.com/entry/2015/02/03/003428
	gettimeofday(&endTime, NULL);
	time_t diffsec = difftime(endTime.tv_sec,startTime.tv_sec);
	suseconds_t diffsub = endTime.tv_usec - startTime.tv_usec;
	double t_span2 = diffsec + diffsub*1.e-6;
// end: sys/time.h */
	double cpusec = (endClock - startClock)/(double)CLOCKS_PER_SEC;

	fprintf(fp_o,"# Input-file = %s\n", argv[1]);
	fprintf(fp_o,"# This file  = %s\n", kappa_file);
	fprintf(fp_o,"# Exe-file   = %s : `%s' was compiled at %s on %s by C-version:%ld\n", argv[0], __FILE__, __TIME__, __DATE__, __STDC_VERSION__);
	fprintf(fp_o,"# LDBL_DIG   = %u, LDBL_EPSILON = %.20LE\n", __LDBL_DIG__,  __LDBL_EPSILON__);
	fprintf(fp_o,"# Exec time  = %1.0f sec, CPU time = %1.6f sec\n", t_span, cpusec);
	fprintf(fp_o,"# Created on %s", ctime(&timer_f));
	fprintf(stderr,"# Created on %s", ctime(&timer_f));
	fclose(fp_o);
}
//======================================================================
//  FnName                                                           
//                                       Last updated on Jul 02, 2018.  
//======================================================================
long double FnName(const char *S, const long double x)
{
	long double Fn;
	long double y = 0.5L*Pi*x;
	if (strncmp(S, "const_odd", 9) == 0) {
		if(x < 1.L && x > 0.L){
			Fn = 1.L;
		} else if(x > -1. && x < 0.) {
			Fn = -1.L;
		} else {
			Fn = 0.L;
		}
	} else if (strncmp(S, "const_even",10) == 0) {
		if(x < 1.L && x > -1.L){
			Fn = 1.L;
		} else {
			Fn = 0.L;
		}
	} else if (strncmp(S, "const_cos",9) == 0) {
		if(x < 1.L && x > -1.L){
			Fn = cos(y);
		} else {
			Fn = 0.L;
		}
	} else if (strncmp(S, "const_sin",9) == 0) {
		if(x < 1.L && x > -1.L){
			Fn = sin(y);
		} else {
			Fn = 0.L;
		}
	} else if (strncmp(S, "cos_cos",7) == 0) {
		if(x < 1.L && x > -1.L){
			Fn = cos(y)*cos(y);
		} else {
			Fn = 0.L;
		}
	} else if (strncmp(S, "cos_sin",7) == 0 || strncmp(S, "sin_cos", 7) == 0) {
		if(x < 1.L && x > -1.L){
			Fn = cos(y)*sin(y);
		} else {
			Fn = 0.L;
		}
	} else if (strncmp(S, "cos2_sin",7) == 0 || strncmp(S, "sin_cos2", 7) == 0) {
		if(x < 1.L && x > -1.L){
			Fn = cos(y)*cos(y)*sin(y);
		} else {
			Fn = 0.L;
		}
	} else {
		fprintf(stderr,"FnName error!  name= %s\n", S);
		exit(1);
	}
	return(Fn);
}
