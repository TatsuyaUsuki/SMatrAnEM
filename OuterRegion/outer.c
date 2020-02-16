//  Start of outer.h
#include "header_macro.h"
#include "constant.h"
int input_file1(FILE *fp_i1, long L[3], char xi_file[]);
int input_file2(FILE *fp_i1, char urxF[], long *cell_num, char outer_file[]);
int input_file3(const char xi_file[], const long L[3], long double u0[], long double u2[], long double dudxi[]);
int input_file4(const char urxF[], const long L[3], long double x[], long double h0[]);
void Calc_main(FILE *fp_o, const long L[3], 
					const long double u0[2*L[0]+1], const long double u2[2*L[2]+1], 
					const long double dudxi[2*L[2]+1], const long double x[2*(2*L[0]+1)*2], 
					const long double h0[2*(2*L[0]+1)], const long cell_num);
//  End of outer.h
//======================================================================
//   main -+-- input_file1, input_file2 -- rm_space                                 
//         +-- input_file3, input_file4
//         +-- Calc_main                                      
//======================================================================
int main(int argc, char **argv)
{
	if(argc != 2) {
		fprintf(stderr,"error: number of files: argc = %d != 2\n", argc);
		exit(1);
	}else if(strncmp(argv[1], "-v", 2) == 0 ) {
		fprintf(stderr,"The '%s' creates 2D grids in outer regions.\n", argv[0]);
		fprintf(stderr,"Version 18.08.04 is compiled at %s on %s\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
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
	char urxF[BUFSIZE], outer_file[BUFSIZE];
	long cell_num;
	if(input_file2(fp_i1, urxF, &cell_num, outer_file) != 0) { 
		fprintf(stderr,"input_file2 error!\n");
		exit(1);
	}
	fclose(fp_i1);
	long double *u0 = calloc((2*L[0]+1),sizeof(long double));//revised on 20180725
	long double *u2 = calloc((2*L[2]+1),sizeof(long double));//revised on 20180725
	long double *dudxi = calloc((2*L[2]+1),sizeof(long double));//revised on 20180725
	if(u0 == NULL || u2 == NULL || dudxi == NULL){fprintf(stderr,"u0, u2 and/or dudxi can not be secured!\n");	exit(EXIT_FAILURE);}
	if(input_file3(xi_file, L, u0, u2, dudxi) != 0) { // inputting parameters from u-mesh data 
		fprintf(stderr,"input_file3 error!\n");
		exit(1);
	}
	long double *x = calloc(2*(2*L[0]+1)*2,sizeof(long double));
	long double *h0 = calloc(2*(2*L[0]+1),sizeof(long double));
	if(x == NULL || h0 == NULL){fprintf(stderr,"x and/or h0 can not be secured!\n");	exit(EXIT_FAILURE);}
	if(input_file4(urxF, L, x, h0) != 0) { // inputting parameters from u-mesh data 
		fprintf(stderr,"input_file4 error!\n");
		exit(1);
	}
	FILE *fp_o;
	fp_o = fopen(outer_file,"w");
	if (fp_o == NULL){
		fprintf(stderr,"open error!:  outer_file\n");
		exit(1);
	}
	time_t timer_s = time(0);
	clock_t startClock = clock();
 /* "sys/time.h" is not ANSI but POSIX.
	struct timeval startTime, endTime;// see http://8ttyan.hatenablog.com/entry/2015/02/03/003428
	gettimeofday(&startTime, NULL);
// end: sys/time.h */
	Calc_main(fp_o, L, u0, u2, dudxi, x, h0, cell_num); // calc. of EM and outputting results to an outputfile
	SAFEFREE(u0);	SAFEFREE(u2);	SAFEFREE(dudxi);	SAFEFREE(x);	SAFEFREE(h0);
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
	fprintf(fp_o,"# This file  = %s, Non-uniform data = %s\n", outer_file, xi_file);
	fprintf(fp_o,"# L0 = %ld, L1 = %ld, L2 = %ld, outer cell num = %ld\n", L[0], L[1], L[2], cell_num);
	fprintf(fp_o,"# Info-file  = %s\n", argv[1]);
	fprintf(fp_o,"# Exe-file   = %s : `%s' was compiled at %s on %s by C-version:%ld\n", argv[0], __FILE__, __TIME__, __DATE__, __STDC_VERSION__);
	fprintf(fp_o,"# LDBL_DIG   = %u, LDBL_EPSILON = %.20LE\n", __LDBL_DIG__,  __LDBL_EPSILON__);
	fprintf(fp_o,"# Exec time  = %1.0f sec, CPU time = %1.6f sec\n", t_span, cpusec);
	fprintf(fp_o,"# Created on %s", ctime(&timer_f));
	fprintf(stderr,"# Created on %s", ctime(&timer_f));
	fclose(fp_o);
}
//======================================================================
//  input_file4                                                         
//                                       Last updated on Jun 28, 2018.  
//======================================================================
int input_file4(const char urxF[], const long L[3], long double x[], long double h0[])
{
	FILE *fp_i1;
	fp_i1 = fopen(urxF,"r");
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = 0;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 2*(2*L[0] + 1)) { 
		long j2p1;
		if(sscanf(buf,"%*[^=] %*[=] %ld", &j2p1) == 1){
			if(fgets(buf, sizeof( buf ), fp_i1) == NULL){fprintf(stderr,"Error1! @ input_file4\n"); exit(1);};// # comment
			if( j2p1 == 0  || j2p1 == 2*L[2] ) {
				long j_shift = 0;
				if(j2p1 == 2*L[2] ){j_shift = 2*L[0] + 1;}
				for(long j0 = 0 ; j0 < 2*L[0] +1 ; j0++){
					if(fgets(buf, sizeof( buf ), fp_i1) == NULL){fprintf(stderr,"Error2! @ input_file4\n"); exit(1);};
					long double u2_dummy, u0_dummy, r2_dummy, r0_dummy, x2_dummy, x0_dummy, h2_dummy, h0_dummy;
					long j_dummy;
					if(sscanf(buf,"%ld, %Lf, %Lf, %Lf, %Lf, %Lf, %Lf, %Lf, %Lf", &j_dummy, 
					&u2_dummy, &u0_dummy, &r2_dummy, &r0_dummy, 
					&x2_dummy, &x0_dummy, &h2_dummy, &h0_dummy) == 9){
						if(j_dummy == j0){
							h0[j0 + j_shift] = h0_dummy;
							x[j0*2 + j_shift*2] = x0_dummy;
							x[j0*2 + 1 + j_shift*2] = x2_dummy;
							j_count++;
						}
					}
				}
			}
		}
	}
	fclose(fp_i1);
	return(j_count - 2*(2*L[0] + 1));
}
//======================================================================
//  input_file3                                                         
//                                       Last updated on May 23, 2018.  
//======================================================================
int input_file3(const char xi_file[], const long L[3], long double u0[], long double u2[], long double dudxi[])
{
	FILE *fp_i1;
	fp_i1 = fopen(xi_file,"r");
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = 0;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < (2*L[0] + 1) + (2*L[1] + 1) + (2*L[2] + 1)) { 
		long double xi_dummy, u_dummy, dudxi_dummy;
		if(strncmp(buf, "#", 1) != 0 && sscanf(buf,"%Lf, %Lf, %Lf", &xi_dummy, &u_dummy, &dudxi_dummy) == 3){
			if( j_count >= 0  && j_count < 2*L[0] + 1 ) {
				u0[j_count] = u_dummy;
			}
			if( j_count >= (2*L[0] + 1) + (2*L[1] + 1)  && j_count < (2*L[0] + 1) + (2*L[1] + 1) + (2*L[2] + 1) ) {
				u2[j_count -  (2*L[0] + 1) - (2*L[1] + 1) ] = u_dummy;
				dudxi[j_count - (2*L[0] + 1) - (2*L[1] + 1) ] = dudxi_dummy;
			}
			j_count++;
		}
	}
	fclose(fp_i1);
	return(j_count - (2*L[0] + 1) - (2*L[1] + 1) - (2*L[2] + 1));
}
//======================================================================
//  input_file1                                                         
//                           Last updated on Jun 16, 2018.  
//======================================================================
void rm_space( char *A );
void rm_comma(char *A);
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
	for(long j = 0 ; j < 3 ; j++) {
		if(L[j] < 0){
			fprintf(stderr,"The L[%ld] was not read @ input_file1!\n",j);
			j_count--;
		}
	}
	rm_comma(xi_file); // added on Jul 06, 2018
	return(j_count - 4);
}
//======================================================================
//  input_file2                                                         
//                          Last updated on Jun 28, 2018.  
//======================================================================
int input_file2(FILE *fp_i1, char urxF[], long *cell_num, char outer_file[])
{
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = 0;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 2) { 
		rm_space(buf);
		if(strncmp(buf, "u2r2x", 5) == 0 && sscanf(buf,"%*[^=] %*[=] %s", urxF) == 1) {
			j_count++;
		}else if(strncmp(buf, "outer", 5) == 0 && sscanf(buf,"%*[^=] %*[=] %ld %*[^=] %*[=] %s", cell_num, outer_file) == 2) {
			j_count++;
		}
	}
	rm_comma(urxF); rm_comma(outer_file);// added on Jul 06, 2018
	return(j_count - 2);
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
//  This program removes comma of end in characters,
//                          Last updated on Jul 06, 2018. 
//======================================================================
void rm_comma(char *A)
{
//	A[BUFSIZE - 1] = '\0';// A[] has to include NULL character.
	int i = 0;
	while(A[i] != '\0' && i < BUFSIZE){
		if(A[i] == ',') {
			A[i] = '\0';
			goto replaced;
		}
		i++;
	}
	replaced:;
}
//======================================================================
//  This program is Calc_main Program of mkCurv.                        
//                                       Last updated on Aug 04, 2018.  
//======================================================================
void Calc_main(FILE *fp_o, const long L[3], 
					const long double u0[2*L[0]+1], const long double u2[2*L[2]+1], 
					const long double dudxi[2*L[2]+1], const long double x[2*(2*L[0]+1)*2], 
					const long double h0[2*(2*L[0]+1)], const long cell_num)
{
	for(long sig = 0 ; sig < 2 ; sig++) {
		long i_shift = (2*L[0]+1)*sig;
		long double p0, p2;
		{
			long double vx0 = x[2*L[0]*2 + i_shift*2] - x[0 + i_shift*2];
			long double vx2 = x[2*L[0]*2+1 + i_shift*2] - x[1 + i_shift*2];
			long double cvx = 1./sqrt(vx0*vx0 + vx2*vx2);
			vx0 *= cvx;
			vx2 *= cvx;
			p0 = -vx2;
			p2 = vx0;
		}
		long double du2dxi2 = dudxi[2*L[2]*sig];
		for(long j0 = 1 ; j0 <= 2*cell_num ; j0++){
			if(sig == 0){fprintf(fp_o,"# Bottom region ");}
			else if(sig == 1){fprintf(fp_o,"# Top region ");}
//			fprintf(fp_o,"2l_2 + 1 = %ld\n", sig*(2*L[2]+1) + (2*sig-1)*j0);
			fprintf(fp_o,"2l_2 + 1 = %ld\n", sig*(2*L[2]) + (2*sig-1)*j0);//revised on 20180725
			fprintf(fp_o,"#2l_0 + 1,            k0*u_2,                     k0*u_0,                     k0*r_2,                     k0*r_0,                     k0*x_2,                     k0*x_0,                        h_2,                        h_0\n");
			long double h2 = 1.;
			for(long l0 = 0 ; l0 < (2*L[0]+1) ; l0++){
				long double x2 = x[l0*2+1 + i_shift*2] + 0.5L*du2dxi2*h2*p2*((long double) (2*sig-1)*j0);
				long double x0 = x[l0*2   + i_shift*2] + 0.5L*du2dxi2*h2*p0*((long double) (2*sig-1)*j0);
				long double au2 = u2[sig*(2*L[2])] + 0.5L*du2dxi2*((long double) (2*sig-1)*j0);
				long double ar2 = h2*au2;
				fprintf(fp_o,"%ld, %.20Le, %.20Le, %.20Le, %.20Le, %.20Le, %.20Le, %.20Le, %.20Le\n", 
					l0, 
					au2, u0[l0], 
					ar2, h0[l0 + i_shift]*u0[l0], 
					x2, x0,
					h2, h0[l0 + i_shift]);
			}
			fprintf(fp_o,"\n");//fprintf(fp_o,"\n\n");// comment out on Jul 02, 2018
		}
		if(sig == 0){fprintf(fp_o,"# Start Top-region -----------------------\n\n");}//modified on Aug 04, 2018
	}
}
