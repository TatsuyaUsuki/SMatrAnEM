//  Start of u2r2x.h
#include "header_macro.h"
#include <sys/time.h>
#include "constant.h"
//#define NDEBUG
//  End of u2r2x.h
//======================================================================
//  Main routine of u2r2x.          
// main --- input_file1, input_file2, input_file3, lip_test, Calc_main
//======================================================================
void input_file1(FILE *fp_i1, long L[3], char non_uniformF[], char kappaF[], char urxF[]);
void input_file2(const char non_uniformF[], const long L[3], long double uxi[]);
void input_file3(const char kappaF[], const long L[3], long double kau[]);
void lip_test(const long L[], const long double uxi[], const long double kau[]);
void Calc_main(FILE *fp_o, const long L[3], const long double uxi[], const long double kau[]);
int main(int argc, char **argv)
{
	if(argc != 2) {
		fprintf(stderr,"error: number of files \n");
		exit(EXIT_FAILURE);
	}else if(strncmp(argv[1], "-v", 2) == 0 ) {
		fprintf(stderr,"The '%s' creates 2D coordinate transformation.\n", argv[0]);
		fprintf(stderr,"Version 18.08.10 is compiled at %s on %s\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(stderr," Source code : '%s'\n Author      : Tatsuya Usuki\n URL         : http://www.smatran.org\n", __FILE__);
		fprintf(stderr," Reference   : 'Wave scattering in frequency domain' as 'Formulation.pdf' on May 20, 2018\n");
		exit(0);//normal end
	}
		fprintf(stderr,"Start of input_file1\n");
	long L[3]; //division number: L0 = L[0], L1 = L[1], L2 = L[2]
	char non_uniformF[BUFSIZE], kappaF[BUFSIZE], urxF[BUFSIZE];//file name.
	{
		FILE *fp_i1;
		fp_i1 = fopen(argv[1],"r");
		if (fp_i1 == NULL){
			fprintf(stderr,"open error!: open input-file!\n");
			exit(EXIT_FAILURE);
		}
		input_file1(fp_i1, L, non_uniformF, kappaF, urxF);
		fclose(fp_i1);
	}	fprintf(stderr,"End of input_file1\n");
	long double *uxi = calloc(2*(2*L[0]+1+2*L[1]+1+2*L[2]+1),sizeof(long double));//revised on 20180725
	if(uxi == NULL){fprintf(stderr,"uxi can not be secured!\n");	exit(EXIT_FAILURE);}//revised on 20180725
	input_file2(non_uniformF, L, uxi);	fprintf(stderr,"End of input_file2\n");
	long double *kau = calloc(2*(2*L[2]+1),sizeof(long double));//revised on 20180725
	if(kau == NULL){fprintf(stderr,"kau can not be secured!\n");	exit(EXIT_FAILURE);}//revised on 20180725
	input_file3(kappaF, L, kau);	fprintf(stderr,"End of input_file3\n");

	FILE *fp_o;
	fp_o = fopen(urxF,"w");
	if (fp_o == NULL){
		fprintf(stderr,"open error!:  output-test\n");
		exit(EXIT_FAILURE);
	}
	time_t timer_s = time(0);
	clock_t startClock = clock();
// /* "sys/time.h" is not ANSI but POSIX.
	struct timeval startTime;// see http://8ttyan.hatenablog.com/entry/2015/02/03/003428
	gettimeofday(&startTime, NULL);// end: sys/time.h */
	
	lip_test(L, uxi, kau);//	fprintf(stderr,"Test 'lip' done\n");
	Calc_main(fp_o, L, uxi, kau);// calc. of EM and outputting results to an outputfile
	SAFEFREE(kau);	SAFEFREE(uxi);
	
	time_t timer_f = time(0);
	double t_span = difftime(timer_f,timer_s);
	clock_t endClock = clock();
// /* "sys/time.h" is not ANSI but POSIX.
	struct timeval endTime;// see http://8ttyan.hatenablog.com/entry/2015/02/03/003428
	gettimeofday(&endTime, NULL);
	time_t diffsec = difftime(endTime.tv_sec,startTime.tv_sec);
	suseconds_t diffsub = endTime.tv_usec - startTime.tv_usec;
	t_span = diffsec + diffsub*1.e-6;// end: sys/time.h */
	double cpusec = (endClock - startClock)/(double)CLOCKS_PER_SEC;
	fprintf(fp_o,"# This file  = %s, Non-uniform data = %s\n", urxF, non_uniformF);
	fprintf(fp_o,"# L0 = %ld, L1 = %ld, L2 = %ld\n", L[0], L[1], L[2]);
	fprintf(fp_o,"# Info-file  = %s\n", argv[1]);
	fprintf(fp_o,"# Exe-file   = %s : `%s' was compiled at %s on %s by C-version:%ld\n", argv[0], __FILE__, __TIME__, __DATE__, __STDC_VERSION__);
	fprintf(fp_o,"# Kappa data = %s\n", kappaF);
	fprintf(fp_o,"# Exec time  = %1.6f sec, CPU time= %1.6f sec, Efficiency = %1.2f %%\n", t_span, cpusec, 100.*cpusec/t_span);
	fprintf(fp_o,"# Created on %s", ctime(&startTime.tv_sec));
	fprintf(stderr,"# Created on %s", ctime(&startTime.tv_sec));
	fclose(fp_o);
}
//======================================================================
//  input_file3                     Last updated on Jun 07, 2018.
//======================================================================
void input_file3(const char kappaF[], const long L[3], long double kau[])
{
	FILE *fp_i;
	fp_i = fopen(kappaF,"r");
	if (fp_i == NULL){
		fprintf(stderr,"open '%s' error @ input_file3!\n", kappaF);
		exit(EXIT_FAILURE);
	}
	char buf[BUFSIZE];	// buffer for fgets
	long j_count = 0;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL && j_count < 2*L[2]+1) {
		if(strncmp(buf, "#", 1) != 0) {
			double R_dummy0, R_dummy1, R_dummy2;
			if(sscanf(buf,"%lf, %lf, %lf, %Lf, %Lf", &R_dummy0, &R_dummy1, &R_dummy2, &kau[2*j_count], &kau[2*j_count + 1]) == 5){
				j_count++;
			}
		}
	}
	fclose(fp_i);

	if(j_count != 2*L[2]+1) { 
		fprintf(stderr,"input_file3 error! # of read data = %ld\n", j_count - (2*L[2]+1));
		exit(EXIT_FAILURE);
	}
}
//======================================================================
//  input_file2                     Last updated on Jun 06, 2018.
//======================================================================
void input_file2(const char non_uniformF[], const long L[3], long double uxi[])
{
	FILE *fp_i;
	fp_i = fopen(non_uniformF,"r");
	if (fp_i == NULL){
		fprintf(stderr,"open '%s' error @ input_file2!\n", non_uniformF);
		exit(EXIT_FAILURE);
	}
	char buf[BUFSIZE];	// buffer for fgets
	long j_count = 0;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL && j_count < 2*L[0]+1+2*L[1]+1+2*L[2]+1) {
		if(strncmp(buf, "#", 1) != 0) {
			double R_dummy;
			if(sscanf(buf,"%lf, %Lf, %Lf", &R_dummy, &uxi[2*j_count], &uxi[2*j_count + 1]) == 3){
				j_count++;
			}
		}
	}
	fclose(fp_i);

	if(j_count != 2*L[0]+1+2*L[1]+1+2*L[2]+1) { 
		fprintf(stderr,"input_file2 error! # of read data = %ld\n", j_count - (2*L[0]+1+2*L[1]+1+2*L[2]+1));
		exit(EXIT_FAILURE);
	}
}
void rm_space( char *A );
//======================================================================
//  input_file1 --- rm_space              Last updated on Jun 16, 2018.
//======================================================================
void input_file1(FILE *fp_i1, long L[3], char non_uniformF[], char kappaF[], char urxF[])
{
	char buf[BUFSIZE];	// buffer for fgets
	for(int j = 0 ; j < 3 ; j++) {L[j] = -1;}
	int j_count = 0;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 3) { 
		rm_space(buf);
		if(strncmp(buf, "xi2u", 4) == 0 && sscanf(buf,"%*[^=] %*[=] %s", non_uniformF) == 1){
			j_count++;//j_count == 1
		}else if(strncmp(buf, "kappa", 5) == 0) {
			if(sscanf(buf,"%*[^=] %*[=] %s", kappaF) == 1){
				j_count++;//j_count == 2
			}
		}else if(strncmp(buf, "u2r2x", 5) == 0) {
			if(sscanf(buf,"%*[^=] %*[=] %s", urxF) == 1){
				j_count++;//j_count == 3
			}
		}
	}
	if(j_count != 3) {fprintf(stderr,"File-names can not read @ input_file1!\n");	exit(EXIT_FAILURE);}
	FILE *fp_ix;
	fp_ix = fopen(non_uniformF,"r");
	if (fp_ix == NULL){
		fprintf(stderr,"%s can not open!\n", non_uniformF);
		exit(EXIT_FAILURE);
	}else{
		while(fgets(buf, sizeof( buf ), fp_ix) != NULL && j_count < 3+3){
			if(sscanf(buf,"%*[^=] %*[=] %ld", &L[j_count - 3]) == 1){	j_count++;}//j_count == 3+3
		}
	}
	fclose(fp_ix);
	if(j_count != 3+3) { 
		fprintf(stderr,"input_file1 error! # of read data = %d\n", j_count);
		exit(EXIT_FAILURE);
	}else{
		for(int j = 0 ; j < 3 ; j++) {
			if(L[j] < 0){
				fprintf(stderr,"The L[%d] was not read @ input_file1!\n",j);
				exit(EXIT_FAILURE);
			}
		}		
	}
}
//======================================================================
//  'rm_space' removes spaces of head in characters,
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
	}//	fprintf(stderr,"Read data after removing spaces =%s\n", A);
}
void check_xh(const double x0[], const double x2[], const double h0[], const double h2[], const long L[], const long double uxi[], const long double kau[], const double r0[], const double r2[]);
void gen_h0h2(double h0[], double h2[], const double r0[], const double r2[], const long L[], const long double uxi[], const long double kau[]);
void gen_x2x0(double x2[], double x0[], const double r0[], const double r2[], const long L[], const long double uxi[], const long double kau[]);
void iter_r2r0(const long l0, double r0[], const long l2, double r2[], const long L[], const long double uxi[], const long double kau[]);
void lip(long *j2p1, double *x, const double u2, const long L[], const long double uxi[], const long D);
double integF(const long l0, const double r0[], const long l2, const double r2[], const long L[], const long double uxi[], const long double kau[]);
double zeta(const double r2, const long L[], const long double uxi[], const long double kau[]);
double u(const char s, const long j, const long l, const long L[],const long double uxi[]);
//======================================================================
//  'Calc_main' calculates coordinates transform.                        
//                                  Last updated on Jul 02, 2018.  
//======================================================================
void Calc_main(FILE *fp_o, const long L[3], const long double uxi[], const long double kau[])
{
	double *r0 = calloc((2*L[0]+1)*(2*L[2]+1),sizeof(double));//revised on 20180725
	double *r2 = calloc((2*L[0]+1)*(2*L[2]+1),sizeof(double));//revised on 20180725
	if(r0 == NULL || r2 == NULL){fprintf(stderr,"r0 and/or r2 can not be secured!\n");	exit(EXIT_FAILURE);}//revised on 20180725
	for(long l2 = 0 ; l2 < (2*L[2]+1) ; l2++){
		double u2 = u('n',2,l2,L,uxi);
		double zetau2 = zeta(u2,L,uxi,kau);
		for(long l0 = 0 ; l0 < (2*L[0]+1) ; l0++){
			r0[l0 + (2*L[0]+1)*l2] = u('n',0,l0,L,uxi)*zetau2;
			r2[l0 + (2*L[0]+1)*l2] = u2;
		}
	}//	fprintf(stderr,"End of Initialization for 'r0' and 'r2'\n");
	for(long l2 = 0 ; l2 < (2*L[2]+1) ; l2++){
		long i2p1_0;
		double x_0;
		lip(&i2p1_0, &x_0, 0., L, uxi,0);	assert(0 < i2p1_0 && i2p1_0 < 2*L[0]+1);
		for(long l0 = i2p1_0 ; l0 < (2*L[0]+1) ; l0++){
			iter_r2r0(l0, r0, l2, r2, L,uxi,kau);
		}
		for(long l0 = i2p1_0 - 1 ; l0 >= 0 ; l0--){
			iter_r2r0(l0, r0, l2, r2, L,uxi,kau);
		}
	}//	fprintf(stderr,"Finish iteration for r0 and r2\n");
	double *h0 = calloc((2*L[0]+1)*(2*L[2]+1),sizeof(double));//revised on 20180725
	double *h2 = calloc((2*L[0]+1)*(2*L[2]+1),sizeof(double));//revised on 20180725
	if(h0 == NULL || h2 == NULL){fprintf(stderr,"h0 and/or h2 can not be secured!\n");	exit(EXIT_FAILURE);}//revised on 20180725
	gen_h0h2(h0, h2, r0, r2, L,uxi,kau);
	
	double *x0 = calloc((2*L[0]+1)*(2*L[2]+1),sizeof(double));
	double *x2 = calloc((2*L[0]+1)*(2*L[2]+1),sizeof(double));
	if(x0 == NULL || x2 == NULL){fprintf(stderr,"x0 and/or x2 can not be secured!\n");	exit(EXIT_FAILURE);}//revised on 20180725
	gen_x2x0(x0, x2, r0, r2, L,uxi,kau);
	check_xh(x0, x2, h0, h2, L,uxi,kau, r0, r2);
	for(long l2 = 0 ; l2 < (2*L[2]+1) ; l2++){
		fprintf(fp_o,"#2l_2 + 1 = %ld\n", l2);
		fprintf(fp_o,"#2l_0 + 1,            k0*u_2,                     k0*u_0,                     k0*r_2,                     k0*r_0,                     k0*x_2,                     k0*x_0,                        h_2,                        h_0\n");
		for(long l0 = 0 ; l0 < (2*L[0]+1) ; l0++){
			fprintf(fp_o,"%ld, %.20e, %.20e, %.20e, %.20e, %.20e, %.20e, %.20e, %.20e\n", 
			l0, u('n',2,l2,L,uxi), u('n',0,l0,L,uxi), 
			r2[l0 + (2*L[0]+1)*l2], r0[l0 + (2*L[0]+1)*l2], 
			x2[l0 + (2*L[0]+1)*l2], x0[l0 + (2*L[0]+1)*l2],
			h2[l0 + (2*L[0]+1)*l2], h0[l0 + (2*L[0]+1)*l2]);
		}
		fprintf(fp_o,"\n"); //fprintf(fp_o,"\n\n"); //comment out on Jul 2, 2018
	}
	SAFEFREE(h2);	SAFEFREE(h0);	SAFEFREE(x2);	SAFEFREE(x0);	SAFEFREE(r2);	SAFEFREE(r0);
}
double kappa(const char s, const double u2, const long L[], const long double uxi[], const long double kau[]);
double Diff(const long l0, const long l2, const long id, const double x[], const long L[], const long double uxi[]);
//======================================================================
//  'check_xh' checks x_0, x_2, h0 and h2 by Section A.2, eqs. (A.9) and (A.10).
//                          Last updated on Jun 19, 2018. 
//    check_xh ---- u, kappa, Diff
//======================================================================
void check_xh(const double x0[], const double x2[], const double h0[], const double h2[], const long L[], const long double uxi[], const long double kau[], const double r0[], const double r2[])
{
	fprintf(stderr,"# Check x-vector, scale factors\n");
	fprintf(stderr,"# 2*xi0, 2*xi2, Ortho_xu, Ortho_ru, h0, test-h0, h2, test-h2s\n");
	double max_err = 0.;	long max_l0 = 0;	long max_l2 = 0;
	double may_err = 0.;	long may_l0 = 0;	long may_l2 = 0;
	for(long l2 = 0 ; l2 < (2*L[2]+1) ; l2++){
		for(long l0 = 0 ; l0 < (2*L[0]+1) ; l0++){
			double mr0kb2 = 1. - ( r0[l0 + (2*L[0]+1)*l2] * kappa('b',r2[l0 + (2*L[0]+1)*l2], L,uxi,kau) );
			mr0kb2 = mr0kb2*mr0kb2;
			
			double dx0du2 = Diff(l0, l2, 2, x0, L, uxi);
			double dx2du2 = Diff(l0, l2, 2, x2, L, uxi);
			double dr0du2 = Diff(l0, l2, 2, r0, L, uxi);
			double dr2du2 = Diff(l0, l2, 2, r2, L, uxi);
			
			double dx0du0 = Diff(l0, l2, 0, x0, L, uxi);
			double dx2du0 = Diff(l0, l2, 0, x2, L, uxi);
			double dr0du0 = Diff(l0, l2, 0, r0, L, uxi);
			double dr2du0 = Diff(l0, l2, 0, r2, L, uxi);
			
			double Ortho_xu = dx0du0*dx0du2 + dx2du0*dx2du2;
			double hh0 = sqrt(dx0du0*dx0du0 + dx2du0*dx2du0);
			double hh2 = sqrt(dx0du2*dx0du2 + dx2du2*dx2du2);
			Ortho_xu /= (hh0*hh2);
			if(fabs(Ortho_xu) > fabs(max_err)){
				max_err = Ortho_xu;
				max_l0 = l0;
				max_l2 = l2;
			}
			double Ortho_ru = (dr0du0*dr0du2 + mr0kb2*dr2du0*dr2du2)/sqrt(dr0du0*dr0du0 + mr0kb2*dr2du0*dr2du0)/sqrt(dr0du2*dr0du2 + mr0kb2*dr2du2*dr2du2);
			if(fabs(Ortho_ru) > fabs(may_err)){
				may_err = Ortho_ru;
				may_l0 = l0;
				may_l2 = l2;
			}
			fprintf(stderr,"%ld, %ld, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e\n",
				l0,l2, Ortho_xu, Ortho_ru, 
				h0[l0 + (2*L[0]+1)*l2], hh0, 
				h2[l0 + (2*L[0]+1)*l2], hh2);
		}
	}
	fprintf(stderr,"#max_err for dx/du = %.4e at (x0= %.4e, x2= %.4e) (2*(l0,l2)=(%ld, %ld))\n",//revised on 20180810
		max_err, x0[max_l0 + (2*L[0]+1)*max_l2], x2[max_l0 + (2*L[0]+1)*max_l2], max_l0, max_l2);
	fprintf(stderr,"#max_err for dr/du = %.4e at (x0= %.4e, x2= %.4e) (2*(l0,l2)=(%ld, %ld))\n",//revised on 20180810
		may_err, x0[may_l0 + (2*L[0]+1)*may_l2], x2[may_l0 + (2*L[0]+1)*may_l2], may_l0, may_l2);
}
//======================================================================
//  'Diff' generates differential operation 'check_xh'
//                          Last updated on Aug 10, 2018. 
//======================================================================
double Diff(const long l0, const long l2, const long id, const double x[], const long L[], const long double uxi[])
{
	long s0, s2;
	if(id == 0){
		s0 = 1;
		s2 = 0;
	}else if(id == 2){
		s0 = 0;
		s2 = 1;
	}else{
		fprintf(stderr,"Error id = %ld @ Diff\n", id);
		exit(EXIT_FAILURE);
	}
	long j = s0*l0 + s2*l2;
	double dxdu;//double du_f, du_b, dxdu;
	if(j == 0){
		double du_f = u('n',id,j+1,L,uxi) - u('n',id,j,L,uxi);
/*		double du_b = u('n',id,j+2,L,uxi) - u('n',id,j,L,uxi);
		dxdu = (du_b*(x[l0+s0 + (2*L[0]+1)*(l2+s2)] - x[l0 + (2*L[0]+1)*(l2)])/du_f
				-du_f*(x[l0+2*s0 + (2*L[0]+1)*(l2+2*s2)] - x[l0 + (2*L[0]+1)*(l2)])/du_b)/(du_b - du_f);*/
		dxdu = (x[l0+s0 + (2*L[0]+1)*(l2+s2)] - x[l0 + (2*L[0]+1)*(l2)])/du_f;//This is better than the above in testing, revised on 20180810
	}else if(j > 0 && j < 2*L[id]){
		double du_f = u('n',id,j+1,L,uxi) - u('n',id,j,L,uxi);
		double du_b = u('n',id,j,L,uxi) - u('n',id,j-1,L,uxi);
		dxdu = (du_b/(du_f+du_b))*(x[l0+s0 + (2*L[0]+1)*(l2+s2)] - x[l0 + (2*L[0]+1)*(l2)])/du_f;
		dxdu += (du_f/(du_f+du_b))*(x[l0 + (2*L[0]+1)*(l2)] - x[l0-s0 + (2*L[0]+1)*(l2-s2)])/du_b;
	}else if(j == 2*L[id]){
		double du_f = u('n',id,j,L,uxi) - u('n',id,j-1,L,uxi);
/*		double du_b = u('n',id,j,L,uxi) - u('n',id,j-2,L,uxi);
		dxdu = (du_b*(x[l0 + (2*L[0]+1)*(l2)] - x[l0-s0 + (2*L[0]+1)*(l2-s2)])/du_f
				-du_f*(x[l0 + (2*L[0]+1)*(l2)] - x[l0-2*s0 + (2*L[0]+1)*(l2-2*s2)])/du_b)/(du_b - du_f);*/
		dxdu = (x[l0 + (2*L[0]+1)*(l2)] - x[l0-s0 + (2*L[0]+1)*(l2-s2)])/du_f;//This is better than the above in testing, revised on 20180810
	}else{	fprintf(stderr,"Error in Diff!\n");	exit(EXIT_FAILURE);}//	assert(du_f > 0.);	assert(du_b > 0.);
	return(dxdu);
}
//======================================================================
//  'gen_h0h2' generates h_0 and h_2 by eq. (2.3).
//                          Last updated on Jun 18, 2018. 
//    gen_h0h2 ---- zeta, u, kappa
//======================================================================
void gen_h0h2(double h0[], double h2[], const double r0[], const double r2[], const long L[], const long double uxi[], const long double kau[])
{
	assert(2*L[2] >= 2);
	for(long l2 = 0 ; l2 < (2*L[2]+1) ; l2++){
		double du2_f;
		double du2_b = 1.;
		if(l2 == 0){
			du2_f = u('n',2,l2+1,L,uxi) - u('n',2,l2,L,uxi);
			du2_b = u('n',2,l2+2,L,uxi) - u('n',2,l2,L,uxi);
		}else if(l2 > 0 && l2 < 2*L[2]){
			du2_f = u('n',2,l2+1,L,uxi) - u('n',2,l2,L,uxi);
			du2_b = u('n',2,l2,L,uxi) - u('n',2,l2-1,L,uxi);
		}else{
			du2_f = u('n',2,l2,L,uxi) - u('n',2,l2-1,L,uxi);
			du2_b = u('n',2,l2,L,uxi) - u('n',2,l2-2,L,uxi);
		}	assert(du2_f > 0.);	assert(du2_b > 0.);
		for(long l0 = 0 ; l0 < (2*L[0]+1) ; l0++){
			double rr2 = r2[l0 + (2*L[0]+1)*l2];
			{
				double rr0 = r0[l0 + (2*L[0]+1)*l2];
				double mr0kb = 1. - rr0*kappa('b',rr2, L,uxi,kau);	assert(mr0kb >= 0.);
				double r0kw = rr0*kappa('w',rr2, L,uxi,kau); 
				double zetar2 = zeta(rr2,L,uxi,kau);
				h2[l0 + (2*L[0]+1)*l2] = sqrt(mr0kb*mr0kb + r0kw*r0kw);
				h0[l0 + (2*L[0]+1)*l2] = mr0kb*zetar2/sqrt(mr0kb*mr0kb + r0kw*r0kw);
			}
			double dr2du2;
			if(l2 == 0){
				double R_f = (r2[l0 + (2*L[0]+1)*(l2+1)] - rr2)/du2_f;
				double R_b = (r2[l0 + (2*L[0]+1)*(l2+2)] - rr2)/du2_b;
				dr2du2 = (du2_b*R_f - du2_f*R_b)/(du2_b - du2_f);
			}else if(l2 > 0 && l2 < 2*L[2]){
				dr2du2 = (du2_b/(du2_f+du2_b))*(r2[l0 + (2*L[0]+1)*(l2+1)] - rr2)/du2_f;
				dr2du2 += (du2_f/(du2_f+du2_b))*(rr2 - r2[l0 + (2*L[0]+1)*(l2-1)])/du2_b;
			}else{
				double R_f = (rr2 - r2[l0 + (2*L[0]+1)*(l2-1)])/du2_f;
				double R_b = (rr2 - r2[l0 + (2*L[0]+1)*(l2-2)])/du2_b;
				dr2du2 = (du2_b*R_f - du2_f*R_b)/(du2_b - du2_f);
			}
			h2[l0 + (2*L[0]+1)*l2] *= (dr2du2);
		}
	}	
}
void integ_r2(double F[], const double f[], const long l0, const double r2[], const long L[], const long double uxi[]);
//======================================================================
//  'gen_x2x0' generates x_2 and x_0 by eq. (2.1).
//                          Last updated on Jun 25, 2018. 
//    gen_x2x0 ---- integ_r2
//======================================================================
void gen_x2x0(double x0[], double x2[], const double r0[], const double r2[], const long L[], const long double uxi[], const long double kau[])
{
	double *kappa_b = calloc((2*L[2]+1),sizeof(double));
	double *fc = calloc((2*L[2]+1),sizeof(double));
	double *gs = calloc((2*L[2]+1),sizeof(double));
	double *dummyM = calloc((2*L[2]+1),sizeof(double));
	if(kappa_b == NULL || fc == NULL || gs == NULL || dummyM == NULL){fprintf(stderr,"kappa_b, fc, gs and/or dummyM can not be secured!\n");	exit(EXIT_FAILURE);}
	for(long l0 = 0 ; l0 < (2*L[0]+1) ; l0++){
		for(long l2 = 0 ; l2 < (2*L[2]+1) ; l2++){ kappa_b[l2] = kappa('b',r2[l0 + (2*L[0]+1)*l2], L,uxi,kau);}
		integ_r2(dummyM, kappa_b, l0,r2,L,uxi);
		for(long l2 = 0 ; l2 < (2*L[2]+1) ; l2++){
			assert(-Pi <= dummyM[l2] && dummyM[l2] <= Pi);
			fc[l2] = cos(dummyM[l2]);
			gs[l2] = sin(dummyM[l2]);
		}
		for(long l2 = 0 ; l2 < (2*L[2]+1) ; l2++){
			x2[l0 + (2*L[0]+1)*l2] = -r0[l0 + (2*L[0]+1)*l2]*gs[l2];
			x0[l0 + (2*L[0]+1)*l2] = r0[l0 + (2*L[0]+1)*l2]*fc[l2];
		}
		integ_r2(dummyM, fc, l0,r2,L,uxi);
		for(long l2 = 0 ; l2 < (2*L[2]+1) ; l2++){	x2[l0 + (2*L[0]+1)*l2] += dummyM[l2];}
		integ_r2(dummyM, gs, l0,r2,L,uxi);
		for(long l2 = 0 ; l2 < (2*L[2]+1) ; l2++){	x0[l0 + (2*L[0]+1)*l2] += dummyM[l2];}
	}
	SAFEFREE(dummyM);	SAFEFREE(gs);	SAFEFREE(fc);	SAFEFREE(kappa_b);
}
//======================================================================
//  'integ_r2' integrates f from 0 to r2 (l2).
//                          Last updated on Jun 19, 2018. 
//    integ_r2 ---- lip
//======================================================================
void integ_r2(double F[], const double f[], const long l0, const double r2[], const long L[], const long double uxi[])
{
	assert(0 <= l0 && l0 < 2*L[0]+1);
	long l2_b0 = -1;
	double x_b0 = -1;
	for(long l2 = 0 ; l2 < 2*L[2] ; l2++){
		double r_b = r2[l0 + (2*L[0]+1)*l2];
		double r_f = r2[l0 + (2*L[0]+1)*(l2+1)];
		if(r_b <= 0. && r_f > 0.){
			l2_b0 = l2;
			x_b0 = r_f/(r_f - r_b);
			goto fix_0;
		}
	}
	fprintf(stderr,"# 'integ_r2' can not define l2_b0!");
	exit(EXIT_FAILURE);
	fix_0:assert(l2_b0 >= 0 && x_b0 > 0. && x_b0 <= 1.);
	for(long l2 = 0 ; l2 < 2*L[2]+1 ; l2++){
		long min_i2p1, max_i2p1;
		double sig, sum;
		if(l2 < l2_b0){
			min_i2p1 = l2;
			max_i2p1 = l2_b0;
			sig = -1.;
			double dr = r2[l0 + (2*L[0]+1)*(l2_b0 + 1)] - r2[l0 + (2*L[0]+1)*(l2_b0)];
			sum = (1. - x_b0)*dr*0.5*((1. + x_b0)*f[l2_b0]+(1. - x_b0)*f[l2_b0+1]);
		}else if(l2_b0 < l2){
			min_i2p1 = l2_b0 + 1;
			max_i2p1 = l2;
			sig = 1.;
			double dr = r2[l0 + (2*L[0]+1)*(l2_b0 + 1)] - r2[l0 + (2*L[0]+1)*(l2_b0)];
			sum = x_b0*dr*0.5*(x_b0*f[l2_b0]+(2. - x_b0)*f[l2_b0+1]);
		}else{
			min_i2p1 = l2_b0;
			max_i2p1 = min_i2p1;
			sig = 0.;
			sum = 0.;
		}
		for(long j2 = min_i2p1 ; j2 < max_i2p1 ; j2++){
			double dr = r2[l0 + (2*L[0]+1)*(j2 + 1)] - r2[l0 + (2*L[0]+1)*(j2)];
			sum += dr*0.5*(f[j2] + f[j2+1]);
		}
		F[l2] = sum*sig;
	}
}
//======================================================================
//  'iter_r2r0' iterates eq. (2.2).
//                          Last updated on Jun 23, 2018. 
//     iter_r2r0 ---- u, integF, zeta
//======================================================================
void iter_r2r0(const long l0, double r0[], const long l2, double r2[], const long L[], const long double uxi[], const long double kau[])
{
	double pre_r2 = r2[l0 + (2*L[0]+1)*l2];
	double pre_r0 = r0[l0 + (2*L[0]+1)*l2];
	double damp = 0.5;
	for(long id = 0 ; id < 5 ; id++){
		r2[l0 + (2*L[0]+1)*l2] = pre_r2;
		r0[l0 + (2*L[0]+1)*l2] = pre_r0;
		for(long i0 = 0 ; i0 < 10000 ; i0++){
			double ini_r2 = r2[l0 + (2*L[0]+1)*l2];
			double ini_r0 = r0[l0 + (2*L[0]+1)*l2];
			
			r2[l0 + (2*L[0]+1)*l2] = u('n',2,l2,L,uxi) - integF(l0, r0, l2, r2, L,uxi,kau);
			r2[l0 + (2*L[0]+1)*l2] *= damp;
			r2[l0 + (2*L[0]+1)*l2] += (1.-damp)*ini_r2;
			
			if(r2[l0 + (2*L[0]+1)*l2] > u('n',2,2*L[2],L,uxi)){r2[l0 + (2*L[0]+1)*l2] = u('n',2,2*L[2],L,uxi);}
			if(r2[l0 + (2*L[0]+1)*l2] < u('n',2,0,L,uxi)){r2[l0 + (2*L[0]+1)*l2] = u('n',2,0,L,uxi);}
			
			r0[l0 + (2*L[0]+1)*l2] = u('n',0,l0,L,uxi)*zeta(r2[l0 + (2*L[0]+1)*l2],L,uxi,kau);
			
			if(ini_r2 == r2[l0 + (2*L[0]+1)*l2] && ini_r0 == r0[l0 + (2*L[0]+1)*l2]) {
				goto convergence;
			}
		}
		damp *= 0.5;
	}
	fprintf(stderr,"# eqs. was not converged at 2*l0=%ld, 2*l2=%ld!", l0,l2);
	exit(EXIT_FAILURE);
	convergence:;
}
double F2D(const double rr0, const double rr2, const long L[], const long double uxi[], const long double kau[]);
//======================================================================
//  'integF' integrates F_2D from 0 to u0 (l0).
//                          Last updated on Jun 13, 2018. 
//     integF ---- lip, F2D
//======================================================================
double integF(const long l0, const double r0[], const long l2, const double r2[], const long L[], const long double uxi[], const long double kau[])
{
	assert(0 <= l0 && l0 < 2*L[0]+1);
	assert(0 <= l2 && l2 < 2*L[2]+1);
	long min_i2p1 = -1;
	long max_i2p1 = -1;
	double min_x = -1.;
	double max_x = -1.;
	double sig = 0.;
	{
		long i2p1_0;
		double x_0;
		lip(&i2p1_0, &x_0, 0., L, uxi,0);
		assert(0. <= x_0 && x_0 <= 1.);
		assert(0 <= i2p1_0 && i2p1_0 < 2*L[0]+1);
		if(uxi[2*l0] >= 0.L){
			min_i2p1 = i2p1_0;
			min_x = x_0;
			max_i2p1 = l0;
			assert(0 <= max_i2p1 && max_i2p1 < 2*L[0]+1);
			max_x = 1.;
			sig = 1.;
		}else{
			min_i2p1 = l0;
			min_x = 1.;
			max_i2p1 = i2p1_0;
			max_x = x_0;
			sig = -1.;
		}
	}

	double sum = 0.;
	if(max_i2p1 > min_i2p1){
		sum = (min_x)*0.5*(
				min_x*F2D(r0[min_i2p1 + (2*L[0]+1)*l2], r2[min_i2p1 + (2*L[0]+1)*l2], L,uxi,kau)*uxi[2*min_i2p1 + 1]
				+ ((1. - min_x) + 1.)*F2D(r0[min_i2p1+1 + (2*L[0]+1)*l2], r2[min_i2p1+1 + (2*L[0]+1)*l2], L,uxi,kau)*uxi[2*(min_i2p1+1) + 1]
			);
		for(long i2p1 = min_i2p1+1 ; i2p1 < max_i2p1 ; i2p1++){
			sum += 0.5*(
				F2D(r0[i2p1 + (2*L[0]+1)*l2], r2[i2p1 + (2*L[0]+1)*l2], L,uxi,kau)*uxi[2*i2p1 + 1]
				+ F2D(r0[i2p1+1 + (2*L[0]+1)*l2], r2[i2p1+1 + (2*L[0]+1)*l2], L,uxi,kau)*uxi[2*(i2p1+1) + 1]
			);
		}
		if(max_i2p1 < 2*L[0]){
			sum += (1.-max_x)*0.5*(
				(1. + max_x)*F2D(r0[max_i2p1 + (2*L[0]+1)*l2], r2[max_i2p1 + (2*L[0]+1)*l2], L,uxi,kau)*uxi[2*max_i2p1 + 1]
				+ (1. - max_x)*F2D(r0[max_i2p1+1 + (2*L[0]+1)*l2], r2[max_i2p1+1 + (2*L[0]+1)*l2], L,uxi,kau)*uxi[2*(max_i2p1+1) + 1]
				);
		}
	}else{
		if(max_i2p1 == min_i2p1 && min_x >= max_x){
			sum = (min_x-max_x)*0.5*((min_x + max_x)*F2D(r0[min_i2p1 + (2*L[0]+1)*l2], r2[min_i2p1 + (2*L[0]+1)*l2], L,uxi,kau)*uxi[2*min_i2p1 + 1]
				+ ((1. - min_x) + (1.- max_x))*F2D(r0[min_i2p1+1 + (2*L[0]+1)*l2], r2[min_i2p1+1 + (2*L[0]+1)*l2], L,uxi,kau)*uxi[2*(min_i2p1+1) + 1]);
		}else{
			fprintf(stderr,"# integF error!: max_i2p1 =%ld, min_i2p1=%ld, min_x=%.2e, max_x=%.2e", max_i2p1, min_i2p1, min_x, max_x);
			exit(EXIT_FAILURE);
		}
	}
	return(sum*sig*0.5);//Note that xi = 0.5*i2p1 on Jun 14, 2018
}
//======================================================================
//  'F2D' calculates F_2D in eq. (2.2).
//                          Last updated on Jun 13, 2018. 
//      F2D ---- kappa
//======================================================================
double F2D(const double rr0, const double rr2, const long L[], const long double uxi[], const long double kau[])
{
	double mr0kb = 1. - rr0*kappa('b',rr2, L,uxi,kau);
	double r0kw = rr0*kappa('w',rr2, L,uxi,kau); 
	return( (r0kw*zeta(rr2, L,uxi,kau))/(mr0kb*mr0kb + r0kw*r0kw) );	
}
//======================================================================
//  'zeta' calculates zeta(r2) in eq.(2.2).
//                          Last updated on Jun 12, 2018. 
//      zeta ---- lip
//======================================================================
double zeta(const double r2, const long L[], const long double uxi[], const long double kau[])
{
	long min_i2p1 = -1;
	long max_i2p1 = -1;
	double min_x = -1.;
	double max_x = -1.;
	double sig = 0.;
	if(r2 >= 0.){
		lip(&max_i2p1, &max_x, r2, L, uxi,2);
		lip(&min_i2p1, &min_x, 0., L, uxi,2);
		sig = 1.;
	}else{
		lip(&max_i2p1, &max_x, 0., L, uxi,2);
		lip(&min_i2p1, &min_x, r2, L, uxi,2);
		sig = -1.;
	}
	long ishift = 2*(2*L[0] + 1 + 2*L[1] + 1);
	double sum = 0.;
	if(max_i2p1 > min_i2p1){//Note that kau[2*i2p1+1] means kappa_w!
		sum = (min_x)*0.5*(min_x*kau[2*min_i2p1+1]*uxi[2*min_i2p1 + 1 + ishift]
				+ ((1. - min_x) + 1.)*kau[2*(min_i2p1+1)+1]*uxi[2*(min_i2p1+1) + 1 + ishift]);
		for(long i2p1 = min_i2p1+1 ; i2p1 < max_i2p1 ; i2p1++){
			sum += 0.5*(kau[2*i2p1+1]*uxi[2*i2p1 + 1 + ishift]
					+ kau[2*(i2p1+1)+1]*uxi[2*(i2p1+1) + 1 + ishift]);
		}
		sum += (1.-max_x)*0.5*((1. + max_x)*kau[2*max_i2p1+1]*uxi[2*max_i2p1 + 1 + ishift]
				+ (1. - max_x)*kau[2*(max_i2p1+1)+1]*uxi[2*(max_i2p1+1) + 1 + ishift]);
	}else{
		if(max_i2p1 == min_i2p1 && min_x >= max_x){
			sum = (min_x-max_x)*0.5*((min_x + max_x)*kau[2*min_i2p1+1]*uxi[2*min_i2p1 + 1 + ishift]
				+ ((1. - min_x) + (1.- max_x))*kau[2*(min_i2p1+1)+1]*uxi[2*(min_i2p1+1) + 1 + ishift]);
		}else{
			fprintf(stderr,"# zeta error!");
			exit(EXIT_FAILURE);
		}
	}
	return(exp(sum*sig*0.5));//Note that xi = 0.5*i2p1 on Jun 14, 2018
}
//======================================================================
//  'lip_test' is test program for 'lip'.
//                          Last updated on Jun 13, 2018. 
//      lip_test ---- lip
//======================================================================
void lip_test(const long L[], const long double uxi[], const long double kau[])
{
	fprintf(stderr,"# Mesh number: L0 = %ld, L1 = %ld, L2 = %ld\n", L[0], L[1], L[2]);
	for(long i0 = 0 ; i0 < 3 ; i0++){
		fprintf(stderr,"# xi_%ld, u_%ld, du_%ld/dxi_%ld\n", i0, i0, i0, i0);
		for(long i1 = 0 ; i1 < 2*L[i0]+1 ; i1++){
			fprintf(stderr,"%.1f, %.20e, %.20e\n", 0.5*((double) i1), u('n',i0,i1,L,uxi), u('d',i0,i1,L,uxi));
		}
	}
	fprintf(stderr,"#\n#Test 'lip'\n");

	fprintf(stderr,"#\n# xi_0, u_0\n");
	for(long i0 = 0 ; i0 < 2*L[0]+1 ; i0++){
		fprintf(stderr,"%.1f, %.16e\n",
			0.5*((double) i0), u('n',0,i0,L,uxi));
	}
	double x = -1.;
	long j2p1 = -1;
	double cu;
	double cfa[6]={0.00028, 0.345678, 0.4999, 0.5001, 0.5657, 0.99977};
	fprintf(stderr,"# xi_0, cu0, lip-u0\n"); 
	for(long i0 = 0 ; i0 < 6 ; i0++){
		cu = (double) (uxi[0] + cfa[i0]*(uxi[2*(2*L[0])] - uxi[0]));
		lip(&j2p1, &x, cu, L, uxi,0);
		fprintf(stderr,"%.1f, %.16e, %.16e\n", 
			0.5*((double) j2p1), cu, (double) (x*uxi[2*j2p1] + (1.-x)*uxi[2*(j2p1+1)]));
	}
	fprintf(stderr,"# xi_2, cu2, lip-u2, kappa_b, kappa_w, zeta\n"); 
	long ishift = 2*(2*L[0] + 1 + 2*L[1] + 1);
	for(long i0 = 0 ; i0 < 6 ; i0++){
		cu = (double) (uxi[ishift] + cfa[i0]*(uxi[2*(2*L[2]) + ishift] - uxi[ishift]));
		lip(&j2p1, &x, cu, L, uxi,2);
		fprintf(stderr,"%.1f, %.16e, %.16e, %.16e, %.16e, %.16e\n", 
			0.5*((double) j2p1), cu, (double) (x*uxi[2*j2p1+ishift] + (1.-x)*uxi[2*(j2p1+1) + ishift]),
			kappa('b',cu,L, uxi,kau), kappa('w',cu,L, uxi,kau), zeta(cu,L,uxi,kau));
	}
	fprintf(stderr,"#\n# xi_2, u_2, kappa_b, kappa_w, zeta\n");
	for(long i0 = 0 ; i0 < 2*L[2]+1 ; i0++){
		fprintf(stderr,"%.1f, %.16e, %.16e, %.16e, %.16e\n",
			0.5*((double) i0), u('n',2,i0,L,uxi),
			kappa('b',u('n',2,i0,L,uxi), L,uxi,kau),
			kappa('w',u('n',2,i0,L,uxi), L,uxi,kau), 
			zeta(u('n',2,i0,L,uxi), L,uxi,kau));
	}
}
//======================================================================
//  'lip' calculates linear interpolate.
//                          Last updated on Jun 13, 2018. 
//======================================================================
void lip(long *j2p1, double *x, const double u2, const long L[], const long double uxi[], const long D)
{
	assert(D == 0 || D == 2);
	long ishift = 0;
	for(long i0 = 0 ; i0 < D ; i0++){ishift += 2*(2*L[i0] + 1);}
///*
	if(((double) uxi[ishift]) > u2 || u2 > ((double) uxi[2*(2*L[D]) + ishift])){
		fprintf(stderr,"#Error! u%ld = %.2e, uxi[ishift] = %.2Le, uxi[2*(2*L[D]) + ishift] = %.2Le", D, u2, uxi[ishift], uxi[2*(2*L[D]) + ishift]);
	}//*/
	assert(((double) uxi[ishift]) <= u2 && u2 <= ((double) uxi[2*(2*L[D]) + ishift]));
	
	for(long i0 = 0 ; i0 < 2*L[D] ; i0++){
		if(((double) uxi[2*i0 + ishift]) == u2 ){
			*j2p1 = i0;
			*x = 1.;
			goto define_x;
		}else if(((double) uxi[2*i0 + ishift]) < u2 && u2 < ((double) uxi[2*(i0 + 1) + ishift])) {
			*j2p1 = i0;
			*x = (uxi[2*(i0 + 1) + ishift] - u2)/( uxi[2*(i0 + 1) + ishift] - uxi[2*i0 + ishift]);
			goto define_x;
		}else if(u2 == ((double) uxi[2*(2*L[D]) + ishift])){
			*j2p1 = 2*L[D]-1;
			*x = 0.;
			goto define_x;
		}
	}
	fprintf(stderr,"# 'lip' can not generate x!");
	fprintf(stderr,"# u2 = %.2e, min uxi = %.2Le, max uxi = %.2Le", u2, uxi[ishift], uxi[2*(2*L[2]) + ishift]);
	exit(EXIT_FAILURE);
	define_x:;
}
//======================================================================
//  'kappa' generates kappa_b and kappa_w in Chapter 2.
//                          Last updated on Jun 12, 2018. 
//      kappa ---- lip
//======================================================================
double kappa(const char s, const double u2, const long L[], const long double uxi[], const long double kau[])
{
	long j2p1;
	double x;
	lip(&j2p1, &x, u2, L, uxi,2);
	long ishift = 0;
	if(s == 'b') {
		ishift = 2*j2p1 ;
	}else if(s == 'w'){
		ishift = 2*j2p1 + 1;
	}else{
		fprintf(stderr,"# '%c' error @ kappa!", s);
		exit(EXIT_FAILURE);
	}
	return(x*kau[ishift] + (1.-x)*kau[ishift+2]);
}
//======================================================================
//  'u' generates u_j by 'n' and du_j/dxi by 'd' in Section 6.2.
//                          Last updated on Aug 10, 2018. 
//======================================================================
double u(const char s, const long j, const long l, const long L[],const long double uxi[])
{
	assert(j >= 0 && j < 3);	assert(l >= 0); assert(l < 2*L[j]+1);
	long ishift = 0;
	for(long i0 = 0 ; i0 < j ; i0++){ishift += 2*(2*L[i0] + 1);}
	if(s == 'n') {
		ishift += 2*l ;
	}else if(s == 'd'){
		ishift += 2*l + 1;
	}else{
		fprintf(stderr,"# '%c' error @ u!", s);
		exit(EXIT_FAILURE);
	}
	assert(ishift < 2*(2*L[0]+1+2*L[1]+1+2*L[2]+1));
	return(uxi[ishift]);
}
