//  Start of Mcalc.h
#include "header_macro.h"
#include "constant.h"
void main_calc(FILE *fp_i2, const int i1, const int L[], const double u01[], const double du01[], const char f_prefix[], const double inv_k0);
void input_file1(FILE *fp_i1, const int L[2], const double inv_k0, double u01[], double du01[]);
void input_file0(FILE *fp_i0, char f_prefix[], char nonuniform[], int L[2], double *inv_k0);
//  End of Mcalc.h
//======================================================================
//   main ---- input_file0, main_calc             
//======================================================================
int main(int argc, char **argv)
{
	{	time_t timer_ini = time(0);	fprintf(stderr,"# Start time of mode calculation = %s\n", ctime(&timer_ini));}
	if(argc != 2) {	fprintf(stderr,"error: number of files \n");	exit(EXIT_FAILURE);}
	else if(strncmp(argv[1], "-v", 2) == 0 || strcmp(argv[1], "--version") == 0 ) {
		fprintf(stderr,"The '%s' creates Poynting vector along z-axis.\n", argv[0]);
		fprintf(stderr,"Version 20.07.25 is compiled at %s on %s.\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(stderr," Source code : '%s'\n Author      : Tatsuya Usuki\n URL         : http://www.smatran.org\n", __FILE__);
		fprintf(stderr," References  : 'Wave scattering in frequency domain' as 'Formulation.pdf' on Aug 14, 2019.\n");
		fprintf(stderr,"There is NO warranty.\n");
		exit(EXIT_SUCCESS);//normal end
	}
//------- begin reading fundamental parameters -------
	FILE *fp_i0;
	fp_i0 = fopen(argv[1],"r");
	if (fp_i0 == NULL){	fprintf(stderr,"open error!: open input-file0!\n");	exit(EXIT_FAILURE);}
	fprintf(stderr,"The input file0: %s\n",argv[1]);
	
	char f_prefix[BUFSIZE], nonuniform[BUFSIZE];
	int L[2];
	double inv_k0;
	input_file0(fp_i0, f_prefix, nonuniform, L, &inv_k0);
	fprintf(stderr,"L[0] = %d, L[1] = %d, inv_k0 = %.5e\n", L[0], L[1], inv_k0);
	if(fclose(fp_i0) != 0) {	fprintf(stderr,"fclose error after input_file0!\n");	exit(EXIT_FAILURE);}
//-------    begin reading nonuniform u0, u1    -------
	FILE *fp_i1;
	fp_i1 = fopen(nonuniform,"r");
	if (fp_i1 == NULL){	fprintf(stderr,"open error!: open input-file1!\n");	exit(EXIT_FAILURE);}
	fprintf(stderr,"The input file1: %s\n",nonuniform);
	double *u01 = calloc(2*(L[0]+L[1]),sizeof(double));
	double *du01 = calloc(2*(L[0]+L[1]),sizeof(double));
	input_file1(fp_i1, L, inv_k0, u01, du01);
	if(fclose(fp_i1) != 0) {	fprintf(stderr,"fclose error after input_file1!\n");	exit(EXIT_FAILURE);}
//-------    begin reading Xi    -------
	for(int i1 = 0 ; i1 < 2 ; i1++){
		char Xi_name[BUFSIZE];
		Xi_name[0] = '\0';
		if(i1 == 0){
			snprintf(Xi_name,sizeof(Xi_name),"%s_bXi.dat", f_prefix);
		}else{
			snprintf(Xi_name,sizeof(Xi_name),"%s_tXi.dat", f_prefix);
		}
		FILE *fp_i2;
		fp_i2 = fopen(Xi_name,"r");
		if(fp_i2 != NULL) {
			main_calc(fp_i2, i1, L, u01, du01, f_prefix, inv_k0);//revised on 20200714
			fprintf(stderr,"The Xi file: %s\n",Xi_name);
			if(fclose(fp_i2) != 0) {	fprintf(stderr,"fclose error after input_file2!\n");	exit(EXIT_FAILURE);}
		}
	}
	SAFEFREE(u01);	SAFEFREE(du01);
}
//======================================================================
//  main_calc   
//                           Last updated on Jul 26, 2020  
//======================================================================
void main_calc(FILE *fp_i2, const int i1, const int L[], const double u01[], const double du01[], const char f_prefix[], const double inv_k0)
{
	char buf[BUFSIZE];	// buffer for fgets
	int mode_n;
	
	double complex *Xi = calloc(4*L[0]*L[1], sizeof(double complex));
	double *Poynting = calloc(2*L[0]*L[1], sizeof(double));
	
	while(fgets(buf, sizeof( buf ), fp_i2) != NULL){
		if(strncmp(buf, "ReXi", 4) == 0 ){
			if(sscanf(buf,"%*[^[] %*[[] %d ", &mode_n)!=1){	fprintf(stderr,"mode_n error!\n");	exit(EXIT_FAILURE);}
			if(mode_n >= 4){goto finished;}
			char out_name[BUFSIZE];
			out_name[0] = '\0';
			if(i1 == 0){	snprintf(out_name,sizeof(out_name),"%s_b%d_Poynting.dat", f_prefix, mode_n);}
			else if(i1 == 1){	snprintf(out_name,sizeof(out_name),"%s_t%d_Poyinting.dat", f_prefix, mode_n);}
			else{	fprintf(stderr,"i1=%d error!\n", i1);	exit(EXIT_FAILURE);}
			for(int j0 = 0 ; j0 < 4*L[0]*L[1] ; j0++){
				if(fgets(buf, sizeof( buf ), fp_i2) == NULL){	fprintf(stderr,"fgets error! @ j0=%d\n", j0);	exit(EXIT_FAILURE);}
				double dummyR, dummyI;
				if(sscanf(buf,"%lf %lf", &dummyR, &dummyI)!=2){	fprintf(stderr,"real and imaginary data error!\n");	exit(EXIT_FAILURE);}
				Xi[j0] = dummyR + I*dummyI;
			}
			double waitA = 0.;	double waitB = 0.;
			for(int j1 = 0 ; j1 < L[1] ; j1++){
				for(int j0 = 0 ; j0 < L[0] ; j0++){
					int Axy = j0 + L[0]*j1;
					int Bxy = Axy + L[0]*L[1];
					Poynting[Axy] = 2*creal(conj(Xi[Axy])*Xi[Axy + 2*L[0]*L[1]]);
					waitA += Poynting[Axy];
					Poynting[Axy] /= du01[2*j0+1]*du01[2*L[0] + 2*j1+0];
					Poynting[Bxy] = 2*creal(conj(Xi[Bxy])*Xi[Bxy + 2*L[0]*L[1]]);
					waitB += Poynting[Bxy];
					Poynting[Bxy] /= du01[2*j0+0]*du01[2*L[0] + 2*j1+1];
				}
			}
			FILE *fp_o;
			fp_o = fopen(out_name,"w");
			fprintf(fp_o,"#mode_n =  %d, -Ey*Hx = %.5e, Ex*Hy = %.5e\n", mode_n, waitA, waitB);
			fprintf(fp_o,"#u_0[m], u_1[m], PoyintingVector_2\n");//	fprintf(stderr,"#mode_n =  %d\n", mode_n);
			for(int j1 = 0 ; j1 < L[1] ; j1++){
				int jj1 = (j1+1)%L[1];
				for(int j0 = 0 ; j0 < L[0] ; j0++){
					int jj0 = (j0+1)%L[0];
					double AveP = (Poynting[j0 + L[0]*j1] + Poynting[j0 + L[0]*j1 + L[0]*L[1]]);
					AveP += (Poynting[j0 + L[0]*jj1] + Poynting[jj0 + L[0]*j1 + L[0]*L[1]]);
					AveP *= 0.5;
					fprintf(fp_o,"%.5e, %.5e, %.5e\n", u01[2*j0+1], u01[2*j1+1 + 2*L[0]], AveP);
				}
				fprintf(fp_o,"\n");
			}
			if(fclose(fp_o) != 0) {	fprintf(stderr,"fclose error after output_file!\n");	exit(EXIT_FAILURE);}
		}
	}
	fprintf(stderr,"Finished error\n");
	finished:;
	SAFEFREE(Xi);	SAFEFREE(Poynting);
}
void rm_space( char *A );
//======================================================================
//  input_file1   ---- rm_space
//                           Last updated on Jul 26, 2020  
//======================================================================
void input_file1(FILE *fp_i1, const int L[2], const double inv_k0, double u01[], double du01[])
{
	char buf[BUFSIZE], Dummy[BUFSIZE];	// buffer for fgets
	int j_count = -2;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 0) { 
		rm_space(buf);
		if(strncmp(buf, "# xi0", 5) == 0){
			for(int i0 = 0 ; i0 < 2*L[0] ; i0++){
				double dummy0, dummy1, dummy2;
				if(fgets(Dummy, sizeof( Dummy ), fp_i1) != NULL && sscanf(Dummy,"%lf, %lf, %lf", &dummy0, &dummy1, &dummy2) == 3){
					u01[i0] = dummy1*inv_k0;
					du01[i0] = dummy2;//	du01[i0] = dummy2*inv_k0;
				}else{	fprintf(stderr,"u0 and du0 cannot be read in input_file1! i0=%d, %s\n", i0, Dummy);	exit(EXIT_FAILURE);}
			}
			j_count++;//-1
		}else if(strncmp(buf, "# xi1", 5) == 0){
			for(int i0 = 2*L[0] ; i0 < 2*(L[0]+L[1]) ; i0++){
				double dummy0, dummy1, dummy2;
				if(fgets(Dummy, sizeof( Dummy ), fp_i1) != NULL && sscanf(Dummy,"%lf, %lf, %lf", &dummy0, &dummy1, &dummy2) == 3){
					u01[i0] = dummy1*inv_k0;
					du01[i0] = dummy2;//	du01[i0] = dummy2*inv_k0;
				}else{	fprintf(stderr,"u1 and du1 cannot be read in input_file1! i0=%d, %s\n", i0, Dummy);	exit(EXIT_FAILURE);}
			}
			j_count++;//0
		}
	}
	if(j_count != 0) {	fprintf(stderr,"u01 and du01 cannot be read in input_file1!");	exit(EXIT_FAILURE);}
}
//======================================================================
//  input_file0 ---- rm_space, ScaleUnit, rm_comma                      
//                           Last updated on Jul 26, 2020  
//======================================================================
void rm_comma( char *A );
double ScaleUnit(char x[]);
void input_file0(FILE *fp_i0, char f_prefix[], char nonuniform[], int L[2], double *inv_k0)
{
	char buf[BUFSIZE], A[BUFSIZE];	// buffer for fgets
	int j_count = -5; int j, L_dummy;
	while(fgets(buf, sizeof( buf ), fp_i0) != NULL && j_count < 0) { 
		rm_space(buf);
		if(strncmp(buf, "inv_k0", 6) == 0 && sscanf(buf,"%*[^=] %*[=] %lf %s", inv_k0, A) == 2){//revised on 20200707
			*inv_k0 *= ScaleUnit(A);
			j_count++;//-4
		}else if(strncmp(buf, "mesh", 4)*strncmp(buf, "Mesh", 4)*strncmp(buf, "MESH", 4) == 0
		 && sscanf(buf,"%*[^=] %*[=] %d %*[^=] %*[=] %d",
		  &j, &L_dummy) == 2){
			if( j >= 0 && j < 2 && L_dummy >= 0) {
				L[j] = L_dummy;
				j_count++;//-3, -2
			}
		}else if(strncmp(buf, "xi2u", 4) == 0 && sscanf(buf,"%*[^=] %*[=] %s", nonuniform) == 1){
			j_count++;//-1
		}else if(strncmp(buf, "Prefix", 6) == 0 && sscanf(buf,"%*[^=] %*[=] %s", f_prefix) == 1){
			rm_comma(f_prefix);
			j_count++;//0
		}
	}
	if(j_count != 0) {	fprintf(stderr,"4 control commands cannot be read in input_file0!");	exit(EXIT_FAILURE);}
}
//======================================================================
//  double ScaleUnit                                                    
//                                       Last updated on Feb 26, 2017.  
//======================================================================
//#include <string.h>
double ScaleUnit(char x[])
{
	double scale;
	if (strncmp(x, "cm", 2) == 0) {
		scale = 1e-2;
	} else if (strncmp(x, "deg", 3) == 0) {
		scale = Pi/180.;
	} else if (strncmp(x, "km", 2) == 0) {
		scale = 1e3;
	} else if (strncmp(x, "m", 1) == 0) {
		if (strncmp(x, "mm", 2) == 0) {
			scale = 1e-3;// order of m, mm, micron is important!
		} else if (strncmp(x, "micron", 6) == 0) {
			scale = 1e-6;// order of m, mm, micron is important!
		} else {
			scale = 1e-0;// order of m, mm, micron is important!
		}
	} else if (strncmp(x, "nm", 2) == 0) {
		scale = 1e-9;
	} else if (strncmp(x, "um", 2) == 0) {
		scale = 1e-6;
	} else if (strncmp(x, "rad", 3) == 0) {
		scale = 1.;
	} else {
		fprintf(stderr,"ScaleUnit error!  x = %s\n", x);
		exit(1);
	}
	return(scale);
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
//  This program removes comma of end in characters.
//                                        Last updated on Jul 06, 2018  
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
