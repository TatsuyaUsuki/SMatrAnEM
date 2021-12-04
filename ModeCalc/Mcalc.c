//  Start of Mcalc.h
#include "header_macro.h"
#include "constant.h"
void input_file0(FILE *fp_i0, int BoundC[], char f_prefix[], char p_cross[], char InPrecision[], char Xi_ctrl[], double *omega0);
void main_calc(char **argv, const char data_name[BUFSIZE], const char out_name[BUFSIZE], const int BoundC[2], const char InPrecision[BUFSIZE], const char Xi_name[BUFSIZE], const double omega0);
//  End of Mcalc.h
//======================================================================
//   main ---- input_file0, main_calc             
//======================================================================
int main(int argc, char **argv)
{
	{	time_t timer_ini = time(0);	fprintf(stderr,"# Start time of mode calculation = %s\n", ctime(&timer_ini));}
	if(argc != 2) {	fprintf(stderr,"error: number of files \n");	exit(EXIT_FAILURE);}
	else if(strncmp(argv[1], "-v", 2) == 0 || strcmp(argv[1], "--version") == 0 ) {
		fprintf(stderr,"The '%s' creates wave-modes.\n", argv[0]);
		fprintf(stderr,"Version 21.08.12 is compiled at %s on %s.\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(stderr," Source code : '%s'\n Author      : Tatsuya Usuki\n URL         : http://www.smatran.org\n", __FILE__);
		fprintf(stderr," References  : 'Wave scattering in frequency domain' as 'Formulation.pdf' on Aug 14, 2019.\n");
		fprintf(stderr,"There is NO warranty.\n");
		exit(EXIT_SUCCESS);//normal end
	}
//------- begin reading fundamental parameters -------
	FILE *fp_i0;
	fp_i0 = fopen(argv[1],"r");
	if (fp_i0 == NULL){	fprintf(stderr,"open error!: open input-file1!\n");	exit(EXIT_FAILURE);}
	fprintf(stderr,"The 1st input file: %s\n",argv[1]);
	
	char f_prefix[BUFSIZE], p_cross[BUFSIZE*2], InPrecision[BUFSIZE], Xi_ctrl[BUFSIZE];
	int BoundC[2];
	double omega0;
	input_file0(fp_i0, BoundC, f_prefix, p_cross, InPrecision, Xi_ctrl, &omega0);//revised on 20200707
	if(fclose(fp_i0) != 0) {	fprintf(stderr,"fclose error after input_file!\n");	exit(EXIT_FAILURE);}
//-------  end reading fundamental parameters  -------
	for(int i1 = 0 ; i1 < 2 ; i1++){
		if(strcmp(p_cross, &p_cross[BUFSIZE]) == 0 && i1 == 1) {goto end_1st;}
		char data_name[BUFSIZE];
		snprintf(data_name,sizeof(data_name),"%s_Med%s.dat", f_prefix, &p_cross[i1*BUFSIZE]);// https://www.ipa.go.jp/security/awareness/vendor/programmingv1/b06_02.html
		char out_name[BUFSIZE], Xi_name[BUFSIZE];
		Xi_name[0] = '\0';
		if(i1 == 0){
			snprintf(out_name,sizeof(out_name),"%s_bMode.dat", f_prefix);
			if(strncmp(Xi_ctrl, "Y", 1) == 0 || strncmp(Xi_ctrl, "y", 1) == 0 ) {
				snprintf(Xi_name,sizeof(Xi_name),"%s_bXi.dat", f_prefix);
			}
		}
		else{
			snprintf(out_name,sizeof(out_name),"%s_tMode.dat", f_prefix);
			if(strncmp(Xi_ctrl, "Y", 1) == 0 || strncmp(Xi_ctrl, "y", 1) == 0 ) {
				snprintf(Xi_name,sizeof(Xi_name),"%s_tXi.dat", f_prefix);
			}
		}
		main_calc(argv, data_name, out_name, BoundC, InPrecision, Xi_name, omega0);//revised on 20200707
	}
	end_1st:;
}
//======================================================================
//  input_file0 ---- rm_space, ScaleUnit, rm_comma                      
//                           Last updated on Jul 07, 2020  
//======================================================================
void rm_space( char *A );
void rm_comma( char *A );
double ScaleUnit(char x[]);
void input_file0(FILE *fp_i0, int BoundC[2], char f_prefix[], char p_cross[], char InPrecision[], char Xi_ctrl[], double *omega0)
{
	char buf[BUFSIZE], A[BUFSIZE];	// buffer for fgets
	double lambda, inv_k0;
	int j_count = -5;
	while(fgets(buf, sizeof( buf ), fp_i0) != NULL && j_count < 0) { 
		rm_space(buf);
		if(strncmp(buf, "inv_k0", 6) == 0 && sscanf(buf,"%*[^=] %*[=] %lf %s", &inv_k0, A) == 2){//revised on 20200707
			inv_k0 *= ScaleUnit(A);
			j_count++;//-4
		}else if(strncmp(buf, "MediumPara", 10) == 0 && sscanf(buf,"%*[^=] %*[=] %lf %s", &lambda, A) == 2){//revised on 20200707
			lambda *= ScaleUnit(A);
			j_count++;//-3
		}else if(strncmp(buf, "Prefix", 6) == 0 && sscanf(buf,"%*[^=] %*[=] %s", f_prefix) == 1){
			rm_comma(f_prefix);
			j_count++;//-2
		}else if(strncmp(buf, "BoundaryCondition", 17) == 0 && sscanf(buf,"%*[^=] %*[=] %d %*[^=] %*[=] %d", &BoundC[0], &BoundC[1]) == 2) {
			j_count++;//-1
		}else if(strncmp(buf, "ModeP", 5) == 0 && sscanf(buf,"%*[^=] %*[=] %s %*[^=] %*[=] %s %*[^=] %*[=] %s", Xi_ctrl, InPrecision, p_cross) == 3) {
			rm_comma(InPrecision);
			rm_comma(p_cross);
			if(sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %*[^=] %*[=] %*[^=] %*[=] %s", &p_cross[BUFSIZE]) != 1){sscanf(p_cross,"%s", &p_cross[BUFSIZE]);}
			else{rm_comma(&p_cross[BUFSIZE]);}
			j_count++;//0
		}
	}
	*omega0 = (Pi2*inv_k0)/(lambda);//revised on 20200707
	if(j_count != 0) {	fprintf(stderr,"3 control commands cannot be read in input_file1!");	exit(EXIT_FAILURE);}
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
//======================================================================
//  main_calc ---- read_para, set_V, calc_Xi                  
//                                       Last updated on Jul 07, 2020.  
//======================================================================
void calc_Xi(const int L[2], const int BoundC[2], const char InPrecision[BUFSIZE], 
	const double omega, long count_num[2], double *norm, double *Xi_Norm, 
	double complex theta_m[6*L[0]*L[1]], double complex Xi_m[4*L[0]*L[1]*4*L[0]*L[1]]);
long set_V(FILE *fp_i, const int L[2], const double omega, double complex Vcomp[]);
void read_para(FILE *fp_i, const char data_name[BUFSIZE], char comment[BUFSIZE], double *du_2, double *k_0, int L[2], double *omega);
void main_calc(char **argv, const char data_name[BUFSIZE], const char out_name[BUFSIZE], const int BoundC[2], const char InPrecision[BUFSIZE], const char Xi_name[BUFSIZE], const double omega0)
{
	time_t timer_s, timer_f;
	timer_s = time(0);
	clock_t startClock = clock();
	FILE *fp_i;
	fp_i = fopen(data_name,"r");	if (fp_i == NULL){	fprintf(stderr,"open error!: open data-file!\n");	exit(EXIT_FAILURE);}
	char comment[BUFSIZE];
	int L[2];
	double du_2, k_0, omega;
	read_para(fp_i, data_name, comment, &du_2, &k_0, L, &omega);
	if(fabs(1. - omega0/omega) > 1.e-10){fprintf(stderr,"omega = %.10e != omega0 = %.10e\n", omega, omega0); omega = omega0;}//revised on 20200707
	double norm, Xi_Norm;
	long count_num[2];//{comp_count, num_real};
	double complex *theta_m = calloc(6*L[0]*L[1],sizeof(double complex));//Note that region from &theta_m[4*L[0]*L[1]] to &theta_m[6*L[0]*L[1]-1] is used as working region.
	count_num[0] = set_V(fp_i, L, omega, theta_m);//set_V(fp_i, L, omega, Vcomp)
	if(fclose(fp_i) != 0) {	fprintf(stderr,"fclose error after set_V!\n");	exit(EXIT_FAILURE);}
	double complex *Xi_m = calloc(4*L[0]*L[1]*4*L[0]*L[1],sizeof(double complex));
	calc_Xi(L, BoundC, InPrecision, omega, count_num, &norm, &Xi_Norm, theta_m, Xi_m);
	if(Xi_name[0] != '\0'){
		FILE *fp_o2;
		timer_f = time(0);
		fp_o2 = fopen(Xi_name,"w");	if (fp_o2 == NULL){	fprintf(stderr,"open error!: open Xi-file!\n");	exit(EXIT_FAILURE);}
		fprintf(fp_o2,"# This file      : ' %s ' on %s", Xi_name, ctime(&timer_f));
		fprintf(fp_o2,"# Reference file : ' %s '\n", out_name);
		for(long i_mod = 0 ; i_mod < 2*L[0]*L[1] ; i_mod++){
			fprintf(fp_o2,"ReXi[%ld] ImXi[%ld] ReXi[%ld] ImXi[%ld]\n", i_mod, i_mod, i_mod+2*L[0]*L[1], i_mod+2*L[0]*L[1]);
			for(long j = 0 ; j < 4*L[0]*L[1] ; j++){
				fprintf(fp_o2,"%.20e %.20e %.20e %.20e\n", creal(Xi_m[j+4*L[0]*L[1]*i_mod]), cimag(Xi_m[j+4*L[0]*L[1]*i_mod]), creal(Xi_m[j+4*L[0]*L[1]*(i_mod+2*L[0]*L[1])]), cimag(Xi_m[j+4*L[0]*L[1]*(i_mod+2*L[0]*L[1])]));
			}
		}
		if(fclose(fp_o2) != 0) {	fprintf(stderr,"fclose error after Xi_data!\n");	exit(EXIT_FAILURE);}
	}
	SAFEFREE(Xi_m);
	FILE *fp_o;
	fp_o = fopen(out_name,"w");	if (fp_o == NULL){	fprintf(stderr,"open error!: open o-file!\n");	exit(EXIT_FAILURE);}
	timer_f = time(0);
	clock_t endClock = clock();
	fprintf(fp_o,"# This file          : ' %s ' on %s", out_name, ctime(&timer_f));
	fprintf(fp_o,"%s\n", comment);
	fprintf(fp_o,"# Exec. file         : ' %s ' when '%s' was compiled at %s on %s by C-version %ld\n", argv[0], __FILE__, __TIME__, __DATE__, __STDC_VERSION__);
	fprintf(fp_o,"# Calc. time         : Exec. time = %1.0f s, CPU time = %1.6f s\n", difftime(timer_f,timer_s), (endClock - startClock)/(double)CLOCKS_PER_SEC);
	fprintf(fp_o,"# Parameter file     : ' %s ' for prefix, boundary conditions and internal precision\n", argv[1]);
	fprintf(fp_o,"# Normalization para.: k_0[m^-1] = %.20e\n", k_0);
	fprintf(fp_o,"# Omega              : k/k_0 = %.20e\n", omega);
	fprintf(fp_o,"# Lattice numbers    : L0 = %d, L1 = %d\n", L[0], L[1]);
	fprintf(fp_o,"# Boundary conditions: B0 = %d, B1 = %d\n", BoundC[0], BoundC[1]);
	fprintf(fp_o,"# Internal precision : %s\n", InPrecision);
	fprintf(fp_o,"# Check of solutions : ||mbb*eigenE/eiSin - eigenH||_F/(||mbb*eigenE/eiSin||_F+||eigenH||_F) = %11.5e\n", 
		norm);
	fprintf(fp_o,"# Check of bi-ortho. : ||(tilde Xi)^+ Xi - E||_F = %11.5e\n", Xi_Norm);
	if(count_num[0] == 0){fprintf(fp_o,"# Matrix info.       : real matrix\n");}
	else{fprintf(fp_o,"# Matrix info.       : complex matrix\n");}
	fprintf(fp_o,"# Eigen-value info.  : real num./2 = %ld, comp. num./2 = %ld\n", count_num[1], 2*L[0]*L[1]-count_num[1]);
	fprintf(fp_o,"# Cell length        : k_0*Du_2 = %.20e\n", du_2);
	for(long j = 0 ; j < 2*L[0]*L[1] ; j++){
		fprintf(fp_o,"theta[%ld]=(%.20e, %.20e), theta[%ld]=(%.20e, %.20e)\n", j, creal(theta_m[j]), cimag(theta_m[j]), j+2*L[0]*L[1], creal(theta_m[j+2*L[0]*L[1]]), cimag(theta_m[j+2*L[0]*L[1]]));
	}
	if(fclose(fp_o) != 0) {	fprintf(stderr,"fclose error after output_data!\n");	exit(EXIT_FAILURE);}
	SAFEFREE(theta_m);
}
void degeneXi(const long N, const long num_real, const double complex theta_m[2*N], double complex Xi_m[2*N*2*N], double complex CWork[]);
double checkXi(const long N, const long num_real, const double complex theta_m[2*N], const double complex Xi_m[2*N*2*N], double complex CWork[]);
void normXi(const long N, const long num_real, const long icolumn, double complex Xi_m[2*N*2*N]);
void sortXi(const long N, const long num_real, double complex theta_m[2*N], double complex Xi_m[2*N*2*N]);
void makeXi(const long N, const double complex eiSin[N], double complex eigenH[N*N], double complex eigenE[N*N],
	double complex Theta_m[2*N], double complex Xi_m[2*N*2*N]);
void checkM(const long N, const double complex maa[N*N], const double complex mbb[N*N]);
double norm_calc(const long N, const double complex A[N*N], const double complex E[N], const double complex V1[N*N], const double complex V2[N*N], double complex CWork[]);
void makeE(const int L[2], const double complex maa[2*L[0]*L[1]*2*L[0]*L[1]], const double complex eiSin[2*L[0]*L[1]],
 const double complex eigenH[2*L[0]*L[1]*2*L[0]*L[1]], double complex eigenE[2*L[0]*L[1]*2*L[0]*L[1]]);
long sortSin(const int L[2], double complex eiSin[2*L[0]*L[1]], double complex eigenH[2*L[0]*L[1]*2*L[0]*L[1]]);
void makeSin(const int L[2], const double complex eiSin2[2*L[0]*L[1]], double complex eiSin[2*L[0]*L[1]]);
void realizeH(const int L[2], const double complex eiSin2[2*L[0]*L[1]], double complex eigenH[2*L[0]*L[1]*2*L[0]*L[1]]);
void ReigenVcalc(const char InPrecision[BUFSIZE], const int L[2], const double complex mbbmaa[2*L[0]*L[1]*2*L[0]*L[1]], double complex eiSin2[2*L[0]*L[1]], double complex eigenH[2*L[0]*L[1]*2*L[0]*L[1]]);
void eigenVcalc(const char InPrecision[BUFSIZE], const int L[2], const double complex mbbmaa[2*L[0]*L[1]*2*L[0]*L[1]], double complex eiSin2[2*L[0]*L[1]], double complex eigenH[2*L[0]*L[1]*2*L[0]*L[1]]);
void set_M(const int BoundC[2], const int L[2], const double complex Vcomp[6*L[0]*L[1]], 
	double complex maa[2*L[0]*L[1]*2*L[0]*L[1]], double complex mbb[2*L[0]*L[1]*2*L[0]*L[1]], double complex CWork[], double DWork[]);
void set_Msub(const int L[2], double complex mbbmaa[2*L[0]*L[1]*2*L[0]*L[1]], double complex maa[2*L[0]*L[1]*2*L[0]*L[1]], double complex mbb[2*L[0]*L[1]*2*L[0]*L[1]]);
//======================================================================
//  calc_Xi ---- set_M, checkM, ReigenVcalc, eigenVcalc, realizeH, 
//               makeSin, sortSin, makeE, norm_calc, makeXi, 
//               sortXi, normXi, chechXi    
//                                        Last updated on Aug 20, 2019  
//======================================================================
void calc_Xi(const int L[2], const int BoundC[2], const char InPrecision[BUFSIZE], 
	const double omega, long count_num[2], double *norm, double *Xi_Norm, 
	double complex theta_m[6*L[0]*L[1]], double complex Xi_m[4*L[0]*L[1]*4*L[0]*L[1]])
{
	long N = 2*L[0]*L[1];
	time_t timer_s, timer_f;	clock_t startClock, endClock;
	//--- maa = &Xi_m[0], mbb = &Xi_m[N*N*1], mbbmaa = &Xi_m[N*N*3], eigenH = &Xi_m[N*N*2], eigenE = &Xi_m[N*N*3] ---
	timer_s = time(0);
	double complex *CWork = calloc(2*N*2*N*2,sizeof(double complex));//Max. setting of CWork for using check_Xi//calloc(2*N*N,sizeof(double complex));
	{
		double complex *maa = &Xi_m[0];	double complex *mbb = &Xi_m[N*N*1];
		{
			double complex *Vcomp = &theta_m[0];
			double *DWork = calloc(N*N,sizeof(double));//revised on 20190114//long NA1 = L[0]*L[1];	long Nshift = NA1*NA1;//calloc(Nshift*4,sizeof(double));
			set_M(BoundC, L, Vcomp, maa, mbb, CWork, DWork);
			SAFEFREE(DWork);
		}
		double complex *eiSin2 = &theta_m[0];	double complex *eigenH = &Xi_m[N*N*2];
		{
			double complex *mbbmaa = &Xi_m[N*N*3];
			set_Msub(L, mbbmaa, maa, mbb);//set_Msub(L, mbbmaa, maa, mbb);
			timer_f = time(0);fprintf(stderr,"set_M time = %1.0f s\n", difftime(timer_f,timer_s));
			//------------------------------------------------------------------------------------------------------------
			timer_s = time(0);
			if(count_num[0] == 0){checkM(N, maa, mbb);};//checkM(N, maa, mbb)
			timer_f = time(0);fprintf(stderr,"checkM time = %1.0f s\n", difftime(timer_f,timer_s));
			//--------------------------------  eiSin2 = &theta_m[0], eiSin = &theta_m[N*2]  --------------------------------
			timer_s = time(0);	startClock = clock();
			if(count_num[0] == 0){ReigenVcalc(InPrecision, L, mbbmaa, eiSin2, eigenH);}
			else{eigenVcalc(InPrecision, L, mbbmaa, eiSin2, eigenH);}
			fprintf(stderr,"||MM*H/(k^2)-H||_F/(||MM*H/(k^2)||_F+||H||_F)= %.5e\n", norm_calc(N, &Xi_m[N*N*3], &theta_m[0], &Xi_m[N*N*2], &Xi_m[N*N*2], CWork));
			endClock = clock();	timer_f = time(0);	fprintf(stderr,"ReigenVcalc or eigenVcalc time = %1.0f s, CPU-time = %1.6f s\n", difftime(timer_f,timer_s), (endClock - startClock)/(double)CLOCKS_PER_SEC);
		}//------------------------------------------------------------------------------------------------------------
		timer_s = time(0);
		if(count_num[0] == 0){realizeH(L, eiSin2, eigenH);};//realizeH(L, eiSin2, eigenH)
		timer_f = time(0);fprintf(stderr,"realizeH time = %1.0f s\n", difftime(timer_f,timer_s));
		//------------------------------------------------------------------------------------------------------------
		timer_s = time(0);
		double complex *eiSin = &theta_m[N*2];
		makeSin(L, eiSin2, eiSin);
		timer_f = time(0);fprintf(stderr,"makeSin time = %1.0f s\n", difftime(timer_f,timer_s));
		//------------------------------------------------------------------------------------------------------------
		timer_s = time(0);
		count_num[1] = sortSin(L, eiSin, eigenH);
		timer_f = time(0);fprintf(stderr,"sortSin time = %1.0f s\n", difftime(timer_f,timer_s));
		{//------------------------------------------------------------------------------------------------------------
			timer_s = time(0);	startClock = clock();
			double complex *eigenE = &Xi_m[N*N*3];
			makeE(L, maa, eiSin, eigenH, eigenE);
			endClock = clock();	timer_f = time(0);fprintf(stderr,"makeE time = %1.0f s, CPU-time = %1.6f s\n", difftime(timer_f,timer_s), (endClock - startClock)/(double)CLOCKS_PER_SEC);
			//------------------------------------------------------------------------------------------------------------
			timer_s = time(0);	startClock = clock();
			*norm = norm_calc(N, mbb, eiSin, eigenE, eigenH, CWork);
			endClock = clock();	timer_f = time(0);	fprintf(stderr,"norm_calc time = %1.0f s, CPU-time = %1.6f s\n", difftime(timer_f,timer_s), (endClock - startClock)/(double)CLOCKS_PER_SEC);
			//------------------------------------------------------------------------------------------------------------
			timer_s = time(0);
			makeXi(N, eiSin, eigenH, eigenE, theta_m, Xi_m);
			timer_f = time(0);fprintf(stderr,"makeXi time = %1.0f s\n", difftime(timer_f,timer_s));
		}
	}//------------------------------------------------------------------------------------------------------------
	timer_s = time(0);
	sortXi(N, count_num[1], theta_m, Xi_m);
	timer_f = time(0);fprintf(stderr,"sortXi time = %1.0f s\n", difftime(timer_f,timer_s));
	//------------------------------------------------------------------------------------------------------------
	timer_s = time(0);
	for(long icolumn = 0 ; icolumn < N ; icolumn++){	normXi(N, count_num[1], icolumn, Xi_m);}
	timer_f = time(0);fprintf(stderr,"normXi time = %1.0f s\n", difftime(timer_f,timer_s));
	//------------------------------------------------------------------------------------------------------------
	timer_s = time(0);
	startClock = clock();
	degeneXi(N, count_num[1], theta_m, Xi_m, CWork);
	endClock = clock();
	timer_f = time(0);
	fprintf(stderr,"degeneXi time = %1.0f s, CPU-time = %1.6f s\n", difftime(timer_f,timer_s), (endClock - startClock)/(double)CLOCKS_PER_SEC);
	//------------------------------------------------------------------------------------------------------------
	timer_s = time(0);
	startClock = clock();
	*Xi_Norm = checkXi(N, count_num[1], theta_m, Xi_m, CWork);
	fprintf(stderr,"||(tilde Xi)^+ Xi - E||_F = %.5e\n", *Xi_Norm);
	SAFEFREE(CWork);
	endClock = clock();
	timer_f = time(0);
	fprintf(stderr,"checkXi time = %1.0f s, CPU-time = %1.6f s\n", difftime(timer_f,timer_s), (endClock - startClock)/(double)CLOCKS_PER_SEC);
}
//======================================================================
//  checkXi ---- makePsi                 Last updated on Aug 20, 2019.  
//======================================================================
extern void zgemm_(char *transa, char *transb, const int *m, const int *n, const int *k,
	const double complex *alpha, const double complex *A, const int *ldA, 
	const double complex *B, const int *ldB,
	const double complex *beta , double complex *C, const int *ldC);
extern double dznrm2_(const int *n, const double complex *x, const int *incx);
void makePsi(const long N, const long num_real, const long mode0, const double complex theta_m[2*N], const double complex Xi_m[2*N*2*N], double complex Psi[2*N*2*N]);
double checkXi(const long N, const long num_real, const double complex theta_m[2*N], const double complex Xi_m[2*N*2*N], double complex CWork[])//revised on 20181011
{
	double norm = 0.;
	double complex *Psi = &CWork[0];//calloc(2*N*2*N,sizeof(double complex));
	for(long mode0 = 0 ; mode0 < N ; mode0++){
		makePsi(N, num_real, mode0, theta_m, Xi_m, Psi);
	}
	/*for(long mode0 = 0 ; mode0 < num_real ; mode0++){
		for(long mode1 = 0 ; mode1 < N ; mode1++){
			Psi[mode1   + 2*N*mode0] = Xi_m[mode1+N + 2*N*mode0];
			Psi[mode1+N + 2*N*mode0] = Xi_m[mode1   + 2*N*mode0];
			Psi[mode1   + 2*N*(mode0+N)] = -Xi_m[mode1+N + 2*N*(mode0+N)];
			Psi[mode1+N + 2*N*(mode0+N)] = -Xi_m[mode1   + 2*N*(mode0+N)];
		}
	}
	for(long mode0 = num_real ; mode0 < N ; mode0++){
		for(long mode1 = 0 ; mode1 < N ; mode1++){
			Psi[mode1   + 2*N*(mode0+N)]= Xi_m[mode1+N + 2*N*(mode0  )];
			Psi[mode1+N + 2*N*(mode0+N)]= Xi_m[mode1   + 2*N*(mode0  )];
			Psi[mode1   + 2*N*(mode0  )]= Xi_m[mode1+N + 2*N*(mode0+N)];
			Psi[mode1+N + 2*N*(mode0  )]= Xi_m[mode1   + 2*N*(mode0+N)];
		}
	}*/
	{// see " http://azalea.s35.xrea.com/blas/index.html "
		int iN = 2*N;
		double complex alpha = 1.;
		double complex beta = 0.;
		double complex *C = &CWork[2*N*2*N];//calloc(iN*iN,sizeof(double complex));
		zgemm_("C", "N", &iN, &iN, &iN, &alpha, &Psi[0], &iN, &Xi_m[0], &iN, &beta, &C[0],&iN);// "H" is an illegal value!
		for(int i = 0 ; i < iN ; i++){C[i+iN*i]-=1.;}
		/*for(int j = 0 ; j < iN ; j++){
			for(int i = 0 ; i < iN ; i++){
				if(fabs(C[i+iN*j]) > sqrt(sqrt((double) iN)*DBL_EPSILON)){
					if(fabs((theta_m[i]-theta_m[j])/(theta_m[i]+theta_m[j])) < sqrt((double) iN)*DBL_EPSILON){
						fprintf(stderr,"Two modes show accidental degeneracy: |Psi[%d]^+ Xi_m[%d] - E| = %.5e,", i, j, fabs(C[i+iN*j]));
						fprintf(stderr," and |(theta_m[%d]-theta_m[%d])/(theta_m[%d]+theta_m[%d])| = %.5e.\n", i, j,i, j, fabs((theta_m[i]-theta_m[j])/(theta_m[i]+theta_m[j])));
					}else{
						fprintf(stderr,"Check orthogonality for two modes: |Psi[%d]^+ Xi_m[%d] - E| = %.5e > 1e-3.\n", i, j, fabs(C[i+iN*j]));
					}
				}// revised on 20190819
			}
		}*/
		int iN2 = iN*iN;//revised on 20190107
		int incx = 1;
		norm = dznrm2_(&iN2,&C[0],&incx);
//		SAFEFREE(C);
	}// http://www.netlib.org/lapack/explore-html/d7/d76/zgemm_8f.html
	/*{
		for(long mode0 = 0 ; mode0 < 2*N ; mode0++){
			for(long mode1 = 0 ; mode1 < 2*N ; mode1++){
				double complex normZ = 0. + I*0.;
				for(long isum = 0 ; isum < 2*N ; isum++){
					normZ += conj(Psi[isum + 2*N*mode0]) * Xi_m[isum + 2*N*mode1];
				}
				normZ -= (mode0 == mode1);//if(mode0 == mode1) {normZ -= 1.;}
				norm += creal(normZ)*creal(normZ) + cimag(normZ)*cimag(normZ);
			}
		}
		norm = sqrt(norm);
	}*/
//	SAFEFREE(Psi);
	return(norm);
}
//======================================================================
//  degeneXi ---- normXi, makePsi        Last updated on Aug 08, 2021.  
//======================================================================
void normXi(const long N, const long num_real, const long icolumn, double complex Xi_m[2*N*2*N]);
void makePsi(const long N, const long num_real, const long mode0, const double complex theta_m[2*N], const double complex Xi_m[2*N*2*N], double complex Psi[2*N*2*N]);
void degeneXi(const long N, const long num_real, const double complex theta_m[2*N], double complex Xi_m[2*N*2*N], double complex CWork[])
{
	double complex *Psi = &CWork[0];
	long N2 = 2*N;
	for(long j = 0 ; j < N-1 ; j++){
		for(long i = j+1 ; i < N ; i++){
			if(fabs((theta_m[j]-theta_m[i])/(theta_m[j]+theta_m[i])) < FLT_EPSILON){
//			if(fabs((theta_m[j]-theta_m[i])/(theta_m[j]+theta_m[i])) < sqrt((double) N2)*DBL_EPSILON){
				//sum_ij = <i|j>, sum_ji = <j|i> for i > j.
				//|i_mod> = |i> - sum_ji*|j>.
				// <i_mod|i_mod> = (<i| - conj(sum_ji)*<j|)(|i> - sum_ji*|j>) = 1 - |sum_ji|^2 - sum_ji*sum_ij + |sum_ji|^2 = 1 - sum_ji*sum_ij.
				double complex sum_ij0; sum_ij0 = 0.;
				double complex sum_ijN; sum_ijN = 0.;
				makePsi(N, num_real, i, theta_m, Xi_m, Psi);
				for(long k = 0 ; k < N2 ; k++){
					sum_ij0 += conj(Psi[k+N2*i])*Xi_m[k+N2*j];
					sum_ijN += conj(Psi[k+N2*(i+N)])*Xi_m[k+N2*(j+N)];
				}
				double complex sum_ji0; sum_ji0 = 0.;
				double complex sum_jiN; sum_jiN = 0.;
				makePsi(N, num_real, j, theta_m, Xi_m, Psi);
				for(long k = 0 ; k < N2 ; k++){
					sum_ji0 += conj(Psi[k+N2*j])*Xi_m[k+N2*i];
					sum_jiN += conj(Psi[k+N2*(j+N)])*Xi_m[k+N2*(i+N)];
				}
				double complex a0, b0;
				{
					double nor0 = 1./sqrt(fabs(1. - sum_ji0*sum_ij0));
					a0 = nor0 + 0.*I;
					b0 = -a0 * sum_ji0;
				}
				double complex aN, bN;
				{
					double norN= 1./sqrt(fabs(1. - sum_jiN*sum_ijN));
					aN = norN + 0.*I;
					bN = -aN * sum_jiN;
				}
				for(long k = 0 ; k < N2 ; k++){
					Xi_m[k+N2*i] = a0*Xi_m[k+N2*i] + b0*Xi_m[k+N2*j];
					Xi_m[k+N2*(i+N)] = aN*Xi_m[k+N2*(i+N)] + bN*Xi_m[k+N2*(j+N)];
				}
				normXi(N, num_real, i, Xi_m);
				double complex check_ji0; check_ji0 = 0.;
				double complex check_jiN; check_jiN = 0.;
				for(long k = 0 ; k < N2 ; k++){
					check_ji0 += conj(Psi[k+N2*j])*Xi_m[k+N2*i];
					check_jiN += conj(Psi[k+N2*(j+N)])*Xi_m[k+N2*(i+N)];
				}
				fprintf(stderr,"Accidental degeneracy: |(theta_m[%ld]-theta_m[%ld])/(theta_m[%ld]+theta_m[%ld])| = %.5e,", j, i, j, i, fabs((theta_m[j]-theta_m[i])/(theta_m[j]+theta_m[i])));
				fprintf(stderr," and |(tilde Xi[%ld])^+ Xi[%ld]| = %.5e --> %.5e\n", j, i, fabs(sum_ji0), fabs(check_ji0));
			}
		}
	}
}
/* Old version on Aug 20, 2019.
{
	double complex *Psi = &CWork[0];
	long N2 = 2*N;
	for(long j = 0 ; j < N-1 ; j++){
		makePsi(N, num_real, j, theta_m, Xi_m, Psi);
		{
			double complex sum_jj0; sum_jj0 = 0.;
			double complex sum_jjN; sum_jjN = 0.;
			for(long k = 0 ; k < N2 ; k++){sum_jj0 += conj(Psi[k+N2*j])*Xi_m[k+N2*j];	sum_jjN += conj(Psi[k+N2*(j+N)])*Xi_m[k+N2*(j+N)];}
			sum_jj0 -= 1.;
			sum_jjN -= 1.;
			if(fabs(sum_jj0) > 0.01 || fabs(sum_jjN) > 0.01){	fprintf(stderr,"Normalization error in degeneXi!: S[%ld]=%.5e, S[%ld]=%.5e\n",j,fabs(sum_jj0),j+N,fabs(sum_jjN));	exit(EXIT_FAILURE);}
		}
		for(long i = j+1 ; i < N ; i++){
			double complex sum_ji0; sum_ji0 = 0.;
			double complex sum_jiN; sum_jiN = 0.;
			for(long k = 0 ; k < N2 ; k++){
				sum_ji0 += conj(Psi[k+N2*j])*Xi_m[k+N2*i];
				sum_jiN += conj(Psi[k+N2*(j+N)])*Xi_m[k+N2*(i+N)];
			}
			if(fabs((theta_m[j]-theta_m[i])/(theta_m[j]+theta_m[i])) < sqrt((double) N2)*DBL_EPSILON){
				fprintf(stderr,"Accidental degeneracy: |(theta_m[%ld]-theta_m[%ld])/(theta_m[%ld]+theta_m[%ld])| = %.5e,", j, i, j, i, fabs((theta_m[j]-theta_m[i])/(theta_m[j]+theta_m[i])));
				fprintf(stderr," and |(tilde Xi[%ld])^+ Xi[%ld]| = %.5e --> zero\n", j, i, fabs(sum_ji0));
				for(long k = 0 ; k < N2 ; k++){
					Xi_m[k+N2*i] -= sum_ji0*Xi_m[k+N2*j];
					Xi_m[k+N2*(i+N)] -= sum_jiN*Xi_m[k+N2*(j+N)];
				}
				normXi(N, num_real, i, Xi_m);
			}
		}
	}
}
*/
//======================================================================
//  makePsi                              Last updated on Aug 20, 2019.  
//======================================================================
void makePsi(const long N, const long num_real, const long mode0, const double complex theta_m[2*N], const double complex Xi_m[2*N*2*N], double complex Psi[2*N*2*N])
{
	if(0 <= mode0 && mode0 < num_real){
		for(long mode1 = 0 ; mode1 < N ; mode1++){
			Psi[mode1   + 2*N*mode0] = Xi_m[mode1+N + 2*N*mode0];
			Psi[mode1+N + 2*N*mode0] = Xi_m[mode1   + 2*N*mode0];
			Psi[mode1   + 2*N*(mode0+N)] = -Xi_m[mode1+N + 2*N*(mode0+N)];
			Psi[mode1+N + 2*N*(mode0+N)] = -Xi_m[mode1   + 2*N*(mode0+N)];
		}
	}else if(num_real <= mode0 && mode0 < N){
		for(long mode1 = 0 ; mode1 < N ; mode1++){
			Psi[mode1   + 2*N*(mode0+N)]= Xi_m[mode1+N + 2*N*(mode0  )];
			Psi[mode1+N + 2*N*(mode0+N)]= Xi_m[mode1   + 2*N*(mode0  )];
			Psi[mode1   + 2*N*(mode0  )]= Xi_m[mode1+N + 2*N*(mode0+N)];
			Psi[mode1+N + 2*N*(mode0  )]= Xi_m[mode1   + 2*N*(mode0+N)];
		}
	}else{	fprintf(stderr,"mode0 error in makePsi!\n");	exit(EXIT_FAILURE);	}
}
//======================================================================
//  normXi                               Last updated on Aug 20, 2019.  
//======================================================================
void normXi(const long N, const long num_real, const long icolumn, double complex Xi_m[2*N*2*N])
{
	if(0 <=icolumn && icolumn < num_real){
		for(long fb = 0 ; fb < 2 ; fb++){
			double complex normZ = 0. + I*0.;
			for(long irow = 0 ; irow < N ; irow++){
				normZ += conj(Xi_m[irow + 2*N*(icolumn+fb*N)]) * Xi_m[irow+N + 2*N*(icolumn+fb*N)];
				normZ += conj(Xi_m[irow+N + 2*N*(icolumn+fb*N)]) * Xi_m[irow + 2*N*(icolumn+fb*N)];
			}
			if(fabs(cimag(normZ)) > DBL_EPSILON) {fprintf(stderr,"Warning of imaginary flow in normXi!: normZ = (%.5e,%.5e) at mode%ld\n", creal(normZ),cimag(normZ),icolumn+fb*N);}//revised on 20181110
			if((0.5 - fb)*creal(normZ) < 0.) {	fprintf(stderr,"flow direction error in normXi! creal(normZ)=%.5e, icolumn/num_real=%ld/%ld, fb=%ld\n", creal(normZ), icolumn, num_real, fb);	exit(EXIT_FAILURE);}//revised on 20200810
			double normR = 1./sqrt(fabs(creal(normZ)));
			for(long irow = 0 ; irow < 2*N ; irow++){
				Xi_m[irow + 2*N*(icolumn+fb*N)] *= normR;
			}
		}
	}else if(num_real <= icolumn && icolumn < N){
			double complex normZ = 0. + I*0.;
			for(long irow = 0 ; irow < N ; irow++){
				normZ += conj(Xi_m[irow + 2*N*(icolumn+N)]) * Xi_m[irow+N + 2*N*(icolumn)];
				normZ += conj(Xi_m[irow+N + 2*N*(icolumn+N)]) * Xi_m[irow + 2*N*(icolumn)];
			}
//			normZ = 1./sqrt(normZ);
			normZ = sqrt(normZ);	normZ = conj(normZ)/(creal(normZ)*creal(normZ)+cimag(normZ)*cimag(normZ));//revised on 20190619
			for(long irow = 0 ; irow < 2*N ; irow++){
				Xi_m[irow + 2*N*(icolumn)] *= normZ;
				Xi_m[irow + 2*N*(icolumn+N)] *= conj(normZ);
			}
	}else{	fprintf(stderr,"icolumn error in normXi!\n");	exit(EXIT_FAILURE);	}
}
//======================================================================
//  sortXi ---- sort_cpair                Last updated on Aug 12, 2021.  
//======================================================================
void sort_cpair(const long N, const long i, const long j, double complex theta_m[2*N], double complex Xi_m[2*N*2*N]);
void sortXi(const long N, const long num_real, double complex theta_m[2*N], double complex Xi_m[2*N*2*N])
{
	for(long icolumn = num_real ; icolumn < N-1 ; icolumn++){
		if(fabs(creal(theta_m[icolumn+N])) >  FLT_EPSILON && fabs(creal(theta_m[icolumn]+theta_m[icolumn+N])) < FLT_EPSILON){
			for(long isearch = icolumn+1 ; isearch < N ; isearch++){
				if(fabs(creal(theta_m[isearch+N]-theta_m[icolumn])) <  FLT_EPSILON
				&& fabs(cimag(theta_m[isearch+N]+theta_m[icolumn])) < FLT_EPSILON){
					sort_cpair(N, icolumn+N, isearch+N, theta_m, Xi_m);
					goto theta_conj;
				}
			}
			fprintf(stderr,"no-pair complex theta_m[%ld] = (%.5e,%.5e)\n", icolumn+N, creal(theta_m[icolumn+N]), cimag(theta_m[icolumn+N]));
			theta_conj:;
		}
	}
	//------- sort theta_m added on 20210818 -------
	for(long icolumn = 0 ; icolumn < num_real-1 ; icolumn++){
		for(long isearch = icolumn+1 ; isearch < num_real ; isearch++){
			if(  creal(theta_m[isearch]) > creal(theta_m[icolumn]) + FLT_EPSILON){
				sort_cpair(N, icolumn, isearch, theta_m, Xi_m);
				sort_cpair(N, icolumn+N, isearch+N, theta_m, Xi_m);
			}
		}
	}
	for(long icolumn = num_real ; icolumn < N-1 ; icolumn++){
		for(long isearch = icolumn+1 ; isearch < N ; isearch++){
			if( cimag(theta_m[isearch]) + FLT_EPSILON < cimag(theta_m[icolumn]) ){
				sort_cpair(N, icolumn, isearch, theta_m, Xi_m);
				sort_cpair(N, icolumn+N, isearch+N, theta_m, Xi_m);
			}
		}
	}
	for(long icolumn = num_real ; icolumn < N ; icolumn++){//revised on 20211107
		if(fabs(theta_m[icolumn]-conj(theta_m[icolumn+N])) > FLT_EPSILON){
			for(long jcolumn = icolumn + 1 ; jcolumn < N ; jcolumn++){
				if(fabs(theta_m[icolumn]-conj(theta_m[jcolumn+N])) <= FLT_EPSILON){	sort_cpair(N, icolumn+N, jcolumn+N, theta_m, Xi_m);	goto checkconj;}
			}
			checkconj:;
		}
	}
	//-------     End of sort theta_m     -------
	//------- check theta_m added on 20210816 -------
	for(long icolumn = 0 ; icolumn < num_real-1 ; icolumn++){
		if(  creal(theta_m[icolumn+1]) > creal(theta_m[icolumn]) + FLT_EPSILON){
			fprintf(stderr,"theta_m[%ld] = (%.5e,%.5e), ", icolumn, creal(theta_m[icolumn]), cimag(theta_m[icolumn]));
			fprintf(stderr,"theta_m[%ld] = (%.5e,%.5e)\n", icolumn+1, creal(theta_m[icolumn+1]), cimag(theta_m[icolumn+1]));	exit(EXIT_FAILURE);
		}
	}
	for(long icolumn = 0 ; icolumn < num_real ; icolumn++){
		if( fabs(cimag(theta_m[icolumn])) > DBL_EPSILON || fabs(creal(theta_m[icolumn+N]) + creal(theta_m[icolumn])) > FLT_EPSILON){
			fprintf(stderr,"theta_m[%ld] = (%.5e,%.5e), ", icolumn, creal(theta_m[icolumn]), cimag(theta_m[icolumn]));
			fprintf(stderr,"theta_m[%ld] = (%.5e,%.5e)\n", icolumn+N, creal(theta_m[icolumn+N]), cimag(theta_m[icolumn+N]));	exit(EXIT_FAILURE);
		}
	}
	for(long icolumn = num_real ; icolumn < N-1 ; icolumn++){
		if( cimag(theta_m[icolumn+1]) + FLT_EPSILON < cimag(theta_m[icolumn]) ){
			fprintf(stderr,"theta_m[%ld] = (%.5e,%.5e), ", icolumn, creal(theta_m[icolumn]), cimag(theta_m[icolumn]));
			fprintf(stderr,"theta_m[%ld] = (%.5e,%.5e)\n", icolumn+1, creal(theta_m[icolumn+1]), cimag(theta_m[icolumn+1]));	exit(EXIT_FAILURE);
		}
	}
	for(long icolumn = num_real ; icolumn < N ; icolumn++){
		if( cimag(theta_m[icolumn]) < DBL_EPSILON || fabs(theta_m[icolumn]-conj(theta_m[icolumn+N])) > FLT_EPSILON){
			fprintf(stderr,"theta_m[%ld] = (%.5e,%.5e), ", icolumn, creal(theta_m[icolumn]), cimag(theta_m[icolumn]));
			fprintf(stderr,"theta_m[%ld] = (%.5e,%.5e)\n", icolumn+N, creal(theta_m[icolumn+N]), cimag(theta_m[icolumn+N]));	exit(EXIT_FAILURE);
		}
	}
	//-------     End of check theta_m     -------
}
void sort_cpair(const long N, const long i, const long j, double complex theta_m[2*N], double complex Xi_m[2*N*2*N])
{
	double complex a = theta_m[i];
	theta_m[i] = theta_m[j];
	theta_m[j] = a;
	for(long irow = 0 ; irow < 2*N ; irow++){
		a = Xi_m[irow + 2*N*i];
		Xi_m[irow + 2*N*i] = Xi_m[irow + 2*N*j];
		Xi_m[irow + 2*N*j] = a;
	}
}
//======================================================================
//  makeXi                             Last updated on Nov 06, 2018.  
//======================================================================
void makeXi(const long N, const double complex eiSin[N], double complex eigenH[N*N], double complex eigenE[N*N],
	double complex theta_m[2*N], double complex Xi_m[2*N*2*N])
{
	for(long j = 0 ; j < N ; j++){
		theta_m[j] = 2.*asin(0.5*eiSin[j]);//revised on 20181106
		theta_m[j+N] = -theta_m[j];
	}
	for(long icolumn = 0 ; icolumn < N ; icolumn++){
		for(long irow = 0 ; irow < N ; irow++){
			Xi_m[irow + 2*N*icolumn] = eigenH[irow + N*icolumn];
			Xi_m[irow+N + 2*N*icolumn] = eigenE[irow + N*icolumn];
		}
	}
	for(long icolumn = 0 ; icolumn < N ; icolumn++){
		for(long irow = 0 ; irow < N ; irow++){
			Xi_m[irow + 2*N*(icolumn+N)] = exp(I*0.5*theta_m[icolumn])*Xi_m[irow + 2*N*icolumn];
			Xi_m[irow + 2*N*icolumn] = exp(-I*0.5*theta_m[icolumn])*Xi_m[irow + 2*N*icolumn];
			
			Xi_m[irow+N + 2*N*(icolumn+N)] = -Xi_m[irow+N + 2*N*icolumn];
		}
	}
}
//======================================================================
//  checkM                             Last updated on Oct 04, 2018.  
//======================================================================
void checkM(const long N, const double complex maa[N*N], const double complex mbb[N*N])
{
	for(long j = 0 ; j < N*N ; j++){
		if(cimag(maa[j]) != 0. || cimag(mbb[j]) != 0.){	fprintf(stderr,"imaginary error of matrices in checkM!\n");	exit(EXIT_FAILURE);}
	}
	for(long icolumn = 0 ; icolumn < N ; icolumn++){
		for(long irow = 0 ; irow < N ; irow++){
			if(creal(maa[irow+N*icolumn]) - creal(maa[icolumn+N*irow]) != 0.){	fprintf(stderr,"transpose error of maa in checkM!\n");	exit(EXIT_FAILURE);}
			if(creal(mbb[irow+N*icolumn]) - creal(mbb[icolumn+N*irow]) != 0.){	fprintf(stderr,"transpose error of mbb in checkM!\n");	exit(EXIT_FAILURE);}
		}
	}
}
//======================================================================
//  norm_calc                            Last updated on Jun 19, 2019.  
//======================================================================
double norm_calc(const long N, const double complex A[N*N], const double complex E[N], const double complex V1[N*N], const double complex V2[N*N], double complex CWork[])
{
	double norm[3] = {0.e0,0.e0,0.e0};
//norm[0]: ||A*V1*diag(E^-1) - V2||_F ,
//norm[1]: ||A*V1*diag(E^-1)||_F ,
//norm[2]: ||V2||_F .
	{// see " http://azalea.s35.xrea.com/blas/index.html "
		int iN = N;
		double complex alpha = 1.;
		double complex beta = 0.;
		double complex *C = &CWork[0];//calloc(iN*iN,sizeof(double complex));
		zgemm_("N", "N", &iN, &iN, &iN, &alpha, &A[0], &iN, &V1[0], &iN, &beta, &C[0],&iN);// "H" is an illegal value!
		for(int icolumn = 0 ; icolumn < iN ; icolumn++){
			for(int irow = 0 ; irow < iN ; irow++){
//				C[irow+iN*icolumn] /= E[icolumn];
				C[irow+iN*icolumn] *= conj(E[icolumn])/(creal(E[icolumn])*creal(E[icolumn])+cimag(E[icolumn])*cimag(E[icolumn]));// revised on 2019.06.19
			}
		}
		int incx = 1;
		int iN2 = iN*iN;
		norm[1] = dznrm2_(&iN2,&C[0],&incx);
		norm[2] = dznrm2_(&iN2,&V1[0],&incx);
		for(int i = 0 ; i < iN2 ; i++){C[i] -= V2[i];}
		norm[0] = dznrm2_(&iN2,&C[0],&incx);
		//SAFEFREE(C);
	}// http://www.netlib.org/lapack/explore-html/d7/d76/zgemm_8f.html
	/*{
		for(long i0 = 0 ; i0 < N ; i0++) {//i0 column
			for(long j0 = 0 ; j0 < N ; j0++) {//j0 row
				double complex z[3] = {0.e0+I*0.e0, 0.e0+I*0.e0, 0.e0+I*0.e0};
				for(long k0 = 0 ; k0 < N ; k0++) {z[1] += A[j0+N*k0]*V1[k0+N*i0];}//A*V1
				z[1] /= E[i0];//A*V1*diag(E^-1)
				z[2] = V2[j0+N*i0];//V2
				z[0] = z[1] - z[2];//A*V1*diag(E^-1) - V2
				for(long l0 = 0 ; l0 < 3 ; l0++) {
					norm[l0] += creal(z[l0])*creal(z[l0])+cimag(z[l0])*cimag(z[l0]);
				}
			}
		}
		for(int i0 = 0 ; i0 < 3 ; i0++) {norm[i0] = sqrt(norm[i0]);}
	}*/
	return(norm[0]/(norm[1]+norm[2]));
}
//======================================================================
//  makeE                                Last updated on Jun 19, 2019.  
//======================================================================
void makeE(const int L[2], const double complex maa[2*L[0]*L[1]*2*L[0]*L[1]], const double complex eiSin[2*L[0]*L[1]],
 const double complex eigenH[2*L[0]*L[1]*2*L[0]*L[1]], double complex eigenE[2*L[0]*L[1]*2*L[0]*L[1]])
{
	{// see " http://azalea.s35.xrea.com/blas/index.html "
		int iN = 2*L[0]*L[1];
		double complex alpha = 1.;
		double complex beta = 0.;
		zgemm_("N", "N", &iN, &iN, &iN, &alpha, &maa[0], &iN, &eigenH[0], &iN, &beta, &eigenE[0],&iN);// "H" is an illegal value!
		for(int icolumn = 0 ; icolumn < iN ; icolumn++){
			for(int irow = 0 ; irow < iN ; irow++){
//				eigenE[irow + iN*icolumn] /= eiSin[icolumn];
				eigenE[irow + iN*icolumn] *= conj(eiSin[icolumn])/(creal(eiSin[icolumn])*creal(eiSin[icolumn])+cimag(eiSin[icolumn])*cimag(eiSin[icolumn]));// revised on 2019.06.19
			}
		}
	}// http://www.netlib.org/lapack/explore-html/d7/d76/zgemm_8f.html
	/*for(long icolumn = 0 ; icolumn < 2*L[0]*L[1] ; icolumn++){
		for(long irow = 0 ; irow < 2*L[0]*L[1] ; irow++){
			eigenE[irow + 2*L[0]*L[1]*icolumn] = 0.+ I*0.;
			for(long isum = 0 ; isum < 2*L[0]*L[1] ; isum++){
				eigenE[irow + 2*L[0]*L[1]*icolumn] += maa[irow + 2*L[0]*L[1]*isum]*eigenH[isum + 2*L[0]*L[1]*icolumn];
			}
			eigenE[irow + 2*L[0]*L[1]*icolumn] /= eiSin[icolumn];
		}
	}*/
}
//======================================================================
//  sortSin                              Last updated on Nov 12, 2018.  
//======================================================================
long sortSin(const int L[2], double complex eiSin[2*L[0]*L[1]], double complex eigenH[2*L[0]*L[1]*2*L[0]*L[1]])
{
	long count_real = 0;
	for(long j = 0 ; j < 2*L[0]*L[1] ; j++){//revised on 20181112
		if(pow(cimag(eiSin[j]),2)/pow(creal(eiSin[j]),2) <= DBL_EPSILON){
			eiSin[j] = creal(eiSin[j]) + I*0.;
			count_real++;
		}else if(pow(creal(eiSin[j]),2)/pow(cimag(eiSin[j]),2) <= DBL_EPSILON){
			eiSin[j] = 0. + I*cimag(eiSin[j]);
		}
	}
	for(long target_j = 0 ; target_j < count_real ; target_j++){
		if(cimag(eiSin[target_j]) != 0.){
			for(long j = target_j + 1 ; j < 2*L[0]*L[1] ; j++){
				if(cimag(eiSin[j]) == 0.){
					double dummy = creal(eiSin[j]);
					eiSin[j] = eiSin[target_j];
					eiSin[target_j] = dummy;
					for(long i = 0 ; i < 2*L[0]*L[1] ; i++){
						double complex cdummy = eigenH[i + 2*L[0]*L[1]*j];
						eigenH[i + 2*L[0]*L[1]*j] = eigenH[i + 2*L[0]*L[1]*target_j];
						eigenH[i + 2*L[0]*L[1]*target_j] = cdummy;
					}
				}
			}
		}
	}
//------------------- real-value sorting ------------------------
	for(long target_j = 0 ; target_j < count_real - 1 ; target_j++){
		for(long j = target_j + 1 ; j < count_real ; j++){
			if(creal(eiSin[j]) > creal(eiSin[target_j])){
				double dummy = creal(eiSin[j]);
				eiSin[j] = eiSin[target_j];
				eiSin[target_j] = dummy;
				for(long i = 0 ; i < 2*L[0]*L[1] ; i++){
					double complex cdummy = eigenH[i + 2*L[0]*L[1]*j];
					eigenH[i + 2*L[0]*L[1]*j] = eigenH[i + 2*L[0]*L[1]*target_j];
					eigenH[i + 2*L[0]*L[1]*target_j] = cdummy;
				}
			}
		}
	}
//------------------- complex-value sorting -----------------------
	for(long target_j = count_real ; target_j < 2*L[0]*L[1] - 1 ; target_j++){
		for(long j = target_j + 1 ; j < 2*L[0]*L[1] ; j++){
			if(cimag(eiSin[j]) < cimag(eiSin[target_j])){
				double complex dummy = eiSin[j];
				eiSin[j] = eiSin[target_j];
				eiSin[target_j] = dummy;
				for(long i = 0 ; i < 2*L[0]*L[1] ; i++){
					double complex cdummy = eigenH[i + 2*L[0]*L[1]*j];
					eigenH[i + 2*L[0]*L[1]*j] = eigenH[i + 2*L[0]*L[1]*target_j];
					eigenH[i + 2*L[0]*L[1]*target_j] = cdummy;
				}
			}
		}
	}
	return(count_real);
}
//======================================================================
//  makeSin                              Last updated on Aug. 10, 2021.  
//======================================================================
void makeSin(const int L[2], const double complex eiSin2[2*L[0]*L[1]], double complex eiSin[2*L[0]*L[1]])
{
	for(long j = 0 ; j < 2*L[0]*L[1] ; j++){
		double complex z = sqrt(eiSin2[j]);
		if(creal(z) < 0.){	fprintf(stderr,"sqrt error in Ecalc!\n");	exit(EXIT_FAILURE);}
		if(cimag(z) < -FLT_EPSILON*fabs(creal(z))){z = -z;}//revised on 20210810	Old: if(cimag(z) < 0.){z = -z;}
		eiSin[j] = z;
	}
}
//======================================================================
//  realizeH                             Last updated on Oct 03, 2018.  
//======================================================================
void realizeH(const int L[2], const double complex eiSin2[2*L[0]*L[1]], double complex eigenH[2*L[0]*L[1]*2*L[0]*L[1]])
{
	for(long j = 0 ; j < 2*L[0]*L[1] ; j++){
		if(cimag(eiSin2[j]) == 0.){
			double max_abs = fabs(eigenH[0 + 2*L[0]*L[1]*j]);
			long i_abs = 0;
			for(long i = 1 ; i < 2*L[0]*L[1] ; i++){
				double abs0 = fabs(eigenH[i + 2*L[0]*L[1]*j]);
				if(abs0 > max_abs){i_abs = i;	max_abs = abs0;}
			}
			if(creal(eigenH[i_abs + 2*L[0]*L[1]*j]) == 0.){
				for(long i = 0 ; i < 2*L[0]*L[1] ; i++){eigenH[i + 2*L[0]*L[1]*j] = cimag(eigenH[i + 2*L[0]*L[1]*j]);}
			}else if(creal(eigenH[i_abs + 2*L[0]*L[1]*j]) != 0. && cimag(eigenH[i_abs + 2*L[0]*L[1]*j]) != 0.){
				double complex norm_H = max_abs/eigenH[i_abs + 2*L[0]*L[1]*j];
				for(long i = 0 ; i < 2*L[0]*L[1] ; i++){eigenH[i + 2*L[0]*L[1]*j] = norm_H*eigenH[i + 2*L[0]*L[1]*j];}
			}
		}
	}
}
//======================================================================
//  ReigenVcalc                           Last updated on Jan 14, 2019  
//======================================================================
//double norm_calc(const long N, const double complex A[N*N], const double complex E[N], const double complex V1[N*N], const double complex V2[N*N]);
long ldgeev_(const long N, long double A[N*N], long double E_real[N], 
            long double E_imag[N], long double V_real[N*N]);
extern void dgeev_(char[], char[], long *N1, double A[], long *N2, double ER[],
                   double EI[], double vdummy[], long *ndummy, double VR[], 
                   long *N3, double work[], long *N4, long *info);
void ReigenVcalc(const char InPrecision[BUFSIZE], const int L[2], const double complex mbbmaa[2*L[0]*L[1]*2*L[0]*L[1]], double complex eiSin2[2*L[0]*L[1]], double complex eigenH[2*L[0]*L[1]*2*L[0]*L[1]])
{
	long NA = 2*L[0]*L[1];
    double *ER;
    ER = malloc(sizeof(double)*NA);
    double *EI;
    EI = malloc(sizeof(double)*NA);
    double *VR;
    VR = malloc(sizeof(double)*NA*NA);
	long info = 0;
	if(strncmp(InPrecision, "D", 1) == 0){
		double *RM;
		RM = malloc(sizeof(double)*NA*NA);
		for(long i = 0 ; i < NA ; i++) {
			for(long j = 0 ; j < NA ; j++) {
				RM[j+i*NA] = creal(mbbmaa[j+i*NA]);
			}
		} 
		double vdummy[1*1];
		long ndummy = 1;
		long n3 = 4 * NA;
		double *work;
		work = malloc(sizeof(double)*n3);
		dgeev_("N", "V", &NA, RM, &NA, ER, EI, vdummy, &ndummy, VR, &NA, work, &n3, 
			&info);
		if(info != 0){fprintf(stderr,"info = %ld in dgeev of ReigenVcalc!\n", info);}
		SAFEFREE(RM);SAFEFREE( work );
	}else if(strncmp(InPrecision, "L", 1) == 0){
		long double *RML;
		RML = malloc(sizeof(long double)*NA*NA);
		for(long i = 0 ; i < NA ; i++) {
			for(long j = 0 ; j < NA ; j++) {
				RML[j+i*NA] = creal(mbbmaa[j+i*NA]);
			}
		} 
		long double *ERL;
		ERL = malloc(sizeof(long double)*NA);
		long double *EIL;
		EIL = malloc(sizeof(long double)*NA);
		long double *VRL;
		VRL = malloc(sizeof(long double)*NA*NA);
		info = ldgeev_(NA, RML, ERL, EIL, VRL);
		SAFEFREE(RML);
		for(long i = 0 ; i < NA ; i++) {
			ER[i] = ERL[i];
			EI[i] = EIL[i];
			for(long j = 0 ; j < NA ; j++) {
				VR[j+i*NA] = VRL[j+i*NA];
			}
		}
		SAFEFREE(ERL);SAFEFREE(EIL);SAFEFREE(VRL);
	}else{	fprintf(stderr,"'%s' is wrong in eigenVcalc!\n", InPrecision);	exit(EXIT_FAILURE);}
	double *VI;
	VI = calloc(NA*NA,sizeof(double));//VI = 0.e0
	for(long i = 0 ; i < NA*NA ; i++) {VI[i] = 0.;}
	for(long i = 0 ; i < NA-1 ; i++) {
		if(fabs(EI[i]) > 0.e0) { 
			for(long j = 0 ; j < NA ; j++) {
				VI[j+i*NA] = VR[j+(i+1)*NA];
				VI[j+(i+1)*NA] = -VI[j+i*NA];
				VR[j+(i+1)*NA] = VR[j+i*NA];
			}
			i++;
		}
	}
	for(long i = 0 ; i < NA ; i++) {eiSin2[i] = ER[i] + I*EI[i];}
	for(long i = 0 ; i < NA*NA ; i++) {eigenH[i] = VR[i] + I*VI[i];}
	SAFEFREE(ER);SAFEFREE(EI);SAFEFREE(VR);SAFEFREE(VI);
//	fprintf(stderr,"||MM*H/(k^2)-H||_F/(||MM*H/(k^2)||_F+||H||_F)= %.5e in ReigenVcalc\n", norm_calc(NA, mbbmaa, eiSin2, eigenH, eigenH));
}
//======================================================================
//  eigenVcalc                            Last updated on Jan 14, 2019  
//======================================================================
//double norm_calc(const long N, const double complex A[N*N], const double complex E[N], const double complex V1[N*N], const double complex V2[N*N]);
long lzgeev_(const long N, long double complex A[N*N], long double complex E[N],
             long double complex V[N*N]);
extern void zgeev_(char[], char[], long *N1, double complex A[], long *N2, 
                   double complex E[], double complex vdummy[], long *ndummy, 
                   double complex V[], long *N3, double complex work[], 
                   long *N4, double Dwork[], long *info);
void eigenVcalc(const char InPrecision[BUFSIZE], const int L[2], const double complex mbbmaa[2*L[0]*L[1]*2*L[0]*L[1]], double complex eiSin2[2*L[0]*L[1]], double complex eigenH[2*L[0]*L[1]*2*L[0]*L[1]])
{
	long NA = 2*L[0]*L[1];
	long info = 0;
	if(strncmp(InPrecision, "D", 1) == 0){
		double complex *ZM;
		ZM = malloc(sizeof(double complex)*NA*NA);
		for(long i = 0 ; i < NA ; i++){
			for(long j = 0 ; j < NA ; j++){ZM[j+i*NA] = mbbmaa[j+i*NA];}
		} 
		double complex vdummy[1*1];
		long ndummy = 1;
		long n3 = 2 * NA;
		double complex *work;
		work = malloc(n3*sizeof(double complex));
		double *Dwork;
		Dwork = malloc(n3*sizeof(double));
		zgeev_("N", "V", &NA, ZM, &NA, eiSin2, vdummy, &ndummy, eigenH, &NA, work, &n3, Dwork, &info);
		if(info != 0){fprintf(stderr,"info = %ld in zgeev of eigenVcalc!\n", info);}
		SAFEFREE(ZM);SAFEFREE(work);SAFEFREE(Dwork);
	}else if(strncmp(InPrecision, "L", 1) == 0){
		long double complex *ZML;
		ZML = malloc(sizeof(long double complex)*NA*NA);
		for(long i = 0 ; i < NA ; i++){
			for(long j = 0 ; j < NA ; j++){ZML[j+i*NA] = mbbmaa[j+i*NA];}
		} 
		long double complex *EL;
		EL = malloc(sizeof(long double complex)*NA);
		long double complex *VL;
		VL = malloc(sizeof(long double complex)*NA*NA);
		info = lzgeev_(NA, ZML, EL, VL);
		for(long j = 0 ; j < NA ; j++){eiSin2[j] = EL[j];}
		for(long i = 0 ; i < NA ; i++){
			for(long j = 0 ; j < NA ; j++){eigenH[j+i*NA] = VL[j+i*NA];}
		} 
		SAFEFREE(ZML);SAFEFREE(EL);SAFEFREE(VL);
	}else{	fprintf(stderr,"'%s' is wrong in eigenVcalc!\n", InPrecision);	exit(EXIT_FAILURE);}
	for(int j = 0 ; j < NA ; j++){//DBL_EPSILON
//		if(creal(eiSin2[j]) != 0. && fabs(creal(eiSin2[j])/cimag(eiSin2[j])) <= ERROR_MIN){eiSin2[j]=0.+I*cimag(eiSin2[j]);}
//		if(cimag(eiSin2[j]) != 0. && fabs(cimag(eiSin2[j])/creal(eiSin2[j])) <= ERROR_MIN){eiSin2[j]=creal(eiSin2[j])+I*0.;}
		if(pow(creal(eiSin2[j]),2)/pow(cimag(eiSin2[j]),2) <= DBL_EPSILON){eiSin2[j]=0.+I*cimag(eiSin2[j]);}//revised on 20181112
		if(pow(cimag(eiSin2[j]),2)/pow(creal(eiSin2[j]),2) <= DBL_EPSILON){eiSin2[j]=creal(eiSin2[j])+I*0.;}//revised on 20181112
	}
//	fprintf(stderr,"||MM*H/(k^2)-H||_F/(||MM*H/(k^2)||_F+||H||_F)= %.5e in eigenVcalc\n", norm_calc(NA, mbbmaa, eiSin2, eigenH, eigenH));//revised on 20181011
}
//======================================================================
//  set_M ---- set_D                     Last updated on Jan 10, 2019.  
//======================================================================
void set_D(const int BoundC[2], const int L[2], double Dup0[], double Dup1[], double Ddown0[], double Ddown1[]);
void set_M(const int BoundC[2], const int L[2], const double complex Vcomp[6*L[0]*L[1]], 
	double complex maa[2*L[0]*L[1]*2*L[0]*L[1]], double complex mbb[2*L[0]*L[1]*2*L[0]*L[1]], double complex CWork[], double DWork[])
{
	long NA1 = L[0]*L[1];
	long Nshift = NA1*NA1;
	double *Dall = &DWork[0];//calloc(Nshift*4,sizeof(double));
	double *Dup0 = &Dall[Nshift*0];
	double *Dup1 = &Dall[Nshift*1];
	double *Ddown0 = &Dall[Nshift*2];
	double *Ddown1 = &Dall[Nshift*3];
	set_D(BoundC, L, Dup0, Dup1, Ddown0, Ddown1);
	long NA2 = 2*NA1;
	for(long j = 0 ; j < NA2*NA2 ; j++){maa[j] = 0. + I*0.;	mbb[j] = 0. + I*0.;}
	for(long irow = 0 ; irow < NA1 ; irow++){
			maa[irow	+ NA2*irow		] = Vcomp[irow + NA1*0];
			maa[irow+NA1	+ NA2*(irow+NA1)	] = Vcomp[irow + NA1*1];
			mbb[irow	+ NA2*irow		] = Vcomp[irow + NA1*3];
			mbb[irow+NA1	+ NA2*(irow+NA1)	] = Vcomp[irow + NA1*4];
	}
	{
		double complex *Call = &CWork[0];//calloc(Nshift*(4+3),sizeof(double complex)); Nshift*(4+3) <= NA2*NA2*12
		double complex *sum0 = &Call[Nshift*0];
		double complex *sum1 = &Call[Nshift*1];
		double complex *sum2 = &Call[Nshift*2];
		double complex *sum3 = &Call[Nshift*3];
		{
			double complex *D1V = &Call[Nshift*4];
			double complex *D0V = &Call[Nshift*5];
			for(long isum = 0 ; isum < NA1 ; isum++){
				for(long irow = 0 ; irow < NA1 ; irow++){
					D1V[irow + NA1*isum] =  Ddown1[irow + NA1*isum]*Vcomp[isum + NA1*2];
					D0V[irow + NA1*isum] =  Ddown0[irow + NA1*isum]*Vcomp[isum + NA1*2];
				}
			}
			double complex *cD = &Call[Nshift*6];
			for(long j = 0 ; j < Nshift ; j++){cD[j] = Dup1[j];}
			int iN = NA1;
			double complex alpha = 1.;
			double complex beta = 0.;
			zgemm_("N", "N", &iN, &iN, &iN, &alpha, &D1V[0], &iN, &cD[0], &iN, &beta, &sum0[0],&iN);// "H" is an illegal value!
			zgemm_("N", "N", &iN, &iN, &iN, &alpha, &D0V[0], &iN, &cD[0], &iN, &beta, &sum2[0],&iN);// "H" is an illegal value!
			for(long j = 0 ; j < Nshift ; j++){cD[j] = Dup0[j];}
			zgemm_("N", "N", &iN, &iN, &iN, &alpha, &D1V[0], &iN, &cD[0], &iN, &beta, &sum1[0],&iN);// "H" is an illegal value!
			zgemm_("N", "N", &iN, &iN, &iN, &alpha, &D0V[0], &iN, &cD[0], &iN, &beta, &sum3[0],&iN);// "H" is an illegal value!
			for(long icolumn = 0 ; icolumn < NA1; icolumn++){
				for(long irow = 0 ; irow < NA1 ; irow++){
					maa[irow	+ NA2*icolumn		] +=  sum0[irow + NA1*icolumn];
					maa[irow	+ NA2*(icolumn+NA1)	]  = -sum1[irow + NA1*icolumn];
					maa[irow+NA1	+ NA2*icolumn		]  = -sum2[irow + NA1*icolumn];
					maa[irow+NA1	+ NA2*(icolumn+NA1)	] +=  sum3[irow + NA1*icolumn];
				}
			}
			for(long isum = 0 ; isum < NA1 ; isum++){
				for(long irow = 0 ; irow < NA1 ; irow++){
					D0V[irow + NA1*isum] =  Dup0[irow + NA1*isum]*Vcomp[isum + NA1*5];
					D1V[irow + NA1*isum] =  Dup1[irow + NA1*isum]*Vcomp[isum + NA1*5];
				}
			}
			for(long j = 0 ; j < Nshift ; j++){cD[j] = Ddown0[j];}
			zgemm_("N", "N", &iN, &iN, &iN, &alpha, &D0V[0], &iN, &cD[0], &iN, &beta, &sum0[0],&iN);// "H" is an illegal value!
			zgemm_("N", "N", &iN, &iN, &iN, &alpha, &D1V[0], &iN, &cD[0], &iN, &beta, &sum2[0],&iN);// "H" is an illegal value!
			for(long j = 0 ; j < Nshift ; j++){cD[j] = Ddown1[j];}
			zgemm_("N", "N", &iN, &iN, &iN, &alpha, &D0V[0], &iN, &cD[0], &iN, &beta, &sum1[0],&iN);// "H" is an illegal value!
			zgemm_("N", "N", &iN, &iN, &iN, &alpha, &D1V[0], &iN, &cD[0], &iN, &beta, &sum3[0],&iN);// "H" is an illegal value!
		}
		for(long icolumn = 0 ; icolumn < NA1; icolumn++){
			for(long irow = 0 ; irow < NA1 ; irow++){
				mbb[irow	+ NA2*icolumn		] += sum0[irow + NA1*icolumn];
				mbb[irow	+ NA2*(icolumn+NA1)	]  = sum1[irow + NA1*icolumn];
				mbb[irow+NA1	+ NA2*icolumn		]  = sum2[irow + NA1*icolumn];
				mbb[irow+NA1	+ NA2*(icolumn+NA1)	] += sum3[irow + NA1*icolumn];
			}
		}
		//SAFEFREE(Call);
	}
	//SAFEFREE(Dall);
}
void set_Msub(const int L[2], double complex mbbmaa[2*L[0]*L[1]*2*L[0]*L[1]], double complex maa[2*L[0]*L[1]*2*L[0]*L[1]], double complex mbb[2*L[0]*L[1]*2*L[0]*L[1]])
{
// see " http://azalea.s35.xrea.com/blas/index.html "
	int iN = 2*L[0]*L[1];
	double complex alpha = 1.;
	double complex beta = 0.;
	zgemm_("N", "N", &iN, &iN, &iN, &alpha, &mbb[0], &iN, &maa[0], &iN, &beta, &mbbmaa[0],&iN);
// http://www.netlib.org/lapack/explore-html/d7/d76/zgemm_8f.html
	/*for(long icolumn = 0 ; icolumn < NA2 ; icolumn++){
		for(long irow = 0 ; irow < NA2 ; irow++){
			double complex sum = 0. + I*0.;
			for(long isum = 0 ; isum < NA2 ; isum++){
				sum += mbb[irow + NA2*isum]*maa[isum + NA2*icolumn];
			}
			mbbmaa[irow + NA2*icolumn] = sum;
		}
	}*/
}
//======================================================================
//  set_D                                Last updated on Oct 06, 2018.  
//======================================================================
void set_D(const int BoundC[2], const int L[2], double Dup0[], double Dup1[], double Ddown0[], double Ddown1[])
{
	for(long j = 0 ; j < L[0]*L[1]*L[0]*L[1] ; j++){	Dup0[j] = 0. + I*0.;	Dup1[j] = 0. + I*0.;}
	long sw[2];
	for(long i1 = 0 ; i1 < L[1] ; i1++){
		sw[1] = (i1 + 1)/L[1];
		long i1a = (i1 + 1)%L[1];
		for(long i0 = 0 ; i0 < L[0] ; i0++){
			sw[0] = (i0 + 1)/L[0];
			long i0a = (i0 + 1)%L[0];
			long ipoint = i0 + L[0]*i1;
			long ipoint_0a = i0a +  L[0]*i1;
			long ipoint_1a = i0  + L[0]*i1a;
			Dup0[ipoint + L[0]*L[1]*ipoint] = -1.;
			Dup0[ipoint + L[0]*L[1]*ipoint_0a] = (double) ((1-sw[0]) + sw[0]*BoundC[0]);
			Dup1[ipoint + L[0]*L[1]*ipoint] = -1.;
			Dup1[ipoint + L[0]*L[1]*ipoint_1a] = (double) ((1-sw[1]) + sw[1]*BoundC[1]);
/*			if(sw[0] != 0 && sw[0] != 1){	fprintf(stderr,"sw[0] error in set_D!\n");	exit(EXIT_FAILURE);}
			if(sw[1] != 0 && sw[1] != 1){	fprintf(stderr,"sw[1] error in set_D!\n");	exit(EXIT_FAILURE);}*/
		}
	}
/*	for(long j = 0 ; j < L[0]*L[1]*L[0]*L[1] ; j++){
		if(Dup0[j] != 0. && Dup0[j] != 1. && Dup0[j] != -1.){	fprintf(stderr,"Dup0 error in set_D!\n");	exit(EXIT_FAILURE);}
		if(Dup1[j] != 0. && Dup1[j] != 1. && Dup1[j] != -1.){	fprintf(stderr,"Dup1 error in set_D!\n");	exit(EXIT_FAILURE);}
	}*/
	for(long icolumn = 0 ; icolumn < L[0]*L[1] ; icolumn++){
		for(long irow = 0 ; irow < L[0]*L[1] ; irow++){
			Ddown0[irow + L[0]*L[1]*icolumn] = -Dup0[icolumn + L[0]*L[1]*irow];
			Ddown1[irow + L[0]*L[1]*icolumn] = -Dup1[icolumn + L[0]*L[1]*irow];
		}
	}
}
//======================================================================
//  set_V                                Last updated on Oct 06, 2018.  
//======================================================================
long set_V(FILE *fp_i, const int L[2], const double omega, double complex Vcomp[6*L[0]*L[1]])
{
	char buf[BUFSIZE];
	for(int i1 = 0 ; i1 < L[1] ; i1++){
		for(int i0 = 0 ; i0 < L[0] ; i0++){
			for(int s0 = 0 ; s0 < 2 ; s0++){
				if(fgets(buf, sizeof( buf ), fp_i) == NULL) {	fprintf(stderr,"Lines can not read in data-file!\n");	exit(EXIT_FAILURE);}
				int dummy0, dummy1, dummys;
				double dr[3], di[3], df[3];
				if(sscanf(buf,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
					&dummy0, &dummy1, &dummys, 
					&dr[0], &di[0], &dr[1], &di[1], &dr[2], &di[2],
					&df[0], &df[1], &df[2]) == 3 + 6 + 3){
					if(s0 != dummys || i0 != dummy0 || i1 != dummy1) {	fprintf(stderr,"%d, %d or %d is not matched data in data-file!\n", i0, i1, s0);	exit(EXIT_FAILURE);}
					for(int j = 0 ; j < 3 ; j++){
						Vcomp[i0 +  L[0]*i1 + L[0]*L[1]*(j + 3*s0)] = (dr[j] + I*di[j])*df[j];
					}
				}else{	fprintf(stderr,"12 parameters cannot be read in data-file!\n");	exit(EXIT_FAILURE);}
			}
		}
	}
	
	for(long j = 0 ; j < 6*L[0]*L[1] ; j++){Vcomp[j] *= omega;}

/*	double max_muepsi[6] = {0.,0.,0.,0.,0.,0.};
	for(int s = 0 ; s < 6 ; s++){
	for(long j = 0 ; j < L[0]*L[1] ; j++){
		if(creal(Vcomp[j+L[0]*L[1]*s]) < 0.) {	fprintf(stderr,"creal(Vcomp[j])<0. in set_V!\n");	exit(EXIT_FAILURE);}
		if(fabs(Vcomp[j+L[0]*L[1]*s]) > max_muepsi[s]){max_muepsi[s]=fabs(Vcomp[j+L[0]*L[1]*s]);}
		}
		fprintf(stderr,"max_omegamuepsi [%d]= %.5e in set_V!\n", s, max_muepsi[s]);
	}*/

	for(long j = 0 ; j < L[0]*L[1] ; j++){
		Vcomp[j + L[0]*L[1]*2] = 1./Vcomp[j + L[0]*L[1]*2];
		Vcomp[j + L[0]*L[1]*5] = 1./Vcomp[j + L[0]*L[1]*5];
	}

/*	double max_Vcomp = 0.;
	for(long j = 0 ; j < 6*L[0]*L[1] ; j++){
		if(fabs(Vcomp[j]) > max_Vcomp){max_Vcomp=fabs(Vcomp[j]);}
	}
	fprintf(stderr,"max_Vcomp = %.5e in set_V!\n", max_Vcomp);*/

	long comp_count = 0;
	for(long j = 0 ; j < 6*L[0]*L[1] ; j++){if(cimag(Vcomp[j]) != 0.){comp_count++;}}
	return(comp_count);
}
//======================================================================
//  read_para                            Last updated on Oct 01, 2018.  
//======================================================================
void read_para(FILE *fp_i, const char data_name[BUFSIZE], char comment[BUFSIZE], double *du_2, double *k_0, int L[2], double *omega)
{
	double lambda;
	char buf[BUFSIZE];
	if(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
		rm_comma(buf);
		snprintf(comment,BUFSIZE,"# Data file          : ' %s ' %s", data_name, &buf[10]);// https://www.ipa.go.jp/security/awareness/vendor/programmingv1/b06_02.html
	}
	int j_count = -5;
	long double u2_0, u2_1;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
		rm_space(buf);
		if(strncmp(buf, "# wavelength", 12) == 0 && sscanf(buf,"%*[^=] %*[=] %lf", &lambda) == 1){
			j_count++;//-4
		}else if(strncmp(buf, "# k_0", 5) == 0 && sscanf(buf,"%*[^=] %*[=]  %Lf %*[;] %*[^;] %*[;]  %Lf %*[^=] %*[=] %lf", &u2_0, &u2_1, k_0) == 3){
			u2_1 = (u2_1 - u2_0);
			j_count++;//-3
		}else if(strncmp(buf, "# l <", 5) == 0 && sscanf(buf,"%*[^<] %*[<] %d %*[^<] %*[<] %d", &L[0], &L[1]) == 2) {
			j_count++;//-2
		}else if(strncmp(buf, "# l m", 5) == 0){j_count++;}//-1,0
		if(j_count == 0){goto start_data;}
	}
	fprintf(stderr,"4 parameters can not read in data-file!");	exit(EXIT_FAILURE);
	start_data:;
	*omega = Pi2/(lambda*(*k_0));
	*du_2 = ((double) u2_1);
}
