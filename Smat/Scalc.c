//  Start of Scalc.h
#include "header_macro.h"
#include "constant.h"
void Out_Smat(char **argv, const time_t timer_ini, const char *f_prefix, const char *p_target, 
	const double omega, const char *S_ctrl, const char *solv, const int BoundC[2], const int RealNum05[2], const double complex Smatrix[]);
void Smat_calc(const int iway, const int NA, const int Noffset, const int Nimag, const double complex F01[], const double complex B01[], double complex Smatrix[]);
void final_calc(const char solv[], const double omega, const int BoundC[2], const int iway, const char f_prefix[], const char p_target[], 
		const int L[4], const int iend, const int Noffset, const int Nimag, 
		double complex F01[], double complex B01[], double complex F2[], double complex B2[], double complex Vmab[], double complex CWork[]);
void scattering_calc(const char solv[], const double omega, const int BoundC[2], const int iway, const char f_prefix[], 
		const int L[4], const int istep, const int Noffset, 
		double complex F01[], double complex B01[], double complex F2[], double complex B2[], double complex Vmab[], double complex CWork[]);
void initial_calc(const char solv[], const double omega, const int BoundC[2], const int iway, const char f_prefix[], const char p_target[], 
		const int L[4], const int istart, const int Noffset, 
		double complex F01[], double complex B01[], double complex F2[], double complex B2[], double complex Vmab[], double complex CWork[]);
double get_RealNum05omega(const char f_prefix[], const char p_target[], const int L[4], int RealNum05[2]);
int num_u2(const int L[4], const char buf[]);
void input_file1(FILE *fp_i1, int L[4]);
void input_file0(FILE *fp_i0, int BoundC[], char f_prefix[], char p_target[], char S_ctrl[], char solv[]);
void main_calc(char **argv, const char data_name[BUFSIZE], const char out_name[BUFSIZE], const int BoundC[2], const char InPrecision[BUFSIZE], const char Xi_name[BUFSIZE]);
//  End of Scalc.h
//======================================================================
//   main ---- input_file0, input_file1, get_RealNum05,                 
//             initial_calc, scattering_calc, final_calc,               
//             Smat_calc, Out_Smat                                      
//======================================================================
int main(int argc, char **argv)
{
	time_t timer_ini = time(0);	fprintf(stderr,"# Start time of mode calculation = %s\n", ctime(&timer_ini));
	if(argc != 2) {fprintf(stderr,"error: number of files \n");	exit(EXIT_FAILURE);}
	else if(strncmp(argv[1], "-v", 2) == 0 || strcmp(argv[1], "--version") == 0 ) {
		fprintf(stderr,"The '%s' creates scattering matrix.\n", argv[0]);
		fprintf(stderr,"Version 19.09.23 is compiled at %s on %s.\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(stderr," Source code : '%s'\n Author      : Tatsuya Usuki\n URL         : http://www.smatran.org\n", __FILE__);
		fprintf(stderr," References  : 'Wave scattering in frequency domain' as 'Formulation.pdf' on May 20, 2018;\n");
		fprintf(stderr,"               './Smat_fig.pdf' on Nov 06, 2018.\n");
		fprintf(stderr,"There is NO warranty.\n");
		exit(EXIT_SUCCESS);//normal end
	}
//------- begin reading fundamental parameters -------
	char f_prefix[BUFSIZE], p_target[BUFSIZE*2], S_ctrl[BUFSIZE], solv[BUFSIZE];
	int BoundC[2];
	{
		FILE *fp_i0;
		fp_i0 = fopen(argv[1],"r");
		if (fp_i0 == NULL){	fprintf(stderr,"open error!: open input-file0!\n");	exit(EXIT_FAILURE);}
		input_file0(fp_i0, BoundC, f_prefix, p_target, S_ctrl, solv);
		if(fclose(fp_i0) != 0) {	fprintf(stderr,"fclose error after input_file0!\n");	exit(EXIT_FAILURE);}
	}
	int L[4];//L[0] = L0, L[1] = L1, L[2] = L2 + 2*outer, L[3] = outer
	{
		char data_name[BUFSIZE];
		snprintf(data_name,sizeof(data_name),"%s_Med%s.dat", f_prefix, p_target);// https://www.ipa.go.jp/security/awareness/vendor/programmingv1/b06_02.html
		FILE *fp_i1;
		fp_i1 = fopen(data_name,"r");
		if (fp_i1 == NULL){	fprintf(stderr,"open error!: open input-file1!\n");	exit(EXIT_FAILURE);}
		input_file1(fp_i1, L);// revised on 20190917
		if(fclose(fp_i1) != 0) {	fprintf(stderr,"fclose error after input_file1!\n");	exit(EXIT_FAILURE);}
	}
//-------  end reading fundamental parameters  -------
	int iB, iT;
	iB = num_u2(L, p_target);	iT = num_u2(L, &p_target[BUFSIZE]);
	if(iB == iT) {	iT = (L[2] - 1) - iB;	fprintf(stderr,"The 2nd target ' %s ' is the same as the 1st, and iT is changed to %d\n", &p_target[BUFSIZE], iT);}
	
	if(iB >= iT){fprintf(stderr,"Target_B = %d for %s is larger than or equal to Target_T = %d for %s\n", iB, p_target, iT, &p_target[BUFSIZE]);	exit(EXIT_FAILURE);}
	else{fprintf(stderr,"Target_B = %d for %s, Target_T = %d for %s\n", iB, p_target, iT, &p_target[BUFSIZE]);}
	
	int RealNum05[2];
	double omega = get_RealNum05omega(f_prefix, p_target, L, RealNum05);	fprintf(stderr,"Omega = k/k0 = %.20e\n", omega);
	double complex *Smatrix = calloc((RealNum05[0]+RealNum05[1])*(RealNum05[0]+RealNum05[1]),sizeof(double complex));
	for(int iway = 1 ; iway >= -(strncmp(S_ctrl, "Y", 1) == 0 || strncmp(S_ctrl, "y", 1) == 0 ) ; iway -= 2){
		int istart = iB*((1+iway)/2)+iT*((1-iway)/2);int iend = iT*((1+iway)/2)+iB*((1-iway)/2);
		if(iway > 0){fprintf(stderr,"Forward calc.: iway = %d\n", iway);}else{fprintf(stderr,"Backward calc.: iway = %d\n", iway);}
		int Noffset = RealNum05[(1-iway)/2];
		int NA = 2*L[0]*L[1];
		size_t n_FB = (Noffset+NA)*Noffset + (Noffset+NA)*NA + NA*Noffset + NA*NA + NA*(NA*2+3);//revised on 20190108
		size_t n_CWork = n_FB + NA*2*NA*2*3;//revised on 20190108
		double complex *FB = calloc(n_CWork,sizeof(double complex));
		double complex *F01 = &FB[0];//(Noffset+NA)*Noffset
		double complex *B01 = &FB[(Noffset+NA)*Noffset];//(Noffset+NA)*NA
		double complex *F2 = &FB[(Noffset+NA)*(Noffset+NA)];//NA*Noffset
		double complex *B2 = &FB[(Noffset+NA)*(Noffset+NA) + NA*Noffset];//NA*NA
		double complex *Vmab = &FB[(Noffset+2*NA)*(Noffset+NA)];//NA*(NA*2+3)
		double complex *CWork = &FB[n_FB];//NA*2*NA*2*3, CWork means Temporary Working memory.
		initial_calc(solv, omega, BoundC, iway, f_prefix, p_target, L, istart, Noffset, F01, B01, F2, B2, Vmab, CWork);
		for(int istep = istart + iway ; iway*istep < iway*iend ; istep += iway){
//			fprintf(stderr,"Temp message: istep = %d\n", istep);//revised on 20190107
			scattering_calc(solv, omega, BoundC, iway, f_prefix, L, istep, Noffset, F01, B01, F2, B2, Vmab, CWork);
		}
		int Nimag = RealNum05[(1+iway)/2];
		final_calc(solv, omega, BoundC, iway, f_prefix, p_target, L, iend, Noffset, Nimag, F01, B01, F2, B2, Vmab, CWork);
		Smat_calc(iway, NA, Noffset, Nimag, F01, B01, Smatrix);
		SAFEFREE(FB);
	}
	Out_Smat(argv, timer_ini, f_prefix, p_target, omega, S_ctrl, solv, BoundC, RealNum05, Smatrix);
	SAFEFREE(Smatrix);
}
//======================================================================
//  num_u2                              Last updated on Sep 23, 2019.  
//======================================================================
int num_u2(const int L[4], const char buf[])
{
	int iI = -1;
	if(strncmp(buf, "B", 1) == 0){sscanf(&buf[1],"%d", &iI);}
	else if(strncmp(buf, "T", 1) == 0){sscanf(&buf[1],"%d", &iI);iI += L[2] - L[3];}
	else{sscanf(buf,"%d", &iI);iI += L[3];}
	return(iI);
}
void check_Smat(const int NS, const double complex Smatrix[], double norm[2]);
//======================================================================
//  Out_Smat ---- check_Smat               Last updated on Nov 12, 2018
//======================================================================
void Out_Smat(char **argv, const time_t timer_ini, const char *f_prefix, const char *p_target, 
	const double omega, const char *S_ctrl, const char *solv, const int BoundC[2],
	const int RealNum05[2], const double complex Smatrix[(RealNum05[0]+RealNum05[1])*(RealNum05[0]+RealNum05[1])])
{
	double norm[2];
	check_Smat(RealNum05[0]+RealNum05[1], Smatrix, norm);
	char out_name[BUFSIZE];snprintf(out_name,sizeof(out_name),"%s_Smatrix.dat", f_prefix);
	FILE *fp_o;
	fp_o = fopen(out_name,"w");	if (fp_o == NULL){	fprintf(stderr,"open error!: open o-file!\n");	exit(EXIT_FAILURE);}
	time_t timer_f = time(0);
	fprintf(fp_o,"# This file          : ' %s ' on %s", out_name, ctime(&timer_f));
	fprintf(fp_o,"# Exec. file         : ' %s ' when '%s' was compiled at %s on %s by C-version %ld\n", argv[0], __FILE__, __TIME__, __DATE__, __STDC_VERSION__);
	fprintf(fp_o,"# Calc. time         : Exec. time = %1.0f s\n", difftime(timer_f,timer_ini));
	fprintf(fp_o,"# Parameter file     : ' %s ' for prefix, boundary conditions and internal precision\n", argv[1]);
	fprintf(fp_o,"# Omega              : k/k_0 = %.20e\n", omega);
	fprintf(fp_o,"# Boundary conditions: B0 = %d, B1 = %d\n", BoundC[0], BoundC[1]);
	fprintf(fp_o,"# Target layers      : Bottom = %s, Top = %s\n", p_target, &p_target[BUFSIZE]);
	fprintf(fp_o,"# Mode numbers       : Bottom = %d, Top = %d\n", RealNum05[0], RealNum05[1]);
	fprintf(fp_o,"# Calc. mode         : both ways = %s, 'zgesv' = %s\n", S_ctrl, solv);
	fprintf(fp_o,"# Smat. residual norm: unitarity = %.16e, transpose = %.16e\n", norm[0], norm[1]);
	fprintf(fp_o,"# m n real(S_mn) imag(S_mn) |S_mn|^2 for m,n < %d\n",RealNum05[0]+RealNum05[1]);
	for(int column = 0 ; column < RealNum05[0]+RealNum05[1] ; column++){
		for(int row = 0 ; row < RealNum05[0]+RealNum05[1] ; row++){
			fprintf(fp_o,"%d %d %.20e %.20e %.20e\n", row, column, creal(Smatrix[row + (RealNum05[0]+RealNum05[1])*column]), cimag(Smatrix[row + (RealNum05[0]+RealNum05[1])*column]), pow(fabs(Smatrix[row + (RealNum05[0]+RealNum05[1])*column]),2));
		}
	}
	if(fclose(fp_o) != 0) {	fprintf(stderr,"fclose error after output_data!\n");	exit(EXIT_FAILURE);}
}
//======================================================================
//  Smat_calc                              Last updated on Nov 05, 2018
//======================================================================
void Smat_calc(const int iway, const int NA, const int Noffset, const int Nimag, 
	const double complex F01[(Noffset+NA)*Noffset], const double complex B01[(Noffset+NA)*NA],
	double complex S[(Noffset+Nimag)*(Noffset + Nimag)])
{
	if(iway == 1){//Noffset = RealNum05[0], Nimag = RealNum05[1]
		for(int column = 0 ; column < Noffset + Nimag; column++){
			for(int row = 0 ; row < Noffset + Nimag ; row++){
				if(column < Noffset){
					S[row + (Noffset+Nimag)*column] = F01[row + (Noffset+NA)*column];
				}else{
					S[row + (Noffset+Nimag)*column] = B01[row + (Noffset+NA)*(column-Noffset)];
				}
			}
		}
	}else if(iway == -1){//Noffset = RealNum05[1], Nimag = RealNum05[0]
		for(int column = Nimag ; column < Noffset + Nimag; column++){
			for(int row = Nimag ; row < Noffset + Nimag ; row++){
				S[row + (Nimag+Noffset)*column] = F01[(row-Nimag) + (Noffset+NA)*(column-Nimag)];
			}
		}
		for(int column = Nimag ; column < Noffset + Nimag; column++){
			for(int row = 0 ; row < Nimag ; row++){
				S[row + (Nimag+Noffset)*column] = F01[(row+Noffset) + (Noffset+NA)*(column-Nimag)];
			}
		}
	}else{
		fprintf(stderr,"iway = %d error in Smat-calc!\n", iway);	exit(EXIT_FAILURE);
	}
}
long set_V(const char data_name[], const int L[2], const double omega, double complex Vcomp[6*L[0]*L[1]]);
void set_M(const int BoundC[2], const int L[2], const double complex Vcomp[6*L[0]*L[1]], 
	double complex maa[2*L[0]*L[1]*2*L[0]*L[1]], double complex mbb[2*L[0]*L[1]*2*L[0]*L[1]], double complex CWork[]);
void matdata_file(const char f_prefix[], const char add_name[], const int Ltot, const int Lout, const int istep, char data_name[]);
void inv_process(const char solv[], const int NA, const int Noffset, double complex F01[], double complex B01[], double complex F2[], double complex B2[], double complex CWork[]);
void set_iniMat(const char Xi_file[], const int iway, const int NA, const int Noffset, double complex F01[], double complex B01[], double complex F2[], double complex B2[], double complex CWork[]);
//======================================================================
//  initial_calc ---- set_iniMat, inv_process, matdata_file, set_V,set_M
//                                        Last updated on Sep 14, 2019  
//======================================================================
void initial_calc(const char solv[], const double omega, const int BoundC[2], const int iway, const char f_prefix[], const char p_target[], 
		const int L[4], const int istart, const int Noffset, 
		double complex F01[(Noffset+2*L[0]*L[1])*Noffset], double complex B01[(Noffset+2*L[0]*L[1])*2*L[0]*L[1]],
		double complex F2[2*L[0]*L[1]*Noffset], double complex B2[2*L[0]*L[1]*2*L[0]*L[1]],
		double complex Vmab[2*L[0]*L[1]*(2*L[0]*L[1]*2+3)], double complex CWork[2*L[0]*L[1]*2*2*L[0]*L[1]*2*3])
{
	char Xi_file[BUFSIZE];
	if(iway > 0 || strcmp(p_target, &p_target[BUFSIZE]) == 0 || strncmp(&p_target[BUFSIZE], "B", 1) == 0){snprintf(Xi_file,sizeof(Xi_file),"%s_bXi.dat", f_prefix);}
	else{snprintf(Xi_file,sizeof(Xi_file),"%s_tXi.dat", f_prefix);}
	int NA = 2*L[0]*L[1];	assert(Noffset <= NA);
	set_iniMat(Xi_file, iway, NA, Noffset, F01, B01, F2, B2, CWork);
	inv_process(solv, NA, Noffset, F01, B01, F2, B2, CWork);
	fprintf(stderr,"\n");
	//transfer to the next calc.
	for(int column = 0 ; column < Noffset ; column++){
		for(int row = 0 ; row < NA ; row++){F2[row + NA*column] = F01[row+Noffset + (Noffset+NA)*column];}
	}
	for(int column = 0 ; column < Noffset ; column++){
		for(int row = 0 ; row < NA ; row++){F01[row+Noffset + (Noffset+NA)*column] = 0. + I*0.;}
	}
	char data_name[BUFSIZE];
	matdata_file(f_prefix, "_Med", L[2], L[3], istart, data_name);// revised on 20190914
	fprintf(stderr,"Xi file: %s, Input file: %s, ", Xi_file, data_name);//fprintf(stderr,"Xi file: %s, Input file: %s, istart = %d, ", Xi_file, data_name, istart);
	double complex *V = &Vmab[NA*NA*2];//double complex Vmab[NA*(NA*2+3)], double complex *V = calloc(6*L[0]*L[1],sizeof(double complex));
	if(set_V(data_name, L, omega, V) != 0){	fprintf(stderr,"V is not real in initial_calc\n");	exit(EXIT_FAILURE);}
	double complex *maa=&Vmab[0];
	double complex *mbb=&Vmab[NA*NA*1];//double complex *maa = calloc(NA*NA*2,sizeof(double complex));
	set_M(BoundC, L, V, maa, mbb, CWork);//set_M(BoundC, L, Vcomp, maa, mbb);
//	SAFEFREE(V);
	for(int column = 0 ; column < NA ; column++){
		for(int row = 0 ; row < NA ; row++){B2[row + NA*column] = B01[row+Noffset + (Noffset+NA)*column] + (I*((double)iway))*maa[row + NA*column + NA*NA*((1+iway)/2)];}//revised on 20181106
	}
//	SAFEFREE(maa);
	for(int column = 0 ; column < NA ; column++){
		for(int row = 0 ; row < NA ; row++){B01[row+Noffset + (Noffset+NA)*column] = 0. + I*0.;}
		B01[column+Noffset + (Noffset+NA)*column] = 1. + I*0.;
	}
	inv_process(solv, NA, Noffset, F01, B01, F2, B2, CWork);
	fprintf(stderr,"\n");
}
//======================================================================
// scattering_calc ---- matdata_file, set_V, set_M, inv_process
//                                         Last updated on Jan 10, 2019 
//======================================================================
void scattering_calc(const char solv[], const double omega, const int BoundC[2], const int iway, const char f_prefix[], 
		const int L[4], const int istep, const int Noffset, 
		double complex F01[(Noffset+2*L[0]*L[1])*Noffset], double complex B01[(Noffset+2*L[0]*L[1])*2*L[0]*L[1]],
		double complex F2[2*L[0]*L[1]*Noffset], double complex B2[2*L[0]*L[1]*2*L[0]*L[1]],
		double complex Vmab[2*L[0]*L[1]*(2*L[0]*L[1]*2+3)], double complex CWork[2*L[0]*L[1]*2*2*L[0]*L[1]*2*3])
{
	int NA = 2*L[0]*L[1];	assert(Noffset <= NA);
	char data_name[BUFSIZE];
	matdata_file(f_prefix, "_Med", L[2], L[3], istep, data_name);
	fprintf(stderr,"%s, ",data_name);
	double complex *V = &Vmab[NA*NA*2];//double complex Vmab[NA*(NA*2+3)], double complex *V = calloc(6*L[0]*L[1],sizeof(double complex));
	if(set_V(data_name, L, omega, V) != 0){	fprintf(stderr,"V is not real in scattering_calc, ");}//	exit(EXIT_FAILURE);} // revised on 20190911
	double complex *maa=&Vmab[0];
	double complex *mbb=&Vmab[NA*NA*1];//double complex *maa = calloc(NA*NA*2,sizeof(double complex));
	set_M(BoundC, L, V, maa, mbb, CWork);
	for(int sub = 0 ; sub < 2 ; sub++){
		for(int column = 0 ; column < Noffset ; column++){
			for(int row = 0 ; row < NA ; row++){F2[row + NA*column] = F01[row+Noffset + (Noffset+NA)*column];}
		}
		for(int column = 0 ; column < Noffset ; column++){
			for(int row = 0 ; row < NA ; row++){F01[row+Noffset + (Noffset+NA)*column] = 0. + I*0.;}
		}
		int i_shift = (1-iway)/2;
		i_shift = (i_shift + sub)%2;
		for(int column = 0 ; column < NA ; column++){
			for(int row = 0 ; row < NA ; row++){B2[row + NA*column] = B01[row+Noffset + (Noffset+NA)*column] + (I*((double)iway))*maa[row + NA*column + NA*NA*i_shift];}//revised on 20181106
		}
		for(int column = 0 ; column < NA ; column++){
			for(int row = 0 ; row < NA ; row++){B01[row+Noffset + (Noffset+NA)*column] = 0. + I*0.;}
			B01[column+Noffset + (Noffset+NA)*column] = 1. + I*0.;
		}
		inv_process(solv, NA, Noffset, F01, B01, F2, B2, CWork);
		if(sub == 0){fprintf(stderr,", ");}
	}
	fprintf(stderr,"\n");
}
void set_finMat(const char Xi_file[], const int iway, const int NA, const int Noffset, const int Nimag, 
	double complex F01[], double complex B01[], double complex F2[], double complex B2[], double complex CWork[]);
//======================================================================
//  final_calc ---- set_finMat, inv_process, matdata_file, set_V, set_M
//                                        Last updated on Sep 11, 2019  
//======================================================================
void final_calc(const char solv[], const double omega, const int BoundC[2], const int iway, const char f_prefix[], const char p_target[], 
		const int L[4], const int iend, const int Noffset, const int Nimag, 
		double complex F01[(Noffset+2*L[0]*L[1])*Noffset], double complex B01[(Noffset+2*L[0]*L[1])*2*L[0]*L[1]],
		double complex F2[2*L[0]*L[1]*Noffset], double complex B2[2*L[0]*L[1]*2*L[0]*L[1]],
		double complex Vmab[2*L[0]*L[1]*(2*L[0]*L[1]*2+3)], double complex CWork[2*L[0]*L[1]*2*2*L[0]*L[1]*2*3])
{
	int NA = 2*L[0]*L[1];	assert(Noffset <= NA);
	for(int column = 0 ; column < Noffset ; column++){
		for(int row = 0 ; row < NA ; row++){F2[row + NA*column] = F01[row+Noffset + (Noffset+NA)*column];}
	}
	for(int column = 0 ; column < Noffset ; column++){
		for(int row = 0 ; row < NA ; row++){F01[row+Noffset + (Noffset+NA)*column] = 0. + I*0.;}
	}
	char data_name[BUFSIZE];
	matdata_file(f_prefix, "_Med", L[2], L[3], iend, data_name);
	double complex *V = &Vmab[NA*NA*2];//double complex Vmab[NA*(NA*2+3)], double complex *V = calloc(3*NA,sizeof(double complex));
	if(set_V(data_name, L, omega, V) != 0){	fprintf(stderr,"V is not real in final_calc\n");	exit(EXIT_FAILURE);} // error message revised on 20190911
	double complex *maa=&Vmab[0];
	double complex *mbb=&Vmab[NA*NA*1];//double complex *maa = calloc(NA*NA*2,sizeof(double complex));
	set_M(BoundC, L, V, maa, mbb, CWork);
	for(int column = 0 ; column < NA ; column++){
		for(int row = 0 ; row < NA ; row++){B2[row + NA*column] = B01[row+Noffset + (Noffset+NA)*column] + (I*((double)iway))*maa[row + NA*column + NA*NA*((1-iway)/2)];}//revised on 20181106
	}
	for(int column = 0 ; column < NA ; column++){
		for(int row = 0 ; row < NA ; row++){B01[row+Noffset + (Noffset+NA)*column] = 0. + I*0.;}
		B01[column+Noffset + (Noffset+NA)*column] = 1. + I*0.;
	}
	inv_process(solv, NA, Noffset, F01, B01, F2, B2, CWork);
	fprintf(stderr,"\n");
	//transfer to the final calc.
	char Xi_file[BUFSIZE];
	if(iway < 0 || strcmp(p_target, &p_target[BUFSIZE]) == 0 || strncmp(&p_target[BUFSIZE], "B", 1) == 0){snprintf(Xi_file,sizeof(Xi_file),"%s_bXi.dat", f_prefix);}
	else{snprintf(Xi_file,sizeof(Xi_file),"%s_tXi.dat", f_prefix);}
	fprintf(stderr,"Input file: %s, Xi file: %s\n", data_name, Xi_file);
	set_finMat(Xi_file, iway, NA, Noffset, Nimag, F01, B01, F2, B2, CWork);
	inv_process(solv, NA, Noffset, F01, B01, F2, B2, CWork);
	fprintf(stderr,"\n");
}
//======================================================================
//  inv_process ---- zgesvx_, zgemm_     Last updated on Jun 19, 2019   
//======================================================================
extern void zgemm_(char *transa, char *transb, 
	const int *m, const int *n, const int *k,
	const double complex *alpha, const double complex *A, const int *ldA, 
	const double complex *B, const int *ldB,
	const double complex *beta , double complex *C, const int *ldC);
extern void zgesv_(int *n, int *nrhs, double complex a[], int *lda, int *ipiv, double complex b[], int *ldb, int *info);
void inv_process(const char solv[], const int NA, const int Noffset, 
		double complex F01[(Noffset+NA)*Noffset], double complex B01[(Noffset+NA)*NA],
		double complex F2[NA*Noffset], double complex B2[NA*NA], double complex CWork[NA*2*NA*2*3])
{
	clock_t startClock, endClock;
	if(strncmp(solv, "Y", 1) == 0 || strncmp(solv, "y", 1) == 0){
		double complex *B01T = &CWork[0];//calloc(NA*(Noffset+NA),sizeof(double complex));
		for(int column = 0 ; column < Noffset+NA ; column++){
			for(int row = 0 ; row < NA ; row++){B01T[row + NA*column] = B01[column + (Noffset+NA)*row];}
		}
		{
			double complex *B2T = &CWork[NA*(Noffset+NA)];//calloc(NA*NA,sizeof(double complex));
			for(int column = 0 ; column < NA ; column++){
				for(int row = 0 ; row < NA ; row++){B2T[row + NA*column] = B2[column + NA*row];}
			}
			int n = NA;
			int nrhs = NA + Noffset;
			int *ipiv = calloc(n,sizeof(int));
			/*long n = NA;
			long nrhs = NA + Noffset;
			long *ipiv = calloc(n,sizeof(long));*/
			int info = -1;	startClock = clock();
			zgesv_(&n, &nrhs, B2T, &n, ipiv, B01T, &n, &info);	endClock = clock();
			SAFEFREE(ipiv);//SAFEFREE(B2T);
			if(info != 0){fprintf(stderr,"info = %d for zgesv in inv_process!\n", info);}
		}
		fprintf(stderr,"zgesv( %1.6f s ), ", (endClock - startClock)/(double)CLOCKS_PER_SEC);
		for(int column = 0 ; column < NA ; column++){
			for(int row = 0 ; row < Noffset+NA ; row++){B01[row + (Noffset+NA)*column] = B01T[column + NA*row];}
		}
		//SAFEFREE(B01T);
	}else{	startClock = clock();
		for(int column0 = NA-1 ; column0 >= 0 ; column0--){
			double maxB = fabs(B2[column0 + NA*column0]);
			int maxc = column0;
			for(int column1 = column0-1 ; column1 >= 0 ; column1--){
				if(maxB < fabs(B2[column0 + NA*column1])) {	maxB = fabs(B2[column0 + NA*column1]);	maxc = column1;}
			}
			if( column0 != maxc){
				double complex dummy;
				for(int row = 0 ; row < Noffset + NA ; row++){
					dummy = B01[row + (Noffset+NA)*column0];
					B01[row + (Noffset+NA)*column0] = B01[row + (Noffset+NA)*maxc];
					B01[row + (Noffset+NA)*maxc] = dummy;
				}
				for(int row = 0 ; row <= column0 ; row++){
					dummy = B2[row + NA*column0];
					B2[row + NA*column0] = B2[row + NA*maxc];
					B2[row + NA*maxc] = dummy;
				}
			}
			//double complex zn = (1. + I*0.)/B2[column0 + NA*column0];
			double complex zn = conj(B2[column0 + NA*column0])/(creal(B2[column0 + NA*column0])*creal(B2[column0 + NA*column0])+cimag(B2[column0 + NA*column0])*cimag(B2[column0 + NA*column0]));// revised on 20190619
			for(int row = 0 ; row < Noffset + NA ; row++){B01[row + (Noffset+NA)*column0] *= zn;}
			for(int row = 0 ; row < column0 ; row++){B2[row + NA*column0] *= zn;}
			for(int column1 = column0-1 ; column1 >= 0 ; column1--){
				zn = B2[column0 + NA*column1];
				for(int row = 0 ; row < Noffset + NA ; row++){B01[row +(Noffset+NA)*column1] -= zn*B01[row + (Noffset+NA)*column0];}
				for(int row = 0 ; row < column0 ; row++){B2[row + NA*column1] -= zn*B2[row + NA*column0];}
			}
		}
		for(int column0 = 0 ; column0 < NA-1 ; column0++){
			for(int column1 = column0+1 ; column1 < NA ; column1++){
				double complex zn = B2[column0 + NA*column1];
				for(int row = 0 ; row < Noffset + NA ; row++){B01[row + (Noffset+NA)*column1] -= zn*B01[row + (Noffset+NA)*column0];}
			}
		}
		endClock = clock();
		fprintf(stderr,"zaxpy_( %1.6f s ), ", (endClock - startClock)/(double)CLOCKS_PER_SEC);
	}
	double complex *Xdummy = &CWork[0];//calloc((Noffset+NA)*Noffset,sizeof(double complex));
	{//B01[(Noffset+NA)*NA]*F2[NA*Noffset] = Xdummy[(Noffset+NA)*Noffset]
		int m = NA + Noffset;
		int n = Noffset;
		int k = NA;
		double complex alpha = 1.;
		double complex beta = 0.;	startClock = clock();
		zgemm_("N", "N", &m, &n, &k, &alpha, &B01[0], &m, &F2[0], &k, &beta, &Xdummy[0], &m);	endClock = clock();
	}	fprintf(stderr,"zgemm( %1.6f s )", (endClock - startClock)/(double)CLOCKS_PER_SEC);
	for(int j = 0 ; j < (Noffset+NA)*Noffset ; j++){F01[j] -= Xdummy[j];}
	//SAFEFREE(Xdummy);
	for(int j = 0 ; j < NA*Noffset ; j++){F2[j] = 0.;}
	for(int j = 0 ; j < NA*NA ; j++){B2[j] = 0.;}
	for(int column = 0 ; column < NA ; column++){B2[column + NA*column] = 1.;}
}
void set_tilPhi(const char Xi_file[], const int iway, const int NA, const int Nimag, double complex tilPhi[]);
void set_Phi(const char Xi_file[], const int iway, const int NA, const int Nimag, double complex Phi[]);
extern double dznrm2_(const int *n, const double complex *x, const int *incx);
//======================================================================
//  set_iniMat                           Last updated on Jan 10, 2019   
//======================================================================
void set_iniMat(const char Xi_file[], const int iway, const int NA, const int Noffset, 
		double complex F01[(Noffset+NA)*Noffset], double complex B01[(Noffset+NA)*NA],
		double complex F2[NA*Noffset], double complex B2[NA*NA], double complex CWork[NA*2*NA*2*3])
{
	for(int column = 0 ; column < NA ; column++){
		for(int row = 0 ; row < Noffset ; row++){
			B01[row + (Noffset+NA)*column] = 0. + I*0.;
		}
	}
	for(int column = 0 ; column < Noffset ; column++){
		for(int row = 0 ; row < Noffset ; row++){
			F01[row + (Noffset+NA)*column] = 0. + I*0.;
		}
		B01[column + (Noffset+NA)*column] = 1. + I*0.;
	}
	{//	size_t n_Work = 2*NA*2*NA*3;//	double complex *CWork = calloc(n_Work,sizeof(double complex));//revised on 20190108
		double complex *Phi = &CWork[0];//double complex *Phi = calloc(2*NA*2*NA,sizeof(double complex));//revised on 20190108
		set_Phi(Xi_file, iway, NA, Noffset, Phi);
		for(int column = 0 ; column < Noffset ; column++){
			for(int row = 0 ; row < NA ; row++){
				F01[(Noffset + row) + (Noffset+NA)*column] = Phi[row + (2*NA)*column];
				F2[row + NA*column] = Phi[row+NA + (2*NA)*column];
			}
		}
		for(int column = 0 ; column < NA ; column++){
			for(int row = 0 ; row < NA ; row++){
				B01[(Noffset + row) + (Noffset+NA)*column] = Phi[row + (2*NA)*(column+NA)];
				B2[row + NA*column] = Phi[row+NA + (2*NA)*(column+NA)];
			}
		}
		{
			double complex *tilPhi = &CWork[2*NA*2*NA];//double complex *tilPhi = calloc(2*NA*2*NA,sizeof(double complex));//revised on 20190108
			set_tilPhi(Xi_file, iway, NA, Noffset, tilPhi);
			double complex *C = &CWork[2*NA*2*NA*2];//double complex *C = calloc(2*NA*2*NA,sizeof(double complex));//revised on 20190108
			{
				int m = 2*NA;
				int n = 2*NA;
				int k = 2*NA;
				double complex alpha = 1.;
				double complex beta = 0.;
				zgemm_("C", "N", &m, &n, &k, &alpha, &tilPhi[0], &m, &Phi[0], &k, &beta, &C[0], &m);
			}
			for(int column = 0 ; column < 2*NA ; column++){C[column + (2*NA)*column] -= 1. + I*0.;}
			int iN2 = (2*NA)*(2*NA);
			int incx = 1;
			fprintf(stderr,"||tilPhi*Phi - 1||_F = %.5e in set_iniMat\n", dznrm2_(&iN2,&C[0],&incx));//revised on 20181112
			//SAFEFREE(tilPhi);SAFEFREE(C);//revised on 20190108
		}
		//SAFEFREE(Work);//SAFEFREE(Phi);//revised on 20190108
	}
}
//======================================================================
//  set_finMat                           Last updated on Jan 10, 2019   
//======================================================================
void set_finMat(const char Xi_file[], const int iway, const int NA, const int Noffset, const int Nimag,
		double complex F01[(Noffset+NA)*Noffset], double complex B01[(Noffset+NA)*NA],
		double complex F2[NA*Noffset], double complex B2[NA*NA], double complex CWork[NA*2*NA*2*3])
{
//	size_t n_Work = 2*NA*(4*NA+2*Noffset);//	double complex *Work = calloc(n_Work,sizeof(double complex));//revised on 20190108
	double complex *tilPhi = &CWork[0];//calloc(2*NA*2*NA,sizeof(double complex));//revised on 20190108
	set_tilPhi(Xi_file, iway, NA, Nimag, tilPhi);
	double complex *F1B1F2B2 = &CWork[2*NA*2*NA];//calloc((2*NA)*(Noffset+NA),sizeof(double complex));//revised on 20190108
	for(int column = 0 ; column < Noffset+NA ; column++){
		for(int row = 0 ; row < 2*NA ; row++){
			if(column < Noffset){
				if(row < NA){
					F1B1F2B2[row + (2*NA)*column] = F01[row+Noffset + (Noffset+NA)*column];
				}else{
					F1B1F2B2[row + (2*NA)*column] = F2[row - NA + NA*column];
				}
			}else{
				if(row < NA){
					F1B1F2B2[row + (2*NA)*column] = B01[row+Noffset + (Noffset+NA)*(column-Noffset)];
				}else{
					F1B1F2B2[row + (2*NA)*column] = B2[row - NA + NA*(column-Noffset)];
				}
			}
			
		}
	}
	double complex *C = &CWork[2*NA*(3*NA+1*Noffset)];//calloc((2*NA)*(Noffset+NA),sizeof(double complex));//revised on 20190108
	clock_t startClock, endClock;
	{//tilPhi[2NA*2NA]*F1B1F2B2[2NA*(Noffset+NA)] = C[2NA*(Noffset+NA)]
		int m = 2*NA;
		int n = Noffset+NA;
		int k = 2*NA;
		double complex alpha = 1.;
		double complex beta = 0.;	startClock = clock();
		zgemm_("C", "N", &m, &n, &k, &alpha, &tilPhi[0], &m, &F1B1F2B2[0], &k, &beta, &C[0], &m);	endClock = clock();
	}	fprintf(stderr,"zgemm( %1.6f s ) for tilPhi*F1B1F2B2, ", (endClock - startClock)/(double)CLOCKS_PER_SEC);
//	SAFEFREE(tilPhi);SAFEFREE(F1B1F2B2);//revised on 20190108
	for(int column = 0 ; column < Noffset+NA ; column++){
		for(int row = 0 ; row < 2*NA ; row++){
			if(column < Noffset){
				if(row < NA){
					F01[row+Noffset + (Noffset+NA)*column] = C[row + (2*NA)*column];
				}else{
					F2[row - NA + NA*column] = C[row + (2*NA)*column];
				}
			}else{
				if(row < NA){
					B01[row+Noffset + (Noffset+NA)*(column-Noffset)] = C[row + (2*NA)*column];
				}else{
					B2[row - NA + NA*(column-Noffset)] = C[row + (2*NA)*column];
				}
			}
			
		}
	}
//	SAFEFREE(Work);//SAFEFREE(C);//revised on 20190108
}
//======================================================================
//  set_Phi                           Last updated on Nov 05, 2018   
//======================================================================
void set_Phi(const char Xi_file[], const int iway, const int NA, const int Nimag, double complex Phi[2*NA*2*NA])
{
	FILE *fp_i;
	fp_i = fopen(Xi_file,"r");
	if (fp_i == NULL){	fprintf(stderr,"open error!: open %s!\n", Xi_file);	exit(EXIT_FAILURE);}
	char buf[BUFSIZE];	// buffer for fgets
	for(int j = 0 ; j < 2 ; j++){
		if(fgets(buf, sizeof( buf ), fp_i) == NULL){fprintf(stderr,"set_tilPhi error!\n");	exit(EXIT_FAILURE);}
	}
	for(int jcolumn = 0 ; jcolumn < NA ; jcolumn++){
		if(fgets(buf, sizeof( buf ), fp_i)== NULL){fprintf(stderr,"set_tilPhi error!\n");	exit(EXIT_FAILURE);}
		int dummyi = -1;
		if(sscanf(buf,"%*[^[] %*[[] %d", &dummyi) == 1){
			if(dummyi != jcolumn){	fprintf(stderr,"dummyi = %d, jcolumn = %d!\n", dummyi, jcolumn);	exit(EXIT_FAILURE);}
		}
		for(int row = 0 ; row < 2*NA ; row++){
			if(fgets(buf, sizeof( buf ), fp_i) == NULL){fprintf(stderr,"set_tilPhi error!\n");	exit(EXIT_FAILURE);}
			double dummy0, dummy1, dummy2, dummy3;
			if(sscanf(buf,"%lf %lf %lf %lf", &dummy0, &dummy1, &dummy2, &dummy3) == 4){
				if(iway == 1){
					Phi[row + (2*NA)*jcolumn] = dummy0 + I*dummy1;
					Phi[row + (2*NA)*(jcolumn+NA)] = dummy2 + I*dummy3;
				}else{
					Phi[(row+NA)%(2*NA) + (2*NA)*jcolumn] = dummy2 + I*dummy3;
					Phi[(row+NA)%(2*NA) + (2*NA)*(jcolumn+NA)] = dummy0 + I*dummy1;
				}
			}else{fprintf(stderr,"From %s, data cannot be read at column = %d, row = %d!\n", Xi_file, jcolumn, row);	fprintf(stderr,"%s\n", buf);	exit(EXIT_FAILURE);}
		}
	}
	if(fclose(fp_i) != 0) {	fprintf(stderr,"fclose error in set_Phi!\n");	exit(EXIT_FAILURE);}
}
//======================================================================
//  set_tilPhi                           Last updated on Nov 05, 2018   
//======================================================================
void set_tilPhi(const char Xi_file[], const int iway, const int NA, const int Nimag, double complex tilPhi[2*NA*2*NA])
{
	FILE *fp_i;
	fp_i = fopen(Xi_file,"r");
	if (fp_i == NULL){	fprintf(stderr,"open error!: open %s!\n", Xi_file);	exit(EXIT_FAILURE);}
	char buf[BUFSIZE];	// buffer for fgets
	for(int j = 0 ; j < 2 ; j++){
		if(fgets(buf, sizeof( buf ), fp_i) == NULL){fprintf(stderr,"set_tilPhi error!\n");	exit(EXIT_FAILURE);}
	}
	for(int jcolumn = 0 ; jcolumn < NA ; jcolumn++){
		if(fgets(buf, sizeof( buf ), fp_i)== NULL){fprintf(stderr,"set_tilPhi error!\n");	exit(EXIT_FAILURE);}
		int dummyi = -1;
		if(sscanf(buf,"%*[^[] %*[[] %d", &dummyi) == 1){
			if(dummyi != jcolumn){	fprintf(stderr,"dummyi = %d, jcolumn = %d!\n", dummyi, jcolumn);	exit(EXIT_FAILURE);}
		}
		for(int row = 0 ; row < 2*NA ; row++){
			if(fgets(buf, sizeof( buf ), fp_i) == NULL){fprintf(stderr,"set_tilPhi error!\n");	exit(EXIT_FAILURE);}
			double dummy0, dummy1, dummy2, dummy3;
			if(sscanf(buf,"%lf %lf %lf %lf", &dummy0, &dummy1, &dummy2, &dummy3) == 4){
				if(iway == 1){
					if(jcolumn < Nimag){
						tilPhi[(row+NA)%(2*NA) + (2*NA)*jcolumn] = dummy0 + I*dummy1;
						tilPhi[(row+NA)%(2*NA) + (2*NA)*(jcolumn+NA)] = -(dummy2 + I*dummy3);
					}else{
						tilPhi[(row+NA)%(2*NA) + (2*NA)*(jcolumn+NA)] = dummy0 + I*dummy1;
						tilPhi[(row+NA)%(2*NA) + (2*NA)*jcolumn] = dummy2 + I*dummy3;
					}
				}else{
					if(jcolumn < Nimag){
						tilPhi[row + (2*NA)*jcolumn] = -(dummy2 + I*dummy3);
						tilPhi[row + (2*NA)*(jcolumn+NA)] = dummy0 + I*dummy1;
					}else{
						tilPhi[row + (2*NA)*jcolumn] = dummy0 + I*dummy1;
						tilPhi[row + (2*NA)*(jcolumn+NA)] = dummy2 + I*dummy3;
					}
				}
			}else{fprintf(stderr,"From %s, data cannot be read at column = %d, row = %d!\n", Xi_file, jcolumn, row);	fprintf(stderr,"%s\n", buf);	exit(EXIT_FAILURE);}
		}
	}
	if(fclose(fp_i) != 0) {	fprintf(stderr,"fclose error in set_tilPhi!\n");	exit(EXIT_FAILURE);}
}
//======================================================================
//  get_RealNum05omega                Last updated on Oct 31, 2018   
//======================================================================
double get_RealNum05omega(const char f_prefix[], const char p_target[], const int L[4], int RealNum05[2])
{
	int sw = 1;
	if(strcmp(p_target, &p_target[BUFSIZE]) == 0) {sw = 0;}//
	char file_name[BUFSIZE];
	snprintf(file_name,sizeof(file_name),"%s_bMode.dat", f_prefix);
	double omega;
	for(int j = 0 ; j <= sw ; j++) {
		FILE *fp_i;
		fp_i = fopen(file_name,"r");
		char buf[BUFSIZE];	// buffer for fgets
		int j_count = -2;
		int dummy = 0;
		double omega_dummy;
		while(fgets(buf, sizeof( buf ), fp_i) != NULL && j_count < 0) { 
			if(strncmp(buf, "# Omega", 7) == 0 && sscanf(buf,"%*[^=] %*[=] %lf", &omega_dummy) == 1){j_count++;}
			if(strncmp(buf, "# Eigen-value", 13) == 0 && sscanf(buf,"%*[^=] %*[=] %d %*[^=] %*[=] %d", &RealNum05[j], &dummy) == 2){j_count++;}
		}
		if(j_count != 0 || RealNum05[j] + dummy != 2*L[0]*L[1]) {fprintf(stderr,"RealNum05 cannot be read in get_RealNum05!\n");	exit(EXIT_FAILURE);}
		if(fclose(fp_i) != 0) {	fprintf(stderr,"fclose error in get_RealNum05!\n");	exit(EXIT_FAILURE);}
		fprintf(stderr,"RealNum05 for %s = %d\n", file_name, RealNum05[j]);
		snprintf(file_name,sizeof(file_name),"%s_tMode.dat", f_prefix);
		if(j == 1 && omega != omega_dummy) {fprintf(stderr,"Top omega = %.5e != Bottom omega = %.5e\n", omega_dummy, omega);	exit(EXIT_FAILURE);}
		omega = omega_dummy;
	}
	if(sw == 0){RealNum05[1] = RealNum05[0];}
	return(omega);
}
//======================================================================
//  matdata_file                         Last updated on Aug 29, 2019   
//======================================================================
void matdata_file(const char f_prefix[], const char add_name[], const int Ltot, const int Lout, const int istep, char data_name[])
{
	if(0 <= istep && istep < Lout){
		snprintf(data_name,BUFSIZE*sizeof(char),"%s%sB%d.dat", f_prefix, add_name, istep);// https://www.ipa.go.jp/security/awareness/vendor/programmingv1/b06_02.html
	}else if(Lout <= istep && istep < Ltot - Lout){
		snprintf(data_name,BUFSIZE*sizeof(char),"%s%s%d.dat", f_prefix, add_name, istep - Lout);
	}else if(Ltot - Lout <= istep && istep < Ltot){
		snprintf(data_name,BUFSIZE*sizeof(char),"%s%sT%d.dat", f_prefix, add_name, istep - (Ltot - Lout));
	}else{fprintf(stderr,"istep = %d error!\n", istep);	exit(EXIT_FAILURE);}
}
//======================================================================
//  input_file1                          Last updated on Sep 17, 2019   
//======================================================================
void input_file1(FILE *fp_i1, int L[4])
{//# l < 30, m < 30, n + 2*outer < 50, media_no < 3, outer = 5,
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = -1;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 0) { 
		if(strncmp(buf, "# l <", 5) == 0 && sscanf(buf,"%*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d %*[^=] %*[=] %d", &L[0], &L[1], &L[2], &L[3]) == 4){
			j_count++;//0
		}
	}
	if(j_count != 0) {	fprintf(stderr,"L[4] cannot be read in input_file1!");	exit(EXIT_FAILURE);}
}
//======================================================================
//  input_file0 ---- rm_space, rm_comma                                 
//                           Last updated on Nov 06, 2018   
//======================================================================
void rm_space( char *A );
void rm_comma( char *A );
void input_file0(FILE *fp_i0, int BoundC[], char f_prefix[], char p_target[], char S_ctrl[], char solv[])
{
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = -4;
	while(fgets(buf, sizeof( buf ), fp_i0) != NULL && j_count < 0) { 
		rm_space(buf);
		if(strncmp(buf, "Prefix", 6) == 0 && sscanf(buf,"%*[^=] %*[=] %s", f_prefix) == 1){
			rm_comma(f_prefix);
			j_count++;//-3
		}else if(strncmp(buf, "BoundaryCondition", 17) == 0 && sscanf(buf,"%*[^=] %*[=] %d %*[^=] %*[=] %d", &BoundC[0], &BoundC[1]) == 2) {
			j_count++;//-2
		}else if(strncmp(buf, "ModeP", 5) == 0 && sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %*[^=] %*[=] %s", p_target) == 1) {
			rm_comma(p_target);
			if(sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %*[^=] %*[=] %*[^=] %*[=] %s", &p_target[BUFSIZE]) != 1){sscanf(p_target,"%s", &p_target[BUFSIZE]);}
			else{rm_comma(&p_target[BUFSIZE]);}
			j_count++;//-1
		}else if(strncmp(buf, "Smatrix", 7) == 0 && sscanf(buf,"%*[^=] %*[=] %s %*[^=] %*[=] %s", S_ctrl, solv) == 2) {
			rm_comma(S_ctrl);rm_comma(solv);
			j_count++;//0
		}
	}
	if(j_count != 0) {	fprintf(stderr,"4 control commands cannot be read in input_file1!");	exit(EXIT_FAILURE);}
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
//  set_M ---- set_D                     Last updated on Jan 10, 2019.  
//======================================================================
void set_D(const int BoundC[2], const int L[2], double Dup0[], double Dup1[], double Ddown0[], double Ddown1[]);
void set_M(const int BoundC[2], const int L[2], const double complex Vcomp[6*L[0]*L[1]], 
	double complex maa[2*L[0]*L[1]*2*L[0]*L[1]], double complex mbb[2*L[0]*L[1]*2*L[0]*L[1]],
	double complex CWork[2*L[0]*L[1]*2*2*L[0]*L[1]*2*3])
{
	long NA1 = L[0]*L[1];
	long Nshift = NA1*NA1;
	double *Dall = calloc(Nshift*4,sizeof(double));
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
	SAFEFREE(Dall);
}
//======================================================================
//  set_D                                Last updated on Oct 06, 2018.  
//======================================================================
void set_D(const int BoundC[2], const int L[2],
	double Dup0[L[0]*L[1]*L[0]*L[1]], double Dup1[L[0]*L[1]*L[0]*L[1]],
	double Ddown0[L[0]*L[1]*L[0]*L[1]], double Ddown1[L[0]*L[1]*L[0]*L[1]])
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
//  set_V                                Last updated on Jun 19, 2019.  
//======================================================================
long set_V(const char data_name[], const int L[2], const double omega, double complex Vcomp[6*L[0]*L[1]])
{
	FILE *fp_i;
	fp_i = fopen(data_name,"r");	if (fp_i == NULL){	fprintf(stderr,"open error! in set_V\n");	exit(EXIT_FAILURE);}
	char buf[BUFSIZE];
	while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
		rm_space(buf);
		if(strncmp(buf, "# l m 1", 7) == 0){goto start_data;}
	}
	fprintf(stderr,"set_V can not read '# l m 1'!");	exit(EXIT_FAILURE);
	start_data:;
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
	if(fclose(fp_i) != 0) {	fprintf(stderr,"fclose error in set_V!\n");	exit(EXIT_FAILURE);}
	for(long j = 0 ; j < 6*L[0]*L[1] ; j++){Vcomp[j] *= omega;}
	for(long j = 0 ; j < L[0]*L[1] ; j++){
//		Vcomp[j + L[0]*L[1]*2] = 1./Vcomp[j + L[0]*L[1]*2];
//		Vcomp[j + L[0]*L[1]*5] = 1./Vcomp[j + L[0]*L[1]*5];
		Vcomp[j + L[0]*L[1]*2] = conj(Vcomp[j + L[0]*L[1]*2])/(creal(Vcomp[j + L[0]*L[1]*2])*creal(Vcomp[j + L[0]*L[1]*2])+cimag(Vcomp[j + L[0]*L[1]*2])*cimag(Vcomp[j + L[0]*L[1]*2]));//revised on 20190619
		Vcomp[j + L[0]*L[1]*5] = conj(Vcomp[j + L[0]*L[1]*5])/(creal(Vcomp[j + L[0]*L[1]*5])*creal(Vcomp[j + L[0]*L[1]*5])+cimag(Vcomp[j + L[0]*L[1]*5])*cimag(Vcomp[j + L[0]*L[1]*5]));//revised on 20190619
	}
	long comp_count = 0;
	for(long j = 0 ; j < 6*L[0]*L[1] ; j++){if(cimag(Vcomp[j]) != 0.){comp_count++;}}
	return(comp_count);
}
//======================================================================
//  check_Smat                             Last updated on Nov 06, 2018 
//======================================================================
void check_Smat(const int NS, const double complex Smatrix[NS*NS], double norm[2])
{
	double complex *C = calloc(NS*NS,sizeof(double complex));
	{
		int m = NS;
		int n = NS;
		int k = NS;
		double complex alpha = 1.;
		double complex beta = 0.;
		zgemm_("C", "N", &m, &n, &k, &alpha, &Smatrix[0], &m, &Smatrix[0], &k, &beta, &C[0], &m);
	}
	for(int row = 0 ; row < NS ; row++){C[row + NS*row] -= 1.;}
	{
		int NS2 = NS*NS;
		int incx = 1;
		norm[0] = dznrm2_(&NS2,&C[0],&incx);
	}
	for(int column = 0 ; column < NS ; column++){
		for(int row = 0 ; row < NS ; row++){
			C[row + NS*column] = Smatrix[row + NS*column] - Smatrix[column + NS*row];
		}
	}	
	{
		int NS2 = NS*NS;
		int incx = 1;
		norm[1] = dznrm2_(&NS2,&C[0],&incx);
	}
	SAFEFREE(C);
}
