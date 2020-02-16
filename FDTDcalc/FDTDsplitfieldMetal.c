//  Start of FDTD.h
#include "header_macro.h"
#include "constant.h"
int S_calc(const int L[4], const int iB, const char f_prefix[], const double complex HB[2*L[0]*L[1]], const double complex EB[2*L[0]*L[1]], double complex SB[4*L[0]*L[1]]);
void main_calc(const int L[4], const double Ah[3*L[0]*L[1]*L[2]*2], const double Bh[3*L[0]*L[1]*L[2]*2], 
				const double Ae[3*L[0]*L[1]*L[2]*2], const double Be[3*L[0]*L[1]*L[2]*2], 
				const int iBc[3], 
				double H_t[3*L[0]*L[1]*L[2]*2], double E_t[3*L[0]*L[1]*L[2]*2], 
				const int iBTIin[4], const double Hi[2*L[0]*L[1]], const double Ei[2*L[0]*L[1]], //int iBTIin[4] = {iB, iT, iI, inci_mode};
				const double complex ExpN2[2], 
				double complex HB[2*L[0]*L[1]], double complex EB[2*L[0]*L[1]], double complex HT[2*L[0]*L[1]], double complex ET[2*L[0]*L[1]],
				const char Cont[L[0]*L[1]*L[2]], const int address_notn[L[0]*L[1]*L[2]]);
double set_env(const int N_2pi, const int time_period[3], const int j_half, const int j_time, const int j_period);
void set_tau(const double omega, const double min_Dt, int *N_2pi, double *invtau_R, double *invtau_I);
void set_incident(const int iB, const int iT, const int iI, const char f_prefix[], const int L[4], const int inci_mode, double complex Hinc[2*L[0]*L[1]], double complex Einc[2*L[0]*L[1]]);
double check_Xi(const int NS, const double complex tilXi[NS*NS], const double complex Xi[NS*NS]);
void set_tilXi(const char Xi_file[], const int NA, const int Nimag, double complex tilXi[2*NA*2*NA]);
void set_Xi(const char Xi_file[], const int NA, const int Nimag, double complex Xi[2*NA*2*NA]);
int get_RealNum(const char file_name[], const int L[4]);
void input_data(const char data_name[], const int L[4], double *min_sqDt, double complex mu[3*L[0]*L[1]], double complex epsi[3*L[0]*L[1]]);
void matdata_file(const char f_prefix[], const int L[4], const int istep, char data_name[]);
int num_u2(const int L[5], const char buf[]);
void input_Lf(const char Yee[], int L[4], double *omega);
void input_file0(FILE *fp_i0, int BoundC[], char f_prefix[], char p_target[], int time_period[3], char inci_point[], int *inci_mode, double *Dt_factor);
void set_sigma(const int L[4], const double complex mu[3*L[0]*L[1]*L[2]], const double complex epsi[3*L[0]*L[1]*L[2]], double sigma[2*(L[0]+L[1]+L[2])]);
void remove_sigma(const int L[4], double complex mu[3*L[0]*L[1]*L[2]], double complex epsi[3*L[0]*L[1]*L[2]], const double sigma[2*(L[0]+L[1]+L[2])]);
void set_address(const int L[4], const double sigma[2*(L[0]+L[1]+L[2])], const double complex mu[3*L[0]*L[1]*L[2]], const double complex epsi[3*L[0]*L[1]*L[2]], char Cont[L[0]*L[1]*L[2]], int *Lnotn, int address_notn[L[0]*L[1]*L[2]]);
//  End of FDTD.h
//======================================================================
//   main ---- input_file0, input_Lf              
//======================================================================
int main(int argc, char **argv)
{
	time_t timer_ini = time(0);	fprintf(stderr,"# Start time of FDTD_SpiltField = %s\n", ctime(&timer_ini));
	if(argc != 2) {fprintf(stderr,"error: number of files \n");	exit(EXIT_FAILURE);}
	else if(strncmp(argv[1], "-v", 2) == 0 || strcmp(argv[1], "--version") == 0 ) {
		fprintf(stderr,"The '%s' creates scattering matrix column vector by FDTD.\n", argv[0]);
		fprintf(stderr,"Version 19.10.29 is compiled at %s on %s.\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(stderr," Source code : '%s'\n Author      : Tatsuya Usuki\n URL         : http://www.smatran.org\n", __FILE__);
		fprintf(stderr," Reference   : 'Finite Difference Time Domain' in 'Formulation.pdf' on Aug 25, 2019.\n");
		fprintf(stderr,"There is NO warranty.\n");
		exit(EXIT_SUCCESS);//normal end
	}
//------- begin reading fundamental parameters -------
	char f_prefix[BUFSIZE], p_target[BUFSIZE*2];
	int BoundC[3], time_period[3];	//BoundC[2] = 1;//	BoundC[2] = 0;
	char inci_point[BUFSIZE];
	int inci_mode = 0;
	double Dt_factor = 1.;
	{
		FILE *fp_i0;
		fp_i0 = fopen(argv[1],"r");
		if (fp_i0 == NULL){	fprintf(stderr,"open error!: open input-file0!\n");	exit(EXIT_FAILURE);}
		input_file0(fp_i0, BoundC, f_prefix, p_target, time_period, inci_point, &inci_mode, &Dt_factor);
		if(fclose(fp_i0) != 0) {	fprintf(stderr,"fclose error after input_file0!\n");	exit(EXIT_FAILURE);}
	}
	int L[4];//L[0] = L0, L[1] = L1, L[2] = L2 + 2*outer, L[3] = outer
	double omega = 0.;
	input_Lf(f_prefix, L, &omega);
	fprintf(stderr,"Ramp-period = %d, cw-period = %d, stop-period = %d\n", time_period[0], time_period[1], time_period[2]);
	fprintf(stderr,"Dt_factor = %.5e, split field impl. does not use stable_factor.\n", Dt_factor);
//------------------- Setting mu, eqsilon and min_Dt -------------------
	double min_Dt = 0.;
	double complex *mu = calloc(3*L[0]*L[1]*L[2], sizeof(double complex));	double complex *epsi = calloc(3*L[0]*L[1]*L[2], sizeof(double complex));
	for(int istep = 0 ; istep < L[2] ; istep++){
		char data_name[BUFSIZE];
		matdata_file(f_prefix, L, istep, data_name);
		double min_sqDt;
		input_data(data_name, L, &min_sqDt, &mu[3*L[0]*L[1]*istep], &epsi[3*L[0]*L[1]*istep]);
		if(istep == 0 || min_Dt > sqrt(min_sqDt)){min_Dt = sqrt(min_sqDt);}
	}
//------------------- Setting sigma -------------------
	double *sigma = calloc(2*(L[0]+L[1]+L[2]), sizeof(double));
	for(int l = 0 ; l < 2*(L[0]+L[1]+L[2]) ; l++){sigma[l] = 0.;}
	set_sigma(L, mu, epsi, sigma);
	remove_sigma(L, mu, epsi, sigma);
//-------------- Setting monitor points and incident point --------------------
	int iB, iT, iI;
	iB = num_u2(L, p_target);
	iT = num_u2(L, &p_target[BUFSIZE]);
	iI = num_u2(L, inci_point);
	if(iB == iT){	iT = (L[2] - 1) - iB;	fprintf(stderr,"Change iT because iT = iB\n");}
	if(iB >= iT){fprintf(stderr,"Target_B = %d for %s is larger than or equal to Target_T = %d for %s\n", iB, p_target, iT, &p_target[BUFSIZE]);	exit(EXIT_FAILURE);}
	int iBTIin[4] = {iB, iT, iI, inci_mode};
	fprintf(stderr,"Start = %d for %s, Incident point = %d for %s, Incident mode = %+d, Target_T = %d for %s\n", iB, p_target, iI, inci_point, inci_mode, iT, &p_target[BUFSIZE]);
//----------------------- Setting incident mode ------------------------------
	double complex *Hinc = calloc(2*L[0]*L[1], sizeof(double complex));	double complex *Einc = calloc(2*L[0]*L[1], sizeof(double complex));
	set_incident(iB, iT, iI, f_prefix, L, inci_mode, Hinc, Einc);
//---------------------- Setting time parameters ------------------------------
	int N_2pi;
	double invtau_R, invtau_I;
	set_tau(omega, min_Dt/Dt_factor, &N_2pi, &invtau_R, &invtau_I);
	double complex *ExpN = calloc(2*N_2pi, sizeof(double complex));
	{double Dphi = Pi2/((double) N_2pi);
		for(int j_time = 0 ; j_time < N_2pi ; j_time++){
			double phi00 = Dphi*((double) j_time);	ExpN[0 + 2*j_time] = exp(I*phi00);
			double phi05 = Dphi*(0.5 + ((double) j_time));	ExpN[1 + 2*j_time] = exp(I*phi05);
		}
	}
//---------------------- Setting equation coefficients ------------------------------
	char *Cont = calloc(L[0]*L[1]*L[2], sizeof(char));
	int Lnotn;
	int *address_notn = calloc(L[0]*L[1]*L[2], sizeof(int));
	set_address(L, sigma, mu, epsi, Cont, &Lnotn, address_notn);
	fprintf(stderr,"L3d = %d, Lnotn = %d\n", L[0]*L[1]*L[2], Lnotn);
	double *Ah = calloc(3*(L[0]*L[1]*L[2]+Lnotn), sizeof(double));	double *Bh = calloc(3*(L[0]*L[1]*L[2]+Lnotn), sizeof(double));
	double *Ae = calloc(3*(L[0]*L[1]*L[2]+Lnotn), sizeof(double));	double *Be = calloc(3*(L[0]*L[1]*L[2]+Lnotn), sizeof(double));
	double sigma6[2*3];
	for(int j2 = 0 ; j2 < L[2] ; j2++){
		for(int js = 0 ; js < 2 ; js++){sigma6[js+2*2] = sigma[2*L[0]+2*L[1]+js+2*j2];}
		for(int j1 = 0 ; j1 < L[1] ; j1++){
			for(int js = 0 ; js < 2 ; js++){sigma6[js+2*1] = sigma[2*L[0]+js+2*j1];}
			for(int j0 = 0 ; j0 < L[0] ; j0++){
				for(int js = 0 ; js < 2 ; js++){sigma6[js+2*0] = sigma[js+2*j0];}
				int l = (j0 + L[0]*j1 + L[0]*L[1]*j2);
				for(int jxyz = 0 ; jxyz < 3 ; jxyz++){
					int l3d = jxyz + 3*l;
					Bh[l3d] = 1./(invtau_R*creal(  mu[l3d]) + invtau_I*(creal(  mu[l3d])*sigma6[0+2*((jxyz+1)%3)] + cimag(  mu[l3d])));// -1+3 = 2
					Ah[l3d] = (invtau_R*creal(  mu[l3d]) - invtau_I*(creal(  mu[l3d])*sigma6[0+2*((jxyz+1)%3)] + cimag(  mu[l3d])))*Bh[l3d];// -1+3 = 2
					Be[l3d] = 1./(invtau_R*creal(epsi[l3d]) + invtau_I*(creal(epsi[l3d])*sigma6[1+2*((jxyz+1)%3)] + cimag(epsi[l3d])));// -1+3 = 2
					Ae[l3d] = (invtau_R*creal(epsi[l3d]) - invtau_I*(creal(epsi[l3d])*sigma6[1+2*((jxyz+1)%3)] + cimag(epsi[l3d])))*Be[l3d];// -1+3 = 2
					if(Cont[l] == 'a'){//if(address_notn[l] != Lnotn){//
						int L3D = jxyz + 3*address_notn[l] + 3*L[0]*L[1]*L[2];
						Bh[L3D] = 1./(invtau_R*creal(  mu[l3d]) + invtau_I*(creal(  mu[l3d])*sigma6[0+2*((jxyz+2)%3)] + cimag(  mu[l3d])));// -1+3 = 2
						Ah[L3D] = (invtau_R*creal(  mu[l3d]) - invtau_I*(creal(  mu[l3d])*sigma6[0+2*((jxyz+2)%3)] + cimag(  mu[l3d])))*Bh[L3D];// -1+3 = 2
						Be[L3D] = 1./(invtau_R*creal(epsi[l3d]) + invtau_I*(creal(epsi[l3d])*sigma6[1+2*((jxyz+2)%3)] + cimag(epsi[l3d])));// -1+3 = 2
						Ae[L3D] = (invtau_R*creal(epsi[l3d]) - invtau_I*(creal(epsi[l3d])*sigma6[1+2*((jxyz+2)%3)] + cimag(epsi[l3d])))*Be[L3D];// -1+3 = 2
					}
					if(Cont[l] == 'e' || Cont[l] == 'b'){// ./check_files/StableForm.pdf
						int L3D = jxyz + 3*address_notn[l] + 3*L[0]*L[1]*L[2];
						Be[L3D] = 1./(0.5*creal(epsi[l3d]) - invtau_R*(0.5*cimag(epsi[l3d])/invtau_I));
						Ae[L3D] = -(0.5*creal(epsi[l3d]) + invtau_R*(0.5*cimag(epsi[l3d])/invtau_I))*Be[L3D];
					}
					if(Cont[l] == 'h' || Cont[l] == 'b'){// ./check_files/StableForm.pdf
						int L3D = jxyz + 3*address_notn[l] + 3*L[0]*L[1]*L[2];
						Bh[L3D] = 1./(0.5*creal(mu[l3d]) - invtau_R*(0.5*cimag(mu[l3d])/invtau_I));
						Ah[L3D] = -(0.5*creal(mu[l3d]) + invtau_R*(0.5*cimag(mu[l3d])/invtau_I))*Bh[L3D];
						
					}
				}
			}
		}
	}
	SAFEFREE(mu);	SAFEFREE(epsi);	SAFEFREE(sigma);
//---------------------- EDTD calculation ------------------------------
	double sum00 = 0.;	double sum05 = 0.;
	double *Hi = calloc(2*L[0]*L[1], sizeof(double));	double  *Ei = calloc(2*L[0]*L[1], sizeof(double));
	double complex *HB = calloc(2*L[0]*L[1], sizeof(double complex));	double complex *EB = calloc(2*L[0]*L[1], sizeof(double complex));
	double complex *HT = calloc(2*L[0]*L[1], sizeof(double complex));	double complex *ET = calloc(2*L[0]*L[1], sizeof(double complex));
	for(int l = 0 ; l < 2*L[0]*L[1] ; l++){	HB[l] = 0.;	EB[l] = 0.;	HT[l] = 0.;	ET[l] = 0.;}
	double *H_t = calloc(3*(L[0]*L[1]*L[2]+Lnotn+1), sizeof(double));	double *E_t = calloc(3*(L[0]*L[1]*L[2]+Lnotn+1), sizeof(double));//Note that 'Lnotn+1' instead of 'Lnoton' is set to represent 'out of special field.'
	for(int l = 0 ; l < 3*(L[0]*L[1]*L[2]+Lnotn+1) ; l++){H_t[l] = 0.;	E_t[l] = 0.;}
	for(int j_period = 0 ; j_period < time_period[2] ; j_period++){
		for(int j_time = 0 ; j_time < N_2pi ; j_time++){
			for(int l = 0 ; l < 2*L[0]*L[1] ; l++){
				Ei[l] = (2.)*set_env(N_2pi, time_period, 1, j_time, j_period) * creal(conj(ExpN[1 + 2*j_time])*Einc[l]);
				Hi[l] = (2.)*set_env(N_2pi, time_period, 0, (j_time+1), j_period) * creal(conj(ExpN[0 + 2*((j_time+1)%N_2pi)])*Hinc[l]);
			}
			{
				double f_env00 = set_env(N_2pi, time_period, 0, (j_time+1), j_period);
				sum00 += f_env00;//	sum00 += f_env00*f_env00;
				double f_env05 = set_env(N_2pi, time_period, 1, j_time, j_period);
				sum05 += f_env05;//	sum05 += f_env05*f_env05;
			}
			main_calc(L, Ah, Bh, Ae, Be, BoundC, H_t, E_t, iBTIin, Hi, Ei, &ExpN[0 + 2*((j_time+1)%N_2pi)], HB, EB, HT, ET, Cont, address_notn);
		}
	}
	fprintf(stderr,"sum00 = %.5e, sum05 = %.5e\n", sum00, sum05);
	fprintf(stderr,"Bc[0] = %d, Bc[1] = %d, Bc[2] = %d\n", BoundC[0], BoundC[1], BoundC[2]);
	
	SAFEFREE(Hi);	SAFEFREE(Ei);
	SAFEFREE(ExpN);
	SAFEFREE(H_t);	SAFEFREE(E_t);
	SAFEFREE(Hinc);	SAFEFREE(Einc);
	SAFEFREE(Ah);	SAFEFREE(Ae);	SAFEFREE(Bh);	SAFEFREE(Be);
	SAFEFREE(Cont);	SAFEFREE(address_notn);
//------------------------- S-parameter calculation ---------------------------
	double complex *SB = calloc(4*L[0]*L[1], sizeof(double complex));
	double complex *ST = calloc(4*L[0]*L[1], sizeof(double complex));
	int N_B = S_calc(L, iB, f_prefix, HB, EB, SB);	int N_T = S_calc(L, iT, f_prefix, HT, ET, ST);
	SAFEFREE(HT);	SAFEFREE(ET);	SAFEFREE(HB);	SAFEFREE(EB);
	fprintf(stderr,"Bottom |S_(m,%+d)|^2 :\n", inci_mode);
	for(int jmode = 0 ; jmode < N_B ; jmode++){fprintf(stderr,"%+d, %.20e, %+d, %.20e\n", jmode + 1, creal(conj(SB[jmode])*SB[jmode]), -(jmode + 1), creal(conj(SB[jmode+2*L[0]*L[1]])*SB[jmode+2*L[0]*L[1]]));}
	fprintf(stderr,"Top |S_(m,%+d)|^2 :\n", inci_mode);
	for(int jmode = 0 ; jmode < N_T ; jmode++){fprintf(stderr,"%+d, %.20e, %+d, %.20e\n", jmode + 1, creal(conj(ST[jmode])*ST[jmode]), -(jmode + 1), creal(conj(ST[jmode+2*L[0]*L[1]])*ST[jmode+2*L[0]*L[1]]));}
	SAFEFREE(ST);	SAFEFREE(SB);
}
void set_address(const int L[4], const double sigma[2*(L[0]+L[1]+L[2])], const double complex mu[3*L[0]*L[1]*L[2]], const double complex epsi[3*L[0]*L[1]*L[2]], char Cont[L[0]*L[1]*L[2]], int *Lnotn, int address_notn[L[0]*L[1]*L[2]])
{
	double sigma6[2*3];//	for(int j6 = 0 ; j6 < 6 ; j6++){sigma6[j6] = 0.;}
	for(int j2 = 0 ; j2 < L[2] ; j2++){
		for(int js = 0 ; js < 2 ; js++){sigma6[js+2*2] = sigma[2*L[0]+2*L[1]+js+2*j2];}
		for(int j1 = 0 ; j1 < L[1] ; j1++){
			for(int js = 0 ; js < 2 ; js++){sigma6[js+2*1] = sigma[2*L[0]+js+2*j1];}
			for(int j0 = 0 ; j0 < L[0] ; j0++){
				for(int js = 0 ; js < 2 ; js++){sigma6[js+2*0] = sigma[js+2*j0];}
				Cont[j0 + L[0]*j1 + L[0]*L[1]*j2] = 'n';// Normal case
				for(int j6 = 0 ; j6 < 6 ; j6++){if(fabs(sigma6[j6]) > ERROR_MIN) {Cont[j0 + L[0]*j1 + L[0]*L[1]*j2] = 'a';}}
				int l = j0 + L[0]*j1 + L[0]*L[1]*j2;
				for(int jxyz = 0 ; jxyz < 3 ; jxyz++){
					if(tan(carg(mu[jxyz + 3*l])) < -ERROR_MIN && tan(carg(epsi[jxyz + 3*l])) < -ERROR_MIN) {
						Cont[j0 + L[0]*j1 + L[0]*L[1]*j2] = 'b';// Both cases
					}else{
						if(tan(carg(mu[jxyz + 3*l])) < -ERROR_MIN) {
							Cont[j0 + L[0]*j1 + L[0]*L[1]*j2] = 'h';// Abnormal mu
						}else if(tan(carg(epsi[jxyz + 3*l])) < -ERROR_MIN) {
							Cont[j0 + L[0]*j1 + L[0]*L[1]*j2] = 'e';// Metallic epsilon
						}
					}
				}
			}
		}
	}
	*Lnotn = 0;
	for(int j2 = 0 ; j2 < L[2] ; j2++){
		for(int j1 = 0 ; j1 < L[1] ; j1++){
			for(int j0 = 0 ; j0 < L[0] ; j0++){	if(Cont[j0 + L[0]*j1 + L[0]*L[1]*j2] != 'n'){	(*Lnotn)++;}}
		}
	}
	int jcount = 0;
	for(int j2 = 0 ; j2 < L[2] ; j2++){
		for(int j1 = 0 ; j1 < L[1] ; j1++){
			for(int j0 = 0 ; j0 < L[0] ; j0++){
				if(Cont[j0 + L[0]*j1 + L[0]*L[1]*j2] != 'n'){
					address_notn[j0 + L[0]*j1 + L[0]*L[1]*j2] = jcount;
					jcount++;
				}else{// An address (*Lnotn) means region out of special treatment.
					address_notn[j0 + L[0]*j1 + L[0]*L[1]*j2] = (*Lnotn);
				}
			}
		}
	}
}
//======================================================================
//  main_calc                              Last updated on Oct 29, 2019 
//                                                     See `PMLform.lyx' 
//======================================================================
void ope_Rmap(const int js, const int L[4], const int iBc[3], const int j[3], const int address_notn[L[0]*L[1]*L[2]], const double H_p[3*L[0]*L[1]*L[2]], double WorkR[3]);
void ope_Rdmap(const int js, const int L[4], const int iBc[3], const int j[3], const int address_notn[L[0]*L[1]*L[2]], const double E_p[3*L[0]*L[1]*L[2]], double WorkR[3]);
void ope_R(const int js, const int L[4], const int iBc[3], const int j[3], const double H_p[3*L[0]*L[1]*L[2]], double WorkR[3]);
void ope_Rd(const int js, const int L[4], const int iBc[3], const int j[3], const double E_p[3*L[0]*L[1]*L[2]], double WorkR[3]);
void main_calc(const int L[4], const double Ah[3*L[0]*L[1]*L[2]*2], const double Bh[3*L[0]*L[1]*L[2]*2], 
				const double Ae[3*L[0]*L[1]*L[2]*2], const double Be[3*L[0]*L[1]*L[2]*2], 
				const int iBc[3], 
				double H_t[3*L[0]*L[1]*L[2]*2], double E_t[3*L[0]*L[1]*L[2]*2], 
				const int iBTIin[4], const double Hi[2*L[0]*L[1]], const double Ei[2*L[0]*L[1]], //int iBTIin[4] = {iB, iT, iI, inci_mode};
				const double complex ExpN2[2], 
				double complex HB[2*L[0]*L[1]], double complex EB[2*L[0]*L[1]], double complex HT[2*L[0]*L[1]], double complex ET[2*L[0]*L[1]],
				const char Cont[L[0]*L[1]*L[2]], const int address_notn[L[0]*L[1]*L[2]])
{
	double WorkR[3];
	int L3D = 3*L[0]*L[1]*L[2];
	double mabsm = (double) (iBTIin[3]/abs(iBTIin[3]));
	int j[3];
	for(j[2] = 0 ; j[2] < L[2] ; j[2]++){// H-field step
		for(j[1] = 0 ; j[1] < L[1] ; j[1]++){
			for(j[0] = 0 ; j[0] < L[0] ; j[0]++){
				int l = j[0] + L[0]*j[1] + L[0]*L[1]*j[2];
				if(Cont[l] == 'n'){
					ope_Rd(2, L, iBc, j, E_t, WorkR);
					if(j[2] == iBTIin[2]){
						WorkR[0] += mabsm*Ei[1 + 2*j[0] + 2*L[0]*j[1]];
						WorkR[1] -= mabsm*Ei[0 + 2*j[0] + 2*L[0]*j[1]];
					}
					for(int jxyz = 0 ; jxyz < 3 ; jxyz++){H_t[jxyz+3*l] = Ah[jxyz+3*l]*H_t[jxyz+3*l] + Bh[jxyz+3*l]*(-WorkR[jxyz]);}
					ope_Rdmap(2, L, iBc, j, address_notn, &E_t[L3D], WorkR);
					for(int jxyz = 0 ; jxyz < 3 ; jxyz++){H_t[jxyz+3*l] += Bh[jxyz+3*l]*(-WorkR[jxyz]);}
				}else if(Cont[l] == 'a'){//See './check_files/PMLForm.pdf'
					for(int mp = 0 ; mp < 2 ; mp++){// mp=0 means "+", mp=1 means "-".
						ope_Rd(mp, L, iBc, j, E_t, WorkR);
						int ja = 3*l;	if(mp == 1){ja = L3D + 3*address_notn[l];}
						for(int jxyz = 0 ; jxyz < 3 ; jxyz++){H_t[ja+jxyz] = Ah[ja+jxyz]*H_t[ja+jxyz] + Bh[ja+jxyz]*(-WorkR[jxyz]);
						}
						ope_Rdmap(mp, L, iBc, j, address_notn, &E_t[L3D], WorkR);
						for(int jxyz = 0 ; jxyz < 3 ; jxyz++){H_t[ja+jxyz] += Bh[ja+jxyz]*(-WorkR[jxyz]);}
					}
				}else if(Cont[l] == 'h' || Cont[l] == 'b'){//See './check_files/StableForm.pdf'
					int ja = L3D + 3*address_notn[l];
					ope_Rd(2, L, iBc, j, E_t, WorkR);
					for(int jxyz = 0 ; jxyz < 3 ; jxyz++){	H_t[jxyz + 3*l] += 0.5*H_t[jxyz + ja];}
					for(int jxyz = 0 ; jxyz < 3 ; jxyz++){	H_t[jxyz + ja] = Ah[jxyz + ja]*H_t[jxyz + ja] + Bh[jxyz + ja]*(-WorkR[jxyz]);}
					for(int jxyz = 0 ; jxyz < 3 ; jxyz++){	H_t[jxyz + 3*l] += 0.5*H_t[jxyz + ja];}
				}
			}
		}
	}
	for(j[2] = 0 ; j[2] < L[2] ; j[2]++){// E-field step
		for(j[1] = 0 ; j[1] < L[1] ; j[1]++){
			for(j[0] = 0 ; j[0] < L[0] ; j[0]++){
				int l = j[0] + L[0]*j[1] + L[0]*L[1]*j[2];
				if(Cont[l] == 'n'){
					ope_R(2, L, iBc, j, H_t, WorkR);
					if(j[2] == iBTIin[2]){
						WorkR[0] += mabsm*Hi[1 + 2*j[0] + 2*L[0]*j[1]];
						WorkR[1] -= mabsm*Hi[0 + 2*j[0] + 2*L[0]*j[1]];
					}
					for(int jxyz = 0 ; jxyz < 3 ; jxyz++){E_t[jxyz+3*l] = Ae[jxyz+3*l]*E_t[jxyz+3*l] + Be[jxyz+3*l]*WorkR[jxyz];}
					ope_Rmap(2, L, iBc, j, address_notn, &H_t[L3D], WorkR);
					for(int jxyz = 0 ; jxyz < 3 ; jxyz++){E_t[jxyz+3*l] += Be[jxyz+3*l]*WorkR[jxyz];}
				}else if(Cont[l] == 'a'){//See './check_files/PMLForm.pdf'
					for(int mp = 0 ; mp < 2 ; mp++){// mp=0 means "+", mp=1 means "-".
						ope_R(mp, L, iBc, j, H_t, WorkR);
						int ja = 3*l;	if(mp == 1){ja = L3D + 3*address_notn[l];}
						for(int jxyz = 0 ; jxyz < 3 ; jxyz++){E_t[ja+jxyz] = Ae[ja+jxyz]*E_t[ja+jxyz] + Be[ja+jxyz]*WorkR[jxyz];}
						ope_Rmap(mp, L, iBc, j, address_notn, &H_t[L3D], WorkR);
						for(int jxyz = 0 ; jxyz < 3 ; jxyz++){E_t[ja+jxyz] += Be[ja+jxyz]*WorkR[jxyz];}
					}
				}else if(Cont[l] == 'e'|| Cont[l] == 'b'){//See './check_files/StableForm.pdf'
					int ja = L3D + 3*address_notn[l];
					ope_Rd(2, L, iBc, j, H_t, WorkR);
					for(int jxyz = 0 ; jxyz < 3 ; jxyz++){	E_t[jxyz + 3*l] += 0.5*E_t[jxyz + ja];}
					for(int jxyz = 0 ; jxyz < 3 ; jxyz++){	E_t[jxyz + ja] = Ae[jxyz + ja]*E_t[jxyz + ja] + Be[jxyz + ja]*WorkR[jxyz];}
					for(int jxyz = 0 ; jxyz < 3 ; jxyz++){	E_t[jxyz + 3*l] += 0.5*E_t[jxyz + ja];}
				}
			}
		}
	}
	for(int j1 = 0 ; j1 < L[1] ; j1++){//Calc. of monitor fields
		for(int j0 = 0 ; j0 < L[0] ; j0++){
			for(int m = 0 ; m < 2 ; m++){if(Cont[j0 + L[0]*j1 + L[0]*L[1]*iBTIin[m]] == 'a'){fprintf(stderr,"iBTIin[%d] error!\n",m);	exit(EXIT_FAILURE);}}
			for(int jv = 0 ; jv < 2 ; jv++){
				HB[jv + 2*j0 + 2*L[0]*j1] += ExpN2[0]*((double complex) H_t[jv + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*iBTIin[0]]);
				EB[jv + 2*j0 + 2*L[0]*j1] += ExpN2[1]*((double complex) E_t[jv + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*iBTIin[0]]);
				HT[jv + 2*j0 + 2*L[0]*j1] += ExpN2[0]*((double complex) H_t[jv + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*iBTIin[1]]);
				ET[jv + 2*j0 + 2*L[0]*j1] += ExpN2[1]*((double complex) E_t[jv + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*iBTIin[1]]);
			}
		}
	}
}
void ope_Rmap(const int js, const int L[4], const int iBc[3], const int j[3], const int address_notn[L[0]*L[1]*L[2]], const double H_p[3*L[0]*L[1]*L[2]], double WorkR[3])
{
	double Bp[3];	int k[3];
	for(int jxyz = 0 ; jxyz < 3 ; jxyz++){
		Bp[jxyz] = 1.;	if(j[jxyz]+1 == L[jxyz]){Bp[jxyz] = ((double) iBc[jxyz]);}
		k[jxyz] = (j[jxyz]+1)%L[jxyz] - j[jxyz];
		for(int kxyz = 0 ; kxyz < jxyz ; kxyz++){
			k[jxyz] *= L[kxyz];
		}
	}
	int l = j[0] + L[0]*j[1] + L[0]*L[1]*j[2];
	for(int jxyz = 0 ; jxyz < 3 ; jxyz++){
		WorkR[jxyz] = 0.;
		if(js != 1){// js=0 means "+"
			WorkR[jxyz] += +Bp[(jxyz+1)%3]*H_p[(jxyz+2)%3 + 3*address_notn[k[(jxyz+1)%3]+l]] - H_p[(jxyz+2)%3 + 3*address_notn[l]];
		}
		if(js != 0){// js=1 means "-"
			WorkR[jxyz] += -Bp[(jxyz+2)%3]*H_p[(jxyz+1)%3 + 3*address_notn[k[(jxyz+2)%3]+l]] + H_p[(jxyz+1)%3 + 3*address_notn[l]];
		}
	}
}
void ope_Rdmap(const int js, const int L[4], const int iBc[3], const int j[3], const int address_notn[L[0]*L[1]*L[2]], const double E_p[3*L[0]*L[1]*L[2]], double WorkR[3])
{
	double Bp[3];	int k[3];
	for(int jxyz = 0 ; jxyz < 3 ; jxyz++){
		Bp[jxyz] = 1.;	if(j[jxyz] == 0){Bp[jxyz] = ((double) iBc[jxyz]);}
		k[jxyz] = (j[jxyz]-1+L[jxyz])%L[jxyz] - j[jxyz];
		for(int kxyz = 0 ; kxyz < jxyz ; kxyz++){
			k[jxyz] *= L[kxyz];
		}
	}
	int l = j[0] + L[0]*j[1] + L[0]*L[1]*j[2];
	for(int jxyz = 0 ; jxyz < 3 ; jxyz++){
		WorkR[jxyz] = 0.;
		if(js != 1){// js=0 means "+"
			WorkR[jxyz] += +E_p[(jxyz+2)%3 + 3*address_notn[l]] - Bp[(jxyz+1)%3]*E_p[(jxyz+2)%3 + 3*address_notn[k[(jxyz+1)%3]+l]];
		}
		if(js != 0){// js=1 means "-"
			WorkR[jxyz] += -E_p[(jxyz+1)%3 + 3*address_notn[l]] + Bp[(jxyz+2)%3]*E_p[(jxyz+1)%3 + 3*address_notn[k[(jxyz+2)%3]+l]];
		}
	}
}
void ope_R(const int js, const int L[4], const int iBc[3], const int j[3], const double H_p[3*L[0]*L[1]*L[2]], double WorkR[3])
{
	double Bp[3];	int k[3];
	for(int jxyz = 0 ; jxyz < 3 ; jxyz++){
		Bp[jxyz] = 1.;	if(j[jxyz]+1 == L[jxyz]){Bp[jxyz] = ((double) iBc[jxyz]);}
		k[jxyz] = (j[jxyz]+1)%L[jxyz] - j[jxyz];
		for(int kxyz = 0 ; kxyz < jxyz ; kxyz++){
			k[jxyz] *= L[kxyz];
		}
	}
	int l = j[0] + L[0]*j[1] + L[0]*L[1]*j[2];
	for(int jxyz = 0 ; jxyz < 3 ; jxyz++){
		WorkR[jxyz] = 0.;
		if(js != 1){// js=0 means "+"
			WorkR[jxyz] += +Bp[(jxyz+1)%3]*H_p[(jxyz+2)%3 + 3*(k[(jxyz+1)%3]+l)] - H_p[(jxyz+2)%3 + 3*l];
		}
		if(js != 0){// js=1 means "-"
			WorkR[jxyz] += -Bp[(jxyz+2)%3]*H_p[(jxyz+1)%3 + 3*(k[(jxyz+2)%3]+l)] + H_p[(jxyz+1)%3 + 3*l];
		}
	}
}
void ope_Rd(const int js, const int L[4], const int iBc[3], const int j[3], const double E_p[3*L[0]*L[1]*L[2]], double WorkR[3])
{
	double Bp[3];	int k[3];
	for(int jxyz = 0 ; jxyz < 3 ; jxyz++){
		Bp[jxyz] = 1.;	if(j[jxyz] == 0){Bp[jxyz] = ((double) iBc[jxyz]);}
		k[jxyz] = (j[jxyz]-1+L[jxyz])%L[jxyz] - j[jxyz];
		for(int kxyz = 0 ; kxyz < jxyz ; kxyz++){
			k[jxyz] *= L[kxyz];
		}
	}
	int l = j[0] + L[0]*j[1] + L[0]*L[1]*j[2];
	for(int jxyz = 0 ; jxyz < 3 ; jxyz++){
		WorkR[jxyz] = 0.;
		if(js != 1){// js=0 means "+"
			WorkR[jxyz] += +E_p[(jxyz+2)%3 + 3*l] - Bp[(jxyz+1)%3]*E_p[(jxyz+2)%3 + 3*(k[(jxyz+1)%3]+l)];
		}
		if(js != 0){// js=1 means "-"
			WorkR[jxyz] += -E_p[(jxyz+1)%3 + 3*l] + Bp[(jxyz+2)%3]*E_p[(jxyz+1)%3 + 3*(k[(jxyz+2)%3]+l)];
		}
	}
}
//======================================================================
//  S_calc  ---- get_RealNum, set_tilXi    Last updated on Sep 30, 2019 
//======================================================================
void set_tilXi(const char Xi_file[], const int NA, const int Nimag, double complex tilXi[2*NA*2*NA]);
int get_RealNum(const char file_name[], const int L[4]);
int S_calc(const int L[4], const int iB, const char f_prefix[], const double complex HB[2*L[0]*L[1]], const double complex EB[2*L[0]*L[1]], double complex SB[4*L[0]*L[1]])
{
	char Xi_file[BUFSIZE], mode_file[BUFSIZE];
	if(0 <= iB && iB < L[3]){
		snprintf(Xi_file,sizeof(Xi_file),"%s_bXi.dat", f_prefix);
		snprintf(mode_file,sizeof(mode_file),"%s_bMode.dat", f_prefix);
/*	}else if(L[3] <= iB && iB < L[2]-L[3]){
		snprintf(Xi_file,sizeof(Xi_file),"%s_tXi.dat", f_prefix);
		snprintf(mode_file,sizeof(mode_file),"%s_tMode.dat", f_prefix);*/
	}else if(L[2]-L[3] <= iB && iB < L[2]){
		snprintf(Xi_file,sizeof(Xi_file),"%s_tXi.dat", f_prefix);
		snprintf(mode_file,sizeof(mode_file),"%s_tMode.dat", f_prefix);
	}else{fprintf(stderr,"iB or iT = %d \n", iB);	exit(EXIT_FAILURE);}
	
	int NA = 2*L[0]*L[1];//	assert(Noffset <= NA);
	int Nimag = get_RealNum(mode_file, L);
	double complex *tilXi = calloc(2*NA*2*NA,sizeof(double complex));
	set_tilXi(Xi_file, NA, Nimag, tilXi);
	for(int column_num = 0 ; column_num < 4*L[0]*L[1] ; column_num++){
		SB[column_num] = 0.;
		for(int irow = 0 ; irow < L[0]*L[1] ; irow++){
			SB[column_num] += conj(tilXi[irow + (2*NA)*column_num])*HB[0 + 2*irow];//revised on 20190930
			SB[column_num] += conj(tilXi[irow+L[0]*L[1] + (2*NA)*column_num])*HB[1 + 2*irow] ;//revised on 20190930
			SB[column_num] += conj(tilXi[irow+NA + (2*NA)*column_num])*(-EB[1 + 2*irow]);//revised on 20190930
			SB[column_num] += conj(tilXi[irow+L[0]*L[1]+NA + (2*NA)*column_num])*EB[0 + 2*irow];//revised on 20190930
		}
	}
	SAFEFREE(tilXi);
	return(Nimag);
}
//======================================================================
//  set_env                                Last updated on Sep 26, 2019 
//======================================================================
double set_env(const int N_2pi, const int time_period[3], const int j_half, const int j_time, const int j_period)
{
	double f_env;
	if(j_period < time_period[0]){
		double Dphi = Pi2/((double) 4*time_period[0]);
		double Dphis = Pi2/((double) 4*time_period[0]*N_2pi);
		double complex ExpP = exp(I*Dphi*((double) j_period));
		double complex ExpPs = exp(I*Dphis*(0.5*((double) j_half) + ((double) j_time)));
		f_env = cimag(ExpP*ExpPs);
	}else if(time_period[0] <= j_period && j_period < time_period[0]+time_period[1]){
		f_env = 1.;
	}else if(time_period[0]+time_period[1]<= j_period && j_period < 2*time_period[0]+time_period[1]){
		double Dphi = Pi2/((double) 4*time_period[0]);
		double Dphis = Pi2/((double) 4*time_period[0]*N_2pi);
		double complex ExpP = exp(I*Dphi*((double) j_period - (time_period[0]+time_period[1])));
		double complex ExpPs = exp(I*Dphis*(0.5*((double) j_half) + ((double) j_time)));
		f_env = creal(ExpP*ExpPs);
	}else if(2*time_period[0]+time_period[1] <= j_period && j_period < time_period[2]){
		f_env = 0.;
	}else{fprintf(stderr,"f_env error!");	exit(EXIT_FAILURE);}
//	double inv_sqrnorm = 1./sqrt((double) 1*(time_period[0]+time_period[1])*N_2pi);
	double inv_norm = 1./((double) 1*(time_period[0]+time_period[1])*N_2pi);
	return(inv_norm*f_env*f_env);
}
//======================================================================
//  set_incident ---- set_Xi, set_tilXi    Last updated on Sep 30, 2019 
//======================================================================
void set_incident(const int iB, const int iT, const int iI, const char f_prefix[], const int L[4], const int inci_mode, double complex Hinc[2*L[0]*L[1]], double complex Einc[2*L[0]*L[1]])
{
	char Xi_file[BUFSIZE], mode_file[BUFSIZE];
	if(abs(iB-iI) < abs(iT-iI)){
		snprintf(Xi_file,sizeof(Xi_file),"%s_bXi.dat", f_prefix);
		snprintf(mode_file,sizeof(mode_file),"%s_bMode.dat", f_prefix);
	}else{
		snprintf(Xi_file,sizeof(Xi_file),"%s_tXi.dat", f_prefix);
		snprintf(mode_file,sizeof(mode_file),"%s_tMode.dat", f_prefix);
	}
	int NA = 2*L[0]*L[1];//	assert(Noffset <= NA);
	double complex *Xi = calloc(2*NA*2*NA,sizeof(double complex));
	int Nimag = get_RealNum(mode_file, L);
	set_Xi(Xi_file, NA, Nimag, Xi);
	double complex *tilXi = calloc(2*NA*2*NA,sizeof(double complex));
	set_tilXi(Xi_file, NA, Nimag, tilXi);
	double tilXiXi = check_Xi(2*NA, tilXi, Xi);
	if(tilXiXi > 1.e-15){fprintf(stderr,"tilXi^+ Xi = %.5e\n", tilXiXi);}
	SAFEFREE(tilXi);
	int column_num = 0;
	if(inci_mode > Nimag){fprintf(stderr,"inci_mode = %d > %d\n", inci_mode, Nimag);	exit(EXIT_FAILURE);}
	if(inci_mode > 0){column_num = inci_mode - 1;}
	else if(inci_mode < 0){column_num = (-inci_mode) - 1 + NA;}
	else{fprintf(stderr,"inci_mode = %d!\n", inci_mode);	exit(EXIT_FAILURE);}
	for(int irow = NA ; irow < 2*NA ; irow++){	if(cimag(Xi[irow + (2*NA)*column_num]) != 0.){fprintf(stderr,"incident wave is not real number!\n");	exit(EXIT_FAILURE);}}
	for(int irow = 0 ; irow < L[0]*L[1] ; irow++){
		Hinc[0 + 2*irow] =  Xi[irow + 0*L[0]*L[1] + (2*NA)*column_num];//Hinc[irow] = Xi[irow + (2*NA)*column_num];
		Hinc[1 + 2*irow] =  Xi[irow + 1*L[0]*L[1] + (2*NA)*column_num];
		Einc[1 + 2*irow] = -Xi[irow + 2*L[0]*L[1] + (2*NA)*column_num];
		Einc[0 + 2*irow] =  Xi[irow + 3*L[0]*L[1] + (2*NA)*column_num];
	}
	SAFEFREE(Xi);
}
//======================================================================
//  set_tau                               Last updated on Sep 24, 2019 
//======================================================================
void set_tau(const double omega, const double min_Dt, int *N_2pi, double *invtau_R, double *invtau_I)
{
	*N_2pi = 0;
	while(omega*min_Dt*(*N_2pi) <= Pi2){(*N_2pi)++;}//{N_2pi += 2;}
	double Dphi = Pi2/((double) *N_2pi);
	double Dt = Dphi/omega;
	*invtau_R = (0.5*Dphi)/(Dt*sin(0.5*Dphi));
	*invtau_I = (0.5*omega)/(cos(0.5*Dphi));
	fprintf(stderr,"omega = %.5e, min_Dt = %.5e, N_2pi = %d, Dt = %.5e, tau_R = %.5e, tau_I = %.5e\n", omega, min_Dt, *N_2pi, Dt, 1./(*invtau_R), 1./(*invtau_I));
}
//======================================================================
//  check_Xi                             Last updated on Sep 23, 2019 
//======================================================================
extern void zgemm_(char *transa, char *transb, const int *m, const int *n, const int *k,
	const double complex *alpha, const double complex *A, const int *ldA, 
	const double complex *B, const int *ldB,
	const double complex *beta , double complex *C, const int *ldC);
extern double dznrm2_(const int *n, const double complex *x, const int *incx);
double check_Xi(const int NS, const double complex tilXi[NS*NS], const double complex Xi[NS*NS])
{
	double complex *C = calloc(NS*NS,sizeof(double complex));
	{
		int m = NS;
		int n = NS;
		int k = NS;
		double complex alpha = 1.;
		double complex beta = 0.;
		zgemm_("C", "N", &m, &n, &k, &alpha, &tilXi[0], &m, &Xi[0], &k, &beta, &C[0], &m);
	}
	double norm;
	for(int row = 0 ; row < NS ; row++){C[row + NS*row] -= 1.;}
	{
		int NS2 = NS*NS;
		int incx = 1;
		norm = dznrm2_(&NS2,&C[0],&incx);
	}
	SAFEFREE(C);
	return(norm);
}
//======================================================================
//  set_Xi                           Last updated on Sep 23, 2019   
//======================================================================
void set_Xi(const char Xi_file[], const int NA, const int Nimag, double complex Xi[2*NA*2*NA])
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
				if(jcolumn < Nimag){
					Xi[row + (2*NA)*jcolumn] = dummy0 + I*dummy1;
					Xi[row + (2*NA)*(jcolumn+NA)] = dummy2 + I*dummy3;
				}else{
					Xi[row + (2*NA)*jcolumn] = dummy0 + I*dummy1;
					Xi[row + (2*NA)*(jcolumn+NA)] = dummy2 + I*dummy3;
				}
/*				if(jcolumn < Nimag){
					tilXi[(row+NA)%(2*NA) + (2*NA)*jcolumn] = dummy0 + I*dummy1;
					tilXi[(row+NA)%(2*NA) + (2*NA)*(jcolumn+NA)] = -(dummy2 + I*dummy3);
				}else{
					tilPhi[(row+NA)%(2*NA) + (2*NA)*(jcolumn+NA)] = dummy0 + I*dummy1;
					tilPhi[(row+NA)%(2*NA) + (2*NA)*jcolumn] = dummy2 + I*dummy3;
				}*/
			}else{fprintf(stderr,"From %s, data cannot be read at column = %d, row = %d!\n", Xi_file, jcolumn, row);	fprintf(stderr,"%s\n", buf);	exit(EXIT_FAILURE);}
		}
	}
	if(fclose(fp_i) != 0) {	fprintf(stderr,"fclose error in set_Phi!\n");	exit(EXIT_FAILURE);}
}
//======================================================================
//  set_tilXi                           Last updated on Nov 05, 2018   
//======================================================================
void set_tilXi(const char Xi_file[], const int NA, const int Nimag, double complex tilXi[2*NA*2*NA])
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
				if(jcolumn < Nimag){
					tilXi[(row+NA)%(2*NA) + (2*NA)*jcolumn] = dummy0 + I*dummy1;
					tilXi[(row+NA)%(2*NA) + (2*NA)*(jcolumn+NA)] = -(dummy2 + I*dummy3);
				}else{
					tilXi[(row+NA)%(2*NA) + (2*NA)*(jcolumn+NA)] = dummy0 + I*dummy1;
					tilXi[(row+NA)%(2*NA) + (2*NA)*jcolumn] = dummy2 + I*dummy3;
				}
			}else{fprintf(stderr,"From %s, data cannot be read at column = %d, row = %d!\n", Xi_file, jcolumn, row);	fprintf(stderr,"%s\n", buf);	exit(EXIT_FAILURE);}
		}
	}
	if(fclose(fp_i) != 0) {	fprintf(stderr,"fclose error in set_tilPhi!\n");	exit(EXIT_FAILURE);}
}
//======================================================================
//  get_RealNum                Last updated on Sep 26, 2019   
//======================================================================
int get_RealNum(const char file_name[], const int L[4])
{
		FILE *fp_i;
		fp_i = fopen(file_name,"r");
		char buf[BUFSIZE];	// buffer for fgets
		int j_count = -1;
		int dummy = 0;
		int N_real = 0;
		/*while(fgets(buf, sizeof( buf ), fp_i) != NULL && j_count < 0) { 
			if(strncmp(buf, "# Eigen-value", 13) == 0 && sscanf(buf,"%*[^=] %*[=] %d %*[^=] %*[=] %d", &N_real, &dummy) == 2){j_count++;}
		}*/
		int theta_count = 0;
		int theta_num0, theta_num1;
		double r_dummy0, i_dummy0, r_dummy1, i_dummy1;
		double complex *theta = calloc(4*L[0]*L[1],sizeof(double complex));
		while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
			if(strncmp(buf, "# Eigen-value", 13) == 0 && sscanf(buf,"%*[^=] %*[=] %d %*[^=] %*[=] %d", &N_real, &dummy) == 2){j_count++;}
			else if(strncmp(buf, "theta[", 6) == 0 && sscanf(buf,"%*[^[] %*[[] %d %*[^(] %*[(] %lf %*[,] %lf %*[^[] %*[[] %d %*[^(] %*[(] %lf %*[,] %lf", 
				&theta_num0, &r_dummy0, &i_dummy0, &theta_num1, &r_dummy1, &i_dummy1) == 6){
				theta[theta_num0] = r_dummy0 + I*i_dummy0;	theta[theta_num1] = r_dummy1 + I*i_dummy1;
				//if(theta_count == 0){fprintf(stderr,"%d, %.5e, %.5e, %d, %.5e, %.5e\n", theta_num0, r_dummy0, i_dummy0, theta_num1, r_dummy1, i_dummy1);}
				if(theta_count != theta_num0){fprintf(stderr,"theta cannot be read in get_RealNum05!\n");	exit(EXIT_FAILURE);}
				theta_count++;
			}
		}
		SAFEFREE(theta);
		if(j_count != 0 || N_real + dummy != 2*L[0]*L[1] || theta_count != 2*L[0]*L[1]) {fprintf(stderr,"RealNum05 cannot be read in get_RealNum05!\n");	exit(EXIT_FAILURE);}
		if(fclose(fp_i) != 0) {	fprintf(stderr,"fclose error in get_RealNum05!\n");	exit(EXIT_FAILURE);}
		return(N_real);
}
//======================================================================
//  input_data                              Last updated on Sep 24, 2019.  
//======================================================================
void input_data(const char data_name[], const int L[4], double *min_sqDt, double complex mu[3*L[0]*L[1]], double complex epsi[3*L[0]*L[1]])
{
	FILE *fp_i;
	fp_i = fopen(data_name,"r");
	if (fp_i == NULL){	fprintf(stderr,"open error for %s\n", data_name);	exit(EXIT_FAILURE);}
	char buf[BUFSIZE];
	int j0 = 0;	int j1 = 0;	int j2 = 0;	int j_count = -1;	int j_count2 = 0;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
		if(strncmp(buf,"# l <",5) == 0 && sscanf(buf,"%*[^<] %*[<] %*[^<] %*[<] %*[^<] %*[<] %*[^<] %*[<] %*[^=] %*[=] %*[^=] %*[=] %*[^=] %*[=] %lf", min_sqDt) == 1){j_count++;}//0
		int jl, jm, jn;
		double Rdummy00, Rdummy11, Rdummy22, Idummy00, Idummy11, Idummy22, fdummy00, fdummy11, fdummy22;
		if(sscanf(buf,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", &jl, &jm, &jn, &Rdummy00, &Idummy00, &Rdummy11, &Idummy11, &Rdummy22, &Idummy22, &fdummy00, &fdummy11, &fdummy22) == 12){
			if(j0 == jl && j1 == jm && j2 == jn){
				if(j2 == 0){
					  mu[0 + 3*j0 + 3*L[0]*j1] = (Rdummy00 + I*Idummy00)*fdummy00;
					  mu[1 + 3*j0 + 3*L[0]*j1] = (Rdummy11 + I*Idummy11)*fdummy11;
					epsi[2 + 3*j0 + 3*L[0]*j1] = (Rdummy22 + I*Idummy22)*fdummy22;
				}else if(j2 == 1){
					epsi[1 + 3*j0 + 3*L[0]*j1] = (Rdummy00 + I*Idummy00)*fdummy00;
					epsi[0 + 3*j0 + 3*L[0]*j1] = (Rdummy11 + I*Idummy11)*fdummy11;
					  mu[2 + 3*j0 + 3*L[0]*j1] = (Rdummy22 + I*Idummy22)*fdummy22;
				}else{fprintf(stderr,"j2 = %d in input_data2!\n", j2);	exit(EXIT_FAILURE);}
				j2++;
				if(j2 == 2){
					j2 = 0;	j0++;
					if(j0 == L[0]){	j0 = 0;	j1++;}
				}
				j_count2++;
			}else{fprintf(stderr,"Data cannot be read in input_data2!\n");	fprintf(stderr,"%s\n", buf);	exit(EXIT_FAILURE);}
		}
	}
	if(fclose(fp_i) != 0){	fprintf(stderr,"fclose error in input_data!\n");	exit(EXIT_FAILURE);}
	else if(j_count < 0){	fprintf(stderr,"min_sqDt can not be read in input_data!\n");	exit(EXIT_FAILURE);}	
	else if(j_count2 != 2*L[0]*L[1]){	fprintf(stderr,"j_count2 = %d != %d in input_data!\n", j_count2, 2*L[0]*L[1]);	exit(EXIT_FAILURE);}	
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
//======================================================================
//  matdata_file                     Last updated on Oct 22, 2018   
//======================================================================
void matdata_file(const char f_prefix[], const int L[4], const int istep, char data_name[])
{
	if(istep < L[3]){snprintf(data_name,BUFSIZE,"%s_MedB%d.dat", f_prefix, istep);}
	else if(L[3] <= istep && istep < L[2] - L[3]) {snprintf(data_name,BUFSIZE,"%s_Med%d.dat", f_prefix, istep - L[3]);}
	else if(L[2] - L[3] <= istep){snprintf(data_name,BUFSIZE,"%s_MedT%d.dat", f_prefix, istep - (L[2] - L[3]));}

}
//======================================================================
//  input_Lf                              Last updated on Sep 17, 2019.  
//======================================================================
void input_Lf(const char Yee[], int L[4], double *omega)
{
	char input_name[BUFSIZE], buf[BUFSIZE];
	snprintf(input_name,sizeof(input_name),"%s_Med0.dat", Yee);
	FILE *fp_i;
	fp_i = fopen(input_name,"r");
	if (fp_i == NULL){	fprintf(stderr,"open error for %s\n", input_name);	exit(EXIT_FAILURE);}
	int j_count = -3;
	double wavelength, k_0;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL && j_count < 0) { 
		if(sscanf(buf,"%*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %*[^=] %*[=] %d", &L[0], &L[1], &L[2], &L[3]) == 4){j_count++;}//-2
		else if(strncmp(buf,"# wavelength",12) == 0 && sscanf(buf,"%*[^=] %*[=] %lf", &wavelength) == 1){j_count++;}//-1
		else if(strncmp(buf,"# k_0",5) == 0 && sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %lf", &k_0) == 1){j_count++;}//0
	}
	if(fclose(fp_i) != 0){	fprintf(stderr,"fclose error after reading the Med0 file!\n");	exit(EXIT_FAILURE);}
	else if(j_count < 0){	fprintf(stderr,"L,M,N,N_outer can not be read at Scatterer @ input_file0!\n");	exit(EXIT_FAILURE);}
	else{	fprintf(stderr,"wavelength = %.5e, k_0 = %.5e\n", wavelength, k_0); *omega = Pi2/(wavelength*k_0);}
	
	snprintf(input_name,sizeof(input_name),"%s_bMode.dat", Yee);
	fp_i = fopen(input_name,"r");
	double omega_dummy;
	j_count = -1;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL && j_count < 0) { 
		if(strncmp(buf, "# Omega", 7) == 0 && sscanf(buf,"%*[^=] %*[=] %lf", &omega_dummy) == 1){j_count++;}
	}
	if(fclose(fp_i) != 0){	fprintf(stderr,"fclose error after reading the bMode file!\n");	exit(EXIT_FAILURE);}
	else if(j_count != 0){	fprintf(stderr,"Omega can not be read in input_Lf!\n");	exit(EXIT_FAILURE);}
	else if(*omega != omega_dummy) {fprintf(stderr,"Omega = %.5e != %.5e in input_Lf!\n", *omega, omega_dummy);	exit(EXIT_FAILURE);}
}
//======================================================================
//  input_file0 ---- rm_space, rm_comma                                 
//                                       Last updated on Oct 01, 2019   
//======================================================================
void rm_space( char *A );
void rm_comma( char *A );
void input_file0(FILE *fp_i0, int BoundC[], char f_prefix[], char p_target[], int time_period[3], char inci_point[], int *inci_mode, double *Dt_factor)
{
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = -4;
	while(fgets(buf, sizeof( buf ), fp_i0) != NULL && j_count < 0) { 
		rm_space(buf);
		if(strncmp(buf, "Prefix", 6) == 0 && sscanf(buf,"%*[^=] %*[=] %s", f_prefix) == 1){
			rm_comma(f_prefix);
			j_count++;//-3
		}else if(strncmp(buf, "BoundaryCondition", 17) == 0 && sscanf(buf,"%*[^=] %*[=] %d %*[^=] %*[=] %d %*[^=] %*[=] %d", &BoundC[0], &BoundC[1], &BoundC[2]) == 3) {//revised on 20191001
			j_count++;//-2
		}else if(strncmp(buf, "ModeP", 5) == 0 && sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %*[^=] %*[=] %s", p_target) == 1) {
			rm_comma(p_target);
			if(sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %*[^=] %*[=] %*[^=] %*[=] %s", &p_target[BUFSIZE]) != 1){sscanf(p_target,"%s", &p_target[BUFSIZE]);}
			else{rm_comma(&p_target[BUFSIZE]);}
			j_count++;//-1
		}else if(strncmp(buf, "FDTD", 4) == 0 && sscanf(buf,"%*[^=] %*[=] %d %*[^=] %*[=] %d %*[^=] %*[=] %d %*[^=] %*[=] %s %*[^=] %*[=] %d %*[^=] %*[=] %lf", &time_period[0], &time_period[1], &time_period[2], inci_point, inci_mode, Dt_factor) == 6) {
			rm_comma(inci_point);
			j_count++;//0
		}
	}
	if(j_count != 0) {	fprintf(stderr,"4 control commands cannot be read in input_file0!");	exit(EXIT_FAILURE);}
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
//  This program removes sigma from mu and eqsi.
//                                        Last updated on Oct 17, 2018  
//======================================================================
void remove_sigma(const int L[4], double complex mu[3*L[0]*L[1]*L[2]], double complex epsi[3*L[0]*L[1]*L[2]], const double sigma[2*(L[0]+L[1]+L[2])])
{
	double complex s[2*3];//	double complex s0[2];	double complex s1[2];	double complex s2[2];
	for(int j2 = 0 ; j2 < L[2] ; j2++){
		for(int k = 0 ; k < 2 ; k++){	s[k+2*2] = 1. + I*sigma[2*L[0]+2*L[1]+k+2*j2];}
		for(int j1 = 0 ; j1 < L[1] ; j1++){
			for(int k = 0 ; k < 2 ; k++){	s[k+2*1] = 1. + I*sigma[2*L[0]+k+2*j1];}
			for(int j0 = 0 ; j0 < L[0] ; j0++){
				for(int k = 0 ; k < 2 ; k++){	s[k+2*0] = 1. + I*sigma[k+2*j0];}
				double complex cnorm[3];
				for(int js = 0 ; js < 3 ; js++){
					int jsp1 = (js+1)%3;	int jsp2 = (js+2)%3;
					cnorm[js] = s[0+2*jsp1]*s[0+2*jsp2]*conj(s[1+2*js])/(creal(s[1+2*js])*creal(s[1+2*js]) + cimag(s[1+2*js])*cimag(s[1+2*js]));//s1[0]*s2[0]/s0[1];
				}
				for(int js = 0 ; js < 3 ; js++){
					mu[js + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] *= conj(cnorm[js])/(creal(cnorm[js])*creal(cnorm[js]) + cimag(cnorm[js])*cimag(cnorm[js]));// /= cnorm[js];
					if(cimag(cnorm[js]) != 0.){if(fabs(carg(mu[js + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2])) > ERROR_MIN ){	fprintf(stderr,"mu error in remove_sigma\n");	exit(EXIT_FAILURE);}}
					if(cimag(cnorm[js]) != 0.){mu[js + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] = creal(mu[js + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2]);}
				}
				for(int js = 0 ; js < 3 ; js++){
					int jsp1 = (js+1)%3;	int jsp2 = (js+2)%3;
					cnorm[js] = s[1+2*jsp1]*s[1+2*jsp2]*conj(s[0+2*js])/(creal(s[0+2*js])*creal(s[0+2*js]) + cimag(s[0+2*js])*cimag(s[0+2*js]));//s1[1]*s2[1]/s0[0];
				}
				for(int js = 0 ; js < 3 ; js++){
					epsi[js + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] *= conj(cnorm[js])/(creal(cnorm[js])*creal(cnorm[js]) + cimag(cnorm[js])*cimag(cnorm[js]));// /= cnorm[js];
					if(cimag(cnorm[js]) != 0.){if(fabs(carg(epsi[js + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2])) > ERROR_MIN ){	fprintf(stderr,"epsi error in remove_sigma\n");	exit(EXIT_FAILURE);}}
					if(cimag(cnorm[js]) != 0.){epsi[js + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] = creal(epsi[js + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2]);}
				}
			}
		}
	}
}
//======================================================================
//  This program creats sigma from mu and eqsi.
//                                        Last updated on Oct 17, 2018  
//======================================================================
void set_sigma(const int L[4], const double complex mu[3*L[0]*L[1]*L[2]], const double complex epsi[3*L[0]*L[1]*L[2]], double sigma[2*(L[0]+L[1]+L[2])])
{
	{	int j2 = L[2]/2 ;
		{	int j1 = L[1]/2;
			for(int j0 = 0 ; j0 < L[0] ; j0++){
				double arg11 = carg( mu[1 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
				if(arg11 != 0.){
					double arg22 = carg( mu[2 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
					double arg00 = -carg( epsi[0 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
					if(fabs(arg00 - arg11) < ERROR_MIN && fabs(arg11 - arg22) < ERROR_MIN &&fabs(arg22 - arg00) < ERROR_MIN){
						sigma[0 + 2*j0] = tan(arg00); 	fprintf(stderr,"l0 = %d.0, sigma = %.5e\n", j0, sigma[0 + 2*j0]);
					}
				}
				arg11 = carg( epsi[1 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
				if(arg11 != 0.){
					double arg22 = carg( epsi[2 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
					double arg00 = -carg( mu[0 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
					if(fabs(arg00 - arg11) < ERROR_MIN && fabs(arg11 - arg22) < ERROR_MIN &&fabs(arg22 - arg00) < ERROR_MIN){
						sigma[1 + 2*j0] = tan(arg00); 	fprintf(stderr,"l0 = %d.5, sigma = %.5e\n", j0, sigma[1 + 2*j0]);
					}
				}
			}
		}
	}
	{	int j2 = L[2]/2 ;
		for(int j1 = 0 ; j1 < L[1] ; j1++){
			{	int j0 = L[0]/2;
				double arg22 = carg( mu[2 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
				if(arg22 != 0.){
					double arg00 = carg( mu[0 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
					double arg11 = -carg( epsi[1 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
					if(fabs(arg00 - arg11) < ERROR_MIN && fabs(arg11 - arg22) < ERROR_MIN &&fabs(arg22 - arg00) < ERROR_MIN){
						sigma[2*L[0] + 0 + 2*j1] = tan(arg00);	fprintf(stderr,"l1 = %d.0, sigma = %.5e\n", j1, sigma[2*L[0] + 0 + 2*j1]);
					}
				}
				arg22 = carg( epsi[2 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
				if(arg22 != 0.){
					double arg00 = carg( epsi[0 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
					double arg11 = -carg( mu[1 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
					if(fabs(arg00 - arg11) < ERROR_MIN && fabs(arg11 - arg22) < ERROR_MIN &&fabs(arg22 - arg00) < ERROR_MIN){
						sigma[2*L[0] + 1 + 2*j1] = tan(arg00);	fprintf(stderr,"l1 = %d.5, sigma = %.5e\n", j1, sigma[2*L[0] + 1 + 2*j1]);
					}
				}
			}
		}
	}
	for(int j2 = 0 ; j2 < L[2] ; j2++){
		{	int j1 = L[1]/2;
			{	int j0 = L[0]/2;
				double arg00 = carg( mu[0 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
				if(arg00 != 0.){
					double arg11 = carg( mu[1 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
					double arg22 = -carg( epsi[2 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
					if(fabs(arg00 - arg11) < ERROR_MIN && fabs(arg11 - arg22) < ERROR_MIN &&fabs(arg22 - arg00) < ERROR_MIN){
						sigma[2*L[0] + 2*L[1] + 0 + 2*j2] = tan(arg00);	fprintf(stderr,"l2 = %d.0, sigma = %.5e\n", j2, sigma[2*L[0] + 2*L[1] + 0 + 2*j2]);
					}
				}
				arg00 = carg( epsi[0 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
				if(arg00 != 0.){
					double arg11 = carg( epsi[1 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
					double arg22 = -carg( mu[2 + 3*j0 + 3*L[0]*j1 + 3*L[0]*L[1]*j2] );
					if(fabs(arg00 - arg11) < ERROR_MIN && fabs(arg11 - arg22) < ERROR_MIN &&fabs(arg22 - arg00) < ERROR_MIN){
						sigma[2*L[0] + 2*L[1] + 1 + 2*j2] = tan(arg00);	fprintf(stderr,"l2 = %d.5, sigma = %.5e\n", j2, sigma[2*L[0] + 2*L[1] + 1 + 2*j2]);
					}
				}
			}
		}
	}
}
