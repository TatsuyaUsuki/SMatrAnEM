//  Start of PMLcalc.h
#include "header_macro.h"
#include "constant.h"
void rename_file(const int l2, const char Yee[], const int L[4]);
void output_file(const int l2, const char Yee[], const int L[4], const double omega, const char *exec);
void input_parameter(const char Yee[], int L[4], double *omega);
void input_filename(FILE *fp_i1, char Yee[]);
//  End of PMLcalc.h
//======================================================================
//   main ---- input_filename, input_parameter              
//======================================================================
int main(int argc, char **argv)
{
	time_t timer_ini = time(0);	fprintf(stderr,"# Start time of setting correction factor = %s\n", ctime(&timer_ini));
	if(argc != 2) {fprintf(stderr,"error: number of files \n");	exit(EXIT_FAILURE);}
	else if(strncmp(argv[1], "-v", 2) == 0 || strcmp(argv[1], "--version") == 0 ) {
		fprintf(stderr,"The '%s' creates correction factor.\n", argv[0]);
		fprintf(stderr,"Version 19.09.19 is compiled at %s on %s.\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(stderr," Source code : '%s'\n Author      : Tatsuya Usuki\n URL         : http://www.smatran.org\n", __FILE__);
		fprintf(stderr," References  : Equations (6.6) and (H.3) in 'Formulation.pdf' on Aug 25, 2019\n");
		fprintf(stderr,"There is NO warranty.\n");
		exit(EXIT_SUCCESS);//normal end
	}
//------- begin reading file names -------
	FILE *fp_i1;
	fp_i1 = fopen(argv[1],"r");
	if (fp_i1 == NULL){	fprintf(stderr,"open error!: open input-file1!\n");	exit(EXIT_FAILURE);}
	fprintf(stderr,"The 1st input file: %s\n",argv[1]);
	char Yee[BUFSIZE];
	input_filename(fp_i1, Yee);
	fprintf(stderr,"Prefix: %s\n", Yee);
	if(fclose(fp_i1) != 0) {	fprintf(stderr,"fclose error after input_file!\n");	exit(EXIT_FAILURE);}
	int L[4];//L[2] = bottom + scatterer + top, L[3] = outer number
	double omega = 0.;
	input_parameter(Yee, L, &omega);
	fprintf(stderr,"L0 = %d, L1 = %d, Lbst = %d, Louter = %d \n", L[0], L[1], L[2], L[3]);
	fprintf(stderr,"omega = k/k_0 = %.5e \n", omega);
//-------  end reading file names  -------
	for(int l2 = 0 ; l2 < L[2] ; l2++){	output_file(l2, Yee, L, omega, argv[0]);}
	for(int l2 = 0 ; l2 < L[2] ; l2++){	rename_file(l2, Yee, L);}
}
//======================================================================
//  rename_file                            Last updated on Sep 19, 2019 
//======================================================================
void matdata_file(const char f_prefix[], const char add_name[], const int Ltot, const int Lout, const int istep, char data_name[]);
void rename_file(const int l2, const char Yee[], const int L[4])
{
	char data_name[BUFSIZE], temp_name[BUFSIZE];
	matdata_file(Yee, "_Med", L[2], L[3], l2, data_name);
	matdata_file(Yee, "_temp", L[2], L[3], l2, temp_name);
	if(remove(data_name)==0){rename(temp_name, data_name);}
}
//======================================================================
//  output_file ---- set_rep_rmu, matdata_file,
//                    set_eta ---- set_epmu, set_DX2
//                                         Last updated on Sep 04, 2019 
//======================================================================
void set_eta(const int L[4], const double omega, const double complex cmed[], const double fx[], const double r_ep[], const double r_mu[], double eta[]);
void set_rep_rmu(const int l2, const char Yee[], const int L[4], double r_ep[], double r_mu[]);
//void matdata_file(const char f_prefix[], const char add_name[], const int Ltot, const int Lout, const int istep, char data_name[]);
void output_file(const int l2, const char Yee[], const int L[4], const double omega, const char *exec)
{
	char data_name[BUFSIZE], temp_name[BUFSIZE];
	matdata_file(Yee, "_Med", L[2], L[3], l2, data_name);
	fprintf(stderr,"%s\n", data_name);
	matdata_file(Yee, "_temp", L[2], L[3], l2, temp_name);
	
	double complex *cmed = calloc(3*2*L[0]*L[1],sizeof(double complex));
	double *fx = calloc(3*2*L[0]*L[1],sizeof(double));
	for(int l1 = 0 ; l1 < L[1] ; l1++){
		for(int l0 = 0 ; l0 < L[0] ; l0++){
			for(int j2 = 0 ; j2 < 2 ; j2++){
				for(int je = 0 ; je < 3; je++){
					cmed[je + 3*j2 + 3*2*l0 + 3*2*L[0]*l1] = 0.+I*0.;
					fx[je + 3*j2 + 3*2*l0 + 3*2*L[0]*l1] = 0.;
				}
			}
		}
	}
	FILE *fp_r, *fp_w;
	fp_r = fopen(data_name,"r");
	fp_w = fopen(temp_name,"w");
	if (fp_r == NULL || fp_w == NULL){	fprintf(stderr,"open error for %s or %s in output_file\n", data_name, temp_name);	exit(EXIT_FAILURE);}
	
	char buf[BUFSIZE];
	while(fgets(buf, BUFSIZE, fp_r) != NULL){
		if(strncmp(buf,"# Created",9) == 0){	fprintf(fp_w,"%s", buf);
		}else if(strncmp(buf,"# l <",5) == 0){
			time_t timer_f = time(0);
			fprintf(fp_w,"%s# Correction factor 'eta' was added by ' %s ' on %s", buf, exec, ctime(&timer_f));
		}else if(strncmp(buf,"# l m 0",7) == 0){
			snprintf(buf, sizeof(buf), "# l m 0 mu_00 mu_11 ep_22 k0*fx_00/(1+eta_100) k0*fx_11/(1+eta_010) k0*fx_22/(1+eta_110) : 'fxjj' is defined as (duj+1dxj+1)*(hj+1)*(duj+2dxj+2)*(hj+2)/(dujdxj*hj)\n");
			fprintf(fp_w,"%s", buf);
		}else if(strncmp(buf,"# l m 1",7) == 0){
			snprintf(buf, sizeof(buf), "# l m 1 ep_11 ep_00 mu_22 k0*fx_11/(1+eta_101) k0*fx_00/(1+eta_011) k0*fx_22/(1+eta_001)\n");
			fprintf(fp_w,"%s", buf);
		}else{
			int l0, l1, j2;
			double re0, im0, re1, im1, re2, im2;
			double fx00, fx11, fx22;
			if(sscanf(buf,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", &l0, &l1, &j2, &re0, &im0, &re1, &im1, &re2, &im2, &fx00, &fx11, &fx22) == 12){
				if(j2 < 0 || j2 >= 2){	fprintf(stderr,"j2 = %d can not be read in output_file\n",j2);	exit(EXIT_FAILURE);}
				else if(l0 < 0 || l0 >= L[0]){	fprintf(stderr,"l0 = %d can not be read in output_file\n",l0);	exit(EXIT_FAILURE);}
				else if(l1 < 0 || l1 >= L[1]){	fprintf(stderr,"l1 = %d can not be read in output_file\n",l1);	exit(EXIT_FAILURE);}
				else{
					cmed[0 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1] = re0 + I*im0;
					cmed[1 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1] = re1 + I*im1;
					cmed[2 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1] = re2 + I*im2;
					fx[0 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1] = fx00;// fx[jdata<3][j2<2][l0<L[0]][l1<L[1]]
					fx[1 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1] = fx11;// f00(1/2, 0 , 0 )[0][0][l0][l1]; f11( 0 ,1/2, 0 )[1][0][l0][l1]; f22(1/2,1/2, 0 )[2][0][l0][l1]; 
					fx[2 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1] = fx22;// f11(1/2, 0 ,1/2)[0][1][l0][l1]; f00( 0 ,1/2,1/2)[1][1][l0][l1]; f22( 0 , 0 ,1/2)[2][1][l0][l1];
				}
			}else{
				if(strncmp(buf,"#",1) == 0){fprintf(fp_w,"%s", buf);}
				else{	fprintf(stderr,"buf error after fget in output_file!\n");	exit(EXIT_FAILURE);}
			}
		}
	}
	if(fclose(fp_r) != 0 ){	fprintf(stderr,"fclose error after reading a data file!\n");	exit(EXIT_FAILURE);}
	{
		double *eta = calloc(3*2*L[0]*L[1],sizeof(double));	for(int jtemp = 0 ; jtemp < 3*2*L[0]*L[1] ; jtemp++){eta[jtemp] = 0.;}
		{
			double *r_ep = calloc(2*2*L[0]*L[1],sizeof(double));
			double *r_mu = calloc(2*2*L[0]*L[1],sizeof(double));
			set_rep_rmu(l2, Yee, L, r_ep, r_mu);
			set_eta(L, omega, cmed, fx, r_ep, r_mu, eta);
			SAFEFREE(r_ep);	SAFEFREE(r_mu);
		}
		for(int l1 = 0 ; l1 < L[1] ; l1++){
			for(int l0 = 0 ; l0 < L[0] ; l0++){
				for(int j2 = 0 ; j2 < 2 ; j2++){
					fprintf(fp_w,"%d %d %d %.20e %.20e %.20e %.20e %.20e %.20e %.20e %.20e %.20e\n", l0, l1, j2, 
						creal(cmed[0 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1]), cimag(cmed[0 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1]), 
						creal(cmed[1 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1]), cimag(cmed[1 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1]), 
						creal(cmed[2 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1]), cimag(cmed[2 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1]), 
						fx[0 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1]/(1. + eta[0 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1]), // fx[jdata<3][j2<2][l0<L[0]][l1<L[1]]
						fx[1 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1]/(1. + eta[1 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1]), // f00(1/2, 0 , 0 )[0][0][l0][l1]; f11( 0 ,1/2, 0 )[1][0][l0][l1]; f22(1/2,1/2, 0 )[2][0][l0][l1]; 
						fx[2 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1]/(1. + eta[2 + 3*j2 + 3*2*l0 + 3*2*L[0]*l1])  // f11(1/2, 0 ,1/2)[0][1][l0][l1]; f00( 0 ,1/2,1/2)[1][1][l0][l1]; f22( 0 , 0 ,1/2)[2][1][l0][l1];
					);
				}
			}
		}
		SAFEFREE(eta);
	}
	SAFEFREE(cmed);	SAFEFREE(fx);
	if(fclose(fp_w) != 0){	fprintf(stderr,"fclose error after reading or writing a data file in output_file!\n");	exit(EXIT_FAILURE);}
}
//======================================================================
//  set_rep_rmu                            Last updated on Sep 19, 2019 
//======================================================================
void set_rep_rmu(const int l2, const char Yee[], const int L[4], double r_ep[], double r_mu[])
{
	if(l2 < 0 || L[2]-1 < l2){	fprintf(stderr,"l2 = %d in set_epmu!\n", l2);	exit(EXIT_FAILURE);}//revised on 20190919
	FILE *fp_r;
	char data_name[BUFSIZE];
	for(int jstep = 0 ; jstep < 2 ; jstep++){
		{
			int l2_shift = l2 + 2*jstep - 1;
			if(l2_shift == -1 || l2_shift == L[2]){	l2_shift = l2;}
			matdata_file(Yee, "_Med", L[2], L[3], l2_shift, data_name);
		}//revised on 20190919
		fp_r = fopen(data_name,"r");
		if (fp_r == NULL){	fprintf(stderr,"open error for %s in set_epmu\n", data_name);	exit(EXIT_FAILURE);}
	
		char buf[BUFSIZE];
		while(fgets(buf, BUFSIZE, fp_r) != NULL){
			int l0, l1, j2;
			double re0, im0, re1, im1, re2, im2;
			double fx00, fx11, fx22;
			if(sscanf(buf,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", &l0, &l1, &j2, &re0, &im0, &re1, &im1, &re2, &im2, &fx00, &fx11, &fx22) == 12){
				if(j2 < 0 || j2 >= 2){	fprintf(stderr,"j2 = %d can not be read in output_file\n",j2);	exit(EXIT_FAILURE);}
				else if(l0 < 0 || l0 >= L[0]){	fprintf(stderr,"l0 = %d can not be read in output_file\n",l0);	exit(EXIT_FAILURE);}
				else if(l1 < 0 || l1 >= L[1]){	fprintf(stderr,"l1 = %d can not be read in output_file\n",l1);	exit(EXIT_FAILURE);}
				else{
					if(jstep == 0 && j2 == 1){
						r_ep[0 + 2*0 + 2*2*l0 + 2*2*L[0]*l1] = re0;
						r_ep[1 + 2*0 + 2*2*l0 + 2*2*L[0]*l1] = re1;
						r_ep[0 + 2*1 + 2*2*l0 + 2*2*L[0]*l1] = fx00;
						r_ep[1 + 2*1 + 2*2*l0 + 2*2*L[0]*l1] = fx11;
					}else if(jstep == 1 && j2 == 0){
						r_mu[0 + 2*0 + 2*2*l0 + 2*2*L[0]*l1] = re0;
						r_mu[1 + 2*0 + 2*2*l0 + 2*2*L[0]*l1] = re1;
						r_mu[0 + 2*1 + 2*2*l0 + 2*2*L[0]*l1] = fx00;
						r_mu[1 + 2*1 + 2*2*l0 + 2*2*L[0]*l1] = fx11;
					}
				}
			}
		}
		if(fclose(fp_r) != 0 ){	fprintf(stderr,"fclose error after reading a data file in set_epmu!\n");	exit(EXIT_FAILURE);}
	}
}
//======================================================================
//  matdata_file                     Last updated on Aug 29, 2019   
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
//  set_eta ---- set_epmu, set_DX2         Last updated on Sep 03, 2019 
//======================================================================
void set_epmu(const int L[4], const int l0, const int l1, const double complex cmed[], const double r_ep[], const double r_mu[], double epmu[6]);
void set_DX2(const int L[4], const int l0, const int l1, const double fx[], const double r_ep[], const double r_mu[], double DX2[]);
void set_eta(const int L[4], const double omega, const double complex cmed[], const double fx[], const double r_ep[], const double r_mu[], double eta[])
{
	for(int l1 = 0 ; l1 < L[1] ; l1++){
		for(int l0 = 0 ; l0 < L[0] ; l0++){
			double DX2[6];
			set_DX2(L, l0, l1, fx, r_ep, r_mu, DX2);
			double epmu[6];
			set_epmu(L, l0, l1, cmed, r_ep, r_mu, epmu);
			for(int j6 = 0 ; j6 < 6; j6++){
				eta[j6 + 3*2*l0 + 3*2*L[0]*l1] = (omega*omega/120.)*epmu[j6]*DX2[j6];//	eta[j6 + 3*2*l0 + 3*2*L[0]*l1] = DX2[j6];
			}
		}
	}
}
//======================================================================
//  set_epmu                               Last updated on Sep 03, 2019 
//======================================================================
void set_epmu(const int L[4], const int l0, const int l1, const double complex cmed[], const double r_ep[], const double r_mu[], double epmu[6])
{
	int l1p = (l1 + 1 + L[1])%L[1];
	int l1m = (l1 - 1 + L[1])%L[1];
	int l0p = (l0 + 1 + L[0])%L[0];
	int l0m = (l0 - 1 + L[0])%L[0];
	for(int j6 = 0 ; j6 < 6; j6++){	epmu[j6] = creal(cmed[j6 + 3*2*l0 + 3*2*L[0]*l1]);}
	{ double ep;//xi - l = (1/2, 0 , 0 )
		ep =  creal(cmed[0 + 3*1 + 3*2*l0 + 3*2*L[0]*l1]) +        r_ep[( 0+0 ) + 2*2*l0 + 2*2*L[0]*l1]; //(1/2, 0 ,1/2) + (1/2, 0 ,-1/2)
		ep += creal(cmed[2 + 3*0 + 3*2*l0 + 3*2*L[0]*l1]) + creal(cmed[2 + 3*0 + 3*2*l0 + 3*2*L[0]*l1m]);//(1/2,1/2, 0 ) + (1/2,-1/2, 0 )
		ep *= 0.25;
		epmu[0] *= ep;
	}{double ep;//xi - l = ( 0 ,1/2, 0 )
		ep =  creal(cmed[1 + 3*1 + 3*2*l0 + 3*2*L[0]*l1]) +        r_ep[( 1+0 ) + 2*2*l0 + 2*2*L[0]*l1]; //( 0 ,1/2,1/2) + ( 0 ,1/2,-1/2)
		ep += creal(cmed[2 + 3*0 + 3*2*l0 + 3*2*L[0]*l1]) + creal(cmed[2 + 3*0 + 3*2*l0m + 3*2*L[0]*l1]);//(1/2,1/2, 0 ) + (-1/2,1/2, 0 )
		ep *= 0.25;
		epmu[1] *= ep;
	}{double mu;//xi - l = (1/2,1/2, 0 )
		mu =  creal(cmed[0 + 3*0 + 3*2*l0 + 3*2*L[0]*l1]) + creal(cmed[0 + 3*0 + 3*2*l0 + 3*2*L[0]*l1p]);//(1/2, 0 , 0 ) + (1/2, 1 , 0 )
		mu += creal(cmed[1 + 3*0 + 3*2*l0 + 3*2*L[0]*l1]) + creal(cmed[1 + 3*0 + 3*2*l0p + 3*2*L[0]*l1]);//( 0 ,1/2, 0 ) + ( 1 ,1/2, 0 )
		mu *= 0.25;
		epmu[2] *= mu;
	}
	{ double mu;//xi - l = (1/2, 0 ,1/2)
		mu =  creal(cmed[0 + 3*0 + 3*2*l0 + 3*2*L[0]*l1]) +        r_mu[( 0+0 ) + 2*2*l0 + 2*2*L[0]*l1]; //(1/2, 0 , 0 ) + (1/2, 0 , 1 )
		mu += creal(cmed[2 + 3*1 + 3*2*l0 + 3*2*L[0]*l1]) + creal(cmed[2 + 3*1 + 3*2*l0p + 3*2*L[0]*l1]);//( 0 , 0 ,1/2) + ( 1 , 0 ,1/2)
		mu *= 0.25;
		epmu[3] *= mu;
	}{double mu;//xi - l = ( 0 ,1/2,1/2)
		mu =  creal(cmed[1 + 3*0 + 3*2*l0 + 3*2*L[0]*l1]) +        r_mu[( 1+0 ) + 2*2*l0 + 2*2*L[0]*l1]; //( 0 ,1/2, 0 ) + ( 0 ,1/2, 1 )
		mu += creal(cmed[2 + 3*1 + 3*2*l0 + 3*2*L[0]*l1]) + creal(cmed[2 + 3*1 + 3*2*l0 + 3*2*L[0]*l1p]);//( 0 , 0 ,1/2) + ( 0 , 1 ,1/2)
		mu *= 0.25;
		epmu[4] *= mu;
	}{double ep;//xi - l = ( 0 , 0 ,1/2)
		ep =  creal(cmed[0 + 3*1 + 3*2*l0 + 3*2*L[0]*l1]) + creal(cmed[0 + 3*1 + 3*2*l0m + 3*2*L[0]*l1]);//(1/2, 0 ,1/2) + (-1/2, 0 ,1/2)
		ep += creal(cmed[1 + 3*1 + 3*2*l0 + 3*2*L[0]*l1]) + creal(cmed[1 + 3*1 + 3*2*l0 + 3*2*L[0]*l1m]);//( 0 ,1/2,1/2) + ( 0 ,-1/2,1/2)
		ep *= 0.25;
		epmu[5] *= ep;
	}
}
//======================================================================
//  set_DX2                                Last updated on Sep 04, 2019 
//======================================================================
void set_DX2(const int L[4], const int l0, const int l1, const double fx[], const double r_ep[], const double r_mu[], double DX2[6])
{
	int l1p = (l1 + 1 + L[1])%L[1];
	int l1m = (l1 - 1 + L[1])%L[1];
	int l0p = (l0 + 1 + L[0])%L[0];
	int l0m = (l0 - 1 + L[0])%L[0];
	double D0, D1, D2; // fx[jdata<3][j2<2][l0<L[0]][l1<L[1]]
	//xi - l = (1/2, 0 , 0 )
	D0 = fx[0 + 3*0 + 3*2*l0 + 3*2*L[0]*l1];//(1/2, 0 , 0 ) ; f00(1/2, 0 , 0 )[0][0][l0][l1];
	D1 = 0.5*(fx[0 + 3*1 + 3*2*l0 + 3*2*L[0]*l1] + r_ep[( 0+2 ) + 2*2*l0 + 2*2*L[0]*l1]); //(1/2, 0 ,1/2) + (1/2, 0 ,-1/2) ; f11(1/2, 0 ,1/2)[0][1][l0][l1];
	D2 = 0.5*(fx[2 + 3*0 + 3*2*l0 + 3*2*L[0]*l1] +   fx[2 + 3*0 + 3*2*l0 + 3*2*L[0]*l1m]);//(1/2,1/2, 0 ) + (1/2,-1/2, 0 ) ; f22(1/2,1/2, 0 )[2][0][l0][l1]; 
	DX2[0] = D1*D2 + D2*D0 + D0*D1;
	
	//xi - l = ( 0 ,1/2, 0 )
	D1 = fx[1 + 3*0 + 3*2*l0 + 3*2*L[0]*l1];//( 0 ,1/2, 0 ) ; f11( 0 ,1/2, 0 )[1][0][l0][l1];
	D2 = 0.5*(fx[2 + 3*0 + 3*2*l0 + 3*2*L[0]*l1] +   fx[2 + 3*0 + 3*2*l0m + 3*2*L[0]*l1]);//(1/2,1/2, 0 ) + (-1/2,1/2, 0 ) ; f22(1/2,1/2, 0 )[2][0][l0][l1];
	D0 = 0.5*(fx[1 + 3*1 + 3*2*l0 + 3*2*L[0]*l1] + r_ep[( 1+2 ) + 2*2*l0 + 2*2*L[0]*l1]); //( 0 ,1/2,1/2) + ( 0 ,1/2,-1/2) ; f00( 0 ,1/2,1/2)[1][1][l0][l1];
	DX2[1] = D1*D2 + D2*D0 + D0*D1;
	
	//xi - l = (1/2,1/2, 0 )
	D2 = fx[2 + 3*0 + 3*2*l0 + 3*2*L[0]*l1];//(1/2,1/2, 0 ) ; f22(1/2,1/2, 0 )[2][0][l0][l1];
	D0 = 0.5*(fx[0 + 3*0 + 3*2*l0 + 3*2*L[0]*l1] + fx[0 + 3*0 + 3*2*l0 + 3*2*L[0]*l1p]); //(1/2, 0 , 0 ) + (1/2, 1 , 0 ) ; f00(1/2, 0 , 0 )[0][0][l0][l1];
	D1 = 0.5*(fx[1 + 3*0 + 3*2*l0 + 3*2*L[0]*l1] + fx[1 + 3*0 + 3*2*l0p + 3*2*L[0]*l1]); //( 0 ,1/2, 0 ) + ( 1 ,1/2, 0 ) ; f11( 0 ,1/2, 0 )[1][0][l0][l1];
	DX2[2] = D1*D2 + D2*D0 + D0*D1;
	
	//xi - l = (1/2, 0 ,1/2) ; f11(1/2, 0 ,1/2)[0][1][l0][l1];
	D1 = fx[0 + 3*1 + 3*2*l0 + 3*2*L[0]*l1];//(1/2, 0 ,1/2) ; f11(1/2, 0 ,1/2)[0][1][l0][l1];
	D2 = 0.5*(fx[2 + 3*1 + 3*2*l0p + 3*2*L[0]*l1] +  fx[2 + 3*1 + 3*2*l0 + 3*2*L[0]*l1]);//( 1 , 0 ,1/2) + ( 0 , 0 ,1/2) ; f22( 0 , 0 ,1/2)[2][1][l0][l1];
	D0 = 0.5*(fx[0 + 3*0 + 3*2*l0 + 3*2*L[0]*l1] + r_mu[( 0+2 ) + 2*2*l0 + 2*2*L[0]*l1]);//(1/2, 0 , 0 ) + (1/2, 0 , 1 ) ; f00(1/2, 0 , 0 )[0][0][l0][l1];
	DX2[3] = D1*D2 + D2*D0 + D0*D1;
	
	//xi - l = ( 0 ,1/2,1/2)
	D0 = fx[1 + 3*1 + 3*2*l0 + 3*2*L[0]*l1];//( 0 ,1/2,1/2) ; f00( 0 ,1/2,1/2)[1][1][l0][l1];
	D1 = 0.5*(fx[1 + 3*0 + 3*2*l0 + 3*2*L[0]*l1] + r_mu[( 1+2 ) + 2*2*l0 + 2*2*L[0]*l1]);//( 0 ,1/2, 0 ) + ( 0 ,1/2, 1 ) ; f11( 0 ,1/2, 0 )[1][0][l0][l1];
	D2 = 0.5*(fx[2 + 3*1 + 3*2*l0 + 3*2*L[0]*l1p] +  fx[2 + 3*1 + 3*2*l0 + 3*2*L[0]*l1]);//( 0 , 1 ,1/2) + ( 0 , 0 ,1/2) ; f22( 0 , 0 ,1/2)[2][1][l0][l1];
	DX2[4] = D1*D2 + D2*D0 + D0*D1;
	
	//xi - l = ( 0 , 0 ,1/2)
	D2 = fx[2 + 3*1 + 3*2*l0 + 3*2*L[0]*l1];//( 0 , 0 ,1/2) ; f22( 0 , 0 ,1/2)[2][1][l0][l1];
	D1 = 0.5*(fx[0 + 3*1 + 3*2*l0 + 3*2*L[0]*l1] + fx[0 + 3*1 + 3*2*l0m + 3*2*L[0]*l1]);//(1/2, 0 ,1/2) + (-1/2, 0 ,1/2) ; f11(1/2, 0 ,1/2)[0][1][l0][l1];
	D0 = 0.5*(fx[1 + 3*1 + 3*2*l0 + 3*2*L[0]*l1] + fx[1 + 3*1 + 3*2*l0 + 3*2*L[0]*l1m]);//( 0 ,1/2,1/2) + ( 0 ,-1/2,1/2) ; f00( 0 ,1/2,1/2)[1][1][l0][l1];
	DX2[5] = D1*D2 + D2*D0 + D0*D1;
}
//======================================================================
//  input_filename ---- rm_space, rm_comma                                 
//                                       Last updated on Aug 27, 2019   
//======================================================================
void rm_space( char *A );
void rm_comma( char *A );
void input_filename(FILE *fp_i1, char Yee[])
{
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = -1;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 0) { 
		rm_space(buf);
		if(strncmp(buf, "Prefix", 6) == 0 && sscanf(buf,"%*[^=] %*[=] %s", Yee) == 1){ //Prefix = Yee
			rm_comma(Yee);
			j_count++;// 0
		}
	}
	if(j_count != 0) {	fprintf(stderr,"1 control commands can not read in input_file1!");	exit(EXIT_FAILURE);}
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
//  input_parameter                  
//                                       Last updated on Sep 11, 2019.  
//======================================================================
void input_parameter(const char Yee[], int L[4], double *omega)
{
	char input_name[BUFSIZE], buf[BUFSIZE], dummy[BUFSIZE];
	snprintf(input_name,sizeof(input_name),"%s_Med0.dat", Yee);
	FILE *fp_i;
	fp_i = fopen(input_name,"r");
	if (fp_i == NULL){	fprintf(stderr,"open error for %s\n", input_name);	exit(EXIT_FAILURE);}
	int j_count = -4;
	double wavelength, k_0;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL && j_count <= 0) { // changed < to <=  on 20190910
		if(strncmp(buf,"# Correction",12) == 0){	fprintf(stderr,"Correction factor has already been added!\n");	exit(EXIT_FAILURE);}// revised on 20190911
		if(j_count < 0){
			if(strncmp(buf,"# info",6) == 0 && sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %s", dummy) == 1){rm_comma(dummy);	j_count++;}
			else if(strncmp(buf,"# wavelength",12) == 0 && sscanf(buf,"%*[^=] %*[=] %lf", &wavelength) == 1){j_count++;}
			else if(strncmp(buf,"# k_0",5) == 0 && sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %lf", &k_0) == 1){j_count++;}
			else if(sscanf(buf,"%*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %*[^=] %*[=] %d", &L[0], &L[1], &L[2], &L[3]) == 4){j_count++;}
		}
	}
	if(fclose(fp_i) != 0){	fprintf(stderr,"fclose error after reading a data file!\n");	exit(EXIT_FAILURE);}
	else if(j_count < 0){	fprintf(stderr,"L,M,N,N_outer can not be read at Scatterer @ input_file0!\n");	exit(EXIT_FAILURE);}
	else{	fprintf(stderr,"wavelength = %.5e, k_0 = %.5e\n", wavelength, k_0); *omega = Pi2/(wavelength*k_0);}
}
