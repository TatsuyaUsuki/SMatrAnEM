//  Start of PMLcalc.h
#include "header_macro.h"
#include "constant.h"
void output_PML(const int l2, const char Yee[], const char exec[], const char PMLfile[], const int L[5], 
				double u0[2*L[0]+1], double u1[2*L[1]+1],
				const int oriPML[L[4]], 
				const double startPML[L[4]], const double thickPML[L[4]], 
				const double sqRef[L[4]], const int powerPML[L[4]], double omega);
void trans_para(const int L[5], const int istaPML[L[4]], const int iendPML[L[4]], 
				const double u0[2*L[0]+1], const double u1[2*L[1]+1], const double u2[2*L[2]+1], 
				int oriPML[L[4]], double startPML[L[4]], double thickPML[L[4]]);
void input_u2(char Yee[], const int L[5], double u2[2*L[2]+1]);
void input_u(const char u0u1file[], const int L[5], double u0[2*L[0]+1], double u1[2*L[1]+1]);
void input_PMLn(const char Yee[], const char PMLfile[], 
				const int L[5], int oriPML[L[4]], 
				int istaPML[L[4]], int iendPML[L[4]], 
				double sqRef[L[4]], int powerPML[L[4]]);
int input_N_str(const char PMLfile[]);
void input_L(const char Yee[], const char PMLfile[], int L[5], double *omega);
void input_filename(FILE *fp_i1, char Yee[], char PMLfile[], char u0u1file[]);
//  End of PMLcalc.h
//======================================================================
//   main ---- input_filename, input_L, input_PML, input_u, output_PML  
//======================================================================
int main(int argc, char **argv)
{
	time_t timer_ini = time(0);	fprintf(stderr,"# Start time of PML setting = %s\n", ctime(&timer_ini));
	if(argc != 2) {fprintf(stderr,"error: number of files \n");	exit(EXIT_FAILURE);}
	else if(strncmp(argv[1], "-v", 2) == 0 || strcmp(argv[1], "--version") == 0 ) {
		fprintf(stderr,"The '%s' creates perfectly matched layer (PML).\n", argv[0]);
		fprintf(stderr,"Version 19.09.18 is compiled at %s on %s.\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(stderr," Source code : '%s'\n Author      : Tatsuya Usuki\n URL         : http://www.smatran.org\n", __FILE__);
		fprintf(stderr," Reference   : Equation (6.6) in Discretization by Yee's lattice as 'Formulation.pdf' on Aug 25, 2019;\n");
		fprintf(stderr,"There is NO warranty.\n");
		exit(EXIT_SUCCESS);//normal end
	}
//------- begin reading file names and parameters -------
	FILE *fp_i1;
	fp_i1 = fopen(argv[1],"r");
	if (fp_i1 == NULL){	fprintf(stderr,"open error!: open input-file1!\n");	exit(EXIT_FAILURE);}
	fprintf(stderr,"The 1st input file: %s\n",argv[1]);
	char Yee[BUFSIZE], PMLfile[BUFSIZE], u0u1file[BUFSIZE];
	input_filename(fp_i1, Yee, PMLfile, u0u1file);
	fprintf(stderr,"Prefix: %s, PMLfile: %s\n", Yee, PMLfile);
	if(fclose(fp_i1) != 0) {	fprintf(stderr,"fclose error after input_file!\n");	exit(EXIT_FAILURE);}
	int L[5];//L[2] = bottom + scatterer + top, L[3] = outer number, L[4] = PML number
	double omega = 0.;
	input_L(Yee, PMLfile, L, &omega);
	fprintf(stderr,"L0 = %d, L1 = %d, Lbst = %d, Louter = %d, PML number = %d \n", L[0], L[1], L[2], L[3], L[4]);
//-------  end reading file names and parameters  -------
	{
		int *oriPML;	oriPML = malloc(sizeof(int)*L[4]);
		int *istaPML;	istaPML = malloc(sizeof(int)*L[4]);
		int *iendPML;	iendPML = malloc(sizeof(int)*L[4]);
		double *sqRef;	sqRef = malloc(sizeof(double)*L[4]);
		int *powerPML;	powerPML = malloc(sizeof(int)*L[4]);
		input_PMLn(Yee, PMLfile, L, oriPML, istaPML, iendPML, sqRef, powerPML);//revised on 201917
		{
			double *u0;	u0 = malloc(sizeof(double)*(2*L[0]+1));
			double *u1;	u1 = malloc(sizeof(double)*(2*L[1]+1));
			input_u(u0u1file, L, u0, u1);
			double *startPML;	startPML = malloc(sizeof(double)*L[4]);
			double *thickPML;	thickPML = malloc(sizeof(double)*L[4]);
			{
				double *u2;	u2 = malloc(sizeof(double)*(2*L[2]+1));
				input_u2(Yee, L, u2);//revised on 201917
//				for(int l2 = 0 ; l2 < 2*L[2]+1 ; l2++){fprintf(stderr,"u2[2*%d+%d]=%.5e[k0^-1]\n", l2/2, l2%2, u2[l2]);}
				trans_para(L, istaPML, iendPML, u0, u1, u2, oriPML, startPML, thickPML);//revised on 201917
				SAFEFREE(u2);
			}
			for(int l2 = 0 ; l2 < L[2] ; l2++){
				output_PML(l2, Yee, argv[0], PMLfile, L, u0, u1, oriPML, startPML, thickPML, sqRef, powerPML, omega);
			}
			SAFEFREE(u0);	SAFEFREE(u1);	SAFEFREE(startPML);	SAFEFREE(thickPML);
		}
		SAFEFREE(oriPML);	SAFEFREE(istaPML);	SAFEFREE(iendPML);	SAFEFREE(sqRef);	SAFEFREE(powerPML);
	}
}
//======================================================================
//  trans_para                           Last updated on Sep 18, 2019.  
//======================================================================
void trans_para(const int L[5], const int istaPML[L[4]], const int iendPML[L[4]], 
				const double u0[2*L[0]+1], const double u1[2*L[1]+1], const double u2[2*L[2]+1], 
				int oriPML[L[4]], double startPML[L[4]], double thickPML[L[4]])
{
	double endPML = 0.;
	fprintf(stderr,"Start trans_para\n");
	for(int j_pml = 0 ; j_pml < L[4] ; j_pml++){
		if(istaPML[j_pml] > iendPML[j_pml]){
			if(oriPML[j_pml] == 1){
				startPML[j_pml] = u0[2*istaPML[j_pml] +2];
				endPML = u0[2*iendPML[j_pml]];
			}else if(oriPML[j_pml] == 2){
				startPML[j_pml] = u1[2*istaPML[j_pml] +2];
				endPML = u1[2*iendPML[j_pml]];
			}else if(oriPML[j_pml] == 3){
				startPML[j_pml] = u2[2*istaPML[j_pml] +2];
				endPML = u2[2*iendPML[j_pml]];
			}
			oriPML[j_pml] *= -1;
		}else{
			if(oriPML[j_pml] == 1){
				startPML[j_pml] = u0[2*istaPML[j_pml]];
				endPML = u0[2*iendPML[j_pml] +2];
			}else if(oriPML[j_pml] == 2){
				startPML[j_pml] = u1[2*istaPML[j_pml]];
				endPML = u1[2*iendPML[j_pml] +2];
			}else if(oriPML[j_pml] == 3){
				startPML[j_pml] = u2[2*istaPML[j_pml]];
				endPML = u2[2*iendPML[j_pml] +2];
			}
		}
		thickPML[j_pml] = fabs(endPML - startPML[j_pml]);
		fprintf(stderr,"oriPML[%d]=%d, startPML[%d]=%.5e[k0^-1], thickPML[%d]=%.5e[k0^-1]\n", j_pml, oriPML[j_pml], j_pml, startPML[j_pml], j_pml, thickPML[j_pml]);
	}
}
//======================================================================
//  input_u2  ----  matdata_file         Last updated on Sep 17, 2019.  
//======================================================================
void matdata_file(const char f_prefix[], const char add_name[], const int Ltot, const int Lout, const int istep, char data_name[]);
void input_u2(char Yee[], const int L[5], double u2[2*L[2]+1])
{
	char data_name[BUFSIZE];
	for(int j2 = 0 ; j2 < L[2] ; j2++){
		matdata_file(Yee, "_Med", L[2], L[3], j2, data_name);
		FILE *fp_i;
		fp_i = fopen(data_name,"r");
		if (fp_i == NULL){	fprintf(stderr,"open error!: open input_u2!\n");	exit(EXIT_FAILURE);}

		char buf[BUFSIZE];
		while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
			if(strncmp(buf,"# k_0 * u_2",11) == 0){
				if(sscanf(buf,"%*[^=] %*[=] %lf %*[;] %lf %*[;] %lf", &u2[2*j2], &u2[2*j2+1], &u2[2*j2+2]) == 3){goto Nextstep;}
				else{fprintf(stderr,"read error1 in input_u2\n");	exit(EXIT_FAILURE);}
			}
		}
		Nextstep:;

		if(fclose(fp_i) != 0) {	fprintf(stderr,"fclose error after input_u2!\n");	exit(EXIT_FAILURE);}
	}
}
//======================================================================
//  output_PML  ----  matdata_file, set_header, out_flag, calc_sigma,
//                    sqDt_CFL           Last updated on Sep 06, 2019.  
//======================================================================
void calc_sigma(FILE *fp_i, FILE *fp_o, const int L[5], 
				const double u0[2*L[0]+1], const double u1[2*L[1]+1], const double u2[2], 
				const int oriPML[L[4]], 
				const double startPML[L[4]], const double thickPML[L[4]], 
				const double sqRef[L[4]], const int powerPML[L[4]], const int out_region, const double omega);
int out_flag(const int Ltot, const int Lout, const int istep);
void set_header(FILE *fp_i, double u2[2], FILE *fp_o, const int N_pml, const char *exec, const char PMLfile[], const double sqDt);
double sqDt_CFL(FILE *fp_i, const int L[5]);
void matdata_file(const char f_prefix[], const char add_name[], const int Ltot, const int Lout, const int istep, char data_name[]);
void output_PML(const int l2, const char Yee[], const char exec[], const char PMLfile[], const int L[5], 
				double u0[2*L[0]+1], double u1[2*L[1]+1],
				const int oriPML[L[4]], 
				const double startPML[L[4]], const double thickPML[L[4]], 
				const double sqRef[L[4]], const int powerPML[L[4]], double omega)
{//	fprintf(stderr,"k/k_0 = oemga = %.5e\n", omega);
	char data_name[BUFSIZE], temp_name[BUFSIZE];
	matdata_file(Yee, "_Med", L[2], L[3], l2, data_name);
	if(l2 == 0){fprintf(stderr,"From %s ... ", data_name);}
	if(l2 == L[2]-1){fprintf(stderr,"To %s\n", data_name);}
	matdata_file(Yee, "_PML", L[2], L[3], l2, temp_name);
	
	FILE *fp_r, *fp_w;
	fp_r = fopen(data_name,"r");
	fp_w = fopen(temp_name,"w");
	if (fp_r == NULL || fp_w == NULL){	fprintf(stderr,"open error for %s or %s in output_PML\n", data_name, temp_name);	exit(EXIT_FAILURE);}
	double sqDt = sqDt_CFL(fp_r, L);
	double u2[2];
	set_header(fp_r, u2, fp_w, L[4], exec, PMLfile, sqDt);//	char buf[BUFSIZE];	if(fgets(buf, sizeof( buf ), fp_r) != NULL) {	fprintf(fp_w,"temporary %s", buf);}
	{
		int out_region = out_flag(L[2], L[3], l2);
		calc_sigma(fp_r, fp_w, L, u0, u1, u2, oriPML, startPML, thickPML, sqRef, powerPML, out_region, omega);
	}
	if(fclose(fp_r) != 0 ||fclose(fp_w) != 0 ){	fprintf(stderr,"fclose error after reading or writing data files!\n");	exit(EXIT_FAILURE);}
	else{if(remove(data_name)==0){rename(temp_name, data_name);}}
}
//======================================================================
//  calc_sigma  ----  set_sigma           Last updated on Sep 09, 2019  
//======================================================================
double set_sigma(const double x, 
				const int oriPML, 
				const double startPML, const double thickPML, 
				const double sqRef, const int powerPML, const double omega);
void calc_sigma(FILE *fp_i, FILE *fp_o, const int L[5], 
				const double u0[2*L[0]+1], const double u1[2*L[1]+1], const double u2[2], 
				const int oriPML[L[4]], 
				const double startPML[L[4]], const double thickPML[L[4]], 
				const double sqRef[L[4]], const int powerPML[L[4]], const int out_region, const double omega)
{
	char buf[BUFSIZE];
	int j_count = 0;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
		if(strncmp(buf,"# ",2) != 0){
			int l0, l1, j2;
			double re[3], im[3];
			double fx00, fx11, fx22;
			if(sscanf(buf,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", &l0, &l1, &j2, &re[0], &im[0], &re[1], &im[1], &re[2], &im[2], &fx00, &fx11, &fx22) == 12){
				if(l0 < 0 || l0 >= L[0] || l1 < 0 || l1 >= L[1] || j2 < 0 || j2 >= 2){fprintf(stderr,"%d, %d, %d error in calc_sigma!\n", l0, l1, j2);	exit(EXIT_FAILURE);}
				double complex ce[3];	for(int jdata = 0 ; jdata < 3 ; jdata++){ ce[jdata] = re[jdata] + I*im[jdata];}
				double x[9];
				for(int jdata = 0 ; jdata < 3 ; jdata++){
					int j0,j1;
					if(jdata == 2){ j0 = 1-j2;	j1 = 1-j2;}
					else{	j0 = 1 - jdata;	j1 = jdata;}// Note that eq.(6.6) and check header comment of data file!
					x[0 + 3*jdata] = u0[2*l0+j0];//x[0 + 3*0] = u0[2*l0+1]; x[0 + 3*1] = u0[2*l0+0]; x[0 + 3*2] = u0[2*l0+1-j2];
					x[1 + 3*jdata] = u1[2*l1+j1];//x[1 + 3*0] = u1[2*l1+0]; x[1 + 3*1] = u1[2*l1+1]; x[1 + 3*2] = u1[2*l1+1-j2];
					x[2 + 3*jdata] = u2[j2];//
				}
				/*fprintf(fp_o,"%d %d %d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %d %.5e\n", l0, l1, j2, 
					x[0 + 3*0], x[1 + 3*0], x[2 + 3*0], 
					x[0 + 3*1], x[1 + 3*1], x[2 + 3*1], 
					x[0 + 3*2], x[1 + 3*2], x[2 + 3*2], out_region, omega); */
				{	double sigma[9]; for(int j = 0 ; j < 9 ; j++){sigma[j] = 0.;}
					for(int jdata = 0 ; jdata < 3 ; jdata++){
						for(int j_pml = 0 ; j_pml < L[4] ; j_pml++){
							int jxyz = abs(oriPML[j_pml]) - 1;
							if(jxyz < 0 || jxyz >= 3){fprintf(stderr,"oriPML[%d] = %d in calc_sigma!\n", j_pml, oriPML[j_pml]);	exit(EXIT_FAILURE);}
							else if(jxyz == 2 || out_region == 0){
								sigma[jxyz+3*jdata] += set_sigma(x[jxyz + 3*jdata], oriPML[j_pml], startPML[j_pml], thickPML[j_pml], sqRef[j_pml], powerPML[j_pml], omega);//20190908 edited
							}
						}
					}
					for(int jdata = 0 ; jdata < 3 ; jdata++){
						int jxyz;
						if(j2 == 1 && jdata < 2){	jxyz = 1-jdata;}else{	jxyz = jdata;}// Note that eq.(6.6) and check header comment of data file!
						ce[jdata] *= (1. + I*sigma[(jxyz+1)%3 +3*jdata])*(1. + I*sigma[(jxyz+2)%3 +3*jdata])/(1. + I*sigma[(jxyz)%3 +3*jdata]);
					}
				}
				fprintf(fp_o,"%d %d %d %.20e %.20e %.20e %.20e %.20e %.20e %.20e %.20e %.20e\n", l0, l1, j2, 
					creal(ce[0]), cimag(ce[0]), 
					creal(ce[1]), cimag(ce[1]), 
					creal(ce[2]), cimag(ce[2]), fx00, fx11, fx22);
				j_count++;
			}
		}
	}
	if(j_count != L[0]*L[1]*2){fprintf(stderr,"%d != L[0]*L[1]*2 in calc_sigma\n", j_count);	exit(EXIT_FAILURE);}
}
//======================================================================
//  set_sigma                             Last updated on Sep 09, 2019  
//======================================================================
double set_sigma(const double x, 
				const int oriPML, 
				const double startPML, const double thickPML, 
				const double sqRef, const int powerPML, const double omega)
{
	double endPML= startPML; if(oriPML > 0){	endPML += thickPML;}else if(oriPML < 0){	endPML -= thickPML;}
	if((x - startPML)*(x - endPML) <= 0.){// sqRef = exp[-4*omega*thickPML*max_sigma/(powerPML+1)], see section 3.3
		double max_sigma = (-0.25*(powerPML+1)/(omega*thickPML))*log(sqRef);
//		double sigma_x = max_sigma * pow(fabs((x - startPML)/(endPML - startPML)), (double) powerPML);
		double sigma_x = max_sigma * pow(fabs((x - startPML)/(endPML - startPML)), powerPML);
		return(sigma_x);
	}else{
		return(0.);
	}
}
//======================================================================
//  out_flag                              Last updated on Sep 06, 2019  
//======================================================================
int out_flag(const int Ltot, const int Lout, const int istep)
{
	int flag;
	if(0 <= istep && istep < Lout){
		flag = 1;//	snprintf(data_name,BUFSIZE*sizeof(char),"%s%sB%d.dat", f_prefix, add_name, istep);
	}else if(Lout <= istep && istep < Ltot - Lout){
		flag = 0;//	snprintf(data_name,BUFSIZE*sizeof(char),"%s%s%d.dat", f_prefix, add_name, istep - Lout);
	}else if(Ltot - Lout <= istep && istep < Ltot){
		flag = 1;//	snprintf(data_name,BUFSIZE*sizeof(char),"%s%sT%d.dat", f_prefix, add_name, istep - (Ltot - Lout));
	}else{fprintf(stderr,"istep = %d error in out_flag!\n", istep);	exit(EXIT_FAILURE);}
	return(flag);
}
//======================================================================
//  sqDt_CFL                              Last updated on Sep 10, 2019  
//======================================================================
double sqDt_CFL(FILE *fp_i, const int L[5])
{
	if (fp_i == NULL){	fprintf(stderr,"open error in set_u2\n");	exit(EXIT_FAILURE);}
	char buf[BUFSIZE];
	int j_count = 0;
	double sqDt = -1.;
	int l0p = -1;	int l1p = -1;
	double mul00 = -1.;	double mul11 = -1.;	double epl22 = 0.;
	double epl11, epl00, mul22;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
		if(strncmp(buf,"# ",2) != 0){
			int l0, l1, j2;
			double re[3], im[3];
			double fx00, fx11, fx22;
			if(sscanf(buf,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", &l0, &l1, &j2, &re[0], &im[0], &re[1], &im[1], &re[2], &im[2], &fx00, &fx11, &fx22) == 12){
				if(l0 < 0 || l0 >= L[0] || l1 < 0 || l1 >= L[1] || j2 < 0 || j2 >= 2){fprintf(stderr,"%d, %d, %d error in sqDt_CFL!\n", l0, l1, j2);	exit(EXIT_FAILURE);}
				j_count++;
			}
			if(j2 == 0){
				mul00 = re[0]*fx00; mul11 = re[1]*fx11; epl22 = re[2]*fx22;
				l0p = l0; l1p = l1;
			}else{
				epl11 = re[0]*fx00; epl00 = re[1]*fx11; mul22 = re[2]*fx22;
				double rdummy = epl00*epl11*epl22*mul00*mul11*mul22/((epl00+epl11+epl22)*(mul00+mul11+mul22));
				if(rdummy > 0.){
					if(sqDt < 0. || sqDt > rdummy){sqDt = rdummy;}
				}else{fprintf(stderr,"rdummy <= 0 at  (%d, %d) error in sqDt_CFL!\n", l0, l1);	exit(EXIT_FAILURE);}
				if(l0 != l0p || l1 != l1p || j2 != 1){fprintf(stderr,"(%d, %d) = (%d, %d) or j2 = %d error in sqDt_CFL!\n", l0, l1, l0p, l1p, j2);	exit(EXIT_FAILURE);}
			}
		}
	}
	if(j_count != L[0]*L[1]*2){fprintf(stderr,"%d != L[0]*L[1]*2 in sqDt_CFL\n", j_count);	exit(EXIT_FAILURE);}
	rewind(fp_i);
	return(sqrt(sqDt));
}
//======================================================================
//  set_header                            Last updated on Sep 10, 2019  
//======================================================================
void set_header(FILE *fp_i, double u2[2], FILE *fp_o, const int N_pml, const char *exec, const char PMLfile[], const double sqDt)
{
	if (fp_i == NULL){	fprintf(stderr,"open error in set_u2\n");	exit(EXIT_FAILURE);}
	char buf[BUFSIZE];
	while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
		if(strncmp(buf,"# k_0 * u_2",11) == 0){
			if(sscanf(buf,"%*[^=] %*[=] %lf %*[;] %lf", &u2[0], &u2[1]) == 2){
				fprintf(fp_o,"%s", buf);
			}
			else{fprintf(stderr,"read error1 in set_u2\n");	exit(EXIT_FAILURE);}
		}else if(strncmp(buf,"# l <",5) == 0){
			char dummy_c[BUFSIZE];
			if(sscanf(buf,"%[^\n]", dummy_c) == 1){
				time_t timer_f = time(0);
				snprintf(buf, sizeof(buf), "%s PML num = %d, min Dt^2 = %.20e\n# PML and CFL data were added from ' %s ' by ' %s ' on %s", dummy_c, N_pml, sqDt, PMLfile, exec, ctime(&timer_f));
				fprintf(fp_o,"%s", buf);
			}else{fprintf(stderr,"read error2 in set_u2\n");	exit(EXIT_FAILURE);}
		}else{	
			if(strncmp(buf,"# ",2) == 0){fprintf(fp_o,"%s", buf);}
			else{goto Nextstep;}
		}
	}
	Nextstep:;
	rewind(fp_i);//	if(fgets(buf, sizeof( buf ), fp_i) != NULL) {	fprintf(fp_o,"%s", buf);}
}
//======================================================================
//  matdata_file                          Last updated on Aug 29, 2019  
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
//  input_u                              Last updated on Sep 05, 2019.  
//======================================================================
void input_u(const char u0u1file[], const int L[5], double u0[2*L[0]+1], double u1[2*L[1]+1])
{
	char buf[BUFSIZE];
	FILE *fp_i;
	fp_i = fopen(u0u1file,"r");
	if (fp_i == NULL){	fprintf(stderr,"open error for %s\n", u0u1file);	exit(EXIT_FAILURE);}
	int j_count = -(2*L[0]+1 + 2*L[1]+1);
	while(fgets(buf, sizeof( buf ), fp_i) != NULL && j_count < 0) { 
		if(strncmp(buf,"# xi0",5) == 0){
			for(int j0 = 0 ; j0 < 2*L[0]+1 ; j0++){
				if(fgets(buf, sizeof( buf ), fp_i) && 
				sscanf(buf,"%*[^,] %*[,] %lf", &u0[j0]) == 1){j_count++;}
			}
		}else if(strncmp(buf,"# xi1",5) == 0){
			for(int j1 = 0 ; j1 < 2*L[1]+1 ; j1++){
				if(fgets(buf, sizeof( buf ), fp_i) && 
				sscanf(buf,"%*[^,] %*[,] %lf", &u1[j1]) == 1){j_count++;}
			}
		}
	}
	if(j_count != 0){	fprintf(stderr,"u0 and u1 can not be read @ input_u!\n");	exit(EXIT_FAILURE);}
}
//======================================================================
//  input_PMLn  ---- input_orientation, input_sqRef                     
//                                       Last updated on Sep 17, 2019.  
//======================================================================
int input_orientation(const int L[5], const char buf[], const int n_pml0, int oriPML[], int istaPML[], int iendPML[]);
int input_sqRef(const char buf[], const int n_pml0, double sqRef[], int powerPML[]);
void rm_space( char *A );
void rm_comma( char *A );
void input_PMLn(const char Yee[], const char PMLfile[], 
				const int L[5], int oriPML[L[4]], 
				int istaPML[L[4]], int iendPML[L[4]], 
				double sqRef[L[4]], int powerPML[L[4]])
{
	int N_pml = L[4];
	char input_name[BUFSIZE], buf[BUFSIZE];
	snprintf(input_name,sizeof(input_name),"%s_Med0.dat", Yee);
	FILE *fp_i;
	fp_i = fopen(input_name,"r");
	if (fp_i == NULL){	fprintf(stderr,"open error for %s\n", input_name);	exit(EXIT_FAILURE);}
	int j_count = -1;
	double k_0;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL && j_count < 0) { 
		if(strncmp(buf,"# k_0",5) == 0 && sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %lf", &k_0) == 1){fprintf(stderr,"k_0 = %.5e\n", k_0);	j_count++;}
	}
	if(fclose(fp_i) != 0){	fprintf(stderr,"fclose error1 in input_PML!\n");	exit(EXIT_FAILURE);}
	else if(j_count < 0){	fprintf(stderr,"k_0 can not be read in input_PML!\n");	exit(EXIT_FAILURE);}
// (int) 0 <= N_pml <= 6    : total number of structures.
// (int) oriPML[N_pml]      : 'x' == +1, 'y' == +2, 'z' == +3.
// (int) istaPML[N_pml]     : start point of the PML,
// (int) iendPML[N_pml]     :  end  point of the PML.
// (double) sqRef[N_pml]    : reflectivity of the PML.
// (int) powerPML[N_pml]    : index of depth dependence for the PML
	fp_i = fopen(PMLfile,"r");
	if (fp_i == NULL){	fprintf(stderr,"open error for %s\n", input_name);	exit(EXIT_FAILURE);}
	char command[16];	// buffer for fgets
	int n_pml0 = 0;
	for(int J0 = 0 ; J0 < N_pml ; J0++) { 
		while(fgets(buf, sizeof( buf ), fp_i) != NULL) {
			rm_space(buf);
			if(strncmp(buf, "#", 1) != 0 && sscanf(buf,"%s", command) != EOF){
				rm_comma(command);
				if (strcmp(command, "begin") == 0) {
					int ori_count = 0;// Initialization
					int sqRef_count = 0;// End of initialization
					while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
						rm_space(buf);
						if(strncmp(buf, "#", 1) != 0 && sscanf(buf,"%s", command) != EOF){
							rm_comma(command);
							if (strncmp(command, "orientation", 11) == 0) {
								ori_count += input_orientation(L, buf, n_pml0, oriPML, istaPML, iendPML);
							}else if (strncmp(command, "squareRef", 9) == 0) {
								sqRef_count += input_sqRef(buf, n_pml0, sqRef, powerPML);
							}else if (strcmp(command, "end") == 0) {
								if(ori_count != 1 || sqRef_count != 1){
									fprintf(stderr,"end error @ input_PML!: ori_count = %d, sqRef_count = %d\n", 
										ori_count, sqRef_count);
									exit(EXIT_FAILURE);
								}
								else{
									if(n_pml0 != J0){
										fprintf(stderr,"Error @ input_PML!: n_pml0 > N_pml, n_pml0 = %d,", n_pml0);
										exit(EXIT_FAILURE);
									}else{
										fprintf(stderr,"oriPML[%d]=%d, istaPML[%d]=%d, iendPML[%d]=%d, ", J0, oriPML[J0], J0, istaPML[J0], J0, iendPML[J0]);
										fprintf(stderr,"sqRef[%d]=%.5e, powerPML[%d]=%d\n", J0, sqRef[J0], J0, powerPML[J0]);
										n_pml0 += 1;
										goto NEXT_step;
									}
								}
							}
						}
					}
				}
			}
		}
		NEXT_step:;
	}
	if(n_pml0 != N_pml){
		fprintf(stderr,"Error @ input_file3!: n_pml0 != N_pml, n_pml0 = %d,", n_pml0);
		exit(1);
	}else if(fclose(fp_i) != 0){	fprintf(stderr,"fclose error2 in input_PML!\n");	exit(EXIT_FAILURE);}
	for(int J0 = 0 ; J0 < N_pml ; J0++) { 
		if(oriPML[J0] == 0 || abs(oriPML[J0]) > 3){	fprintf(stderr,"oriPML[%d] == oriPML[%d]\n", J0, oriPML[J0]);	exit(EXIT_FAILURE);}
	}
}
//======================================================================
//  input_orientation  ---- num_u2       Last updated on Sep 17, 2019.  
//======================================================================
int num_u2(const int L[5], const char buf[]);
int input_orientation(const int L[5], const char buf[], const int n_pml0, int oriPML[], int istaPML[], int iendPML[])
{
	char para[BUFSIZE], A[BUFSIZE], B[BUFSIZE];
	if(sscanf(buf,"%*[^=] %*[=] %s %*[^=] %*[=] %s %*[^=] %*[=] %s", para, A, B) == 3){
		rm_comma(para);	rm_comma(A);	rm_comma(B);
		if(strncmp(para,"z",1) == 0 ){//revised on 20190917
			oriPML[n_pml0] = 3;
			istaPML[n_pml0] = num_u2(L, A);
			iendPML[n_pml0] = num_u2(L, B);
		}else if(strncmp(para,"y",1) == 0 ){//revised on 20190917
			oriPML[n_pml0] = 2;
			sscanf(A,"%d", &istaPML[n_pml0]);
			sscanf(B,"%d", &iendPML[n_pml0]);
		}else if(strncmp(para,"x",1) == 0 ){//revised on 20190917
			oriPML[n_pml0] = 1;
			sscanf(A,"%d", &istaPML[n_pml0]);
			sscanf(B,"%d", &iendPML[n_pml0]);
		}else{
			oriPML[n_pml0] = 0;
		}
		return(1);
	}else{
		return(0);
	}
}
//======================================================================
//  num_u2                              Last updated on Sep 17, 2019.  
//======================================================================
int num_u2(const int L[5], const char buf[])
{
	int iI = -1;
	if(strncmp(buf, "B", 1) == 0){sscanf(&buf[1],"%d", &iI);}
	else if(strncmp(buf, "T", 1) == 0){sscanf(&buf[1],"%d", &iI);iI += L[2] - L[3];}
	else{sscanf(buf,"%d", &iI);iI += L[3];}
	return(iI);
}
//======================================================================
//  input_sqRef                          Last updated on Sep 05, 2019.  
//======================================================================
int input_sqRef(const char buf[], const int n_pml0, double sqRef[], int powerPML[])
{
	if(sscanf(buf,"%*[^=] %*[=] %lf %*[^=] %*[=] %d", 
		&sqRef[n_pml0], &powerPML[n_pml0]) == 2){
//		fprintf(stderr,"sqRef[%d]=%.5e, powerPML[%d]=%d\n", n_pml0, sqRef[n_pml0], n_pml0, powerPML[n_pml0]);
		return(1);
	}else{
//		fprintf(stderr,"input_start error! `%s'\n", buf);	exit(EXIT_FAILURE);
		return(0);
	}
}
//======================================================================
//  input_L ---- input_N_str                                    
//                                       Last updated on Sep 18, 2019.  
//======================================================================
int input_N_str(const char PMLfile[]);
void input_L(const char Yee[], const char PMLfile[], int L[5], double *omega)
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
	if(fclose(fp_i) != 0){	fprintf(stderr,"fclose error after reading a data file!\n");	exit(EXIT_FAILURE);}
	else if(j_count < 0){	fprintf(stderr,"L,M,N,N_outer can not be read at Scatterer @ input_file0!\n");	exit(EXIT_FAILURE);}
	else{ *omega = Pi2/(wavelength*k_0);	fprintf(stderr,"wavelength = %.5e, k_0 = %.5e, omega = %.5e\n", wavelength, k_0, *omega);}
	L[4] = input_N_str(PMLfile);
}
//======================================================================
//  input_filename ---- rm_space, rm_comma                                 
//                                       Last updated on Sep 05, 2019   
//======================================================================
//void rm_space( char *A );
//void rm_comma( char *A );
void input_filename(FILE *fp_i1, char Yee[], char PMLfile[], char u0u1file[])
{
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = -3;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 0) { 
		rm_space(buf);
		if(strncmp(buf, "xi2u", 4) == 0 && sscanf(buf,"%*[^=] %*[=] %s", u0u1file) == 1){ 
			rm_comma(u0u1file);
			j_count++;// -2
		}else if(strncmp(buf, "Prefix", 6) == 0 && sscanf(buf,"%*[^=] %*[=] %s", Yee) == 1){ //Prefix = Yee
			rm_comma(Yee);
			j_count++;// -1
		}else if(strncmp(buf, "PML", 3) == 0 && sscanf(buf,"%*[^=] %*[=] %s", PMLfile) == 1){ 
			rm_comma(PMLfile);
			j_count++;// 0
		}
	}
	if(j_count != 0) {	fprintf(stderr,"Control commands can not read in input_filename!");	exit(EXIT_FAILURE);}
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
//  input_N_str                                                         
//                                        Last updated on Sep 18, 2014  
//======================================================================
int input_N_str(const char PMLfile[])
{
	FILE *fp_i3;
	fp_i3 = fopen(PMLfile,"r");
	if (fp_i3 == NULL){	fprintf(stderr,"open error in input_N_str!\n");	exit(EXIT_FAILURE);}
	char buf[16], command[8];	// buffer for fgets, command for seeking "begin" and "end".
	int n_str0 = 0;
	while(fgets(buf, sizeof( buf ), fp_i3) != NULL) { 
		if(sscanf(buf,"%s", command) != EOF){
			if (strcmp(command, "begin") == 0) {
				while(fgets(buf, sizeof( buf ), fp_i3) != NULL) { 
					if(sscanf(buf,"%s", command) != EOF){
						if ( strcmp(command, "end") == 0) {
							n_str0 += 1;
							goto NEXT_step;
						}
					}
				}
			}
		}
		NEXT_step:;
	}
	if(fclose(fp_i3) != 0){	fprintf(stderr,"fclose error in input_N_str!\n");	exit(EXIT_FAILURE);}
	return(n_str0);
}
