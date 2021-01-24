//  Start of Med3D.h
#include "header_macro.h"
#include "constant.h"
int input_N_table(FILE *fp_i2);
void input_file2(FILE *fp_i2, const int N_table, const double env_para[2], double complex ep_omega[N_table], double complex mu_omega[N_table]);
void output_file(char **argv, const char Yee[], const int N_table, const double complex ep_omega[N_table], const double complex mu_omega[N_table], const char MediumTable[], const double env_para[2], const char factor_file[]);
int input_file1(FILE *fp_i1, char Yee[], char MediumTable[], double env_para[2], char factor_file[]);
double complex epsi(const double Lambda, const double Temperature, char material[], const double charge_dens_cm3[2]);
//  End of Med3D.h
//======================================================================
//   main ---- input_file1, input_file2, output_file       
//======================================================================
int main(int argc, char **argv)
{
	{
		time_t timer_ini = time(0);
		fprintf(stderr,"# Start time of Yee's lattice creation = %s\n", ctime(&timer_ini));
	}
	if(argc != 2) {
		fprintf(stderr,"error: number of files \n");
		exit(1);
	}else if(strncmp(argv[1], "-v", 2) == 0 || strcmp(argv[1], "--version") == 0 ) {
		fprintf(stderr,"The '%s' creates permititivity averaging in quadruple Yee's cell.\n", argv[0]);
		fprintf(stderr,"Version 21.01.24 is compiled at %s on %s.\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(stderr," Source code : '%s'\n Author      : Tatsuya Usuki\n URL         : http://www.smatran.org\n", __FILE__);
		fprintf(stderr," References  : 'Wave scattering in frequency domain' as 'Formulation.pdf' on Nov 01, 2020;\n");
		fprintf(stderr,"There is NO warranty.\n");
		exit(0);//normal end
	}
//------ begin reading sub-mesh information and output file name -------
	FILE *fp_i1;
	fp_i1 = fopen(argv[1],"r");
	if (fp_i1 == NULL){
		fprintf(stderr,"open error!: open input-file1!\n");
		exit(1);
	}
	fprintf(stderr,"The 1st input file: %s\n",argv[1]);
	char Yee[BUFSIZE], MediumTable[BUFSIZE], factor_file[BUFSIZE];
	double env_para[2];// env_para[0]: wavelength, env_para[1]: temperature
//	int factor_no;
	if(input_file1(fp_i1, Yee, MediumTable, env_para, factor_file) != 0) { 
		fprintf(stderr,"input_file1 error!\n");
		exit(1);
	}
	fprintf(stderr,"Yee's lattice output file `%s*.dat'\n", Yee);
	fprintf(stderr,"factor_file = %s\n", factor_file);
	if(fclose(fp_i1) != 0) {
		fprintf(stderr,"fclose error after input_file!\n");
		exit(1);
	}
//------  end reading sub-mesh information and output file name  -------

//---------- begin reading 3D-structure information ------------
	FILE *fp_i2;
	fp_i2 = fopen(MediumTable,"r");
	if (fp_i2 == NULL){
		fprintf(stderr,"open error!: open medium-table (the 1st input file)!\n");
		exit(1);
	}
	fprintf(stderr,"\nMedium-table (the 2nd input file): %s\n",MediumTable);
	
	int N_table = input_N_table(fp_i2); //maximum of structure number
	rewind(fp_i2);
	
	fprintf(stderr,"N_table = %d\n", N_table);
	double complex ep_omega[N_table], mu_omega[N_table];
	
	input_file2(fp_i2, N_table, env_para, ep_omega, mu_omega);
	if(fclose(fp_i2) != 0) {	fprintf(stderr,"fclose error after input_file2\n");	exit(1);}
	
	for(int j0 = 0 ; j0 < N_table ; j0++) { 
		fprintf(stderr,"media_no = %d, ",j0);
		fprintf(stderr,"epsilon = (%.20E, %.20E), ", creal(ep_omega[j0]), cimag(ep_omega[j0]));
		fprintf(stderr,"mu = (%.20E, %.20E)\n", creal(mu_omega[j0]), cimag(mu_omega[j0]));
	}
//----------  end reading 3D-structure information  ------------
	output_file(argv, Yee, N_table, ep_omega, mu_omega, MediumTable, env_para, factor_file);
	{
		time_t timer_f = time(0);
		fprintf(stderr,"\n# Finish time of Yee's lattice creation = %s", ctime(&timer_f));
	}
}
//======================================================================
//  input_file1 ---- rm_space, ScaleUnit, rm_comma                      
//                                       Last updated on Jul 19, 2018.  
//======================================================================
double ScaleUnit(char x[]);// unit variation: km, m, cm, mm, micron, um, nm, deg, degree, rad, radian.
double TempUnit(double x, char s[]);
void rm_space( char *A );
void rm_comma( char *A );
int input_file1(FILE *fp_i1, char Yee[], char MediumTable[], double env_para[2], char factor_file[])
{
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = 0;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 2) { 
		rm_space(buf);
		double R_dummy;
		char A[BUFSIZE], B[BUFSIZE];
		if(strncmp(buf, "Prefix", 6) == 0 && sscanf(buf,"%*[^=] %*[=] %s", Yee) == 1){//Prefix = Yee
			j_count++;// 1
		}else if(strncmp(buf, "MediumPara", 6) == 0
		&& sscanf(buf,"%*[^=] %*[=] %lf %s %*[^=] %*[=] %lf %s %*[^=] %*[=] %s %*[^=] %*[=] %s", &env_para[0], A, &R_dummy, B, MediumTable, factor_file) == 6){
			env_para[0] *= ScaleUnit(A);
			env_para[1] = TempUnit(R_dummy, B);
			j_count++;// 2
		}
	}
	rm_comma(Yee);	rm_comma(MediumTable);
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
//  output_file ---- input_data_file0, gen_f, gen_dSi, input_data_file  
//          matdata_file, get_file_info, muep_cell2cell, media_ave2000  
//                                       Last updated on Dec 06, 2020.  
//======================================================================
void input_data_file0(const char Yee[], int L[5]);
void input_data_file(const char Yee[], const int n_w, const int L[5], 
			const double complex ep_omega[], const double complex mu_omega[], double complex ep_cell[8*L[0]*L[1]*2], 
			double complex mu_cell[8*L[0]*L[1]*2], char input_name[], const int ave_sw);
void get_file_info(const char input_name[], char xi_file[], char str_file[], double head_u2[4]);
void matdata_file(const char f_prefix[], const char add_name[], const int Ltot, const int Lout, const int istep, char data_name[]);
void muep_cell2cell(double complex cell[], const double complex mu_cell[], const int L_u, const int L_v, const int x_shift, const int y_shift, const int z_shift, const int L[2], const int factor_no, const int ave_sw);
double complex media_ave2000(const double complex cell[], const int direction, const int sum_rad2, const double dSi[], const int ave_sw);
void gen_f(const char factor_file[], int *factor_no, int *ave_sw);
void gen_dSi(const char factor_file[], const int sum_rad, double dSi[]);
void output_file(char **argv, const char Yee[], const int N_table, const double complex ep_omega[N_table], const double complex mu_omega[N_table], const char MediumTable[], const double env_para[2], const char factor_file[])
{
// This 'output_file' creates a set of 3D files for optical media.
// Yee[]       : an output file name (char).
	int L[5];//[4] -> [5] revised on 20180712
	input_data_file0(Yee, L);
	if(L[3] > N_table) {
		fprintf(stderr,"N_med > N_table!\n");
		exit(1);
	}
	fprintf(stderr,"From a data-file, L = %d, M = %d, N_bottom + N + N_top = %d, N_med = %d, CellNum = %d\n", L[0], L[1], L[2], L[3], L[4]);
	int factor_no, ave_sw;
	gen_f(factor_file, &factor_no, &ave_sw);
	int sum_rad2 = 2*(2*factor_no + 1);//	int sum_rad = 2*factor_no + 1; 
	{	fprintf(stderr,"sum_rad = 2*factor_no + 1 = %d, ave_sw = %d\n", 2*factor_no + 1, ave_sw);}
	double *dSi = malloc(sum_rad2*sizeof(double));//	double *dSi = calloc(2*sum_rad,sizeof(double));
	gen_dSi(factor_file, 2*factor_no + 1, dSi);

	int shift16L0L1 = 8*L[0]*L[1];//8*L[0]*L[1]*2;
	if(ave_sw != 0){	shift16L0L1 *= 2;}
	double complex *ep_cell = malloc(2*(factor_no+1)*shift16L0L1*sizeof(double complex));
	double complex *mu_cell = malloc(2*(factor_no+1)*shift16L0L1*sizeof(double complex));
	char input_name[BUFSIZE];//revised on 20180718
	for(int jj = 0 ; jj < (factor_no + 1); jj++){
		int iw = jj;	if(iw > L[2] - 1){	iw = L[2] - 1;}
		input_data_file(Yee, iw, L, ep_omega, mu_omega, &ep_cell[(jj + factor_no + 1)*shift16L0L1], &mu_cell[(jj + factor_no + 1)*shift16L0L1], input_name, ave_sw);
	}

	for(int jj = 0 ; jj < (factor_no + 1); jj++){
		for(int ll = 0 ; ll < shift16L0L1 ; ll++){
			ep_cell[ll+jj*shift16L0L1] = ep_cell[ll+(factor_no + 1)*shift16L0L1];
			mu_cell[ll+jj*shift16L0L1] = mu_cell[ll+(factor_no + 1)*shift16L0L1];
		}
	}
	double complex *cell;
	if(ave_sw == 0){
		cell = (double complex *)malloc((sum_rad2)*(sum_rad2)*(sum_rad2)*sizeof(double complex));
	}else{
		cell = (double complex *)malloc(2*(sum_rad2)*(sum_rad2)*(sum_rad2)*sizeof(double complex));
	}
	for(int jw0 = 0 ; jw0 < L[2] ; jw0++) { 
		char output_name[BUFSIZE];
		matdata_file(Yee, "_Med", L[2], L[4], jw0, output_name);
		FILE *fp_o;
		fp_o = fopen(output_name,"w");
		if (fp_o == NULL){	fprintf(stderr,"open error for %s\n", output_name);	exit(1);}
		time_t timer_f = time(0);
		fprintf(fp_o,"# Created on %s", ctime(&timer_f));
		fprintf(fp_o,"# exec = %s : `%s' was compiled at %s on %s by C-version:%ld\n", argv[0], __FILE__, __TIME__, __DATE__, __STDC_VERSION__);
		{
			int read_jw = jw0 + factor_no;
			if(jw0 + factor_no >= L[2]){	read_jw = L[2]-1;}
			input_data_file(Yee, read_jw, L, ep_omega, mu_omega, &ep_cell[(2*factor_no + 1)*shift16L0L1], &mu_cell[(2*factor_no + 1)*shift16L0L1], input_name, ave_sw);
		}
		{
			char xi_file[BUFSIZE], str_file[BUFSIZE]; double head_u2[4];//revised on 20180801
			get_file_info(input_name, xi_file, str_file, head_u2);//revised on 20180718
			fprintf(fp_o,"# info = %s, u-xi = %s, str. = %s, s-map = %s, med. = %s, this = %s\n", 
				argv[1], xi_file, str_file, input_name, MediumTable, output_name);//revised on 20180718
			fprintf(fp_o,"# wavelength [m] = %.20e, temp. [K] = %.20e, sinc_rad = %d\n", env_para[0], env_para[1], factor_no);
			fprintf(fp_o,"# k_0 * u_2  = %.20e ; %.20e ; %.20e, k_0[m^-1] = %.20e\n", head_u2[0], head_u2[1], head_u2[2], head_u2[3]);//revised on 20180801
		}
		fprintf(fp_o,"# l < %d, m < %d, n + 2*outer < %d, media_no < %d, outer = %d,\n", L[0], L[1], L[2], L[3], L[4]);
		fprintf(fp_o,"# l m 0 mu_u mu_v ep_w\n");
		fprintf(fp_o,"# l m 1 ep_v ep_u mu_w\n");
		for(int l_v = 0 ; l_v < L[1] ; l_v++){
			for(int l_u = 0 ; l_u < L[0] ; l_u++){
				{
					fprintf(fp_o,"%d %d 0 ", l_u, l_v);
					muep_cell2cell(cell, mu_cell, l_u, l_v, 2, 1, 1, L, factor_no, ave_sw);//mu_u: x_shift = 2, y_shift = 1, z_shift = 1
					double complex mu_u = media_ave2000(cell,0, sum_rad2, dSi, ave_sw);
//
					muep_cell2cell(cell, mu_cell, l_u, l_v, 1, 2, 1, L, factor_no, ave_sw);//mu_v: x_shift = 1, y_shift = 2, z_shift = 1
					double complex mu_v = media_ave2000(cell,1, sum_rad2, dSi, ave_sw);
//
					muep_cell2cell(cell, ep_cell, l_u, l_v, 2, 2, 1, L, factor_no, ave_sw);//ep_w: x_shift = 2, y_shift = 2, z_shift = 1
					double complex ep_w = media_ave2000(cell,2, sum_rad2, dSi, ave_sw);
//
					fprintf(fp_o,"%.20E %.20E %.20E %.20E %.20E %.20E\n",
						creal(mu_u),cimag(mu_v),creal(mu_v),cimag(mu_v),creal(ep_w),cimag(ep_w));
				}
				
				{
					fprintf(fp_o,"%d %d 1 ", l_u, l_v);
					muep_cell2cell(cell, ep_cell, l_u, l_v, 2, 1, 2, L, factor_no, ave_sw);//ep_v: x_shift = 2, y_shift = 1, z_shift = 2
					double complex ep_v = media_ave2000(cell,1, sum_rad2, dSi, ave_sw);
//
					muep_cell2cell(cell, ep_cell, l_u, l_v, 1, 2, 2, L, factor_no, ave_sw);//ep_u: x_shift = 1, y_shift = 2, z_shift = 2
					double complex ep_u = media_ave2000(cell,0, sum_rad2, dSi, ave_sw);
//
					muep_cell2cell(cell, mu_cell, l_u, l_v, 1, 1, 2, L, factor_no, ave_sw);//mu_w: x_shift = 1, y_shift = 1, z_shift = 2
					double complex mu_w = media_ave2000(cell,2, sum_rad2, dSi, ave_sw);
//
					fprintf(fp_o,"%.20E %.20E %.20E %.20E %.20E %.20E\n",
						creal(ep_v),cimag(ep_v),creal(ep_u),cimag(ep_u),creal(mu_w),cimag(mu_w));
				}
			}
		}
		if(fclose(fp_o) != 0) {	fprintf(stderr,"fclose error after writing a data file!\n");	exit(1);}
//		for(int ll = 0 ; ll < shift16L0L1 ; ll++){
//			for(int ls = 0 ; ls < sum_rad ; ls++){
		for(int ls = 0 ; ls < 2*factor_no + 1 ; ls++){
			for(int ll = 0 ; ll < shift16L0L1 ; ll++){
				ep_cell[ll + shift16L0L1*ls] = ep_cell[ll + shift16L0L1*(ls+1)];
				mu_cell[ll + shift16L0L1*ls] = mu_cell[ll + shift16L0L1*(ls+1)];
			}
		}
	}
	SAFEFREE(ep_cell);	SAFEFREE(mu_cell);	SAFEFREE(cell);	SAFEFREE(dSi);
}
//======================================================================
//  media_ave2000 ---- normal_media, media_ave  
//                                       Last updated on Dec 28, 2020.  
//======================================================================
void normal_media(const int sum_rad2, const double complex cell[], double normal_sqrV[3]);
double complex media_ave(const double complex cell[], const int sum_rad, const double dSi[]);
double complex media_ave2000(const double complex cell[], const int direction, const int sum_rad2, const double dSi[], const int ave_sw)
{
	double complex media_out = 0.;
	if(ave_sw == 4){//Case 4: geometric mean along normal vector only in each sub-cell
			double complex media2;
			{
				double complex media1 = media_ave(cell, sum_rad2, dSi);
				double complex cmedia_inv = media_ave(&cell[sum_rad2*sum_rad2*sum_rad2], sum_rad2, dSi);
				double complex media3 = conj(cmedia_inv)/(creal(cmedia_inv)*creal(cmedia_inv)+cimag(cmedia_inv)*cimag(cmedia_inv));
				media2 = sqrt(media1*media3);
				if(cos(carg(media1+media3)) < 0.) {media2 = -media2;}
			}
			media_out = media2;
	}else if(ave_sw == 3){//Case 3: harmonic mean for global cells
			double complex media2;
			{
				double complex cmedia_inv= media_ave(&cell[sum_rad2*sum_rad2*sum_rad2], sum_rad2, dSi);
				media2 = conj(cmedia_inv)/(creal(cmedia_inv)*creal(cmedia_inv)+cimag(cmedia_inv)*cimag(cmedia_inv));
			}
			media_out = media2;
	}else if(ave_sw == 2){//Case 2: geometric mean along normal vector only in each sub-cell
		double sqrV[3];
		normal_media(sum_rad2, cell, sqrV);
		if(sqrV[direction] > 1. || sqrV[direction] < 0. || direction > 2 || direction <0){
			fprintf(stderr,"Error @ media_ave2000! sqrV[%d] = %.5e\n", direction, sqrV[direction]);	exit(1);
		}else{
			double complex media1, media2;
			{//revised on 20201221
				media1 = media_ave(cell, sum_rad2, dSi);
				double complex cmedia_inv = media_ave(&cell[sum_rad2*sum_rad2*sum_rad2], sum_rad2, dSi);
				double complex media3 = conj(cmedia_inv)/(creal(cmedia_inv)*creal(cmedia_inv)+cimag(cmedia_inv)*cimag(cmedia_inv));
				media2 = sqrt(media1*media3);
				if(cos(carg(media1+media3)) < 0.) {media2 = -media2;}
			}
			media_out = (1.-sqrV[direction])*media1 + sqrV[direction]*media2;
		}
	}else if(ave_sw == 1){//Case 1: harmonic mean along normal vector for global cells
		double sqrV[3];
		normal_media(sum_rad2, cell, sqrV);
		if(sqrV[direction] > 1. || sqrV[direction] < 0. || direction > 2 || direction <0){
			fprintf(stderr,"Error @ media_ave2000! sqrV[%d] = %.5e\n", direction, sqrV[direction]);	exit(1);
		}else{
			double complex media1, media2;
			{
				media1 = media_ave(cell, sum_rad2, dSi);
				double complex cmedia_inv= media_ave(&cell[sum_rad2*sum_rad2*sum_rad2], sum_rad2, dSi);
				media2 = conj(cmedia_inv)/(creal(cmedia_inv)*creal(cmedia_inv)+cimag(cmedia_inv)*cimag(cmedia_inv));
			}
			media_out = (1.-sqrV[direction])*media1 + sqrV[direction]*media2;
		}
	}else if(ave_sw == 0){//Case 0: arithmetic mean
		media_out = media_ave(cell, sum_rad2, dSi);
	}else if(ave_sw == -1){//Case -1: harmonic mean along normal vector only in each sub-cell added on 20201225
		double sqrV[3];
		normal_media(sum_rad2, cell, sqrV);
		if(sqrV[direction] > 1. || sqrV[direction] < 0. || direction > 2 || direction <0){
			fprintf(stderr,"Error @ media_ave2000! sqrV[%d] = %.5e\n", direction, sqrV[direction]);	exit(1);
		}else{
			double complex media1, media2;
			{
				media1 = media_ave(cell, sum_rad2, dSi);
				media2= media_ave(&cell[sum_rad2*sum_rad2*sum_rad2], sum_rad2, dSi);
			}
			media_out = (1.-sqrV[direction])*media1 + sqrV[direction]*media2;
		}
	}else if(ave_sw == -2){//Case -2: geometric mean along normal vector only in each sub-cell added on 20201225
		double sqrV[3];
		normal_media(sum_rad2, cell, sqrV);
		if(sqrV[direction] > 1. || sqrV[direction] < 0. || direction > 2 || direction <0){
			fprintf(stderr,"Error @ media_ave2000! sqrV[%d] = %.5e\n", direction, sqrV[direction]);	exit(1);
		}else{
			double complex media1, media2;
			{
				media1 = media_ave(cell, sum_rad2, dSi);
				double complex media3 = media_ave(&cell[sum_rad2*sum_rad2*sum_rad2], sum_rad2, dSi);
				media2 = sqrt(media1*media3);
				if(cos(carg(media1+media3)) < 0.) {media2 = -media2;}
			}
			media_out = (1.-sqrV[direction])*media1 + sqrV[direction]*media2;
		}
	}else if(ave_sw == -3){//Case -3: harmonic mean only in each sub-cell
		media_out = media_ave(&cell[sum_rad2*sum_rad2*sum_rad2], sum_rad2, dSi);
	}else if(ave_sw == -4){//Case -4: arithmetic mean mixed with harmonic mean only in each sub-cell
			double complex media1, media2;
			{
				media1 = media_ave(cell, sum_rad2, dSi);
				media2 = media_ave(&cell[sum_rad2*sum_rad2*sum_rad2], sum_rad2, dSi);
			}
			media_out = (0.9)*media1 + (0.1)*media2;
	}else if(ave_sw == -5){//Case -5: geometric mean
			double complex media1, media2;
			{
				media1 = media_ave(cell, sum_rad2, dSi);
				media2 = media_ave(&cell[sum_rad2*sum_rad2*sum_rad2], sum_rad2, dSi);
			}
			media_out = sqrt(media1*media2);
			if(cos(carg(media1+media2)) < 0.) {media_out = -media_out;}
	}else if(ave_sw < -5 || ave_sw > 4){
		fprintf(stderr,"Error @ media_ave2000! ave_sw = %d\n", ave_sw);	exit(1);
	}
	return(media_out);
}
double complex media_ave(const double complex cell[], const int sum_rad2, const double dSi[])
{
	double complex mu_u  =  0.;
	for(int ls_w = 0 ; ls_w < sum_rad2 ; ls_w++){
		for(int ls_v = 0 ; ls_v < sum_rad2 ; ls_v++){
			for(int ls_u = 0 ; ls_u < sum_rad2 ; ls_u++){
				mu_u += dSi[ls_u]*dSi[ls_v]*dSi[ls_w]*cell[ls_u+(sum_rad2)*ls_v+(sum_rad2)*(sum_rad2)*ls_w];
			}
		}
	}
	return(mu_u);
}
void normal_media(const int sum_rad2, const double complex cell[], double normal_sqrV[3])
{
	int sum_rad = sum_rad2/2;
	{
		double mu_u  =  0.;
		for(int ls_u = sum_rad - 1 ; ls_u <= sum_rad ; ls_u++){
			double complex sub_sum = 0.;
			for(int ls_w = sum_rad - 1 ; ls_w <= sum_rad ; ls_w++){
				for(int ls_v = sum_rad - 1 ; ls_v <= sum_rad ; ls_v++){
					sub_sum += cell[ls_u+sum_rad2*ls_v+sum_rad2*sum_rad2*ls_w];
				}
			}
			sub_sum *= 0.25;
			mu_u += ((double complex) ((2*ls_u) - sum_rad2 + 1))*sub_sum;
		}
		normal_sqrV[0] = creal(mu_u)*creal(mu_u) + cimag(mu_u)*cimag(mu_u);
	}{
		double mu_v  =  0.;
		for(int ls_v = sum_rad - 1 ; ls_v <= sum_rad ; ls_v++){
			double complex sub_sum = 0.;
			for(int ls_w = sum_rad - 1 ; ls_w <= sum_rad ; ls_w++){
				for(int ls_u = sum_rad - 1 ; ls_u <= sum_rad ; ls_u++){
					sub_sum += cell[ls_u+sum_rad2*ls_v+sum_rad2*sum_rad2*ls_w];
				}
			}
			sub_sum *= 0.25;
			mu_v += ((double complex) ((2*ls_v) - sum_rad2 + 1))*sub_sum;
		}
		normal_sqrV[1] = creal(mu_v)*creal(mu_v) + cimag(mu_v)*cimag(mu_v);
	}{
		double mu_w  =  0.;
		for(int ls_w = sum_rad - 1 ; ls_w <= sum_rad ; ls_w++){
			double complex sub_sum = 0.;
			for(int ls_v = sum_rad - 1 ; ls_v <= sum_rad ; ls_v++){
				for(int ls_u = sum_rad - 1 ; ls_u <= sum_rad ; ls_u++){
					sub_sum += cell[ls_u+sum_rad2*ls_v+sum_rad2*sum_rad2*ls_w];
				}
			}
			sub_sum *= 0.25;
			mu_w += ((double complex) ((2*ls_w) - sum_rad2 + 1))*sub_sum;
		}
		normal_sqrV[2] = creal(mu_w)*creal(mu_w) + cimag(mu_w)*cimag(mu_w);
	}{
		double norm2 = (normal_sqrV[0] + normal_sqrV[1] + normal_sqrV[2]);
		if(norm2 > DBL_EPSILON){
			norm2 = 1./norm2;
			for(int j = 0 ; j < 3 ; j++){normal_sqrV[j] *= norm2;}
		}
	}
}
//======================================================================
//  gen_f                           Last updated on Dec 04, 2020   
//======================================================================
void gen_f(const char factor_file[], int *factor_no, int *ave_sw)
{
	FILE *fp_i;
	fp_i = fopen(factor_file,"r");
	if (fp_i != NULL){
		fprintf(stderr,"Inport from %s\n", factor_file);
		char buf[BUFSIZE];
		while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
			if(sscanf(buf,"%*[^=] %*[=] %d %*[^=] %*[=] %d", factor_no, ave_sw) == 2){	goto EndOfSmoothingData;
			}
		}
		fprintf(stderr,"Cannot get factor_no and/or ave_sw in gen_f\n");	exit(1);
		EndOfSmoothingData:;
		if(fclose(fp_i) != 0) {	fprintf(stderr,"fclose error of %s\n", factor_file);	exit(1);}
	}else{fprintf(stderr,"No file as %s in gen_f()\n", factor_file);}
}
//======================================================================
//  gen_dSi                         Last updated on Jan 05, 2021   
//======================================================================
void gen_dSi(const char factor_file[], const int sum_rad, double dSi[])
{
	{
		double *eta = malloc(sum_rad*sizeof(double));
		eta[0] = 0.5;
		for(int jj = 1; jj < sum_rad; jj++){	eta[jj] = 0.;}//	double eta[13]={0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};//test_sample
		FILE *fp_i;
		fp_i = fopen(factor_file,"r");
		if (fp_i != NULL){
			fprintf(stderr,"Inport from %s\n", factor_file);
			char buf[BUFSIZE];
			int Ndummy; double Rdummy;
			while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
				if(strncmp(buf, "end",3) == 0||strncmp(buf, "End",3) == 0||strncmp(buf, "END",3) == 0){	goto EndOfSmoothingData;
				}else if(sscanf(buf,"%d %*[,] %lf", &Ndummy, &Rdummy) == 2){
//					if(Ndummy >= 0 && Ndummy < 13){eta[Ndummy] = Rdummy;	fprintf(stderr,"eta[%d] = %.5e\n", Ndummy,eta[Ndummy]);}
					if(Ndummy >= 0 && Ndummy < sum_rad){eta[Ndummy] = Rdummy;	fprintf(stderr,"eta[%d] = %.5e\n", Ndummy,eta[Ndummy]);}//revised on 20210113
				}
			}
			fprintf(stderr,"Warning: Cannot get end point in %s\n",factor_file);//	exit(1);
			EndOfSmoothingData:;
			if(fclose(fp_i) != 0) {	fprintf(stderr,"fclose error of %s\n", factor_file);	exit(1);}
		}else{fprintf(stderr,"No file as %s in gen_dSi()\n", factor_file);}
		
		//if(sum_rad > 13) {	fprintf(stderr,"sum_rad > 13, sum_rad = %d\n", sum_rad);	exit(1);}
		for(int jj = 0; jj < sum_rad; jj++){
			dSi[sum_rad -1 -jj] = eta[jj];
		}
		SAFEFREE(eta);
	}
	for(int jj = 0; jj < sum_rad; jj++){
		dSi[2*sum_rad -1 -jj] = dSi[jj];
	}
	{
		double normal_sum = 0.0;
		for(int jj = 0; jj < 2*sum_rad; jj++){	normal_sum += dSi[jj];}
		for(int jj = 0; jj < 2*sum_rad; jj++){	dSi[jj] /= normal_sum;}
	}
	{
		double sum = 0.5;
		for(int jj = 0; jj < sum_rad; jj++){	sum -= dSi[jj];}
		if(fabs(sum) >  2*sum_rad*DBL_EPSILON) {	fprintf(stderr,"sum != 0., sum = %.20E\n", sum);	exit(1);}
	}
}
//======================================================================
//  muep_cell2cell                   Last updated on Dec 06, 2020   
//======================================================================
// factor_no = 2, sum_rad = 2*factor_no + 1 = 5 
// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11 for mu_cell, ep_cell;  < 2*(sum_rad + 1) = 12
// 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 for above/2;  < (sum_rad + 1) = 6
//-3,-3,-2,-2,-1,-1, 0, 0,+1,+1,+2,+2 for shift of L_u, L_v;  -(sum_rad + 1)/2 <=   < (sum_rad + 1)/2 = factor_no + 1 = 3
//    0, 1, 2, 3,{4, 5},6, 7, 8, 9 for ls_u,v,w of cell
//    0, 0, 1, 1,{2, 2},3, 3, 4, 4 for ls_u,v,w/2
//
// factor_no = 1, sum_rad = 2*factor_no + 1 = 3 
// 0, 1, 2, 3, 4, 5, 6, 7 for mu_cell, ep_cell;  < 2*(sum_rad + 1) = 8
// 0, 0, 1, 1, 2, 2, 3, 3 for above/2;  < (sum_rad + 1) = 4
// n, n, 0, 0, 1, 1, 2, 2 for mu_cellJ, ep_cellJ
//-2,-2,-1,-1, 0, 0,+1,+1 for shift of L_u, L_v;  -(sum_rad + 1)/2 <=   < (sum_rad + 1)/2 = factor_no + 1 = 2
//    0, 1,{2, 3},4, 5 for ls_u,v,w of cell
//    0, 0,{1, 1},2, 2 for ls_u,v,w/2
//
//mu_u: x_shift = 2, y_shift = 1, z_shift = 1
//mu_v: x_shift = 1, y_shift = 2, z_shift = 1
//ep_w: x_shift = 2, y_shift = 2, z_shift = 1
//ep_v: x_shift = 2, y_shift = 1, z_shift = 2
//ep_u: x_shift = 1, y_shift = 2, z_shift = 2
//mu_w: x_shift = 1, y_shift = 1, z_shift = 2
void muep_cell2cell(double complex cell[], const double complex mu_cell[], const int L_u, const int L_v, const int x_shift, const int y_shift, const int z_shift, const int L[2], const int factor_no, const int ave_sw)
{
	int sum_rad2 = (2*factor_no + 1)*2; 
	for(int ls_w = 0 ; ls_w < sum_rad2 ; ls_w++){//(int ls_w = 0 ; ls_w < 6 ; ls_w++){
		int l_w, s_w;
		{int lshift_w = ls_w + z_shift;	l_w = lshift_w/2;	s_w = lshift_w%2;}//1 <= lshift_w <= 5+2; 0 <= l_w < 4
		for(int ls_v = 0 ; ls_v < sum_rad2 ; ls_v++){//for(int ls_v = 0 ; ls_v < 6 ; ls_v++){
			int l_v, s_v;
			{int lshift_v = ls_v + y_shift;	l_v = lshift_v/2;	s_v = lshift_v%2;}
			int Ll_v = (L_v + l_v - (factor_no + 1));//(L_v + l_v - 2 + 2*L[1])%L[1];
			while(Ll_v < 0){	Ll_v += L[1];}
			Ll_v %= L[1];//Ll_v = Ll_v%L[1];
			for(int ls_u = 0 ; ls_u < sum_rad2 ; ls_u++){//for(int ls_u = 0 ; ls_u < 6 ; ls_u++){
				int l_u, s_u;
				{int lshift_u = ls_u + x_shift;	l_u = lshift_u/2;	s_u = lshift_u%2;}
				int Ll_u = (L_u + l_u - (factor_no + 1));//(L_u + l_u - 2 + 2*L[0])%L[0];
				while(Ll_u < 0){	Ll_u += L[0];}
				Ll_u %= L[0];//Ll_u = Ll_u%L[0];
				if(ave_sw == 0){
					cell[ls_u+(sum_rad2)*ls_v+(sum_rad2)*(sum_rad2)*ls_w] = mu_cell[s_u+2*s_v+4*s_w + 8*Ll_u + 8*L[0]*Ll_v + 8*L[0]*L[1]*l_w ];
				}else{
					cell[ls_u+(sum_rad2)*ls_v+(sum_rad2)*(sum_rad2)*ls_w                                   ] = mu_cell[s_u+2*s_v+4*s_w + 8*Ll_u + 8*L[0]*Ll_v               + 16*L[0]*L[1]*l_w ];
					cell[ls_u+(sum_rad2)*ls_v+(sum_rad2)*(sum_rad2)*ls_w + (sum_rad2)*(sum_rad2)*(sum_rad2)] = mu_cell[s_u+2*s_v+4*s_w + 8*Ll_u + 8*L[0]*Ll_v + 8*L[0]*L[1] + 16*L[0]*L[1]*l_w ];
				}
			}
		}
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
//  get_file_info                    Last updated on Aug 01, 2018   
//======================================================================
void get_file_info(const char input_name[], char xi_file[], char str_file[], double head_u2[4])//revised on 20180801
{
	FILE *fp_i;
	fp_i = fopen(input_name,"r");
	if (fp_i == NULL){
		fprintf(stderr,"open error for %s in get_file_info\n", input_name);
		exit(1);
	}
	char buf[BUFSIZE], dummy[BUFSIZE];//
	int i_count = 0;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL && i_count < 2) { 
		if(strncmp(buf,"# info",6) == 0) {//revised on 20180801
			if(sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %s %*[^=] %*[=] %s %*[^=] %*[=] %s", xi_file, str_file, dummy) != 3){fprintf(stderr,"File names can not read in get_file_info!\n");	exit(1);}
			rm_comma(xi_file); rm_comma(str_file); rm_comma(dummy);
			if(strcmp(input_name, dummy) != 0){fprintf(stderr,"Error! %s != %s in get_file_info!\n", input_name, dummy);	exit(1);}
			i_count++;
		}else if(strncmp(buf,"# k_0 * u_2",11) == 0){//revised on 20180801
			if(sscanf(buf,"%*[^=] %*[=] %lf %*[;] %lf %*[;] %lf %*[^=] %*[=] %lf", &head_u2[0], &head_u2[1], &head_u2[2], &head_u2[3]) != 4){fprintf(stderr,"Parameters can not read in get_file_info!\n");	exit(1);}
			i_count++;
		}
	}
	if(fclose(fp_i) != 0){	fprintf(stderr,"fclose error in get_file_info\n");	exit(1);}
	if(i_count != 2){	fprintf(stderr,"Error! i_count = %d != 2 in get_file_info!\n", i_count);	exit(1);}
}
//======================================================================
//  input_data_file                  Last updated on Dec 04, 2020   
//======================================================================
void input_data_file(const char Yee[], const int n_w, const int L[5], 
			const double complex ep_omega[], const double complex mu_omega[], double complex ep_cell[8*L[0]*L[1]], 
			double complex mu_cell[8*L[0]*L[1]], char input_name[], const int ave_sw)
{
	FILE *fp_i;
	matdata_file(Yee, "", L[2], L[4], n_w, input_name);

	fp_i = fopen(input_name,"r");
	if (fp_i == NULL){
		fprintf(stderr,"open error for %s in input_data_file\n", input_name);
		exit(1);
	}
	char buf[BUFSIZE];
	int dummy[4];
	while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
		if(sscanf(buf,"%*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d", &dummy[0], &dummy[1], &dummy[2], &dummy[3]) == 4) goto End_of_LMN;
	}
	fprintf(stderr,"No L,M,N,N_med are read @ input_file1!\n");
	exit(1);
	End_of_LMN:;
	if(L[0] != dummy[0] || L[1] != dummy[1] || L[3] != dummy[3]){//revised on 20180712
		fprintf(stderr,"L,M,N,N_med are incorrect in %s!\n", input_name);
		exit(1);
	}
	
	for(int l_v = 0 ; l_v < L[1] ; l_v++){
		for(int l_u = 0 ; l_u < L[0] ; l_u++){
			while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
				int m_u, m_v;
				if(sscanf(buf,"%d %d", &m_u, &m_v) == 2) {
					if(m_u == l_u && m_v == l_v) goto Next_Step_0;
				}
			}
			fprintf(stderr,"l_u and l_v were not read @ input_data_file!\n");
			exit(1);
			Next_Step_0:;
			for(int k_w = 0 ; k_w < 2 ; k_w++){
				for(int k_v = 0 ; k_v < 2 ; k_v++){
					for(int k_u = 0 ; k_u < 2 ; k_u++){
						int N_med;
						while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
							int m_u, m_v, m_w;
							if(sscanf(buf,"%d %d %d %d", &m_u, &m_v, &m_w, &N_med) == 4) {
								if(m_u == k_u && m_v == k_v && m_w == k_w) goto Next_Step_1;
							}
						}
						fprintf(stderr,"k_u, k_v, k_w and N_med were not read @ input_data_file!\n");
						exit(1);
						Next_Step_1:;
						int j_med;
						double weight_med;
						double sum_weight_med = 0.;
						if(k_u+2*k_v + 4*k_w + 8*l_u + 8*L[0]*l_v > 8*L[0]*L[1]-1){fprintf(stderr,"Loop error in input_data_file\n"); exit(1);}
						ep_cell[k_u+2*k_v + 4*k_w + 8*l_u + 8*L[0]*l_v] = 0.;
						mu_cell[k_u+2*k_v + 4*k_w + 8*l_u + 8*L[0]*l_v] = 0.;
						if(ave_sw != 0){
							ep_cell[k_u+2*k_v + 4*k_w + 8*l_u + 8*L[0]*l_v + 8*L[0]*L[1]] = 0.;
							mu_cell[k_u+2*k_v + 4*k_w + 8*l_u + 8*L[0]*l_v + 8*L[0]*L[1]] = 0.;
						}
						for(int med_no = 0 ; med_no < N_med ; med_no++){
							while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
								if(sscanf(buf,"%d %lf", &j_med, &weight_med) == 2) goto Next_Step_2;
							}
							fprintf(stderr,"j_med and weight_med were not read @ input_data_file!\n");
							exit(1);
							Next_Step_2:;
							sum_weight_med += weight_med;
							if(j_med >= L[3]) {
								fprintf(stderr,"j_med >= N_med! j_med = %d, N_med = %d\n", j_med, N_med);
								exit(1);
							}
							{
								double complex ep = ep_omega[j_med];
								ep_cell[k_u+2*k_v + 4*k_w + 8*l_u + 8*L[0]*l_v] += weight_med*ep;
								double complex mu = mu_omega[j_med];
								mu_cell[k_u+2*k_v + 4*k_w + 8*l_u + 8*L[0]*l_v] += weight_med*mu;
								if(ave_sw != 0){//revised on 20201204
									ep_cell[k_u+2*k_v + 4*k_w + 8*l_u + 8*L[0]*l_v + 8*L[0]*L[1]] += weight_med*conj(ep)/(creal(ep)*creal(ep)+cimag(ep)*cimag(ep));
									mu_cell[k_u+2*k_v + 4*k_w + 8*l_u + 8*L[0]*l_v + 8*L[0]*L[1]] += weight_med*conj(mu)/(creal(mu)*creal(mu)+cimag(mu)*cimag(mu));
								}
							}
						}
//						if(ave_sw != 0 && ave_sw != 1){//added on 20201218. Note that a case of ave_sw = 1 sets harmonic mean for global cells
						if(ave_sw < 0){//added on 20201218. Note that a case of ave_sw = 1 sets harmonic mean for global cells
							double complex ep = ep_cell[k_u+2*k_v + 4*k_w + 8*l_u + 8*L[0]*l_v + 8*L[0]*L[1]];
							ep_cell[k_u+2*k_v + 4*k_w + 8*l_u + 8*L[0]*l_v + 8*L[0]*L[1]] = conj(ep)/(creal(ep)*creal(ep)+cimag(ep)*cimag(ep));
							double complex mu = mu_cell[k_u+2*k_v + 4*k_w + 8*l_u + 8*L[0]*l_v + 8*L[0]*L[1]];
							mu_cell[k_u+2*k_v + 4*k_w + 8*l_u + 8*L[0]*l_v + 8*L[0]*L[1]] = conj(mu)/(creal(mu)*creal(mu)+cimag(mu)*cimag(mu));
						}
						if(sum_weight_med != 1.) {
							fprintf(stderr,"sum_weight_med = %g\n", sum_weight_med);
							exit(1);
						}
					}
				}
			}
		}
	}

	if(fclose(fp_i) != 0) {
		fprintf(stderr,"fclose error after reading a data file!\n");
		exit(1);
	}
}
//======================================================================
//  input_data_file0                 Last updated on Jul 23, 2018   
//======================================================================
void input_data_file0(const char Yee[], int L[5])//revised on 20180712
{
	FILE *fp_i;
	char input_nameB[BUFSIZE], input_name[BUFSIZE], input_nameT[BUFSIZE], buf[BUFSIZE];
	snprintf(input_nameB,sizeof(input_nameB),"%sB0.dat", Yee);//revised on 20180723 https://www.ipa.go.jp/security/awareness/vendor/programmingv1/b06_02.html
	snprintf(input_name,sizeof(input_name),"%s0.dat", Yee);//revised on 20180723 https://www.ipa.go.jp/security/awareness/vendor/programmingv1/b06_02.html
	snprintf(input_nameT,sizeof(input_nameT),"%sT0.dat", Yee);//revised on 20180723 https://www.ipa.go.jp/security/awareness/vendor/programmingv1/b06_02.html
	int Ld[5];
//------------------------------- Bottom -------------------------------
	fp_i = fopen(input_nameB,"r");
	if (fp_i == NULL){
		fprintf(stderr,"open error for %s\n", input_nameB);
		exit(1);
	}
	while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
		if(sscanf(buf,"%*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d", &L[0], &L[1], &L[4], &L[3]) == 4) goto End_of_LMNB;
	}
	fprintf(stderr,"No L,M,N,N_med are read at Bottom @ input_file1!\n");
	exit(1);
	End_of_LMNB:;
	if(fclose(fp_i) != 0) {
		fprintf(stderr,"fclose error after reading a data file!\n");
		exit(1);
	}
//----------------------------- Scatterer  -----------------------------
	fp_i = fopen(input_name,"r");
	if (fp_i == NULL){
		fprintf(stderr,"open error for %s\n", input_name);
		exit(1);
	}
	while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
		if(sscanf(buf,"%*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d", &Ld[0], &Ld[1], &L[2], &Ld[3]) == 4) goto End_of_LMN;
	}
	fprintf(stderr,"No L,M,N,N_med are read at Scatterer @ input_file1!\n");
	exit(1);
	End_of_LMN:;
	if(Ld[0] != L[0] || Ld[1] != L[1] || Ld[3] != L[3]){fprintf(stderr,"Data have no consistency between Bottom and Scatterer!\n"); exit(1);}
	if(fclose(fp_i) != 0) {
		fprintf(stderr,"fclose error after reading a data file!\n");
		exit(1);
	}
//-------------------------------  Top   -------------------------------
	fp_i = fopen(input_nameT,"r");
	if (fp_i == NULL){
		fprintf(stderr,"open error for %s\n", input_nameT);
		exit(1);
	}
	while(fgets(buf, sizeof( buf ), fp_i) != NULL) { 
		if(sscanf(buf,"%*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d", &Ld[0], &Ld[1], &Ld[4], &Ld[3]) == 4) goto End_of_LMNT;
	}
	fprintf(stderr,"No L,M,N,N_med are read at Top @ input_file1!\n");
	exit(1);
	End_of_LMNT:;
	if(Ld[0] != L[0] || Ld[1] != L[1] || Ld[3] != L[3] || Ld[4] != L[4]){fprintf(stderr,"Data have no consistency between Bottom and Scatterer!\n"); exit(1);}
	if(fclose(fp_i) != 0) {
		fprintf(stderr,"fclose error after reading a data file!\n");
		exit(1);
	}
//----------------------------- Scatterer  -----------------------------
	L[2] += 2*L[4];
}
//======================================================================
//  double ScaleUnit, TempUnit                                                    
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
double TempUnit(double x, char s_char[])
{
	double offset_temp = -1.;
	double scale = -1.;
	if (strncmp(s_char, "K", 1) == 0) {
		scale = 1e-0;
		offset_temp = 0.;
	}// Kelvin
	if (strncmp(s_char, "degC", 4) == 0 || strncmp(s_char, "deg.C", 5) == 0 || strncmp(s_char, "deg_C", 5) == 0) {
		scale = 1.;
		offset_temp = 273.15;
	}// Celsius
	if (strncmp(s_char, "degF", 4) == 0 || strncmp(s_char, "deg.F", 5) == 0 || strncmp(s_char, "deg_F", 5) == 0) {
		scale = 5./9.;
		offset_temp = 459.67;
	}// Fahrenheit
	if (scale < 0. || offset_temp < -0.1) {
		fprintf(stderr,"TempUnit error!: %s\n", s_char);
		exit(1);
	}
	return(scale*(x + offset_temp));
}
//======================================================================
//  int input_N_table                                                   
//                                       Last updated on May 23, 2020.  
//======================================================================
int input_N_table(FILE *fp_i2)
{
	char buf[BUFSIZE], command[16]; // buffer for fgets.
	int n_table0 = 0;
	int n_table1 = -10;
	while(fgets(buf, sizeof( buf ), fp_i2) != NULL) {
		char number_int[16];
		sscanf(buf,"%s", number_int);
		if (strcmp(number_int, "end") == 0||strcmp(number_int, "End") == 0||strcmp(number_int, "END") == 0) break;
		if(sscanf(buf,"%d %s", &n_table1, command) == 2){
			if(n_table1 == n_table0) {n_table0 += 1;}
		}
	}
	return(n_table0);
}
//======================================================================
//  void input_file2                                                    
//                                       Last updated on Oct 22, 2014.  
//======================================================================
void input_file2(FILE *fp_i2, const int N_table, const double env_para[2], double complex ep_omega[N_table], double complex mu_omega[N_table])
{
	char buf[BUFSIZE], command[16];	// buffer for fgets.
	int n_table0 = 0;
	int n_table1 = -10;
	while(fgets(buf, sizeof( buf ), fp_i2) != NULL) { 
		char number_int[16];
		sscanf(buf,"%s", number_int);
		if (strcmp(number_int, "end") == 0||strcmp(number_int, "End") == 0||strcmp(number_int, "END") == 0) break;
		if(sscanf(buf,"%d %s", &n_table1, command) == 2){
			if(n_table1 == n_table0){
				int dummy;
				mu_omega[n_table1] = 1.;
				if (strcmp(command, "name") == 0) {
					char material[BUFSIZE];
					double charge_dens_cm3[2];
					if(sscanf(buf,"%d %s %s %lf %lf", &dummy, command, material, &charge_dens_cm3[0], &charge_dens_cm3[1]) < 3){
						fprintf(stderr,"command, material error!: buf = %s", buf);
						exit(1);
					}
					ep_omega[n_table1] = epsi(env_para[0], env_para[1], material, charge_dens_cm3);
					goto EPSI;
				}
				if (strcmp(command, "index") == 0) {
					double Re_n, Im_n;
					if(sscanf(buf,"%d %s %*[^=] %*[=] %lf %*[^=] %*[=] %lf", &dummy, command, &Re_n, &Im_n) != 4){
						fprintf(stderr,"command, Re_n, Im_n error!: buf = %s", buf);
						exit(1);
					}
					ep_omega[n_table1] = (Re_n + I*Im_n)*(Re_n + I*Im_n);
					goto EPSI;
				}
				if (strcmp(command, "epsilon")== 0) {
					double Re_epsi, Im_epsi;
					if(sscanf(buf,"%d %s %*[^=] %*[=] %lf %*[^=] %*[=] %lf", &dummy, command, &Re_epsi, &Im_epsi) != 4){
						fprintf(stderr,"command, Re_epsi, Im_epsi error! `%s'\n", buf);
						exit(1);
					}
					ep_omega[n_table1] = Re_epsi + I*Im_epsi;
					goto EPSI;
				}
				if (strcmp(command, "epsilon+mu")== 0) {
					double Re_epsi, Im_epsi, Re_mu, Im_mu;
					if(sscanf(buf,"%d %s %*[^=] %*[=] %lf %*[^=] %*[=] %lf %*[^=] %*[=] %lf %*[^=] %*[=] %lf", &dummy, command, &Re_epsi, &Im_epsi, &Re_mu, &Im_mu) != 6){
						fprintf(stderr,"command, epsilon+mu error! `%s'\n", buf);
						exit(1);
					}
					ep_omega[n_table1] = Re_epsi + I*Im_epsi;
					mu_omega[n_table1] = Re_mu + I*Im_mu;
					goto EPSI;
				}
				fprintf(stderr,"cannot get epsilon! `%s'\n", buf);
				exit(1);
				EPSI:;
				n_table0 += 1;
			}
		}
		if(n_table0 == N_table) goto Step_end;
	}
	fprintf(stderr,"Never find final step!  n_table0 = %d, N_table = %d\n", n_table0, N_table);
	Step_end:;
}

