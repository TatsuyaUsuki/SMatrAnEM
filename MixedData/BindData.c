//  Start of Med3D.h
#include "header_macro.h"
#include "constant.h"
void output_file(const int l2, const int L[4], double du0dxi0[2*L[0]+1], double du1dxi1[2*L[1]+1], double du2dxi2[2*L[2]+1], const char *urxF, const char *outerF, const char *Yee, const char *exec);
void get_dudxi(const char *non_uniformF, const int L[4], double du0dxi0[], double du1dxi1[], double du2dxi2[]);
void input_data_file0(const char Yee[], const char non_uniformF[], int L[4]);
void input_file1(FILE *fp_i1, char non_uniformF[], char urxF[], char outerF[], char Yee[]);
//  End of Med3D.h
//======================================================================
//   main ---- input_file1_count, input_file1, output_file              
//======================================================================
int main(int argc, char **argv)
{
	{	time_t timer_ini = time(0);	fprintf(stderr,"# Start time of Yee's lattice creation = %s\n", ctime(&timer_ini));}
	if(argc != 2) {	fprintf(stderr,"error: number of files \n");	exit(EXIT_FAILURE);}
	else if(strncmp(argv[1], "-v", 2) == 0 || strcmp(argv[1], "--version") == 0 ) {
		fprintf(stderr,"The '%s' binds medium-data and scale-factor data.\n", argv[0]);
		fprintf(stderr,"Version 19.08.29 is compiled at %s on %s.\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(stderr," Source code : '%s'\n Author      : Tatsuya Usuki\n URL         : http://www.smatran.org\n", __FILE__);
		fprintf(stderr," References  : 'Wave scattering in frequency domain' as 'Formulation.pdf' on May 20, 2018.\n");
		fprintf(stderr,"There is NO warranty.\n");
		exit(EXIT_SUCCESS);//normal end
	}
//------- begin reading file names -------
	FILE *fp_i1;
	fp_i1 = fopen(argv[1],"r");
	if (fp_i1 == NULL){	fprintf(stderr,"open error!: open input-file1!\n");	exit(EXIT_FAILURE);}
	fprintf(stderr,"The 1st input file: %s\n",argv[1]);
	
	char non_uniformF[BUFSIZE], urxF[BUFSIZE], outerF[BUFSIZE], Yee[BUFSIZE];
	input_file1(fp_i1, non_uniformF, urxF, outerF, Yee);
	fprintf(stderr,"non_uniformF: %s\n", non_uniformF);
	fprintf(stderr,"urxF        : %s\n", urxF);
	fprintf(stderr,"outerF      : %s\n", outerF);
	fprintf(stderr,"Prefix      : %s\n", Yee);
	if(fclose(fp_i1) != 0) {	fprintf(stderr,"fclose error after input_file!\n");	exit(EXIT_FAILURE);}
	int L[4];//L[2] = bottom + scatterer + top, L[3] = outer number
	input_data_file0(Yee, non_uniformF, L);
//-------  end reading file names  -------
	double *du0dxi0 = calloc(2*L[0]+1,sizeof(double));
	double *du1dxi1 = calloc(2*L[1]+1,sizeof(double));
	double *du2dxi2 = calloc(2*L[2]+1,sizeof(double));
	if(du0dxi0 == NULL || du1dxi1 == NULL || du2dxi2 == NULL){fprintf(stderr,"du0dxi0, du1dxi1 and/or du2dxi2 can not be secured!\n");	exit(EXIT_FAILURE);}
	
	get_dudxi(non_uniformF, L, du0dxi0, du1dxi1, du2dxi2);
	for(int l2 = 0 ; l2 < L[2] ; l2++){	output_file(l2, L, du0dxi0, du1dxi1, du2dxi2, urxF, outerF, Yee, argv[0]);}
	SAFEFREE(du0dxi0);	SAFEFREE(du1dxi1);	SAFEFREE(du2dxi2);
}
//======================================================================
//  input_file1 ---- rm_space, rm_comma                                 
//                                       Last updated on Oct 01, 2018   
//======================================================================
void rm_space( char *A );
void rm_comma( char *A );
void input_file1(FILE *fp_i1, char non_uniformF[], char urxF[], char outerF[], char Yee[])
{
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = -4;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 0) { 
		rm_space(buf);
		if(strncmp(buf, "xi2u", 4) == 0 && sscanf(buf,"%*[^=] %*[=] %s", non_uniformF) == 1){
			rm_comma(non_uniformF);
			j_count++;//-3
		}else if(strncmp(buf, "u2r2x", 5) == 0 && sscanf(buf,"%*[^=] %*[=] %s", urxF) == 1) {
			rm_comma(urxF);
			j_count++;//-2
		}else if(strncmp(buf, "outer", 5) == 0 && sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %s", outerF) == 1) {
			rm_comma(outerF);
			j_count++;//-1
		}else if(strncmp(buf, "Prefix", 6) == 0 && sscanf(buf,"%*[^=] %*[=] %s", Yee) == 1){//Prefix = Yee
			rm_comma(Yee);
			j_count++;// 0
		}
	}
	if(j_count != 0) {	fprintf(stderr,"4 control commands can not read in input_file1!");	exit(EXIT_FAILURE);}//revised on 20181001
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
//  input_data_file0                  
//                                       Last updated on Aug 01, 2018.  
//======================================================================
void input_data_file0(const char Yee[], const char non_uniformF[], int L[4])
{
	char input_nameB[BUFSIZE], input_name[BUFSIZE], input_nameT[BUFSIZE], buf[BUFSIZE], dummy[BUFSIZE];
	snprintf(input_nameB,sizeof(input_nameB),"%s_MedB0.dat", Yee);
	snprintf(input_name,sizeof(input_name),"%s_Med0.dat", Yee);
	snprintf(input_nameT,sizeof(input_nameT),"%s_MedT0.dat", Yee);// https://www.ipa.go.jp/security/awareness/vendor/programmingv1/b06_02.html
//------------------------------- Bottom -------------------------------
	FILE *fp_i;
	fp_i = fopen(input_nameB,"r");
	if (fp_i == NULL){	fprintf(stderr,"open error for %s\n", input_nameB);	exit(EXIT_FAILURE);}
	int j_count = -2;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL && j_count < 0) {
		if(strncmp(buf,"# info",6) == 0 && sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %s", dummy) == 1){rm_comma(dummy);	j_count++;}
		else if(sscanf(buf,"%*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %*[^=] %*[=] %d", &L[0], &L[1], &L[2], &L[3]) == 4){j_count++;}
	}
	if(fclose(fp_i) != 0){	fprintf(stderr,"fclose error after reading a data file!\n");	exit(EXIT_FAILURE);}
	else if(j_count < 0){	fprintf(stderr,"L,M,N,N_outer can not be read at Bottom @ input_file0!\n");	exit(EXIT_FAILURE);}
	else if(strcmp(dummy, non_uniformF) != 0){	fprintf(stderr,"Data have no consistency between file names in Bottom!\n"); exit(EXIT_FAILURE);}
//----------------------------- Scatterer  -----------------------------
	int Ld[4];
	fp_i = fopen(input_name,"r");
	if (fp_i == NULL){	fprintf(stderr,"open error for %s\n", input_name);	exit(EXIT_FAILURE);}
	j_count = -2;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL && j_count < 0) { 
		if(strncmp(buf,"# info",6) == 0 && sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %s", dummy) == 1){rm_comma(dummy);	j_count++;}
		else if(sscanf(buf,"%*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %*[^=] %*[=] %d", &Ld[0], &Ld[1], &Ld[2], &Ld[3]) == 4){j_count++;}
	}
	if(fclose(fp_i) != 0){	fprintf(stderr,"fclose error after reading a data file!\n");	exit(EXIT_FAILURE);}
	else if(j_count < 0){	fprintf(stderr,"L,M,N,N_outer can not be read at Scatterer @ input_file0!\n");	exit(EXIT_FAILURE);}
	else if(strcmp(dummy, non_uniformF) != 0){	fprintf(stderr,"Data have no consistency between file names in Scatterer!\n"); exit(EXIT_FAILURE);}
	else if(Ld[0] != L[0] || Ld[1] != L[1] || Ld[2] != L[2] || Ld[3] != L[3]){fprintf(stderr,"Data have no consistency between Bottom and Scatterer!\n"); exit(EXIT_FAILURE);}
//-------------------------------  Top   -------------------------------
	fp_i = fopen(input_nameT,"r");
	if (fp_i == NULL){	fprintf(stderr,"open error for %s\n", input_nameT);	exit(EXIT_FAILURE);}
	j_count = -2;
	while(fgets(buf, sizeof( buf ), fp_i) != NULL && j_count < 0) { 
		if(strncmp(buf,"# info",6) == 0 && sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %s", dummy) == 1){rm_comma(dummy);	j_count++;}
		else if(sscanf(buf,"%*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %*[^=] %*[=] %d", &Ld[0], &Ld[1], &Ld[2], &Ld[3]) == 4){j_count++;}
	}
	if(fclose(fp_i) != 0){	fprintf(stderr,"fclose error after reading a data file!\n");	exit(EXIT_FAILURE);}
	else if(j_count < 0){	fprintf(stderr,"L,M,N,N_outer can not be read at Top @ input_file0!\n");	exit(EXIT_FAILURE);}
	else if(strcmp(dummy, non_uniformF) != 0){	fprintf(stderr,"Data have no consistency between file names in Top!\n"); exit(EXIT_FAILURE);}
	else if(Ld[0] != L[0] || Ld[1] != L[1] || Ld[2] != L[2] || Ld[3] != L[3]){fprintf(stderr,"Data have no consistency between Top and Scatterer!\n"); exit(EXIT_FAILURE);}
}
//======================================================================
//  get_dudxi ---- get_head              Last updated on Aug 02, 2018.  
//======================================================================
void get_head(FILE *fp_i, char *buf);
void get_dudxi(const char *non_uniformF, const int L[4], double du0dxi0[2*L[0]+1], double du1dxi1[2*L[1]+1], double du2dxi2[2*L[2]+1])
{
	FILE *fp_i;
	fp_i = fopen(non_uniformF,"r");
	if (fp_i == NULL){	fprintf(stderr,"open error for %s in get_dudxi\n", non_uniformF);	exit(EXIT_FAILURE);}
	char buf[BUFSIZE];
	double dummy_a, dummy_b;
	get_head(fp_i, buf);
	for(int i1 = 0 ; i1 < 2*L[0] + 1 ; i1++){
		if(sscanf(buf,"%lf, %lf, %lf", &dummy_a, &dummy_b, &du0dxi0[i1]) != 3){	fprintf(stderr,"data_0 can not be read at %d!\n", i1);	exit(EXIT_FAILURE);}
		if(fgets(buf, BUFSIZE, fp_i) == NULL){	fprintf(stderr,"%s can not have data_0!\n", non_uniformF);	exit(EXIT_FAILURE);}
	}
	get_head(fp_i, buf);
	for(int i1 = 0 ; i1 < 2*L[1] + 1 ; i1++){
		if(sscanf(buf,"%lf, %lf, %lf", &dummy_a, &dummy_b, &du1dxi1[i1]) != 3){	fprintf(stderr,"data_1 can not be read at %d!\n", i1);	exit(EXIT_FAILURE);}
		if(fgets(buf, BUFSIZE, fp_i) == NULL){	fprintf(stderr,"%s can not have data_1!\n", non_uniformF);	exit(EXIT_FAILURE);}
	}
	get_head(fp_i, buf);
	for(int i1 = 2*L[3] ; i1 < 2*(L[2]-L[3]) + 1 ; i1++){
		if(sscanf(buf,"%lf, %lf, %lf", &dummy_a, &dummy_b, &du2dxi2[i1]) != 3){	fprintf(stderr,"data_2 can not be read at %d!\n", i1);	exit(EXIT_FAILURE);}
		if(fgets(buf, BUFSIZE, fp_i) == NULL){	fprintf(stderr,"%s cdoes not have any footer!\n", non_uniformF);	exit(EXIT_FAILURE);}
	}
	for(int i1 = 0 ; i1 < 2*L[3] ; i1++){
		du2dxi2[i1] = du2dxi2[2*L[3]];
	}
	for(int i1 = 2*(L[2]-L[3]) + 1 ; i1 < 2*L[2] + 1 ; i1++){
		du2dxi2[i1] = du2dxi2[2*(L[2]-L[3])];
	}
}
void get_head(FILE *fp_i, char *buf)
{
	while(fgets(buf, BUFSIZE, fp_i) != NULL) {	if(strncmp(buf,"#",1) != 0){goto head_of_data;}}
	fprintf(stderr,"head of data can not be read!\n");	exit(EXIT_FAILURE);
	head_of_data:;
}
//======================================================================
// output_file ---- matdata_file, get_h0h2 Last updated on Aug 03, 2018.  
//======================================================================
void get_h0h2(const int l2, const int L[4], const char *u2r2x_name, double h0[], double h2[]);
void matdata_file(const char f_prefix[], const char add_name[], const int Ltot, const int Lout, const int istep, char data_name[]);
void output_file(const int l2, const int L[4], double du0dxi0[2*L[0]+1], double du1dxi1[2*L[1]+1], double du2dxi2[2*L[2]+1], const char *urxF, const char *outerF, const char *Yee, const char *exec)
{
	double *h0 = calloc(2*(2*L[0]+1),sizeof(double));
	double *h2 = calloc(2*(2*L[0]+1),sizeof(double));
	if(h0 == NULL || h2 == NULL){fprintf(stderr,"h0 and/or h2 can not be secured in output_file!\n");	exit(EXIT_FAILURE);}
	if(L[3] <= l2 && l2 < L[2] - L[3]){	get_h0h2(l2, L, urxF, h0, h2);}
	else{	get_h0h2(l2, L, outerF, h0, h2);}
	if(l2 == L[2] - L[3]){	get_h0h2(l2, L, urxF, h0, h2);}//Because of mis-match between outer-region and Top-region
	double *fx = calloc(3*(2*L[0]+1)*(2*L[1]+1)*2,sizeof(double));
	for(int j2 = 0 ; j2 < 2 ; j2++){
		for(int j1 = 0 ; j1 < 2*L[1]+1 ; j1++){
			for(int j0 = 0 ; j0 < 2*L[0]+1 ; j0++){
				double h1 = 1.;
				double H[3];
				if(h0[j0+(2*L[0]+1)*j2] == 0. || h2[j0+(2*L[0]+1)*j2] == 0.){
					fprintf(stderr,"h0=%.5e and/or h2=%.5e can not be read at (%d,2*%d+%d) in output_file!\n",h0[j0+(2*L[0]+1)*j2],h2[j0+(2*L[0]+1)*j2],j0,l2,j2);
					exit(EXIT_FAILURE);
				}
				H[0] = (du0dxi0[j0]*h0[j0+(2*L[0]+1)*j2]);
				H[1] = (du1dxi1[j1]*h1);
				H[2] = (du2dxi2[2*l2+j2]*h2[j0+(2*L[0]+1)*j2]);
				for(int jj = 0 ; jj < 3 ; jj++){
//					fx[jj+3*j0+3*(2*L[0]+1)*j1+3*(2*L[0]+1)*(2*L[1]+1)*j2] = H[jj]; tested on 20180831
					fx[jj+3*j0+3*(2*L[0]+1)*j1+3*(2*L[0]+1)*(2*L[1]+1)*j2] = H[(jj+1)%3]*H[(jj+2)%3]/H[jj];
				}
			}
		}
	}
	SAFEFREE(h0);	SAFEFREE(h2);
	
	char file_name[BUFSIZE], temp_name[BUFSIZE];
	matdata_file(Yee, "_Med", L[2], L[3], l2, file_name);	matdata_file(Yee, "_temp", L[2], L[3], l2, temp_name);//revised on 20190829
	/*if(0 <= l2 && l2 < L[3]){
		snprintf(file_name, sizeof(file_name), "%s_MedB%d.dat",Yee,l2);
		snprintf(temp_name, sizeof(temp_name), "%s_tempB%d.dat",Yee,l2);
	}
	else if(L[3] <= l2 && l2 < L[2] - L[3]){
		snprintf(file_name, sizeof(file_name), "%s_Med%d.dat",Yee,l2 - L[3]);
		snprintf(temp_name, sizeof(temp_name), "%s_temp%d.dat",Yee,l2 - L[3]);
	}
	else if(L[2] - L[3] <= l2 && l2 < L[2]){
		snprintf(file_name, sizeof(file_name), "%s_MedT%d.dat",Yee,l2 - (L[2] - L[3]));
		snprintf(temp_name, sizeof(temp_name), "%s_tempT%d.dat",Yee,l2 - (L[2] - L[3]));
	}*/
	FILE *fp_r, *fp_w;
	fp_r = fopen(file_name,"r");
	fp_w = fopen(temp_name,"w");
	if (fp_r == NULL || fp_w == NULL){	fprintf(stderr,"open error for %s or %s in output_file\n", file_name, temp_name);	exit(EXIT_FAILURE);}
	char buf[BUFSIZE];
	while(fgets(buf, BUFSIZE, fp_r) != NULL){
		char dummy_c[BUFSIZE];
		if(strncmp(buf,"# Created",9) == 0){
			if(sscanf(buf,"%[^\n]", dummy_c) == 1){
				rm_comma(dummy_c);
				time_t timer_f = time(0);
				snprintf(buf, sizeof(buf), "%s, fxjj was added from ' %s ' and ' %s ' by ' %s ' on %s", dummy_c, urxF, outerF, exec, ctime(&timer_f));
				fprintf(fp_w,"%s", buf);
			}
		}else if(strncmp(buf,"# l m 0",7) == 0){
			snprintf(buf, sizeof(buf), "# l m 0 mu_00 mu_11 ep_22 k0*fx_00 k0*fx_11 k0*fx_22 : 'fxjj' is defined as (duj+1dxj+1)*(hj+1)*(duj+2dxj+2)*(hj+2)/(dujdxj*hj)\n");//revised on 20180917//# l m 0 mu_u mu_v ep_w
			fprintf(fp_w,"%s", buf);
		}else if(strncmp(buf,"# l m 1",7) == 0){
			snprintf(buf, sizeof(buf), "# l m 1 ep_11 ep_00 mu_22 k0*fx_11 k0*fx_00 k0*fx_22\n");//revised on 20180917//# l m 1 ep_v ep_u mu_w
			fprintf(fp_w,"%s", buf);
		}else{
			if(sscanf(buf,"%[^\n]", dummy_c) == 1){
				rm_comma(dummy_c);
				int l0, l1, j2;
				double re0, im0, re1, im1, re2, im2;
				if(sscanf(dummy_c,"%d %d %d %lf %lf %lf %lf %lf %lf", &l0, &l1, &j2, &re0, &im0, &re1, &im1, &re2, &im2) == 9){
					if(l0 < 0 || l0 >= L[0]){	fprintf(stderr,"l0 = %d can not be read in output_file\n",l0);	exit(EXIT_FAILURE);}
					if(l1 < 0 || l1 >= L[1]){	fprintf(stderr,"l1 = %d can not be read in output_file\n",l1);	exit(EXIT_FAILURE);}
					snprintf(dummy_c, sizeof(dummy_c), "%d %d %d %.20e %.20e %.20e %.20e %.20e %.20e", l0, l1, j2, re0, im0, re1, im1, re2, im2);
					if(j2 == 0){
						snprintf(buf, sizeof(buf), "%s %.20e %.20e %.20e\n", dummy_c,
						fx[0+3*(2*l0+1)+3*(2*L[0]+1)*(2*l1)+3*(2*L[0]+1)*(2*L[1]+1)*j2],
							fx[1+3*(2*l0)+3*(2*L[0]+1)*(2*l1+1)+3*(2*L[0]+1)*(2*L[1]+1)*j2],
								fx[2+3*(2*l0+1)+3*(2*L[0]+1)*(2*l1+1)+3*(2*L[0]+1)*(2*L[1]+1)*j2]);
					}else if(j2 == 1){
						snprintf(buf, sizeof(buf), "%s %.20e %.20e %.20e\n", dummy_c,
						fx[1+3*(2*l0+1)+3*(2*L[0]+1)*(2*l1)+3*(2*L[0]+1)*(2*L[1]+1)*j2],
							fx[0+3*(2*l0)+3*(2*L[0]+1)*(2*l1+1)+3*(2*L[0]+1)*(2*L[1]+1)*j2],
								fx[2+3*(2*l0)+3*(2*L[0]+1)*(2*l1)+3*(2*L[0]+1)*(2*L[1]+1)*j2]);
					}else{	fprintf(stderr,"j2 can not be read in output_file\n");	exit(EXIT_FAILURE);}
				}
			}
			fprintf(fp_w,"%s", buf);
		}
	}
	SAFEFREE(fx);
	if(fclose(fp_r) == 0 && fclose(fp_w) == 0){	remove(file_name); rename(temp_name, file_name);}
	else{	fprintf(stderr,"fclose error after reading or writing a data file in output_file!\n");	exit(EXIT_FAILURE);}
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
void get_h0h2(const int l2, const int L[4], const char *u2r2x_name, double h0[], double h2[])
{
	FILE *fp_i;
	fp_i = fopen(u2r2x_name,"r");
	if (fp_i == NULL){	fprintf(stderr,"open error for %s in get_h0h2\n", u2r2x_name);	exit(EXIT_FAILURE);}
	char buf[BUFSIZE];
	while(fgets(buf, BUFSIZE, fp_i) != NULL){
		if(strncmp(buf,"# This",6) == 0){	goto EndData;}
		int dummy_i;
		if(sscanf(buf,"%*[^=] %*[=] %d", &dummy_i) == 1){
			if(fgets(buf, BUFSIZE, fp_i) == NULL){	fprintf(stderr,"End of file in get_h0h2\n");	exit(EXIT_FAILURE);}
			if((dummy_i == 2*(l2 - L[3]) || dummy_i == 2*(l2 - L[3]) + 1) && strncmp(buf,"#2l_0 + 1",9) == 0){
				int s_i = dummy_i - 2*(l2 - L[3]);
				for(int i0 = 0 ; i0 < (2*L[0]+1) ; i0++){
					int dl;
					if(fgets(buf, BUFSIZE, fp_i) != NULL){
						double du2, du0, dr2, dr0, dx2, dx0;
						if(sscanf(buf,"%d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf", 
						&dl, &du2, &du0, &dr2, &dr0, &dx2, &dx0, &h2[i0+s_i*(2*L[0]+1)], &h0[i0+s_i*(2*L[0]+1)]) != 9){
							fprintf(stderr,"Data can not be read at 2*l0+1 = %d for %s in get_h0h2\n", dummy_i, u2r2x_name);	exit(EXIT_FAILURE);
						}else{
							if(i0 != dl){fprintf(stderr,"%d != %d in get_h0h2\n", i0, dl);	exit(EXIT_FAILURE);}
							if(h2[i0+s_i*(2*L[0]+1)] == 0. || h0[i0+s_i*(2*L[0]+1)] == 0.){fprintf(stderr,"h2 and/or h0 can not be read in get_h0h2\n");	exit(EXIT_FAILURE);}
							if(dl > 2*L[0] || dl <0){	fprintf(stderr,"2l0+1=%d in get_h0h2",dl);	exit(EXIT_FAILURE);}
						}
					}
				}
			}
		}
	}
	fprintf(stderr,"Footer can not be detected in get_h0h2");	exit(EXIT_FAILURE);
	EndData:;
	if(fclose(fp_i) != 0){	fprintf(stderr,"fclose error after reading a data file in get_h0h2!\n");	exit(EXIT_FAILURE);}
}
