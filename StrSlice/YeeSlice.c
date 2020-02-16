//  Start of Med3D.h
#include "header_macro.h"
#include "constant.h"
void output_file(const char *Yee, const int s_axis, const int s_2lp1, const char *s_name);
int input_file1_count(FILE *fp_i1);
void input_file1(FILE *fp_i1, const int s_num, char Yee[], int s_axis[], int s_2lp1[], char s_name[]);
//  End of Med3D.h
//======================================================================
//   main ---- input_file1_count, input_file1, output_file
//======================================================================
int main(int argc, char **argv)
{
	{	time_t timer_ini = time(0);	fprintf(stderr,"# Start time of Yee's lattice creation = %s\n", ctime(&timer_ini));}
	if(argc != 2) {	fprintf(stderr,"error: number of files \n");	exit(EXIT_FAILURE);}
	else if(strncmp(argv[1], "-v", 2) == 0 || strcmp(argv[1], "--version") == 0 ) {
		fprintf(stderr,"The '%s' slices Med3D data to cross-sections.\n", argv[0]);
		fprintf(stderr,"Version 18.09.17 is compiled at %s on %s.\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(stderr," Source code : '%s'\n Author      : Tatsuya Usuki\n URL         : http://www.smatran.org\n", __FILE__);
		fprintf(stderr," References  : 'Wave scattering in frequency domain' as 'Formulation.pdf' on May 20, 2018.\n");
		fprintf(stderr,"There is NO warranty.\n");
		exit(EXIT_SUCCESS);//normal end
	}
//------ begin reading sub-mesh information and output file name -------
	FILE *fp_i1;
	fp_i1 = fopen(argv[1],"r");
	if (fp_i1 == NULL){	fprintf(stderr,"open error!: open input-file1!\n");	exit(EXIT_FAILURE);}
	fprintf(stderr,"The 1st input file: %s\n",argv[1]);
	int s_num = input_file1_count(fp_i1);
	
	rewind(fp_i1);
	
	char *s_name = calloc(BUFSIZE*s_num,sizeof(char));
	int *s_axis = calloc(s_num,sizeof(int));
	int *s_2lp1 = calloc(s_num,sizeof(int));
	if(s_name == NULL || s_axis == NULL || s_2lp1 == NULL){fprintf(stderr,"s_name, s_axis and/or s_2lp1 can not be secured!\n");	exit(EXIT_FAILURE);}
	
	char Yee[BUFSIZE];
	input_file1(fp_i1, s_num, Yee, s_axis, s_2lp1, s_name);
	fprintf(stderr,"Prefix: %s\n", Yee);
	if(fclose(fp_i1) != 0) {	fprintf(stderr,"fclose error after input_file!\n");	exit(EXIT_FAILURE);}
//--------------------- begin slicing Yee lattice  ---------------------
	for(int i1 = 0 ; i1 < s_num ; i1++){
		fprintf(stderr,"\ncount=%d,	", i1);
		output_file(Yee, s_axis[i1], s_2lp1[i1], &s_name[BUFSIZE*i1]);//revised on 20180801
	}
	{	time_t timer_f = time(0);	fprintf(stderr,"\n# Finish time of Yee-slice = %s", ctime(&timer_f));}
	SAFEFREE(s_name);	SAFEFREE(s_axis);	SAFEFREE(s_2lp1);
}
//======================================================================
//  input_file1 ---- rm_space, rm_comma                      
//                                       Last updated on Jul 23, 2018.  
//======================================================================
void rm_space( char *A );
void rm_comma( char *A );
void input_file1(FILE *fp_i1, const int s_num, char Yee[], int s_axis[], int s_2lp1[], char s_name[])
{
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = -1;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < s_num) { 
		rm_space(buf);
		
		if(strncmp(buf, "Prefix", 6) == 0 && sscanf(buf,"%*[^=] %*[=] %s", Yee) == 1){	rm_comma(Yee);	j_count++;}// 0
		else if((strncmp(buf, "Slice", 5) == 0 || strncmp(buf, "slice", 5) == 0 || strncmp(buf, "SLICE", 5) == 0) && j_count > -1 && j_count < s_num){
			if(buf[5] == '0'){s_axis[j_count] = 0;}
			else if(buf[5] == '1'){s_axis[j_count] = 1;}
			else if(buf[5] == '2'){s_axis[j_count] = 2;}
			else{fprintf(stderr,"s_axis can not be read in input_file1\n");	exit(EXIT_FAILURE);}
			int l0, s0; 
			char A[BUFSIZE];
			if(sscanf(buf,"%*[^=] %*[=] %d %*[^=] %*[=] %d %*[^=] %*[=] %s", &l0, &s0, A) == 3){
				s_2lp1[j_count] = 2*l0 + s0;
				rm_comma(A);
				strncpy(&s_name[BUFSIZE*j_count], A, BUFSIZE-1);//revised on 20180801
			}else{fprintf(stderr,"s_2lp1 and s_name can not be read\n");	exit(EXIT_FAILURE);}
			j_count++;// 0+s_num
		}
	}
	if(j_count != s_num) {	fprintf(stderr,"j_count = %d != s_num = %d", j_count, s_num);	exit(EXIT_FAILURE);}
}
//======================================================================
//  input_file1_count ---- rm_space                                     
//                                       Last updated on Jul 23, 2018.  
//======================================================================
int input_file1_count(FILE *fp_i1)
{
	char buf[BUFSIZE];	// buffer for fgets
	int s_num = 0;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL) { 
		rm_space(buf);
		if(strncmp(buf, "Slice", 5) == 0 || strncmp(buf, "slice", 5) == 0 || strncmp(buf, "SLICE", 5) == 0){
			s_num++;
		}
	}	
	return(s_num);
}
//======================================================================
//  This program removes spaces of head in characters,
//       it needs #include <ctype.h>      Last updated on Jun 05, 2018. 
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
//  This program outputs slice-data from Yee lattice.
//  output_file ---- input_file2, Case_u0, Case_u1, Case_u2             
//                                        Last updated on Jul 31, 2018  
//======================================================================
void Case_u1(const char *Yee, const char *non_uniform, const int L[], const int s_2lp1, const char *s_name, const char *u2r2x_name, const char *outer_name);
void Case_u0(const char *Yee, const char *non_uniform, const int L[], const int s_2lp1, const char *s_name, const char *u2r2x_name, const char *outer_name);
void Case_u2(const char *Yee, const char *non_uniform, const int L[], const int s_2lp1, const char *s_name, const char *u2r2x_name, const char *outer_name);
void input_file2(const char Yee[], int L[4], char non_uniform[BUFSIZE], char u2r2x_name[BUFSIZE], char outer_name[BUFSIZE]);
void output_file(const char *Yee, const int s_axis, const int s_2lp1, const char *s_name)
{//https://www.ipa.go.jp/security/awareness/vendor/programmingv2/cc01.html
	fprintf(stderr,"axis=%d,	 2l+1=%d,	name=%s\n", s_axis, s_2lp1, s_name);
	char non_uniform[BUFSIZE], u2r2x_name[BUFSIZE], outer_name[BUFSIZE];	// buffer for fgets
	int L[4];//L[2] is redefined as L_2 + 2*outer_number; L[3] is outer_number.
	input_file2(Yee, L, non_uniform, u2r2x_name, outer_name);
	fprintf(stderr,"L0=%d, L1=%d, L2+2*outer=%d, outer=%d\n",L[0],L[1],L[2],L[3]);
	fprintf(stderr,"non_uniform=%s, u2r2x_name=%s, outer_name=%s\n",non_uniform,u2r2x_name,outer_name);
	if(s_axis == 2){	Case_u2(Yee, non_uniform, L, s_2lp1, s_name, u2r2x_name, outer_name);}
	else if(s_axis == 0){	Case_u0(Yee, non_uniform, L, s_2lp1, s_name, u2r2x_name, outer_name);}
	else if(s_axis == 1){	Case_u1(Yee, non_uniform, L, s_2lp1, s_name, u2r2x_name, outer_name);}
}
//======================================================================
//  input_file2 ---- rm_space, rm_comma                      
//                                       Last updated on Jul 31, 2018.  
//======================================================================
void input_file2(const char Yee[], int L[4], char non_uniform[BUFSIZE], char u2r2x_name[BUFSIZE], char outer_name[BUFSIZE])
{
	FILE *fp_i1;
	char file_name[BUFSIZE];
	snprintf(file_name, sizeof(file_name), "%s_MedB0.dat",Yee);// https://www.ipa.go.jp/security/awareness/vendor/programmingv1/b06_02.html
	fp_i1 = fopen(file_name,"r");
	if (fp_i1 == NULL){fprintf(stderr,"%s can not be read in output_file!\n",file_name);	exit(EXIT_FAILURE);}
	
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = -2;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 0) {
		char ini_file[BUFSIZE];
		if(strncmp(buf,"# info",6) == 0 && sscanf(buf,"%*[^=] %*[=] %s %*[^=] %*[=] %s", ini_file, non_uniform) == 2){
			rm_comma(ini_file); rm_comma(non_uniform);
			FILE *fp_i2 = fopen(ini_file,"r");
			int j_count2 = -2;
			while(fgets(buf, sizeof( buf ), fp_i2) != NULL && j_count2 < 0) {
				rm_space(buf);
				if(strncmp(buf,"u2r2x",5) == 0 && sscanf(buf,"%*[^=] %*[=] %s ", u2r2x_name) == 1){j_count2++;}
				else if(strncmp(buf,"outer",5) == 0 && sscanf(buf,"%*[^=] %*[=] %*[^=] %*[=] %s", outer_name) == 1){j_count2++;}
			}
			if(j_count2 < 0){fprintf(stderr,"file names can not be read in input_file2!\n");	exit(EXIT_FAILURE);}
			if(fclose(fp_i2) != 0) {fprintf(stderr,"fclose error2 in input_file2!\n");	exit(EXIT_FAILURE);}
			rm_comma(u2r2x_name); rm_comma(outer_name);
			j_count++;
		}else if(sscanf(buf,"%*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %d %*[^<] %*[<] %*[^=] %*[=] %d", &L[0],&L[1],&L[2],&L[3]) == 4){//Note L[3] means outer number not medium number!
			j_count++;
		}
	}
	if(fclose(fp_i1) != 0){fprintf(stderr,"fclose error1 in input_file2!\n");	exit(EXIT_FAILURE);}
	if(j_count < 0){fprintf(stderr,"Parameters can not be read in input_file2!\n");	exit(EXIT_FAILURE);}
}
//======================================================================
//  Case_u2 ---- get_u, get_u0r0         Last updated on Aug 29, 2018.  
//======================================================================
void get_u(const char *non_uniform, const int L[], const int axis_num, double u1[]);
void get_u0r0(const int L[], const int axis_num, const int s_2lp1, const char *file_name, double u0V[], double r0V[]);
void Case_u2(const char *Yee, const char *non_uniform, const int L[], const int s_2lp1, const char *s_name, const char *u2r2x_name, const char *outer_name)
{
	int l2 = s_2lp1/2;
	int s2 = s_2lp1 - 2*l2;
	char file_name[BUFSIZE];
	
	if(0 <= l2 && l2 < L[3]){	snprintf(file_name, sizeof(file_name), "%s_MedB%d.dat",Yee,l2);}
	else if(L[3] <= l2 && l2 < L[2] - L[3]){	snprintf(file_name, sizeof(file_name), "%s_Med%d.dat",Yee,l2 - L[3]);}
	else if(L[2] - L[3] <= l2 && l2 < L[2]){	snprintf(file_name, sizeof(file_name), "%s_MedT%d.dat",Yee,l2 - (L[2] - L[3]));}
	
	fprintf(stderr,"file name = %s, l2 = %d, s2 = %d in output_file\n",file_name, l2, s2);
	double *u0 = calloc((2*L[0]+1),sizeof(double));
	double *r0 = calloc((2*L[0]+1),sizeof(double));
	if(u0 == NULL || r0 == NULL){fprintf(stderr,"u0 and/or r0 can not be secured!\n");	exit(EXIT_FAILURE);}
	if(L[3] <= l2 && l2 < L[2] - L[3]){	get_u0r0(L,0, s_2lp1, u2r2x_name, u0, r0);}
	else{	get_u0r0(L,0, s_2lp1, outer_name, u0, r0);}
	
	FILE *fp_i1, *fpo_a, *fpo_b, *fpo_c;
	if ((fp_i1 = fopen(file_name,"r")) == NULL){fprintf(stderr,"%s can not be read in Case_u2!\n",file_name);	exit(EXIT_FAILURE);}
	char s_name_a[BUFSIZE], s_name_b[BUFSIZE], s_name_c[BUFSIZE];
	snprintf(s_name_a, sizeof(s_name_a), "%s_0.dat", s_name);	snprintf(s_name_b, sizeof(s_name_b), "%s_1.dat", s_name);	snprintf(s_name_c, sizeof(s_name_c), "%s_2.dat", s_name);
	if ((fpo_a  = fopen(s_name_a,"w")) == NULL){fprintf(stderr,"%s can not be written in Case_u2!\n",s_name_a);	exit(EXIT_FAILURE);}
	if ((fpo_b  = fopen(s_name_b,"w")) == NULL){fprintf(stderr,"%s can not be written in Case_u2!\n",s_name_b);	exit(EXIT_FAILURE);}
	if ((fpo_c  = fopen(s_name_c,"w")) == NULL){fprintf(stderr,"%s can not be written in Case_u2!\n",s_name_c);	exit(EXIT_FAILURE);}
	fprintf(fpo_a,"# Slice point of u2-axis: cell number = %d (with B and T outers = %d * 2), sub-cell number = %d \n",l2,L[3],s2);
	fprintf(fpo_b,"# Slice point of u2-axis: cell number = %d (with B and T outers = %d * 2), sub-cell number = %d \n",l2,L[3],s2);
	fprintf(fpo_c,"# Slice point of u2-axis: cell number = %d (with B and T outers = %d * 2), sub-cell number = %d \n",l2,L[3],s2);
	char buf[BUFSIZE];	// buffer for fgets
	{
		int j_count = -1;
		while(j_count < 0) {
			if(fgets(buf, sizeof( buf ), fp_i1) == NULL){fprintf(stderr,"comment1 read error!: %s\n",buf);	exit(EXIT_FAILURE);}
			if(strncmp(buf,"# k_0 * u_2",11) == 0){
				double us0, us1;
				if(sscanf(buf,"%*[^=] %*[=] %lf %*[;] %lf", &us0, &us1)!=2){fprintf(stderr,"# k_0 * u_2 read error!@case2: %s\n",buf);	exit(EXIT_FAILURE);};//revised on 20180829
				if(s2 == 0){fprintf(fpo_a,"# k0*u2 = %.5e\n",us0);	fprintf(fpo_b,"# k0*u2 = %.5e\n",us0);	fprintf(fpo_c,"# k0*u2 = %.5e\n",us0);}
				else if(s2 == 1){fprintf(fpo_a,"# k0*u2 = %.5e\n",us1);	fprintf(fpo_b,"# k0*u2 = %.5e\n",us1);	fprintf(fpo_c,"# k0*u2 = %.5e\n",us1);}
				j_count++;
			}
		}
		if(j_count < 0){fprintf(stderr,"Parameters can not be read in input_file2!\n");	exit(EXIT_FAILURE);}
	}
	while(strncmp(buf,"0 0 0",5) != 0) {
		if(fgets(buf, sizeof( buf ), fp_i1) == NULL){fprintf(stderr,"0 0 0 read error!: %s\n",buf);	exit(EXIT_FAILURE);}
	}
	int max_mls = (L[1]-1) + (L[0]-1) + 1;
	if(s2 == 0){
		fprintf(fpo_a,"#   k0*u0(1),        k0*r0,     k0*u1(0),        |mu|,      Arg(mu),     k0*fx00\n");//revised on 20180917
		fprintf(fpo_b,"#   k0*u0(0),        k0*r0,     k0*u1(1),        |mu|,      Arg(mu),     k0*fx11\n");//revised on 20180917
		fprintf(fpo_c,"#   k0*u0(1),        k0*r0,     k0*u1(1),      |epsi|,    Arg(epsi),     k0*fx22\n");//revised on 20180917
	}else if(s2 == 1){
		fprintf(fpo_a,"#   k0*u0(1),        k0*r0,     k0*u1(0),      |epsi|,    Arg(epsi),     k0*fx11\n");//revised on 20180917
		fprintf(fpo_b,"#   k0*u0(0),        k0*r0,     k0*u1(1),      |epsi|,    Arg(epsi),     k0*fx00\n");//revised on 20180917
		fprintf(fpo_c,"#   k0*u0(0),        k0*r0,     k0*u1(0),        |mu|,      Arg(mu),     k0*fx22\n");//revised on 20180917
	}
	double *u1 = calloc(2*L[1]+1,sizeof(double));
	get_u(non_uniform, L, 1,u1);
	for(int ma = 0 ; ma < L[1] ; ma++){
		for(int la = 0 ; la < L[0] ; la++){
			for(int sa = 0 ; sa < 2 ; sa++){
				int dl, dm, ds;
				double dr0, di0, dr1, di1, dr2, di2, fx_a = 0., fx_b = 0., fx_c = 0.;
				if(sscanf(buf,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf ", 
				&dl, &dm, &ds, &dr0, &di0, &dr1, &di1, &dr2, &di2, &fx_a, &fx_b, &fx_c) != 12){
					if(sscanf(buf,"%d %d %d %lf %lf %lf %lf %lf %lf ", 
					&dl, &dm, &ds, &dr0, &di0, &dr1, &di1, &dr2, &di2) != 9){	fprintf(stderr,"epsi read error!: %s\n",buf);	exit(EXIT_FAILURE);}
				}
				if(dl != la || dm != ma || ds != sa){fprintf(stderr,"int parameters can not be matched to loop!\n");	exit(EXIT_FAILURE);}
				if(fgets(buf, sizeof( buf ), fp_i1) == NULL && sa + la + ma != max_mls){fprintf(stderr,"End of file error! @(%d,%d,%d)\n",la,ma,sa);	exit(EXIT_FAILURE);}
				if(s2 == sa){//# l m 0 mu_00 mu_11 ep_22 fx_00 fx_11 fx_22; //# l m 1 ep_11 ep_00 mu_22 fx_11 fx_00 fx_22
					fprintf(fpo_a,"%+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",u0[ 2*la+1  ], r0[ 2*la+1  ], u1[ 2*ma+0  ], sqrt(dr0*dr0 + di0*di0), atan2(di0,dr0), fx_a);
					fprintf(fpo_b,"%+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",u0[ 2*la+0  ], r0[ 2*la+0  ], u1[ 2*ma+1  ], sqrt(dr1*dr1 + di1*di1), atan2(di1,dr1), fx_b);
					fprintf(fpo_c,"%+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",u0[2*la+1-s2], r0[2*la+1-s2], u1[2*ma+1-s2], sqrt(dr2*dr2 + di2*di2), atan2(di2,dr2), fx_c);
				}
			}
		}
		fprintf(fpo_a,"\n");	fprintf(fpo_b,"\n");	fprintf(fpo_c,"\n");
	}
	SAFEFREE(u1);	SAFEFREE(u0);	SAFEFREE(r0);
	if(fclose(fp_i1) != 0 || fclose(fpo_a) != 0 || fclose(fpo_b) != 0 || fclose(fpo_c) != 0){fprintf(stderr,"fclose error input or output in Case_u2!\n");	exit(EXIT_FAILURE);}
}
//======================================================================
//  get_u                                Last updated on Jul 31, 2018.  
//======================================================================
void get_u(const char *non_uniform, const int L[], const int axis_num, double u1[])
{//2*L[1]+1
	FILE *fp_i;
	if ((fp_i = fopen(non_uniform,"r")) == NULL){fprintf(stderr,"%s can not be read in output_file!\n",non_uniform);	exit(EXIT_FAILURE);}
	char buf[BUFSIZE];
	while(fgets(buf, sizeof( buf ), fp_i) != NULL){
		if(axis_num == 0 && strncmp(buf,"# xi0",5) == 0) {goto HeadData;}
		if(axis_num == 1 && strncmp(buf,"# xi1",5) == 0) {goto HeadData;}
	}
	fprintf(stderr,"# xi0 or xi1 can not be read in get_u!\n");	exit(EXIT_FAILURE);
	HeadData:;
	int L1x2p1 = 2*L[axis_num]+1;
	for(int i1 = 0 ; i1 < L1x2p1 ; i1++){
		if(fgets(buf, sizeof( buf ), fp_i) != NULL){
			double dummy;
			if(sscanf(buf,"%lf, %lf", 	&dummy, &u1[i1]) != 2){fprintf(stderr,"u1 can not be get in output_file!\n");	exit(EXIT_FAILURE);}
		}else{fprintf(stderr,"buf can not be get in output_file!\n");	exit(EXIT_FAILURE);}
	}
	if(fclose(fp_i) != 0){fprintf(stderr,"fclose error input in output_file!\n");	exit(EXIT_FAILURE);}
}
//======================================================================
//  get_u0r0                             Last updated on Jul 31, 2018.  
//======================================================================
void get_u0r0(const int L[], const int axis_num, const int s_2lp1, const char *file_name, double uV[], double rV[])
{
	int l2x2p1 = s_2lp1 - (2*L[3]);
	FILE *fp_i;
	if ((fp_i = fopen(file_name,"r")) == NULL){fprintf(stderr,"%s can not be read in get_u0r0!\n",file_name);	exit(EXIT_FAILURE);}
	char buf[BUFSIZE];	// buffer for fgets
	while(fgets(buf, sizeof( buf ), fp_i) != NULL) {
		int dummy;
		if(sscanf(buf,"%*[^=] %*[=] %d", &dummy) == 1) {
			if(dummy == l2x2p1 && fgets(buf, sizeof( buf ), fp_i) != NULL){
				if(strncmp(buf,"#2l_0",5) == 0){	goto HeadData;}
			}
		}
	}
	fprintf(stderr,"Data head can not be read in get_u0r0!: 2l2+1=%d, file_name=%s\n",l2x2p1,file_name);
	fprintf(stderr,"L[0]=%d, L[1]=%d, L[2]=%d, L[3]=%d\n",L[0],L[1],L[2],L[3]);
	exit(1);
	HeadData:;
//#2l_0 + 1,            k0*u_2,                     k0*u_0,                     k0*r_2,                     k0*r_0
	for(int i0 = 0 ; i0 < L[0]*2 + 1 ; i0++){
		int dummy;
		double u2, r2, u0, r0;
		if(fgets(buf, sizeof( buf ), fp_i) != NULL){
			if(sscanf(buf,"%d, %lf, %lf, %lf, %lf", &dummy, &u2, &u0, &r2, &r0) == 5) {
				if(dummy != i0){fprintf(stderr,"Error! %d != i0 = %d in get_file_info!\n", dummy, i0);	exit(1);}
			}else{fprintf(stderr,"5 data can not be read in get_u0r0!\n");	exit(1);}
		}else{fprintf(stderr,"Data can not be read in get_u0r0!\n");	exit(1);}
		if(axis_num == 0){
			uV[i0] = u0;
			rV[i0] = r0;
		}else if(axis_num == 2){
			uV[i0] = u2;
			rV[i0] = r2;
		}else{fprintf(stderr,"Error! axis num = %d\n", axis_num);	exit(1);}
	}
	if(fclose(fp_i) != 0){fprintf(stderr,"fclose error input in get_u0r0!\n");	exit(EXIT_FAILURE);}
}
//======================================================================
//  Case_u0 ---- get_u0r0, get_u         Last updated on Aug 07, 2018.  
//======================================================================
void Case_u0(const char *Yee, const char *non_uniform, const int L[], const int s_2lp1, const char *s_name, const char *u2r2x_name, const char *outer_name)
{
	if( s_2lp1 < 0 || s_2lp1 >= 2*L[0]+1){fprintf(stderr,"Error in Case_u0: s_2lp1 = %d\n", s_2lp1);	exit(EXIT_FAILURE);}
	double *u2 = calloc((2*L[2]+1),sizeof(double));
	double *r2 = calloc((2*L[2]+1),sizeof(double));
	for(int l2x2 = 0 ; l2x2 < 2*L[2]+1; l2x2++){//l2x2 = 2*L[2] is not used. commented on 20180731 
		double *uV2 = calloc((2*L[0]+1),sizeof(double));
		double *rV2 = calloc((2*L[0]+1),sizeof(double));
		if(uV2 == NULL || rV2 == NULL){fprintf(stderr,"u0 and/or r0 can not be secured!\n");	exit(EXIT_FAILURE);}
		if(2*L[3] <= l2x2 && l2x2 <= 2*L[2] - 2*L[3]){
			get_u0r0(L,2, l2x2, u2r2x_name, uV2, rV2);
		}else{
			get_u0r0(L,2, l2x2, outer_name, uV2, rV2);
		}
		u2[l2x2] = uV2[s_2lp1];
		r2[l2x2] = rV2[s_2lp1];
		if(uV2 != NULL){free(uV2);	uV2 = NULL;}	// To avoid "dangling pointer" and "double-free"
		if(rV2 != NULL){free(rV2);	rV2 = NULL;}	// To avoid "dangling pointer" and "double-free"
	}
	int l0 = s_2lp1/2;
	int s0 = s_2lp1 - 2*l0;
	FILE *fpo_a, *fpo_b, *fpo_c;
	char s_name_a[BUFSIZE], s_name_b[BUFSIZE], s_name_c[BUFSIZE];
	snprintf(s_name_a, sizeof(s_name_a), "%s_0.dat", s_name);	snprintf(s_name_b, sizeof(s_name_b), "%s_1.dat", s_name);	snprintf(s_name_c, sizeof(s_name_c), "%s_2.dat", s_name);
	if ((fpo_a  = fopen(s_name_a,"w")) == NULL){fprintf(stderr,"%s can not be written in Case_u0!\n",s_name_a);	exit(EXIT_FAILURE);}
	if ((fpo_b  = fopen(s_name_b,"w")) == NULL){fprintf(stderr,"%s can not be written in Case_u0!\n",s_name_b);	exit(EXIT_FAILURE);}
	if ((fpo_c  = fopen(s_name_c,"w")) == NULL){fprintf(stderr,"%s can not be written in Case_u0!\n",s_name_c);	exit(EXIT_FAILURE);}
	fprintf(fpo_a,"# Slice point of u0-axis: cell number(l_0) = %d, sub-cell number = %d\n",l0 ,s0);
	fprintf(fpo_b,"# Slice point of u0-axis: cell number(l_0) = %d, sub-cell number = %d\n",l0 ,s0);
	fprintf(fpo_c,"# Slice point of u0-axis: cell number(l_0) = %d, sub-cell number = %d\n",l0 ,s0);
	double *u0 = calloc(2*L[0]+1,sizeof(double));
	get_u(non_uniform, L, 0,u0);
	fprintf(fpo_a,"# k0*u0 = %.5e\n", u0[s_2lp1]);
	fprintf(fpo_b,"# k0*u0 = %.5e\n", u0[s_2lp1]);
	fprintf(fpo_c,"# k0*u0 = %.5e\n", u0[s_2lp1]);
	SAFEFREE(u0);
	if(s0 == 0){
		fprintf(fpo_a,"#   k0*u1(1),     k0*u2(0),        k0*r2,        |mu|,      Arg(mu),     k0*fx11\n");//revised on 20180917
		fprintf(fpo_b,"#   k0*u1(0),     k0*u2(1),        k0*r2,        |mu|,      Arg(mu),     k0*fx22\n");//revised on 20180917
		fprintf(fpo_c,"#   k0*u1(1),     k0*u2(1),        k0*r2,      |epsi|,    Arg(epsi),     k0*fx00\n");//revised on 20180917
	}else if(s0 == 1){
		fprintf(fpo_a,"#   k0*u1(1),     k0*u2(0),        k0*r2,      |epsi|,    Arg(epsi),     k0*fx22\n");//revised on 20180917
		fprintf(fpo_b,"#   k0*u1(0),     k0*u2(1),        k0*r2,      |epsi|,    Arg(epsi),     k0*fx11\n");//revised on 20180917
		fprintf(fpo_c,"#   k0*u1(0),     k0*u2(0),        k0*r2,        |mu|,      Arg(mu),     k0*fx00\n");//revised on 20180917
	}
	double *u1 = calloc(2*L[1]+1,sizeof(double));
	get_u(non_uniform, L, 1,u1);
	for(int l2 = 0 ; l2 < L[2]; l2++){
		char file_name[BUFSIZE];
		if(0 <= l2 && l2 < L[3]){	snprintf(file_name, sizeof(file_name), "%s_MedB%d.dat",Yee,l2);}
		else if(L[3] <= l2 && l2 < L[2] - L[3]){snprintf(file_name, sizeof(file_name), "%s_Med%d.dat",Yee,l2 - L[3]);}
		else if(L[2] - L[3] <= l2 && l2 < L[2]){snprintf(file_name, sizeof(file_name), "%s_MedT%d.dat",Yee,l2 - (L[2] - L[3]));}
		FILE *fp_i1;
		if ((fp_i1 = fopen(file_name,"r")) == NULL){fprintf(stderr,"%s can not be read in Case_u0!\n",file_name);	exit(EXIT_FAILURE);}
		char buf[BUFSIZE];	// buffer for fgets
		{
			int j_count = -1;
			while(j_count < 0) {
				if(fgets(buf, sizeof( buf ), fp_i1) == NULL){fprintf(stderr,"comment1 read error!: %s\n",buf);	exit(EXIT_FAILURE);}
				if(strncmp(buf,"# k_0 * u_2",11) == 0){/*fprintf(fp_o,"%s",buf);*/	j_count++;}
			}
			if(j_count < 0){fprintf(stderr,"Parameters can not be read in input_file2!\n");	exit(EXIT_FAILURE);}
		}
		while(strncmp(buf,"0 0 0",5) != 0) {
			if(fgets(buf, sizeof( buf ), fp_i1) == NULL){fprintf(stderr,"0 0 0 read error!: %s\n",buf);	exit(EXIT_FAILURE);}
		}
		int max_mls = (L[1]-1) + (L[0]-1) + 1; 
		for(int ma = 0 ; ma < L[1] ; ma++){
			for(int la = 0 ; la < L[0] ; la++){
				for(int sa = 0 ; sa < 2 ; sa++){
					int dl, dm, ds;
					double dr0, di0, dr1, di1, dr2, di2, fx_a = 0., fx_b = 0., fx_c = 0.;
					if(sscanf(buf,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf ", 
					&dl, &dm, &ds, &dr0, &di0, &dr1, &di1, &dr2, &di2, &fx_a, &fx_b, &fx_c) != 12){
						if(sscanf(buf,"%d %d %d %lf %lf %lf %lf %lf %lf ", 
						&dl, &dm, &ds, &dr0, &di0, &dr1, &di1, &dr2, &di2) != 9){	fprintf(stderr,"epsi read error!: %s\n",buf);	exit(EXIT_FAILURE);}
					}
					if(dl != la || dm != ma || ds != sa){fprintf(stderr,"int parameters can not be matched to loop!\n");	exit(EXIT_FAILURE);}
					if(fgets(buf, sizeof( buf ), fp_i1) == NULL && sa + la + ma != max_mls){fprintf(stderr,"End of file error! @(%d,%d,%d)\n",la,ma,sa);	exit(EXIT_FAILURE);}
					if(la == l0){//# l m 0 mu_00 mu_11 ep_22 fx_00 fx_11 fx_22; # l m 1 ep_11 ep_00 mu_22 fx_11 fx_00 fx_22
						if(s0 == 0){
							if(sa == 0){fprintf(fpo_a,"%+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",//# l m 0 ***** mu_11 ***** +++++ fx_11 +++++
												u1[2*ma+1], u2[2*l2+0], r2[2*l2+0], sqrt(dr1*dr1 + di1*di1), atan2(di1,dr1), fx_b);}
							if(sa == 1){fprintf(fpo_b,"%+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",//# l m 1 ***** ***** mu_22 +++++ +++++ fx_22
												u1[2*ma+0], u2[2*l2+1], r2[2*l2+1], sqrt(dr2*dr2 + di2*di2), atan2(di2,dr2), fx_c);}
							if(sa == 1){fprintf(fpo_c,"%+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",//# l m 1 ***** ep_00 ***** +++++ fx_00 +++++
												u1[2*ma+1], u2[2*l2+1], r2[2*l2+1], sqrt(dr1*dr1 + di1*di1), atan2(di1,dr1), fx_b);
							}
						}else if(s0 == 1){
							if(sa == 0){fprintf(fpo_a,"%+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",//# l m 0 ***** ***** ep_22 +++++ +++++ fx_22
												u1[2*ma+1], u2[2*l2+0], r2[2*l2+0], sqrt(dr2*dr2 + di2*di2), atan2(di2,dr2), fx_c);}
							if(sa == 1){fprintf(fpo_b,"%+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",//# l m 1 ep_11 ***** ***** fx_11 +++++ +++++
												u1[2*ma+0], u2[2*l2+1], r2[2*l2+1], sqrt(dr0*dr0 + di0*di0), atan2(di0,dr0), fx_a);}
							if(sa == 0){fprintf(fpo_c,"%+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",//# l m 0 mu_00 ***** ***** fx_00 +++++ +++++
												u1[2*ma+0], u2[2*l2+0], r2[2*l2+0], sqrt(dr0*dr0 + di0*di0), atan2(di0,dr0), fx_a);}
						}
					}
				}
			}
		}
		if(fclose(fp_i1) != 0){fprintf(stderr,"fclose error input in Case_u0!\n");	exit(EXIT_FAILURE);}
		fprintf(fpo_a,"\n");	fprintf(fpo_b,"\n");	fprintf(fpo_c,"\n");
	}
	SAFEFREE(u1);	SAFEFREE(u2);	SAFEFREE(r2);
	if(fclose(fpo_a) != 0 || fclose(fpo_b) != 0 || fclose(fpo_c) != 0){fprintf(stderr,"fclose error output in Case_u0!\n");	exit(EXIT_FAILURE);}
}
//======================================================================
//  Case_u1 ---- get_u2x2, get_u         Last updated on Aug 08, 2018.  
//======================================================================
void get_u2x2(const int L[], const int s_2lp1, const char *file_name, double uV2[], double rV2[], double uV0[], double rV0[]);
void Case_u1(const char *Yee, const char *non_uniform, const int L[], const int s_2lp1, const char *s_name, const char *u2r2x_name, const char *outer_name)
{
	if( s_2lp1 < 0 || s_2lp1 >= 2*L[1]+1){fprintf(stderr,"Error in Case_u1: s_2lp1 = %d\n", s_2lp1);	exit(EXIT_FAILURE);}
	double *u2 = calloc((2*L[2]+1)*(2*L[0]+1),sizeof(double));
	double *u0 = calloc((2*L[2]+1)*(2*L[0]+1),sizeof(double));
	double *x2 = calloc((2*L[2]+1)*(2*L[0]+1),sizeof(double));
	double *x0 = calloc((2*L[2]+1)*(2*L[0]+1),sizeof(double));
	for(int l2x2 = 0 ; l2x2 < 2*L[2]+1; l2x2++){
		int shift = l2x2*(2*L[0]+1);
		if(2*L[3] <= l2x2 && l2x2 <= 2*L[2] - 2*L[3]){
			get_u2x2(L, l2x2, u2r2x_name, &u2[shift], &x2[shift], &u0[shift], &x0[shift]);//revised on 20180801
		}else{
			get_u2x2(L, l2x2, outer_name, &u2[shift], &x2[shift], &u0[shift], &x0[shift]);//revised on 20180801
		}
	}
	
	FILE *fpo_a, *fpo_b, *fpo_c;
	char s_name_a[BUFSIZE], s_name_b[BUFSIZE], s_name_c[BUFSIZE];
	snprintf(s_name_a, sizeof(s_name_a), "%s_0.dat", s_name);	snprintf(s_name_b, sizeof(s_name_b), "%s_1.dat", s_name);	snprintf(s_name_c, sizeof(s_name_c), "%s_2.dat", s_name);
	if ((fpo_a  = fopen(s_name_a,"w")) == NULL){fprintf(stderr,"%s can not be written in Case_u1!\n",s_name_a);	exit(EXIT_FAILURE);}
	if ((fpo_b  = fopen(s_name_b,"w")) == NULL){fprintf(stderr,"%s can not be written in Case_u1!\n",s_name_b);	exit(EXIT_FAILURE);}
	if ((fpo_c  = fopen(s_name_c,"w")) == NULL){fprintf(stderr,"%s can not be written in Case_u1!\n",s_name_c);	exit(EXIT_FAILURE);}
	
	
	
	
	int l1 = s_2lp1/2;
	int s1 = s_2lp1 - 2*l1;
	fprintf(fpo_a,"# Slice point of u1-axis: cell number(l_1) = %d, sub-cell number = %d\n",l1 ,s1);
	fprintf(fpo_b,"# Slice point of u1-axis: cell number(l_1) = %d, sub-cell number = %d\n",l1 ,s1);
	fprintf(fpo_c,"# Slice point of u1-axis: cell number(l_1) = %d, sub-cell number = %d\n",l1 ,s1);
	double *u1 = calloc(2*L[1]+1,sizeof(double));
	get_u(non_uniform, L, 1,u1);
	fprintf(fpo_a,"# k0*u1 = %.5e\n", u1[s_2lp1]);
	fprintf(fpo_b,"# k0*u1 = %.5e\n", u1[s_2lp1]);
	fprintf(fpo_c,"# k0*u1 = %.5e\n", u1[s_2lp1]);
	SAFEFREE(u1);
	if(s1 == 0){//     #   k0*u2(1),        k0*x2,     k0*u0(0),        k0*x0,        |mu|,        k0*x0,        k0*x0
		fprintf(fpo_a,"#   k0*u2(1),        k0*x2,     k0*u0(0),        k0*x0,        |mu|,      Arg(mu),     k0*fx22\n");//revised on 20180917
		fprintf(fpo_b,"#   k0*u2(0),        k0*x2,     k0*u0(1),        k0*x0,        |mu|,      Arg(mu),     k0*fx00\n");//revised on 20180917
		fprintf(fpo_c,"#   k0*u2(1),        k0*x2,     k0*u0(1),        k0*x0,      |epsi|,    Arg(epsi),     k0*fx11\n");//revised on 20180917
	}else if(s1 == 1){
		fprintf(fpo_a,"#   k0*u2(1),        k0*x2,     k0*u0(0),        k0*x0,      |epsi|,    Arg(epsi),     k0*fx00\n");//revised on 20180917
		fprintf(fpo_b,"#   k0*u2(0),        k0*x2,     k0*u0(1),        k0*x0,      |epsi|,    Arg(epsi),     k0*fx22\n");//revised on 20180917
		fprintf(fpo_c,"#   k0*u2(0),        k0*x2,     k0*u0(0),        k0*x0,        |mu|,      Arg(mu),     k0*fx11\n");//revised on 20180917
	}
	for(int l2 = 0 ; l2 < L[2]; l2++){
		char file_name[BUFSIZE];
		if(0 <= l2 && l2 < L[3]){	snprintf(file_name, sizeof(file_name), "%s_MedB%d.dat",Yee,l2);}
		else if(L[3] <= l2 && l2 < L[2] - L[3]){snprintf(file_name, sizeof(file_name), "%s_Med%d.dat",Yee,l2 - L[3]);}
		else if(L[2] - L[3] <= l2 && l2 < L[2]){snprintf(file_name, sizeof(file_name), "%s_MedT%d.dat",Yee,l2 - (L[2] - L[3]));}
		FILE *fp_i1;
		if ((fp_i1 = fopen(file_name,"r")) == NULL){fprintf(stderr,"%s can not be read in Case_u1!\n",file_name);	exit(EXIT_FAILURE);}
		char buf[BUFSIZE];	// buffer for fgets
		{
			int j_count = -1;
			while(j_count < 0) {
				if(fgets(buf, sizeof( buf ), fp_i1) == NULL){fprintf(stderr,"comment1 read error!: %s\n",buf);	exit(EXIT_FAILURE);}
				if(strncmp(buf,"# k_0 * u_2",11) == 0){/*fprintf(fp_o,"%s",buf);*/	j_count++;}
			}
			if(j_count < 0){fprintf(stderr,"Parameters can not be read in Case_u1!\n");	exit(EXIT_FAILURE);}
		}
		while(strncmp(buf,"0 0 0",5) != 0) {
			if(fgets(buf, sizeof( buf ), fp_i1) == NULL){fprintf(stderr,"0 0 0 read error!: %s\n",buf);	exit(EXIT_FAILURE);}
		}
		int max_mls = (L[1]-1) + (L[0]-1) + 1; 
		for(int ma = 0 ; ma < L[1] ; ma++){
			for(int la = 0 ; la < L[0] ; la++){
				int mat00 = 2*la+0+(2*l2+0)*(2*L[0]+1);
				int mat01 = 2*la+0+(2*l2+1)*(2*L[0]+1);
				int mat10 = 2*la+1+(2*l2+0)*(2*L[0]+1);
				int mat11 = 2*la+1+(2*l2+1)*(2*L[0]+1);
				for(int sa = 0 ; sa < 2 ; sa++){
					int dl, dm, ds;
					double dr0, di0, dr1, di1, dr2, di2, fx_a = 0., fx_b = 0., fx_c = 0.;
					if(sscanf(buf,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf ", 
					&dl, &dm, &ds, &dr0, &di0, &dr1, &di1, &dr2, &di2, &fx_a, &fx_b, &fx_c) != 12){
						if(sscanf(buf,"%d %d %d %lf %lf %lf %lf %lf %lf ", 
						&dl, &dm, &ds, &dr0, &di0, &dr1, &di1, &dr2, &di2) != 9){	fprintf(stderr,"epsi read error!: %s\n",buf);	exit(EXIT_FAILURE);}
					}
					if(dl != la || dm != ma || ds != sa){fprintf(stderr,"int parameters can not be matched to loop!\n");	exit(EXIT_FAILURE);}
					if(fgets(buf, sizeof( buf ), fp_i1) == NULL && sa + la + ma != max_mls){fprintf(stderr,"End of file error! @(%d,%d,%d)\n",la,ma,sa);	exit(EXIT_FAILURE);}
					if(ma == l1){//# l m 0 mu_00 mu_11 ep_22 fx_00 fx_11 fx_22; # l m 1 ep_11 ep_00 mu_22 fx_11 fx_00 fx_22
						if(s1 == 0){
							if(sa == 1){fprintf(fpo_a,"%+.5e, %+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",//# l m 1 ***** ***** mu_22 +++++ +++++ fx_22
												u2[mat01], x2[mat01], u0[mat01], x0[mat01], sqrt(dr2*dr2 + di2*di2), atan2(di2,dr2), fx_c);}
							if(sa == 0){fprintf(fpo_b,"%+.5e, %+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",//# l m 0 mu_00 ***** ***** fx_00 +++++ +++++
												u2[mat10], x2[mat10], u0[mat10], x0[mat10], sqrt(dr0*dr0 + di0*di0), atan2(di0,dr0), fx_a);}
							if(sa == 1){fprintf(fpo_c,"%+.5e, %+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",//# l m 1 ep_11 ***** ***** fx_11 +++++ +++++
												u2[mat11], x2[mat11], u0[mat11], x0[mat11], sqrt(dr0*dr0 + di0*di0), atan2(di0,dr0), fx_a);}
						}else if(s1 == 1){
							if(sa == 1){fprintf(fpo_a,"%+.5e, %+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",//# l m 1 ***** ep_00 ***** +++++ fx_00 +++++
												u2[mat01], x2[mat01], u0[mat01], x0[mat01], sqrt(dr1*dr1 + di1*di1), atan2(di1,dr1), fx_b);}
							if(sa == 0){fprintf(fpo_b,"%+.5e, %+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",//# l m 0 ***** ***** ep_22 +++++ +++++ fx_22
												u2[mat10], x2[mat10], u0[mat10], x0[mat10], sqrt(dr2*dr2 + di2*di2), atan2(di2,dr2), fx_c);}
							if(sa == 0){fprintf(fpo_c,"%+.5e, %+.5e, %+.5e, %+.5e, %.5e, %+.5e, %.5e\n",//# l m 0 ***** mu_11 ***** +++++ fx_11 +++++
												u2[mat00], x2[mat00], u0[mat00], x0[mat00], sqrt(dr1*dr1 + di1*di1), atan2(di1,dr1), fx_b);}
						}
					}
				}
			}
		}
		if(fclose(fp_i1) != 0){fprintf(stderr,"fclose error input in Case_u1!\n");	exit(EXIT_FAILURE);}
		fprintf(fpo_a,"\n");	fprintf(fpo_b,"\n");	fprintf(fpo_c,"\n");
	}
	if(fclose(fpo_a) != 0 || fclose(fpo_b) != 0 || fclose(fpo_c) != 0){fprintf(stderr,"fclose error input or output in Case_u1!\n");	exit(EXIT_FAILURE);}
	SAFEFREE(u2);	SAFEFREE(x2);	SAFEFREE(u0);	SAFEFREE(x0);
}
//======================================================================
//  get_u2x2                             Last updated on Jul 31, 2018.  
//======================================================================
void get_u2x2(const int L[], const int s_2lp1, const char *file_name, double u2[], double x2[], double u0[], double x0[])
{
	int l2x2p1 = s_2lp1 - (2*L[3]);
	FILE *fp_i;
	if ((fp_i = fopen(file_name,"r")) == NULL){fprintf(stderr,"%s can not be read in get_u0r0!\n",file_name);	exit(EXIT_FAILURE);}
	char buf[BUFSIZE];	// buffer for fgets
	while(fgets(buf, sizeof( buf ), fp_i) != NULL) {
		int dummy;
		if(sscanf(buf,"%*[^=] %*[=] %d", &dummy) == 1) {
			if(dummy == l2x2p1 && fgets(buf, sizeof( buf ), fp_i) != NULL){
				if(strncmp(buf,"#2l_0",5) == 0){	goto HeadData;}
			}
		}
	}
	fprintf(stderr,"Data head can not be read in get_u2r2!: 2l2+1=%d, file_name=%s\n",l2x2p1,file_name);
	fprintf(stderr,"L[0]=%d, L[1]=%d, L[2]=%d, L[3]=%d\n",L[0],L[1],L[2],L[3]);
	exit(1);
	HeadData:;
//#2l_0 + 1,            k0*u_2,                     k0*u_0,                     k0*r_2,                     k0*r_0,                     k0*x_2,                     k0*x_0,                        h_2,                        h_0
	for(int i0 = 0 ; i0 < L[0]*2 + 1 ; i0++){
		int dummy;
		double r2, r0;
		if(fgets(buf, sizeof( buf ), fp_i) != NULL){
			if(sscanf(buf,"%d, %lf, %lf, %lf, %lf, %lf, %lf", &dummy, &u2[i0], &u0[i0], &r2, &r0, &x2[i0], &x0[i0]) == 7) {
				if(dummy != i0){fprintf(stderr,"Error! %d != i0 = %d in get_file_info!\n", dummy, i0);	exit(1);}
			}else{fprintf(stderr,"5 data can not be read in get_u0r0!\n");	exit(1);}
		}else{fprintf(stderr,"Data can not be read in get_u0r0!\n");	exit(1);}
	}
	if(fclose(fp_i) != 0){fprintf(stderr,"fclose error input in get_u0r0!\n");	exit(EXIT_FAILURE);}
}
