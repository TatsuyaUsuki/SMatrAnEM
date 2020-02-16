//  Start of StrMap.h
#include "header_macro.h"
#include "constant.h"
int input_file1(FILE *fp_i1, int L[3], char xi_file[], double *k_0, char str_file[], int *RecDep, char Yee[], int *cell_num, char Region[]);
void input_file2(FILE *fp_i1, int L, double u0[L*2+1], double du0dxi0[L*2+1]);
int input_N_str(FILE *fp_i3);
void input_file3(FILE *fp_i3, const double k_0, const int N_str, 
				int media_no[N_str], int shape_str[N_str], double shape_length[N_str*6], 
				double origin_str[N_str*3], double array_step[N_str*9], int n_step[N_str*6]);
void output_file(char **argv, const char Yee[], 
			const int RecDep, const char str_file[], const int L[3], 
			const double u0[L[0]*2+1], const double u1[L[1]*2+1], const double u2[L[2]*2+1], 
			const int N_str, const int media_no[N_str], const int shape_str[N_str], 
			const double shape_length[N_str*6], const double origin_str[N_str*3], 
			const double array_step[N_str*9], const int n_step[N_str*6], const char xi_file[], const double k_0);
//  End of StrMap.h
//======================================================================
//   main ---- input_file1, input_file2, input_file3, output_file       
//======================================================================
int main(int argc, char **argv)
{
	if(argc != 2) {
		fprintf(stderr,"error: number of files \n");
		exit(1);
	}else if(strncmp(argv[1], "-v", 2) == 0 || strcmp(argv[1], "--version") == 0 ) {
		fprintf(stderr,"The '%s' creates media-infomation in each Yee's cell.\n", argv[0]);
		fprintf(stderr,"Version 18.08.03 is compiled at %s on %s.\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(stderr," Source code : '%s'\n Author      : Tatsuya Usuki\n URL         : http://www.smatran.org\n", __FILE__);
		fprintf(stderr," References  : 'Wave scattering in frequency domain' as 'Formulation.pdf' on May 20, 2018;\n");
		fprintf(stderr,"                 Section 8.2 in 'Formulation for SMatrAn' as 'manual.pdf' at Jan 16, 2017.\n");
		fprintf(stderr,"There is NO warranty.\n");
		exit(0);//normal end
	}
//---------- begin reading Yee lattice information ------------
	time_t timer_s = time(0);	clock_t startClock = clock();
	FILE *fp_i1;
	fp_i1 = fopen(argv[1],"r");
	if (fp_i1 == NULL){
		fprintf(stderr,"open error!: open input-file1!\n");
		exit(1);
	}
	fprintf(stderr,"Mesh data file (the 1st input file): %s\n",argv[1]);
	int L[3], RecDep, cell_num; //division number: L0 = L[0], L1 = L[1], L2 = L[2], Recursive depth
	char xi_file[BUFSIZE], str_file[BUFSIZE], Yee[BUFSIZE], Region[BUFSIZE];
	double k_0;
	if(input_file1(fp_i1, L, xi_file, &k_0, str_file, &RecDep, Yee, &cell_num, Region) != 0) { 
		fprintf(stderr,"input_file1 error!\n");
		exit(1);
	}
	fprintf(stderr,"k0 [m^-1] = %.5e\n", k_0);
	fprintf(stderr,"Yee's lattice output file `%s*.dat'\n", Yee);
	fprintf(stderr,"Recursive depth = %d\n", RecDep);
	if(fclose(fp_i1) != 0) {
		fprintf(stderr,"fclose error after input_file!\n");
		exit(1);
	}
	FILE *fp_i2;
	fp_i2 = fopen(xi_file,"r");
	if (fp_i2 == NULL){
		fprintf(stderr,"open error!: open xi_file!\n");
		exit(1);
	}
	double *u0 = calloc((L[0]*2+1),sizeof(double));
	double *du0dxi0 = calloc((L[0]*2+1),sizeof(double));
	if(u0 == NULL || du0dxi0 == NULL){fprintf(stderr,"u0 and/or du0dxi0 can not be secured!\n");	exit(EXIT_FAILURE);}
	input_file2(fp_i2, L[0], u0, du0dxi0);
	
	double *u1 = calloc((L[1]*2+1),sizeof(double));
	double *du1dxi1 = calloc((L[1]*2+1),sizeof(double));
	if(u1 == NULL || du1dxi1 == NULL){fprintf(stderr,"u1 and/or du1dxi1 can not be secured!\n");	exit(EXIT_FAILURE);}
	input_file2(fp_i2, L[1], u1, du1dxi1);
	
	double *u2 = calloc((L[2]*2+1),sizeof(double));
	double *du2dxi2 = calloc((L[2]*2+1),sizeof(double));
	if(u2 == NULL || du2dxi2 == NULL){fprintf(stderr,"u2 and/or du2dxi2 can not be secured!\n");	exit(EXIT_FAILURE);}
	input_file2(fp_i2, L[2], u2, du2dxi2);
	if(fclose(fp_i2) != 0) {
		fprintf(stderr,"fclose error after xi_file!\n");
		exit(1);
	}
	for(int j0 = 0 ; j0 < 3 ; j0++) { 
		if (j0 == 0) fprintf(stderr,"\nL, k0*u_0, k0*du0dxi0\n");
		if (j0 == 1) fprintf(stderr,"\nM, k0*u_1, k0*du1dxi1\n");
		if (j0 == 2) fprintf(stderr,"\nN, k0*u_2, k0*du2dxi2\n");
		for(int l0 = 0 ; l0 < 2*L[j0]+1 ; l0++) { 
			if (j0 == 0) fprintf(stderr,"%d, %.5E, %.5E\n", l0, u0[l0], du0dxi0[l0]);
			if (j0 == 1) fprintf(stderr,"%d, %.5E, %.5E\n", l0, u1[l0], du1dxi1[l0]);
			if (j0 == 2) fprintf(stderr,"%d, %.5E, %.5E\n", l0, u2[l0], du2dxi2[l0]);
		}
	}
	SAFEFREE(du0dxi0);	SAFEFREE(du1dxi1);
//----------  end reading Yee lattive information  ------------
//---------- begin reading 3D-structure information ------------
	FILE *fp_i3;
	fp_i3 = fopen(str_file,"r");
	if (fp_i3 == NULL){
		fprintf(stderr,"open error!: open str_file!\n");
		exit(1);
	}
	fprintf(stderr,"\n3D-structure data file: %s\n",str_file);
	
	int N_str = input_N_str(fp_i3); //total number of structures
	rewind(fp_i3);
	
	fprintf(stderr,"N_str = %d", N_str);
	int shape_str[N_str], media_no[N_str];
	double shape_length[N_str*6], origin_str[N_str*3];
	double array_step[N_str*9];
	int n_step[N_str*6];
	
	input_file3(fp_i3, k_0, N_str, media_no, shape_str, shape_length, origin_str, array_step, n_step);
	if(fclose(fp_i3) != 0) {
		fprintf(stderr,"\nfclose error after input_file3\n");
		exit(1);
	}
	int Max_media_no = media_no[0];
	for(int j0 = 1 ; j0 < N_str ; j0++) {
		if(Max_media_no < media_no[j0]) Max_media_no = media_no[j0];
	}
	fprintf(stderr,", media_no < %d\n", Max_media_no + 1);
	for(int j0 = 0 ; j0 < N_str ; j0++) { 
		fprintf(stderr,"\nno_str = %d\n", j0);
		fprintf(stderr,"media_no = %d\n", media_no[j0]);
		if (shape_str[j0] < 0 || shape_str[j0] > 4) fprintf(stderr,"shape_str error!\n");
		if (shape_str[j0] == 0) fprintf(stderr,"shape_str = block\n");
		if (shape_str[j0] == 1) fprintf(stderr,"shape_str = cylinder\n");
		if (shape_str[j0] == 2) fprintf(stderr,"shape_str = prism\n");
		if (shape_str[j0] == 3) fprintf(stderr,"shape_str = pyramid\n");
		if (shape_str[j0] == 4) fprintf(stderr,"shape_str = spheroid\n");
		fprintf(stderr,"shape_length = (%.5E, %.5E, %.5E)\n", shape_length[j0*6+0], shape_length[j0*6+1], shape_length[j0*6+2]);
		fprintf(stderr,"shape_angle = (%.5E, %.5E, %.5E)\n", shape_length[j0*6+3], shape_length[j0*6+4], shape_length[j0*6+5]);
		fprintf(stderr,"origin_str = (%.5E, %.5E, %.5E)\n", origin_str[j0*3+0], origin_str[j0*3+1], origin_str[j0*3+2]);
		for(int i0 = 0 ; i0 < 3 ; i0++){
			fprintf(stderr,"Step %d:\n", i0+1);
			fprintf(stderr,"array_step = (%.5E, %.5E, %.5E)\n", 
				array_step[j0*9+0+i0*3], array_step[j0*9+1+i0*3], array_step[j0*9+2+i0*3]);
			fprintf(stderr,"n_step start = %d, finish = %d\n", n_step[j0*6+0+i0*2], n_step[j0*6+1+i0*2]);
		}
	}
//----------  end reading 3D-structure information  ------------
	if(strcmp(Region, "all") == 0){
		output_file(argv, Yee, RecDep, str_file, L, 
			u0, u1, u2, 
			N_str, media_no, 
			shape_str, shape_length, 
			origin_str, array_step, n_step, xi_file, k_0);//revised on 20180801
	}
	if(strcmp(Region, "outer") == 0 || strcmp(Region, "all") == 0){
		int Ld[3]={L[0],L[1],cell_num};
		char Yeed[BUFSIZE];
		// Prepare Bottom region.
		snprintf(Yeed,sizeof(Yeed),"%sB", Yee);//revised on 20180723 https://www.ipa.go.jp/security/awareness/vendor/programmingv1/b06_02.html
		double *u2d = calloc((cell_num*2+1),sizeof(double));
		for(int i2 = 0 ; i2 < cell_num*2+1 ; i2++){
			u2d[cell_num*2 - i2] = u2[0] - 0.5*((double) i2)*du2dxi2[0];
//			u2d[i2] = u2[0] - 0.5*((double) i2)*du2dxi2[0]; comment-out on 20180711
		}
		for(int i_step = 0 ; i_step < 2 ; i_step++){
			output_file(argv, Yeed, RecDep, str_file, Ld, 
				u0, u1, u2d, 
				N_str, media_no, 
				shape_str, shape_length, 
				origin_str, array_step, n_step, xi_file, k_0);
			// Prepare Top region.
			snprintf(Yeed,sizeof(Yeed),"%sT", Yee);//revised on 20180723 https://www.ipa.go.jp/security/awareness/vendor/programmingv1/b06_02.html
			for(int i2 = 0 ; i2 < cell_num*2+1 ; i2++){
				u2d[i2] = u2[L[2]*2] + 0.5*((double) i2)*du2dxi2[L[2]*2];
			}
		}
		SAFEFREE(u2d);
	}else{
		fprintf(stderr,"'Region' is not 'all' or 'outer'!\n");
		exit(1);
	}
	SAFEFREE(u0);	SAFEFREE(u1);	SAFEFREE(u2);	SAFEFREE(du2dxi2);
	{
		time_t timer_f = time(0);
		fprintf(stderr,"# Finish time of Yee's lattice creation = %s\n", ctime(&timer_f));
		clock_t endClock = clock();
		double cpusec = (endClock - startClock)/(double)CLOCKS_PER_SEC;
		fprintf(stderr,"# Exec time  = %1.0f sec, CPU time = %1.6f sec\n", difftime(timer_f,timer_s), cpusec);
	}
}
//======================================================================
//  input_file1 ---- rm_space, ScaleUnit, rm_comma                      
//                                       Last updated on Jul 19, 2018.  
//======================================================================
double ScaleUnit(char x[]);// unit variation: km, m, cm, mm, micron, um, nm, deg, degree, rad, radian.
void rm_space( char *A );
void rm_comma( char *A );
int input_file1(FILE *fp_i1, int L[3], char xi_file[], double *k_0, char str_file[], int *RecDep, char Yee[], int *cell_num, char Region[])
{
	char buf[BUFSIZE];	// buffer for fgets
	for(int j = 0 ; j < 3 ; j++) {L[j] = -1;}
	int j_count = 0;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 4) { 
		rm_space(buf);
		if(strncmp(buf, "xi2u", 4) == 0 && sscanf(buf,"%*[^=] %*[=] %s", xi_file) == 1){
			j_count++;// 1
		}else if(strncmp(buf, "outer", 5) == 0 && sscanf(buf,"%*[^=] %*[=] %d", cell_num) == 1){
			j_count++;// 2
		}else if(strncmp(buf, "StrMap", 6) == 0 && sscanf(buf,"%*[^=] %*[=] %s %*[^=] %*[=] %d %*[^=] %*[=] %s", 
			str_file, RecDep, Region) == 3){//InputFile = ./3d_str.txt, RecursiveDepth = 6,  Region = all / outer
			j_count++;// 3
		}else if(strncmp(buf, "Prefix", 6) == 0 && sscanf(buf,"%*[^=] %*[=] %s", Yee) == 1){//Prefix name = ./Data/joint
			j_count++;// 4
		}
	}
	rm_comma(xi_file);	rm_comma(str_file);	rm_comma(Yee);	rm_comma(Region);
	if(j_count != 4) {fprintf(stderr,"File-names can not read @ input_file1!\n");	exit(1);}
	else{j_count = 0;}
	FILE *fp_ix;
	fp_ix = fopen(xi_file,"r");
	if (fp_ix == NULL){
		fprintf(stderr,"%s can not open!\n", xi_file);
		exit(1);
	}else{
		while(fgets(buf, sizeof( buf ), fp_ix) != NULL){
			if(j_count < 3){
				if(sscanf(buf,"%*[^=] %*[=] %d", &L[j_count]) == 1){	j_count++;}//j_count == 3
			}else if(j_count == 3){
				char dummy[BUFSIZE];
				if(sscanf(buf,"%*[#] %s %*[^=] %*[=] %lf", dummy, k_0) == 2){//revised on 20180718
					if(strncmp(dummy,"k0",2) == 0){
						j_count++;
					}else{
						fprintf(stderr,"%s is not 'k0'!\n", dummy);
					}
				}//j_count == 3+1
			}
		}
	}
	fclose(fp_ix);
	for(int j = 0 ; j < 3 ; j++) {
		if(L[j] < 0){
			fprintf(stderr,"The L[%d] was not read @ input_file1!\n",j);
			j_count--;
		}
	}
	return(j_count - (3+1));
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
//  input_file2                                            
//                                        Last updated on Jul 05, 2018  
//======================================================================
void input_file2(FILE *fp_i2, const int L, double u0[L*2+1], double du0dxi0[L*2+1])
{
	char buf[BUFSIZE];	// buffer for fgets
	int i0 = 0;
	while(fgets(buf, sizeof( buf ), fp_i2) != NULL) { 
		double dummy;
		if(strncmp(buf, "#", 1) != 0 && sscanf(buf,"%lf  %*[,] %lf %*[,] %lf", 
			&dummy, &u0[i0], &du0dxi0[i0]) == 3){
			if(dummy == 0.5*((double) i0)) i0++;
			if(i0 == L*2+1) goto End_of_x_data;
		}
	}
	fprintf(stderr,"Full %d-data were not read @ input_fileM!\n",L);
	exit(1);
	End_of_x_data:;
}
//======================================================================
//  input_N_str                                                         
//                                        Last updated on Sep 08, 2014  
//======================================================================
int input_N_str(FILE *fp_i3)
{
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
	return(n_str0);
}
//======================================================================
// void input_file3 --+-- int input_media , rm_space, rm_comma          
//                    |                                                 
//                    +-- int input_shape ------+                       
//                    |                         +-- double ScaleUnit    
//                    +-- int input_origin -----+                       
//                    |                         |                       
//                    +-- int input_array ------+                       
//                                        Last updated on Jul 09, 2018  
//======================================================================
int input_media(const char buf[], const int n_str0, int media_no[]);
int input_shape(const char buf[], const double k_0, const int n_str0, int shape_str[], double shape_length[]);
int input_origin(const char buf[], const double k_0, const int n_str0, double origin_str[]);
int input_array(const char buf[], const double k_0, const int n_str0, const int k_array, double array_step[], int n_step[]);
void input_file3(FILE *fp_i3, const double k_0, const int N_str, 
				int media_no[N_str], int shape_str[N_str], double shape_length[N_str*6], 
				double origin_str[N_str*3], double array_step[N_str*9], int n_step[N_str*6])
{
// (int) N_str                    : total number of structures.
// (int) media_no[N_str]          : media number of a structure.
// (int) shape_str[N_str]         : 'block' == 0, 'cylinder' == 1, 
//                                  'prism' == 2,  'pyramid' == 3, 'spheroid' == 4.
// (double) shaple_length[N_str*6]: length parameter [meter],
//                                  shape_length[][0] = A, 
//                                  shape_length[][1] = B, 
//                                  shape_length[][2] = C.
//                                  structure rotation [radian], roll -> pitch -> yaw,
//                                  shape_length[][3] = theta_x ( roll angle: theta), 
//                                  shape_length[][4] = theta_y (pitch angle:  phi ), 
//                                  shape_length[][5] = theta_z ( yaw  angle:  psi ).
// (double) origin_str[N_str*3]   : an origin [meter] of (u0, u1, u2) coordinates for a structure.
// (double) array_step[N_str*9]   : step lengths [meter] for (1st stage, 2nd stage, 3rd stage),
//                                  x-direction (array_step[][0], array_step[][3], array_step[][6]), 
//                                  y-direction (array_step[][1], array_step[][4], array_step[][7]), 
//                                  z-direction (array_step[][2], array_step[][5], array_step[][8]). 
// (int) n_step[N_str*6]          : step number for (1st stage, 2nd stage, 3rd stage),
//                                  initial number (n_step[][0], n_step[][2], n_step[][4]), 
//                                   final  number (n_step[][1], n_step[][3], n_step[][5]). 
	char buf[BUFSIZE], command[16];	// buffer for fgets
	int n_str0 = 0;
	for(int J0 = 0 ; J0 < N_str ; J0++) { 
		while(fgets(buf, sizeof( buf ), fp_i3) != NULL) {
			rm_space(buf);
			if(strncmp(buf, "#", 1) != 0 && sscanf(buf,"%s", command) != EOF){
				rm_comma(command);
				if (strcmp(command, "begin") == 0) {
					int media_count = 0;// Initialization
					int shape_count = 0;
					int origin_count = 0;
					int array_count = 0;
					for(int k_array = 0 ; k_array < 3 ; k_array++) { 
						array_step[n_str0*9+  k_array*3] = 0.;
						array_step[n_str0*9+1+k_array*3] = 0.;
						array_step[n_str0*9+2+k_array*3] = 0.;
						n_step[n_str0*6+  k_array*2] = 0;
						n_step[n_str0*6+1+k_array*2] = 0;
					}// End of initialization
					while(fgets(buf, sizeof( buf ), fp_i3) != NULL) { 
						rm_space(buf);
						if(strncmp(buf, "#", 1) != 0 && sscanf(buf,"%s", command) != EOF){
							rm_comma(command);
							if (strcmp(command, "media") == 0) {
								media_count += input_media(buf, n_str0, media_no);
							}else if (strcmp(command, "shape") == 0) {
								shape_count += input_shape(buf, k_0, n_str0, shape_str, shape_length);
							}else if (strcmp(command, "origin") == 0) {
								origin_count += input_origin(buf, k_0, n_str0, origin_str);
							}else if (strcmp(command, "array") == 0) {
								array_count += input_array(buf, k_0, n_str0, array_count, array_step, n_step);
							}else if (strcmp(command, "end") == 0) {
								if(media_count != 1 || shape_count != 1 || origin_count != 1 || array_count > 3){
									fprintf(stderr,"end error @ input_file3!: media_count = %d, shape_count = %d,", 
										media_count, shape_count);
									fprintf(stderr," origin_count = %d, array_count= %d\n", origin_count, array_count);
									exit(1);
								}
								else{
									n_str0 += 1;
									if(n_str0 > N_str){
										fprintf(stderr,"Error @ input_file3!: n_str0 > N_str, n_str0 = %d,", n_str0);
										exit(1);
									}
									goto NEXT_step;
								}
							}
						}
					}
				}
			}
		}
		NEXT_step:;
	}
	if(n_str0 != N_str){
		fprintf(stderr,"Error @ input_file3!: n_str0 != N_str, n_str0 = %d,", n_str0);
		exit(1);
	}
}
int input_media(const char buf[], const int n_str0, int media_no[])
{
// (int) media_no[N_str]     	  : media number.
	if(sscanf(buf,"%*[^=] %*[=] %d", &media_no[n_str0]) != 1){
		fprintf(stderr,"Error in input_media that cannot read media_no @ n_str0 = %d\n", n_str0);
		return(0);
	}else{
		return(1);
	}
}
int input_shape(const char buf[], const double k_0, const int n_str0, int shape_str[], double shape_length[])
{
// (int) shape_str[N_str]         : 'block' == 0, 'cylinder' == 1, 
//                                  'prism' == 2,  'pyramid' == 3, 'spheroid' == 4.
// (double) shaple_length[N_str*6]: shape_length[][0] = A, 
//                                  shape_length[][1] = B, 
//                                  shape_length[][2] = C.
//                                  structure rotation [radian], roll -> pitch -> yaw,
//                                  shape_length[][3] = theta_x ( roll angle: theta), 
//                                  shape_length[][4] = theta_y (pitch angle:  phi ), 
//                                  shape_length[][5] = theta_z ( yaw  angle:  psi ).
	int count = 0;
	char command[BUFSIZE], shape[10], A[BUFSIZE], B[BUFSIZE], C[BUFSIZE], x[BUFSIZE], y[BUFSIZE], z[BUFSIZE];
	if(sscanf(buf,"%s %s %*[^=] %*[=] %lf %s %*[^=] %*[=] %lf %s %*[^=] %*[=] %lf %s %*[^=] %*[=] %lf %s %*[^=] %*[=] %lf %s %*[^=] %*[=] %lf %s", 
		command, shape, 
		&shape_length[n_str0*6+0],A, &shape_length[n_str0*6+1],B, &shape_length[n_str0*6+2],C, 
		&shape_length[n_str0*6+3],x, &shape_length[n_str0*6+4],y, &shape_length[n_str0*6+5],z) != 14){
		fprintf(stderr,"command, shape error! `%s'\n", buf);
		exit(1);
	}
	shape_str[n_str0] = -1;
	if (strcmp(shape, "block") == 0)  shape_str[n_str0] = 0;
	if (strcmp(shape, "cylinder") == 0)  shape_str[n_str0] = 1;
	if (strcmp(shape, "prism") == 0)  shape_str[n_str0] = 2;
	if (strcmp(shape, "pyramid") == 0)  shape_str[n_str0] = 3;
	if (strcmp(shape, "spheroid") == 0)  shape_str[n_str0] = 4;
	shape_length[n_str0*6+0] *= ScaleUnit(A)*k_0;
	shape_length[n_str0*6+1] *= ScaleUnit(B)*k_0;
	shape_length[n_str0*6+2] *= ScaleUnit(C)*k_0;
	shape_length[n_str0*6+3] *= ScaleUnit(x)*k_0;
	shape_length[n_str0*6+4] *= ScaleUnit(y)*k_0;
	shape_length[n_str0*6+5] *= ScaleUnit(z)*k_0;
	if(shape_length[n_str0*6+0] < 0. || shape_length[n_str0*6+1] < 0. || shape_length[n_str0*6+2] < 0.) {
		fprintf(stderr,"command, shape error!  Parameter is negative. `%s'\n", buf);
		exit(1);
	}
	count = 1;
	return(count);
}
int input_origin(const char buf[], const double k_0, const int n_str0, double origin_str[])
{
// origin_str[N_str*3]   :  an origin of (u0, u1, u2) coordinates for a structure.(meter, double).
	int count = 0;
	char command[BUFSIZE], x[BUFSIZE], y[BUFSIZE], z[BUFSIZE];
	if(sscanf(buf,"%s %*[^=] %*[=] %lf %s %*[^=] %*[=] %lf %s %*[^=] %*[=] %lf %s", command, &origin_str[n_str0*3+0],x, &origin_str[n_str0*3+1],y, &origin_str[n_str0*3+2],z) != 7){
		fprintf(stderr,"command, shape error! `%s'\n", buf);
		exit(1);
	}
	origin_str[n_str0*3+0] *= ScaleUnit(x)*k_0;
	origin_str[n_str0*3+1] *= ScaleUnit(y)*k_0;
	origin_str[n_str0*3+2] *= ScaleUnit(z)*k_0;
	count = 1;
	return(count);
}
int input_array(const char buf[], const double k_0, const int n_str0, const int k_array, double array_step[], int n_step[])
{
// (double) array_step[N_str*9]   : step lengths [meter] for (1st stage, 2nd stage, 3rd stage),
//                                  x-direction (array_step[][0], array_step[][3], array_step[][6]), 
//                                  y-direction (array_step[][1], array_step[][4], array_step[][7]), 
//                                  z-direction (array_step[][2], array_step[][5], array_step[][8]). 
// (int) n_step[N_str*6]          : step number for (1st stage, 2nd stage, 3rd stage),
//                                  initial number (n_step[][0], n_step[][2], n_step[][4]), 
//                                   final  number (n_step[][1], n_step[][3], n_step[][5]). 
	int count = 0;
	char command[BUFSIZE], x[BUFSIZE], y[BUFSIZE], z[BUFSIZE];
	if(sscanf(buf,"%s %*[^=] %*[=] %lf %s %*[^=] %*[=] %lf %s %*[^=] %*[=] %lf %s  %*[^=] %*[=] %d  %*[^=] %*[=] %d", 
		command, &array_step[n_str0*9+k_array*3], x, &array_step[n_str0*9+1+k_array*3], y, &array_step[n_str0*9+2+k_array*3], 
		z, &n_step[n_str0*6+k_array*2], &n_step[n_str0*6+1+k_array*2]) != 9){
		fprintf(stderr,"command, array error! `%s'\n", buf);
		exit(1);
	}
	array_step[n_str0*9+  k_array*3] *= ScaleUnit(x)*k_0;
	array_step[n_str0*9+1+k_array*3] *= ScaleUnit(y)*k_0;
	array_step[n_str0*9+2+k_array*3] *= ScaleUnit(z)*k_0;
	count = 1;
	return(count);
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
//     void output_file --+----------------+                            
//                        |                +-- media_f                  
//                        +-- est_volume --+                            
//                          (recursive use)                             
//                                        Last updated on Jul 11, 2018  
//======================================================================
int media_f(const double X[3], const int N_str, 
			const int media_no[N_str], const int shape_str[N_str], 
			const double shape_length[N_str*6], const double origin_str[N_str*3], 
			const double array_step[N_str*9], const int n_step[N_str*6]);
void est_volume(int media_p[8*8], double media_r[8*8], int D_depth[2],
			const double XXX[9], const int N_str, 
			const int media_no[N_str], const int shape_str[N_str],
			const double shape_length[N_str*6], const double origin_str[N_str*3], 
			const double array_step[N_str*9], const int n_step[N_str*6]);
void output_file(char **argv, const char Yee[], 
			const int RecDep, const char str_file[], const int L[3], 
			const double u0[L[0]*2+1], const double u1[L[1]*2+1], const double u2[L[2]*2+1], 
			const int N_str, const int media_no[N_str], const int shape_str[N_str], 
			const double shape_length[N_str*6], const double origin_str[N_str*3], 
			const double array_step[N_str*9], const int n_step[N_str*6], const char xi_file[], const double k_0)
{
// This 'output_file' creates a set of 3D files for optical media.
// (char) Yee[]                   : an output file name.
// (int) RecDep,           : RecDep, means depth of recursion for est_volume.
// (int) L[3]                     : space size of periodic (u, v, w) coordinates.
// (double) u0[L[0]*2+1],
//          u1[L[1]*2+1],
//          u2[L[2]*2+1]          : sub-lattice position [meter] in (u0, u1, u2) coordinates.
// (int) N_str                    : structure quantity.
// (int) media_no[N_str]     	  : media number.
// (int) shape_str[N_str]         : 'block' == 0, 'cylinder' == 1, 
//                                  'prism' == 2,  'pyramid' == 3, 'spheroid' == 4.
// (double) shaple_length[N_str*6]: shape_length[][0] = A, 
//                                  shape_length[][1] = B, 
//                                  shape_length[][2] = C.
//                                  structure rotation [radian], roll -> pitch -> yaw,
//                                  shape_length[][3] = theta_x ( roll angle: theta), 
//                                  shape_length[][4] = theta_y (pitch angle:  phi ), 
//                                  shape_length[][5] = theta_z ( yaw  angle:  psi ).
// (double) origin_str[N_str*3]   : an origin [meter] of (u0, u1, u2) coordinates for a structure.
// (double) array_step[N_str*9]   : step lengths [meter] for (1st stage, 2nd stage, 3rd stage),
//                                  x-direction (array_step[][0], array_step[][3], array_step[][6]), 
//                                  y-direction (array_step[][1], array_step[][4], array_step[][7]), 
//                                  z-direction (array_step[][2], array_step[][5], array_step[][8]). 
// (int) n_step[N_str*6]          : step number for (1st stage, 2nd stage, 3rd stage),
//                                  initial number (n_step[][0], n_step[][2], n_step[][4]), 
//                                   final  number (n_step[][1], n_step[][3], n_step[][5]). 
// (double) k_0                   : normalization constant [meter^-1].
	int Max_media_no = media_no[0];
	for(int j0 = 1 ; j0 < N_str ; j0++) {
		if(Max_media_no < media_no[j0]) Max_media_no = media_no[j0];
	}
	FILE *fp_o;
	int L_max = L[0]+1;
	int M_max = L[1]+1;
	int *sheet_uvw = calloc(L_max*M_max*2,sizeof(int));//revised on 20180725
	char buf[BUFSIZE];
	for(int jw0 = 0 ; jw0 < L[2] ; jw0++) { 
		snprintf(buf,sizeof(buf),"%s%d.dat", Yee, jw0);//revised on 20180723 https://www.ipa.go.jp/security/awareness/vendor/programmingv1/b06_02.html
		fp_o = fopen(buf,"w");
		if (fp_o == NULL){
			fprintf(stderr,"open error for %s\n", buf);
			exit(1);
		}
		time_t timer_s = time(0);
		clock_t startClock = clock();
		fprintf(stderr,"# Created on %s", ctime(&timer_s));
		fprintf(fp_o,"# Created on %s", ctime(&timer_s));
		fprintf(fp_o,"# exec = %s : `%s' was compiled at %s on %s by C-version:%ld\n", argv[0], __FILE__, __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(fp_o,"# info = %s, u-xi = %s, str. = %s, this = %s\n", argv[1], xi_file, str_file, buf);//revised on 20180718
		fprintf(fp_o,"# k_0 * u_2  = %.20e ; %.20e ; %.20e, k_0 = %.20e m^-1\n", u2[jw0*2+0], u2[jw0*2+1], u2[jw0*2+2], k_0);// revised on 20180801
		fprintf(fp_o,"# j_u0 < %d, j_u1 < %d, j_u2 < %d, media_no < %d\n", L[0], L[1], L[2], Max_media_no+1);
		fprintf(fp_o,"# bi_u0 bi_u1 bi_u2 media_max\n");// u0, u1 and u2 mean u, v and w, respectively.
		fprintf(fp_o,"# media_no. ratio\n");
		for(int jv0 = 0 ; jv0 <= L[1] ; jv0++) { 
			for(int ju0 = 0 ; ju0 <= L[0] ; ju0++) {
				if(jw0 == 0){
					double X[3] = {u0[ju0*2+0], u1[jv0*2+0], u2[jw0*2+0]};
					sheet_uvw[ju0+L_max*(jv0+M_max*0)] = 
						media_f(X, N_str, media_no, shape_str, shape_length, origin_str, array_step, n_step);
				}else{
					sheet_uvw[ju0+L_max*(jv0+M_max*0)] = sheet_uvw[ju0+L_max*(jv0+M_max*1)];
				}
				double X[3] = {u0[ju0*2+0], u1[jv0*2+0], u2[jw0*2+2]};
				sheet_uvw[ju0+L_max*(jv0+M_max*1)] = 
					media_f(X, N_str, media_no, shape_str, shape_length, origin_str, array_step, n_step);
			}
		}
		int media_p[8*8];
		double media_r[8*8];
		for(int jv0 = 0 ; jv0 < L[1] ; jv0++) { 
			for(int ju0 = 0 ; ju0 < L[0] ; ju0++) {
				for(int kw0 = 0 ; kw0 < 2 ; kw0++) { 
					for(int kv0 = 0 ; kv0 < 2 ; kv0++) { 
						for(int ku0 = 0 ; ku0 < 2 ; ku0++) { 
							media_p[0+8*(ku0+2*(kv0+2*kw0))] = sheet_uvw[(ju0+ku0) + L_max*((jv0+kv0) + M_max*kw0)];
							media_r[0+8*(ku0+2*(kv0+2*kw0))] = 1.e0;
							for(int i0 = 1 ; i0 < 8 ; i0++) {
								media_p[i0+8*(ku0+2*(kv0+2*kw0))] = -1;//// This negative value means void of media when 1 <= i0 <= 7.
								media_r[i0+8*(ku0+2*(kv0+2*kw0))] = 0.;
							}
						}
					}
				}
				int check_i = -1;
				for(int i0 = 1 ; i0 < 8 ; i0++){
					if(media_p[0+8*i0] != media_p[0+8*0]) check_i = 1;
				}
				fprintf(fp_o,"%d %d\n", ju0, jv0);
				if(check_i < 0){
					for(int kw0 = 0 ; kw0 < 2 ; kw0++) { 
						for(int kv0 = 0 ; kv0 < 2 ; kv0++) { 
							for(int ku0 = 0 ; ku0 < 2 ; ku0++) { 
								fprintf(fp_o,"%d %d %d %d\n", ku0,kv0,kw0, 1);
								fprintf(fp_o,"%d %.20E\n", media_p[0+8*(ku0+2*(kv0+2*kw0))], media_r[0+8*(ku0+2*(kv0+2*kw0))]);
							}
						}
					}
				}else{
					double XXX[9] = {u0[ju0*2+0], u1[jv0*2+0], u2[jw0*2+0],// Yee cell origin
										u0[ju0*2+1], u1[jv0*2+1], u2[jw0*2+1],// Yee cell center
										u0[ju0*2+2], u1[jv0*2+2], u2[jw0*2+2]};// Next cell origin
					int D_depth[2] = {0,RecDep};					
					est_volume(media_p, media_r, D_depth, XXX, N_str, 
									media_no, shape_str, shape_length, origin_str, array_step, n_step);
					for(int kw0 = 0 ; kw0 < 2 ; kw0++) { 
						for(int kv0 = 0 ; kv0 < 2 ; kv0++) { 
							for(int ku0 = 0 ; ku0 < 2 ; ku0++) {
								int i0 = 0;
								while(media_p[i0+8*(ku0+2*(kv0+2*kw0))] >= 0 && i0 < 9) i0++;// If media_p is not void, 8 kinds of media is maximum case.
								fprintf(fp_o,"%d %d %d %d\n", ku0,kv0,kw0, i0);
								for(int j0 = 0 ; j0 < i0 ; j0++) fprintf(fp_o,"%d %.20E\n", 
									media_p[j0+8*(ku0+2*(kv0+2*kw0))], media_r[j0+8*(ku0+2*(kv0+2*kw0))]);
								double sum = 0.;// Check sum
								for(int j0 = 0 ; j0 < i0 ; j0++) sum += media_r[j0+8*(ku0+2*(kv0+2*kw0))];
								if(sum != 1.) fprintf(fp_o,"volume error: sum = %.20E\n", sum);
							}
						}
					}
				}
			}
		}
		time_t timer_f = time(0);
		double t_span = difftime(timer_f,timer_s);
		clock_t endClock = clock();
		double cpusec = (endClock - startClock)/(double)CLOCKS_PER_SEC;
		fprintf(fp_o,"# Exec time  = %1.0f sec, CPU time = %1.6f sec\n", t_span, cpusec);
		if(fclose(fp_o) != 0) {
			fprintf(stderr,"fclose error for %s\n", buf);
			exit(1);
		}
	}
	SAFEFREE(sheet_uvw);
}
//======================================================================
//  void est_volume: This is recursively used.                          
//  Dey-Mittra averaging is coded to this function.                     
//  See: S. Dey & R. Mittra, "A conformal finite-difference time-domain 
//  technique for modeling cylindrical dielectric resonators,"          
//  IEEE Transactions on Microwave Theory and Techniques,               
//  vol. 47, pp.1737-1739, 1999.                                        
//                                        Last updated on Aug 15, 2017  
//======================================================================
void est_volume(int media_p[8*8], double media_r[8*8], int D_depth[2],
		const double XXX[9], const int N_str, 
		const int media_no[N_str], const int shape_str[N_str], 
		const double shape_length[N_str*6], const double origin_str[N_str*3], 
		const double array_step[N_str*9], const int n_step[N_str*6])
{
// (int) N_str                    : total number of structures.
// (int) media_no[N_str]          : media number of a structure.
// (int) shape_str[N_str]         : 'block' == 0, 'cylinder' == 1, 
//                                  'prism' == 2,  'pyramid' == 3, 'spheroid' == 4.
// (double) shaple_length[N_str*6]: length parameter [meter],
//                                  shape_length[][0] = A, 
//                                  shape_length[][1] = B, 
//                                  shape_length[][2] = C.
//                                  structure rotation [radian], roll -> pitch -> yaw,
//                                  shape_length[][3] = theta_x ( roll angle: theta), 
//                                  shape_length[][4] = theta_y (pitch angle:  phi ), 
//                                  shape_length[][5] = theta_z ( yaw  angle:  psi ).
// (double) origin_str[N_str*3]   : an origin [meter] of (u0, u1, u2) coordinates for a structure.
// (double) array_step[N_str*9]   : step lengths [meter] for (1st stage, 2nd stage, 3rd stage),
//                                  x-direction (array_step[][0], array_step[][3], array_step[][6]), 
//                                  y-direction (array_step[][1], array_step[][4], array_step[][7]), 
//                                  z-direction (array_step[][2], array_step[][5], array_step[][8]). 
// (int) n_step[N_str*6]          : step number for (1st stage, 2nd stage, 3rd stage),
//                                  initial number (n_step[][0], n_step[][2], n_step[][4]), 
//                                   final  number (n_step[][1], n_step[][3], n_step[][5]). 
	D_depth[0]++;//Then D_depth[0] starts from 1.
	if( D_depth[0] <= D_depth[1] ){//20141028 revised, D_depth[1] = RecDep,
		int media_p3[3*3*3]; // media numbers at 27 points in a devided Yee cell for recursive operation
		for(int kw0 = 0 ; kw0 < 2 ; kw0++) {
			for(int kv0 = 0 ; kv0 < 2 ; kv0++) { 
				for(int ku0 = 0 ; ku0 < 2 ; ku0++) {
					media_p3[(ku0*2)+3*((kv0*2)+3*(kw0*2))] = media_p[0+8*(ku0+2*(kv0+2*kw0))];
				}// copy media numbers at 8 points in 27 points of a devided Yee cell for recursive operation
			}
		}
		for(int kw0 = 0 ; kw0 < 3 ; kw0++) { 
			for(int kv0 = 0 ; kv0 < 3 ; kv0++) { 
				for(int ku0 = 0 ; ku0 < 3 ; ku0++) {
					if(ku0 == 1 || kv0 == 1 || kw0 == 1){//20170814 revised: if(ku0*(ku0-2) != 0 || kv0*(kv0-2) != 0 || kw0*(kw0-2) != 0){
						double X[3];// one of 27 points in a devided Yee cell for recursive operation
						X[0] = XXX[0+3*ku0];// x(ku0,kv0,kw0)
						X[1] = XXX[1+3*kv0];// y(ku0,kv0,kw0)
						X[2] = XXX[2+3*kw0];// z(ku0,kv0,kw0)
						media_p3[ku0+3*(kv0+3*kw0)] = // media numbers at 19 points in a devided Yee cell for recursive operation
							media_f(X, N_str, media_no, shape_str, shape_length, origin_str, array_step, n_step);
					}// select 19 points in 27 points of a devided Yee cell for recursive operation
				}
			}
		}
		for(int kw0 = 0 ; kw0 < 2 ; kw0++) { 
			for(int kv0 = 0 ; kv0 < 2 ; kv0++) { 
				for(int ku0 = 0 ; ku0 < 2 ; ku0++) {// 8 points in 27 points of a devided Yee cell for recursive operation
					int media_p2[2*2*2];
					for(int i0 = 0 ; i0 < 8 ; i0++) media_p2[i0] = -1;// media_p2 is initialized by negative value.
					for(int LW = 0 ; LW < 2 ; LW++) { 
						for(int LV = 0 ; LV < 2 ; LV++) { 
							for(int LU = 0 ; LU < 2 ; LU++) {// copy media_p3 to media_p2 at 8 points for each of 8 points
								media_p2[LU+2*(LV+2*LW)] = media_p3[(ku0+LU)+3*((kv0+LV)+3*(kw0+LW))];
							}
						}
					}
					int check_i = -1;
					for(int i0 = 1 ; i0 < 8 ; i0++){
						if(media_p2[i0] != media_p2[0]) check_i = 1;
					}
					if(check_i < 0){// all media_p2[i0] have same value in sub cell (ku0,kv0,kw0).
						media_p[0+8*(ku0+2*(kv0+2*kw0))] = media_p2[0];
						media_r[0+8*(ku0+2*(kv0+2*kw0))] = 1.;
						for(int i0 = 1 ; i0 < 8 ; i0++){
							media_p[i0+8*(ku0+2*(kv0+2*kw0))] = -1;// This negative value means void of media
							media_r[i0+8*(ku0+2*(kv0+2*kw0))] = 0.;
						}
					}else{// If all media_p2[i0] are not same value, this scope prepares recursive operation.
						int media_pD[8*8];//media_p for sub region
						double media_rD[8*8];//media_r for sub region
						for(int i0 = 0 ; i0 < 8 ; i0++) {
							media_pD[0+8*i0] = media_p2[i0];
							media_rD[0+8*i0] = 1.e0;
							for(int j0 = 1 ; j0 < 8 ; j0++) {
								media_pD[j0+8*i0] = -1;
								media_rD[j0+8*i0] = 0.;
							}
						}
						{//recursive region
							double XXXD[9] = {XXX[0+3*ku0],     XXX[1+3*kv0],     XXX[2+3*kw0], 
												0,                0,                0, // This line will be defined by average. 
												XXX[0+3*(ku0+1)], XXX[1+3*(kv0+1)], XXX[2+3*(kw0+1)]};//20141028 revised
							for(int i0 = 0 ; i0 < 3 ; i0++) XXXD[i0+3*1] = ldexp(XXXD[i0+3*0]+XXXD[i0+3*2],-1);// XXXD[i0+3*1] = 0.5*(XXXD[i0+3*0]+XXXD[i0+3*2]);
							int D_depthD[2] = {D_depth[0],D_depth[1]};					
							est_volume(media_pD, media_rD, D_depthD, XXXD, N_str, media_no, 
										shape_str, shape_length, origin_str, array_step, n_step);
						}//end of recursive region
						for(int i0 = 1 ; i0 < 8 ; i0++) {//shift i0-th column elements to 0-th column for media_pD and media_rD
							for(int k0 = 0 ; k0 < 8 ; k0++) {
								if(media_pD[k0+8*i0] >= 0){
									for(int j0 = 0 ; j0 < 8 ; j0++) {
										if(media_pD[j0+8*0] == media_pD[k0+8*i0]){
											media_rD[j0+8*0] += media_rD[k0+8*i0];
											goto next_search;
										}
										if(media_pD[j0+8*0] < 0){
											media_pD[j0+8*0] = media_pD[k0+8*i0];
											media_rD[j0+8*0] = media_rD[k0+8*i0];
											goto next_search;
										}
									}
								}
								next_search:;
								media_rD[k0+8*i0] = 0.;
								media_pD[k0+8*i0] = -1;
							}
						}//end of shift i0-th column elements to 0-th column for media_pD and media_rD
						for(int j0 = 0 ; j0 < 8 ; j0++) {
							media_p[j0+8*(ku0+2*(kv0+2*kw0))] = media_pD[j0+8*0];
							media_r[j0+8*(ku0+2*(kv0+2*kw0))] = ldexp(media_rD[j0+8*0],-3);//normalizing volume: 0.125*media_rD[j0+8*0];
						}//Copy sub-cell data to cell data
					}
				}
			}
		}
	}
}
//======================================================================
//       int media_f ---- void rotation, int shap_region                
//                                        Last updated on Aug 15, 2017  
//======================================================================
void rotation(const double theta, 
				const double x0, const double y0, 
				double *x1, double *y1);
int shape_region(const int type, 
				const double A, const double B, const double C, 
				const double u0, const double u1, const double u2);
//======================================================================
int media_f(const double X[3], const int N_str, 
		const int media_no[N_str], const int shape_str[N_str], 
		const double shape_length[N_str*6], const double origin_str[N_str*3], 
		const double array_step[N_str*9], const int n_step[N_str*6])
{
// This 'media_f' gives a media number at a point.
// (double) X[3]                  : (X[0],X[1],X[2])[meter] is a point in 3D structure.
// (int) N_str                    : total number of structures.
// (int) media_no[N_str]          : media number of a structure.
// (int) shape_str[N_str]         : 'block' == 0, 'cylinder' == 1, 
//                                  'prism' == 2,  'pyramid' == 3, 'spheroid' == 4.
// (double) shaple_length[N_str*6]: length parameter [meter],
//                                  shape_length[][0] = A, 
//                                  shape_length[][1] = B, 
//                                  shape_length[][2] = C.
//                                  structure rotation [radian], roll -> pitch -> yaw,
//                                  shape_length[][3] = theta_x ( roll angle: theta), 
//                                  shape_length[][4] = theta_y (pitch angle:  phi ), 
//                                  shape_length[][5] = theta_z ( yaw  angle:  psi ).
// (double) origin_str[N_str*3]   : an origin [meter] of (u0, u1, u2) coordinates for a structure.
// (double) array_step[N_str*9]   : step lengths [meter] for (1st stage, 2nd stage, 3rd stage),
//                                  x-direction (array_step[][0], array_step[][3], array_step[][6]), 
//                                  y-direction (array_step[][1], array_step[][4], array_step[][7]), 
//                                  z-direction (array_step[][2], array_step[][5], array_step[][8]). 
// (int) n_step[N_str*6]          : step number for (1st stage, 2nd stage, 3rd stage),
//                                  initial number (n_step[][0], n_step[][2], n_step[][4]), 
//                                   final  number (n_step[][1], n_step[][3], n_step[][5]). 
	int med = 0;// media number
	for(int J = N_str-1 ; J >= 0 ; J--){//revised on 20170815: for(int J = 0 ; J < N_str ; J++){
		for(int L2 = n_step[J*6+4] ; L2 <= n_step[J*6+5]; L2++){
			for(int L1 = n_step[J*6+2] ; L1 <= n_step[J*6+3]; L1++){
				for(int L0 = n_step[J*6+0] ; L0 <= n_step[J*6+1]; L0++){
					double xx = origin_str[J*3+0] + L0*array_step[J*9+0+0*3] + L1*array_step[J*9+0+1*3] + L2*array_step[J*9+0+2*3];
					double yy = origin_str[J*3+1] + L0*array_step[J*9+1+0*3] + L1*array_step[J*9+1+1*3] + L2*array_step[J*9+1+2*3];
					double zz = origin_str[J*3+2] + L0*array_step[J*9+2+0*3] + L1*array_step[J*9+2+1*3] + L2*array_step[J*9+2+2*3];
					// begin Cardan angles rotation (object rotation) X(roll) -> Y(pitch) -> Z(yaw) (X[0]-xx,X[1]-yy,X[2]-zz) --> (x1,y1,z1)
					// following calc. is inverse rotation (world rotation) Z(yaw) -> Y(pitch) -> X(roll)
					double x0, y0;// begin yaw (z-axis) angle rotation
					rotation(-shape_length[J*6+5], X[0] - xx, X[1] - yy, &x0, &y0);// end yaw angle rotation (-
					double z0, x1;// begin pitch (y-axis) angle rotation
					rotation(-shape_length[J*6+4], X[2] - zz, x0, &z0, &x1);// end pitch angle rotation
					double y1, z1;// begin roll (x-axis) angle rotation
					rotation(-shape_length[J*6+3], y0, z0, &y1, &z1);// end roll angle rotation
					//  end  Cardan angles rotation ZYX (X[0]-xx,X[1]-yy,X[2]-zz) --> (x1,y1,z1)
					if(shape_region(shape_str[J], shape_length[J*6+0],shape_length[J*6+1],shape_length[J*6+2], x1,y1,z1) == 1) {
						med = media_no[J];
						goto search_end;//added on 20170815
					}
				}
			}
		}
	}
	search_end:;//added on 20170815
	return(med);
}
//======================================================================
//  void rotation                                                       
//                                        Last updated on Sep 06, 2014  
//======================================================================
void rotation(const double theta, 
				const double x0, const double y0, 
				double *x1, double *y1)
{
// This `rotation' is the rotating operator of coordinates.
// theta is rotating angle (radian: double).
// (x0, y0) is an input coordinate.
// (x1, y1) is a rotated coordinate.
//
	double s = sin(theta);
	double c = cos(theta);
	*x1 = c*x0 - s*y0;
	*y1 = s*x0 + c*y0;
}
//======================================================================
//  int shape_region                                                    
//                                        Last updated on Aug 15, 2017  
//======================================================================
int shape_region(const int type, 
				const double A, const double B, const double C, 
				const double x, const double y, const double z)
{
// This 'shape_region' gives 1 (0) when (x,y,z) is in (out) a shape.
// type: 'block' == 0, 'cylinder' == 1, 'prism' == 2, 'pyramid' == 3, 'spheroid' == 4 (int).
//
	int tag = 0;
	if (type == 0) {//block
		if(A*A >= x*x && B*B >= y*y && C*C >= z*z) {tag = 1;}
	} else if (type == 1) {//cylinder
		if(x*x/(A*A) + y*y/(B*B) <= 1. && C*C >= z*z) {tag = 1;}
	} else if (type == 2) {//prism
		if(x >= 0. && y >= 0. && x/A + y/B <= 1. && C*C >= z*z) {tag = 1;}
	} else if (type == 3) {//pyramid
		if(x >= 0. && y >= 0. && z >= 0. && x/A + y/B + z/C <= 1.) {tag = 1;}
	} else if (type == 4) {//spheroid
		if(x*x/(A*A) + y*y/(B*B) + z*z/(C*C) <= 1.) {tag = 1;}
	} else {
		fprintf(stderr,"shape_region error!\n");
		tag = -1;//20170815 revised
		exit(1);
	}
	return(tag);
}
//////////////////////  block, cylinder, spheroid  /////////////////////
//                                                                      
//                                                                      
//             (-A,-B,+C) --                    -- (-A,+B,+C)           
//                      /|                       /|                     
//                                                                      
//                                                                      
//                                                                      
//                                                                      
//                                                                      
//                 /                        /                           
//      (+A,-B,+C) --                    -- (+A,+B,+C)                  
//                |                        |                            
//                                                                      
//                               |/                                     
//                             -----(0,0,0)                             
//                              /|                                      
//                                                                      
//                       |                        |                     
//             (-A,-B,-C) --                    -- (-A,+B,-C)           
//                      /                        /                      
//                                                                      
//                                                                      
//                                                                      
//                                                                      
//                                                                      
//                |/                       |/                           
//      (+A,-B,-C) --                    -- (+A,+B,-C)                  
//                                                                      
//                                                                      
////////////////////////////////////////////////////////////////////////
//
///////////////////////////////   prism   //////////////////////////////
//                                                                      
//                                                                      
//                       (0,0,+C) --           -- (0,+B,+C)             
//                              /|           *   |                      
//                                       *                              
//                                   *                                  
//                            /  *                                      
//                 (+A,0,+C) *                                          
//                           |                                          
//                               |/                                     
//                             -----(0,0,0)                             
//                              /|                                      
//                                                                      
//                                                                      
//                                                                      
//                                                                      
//                                                                      
//                               |               |                      
//                       (0,0,-C) --           -- (0,+B,-C)             
//                              /            *                          
//                                       *                              
//                                   *                                  
//                           |/  *                                      
//                 (+A,0,-C) *                                          
//                                                                      
//                                                                      
////////////////////////////////////////////////////////////////////////
//
//////////////////////////////   pyramid   /////////////////////////////
//                                                                      
//                                                                      
//                    (0,0,+C)                                          
//                            *                                         
//                            |   *                                     
//                           *|        *                                
//                            |             *                           
//                          *  ------------------* (0,+B,0)             
//                           /(0,0,0)       *                           
//                         */         *                                 
//                         /    *                                       
//               (+A,0,0) *                                             
//                                                                      
//                                                                      
////////////////////////////////////////////////////////////////////////
