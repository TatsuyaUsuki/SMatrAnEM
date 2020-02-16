#include "header_macro.h"
#include "constant.h"
struct SysPara{
	long double k0;
	char name[BUFSIZE];
	long L[3];
	long double Lsys[3];
	long K[3];
	long double min_dudxi[3];
};
//----------------------------------------------------------------------
//     main -+-- input_file0, input_file1 -- rm_space, ScaleUnit        
//           +-- Calc_main -- Calc_F -- Generate_G -- LogMulCal_ini 
//----------------------------------------------------------------------
void input_file0(FILE *fp_i1, struct SysPara *Mesh);
void input_file1(FILE *fp_i1, struct SysPara *Mesh);
void Calc_main(char **argv, const struct SysPara *Mesh);
//======================================================================
//  Main routine of xi2u.                                               
//======================================================================
int main(int argc, char **argv)
{
	if(argc != 2) {
		fprintf(stderr,"error: number of files \n");
		exit(1);
	}else if(strncmp(argv[1], "-v", 2) == 0 ) {// This scope was added on Jul 02, 2018
		fprintf(stderr,"The '%s' creates non-uniform 3D-lattice.\n", argv[0]);
		fprintf(stderr,"Version 18.08.03 is compiled at %s on %s\n C-version   : %ld\n", __TIME__, __DATE__, __STDC_VERSION__);
		fprintf(stderr," Source code : '%s'\n Author      : Tatsuya Usuki\n URL         : http://www.smatran.org\n", __FILE__);
		fprintf(stderr," Reference   : 'Wave scattering in frequency domain' as 'Formulation.pdf' on May 20, 2018\n");
		exit(0);//normal end
	}
	
	FILE *fp_i1;
	fp_i1 = fopen(argv[1],"r");
	if (fp_i1 == NULL){
		fprintf(stderr,"open error!: open input-file!\n");
		exit(1);
	}

	struct SysPara Mesh;
	{
		Mesh.k0 = 0.L;
		input_file0(fp_i1, &Mesh); // inputting parameters from an inputfile 1
		rewind(fp_i1);
		input_file1(fp_i1, &Mesh); // inputting parameters from an inputfile 1
	}
	fclose(fp_i1);
	
	Calc_main(argv, &Mesh); // calc. of EM and outputting results to an outputfile
	
}
void rm_space( char *A );
long double ScaleUnit(char x[]);
//======================================================================
//  input_file0, input_file1                                                        
//                           Last updated on Jul 25, 2018.  
//======================================================================
void input_file0(FILE *fp_i1, struct SysPara *Mesh)
{
	char buf[BUFSIZE];	// buffer for fgets
	int j_count = -1;
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 0) { 
		rm_space(buf);
		long double R_dummy;
		char A[BUFSIZE];
		if(strncmp(buf, "inv_k0", 6)*strncmp(buf, "Inv_k0", 6)*strncmp(buf, "INV_k0", 6) == 0
		&& sscanf(buf,"%*[^=] %*[=] %Lf %s", &R_dummy, A) == 2){
			if(R_dummy > 0.L) {
				Mesh->k0 = 1.L/(R_dummy*ScaleUnit(A));
				j_count++;
			}
		}
	}
	if(j_count != 0){
		fprintf(stderr,"The inv_k0 line can not be read @ input_file0!\n");
		exit(1);
	}
}
void input_file1(FILE *fp_i1, struct SysPara *Mesh)
{
	for(long j = 0 ; j < 3 ; j++) {Mesh->L[j] = -1;}
	int j_count = -4;
	char buf[BUFSIZE];	// buffer for fgets
	while(fgets(buf, sizeof( buf ), fp_i1) != NULL && j_count < 0) { 
		rm_space(buf);
		long j, L_dummy, K_dummy;
		long double R_dummy;
		char A[BUFSIZE];
		if(strncmp(buf, "mesh", 4)*strncmp(buf, "Mesh", 4)*strncmp(buf, "MESH", 4) == 0
		 && sscanf(buf,"%*[^=] %*[=] %ld %*[^=] %*[=] %ld %*[^=] %*[=] %Lf %s %*[^=] %*[=] %ld",
		  &j, &L_dummy, &R_dummy, A, &K_dummy) == 5){
			if( j >= 0 && j < 3 && L_dummy >= 0 && R_dummy > 0) {
				Mesh->L[j] = L_dummy;
				Mesh->Lsys[j] = R_dummy*(Mesh->k0*ScaleUnit(A));
				if(K_dummy <= 0){
					Mesh->K[j] = 0;
					Mesh->min_dudxi[j] = 0.L;
					j_count++;
					fprintf(stderr,"The min du%ld/dxi%ld was not used, since K%ld <= 0 @ input_file1.\n",j,j,j);
				}else{
					long double dudxi;
					char B[BUFSIZE];
					if(sscanf(buf,"%*[^=] %*[=] %ld %*[^=] %*[=] %ld %*[^=] %*[=] %Lf %s %*[^=] %*[=] %ld %*[^=] %*[=] %Lf %s",
							&j, &L_dummy, &R_dummy, A, &K_dummy, &dudxi, B) == 7){
						if(dudxi > 0.L) {
							Mesh->K[j] = K_dummy;
							Mesh->min_dudxi[j] = dudxi*(Mesh->k0*ScaleUnit(B));
							j_count++;
						}
					}
				}
			}
		}else if(strncmp(buf, "xi2u", 4) == 0 && sscanf(buf,"%*[^=] %*[=] %s", Mesh->name) == 1){
			j_count++;
		}
	}
	for(long j = 0 ; j < 3 ; j++) {
		if(Mesh->L[j] < 0){
			fprintf(stderr,"The L[%ld] was not read @ input_file1!\n",j);
			j_count--;
		}
	}
	if(j_count != 0){
		fprintf(stderr,"The %d lines was only read @ input_file1!\n",j_count);
		exit(1);
	}
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
//  ScaleUnit                                                           
//                              Last updated on May 17, 2018.  
//======================================================================
long double ScaleUnit(char x[])
{
	long double scale;
	if (strncmp(x, "cm", 2) == 0) {
		scale = 1.e-2L;
	} else if (strncmp(x, "deg", 3) == 0) {
		scale = Pi/180.L;
	} else if (strncmp(x, "km", 2) == 0) {
		scale = 1.e3L;
	} else if (strncmp(x, "m", 1) == 0) {
		if (strncmp(x, "mm", 2) == 0) {
			scale = 1.e-3L;// order of m, mm, micron is important!
		} else if (strncmp(x, "micron", 6) == 0) {
			scale = 1.e-6L;// order of m, mm, micron is important!
		} else {
			scale = 1.e-0L;// order of m, mm, micron is important!
		}
	} else if (strncmp(x, "nm", 2) == 0) {
		scale = 1.e-9L;
	} else if (strncmp(x, "um", 2) == 0) {
		scale = 1.e-6L;
	} else if (strncmp(x, "rad", 3) == 0) {
		scale = 1.L;
	} else {
		fprintf(stderr,"ScaleUnit error!  x = %s\n", x);
		exit(1);
	}
	return(scale);
}
//======================================================================
//  This program is Calc_main Program of xi2u.                        
//                                Last updated on Jun 22, 2018.  
//======================================================================
void Calc_F(const long L, const long double xi[L], const long N_max, long double F[L], long double IntegF[L]);
void Calc_main(char **argv, const struct SysPara *Mesh)
{
	FILE *fp_o2;
	fp_o2 = fopen(Mesh->name,"w");
	if (fp_o2 == NULL){
		fprintf(stderr,"open error!:  Mesh->name @ Calc_main\n");
		exit(1);
	}
	
	time_t timer_s = time(0);
	clock_t startClock = clock();

	for(long i0 = 0 ; i0 < 3 ; i0++){
		if(Mesh->K[i0] <= 0){
			fprintf(fp_o2,"#\n# Mesh number: L%ld = %ld, k0*Lsys%ld = %.2LE, K%ld = %ld\n", i0, Mesh->L[i0], i0, Mesh->Lsys[i0], i0, Mesh->K[i0]);
		}else{
			fprintf(fp_o2,"#\n# Mesh number: L%ld = %ld, k0*Lsys%ld = %.2LE, K%ld = %ld, k0*Min(du%ld/dxi%ld) = %.2LE\n", i0, Mesh->L[i0], i0, Mesh->Lsys[i0], i0, Mesh->K[i0], i0, i0, Mesh->min_dudxi[i0]);
		}
		fprintf(fp_o2,"# xi%ld,                  k0*u%ld,                k0*du%lddxi%ld\n",i0,i0,i0,i0);
		long double *xipL = calloc(2*Mesh->L[i0],sizeof(long double));//revised on 20180725
		long double *F = calloc(2*Mesh->L[i0],sizeof(long double));
		long double *IntegF = calloc(2*Mesh->L[i0],sizeof(long double));
		if(xipL == NULL || F == NULL || IntegF == NULL){fprintf(stderr,"xipL, F and/or IntegF can not be secured!\n");	exit(EXIT_FAILURE);}
		for(long J0 = 0 ; J0 < Mesh->L[i0] ; J0++) { 
			xipL[J0] = ((long double) J0)/((long double) Mesh->L[i0]);
			xipL[J0 + Mesh->L[i0]] = (0.5L + ((long double) J0))/((long double) Mesh->L[i0]);
		}
		
		Calc_F(2*Mesh->L[i0], xipL, Mesh->K[i0], F, IntegF);
		long double dudxi_last = 0.L;
		for(long J0 = 0 ; J0 < Mesh->L[i0] ; J0++) { 
			long double diff_u = Mesh->Lsys[i0];
			long double Cj = diff_u - ((long double) Mesh->L[i0])*Mesh->min_dudxi[i0];// see eq. (F.4)
//			long double dudxi0 = Mesh->min_dudxi[i0] + (diff_u/((long double) Mesh->L[i0]) - Mesh->min_dudxi[i0])*F[J0];
			long double dudxi0 = Mesh->min_dudxi[i0] + (Cj/((long double) Mesh->L[i0]))*F[J0];// see eq. (F.4)
			if(J0 == 0) {dudxi_last = dudxi0;}
//			long double dudxi1 = Mesh->min_dudxi[i0] + (diff_u/((long double) Mesh->L[i0]) - Mesh->min_dudxi[i0])*F[J0 + Mesh->L[i0]];
			long double dudxi1 = Mesh->min_dudxi[i0] + (Cj/((long double) Mesh->L[i0]))*F[J0 + Mesh->L[i0]];// see eq. (F.4)
			long double u0 = -0.5L*Mesh->Lsys[i0];
			u0 += ((long double) Mesh->L[i0])*xipL[J0]*Mesh->min_dudxi[i0] + Cj*IntegF[J0];
			long double u1 = -0.5L*Mesh->Lsys[i0];
			u1 += ((long double) Mesh->L[i0])*xipL[J0 + Mesh->L[i0]]*Mesh->min_dudxi[i0] + Cj*IntegF[J0 + Mesh->L[i0]];
			fprintf(fp_o2,"%.1f, %.20LE, %.20LE\n", (double)J0, u0, dudxi0);
			fprintf(fp_o2,"%.1f, %.20LE, %.20LE\n", 0.5+(double)J0, u1, dudxi1);
		}
		fprintf(fp_o2,"%.1f, %.20LE, %.20LE\n", (double)Mesh->L[i0], 0.5L*Mesh->Lsys[i0], dudxi_last);
		SAFEFREE(xipL);		SAFEFREE(F);	SAFEFREE(IntegF);
	}
	
	time_t timer_f = time(0);
	double t_span = difftime(timer_f,timer_s);
	clock_t endClock = clock();
	double cpusec = (endClock - startClock)/(double)CLOCKS_PER_SEC;
	fprintf(fp_o2,"# k0 [m^-1] = %.20LE\n", Mesh->k0);//revised on 20180716
	fprintf(fp_o2,"# Input-file = %s\n", argv[1]);
	fprintf(fp_o2,"# This file  = %s\n", Mesh->name);
	fprintf(fp_o2,"# Exe-file   = %s : `%s' was compiled at %s on %s by C-version:%ld\n", argv[0], __FILE__, __TIME__, __DATE__, __STDC_VERSION__);
	fprintf(fp_o2,"# See `Appendix F of Formulation.pdf on 2018/05/20'\n");
	fprintf(fp_o2,"# k0*Lsys0 = %.2LE, k0*Lsys1 = %.2LE, k0*Lsys2 = %.2LE.\n", Mesh->Lsys[0], Mesh->Lsys[1], Mesh->Lsys[2]);
	fprintf(fp_o2,"# K0 = %ld, K1 = %ld, K2 = %ld. ( K: F(xi, K) in Eq. (F.1).)\n", Mesh->K[0], Mesh->K[1], Mesh->K[2]);
	fprintf(fp_o2,"# k0*Min(du0/dxi0) = %.2LE, k0*Min(du1/dxi1) = %.2LE, k0*Min(du2/dxi2) = %.2LE.\n", Mesh->min_dudxi[0], Mesh->min_dudxi[1], Mesh->min_dudxi[2]);
	fprintf(fp_o2,"# <<<<< Note that Min(duj/dxij) is unused when K_j = 0. >>>>>\n");
	fprintf(fp_o2,"# LDBL_DIG   = %u, LDBL_EPSILON = %.20LE\n", __LDBL_DIG__,  __LDBL_EPSILON__);
	fprintf(fp_o2,"# Exec time  = %1.0f sec, CPU time = %1.6f sec\n", t_span, cpusec);
	fprintf(fp_o2,"# Created on %s", ctime(&timer_f));
	fprintf(stderr,"# Created on %s", ctime(&timer_f));
//	fprintf(stderr,"# Check of isspace(null) = %d\n",isspace('\0'));
	fclose(fp_o2);
}
//======================================================================
//
//                Calc_F -- LogMulCal_ini -- Calc_F.  Feb/14/2011
//                                   Last updated on Jul 08, 2012
//       This was modified to long-double version on May 17, 2018
//
//======================================================================
void Generate_G(const long N_max, long double G[N_max+1]);
long double LogMulCal_ini(const long N1);
void Calc_F(const long L, const long double xi[L], const long N_max, long double F[L], long double IntegF[L])
{
	if(N_max == 0){
		for(long ix = 0 ; ix < L ; ix++) {
			F[ix] = 1.;
			IntegF[ix] = xi[ix];
		}
	}
	else{
		long double G[N_max + 1];
		Generate_G(N_max, G);
		for(long ix = 0 ; ix < L ; ix++) {
//			F[ix] = exp(G[0]+((long double) N_max)*log(pow(cos(Pi*xi[ix]),2)));// Eq. (26) in "WageguideModeNum_v3.pdf"
			long double ex_in = G[0]+((long double) N_max*2)*log(fabs(cos(Pi*xi[ix])));//F is calculated by two step, 
			F[ix] = exp(ex_in);// because gcc 4.8.2 for ubuntu can not compile the above command(commented out), although Mingw 4.8.1 is OK.
			IntegF[ix] = 0.;
			for(long J = N_max ; J >= 1 ; J--){
				IntegF[ix] += G[J]*sin(Pi2*((long double) J)*xi[ix])/(Pi2*((long double) J));// Eq. (27) in "WageguideModeNum_v3.pdf"
			}
			IntegF[ix] += xi[ix];
		}
	}
}
void Generate_G(const long N_max, long double G[N_max+1])
{
	G[0] = (((long double) 2*N_max -1 )*log(2.) + LogMulCal_ini(N_max)); // G[0] is logarithm of G.
	
	for(long M = 1 ; M <= N_max ; M++){
		G[M] = 2.;
		for(long I1 = 1 ; I1 <= M ; I1++){
			G[M] *= ((long double) N_max - M + I1)/((long double) N_max + I1);
		}
	}
}
long double LogMulCal_ini(const long N1)
{
	long double Mul = 0.;
	
	for(long J = 1 ; J < N1 ; J++){
		Mul += log(((long double) J)/((long double) N1+J));
	}
	return(Mul);
}

