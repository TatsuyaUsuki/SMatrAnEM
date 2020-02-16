/* These programs set relative permittivities for materials Jan/20/2012   */
/*                                               revised at Jul/11/2018   */
/*     The programs were coded by Tatsuya Usuki.                          */
/*         Programming language is C99.                                   */
/*         See document file "WageguideModeNum_v3.pdf."                   */
#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define Pi2 6.283185307179586476925286
//#include "constant.h"
//
// double complex epsi ---- double complex epsi_Si                        
//                      |                                                 
//                       -- double complex epsi_SiO2                      
//                      |                                                 
//                       -- double complex epsi_dSi                       
//                      |                                                 
//                       -- double complex epsi_Al     --                 
//                      |                                |                
//                       -- double complex epsi_Cu     --|                
//                      |                                |                
//                       -- double complex epsi_Au     --|                
//                      |                                |-- double high_F
//                       -- double complex epsi_Ag     --|                
//                      |                                |-- double low_F 
//                       -- double complex epsi_Pt     --|                
//                      |                                |                
//                       -- double complex epsi_Ti     --|                
//                      |                                |                
//                       -- double complex epsi_W      --|                
//                      |                                |                
//                       -- double complex epsi_TiN    --                 
//
double complex epsi_Si(const double Lambda, const double Temperature);
double complex epsi_SiO2(const double Lambda, const double Temperature);// 20 deg.C <= Temperature <= 30 deg.C
double complex epsi_dSi(const double Lambda, const double Temperature, const double charge_dens_cm3[2]);// From O-band to U-band,  Temperature: 77,195,273,373 K, charge_dens_cm3[0]: electron charge density, charge_dens_cm3[1]: hole charge density
double complex epsi_Al(const double Lambda);//  1.24 um <= lambda <= 2.07 um
double complex epsi_Cu(const double Lambda);//  1.25 um <= lambda <= 2.0  um
double complex epsi_Au(const double Lambda);//  1.24 um <= lambda <= 1.65 um
double complex epsi_Ag(const double Lambda);//  1.24 um <= lambda <= 2.48 um
double complex epsi_Pt(const double Lambda);//  1.24 um <= lambda <= 1.77 um
double complex epsi_Ti(const double Lambda);//  1.22 um <= lambda <= 1.61 um
double complex epsi_W(const double Lambda); //  1.24 um <= lambda <= 1.77 um
double complex epsi_TiN(const double Lambda);//  826.6 nm <= lambda <= 2479.7 nm
//
double complex epsi(const double Lambda, const double Temperature, const char material[], const double charge_dens_cm3[2])
{
// material variation: SiO2, Si, Al, Cu, Au, Ag, Pt, Ti, W, TiN, Vacuum.
	double complex epsi = - I*10.;
	if (strcmp(material, "SiO2") == 0) {epsi = epsi_SiO2(Lambda, Temperature);}
	else if (strcmp(material, "d-Si") == 0) {epsi = epsi_dSi(Lambda, Temperature, charge_dens_cm3);}
	else if (strcmp(material, "Si") == 0) {epsi = epsi_Si(Lambda, Temperature);}
	else if (strcmp(material, "Al") == 0) {epsi = epsi_Al(Lambda);}
	else if (strcmp(material, "Cu") == 0) {epsi = epsi_Cu(Lambda);}
	else if (strcmp(material, "Au") == 0) {epsi = epsi_Au(Lambda);}
	else if (strcmp(material, "Ag") == 0) {epsi = epsi_Ag(Lambda);}
	else if (strcmp(material, "Pt") == 0) {epsi = epsi_Pt(Lambda);}
	else if (strcmp(material, "Ti") == 0) {epsi = epsi_Ti(Lambda);}
	else if (strcmp(material, "W") == 0) {epsi = epsi_W(Lambda);}
	else if (strcmp(material, "TiN") == 0) {epsi = epsi_TiN(Lambda);}
	else if (strcmp(material, "Vacuum") == 0 || strcmp(material, "vacuum") == 0) {epsi = 1. + I*0.;}
	else {
		fprintf(stderr,"epsi error!:  material = %s\n", material);
		exit(1);
	}
	return(epsi);
}
double complex epsi_SiO2(const double Lambda, const double Temperature)
{
	// SiO2
	// Original: I. H. Malitson, JOSA v55 p1205, 1965.
	// L. Tong et al., Opex v12 p1025, 2004.
	double lambda2 = pow(Lambda*1.e6,2.);// lambda [um]
	double Re_epsi_SiO2  = 0.6961663/(lambda2 - pow(0.0684043,2));
	Re_epsi_SiO2 += 0.4079426/(lambda2 - pow(0.1162414,2));
	Re_epsi_SiO2 += 0.8974794/(lambda2 - pow(9.896161,2));
	Re_epsi_SiO2 *= lambda2;
	Re_epsi_SiO2 += 1.;
	//Section V and Fig.3 in Malitson.
	// 20 deg C <= T <= 30 deg C
	double dndT = (11.5 - 0.7*(Lambda*1.e6 -2.)*(Lambda*1.e6 -2.))*1.e-6;
	dndT *= Temperature - (20.+273.15);
	Re_epsi_SiO2 += dndT*(2.*sqrt(Re_epsi_SiO2) + dndT);

	return(Re_epsi_SiO2+I*0.);
}
double Integral_beta(double b[7], double T1);
double complex epsi_Si(const double Lambda, const double Temperature)// From O-band to U-band,  Temperature: 77,195,273,373 K
{
	// Silicon
	// L. Zhang et al., Nanophotonics vol.3, p.247 (2014).
	// V. Passaro and F. de Passaro, Optical and Quantum Electronics, vol.38, pp.877-888 (2006).
	// (typo) R. M. Osgood et al., Adv. Opt. Photon. vol. 1, pp.162-235 (2009).  
	// Note incorrect formula of Si for L. Tong et al., Opex vol. 12, p.1025 (2004).
	double lambda = Lambda*1.e6;// lambda [um]
	double lambda2 = lambda*lambda;
	double lambda2_2 = 1.1071*1.1071;
	double Re_epsi_Si = 11.6858 + 0.939816/lambda2 + ((8.10461e-3)*lambda2_2)/(lambda2 - lambda2_2);
	if (lambda < 0.5*(1.15+1.31) || lambda > 0.5*(1.53+2.39)){
		fprintf(stderr,"epsi_Si error!:  lambda = %.5E um\n", lambda);
		exit(1);
	}

	// Following formulae from McCaulley et al., PRB v49 p7408, 1994. 
	double b_1p31[7] = {-3.6137e-5, 8.65085e-7, -3.83712e-9, 1.00556e-11, -1.49840e-14, 1.18078e-17,-3.80552e-21};
	double b_1p53[7] = {-3.7239e-5, 8.61435e-7, -3.72908e-9, 0.92278e-11, -1.27065e-14, 0.91077e-17,-2.64153e-21};
	
	double w_1p31 = (1.53 - lambda)/(1.53 - 1.31);
	double w_1p53 = (lambda - 1.31)/(1.53 - 1.31);
	
	double Beta = w_1p31*Integral_beta(b_1p31, Temperature) + w_1p53*Integral_beta(b_1p53, Temperature);
//	fprintf(stderr,"Beta = %0.5E, Integral_beta(b_1p31, Temperature) = %0.5E\n", Beta, Integral_beta(b_1p31, Temperature)); 
	Re_epsi_Si *= exp(2.*Beta);
	
	return(Re_epsi_Si+I*0.);
}
double Integral_beta(double b[7], double T1)
{
	double T0 = 293.;//see Appendix A of L.Zhang et al., Nanophotonics vol.3, p.247 (2014).
	double Beta = 0.;
//	fprintf(stderr,"b[0] = %0.5E, T1 = %0.5E\n", b[0], T1); 
	double TT0 = T0;
	double TT1 = T1;
	for(int i0 = 0 ; i0 < 7 ; i0++){
		Beta += b[i0]*(TT1 - TT0)/((double) (i0+1) );
//		fprintf(stderr,"i0 = %d, Beta = %0.5E, TT1 = %0.5E\n", i0, Beta, TT1); 
		TT1 *= T1;
		TT0 *= T0;
	}
	return(Beta);
}	
double complex epsi_dSi(const double Lambda, const double Temperature, const double charge_dens_cm3[2])// From O-band to U-band,  Temperature: 77,195,273,373 K, charge_dens_cm3[0]: electron charge density, charge_dens_cm3[1]: hole charge density
{
	// L. Tong et al., Opex v12 p1025, 2004.
	// Silicon
	double Re_epsi_Si = creal(epsi_Si(Lambda, Temperature));
	double n_Si = sqrt(Re_epsi_Si);

	// see G.T.Reed et al., nature photonics vol.4, p.518 (2010);
	//     R.A. Soref IEEE v. QE23, p.123 (1987).
	if (charge_dens_cm3[0] < 1e15 && charge_dens_cm3[1] < 1e15){
		fprintf(stderr,"epsi_dSi error!:  charge_dens_cm3[0]=%0.5E, charge_dens_cm3[0]=%0.5E\n", charge_dens_cm3[0], charge_dens_cm3[1]);
		exit(1);
	}
	if (charge_dens_cm3[0] > 1e23 || charge_dens_cm3[1] > 1e23){
		fprintf(stderr,"epsi_dSi error!:  charge_dens_cm3[0]=%0.5E, charge_dens_cm3[0]=%0.5E\n", charge_dens_cm3[0], charge_dens_cm3[1]);
		exit(1);
	}
	
	double Dn_1p55 = -(8.8e-22)*charge_dens_cm3[0] - (8.5e-18)*pow(charge_dens_cm3[1],0.8);
	double Dn_1p3 = -(6.2e-22)*charge_dens_cm3[0] - (6.0e-18)*pow(charge_dens_cm3[1],0.8);

	double Dalpha_1p55 = (8.5e-18)*charge_dens_cm3[0] + (6.0e-18)*charge_dens_cm3[1];// Dalpha[cm^-1]
	double Dalpha_1p3 = (6.0e-18)*charge_dens_cm3[0] + (4.0e-18)*charge_dens_cm3[1]; // Dalpha[cm^-1]

	double Dk_1p55 = Dalpha_1p55*(1.55e-4)/(2.*Pi2);// k = alpha*lambda/(4*pi) from R.A. Soref IEEE v. QE23, p.123 (1987).
	double Dk_1p3 = Dalpha_1p3*(1.3e-4)/(2.*Pi2);// k = alpha*lambda/(4*pi) from R.A. Soref IEEE v. QE23, p.123 (1987).

	double Dn, Dk;
	double lambda = Lambda*1.e6;// lambda [um]
	if (lambda < 1.2599 || lambda > 1.6751){
		fprintf(stderr,"epsi_dSi error!:  lambda < O-band (1.26 um) or > U-band (1.675 um)\n");
		exit(1);
	}
	else{
		double w_long_wavelength = (lambda - 1.3)/(1.55 - 1.3);
		double w_short_wavelength = (1.55 - lambda)/(1.55 - 1.3);
		Dn = Dn_1p55*w_long_wavelength + Dn_1p3*w_short_wavelength;
		Dk = Dk_1p55*w_long_wavelength + Dk_1p3*w_short_wavelength;
	}
	return((n_Si+Dn+I*Dk)*(n_Si+Dn+I*Dk));
}
double high_F(const double x, const double long_lambda, const double short_lambda);
double low_F(const double x, const double long_lambda, const double short_lambda);
double complex epsi_Al(const double Lambda)
{
	// M.A. Ordal et al., Appl. Opt. v22 p1099, 1983.
	// Table 1.(E. Shiles et al., PRB v22 p1612, 1980)
	// Aluminum
	double lambda = Lambda*1.e6;// lambda [um]
	double Re_epsi_Al = 1.;
	double Im_epsi_Al = -1.;
	if (lambda < 1.24 || lambda > 2.07){
		fprintf(stderr,"epsi_Al error!:  lambda < 1.24 um or > 2.07 um\n");
		exit(1);
	}
	else{
		if (lambda <= 2.07 && lambda >= 1.55){
			Re_epsi_Al = -4.53e2*low_F(lambda, 2.07, 1.55) - 2.52e2*high_F(lambda, 2.07, 1.55);
			Im_epsi_Al =  9.73e1*low_F(lambda, 2.07, 1.55) + 4.61e1*high_F(lambda, 2.07, 1.55);
		}
		if (lambda <= 1.55 && lambda >= 1.24){
			Re_epsi_Al = -2.52e2*low_F(lambda, 1.55, 1.24) - 1.54e2*high_F(lambda, 1.55, 1.24);
			Im_epsi_Al =  4.61e1*low_F(lambda, 1.55, 1.24) + 3.02e1*high_F(lambda, 1.55, 1.24);
		}
		if (Re_epsi_Al >= 0. || Im_epsi_Al <= 0.){
			fprintf(stderr,"epsi_Al error!:  program error\n");
			exit(1);
		}
	}
	return(Re_epsi_Al + I*Im_epsi_Al);
}
double high_F(const double x, const double long_lambda, const double short_lambda)
{
	if (x <= long_lambda && x >= short_lambda){
		double a_highF = (long_lambda - x)/(long_lambda - short_lambda);
		return(a_highF);
	}
	else{
		fprintf(stderr,"high_F error!:  x > long_lambda or x < short_lambda\n");
		exit(1);
	}
}
double low_F(const double x, const double long_lambda, const double short_lambda)
{
	if (x <= long_lambda && x >= short_lambda){
		double a_lowF  = (x - short_lambda)/(long_lambda - short_lambda);
		return(a_lowF);
	}
	else{
		fprintf(stderr,"low_F error!:  x > long_lambda or x < short_lambda\n");
		exit(1);
	}
}
double complex epsi_Cu(const double Lambda)
{
	// M.A. Ordal et al., Appl. Opt. v22 p1099, 1983.
	// Table 2.(B. Dold et al., Optik v22 p435, 1965)
	// Copper
	double lambda = Lambda*1.e6;// lambda [um]
	double Re_epsi_Cu = 1.;
	double Im_epsi_Cu = -1.;
	if (lambda < 1.25 || lambda > 2.0){
		fprintf(stderr,"epsi_Cu error!:  lambda < 1.25 um or > 2.0 um\n");
		exit(1);
	}
	else{
		if (lambda <= 2.0 && lambda >= 1.5){
			Re_epsi_Cu = -1.12e2*low_F(lambda, 2.0, 1.5) - 6.37e1*high_F(lambda, 2.0, 1.5);
			Im_epsi_Cu =  1.80e1*low_F(lambda, 2.0, 1.5) + 9.28e0*high_F(lambda, 2.0, 1.5);
		}
		if (lambda <= 1.5 && lambda >= 1.25){
			Re_epsi_Cu = -6.37e1*low_F(lambda, 1.5, 1.25) - 4.46e1*high_F(lambda, 1.5, 1.25);
			Im_epsi_Cu =  9.28e0*low_F(lambda, 1.5, 1.25) + 6.57e0*high_F(lambda, 1.5, 1.25);
		}
		if (Re_epsi_Cu >= 0. || Im_epsi_Cu <= 0.){
			fprintf(stderr,"epsi_Cu error!:  program error\n");
			exit(1);
		}
	}
	return(Re_epsi_Cu + I*Im_epsi_Cu);
}
double complex epsi_Au(const double Lambda)
{
	// M.A. Ordal et al., Appl. Opt. v22 p1099, 1983.
	// Table 3.(J.H. Weaver et al. "Physics Data, Optical Properties of Metals," Fach-Information Zentrum, Kalsrube, FOR, 1981)
	// Gold
	double lambda = Lambda*1.e6;// lambda [um]
	double Re_epsi_Au = 1.;
	double Im_epsi_Au = -1.;
	if (lambda < 1.24 || lambda > 1.65){
		fprintf(stderr,"epsi_Au error!:  lambda < 1.24 um or > 1.65 um\n");
		exit(1);
	}
	else{
		if (lambda <= 1.65 && lambda >= 1.55){
			Re_epsi_Au = -1.19e2*low_F(lambda, 1.65, 1.55) - 1.04e2*high_F(lambda, 1.65, 1.55);
			Im_epsi_Au =  4.15e0*low_F(lambda, 1.65, 1.55) + 3.68e0*high_F(lambda, 1.65, 1.55);
		}
		if (lambda <= 1.55 && lambda >= 1.46){
			Re_epsi_Au = -1.04e2*low_F(lambda, 1.55, 1.46) - 9.16e1*high_F(lambda, 1.55, 1.46);
			Im_epsi_Au =  3.68e0*low_F(lambda, 1.55, 1.46) + 3.06e0*high_F(lambda, 1.55, 1.46);
		}
		if (lambda <= 1.46 && lambda >= 1.38){
			Re_epsi_Au = -9.16e1*low_F(lambda, 1.46, 1.38) - 8.12e1*high_F(lambda, 1.46, 1.38);
			Im_epsi_Au =  3.06e0*low_F(lambda, 1.46, 1.38) + 2.70e0*high_F(lambda, 1.46, 1.38);
		}
		if (lambda <= 1.38 && lambda >= 1.31){
			Re_epsi_Au = -8.12e1*low_F(lambda, 1.38, 1.31) - 7.21e1*high_F(lambda, 1.38, 1.31);
			Im_epsi_Au =  2.70e0*low_F(lambda, 1.38, 1.31) + 2.38e0*high_F(lambda, 1.38, 1.31);
		}
		if (lambda <= 1.31 && lambda >= 1.24){
			Re_epsi_Au = -7.21e1*low_F(lambda, 1.31, 1.24) - 6.45e1*high_F(lambda, 1.31, 1.24);
			Im_epsi_Au =  2.38e0*low_F(lambda, 1.31, 1.24) + 2.09e0*high_F(lambda, 1.31, 1.24);
		}
		if (Re_epsi_Au >= 0. || Im_epsi_Au <= 0.){
			fprintf(stderr,"epsi_Au error!:  program error\n");
			exit(1);
		}
	}
	return(Re_epsi_Au + I*Im_epsi_Au);
}
double complex epsi_Ag(const double Lambda)
{
	// M.A. Ordal et al., Appl. Opt. v22 p1099, 1983.
	// Table 5.(H.J. Hageman et al., J.Opt.Soc.Am. v65 p742, 1975)
	// Silver
	double lambda = Lambda*1.e6;// lambda [um]
	double Re_epsi_Ag = 1.;
	double Im_epsi_Ag = -1.;
	if (lambda < 1.24 || lambda > 2.48){
		fprintf(stderr,"epsi_Ag error!:  lambda < 1.24 um or > 2.48 um\n");
		exit(1);
	}
	else{
		if (lambda <= 2.48 && lambda >= 1.24){
			Re_epsi_Ag = -3.35e2*low_F(lambda, 2.48, 1.24) - 8.15e1*high_F(lambda, 2.48, 1.24);
			Im_epsi_Ag =  2.45e1*low_F(lambda, 2.48, 1.24) + 5.06e0*high_F(lambda, 2.48, 1.24);
		}
		if (Re_epsi_Ag >= 0. || Im_epsi_Ag <= 0.){
			fprintf(stderr,"epsi_Ag error!:  program error\n");
			exit(1);
		}
	}
	return(Re_epsi_Ag + I*Im_epsi_Ag);
}
double complex epsi_Pt(const double Lambda)
{
	// M.A. Ordal et al., Appl. Opt. v22 p1099, 1983.
	// Table 10.(J.H. Weaver et al., PRB v11 p1416, 1975)
	// Platinum
	double lambda = Lambda*1.e6;// lambda [um]
	double Re_epsi_Pt = 1.;
	double Im_epsi_Pt = -1.;
	if (lambda < 1.24 || lambda > 1.77){
		fprintf(stderr,"epsi_Pt error!:  lambda < 1.24 um or > 1.77 um\n");
		exit(1);
	}
	else{
		if (lambda <= 1.77 && lambda >= 1.55){
			Re_epsi_Pt = -1.40e1*low_F(lambda, 1.77, 1.55) - 2.14e1*high_F(lambda, 1.77, 1.55);
			Im_epsi_Pt =  7.80e1*low_F(lambda, 1.77, 1.55) + 7.48e1*high_F(lambda, 1.77, 1.55);
		}
		if (lambda <= 1.55 && lambda >= 1.24){
			Re_epsi_Pt = -2.14e1*low_F(lambda, 1.55, 1.24) - 2.58e1*high_F(lambda, 1.55, 1.24);
			Im_epsi_Pt =  7.48e1*low_F(lambda, 1.55, 1.24) + 5.63e1*high_F(lambda, 1.55, 1.24);
		}
		if (Re_epsi_Pt >= 0. || Im_epsi_Pt <= 0.){
			fprintf(stderr,"epsi_Pt error!:  program error\n");
			exit(1);
		}
	}
	return(Re_epsi_Pt + I*Im_epsi_Pt);
}
double complex epsi_Ti(const double Lambda)
{
	// M.A. Ordal et al., Appl. Opt. v22 p1099, 1983.
	// Table 11.(P.B. Johnson et al., PRB v9 p5056, 1974)
	// Titanium
	double lambda = Lambda*1.e6;// lambda [um]
	double Re_epsi_Ti = 1.;
	double Im_epsi_Ti = -1.;
	if (lambda < 1.22 || lambda > 1.61){
		fprintf(stderr,"epsi_Ti error!:  lambda < 1.22 um or > 1.61 um\n");
		exit(1);
	}
	else{
		if (lambda <= 1.61 && lambda >= 1.39){
			Re_epsi_Ti = -8.47e0*low_F(lambda, 1.61, 1.39) - 5.63e0*high_F(lambda, 1.61, 1.39);
			Im_epsi_Ti =  3.47e1*low_F(lambda, 1.61, 1.39) + 3.21e1*high_F(lambda, 1.61, 1.39);
		}
		if (lambda <= 1.39 && lambda >= 1.22){
			Re_epsi_Ti = -5.63e0*low_F(lambda, 1.39, 1.22) - 4.12e0*high_F(lambda, 1.39, 1.22);
			Im_epsi_Ti =  3.21e1*low_F(lambda, 1.39, 1.22) + 3.00e1*high_F(lambda, 1.39, 1.22);
		}
		if (Re_epsi_Ti >= 0. || Im_epsi_Ti <= 0.){
			fprintf(stderr,"epsi_Ti error!:  program error\n");
			exit(1);
		}
	}
	return(Re_epsi_Ti + I*Im_epsi_Ti);
}
double complex epsi_W(const double Lambda)
{
	// M.A. Ordal et al., Appl. Opt. v22 p1099, 1983.
	// Table 12.(J.H. Weaver et al., PRB v12 p1293, 1975)
	// Tungsten
	double lambda = Lambda*1.e6;// lambda [um]
	double Re_epsi_W = 1.;
	double Im_epsi_W = -1.;
	if (lambda < 1.24 || lambda > 1.77){
		fprintf(stderr,"epsi_W error!:  lambda < 1.24 um or > 1.77 um\n");
		exit(1);
	}
	else{
		if (lambda <= 1.77 && lambda >= 1.51){
			Re_epsi_W = -3.50e1*low_F(lambda, 1.77, 1.51) - 1.57e1*high_F(lambda, 1.77, 1.51);
			Im_epsi_W =  1.95e1*low_F(lambda, 1.77, 1.51) + 2.18e1*high_F(lambda, 1.77, 1.51);
		}
		if (lambda <= 1.51 && lambda >= 1.38){
			Re_epsi_W = -1.57e1*low_F(lambda, 1.51, 1.38) - 1.00e1*high_F(lambda, 1.51, 1.38);
			Im_epsi_W =  2.18e1*low_F(lambda, 1.51, 1.38) + 2.76e1*high_F(lambda, 1.51, 1.38);
		}
		if (lambda <= 1.38 && lambda >= 1.24){
			Re_epsi_W = -1.00e1*low_F(lambda, 1.38, 1.24) - 8.80e0*high_F(lambda, 1.38, 1.24);
			Im_epsi_W =  2.76e1*low_F(lambda, 1.38, 1.24) + 2.71e1*high_F(lambda, 1.38, 1.24);
		}
		if (Re_epsi_W >= 0. || Im_epsi_W <= 0.){
			fprintf(stderr,"epsi_W error!:  program error\n");
			exit(1);
		}
	}
	return(Re_epsi_W + I*Im_epsi_W);
}
double complex epsi_TiN(const double Lambda)
{
	// http://www.filmetricsinc.jp/refractive-index-database/TiN/Titanium-nitride
	// Table 11.(P.B. Johnson et al., PRB v9 p5056, 1974)
	// Titanium
	double lambda = Lambda*1.e9;// lambda [nm]
	double Re_n_TiN = -1.;
	double Im_n_TiN = -1.;
	if (lambda < 826.6 || lambda > 2479.7){
		fprintf(stderr,"epsi_TiN error!:  lambda < 826.6 nm or > 2479.7 nm\n");
		exit(1);
	}
	else{
		if (lambda <= 2479.7 && lambda >= 1239.9){
			Re_n_TiN = 4.52*low_F(lambda, 2479.7, 1239.9) + 2.69*high_F(lambda, 2479.7, 1239.9);
			Im_n_TiN = 8.25*low_F(lambda, 2479.7, 1239.9) + 5.04*high_F(lambda, 2479.7, 1239.9);
		}
		if (lambda <= 1239.9 && lambda >= 826.6){
			Re_n_TiN = 2.69*low_F(lambda, 1239.9, 826.6) + 1.82*high_F(lambda, 1239.9, 826.6);
			Im_n_TiN = 5.04*low_F(lambda, 1239.9, 826.6) + 3.81*high_F(lambda, 1239.9, 826.6);
		}
		if (Re_n_TiN <= 0. || Im_n_TiN <= 0.){
			fprintf(stderr,"epsi_TiN error!:  program error\n");
			exit(1);
		}
	}
	return(pow(Re_n_TiN + I*Im_n_TiN,2));
}
