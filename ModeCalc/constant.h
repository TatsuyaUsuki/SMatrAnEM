//This header file defines numerical and physical constants for simulation. 
//                                           Last updated on 2018/08/02.
#ifndef CONSTANT_H
#define CONSTANT_H

#define BUFSIZE 1024 // buffer for fgets
//#define ERROR_MIN 1.e-10 // for zgeev_
//#define ERROR_MIN 1.e-13 // for lzgeev_
//#define ERROR_MIN 0. // for dgeev or ldgeev

//#define Pi (acosl(-1.L)) // circular constant
#define Pi (3.1415926535897932384626433832795028841971L) // circular constant
//#define Pi2 (2.L*acosl(-1.L))
#define Pi2 (2.L*3.1415926535897932384626433832795028841971L)

//Physical constants
//See http://physics.nist.gov/cuu/Constants/
//#define BOLTZMANN_CONST 	(1.38064852e-23L) //k_B: Boltzmann constant [J/K]
//#define PLANCK_BAR (1.054571800e-34L) //h^bar: Planck constant over 2*Pi [J*s]
//#define ELECTRON_MASS 	(9.10938356e-31L) //m_e free electron mass [kg]
//#define ELECTRON_CHARGE 	(1.6021766208e-19L) //q_e electron charge [C]

//Exact value in physical constants
#define LIGHT_VELOCITY (2.99792458e8L) //c_0: light velocity in vacuum [m/s] 
//#define MAGNETIC_CONST (0.4L*3.1415926535897932384626433832795028841971e-6L) //mu_0: magnetic constant [N*A^-2]
//#define ELECTRIC_CONST (1.e-10L/(0.4L*3.1415926535897932384626433832795028841971L*2.99792458L*2.99792458L)) //epsilon_0: electric constant [F/m]

#endif /* CONSTANT_H */
