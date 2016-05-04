/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	Global constants file
*	No need to change these
*/
#include  "udf.h"

#ifndef _TUTMAM_CONSTANTS_H
#define _TUTMAM_CONSTANTS_H

	/* 1/3 */
	#define TUTMAM_13 0.333333333333333
	
	/* 2/3 */
	#define TUTMAM_23 0.666666666666667
	
	/* pi/6 */
	#define TUTMAM_PI 3.141592653589793

	/* pi/6 */
	#define TUTMAM_PI6 0.523598775598299
	
	/* pi/4 */
	#define TUTMAM_PI4 0.785398163397448
	
	/* pi/8 */
	#define TUTMAM_PI8 0.392699081698724
	
	/* 2*pi */
	#define TUTMAM_2PI 6.283185307179586
	
	/* 3*pi */
	#define TUTMAM_3PI 9.424777960769379
	
	/* sqrt(pi) */
	#define TUTMAM_SQRTPI 1.772453850905516
	
	/* (pi/6)^(-1/3) */
	#define TUTMAM_PI6M13 1.2407009817988
	
	/* 2/3*(pi/6)^(-1/3) */
	#define TUTMAM_23PI6M13 0.827133987865867
	
	/* (pi/6)^(2/3) */
	#define TUTMAM_PI623 0.649629514953459
	
	/* Boltzmann constant (J/K) */
	#define TUTMAM_BOLTZ 1.3806488e-23
	
	/* Universal gas constant (J/K mol) */
	#define TUTMAM_R 8.3145
	
	/* Avogadro constant (1/mol) */
	#define TUTMAM_AVOGADRO (TUTMAM_R/TUTMAM_BOLTZ)
	
#endif
