/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	Header-file for tutmam_alpha_d_interpolation.c
*	Do not change this.
*/
#include  "udf.h"
#include  "tutmam_settings.h"

#ifndef _TUTMAM_ALPHA_D_INTERPOLATION_H
#define _TUTMAM_ALPHA_D_INTERPOLATION_H

	#if powerLawDistribution > -1
	
		/* 460, 2425, 12300, 49600 */
		#define TUTMAM_INTERP_VEC_LENGTH 2425
		extern real interpAVec[TUTMAM_INTERP_VEC_LENGTH];
		extern real interpBVec[TUTMAM_INTERP_VEC_LENGTH];
		extern real interpAlphaVec[TUTMAM_INTERP_VEC_LENGTH];
		extern real interpDVec[TUTMAM_INTERP_VEC_LENGTH];
	#endif


	/* Finds alpha and d from interpolation table */
	void interpolateAlphaAndD(real *alphaAndD, real A, real B);

#endif
