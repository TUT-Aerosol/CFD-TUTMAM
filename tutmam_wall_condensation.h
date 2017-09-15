/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	Header-file for tutmam_wall_condensation.c
*	Do not change this.
*/
#include  "udf.h"

#ifndef _TUTMAM_WALL_CONDENSATION_H
#define _TUTMAM_WALL_CONDENSATION_H

	/* mass fraction of iSpecies taking saturation into account */
	real mass_fraction_wall_condensation(cell_t c, Thread *t, int iSpecies);
	
	/* mass fraction of iSpecies taking activity into account */
	real mass_fraction_activity_wall_condensation(cell_t c, Thread *t, face_t f, Thread *thread, int iSpecies);

#endif
