/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	Header-file for tutmam_diffusion.c
*	Do not change this.
*/
#include  "udf.h"

#ifndef _TUTMAM_DIFFUSION_H
#define _TUTMAM_DIFFUSION_H

	/* Diffusion coefficient: size-inpedependent part (m^3/s) */
	real diff_indep(cell_t c, Thread *t);
	
	/* Diffusion coefficient: size-dependent part (1/m) */
	real diff_dep(cell_t c, Thread *t, real dp);
	
	/* k-th moment averaged diffusion coefficient */
	/* integrated numerically by Gauss-Hermite quadrature method */
	real diffusion_coefficient_uds(cell_t c, Thread *t, int i);
	
	/* k-th moment averaged diffusion coefficient */
	/* calculated by analytical solution of integrals */
	/* using polynomial-parameterized slip correction coefficient */
	real diffusion_coefficient_uds_parametrisation(cell_t c, Thread *t, int i);
	
	/* k-th moment averaged diffusion coefficient */
	/* calculated by analytical solution of integrals */
	/* using constant slip correction coefficient */
	real diffusion_coefficient_uds_constant_slip_correction(cell_t c, Thread *t, int i);

#endif
