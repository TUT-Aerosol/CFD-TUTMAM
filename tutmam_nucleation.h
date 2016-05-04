/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	Header-file for tutmam_nucleation.c
*	Do not change this.
*/
#include  "udf.h"

#ifndef _TUTMAM_NUCLEATION_H
#define _TUTMAM_NUCLEATION_H

	/* Sulfuric acid-water nucleation rate from Vehkamaki et al. (2002,2003) parametrisation */
	void nucleation_rate_vehkamaki(real *sourceTerm, cell_t c, Thread *t);
	
	/* Olin nucleation rate stored to sourceTerm */
	void nucleation_rate_olin(real *sourceTerm, cell_t c, Thread *t);
	
	/* Calculating nucleation rate source terms and storing them to sourceTerm vector */
	void calculate_nucleation_rate(real *sourceTerm, cell_t c, Thread *t);
	
	/* Making under-relaxation to sourceTerm vector */
	void under_relax_nucleation_rate(real *sourceTerm, real uRFNucleation, cell_t c, Thread *t);
	
	/* Storing nucleation rate source terms to UDMs */
	void store_nucleation_rate(real *sourceTerm, cell_t c, Thread *t);
	
	
#endif
