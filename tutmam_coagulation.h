/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	Header-file for tutmam_coagulation.c
*	Do not change this.
*/
#include  "udf.h"

#ifndef _TUTMAM_COAGULATION_H
#define _TUTMAM_COAGULATION_H

	/* Fuchs-Sutugin correction for coagulation (Fuchs & Sutugin, 1970) */
	real fuchs_sutugin_for_coagulation(real dpi, real dpj, real diffPI, real diffPJ, real temp, real rhoPI, real rhoPJ);
	
	/* Dahneke correction for coagulation (Dahneke, 1983) */
	real dahneke_for_coagulation(real dpi, real dpj, real diffPI, real diffPJ, real temp, real rhoPI, real rhoPJ);
	
	/* Knudsen number for coagulation */
	real knudsen_number_for_coagulation(real dpi, real dpj, real diffPI, real diffPJ, real temp, real rhoPI, real rhoPJ);
	
	/* Coagulation coefficient divided by 2pi (m^3/s) */
	real coagulation_coefficient(real dpi, real dpj, real diffPI, real diffPJ, real fuchsSutugin);
	
	/* Calculating intramodal coagulation rate source terms and storing them to sourceTerm vector */
	void calculate_intra_coagulation_rate(real *sourceTerm, cell_t c, Thread *t, int j);
	
	/* Use robustness coagulation model */
	void use_robustness_coagulation_model(real *sourceTerm, cell_t c, Thread *t, int j);
	
	/* Making under-relaxation to sourceTerm vector */
	void under_relax_intra_coagulation_rate(real *sourceTerm, real uRFCoagulation, cell_t c, Thread *t, int j);
	
	/* Storing intramodal coagulation rate source terms to UDMs */
	void store_intra_coagulation_rate(const real *sourceTerm, cell_t c, Thread *t, int j);
	
	/* Returning intermodal coagulation rate source term */
	/* For cases with loss terms only */
	real inter_coagulation_rate_L(cell_t c, Thread *t, int iUds, int j, int jOther);
	
	/* Returning intermodal coagulation rate source term */
	/* For cases with loss and gain terms */
	real inter_coagulation_rate_LG(cell_t c, Thread *t, int iUds, int j, int jOther);
	
#endif
