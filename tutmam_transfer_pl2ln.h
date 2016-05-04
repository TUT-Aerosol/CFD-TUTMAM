/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	Header-file for tutmam_transfer_pl2ln.c
*	Do not change this.
*/
#include  "udf.h"

#ifndef _TUTMAM_TRANSFER_PL2LN_H
#define _TUTMAM_TRANSFER_PL2LN_H

	/* Calculating self-coagulational rate source terms and storing them to sourceTerm vector */
	void calculate_self_coagulational_transfer_rate(real *sourceTerm, cell_t c, Thread *t);
	
	/* Calculating intramodal condensation rate source terms and adding them to sourceTerm vector */
	void calculate_inter_modal_condensation_rate(real *sourceTerm, cell_t c, Thread *t);
	
	/* Storing transfer rate source terms to UDMs */
	void store_transfer_rate(real *sourceTerm, cell_t c, Thread *t);
	
#endif
