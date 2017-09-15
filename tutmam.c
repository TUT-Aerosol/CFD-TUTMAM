/*
*	TUT Modal Aerosol Model for CFD (CFD-TUTMAM)
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	This is the main file for operation of CFD-TUTMAM in Fluent.
*	Do not change this.
*/

#include  "udf.h"
#include  "tutmam_settings.h"
#include  "tutmam_func.h"

/* This initializes UDM and UDS values defined in tutmam_settings.h */
/* Initialization can be done manually from Fluent's GUI. */
DEFINE_ON_DEMAND(tutmam_initialize_manually) 
{
	Domain *mixtureDomain;
    mixtureDomain = Get_Domain(1); /* in single-phase model, the mixture is in domain 1 */
	
	tutmam_initialize(mixtureDomain);
}

/* This initializes UDM and UDS values defined in tutmam_settings.h */
/* Initialization is done automatically as the domain is initialized from Fluent, */
/* if this is hooked into Fluent's initialization hook. */
DEFINE_INIT(tutmam_initialize_automatically,d) 
{
	tutmam_initialize(d);
}

/* This saves distribution variables to UDMs after every iteration. */
/* This must be hooked to Fluent's Execute at end hook. */
DEFINE_EXECUTE_AT_END(tutmam_save_udm_variables_automatically) 
{	
	tutmam_save_udm_variables();
}

/* This shows what UDS and USM slots are mapped. */
/* This can be found in Fluent's execute on demand function list */
/* and in CFD-TUTMAM GUI. */
DEFINE_ON_DEMAND(show_uds_udm_mapping_manually) 
{
#if RP_HOST || !PARALLEL /* serial solver or host process of parallel solver */
	show_uds_udm_mapping();
#endif
}

DEFINE_ON_DEMAND(tutmam_save_udm_variables_manually) 
{	
	tutmam_save_udm_variables();
}

/* These transfer the settings from CFD-TUTMAM GUI to these C-codes. */
/* Settings transfer can be done manually from Fluent's GUI, */
/* but CFD-TUTMAM GUI does this automatically after the settings are changed. */

DEFINE_ON_DEMAND(tutmam_transfer_settings_manually_general_settings) 
{
	tutmam_transfer_settings_general_settings();
}

DEFINE_ON_DEMAND(tutmam_transfer_settings_manually_aerosol_distribution) 
{
	tutmam_transfer_settings_aerosol_distribution();
}

DEFINE_ON_DEMAND(tutmam_transfer_settings_manually_aerosol_process_control) 
{
	tutmam_transfer_settings_aerosol_process_control();
}

DEFINE_ON_DEMAND(tutmam_transfer_settings_manually_diffusion) 
{
	tutmam_transfer_settings_diffusion();
}

DEFINE_ON_DEMAND(tutmam_transfer_settings_manually_nucleation) 
{
	tutmam_transfer_settings_nucleation();
}

DEFINE_ON_DEMAND(tutmam_transfer_settings_manually_condensation) 
{
	tutmam_transfer_settings_condensation();
}

DEFINE_ON_DEMAND(tutmam_transfer_settings_manually_coagulation) 
{
	tutmam_transfer_settings_coagulation();
}

/* This transfers the settings from CFD-TUTMAM GUI to these C-codes automatically, */
/* when the compiled tutmamudf-library is loaded. */
/* No need to hook it to Fluent. */
DEFINE_EXECUTE_ON_LOADING(tutmam_transfer_settings_automatically,tutmamudf) 
{
	tutmam_transfer_settings_general_settings();
	tutmam_transfer_settings_aerosol_distribution();
	tutmam_transfer_settings_aerosol_process_control();
	tutmam_transfer_settings_diffusion();
	tutmam_transfer_settings_nucleation();
	tutmam_transfer_settings_condensation();
	tutmam_transfer_settings_coagulation();
	
	tutmam_construct_settings_matrices();
	
#if RP_HOST || !PARALLEL
	Message("\n************************************************************\n      TUT Modal Aerosol Model for CFD 2.0 libraries loaded      \n************************************************************\n");
#endif
  
}

