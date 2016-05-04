/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	Header-file for tutmam_func.c
*	Do not change this.
*/
#include  "udf.h"
#include  "tutmam_settings.h"

#ifndef _TUTMAM_FUNC_H
#define _TUTMAM_FUNC_H

	/* these macros are used in particle distribution parameter saving and acquiring */
	#define C_NTOT(c,t,j) C_UDMI(c,t,iUdmNumberConc+j*nUdmPerMode)
	#define C_CMD(c,t,j) C_UDMI(c,t,iUdmCmd+j*nUdmPerMode)
	#define C_LN2S(c,t,j) C_UDMI(c,t,iUdmLn2s+j*nUdmPerMode)
	
	/* these macros are used in UDS saving and acquiring */
	#define C_UDS_M0(c,t,j) C_UDSI(c,t,iUdsM0+j*nUdsPerMode)
	#define C_UDS_M23(c,t,j) C_UDSI(c,t,iUdsM23+j*nUdsPerMode)
	#define C_UDS_M1(c,t,i,j) C_UDSI(c,t,iUdsM1First+j*nUdsPerMode+i)
	
	/* this macro is used to store and acquire particle temperature (K) */
	#define C_T_PARTICLE(c,t,j) C_UDMI(c,t,iUdmParticleTemperature+j*nUdmPerMode)
	
	/* this macro is used to store and acquire particle density (kg/m^3) */
	#define C_R_PARTICLE(c,t,j) C_UDMI(c,t,iUdmParticleDensity+j*nUdmPerMode)
	
	/* this macro calculates total pressure (Pa) */
	#define C_P_TOT(c,t) (C_P(c,t)+RP_Get_Real("operating-pressure"))
	
	/* this macro is used to store and acquire water equilibrium factor, kappa */
	#define C_WATER_KAPPA(c,t,j) C_UDMI(c,t,iWaterKappaUdm+j*nUdmPerMode)
	
	/* x^3 */
	#define CBC(x) ((x)*(x)*(x))

	/* mole fraction of iSpecies in fluid */
	real C_XI(cell_t c, Thread *t, int iSpecies);
	
	/* mass fraction of iSpecies in particle */
	real C_YI_PARTICLE(cell_t c, Thread *t, int iSpecies, int j);
		
	/* mole fraction of iSpecies in particle */
	real C_XI_PARTICLE(cell_t c, Thread *t, int iSpecies, int j);
	
	/* Mass fraction of ph phase in particle */
	real C_YPH_PARTICLE(cell_t c, Thread *t, int ph, int j);
	
	/* Volume fraction of ph phase in particle */
	real C_VPH_PARTICLE(cell_t c, Thread *t, int ph, int j);
	
	/* Mass fraction of iSpecies in ph phase */
	real C_YI_PHASE(cell_t c, Thread *t, int iSpecies, int ph, int j);
	
	/* Mole fraction of iSpecies in ph phase */
	real C_XI_PHASE(cell_t c, Thread *t, int iSpecies, int ph, int j);
	
	/* Condensed mass fraction of diesel exhaust hydrocarbons */
	real C_YHC_CONDENSED(cell_t c, Thread *t);
	
	/* Particle density (kg/m^3) */
	real particle_density(cell_t c, Thread *t, int j);
	
	/* Particle temperature function (K) */
	real particle_temperature(cell_t c, Thread *t, int j);
	
	/* mean free path of air (Willeke 1976, J. Aerosol Sci. 7, 381-387) */
	real gas_mean_free_path(real temp, real pressure);
	
	/* particle Knudsen number */
	real knudsen_number(real dp, real temp, real pressure);
	
	/* particle slip correction coefficent (Allen & Raabe 1985, Aerosol Sci. Technol. 4, 269-286) */
	real slip_correction_coefficient(real temp, real dp, real pressure);
	
	/* limits value to lowerLimit */
	real tutmam_lower_limit(real, real);
	
	/* limits value to upperLimit */
	real tutmam_upper_limit(real, real);
	
	/* limits value to lowerLimit and upperLimit */
	real tutmam_limits(real, real, real);
	
	/* returns only a positive value */
	real tutmam_positive(real);
	
	/* returns only a negative value */
	real tutmam_negative(real value);
	
	/* abscissas for Gauss-Hermite quadrature numerical integration (efunda.com) */
	real gauss_hermite_abscissas(int i);
	
	/* weights for Gauss-Hermite quadrature numerical integration (efunda.com) */
	real gauss_hermite_weights(int i);
	
	/* abscissas for Gauss-Olin quadrature numerical integration */
	real gauss_olin_abscissas_dp(int i, real D2);
	
	/* weights for Gauss-Olin quadrature numerical integration */
	real gauss_olin_weights(int i, real cc);
	
	/* abscissas for Gauss-Olin quadrature numerical integration, D1 replaced with dCoag */
	real gauss_olin_abscissas_dp_dCoag(int i, real dCoag, real D2);
	
	/* factor for Gauss-Olin quadrature */
	real quadrature_factor_olin(real cc);
	
	/* k-th order of the moment is defined by UDS id */
	real k_moment(int iUds);
	
	/* Calculates ln^2(GSD) from moment concentrations. */
	real calculateLn2s(real numberConc,real surfaceConc,real massConc);
	
	/* Calculates count median diameter (m) from moment concentrations. */
	real calculateCmd(real numberConc,real surfaceConc,real massConc,real rhoParticle);
	
	/* Calculates alpha from moment concentrations. */
	void calculateAlphaAndD2(real *alphaAndD2, real numberConc,real surfaceConc,real massConc,real rhoParticle);
	
	/* Calculates CMD of power law distribution from D1,D2, and alpha */
	real calculateCmdOfPowerLaw(real D2, real alpha);
		
	/* calculates a multiplier for moment averaged D_p for power law distribution */
	real momentAveragingPowerLaw(real alpha, real d, real k, real m);
	
	/* Initialization function that sets initial values for UDSs and UDMs, when initialization is done in Fluent. */
	void tutmam_initialize(Domain *mixtureDomain);
	
	/* The following functions are settings transfer functions that sets the values of */
	/* C-code varibles to the values of RP-variables. */
	/* In parallel solver, C->RP variable transfer is done in host process, */
	/* but C->C variable transfer is done from host to node in node processes. */
	
	/* General settings -settings transfer function */
	void tutmam_transfer_settings_general_settings();
	
	/* Diffusion -settings transfer function */
	void tutmam_transfer_settings_diffusion();
	
	/* Aerosol distribution -settings transfer function */
	void tutmam_transfer_settings_aerosol_distribution();
	
	/* Nucleation -settings transfer function */
	void tutmam_transfer_settings_nucleation();
	
	/* Condensation -settings transfer function */
	void tutmam_transfer_settings_condensation();
	
	/* Aerosol process control -settings transfer function */
	void tutmam_transfer_settings_aerosol_process_control();
	
	/* Coagulation -settings transfer function */
	void tutmam_transfer_settings_coagulation();
	
	/* this converts a string to a real vector */
	void string_to_vector(const char *string, real *vector);
	
	/* this converts a string to an integer vector */
	void string_to_int_vector(const char *string, int *vector, int length);
	
	void message_iCondensingSpecies_vector(int *iCondensingSpeciesVector);
	
	void message_condensationDirectionVector_vector(int *condensationDirectionVector);
	
	/* prints a message of URF vector values and checks that the values are between 0 and 1 */
	void message_URF_vector(real *vector);
	
	/* prints a message of Olin nucleation law parameters */
	void message_Olin_nucleation_parameters_vector(int *iNucleatingSpeciesVector, real *nucleationExponentsVector, real *nucleationSatVapPresExponentsVector, real *nMolecClusterVectorVector);
	
	/* This function saves the particle distribution parameters to UDMs. */
	/* This is called after each iteration. */
	void tutmam_save_udm_variables();

	/* This shows what UDS and UDM slots are mapped. */
	void show_uds_udm_mapping();
	
	/* This construct matrices containing some settings from settings vectors */
	void tutmam_construct_settings_matrices();

#endif
