/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	Header-file for tutmam_condensation.c
*	Do not change this.
*/
#include  "udf.h"

#ifndef _TUTMAM_CONDENSATION_H
#define _TUTMAM_CONDENSATION_H

	/* Fuchs-Sutugin correction for heat transfer of a particle (Fuchs & Sutugin, 1970) */
	real fuchs_sutugin_for_heat(real dp, real temp, real pressure);
	
	/* Fuchs-Sutugin correction for mass transfer of a particle (Fuchs & Sutugin, 1970) */
	real fuchs_sutugin_for_mass(real dp, real temp, real pressure, real iSpecies, real diffCoeffGas, real diffCoeffGasDry, real diffCoeffParticle, real dMolecule, real mParticle);
	
	/* Knudsen number for mass transfer of species iSpecies */
	/* Single-component approximation for mean free path */
	real knudsen_number_for_mass(real dp, real temp, real pressure, int iSpecies, real diffCoeffGas, real diffCoeffGasDry, real diffCoeffParticle, real dMolecule, real mParticle);
			
	/* Single particle mass growth rate of species iSpecies (size-independent part) (s^2/m^2) */
	real mass_growth_rate_indep(real temp, int iSpecies);
	
	/* Single particle mass growth rate of species iSpecies (size-dependent part) (kg m^2/s^3) */
	real mass_growth_rate_dep(real dp, real fuchsSutuginForMass, real dMolecule, real diffCoeffGas, real diffCoeffParticle, real moleFractionInFluid, real moleFractionInPhase, real actCoeff, real saturationVaporPressure, real pressure, real kelvinFactor);	
	
	/* Kelvin diameter (m) */
	real kelvin_diameter(real temp, int iSpecies, int j, const real *moleFractionsInPhase, const real *massFractionsInPhase);
	
	/* Calculating nucleation rate source terms and storing them to sourceTerm vector */
	void calculate_condensation_rate(real *sourceTerm, cell_t c, Thread *t, int j);

	/* Making under-relaxation to sourceTerm vector */
	void under_relax_condensation_rate(real *sourceTerm, real uRFCondensationM23, const real *uRFCondensationM1, cell_t c, Thread *t, int j);
	
	/* Storing condensation rate source terms to UDMs */
	void store_condensation_rate(const real *sourceTerm, cell_t c, Thread *t, int j);
	
	/* calculates relative humidity */
	real tutmam_rh(cell_t c, Thread *t, int iSpeciesWater);
	
	/* calculates water equilibrium relative humidity */
	real rh_eq(cell_t c, Thread *t, int iSpeciesWater, int j, const real *M1, real ntot, real ln2s);
	
	/* calculates water equilibrium factor, kappa */
	real water_kappa(cell_t c, Thread *t, int j);
	
	/* calculates water condensation multiplier */
	real condensation_multiplier(cell_t c, Thread *t, int j);
	
#endif
