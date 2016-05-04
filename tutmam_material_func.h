/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	Header-file for tutmam_material_func.c
*	Do not change this.
*/
#include  "udf.h"

#ifndef _TUTMAM_MATERIAL_FUNC_H
#define _TUTMAM_MATERIAL_FUNC_H

	
	/* Density of iSpecies in particle phase (kg/m^3) */
	real particle_bulk_density(int iSpecies);
	
	/* Density of particle phase ph (kg/m^3) */
	real particle_phase_density(real temp, int ph, const real *massFractionsInPhase);
	
	/* Density of sulfuric acid-water solution (kg/m^3) */
	real sulfuric_acid_water_density(real temp, const real *massFractionsInPhase);
	
	/* Diffusion coefficient of a gas iFluentSpecies in air (m^2/s) */
	real diffusion_coefficient_gas(real temp, real pressure, int iFluentSpecies);
	
	/* Latent heat of condensation of iSpecies (J/kg) */
	real latent_heat(int iSpecies);
	
	/* Phase activity. 1 or 2 phases only. */
	real phase_activity(real volumeFractionInParticle);
	
	/* Activity coefficient of iSpecies */
	real activity_coefficient(real temp, int iSpecies, const real *moleFractionsInPhase);
	
	/* activity coefficient of iSpecies in sulfuric acid-water solution */
	real activity_coefficient_in_sulfuric_acid_water_solution(real temp, int iSpecies, const real *moleFractionsInPhase);
	
	/* Saturation vapor pressure of iSpecies (Pa) */
	real saturation_vapor_pressure(real temp, int iSpecies);
	
	/* Diameter of vapor molecule (m) */
	real molecule_diameter(int iSpecies);
	
	/* Surface tension (N/m) */
	real surface_tension(real temp, int ph, const real *moleFractionsInPhase);
	
	/* Surface tension of sulfuric acid-water solution (N/m) */
	real sulfuric_acid_water_surface_tension(real temp, const real *moleFractionsInPhase);
	
	/* Equilibrium mole fraction of water in particle */
	real eq_mole_fraction_water(real temp, real rh, real dp);
	
	/* Fraction of condensing HCs in diesel exhaust */
	real diesel_exhaust_hc_fraction(real temp, real partialPressureHC);

#endif
