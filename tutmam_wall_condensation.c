/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	This file defines functions for wall condensation
*	You should not need to change this
*/
#include  "udf.h"
#include  "tutmam_wall_condensation.h"
#include  "tutmam_condensation.h"
#include  "tutmam_constants.h"
#include  "tutmam_settings.h"
#include  "tutmam_func.h"
#include  "tutmam_material_func.h"

/* mass fraction of iSpecies taking saturation into account */
real mass_fraction_wall_condensation(cell_t c, Thread *t, int iSpecies) {
	real pressure;
	real temp;
	real partialPressure;
	real saturationPressure;
	real fractionHC;
	real massFraction = 0.0;
	
	pressure = C_P_TOT(c,t);
	temp = C_T(c,t);
	partialPressure = pressure*C_XI(c,t,iSpecies);
	saturationPressure = saturation_vapor_pressure(temp,iSpecies);
	
	if (partialPressure < saturationPressure) {
		massFraction = C_YI(c,t,fluentSpeciesIdVector[iSpecies]);
	}
	
	if (dieselExhaustHCFractionModel == 1 && iSpecies == dieselExhaustHC) {
		fractionHC = tutmam_limits(0.0,diesel_exhaust_hc_fraction(temp,partialPressure) - C_YHC_CONDENSED(c,t),1.0);
		massFraction = C_YI(c,t,fluentSpeciesIdVector[iSpecies])*(1.0-fractionHC);
	}
	
	return massFraction;
}

/* mass fraction of iSpecies taking activity into account */
real mass_fraction_activity_wall_condensation(cell_t c, Thread *t, face_t f, Thread *thread, int iSpecies) {
	real act;
	real T;
	real rh;
	real xiwall;
	real yiwall = 0.0;
	real moleFractionsInPhase[nTutmamSpecies] = { 0.0 };
	int i = 0;
	int iFluentSpeciesWater = -1;
	int iSpeciesWater = -1;
	
	/* find Fluent species ID for water */
	iFluentSpeciesWater = SV_SpeciesIndex("h2o");
	if (waterEq == 0 || iFluentSpeciesWater == -1) {
		Error("Water equilibrium calculation is off or h2o not found.\nActivity based wall condensation cannot be used!\n");
		return 0.0;
	}
	
	iSpeciesWater = tutmamSpeciesIdVector[iFluentSpeciesWater];
	rh = tutmam_rh(c,t,iSpeciesWater);
	T = C_T(c,t);
	
	if (iSpecies == iSpeciesWater) {
		yiwall = molarMassVector[iSpecies]*rh*saturation_vapor_pressure(F_T(f,thread),iSpecies)/(TUTMAM_R*T*C_R(c,t));
		
	} else {
		xiwall = tutmam_limits(0.0,1.0-exp(0.01429*exp(0.007284*tutmam_limits(200.0,T,600.0))*log(tutmam_limits(0.0001,rh,1.0))),1.0);
		
		moleFractionsInPhase[0] = xiwall;
		moleFractionsInPhase[1] = 1.0 - xiwall;
		
		act = activity_coefficient(T,iSpeciesWater,moleFractionsInPhase)*(1.0 - xiwall);

		while (fabs(rh/act - 1.0) > 0.01 && i < 20) {
			xiwall = xiwall/2.0 * (1.0 + log(rh)/log(act));
			xiwall = tutmam_limits(0.0,xiwall,1.0);
			
			moleFractionsInPhase[0] = xiwall;
			moleFractionsInPhase[1] = 1.0 - xiwall;
			
			act = activity_coefficient(T,iSpeciesWater,moleFractionsInPhase)*(1.0 - xiwall);
			
			i += 1;
		}

		if (i == 20) {
		  xiwall = tutmam_limits(0.0,1.0-exp(0.01429*exp(0.007284*tutmam_limits(200.0,T,600.0))*log(tutmam_limits(0.0001,rh,1.0))),1.0);
		}

		moleFractionsInPhase[0] = xiwall;
		moleFractionsInPhase[1] = 1.0 - xiwall;
		
		yiwall = molarMassVector[iSpecies]*moleFractionsInPhase[iSpecies]*activity_coefficient(T,iSpecies,moleFractionsInPhase)*saturation_vapor_pressure(F_T(f,thread),iSpecies)/(TUTMAM_R*T*C_R(c,t));
	}
	
	return yiwall;
}

/* udf for boundary condition */
DEFINE_PROFILE(wall_condensation, thread, iFluentSpecies) 
{
	int iSpecies;
	cell_t c0;
	Thread *t0 = NULL;
	face_t f;
	real massFraction = 0.0;
	
	if (iFluentSpecies < 0 || iFluentSpecies >= nFluentSpecies) {
		Error("Wall condensation error: iFluentSpecies=%d",iFluentSpecies);
	}
	
	iSpecies = tutmamSpeciesIdVector[iFluentSpecies];

	
	begin_f_loop(f,thread)
	{
		c0 = F_C0(f,thread);
		t0 = F_C0_THREAD(f,thread);
		
		if (wallCondensationLawVector[iSpecies] == 0) {
			massFraction = C_YI(c0,t0,iFluentSpecies);
			
		} else if (wallCondensationLawVector[iSpecies] == 1) {
			massFraction = mass_fraction_wall_condensation(c0,t0,iSpecies);
			
		} else if (wallCondensationLawVector[iSpecies] == 2) {
			massFraction = tutmam_upper_limit(mass_fraction_activity_wall_condensation(c0,t0,f,thread,iSpecies),C_YI(c0,t0,iFluentSpecies));
		}

		F_PROFILE(f,thread,iFluentSpecies) = massFraction; /* wallCondensationLawVector[iSpecies] == 3 (full condensation) */
	}
	end_f_loop(f,thread)
}
	
