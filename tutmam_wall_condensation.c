/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	This file defines functions for wall condensation
*	You should not need to change this
*/
#include  "udf.h"
#include  "tutmam_wall_condensation.h"
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
		}

		F_PROFILE(f,thread,iFluentSpecies) = massFraction;
	}
	end_f_loop(f,thread)
}
	