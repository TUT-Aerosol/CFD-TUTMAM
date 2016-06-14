/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	This file defines particle material properties and functions.
*	You can change these to correspond your case.
*/
#include  "udf.h"
#include  "tutmam_material_func.h"
#include  "tutmam_settings.h"
#include  "tutmam_func.h"

/* Density of iSpecies in particle phase (kg/m^3) */
/* You can change the function definition. */
real particle_bulk_density(int iSpecies) {
	real rho;
	
	switch(iSpecies) {
		case 0 :
			rho = 1840.; /* h2so4 */
			break;
			
		case 1 :
			rho = 1000.; /* h2o */
			break;
			
			
		default :
			Error("Species ID not found in particle_bulk_density function\n");
			rho = 1000.;
	}
	
	return rho;
}

/* Density of particle phase ph (kg/m^3) */
/* You can change the function definition. */
real particle_phase_density(real temp, int ph, const real *massFractionsInPhase) {
	real rho;
	
	switch(ph) {
		case 0 :
			rho = sulfuric_acid_water_density(temp,massFractionsInPhase);
			break;
			

		default :
			Error("Phase ID not found in particle_phase_density function\n");
			rho = 1000.;
	}
	
	return rho;
}

/* Density of sulfuric acid-water solution (kg/m^3) */
/* You can change the function definition. */
real sulfuric_acid_water_density(real temp, const real *massFractionsInPhase) {
	real rho;		/* density */
	real aa,bee,cee;/* parameters */
	real massFractionSulfuricAcidInPhase,massFractionSulfuricAcidInPhase2,massFractionSulfuricAcidInPhase3,massFractionSulfuricAcidInPhase4,massFractionSulfuricAcidInPhase5,massFractionSulfuricAcidInPhase6;
	int iSpeciesSulfuricAcid;
	int iFluentSpeciesSulfuricAcid = SV_SpeciesIndex("h2so4"); /* find sulfuric acid ID */
	
	if (iFluentSpeciesSulfuricAcid == -1) {
		Error("h2so4 not found\n");
		return 1000.0;
	}
	iSpeciesSulfuricAcid = tutmamSpeciesIdVector[iFluentSpeciesSulfuricAcid];
	
	massFractionSulfuricAcidInPhase = massFractionsInPhase[iSpeciesSulfuricAcid];
	massFractionSulfuricAcidInPhase2 = massFractionSulfuricAcidInPhase*massFractionSulfuricAcidInPhase;
	massFractionSulfuricAcidInPhase3 = massFractionSulfuricAcidInPhase*massFractionSulfuricAcidInPhase2;
	massFractionSulfuricAcidInPhase4 = massFractionSulfuricAcidInPhase*massFractionSulfuricAcidInPhase3;
	massFractionSulfuricAcidInPhase5 = massFractionSulfuricAcidInPhase*massFractionSulfuricAcidInPhase4;
	massFractionSulfuricAcidInPhase6 = massFractionSulfuricAcidInPhase*massFractionSulfuricAcidInPhase5;
	aa = 0.7681724+2.184714*massFractionSulfuricAcidInPhase+7.163002*massFractionSulfuricAcidInPhase2-44.31447*massFractionSulfuricAcidInPhase3+88.75606*massFractionSulfuricAcidInPhase4-75.73729*massFractionSulfuricAcidInPhase5+23.43228*massFractionSulfuricAcidInPhase6;
	bee = 1.808225e-3-9.294656e-3*massFractionSulfuricAcidInPhase-0.03742148*massFractionSulfuricAcidInPhase2+0.2565321*massFractionSulfuricAcidInPhase3-0.5362872*massFractionSulfuricAcidInPhase4+0.4857736*massFractionSulfuricAcidInPhase5-0.1629592*massFractionSulfuricAcidInPhase6;
	cee = -3.478524e-6+1.335867e-5*massFractionSulfuricAcidInPhase+5.195706e-5*massFractionSulfuricAcidInPhase2-3.717636e-4*massFractionSulfuricAcidInPhase3+7.990811e-4*massFractionSulfuricAcidInPhase4-7.45806e-4*massFractionSulfuricAcidInPhase5+2.58139e-4*massFractionSulfuricAcidInPhase6;

	rho = aa+bee*temp+cee*temp*temp;
	rho = tutmam_limits(600.0,rho*1000.0,2000.0);
	
	return rho;
}

/* Diffusion coefficient of a gas iFluentSpecies in air (m^2/s) */
/* You can change the function definition. */
/* You can use this also in Fluent by choosing diff_gas for UDF in Fluent materials panel. */
real diffusion_coefficient_gas(real temp, real pressure, int iFluentSpecies, real relativeHumidity) {
	real diffCoeff;		/* laminar diffusion coefficient (m^2/s) */

	switch(iFluentSpecies) {
		case 0 :
			/* h2so4, no rh dependence */
			/* diffCoeff = -1.8832e-6 + 2.3024e-8*temp + 2.2366e-11*temp*temp;  */
			
			/* h2so4, rh dependence */
			diffCoeff = 1.8e-9*pow(temp,1.5)/(1.0+0.2876*pow(tutmam_limits(0.0,relativeHumidity,3.0),0.5643)); 
			break;
			
		case 1 :
			diffCoeff = -6.8158e-6 + 8.333e-8*temp + 8.0949e-11*temp*temp; /* h2o */
			break;
			
		case 2 :
			diffCoeff = -6.4141e-6 + 7.1336e-8*temp + 6.9298e-11*temp*temp; /* air */
			break;

		default :
			Error("Species ID not found in diffusion_coefficient_gas function\n");
			diffCoeff = 1.;
	}
	
	return 101325.0/pressure*diffCoeff;
}

/* Latent heat of condensation of iSpecies (J/kg) */
/* You can change the function definition. */
real latent_heat(int iSpecies) {
	real latentHeat;
	
	switch(iSpecies) {
		case 0 :
			latentHeat = 570863.; /* h2so4 */
			break;
			
		case 1 :
			latentHeat = 2257000.; /* h2o */
			break;

			
		default :
			Error("Species ID not found in latent_heat function\n");
			latentHeat = 1.;
	}
	
	return latentHeat;
}

/* Phase activity. 1 or 2 phases only. */
/* You can change the function definition. */
real phase_activity(real volumeFractionInParticle) {
	real phaseActivity;
	
	if (nTutmamPhases > 2) {
		Message("More than 2 phases! Disable phase activity model.");
		return 1.0;
	}
	
	if (volumeFractionInParticle < 0.5) {
		phaseActivity = 0.2371*volumeFractionInParticle+0.5395*sqrt(volumeFractionInParticle);
		
	} else {
		phaseActivity = 1.0-0.2371*(1.0-volumeFractionInParticle)-0.5395*sqrt(1.0-volumeFractionInParticle);
	}

	
	return tutmam_limits(0.0,phaseActivity,1.0);
}

/* Activity coefficient of iSpecies */
/* You can change the function definition. */
real activity_coefficient(real temp, int iSpecies, const real *moleFractionsInPhase) {
	real actCoeff = 0.0;
		
	switch(iSpecies) {
		case 0 :
		case 1 :
			actCoeff = activity_coefficient_in_sulfuric_acid_water_solution(temp,iSpecies,moleFractionsInPhase);
			break;


		default :
			Error("Species ID not found in activity_coefficient function\n");
	}
	
	return actCoeff;
}

/* activity coefficient of iSpecies in sulfuric acid-water solution */
real activity_coefficient_in_sulfuric_acid_water_solution(real temp, int iSpecies, const real *moleFractionsInPhase) {
	real actCoeff;
	real moleFractionSulfuricAcid;
	real moleFractionWater;
	real acta;
	real actb = 0.527;
	real logact;
	
	int iSpeciesSulfuricAcid;
	int iFluentSpeciesSulfuricAcid = SV_SpeciesIndex("h2so4"); /* find sulfuric acid ID */
	
	if (iFluentSpeciesSulfuricAcid == -1) {
		Error("h2so4 not found\n");
		return 1.0;
	}
	iSpeciesSulfuricAcid = tutmamSpeciesIdVector[iFluentSpeciesSulfuricAcid];

	moleFractionSulfuricAcid = moleFractionsInPhase[iSpeciesSulfuricAcid];
	moleFractionWater = 1.0-moleFractionSulfuricAcid;
	
	if (iSpecies == 0) { /* h2so4 */
		acta = 5672.-4.074e6/temp+4.421e8/SQR(temp);
		actb = 1.0/actb;
		logact = acta*SQR(moleFractionWater)/(SQR(moleFractionWater+actb*moleFractionSulfuricAcid)*temp);
		
	} else {	/* h2o */
		acta = 2989.-2.147e6/temp+2.33e8/SQR(temp);
		logact = acta*SQR(moleFractionSulfuricAcid)/(SQR(moleFractionSulfuricAcid+actb*moleFractionWater)*temp);
	}
	
	actCoeff = pow(10.0,logact);
	return actCoeff;
}

/* Saturation vapor pressure of iSpecies (Pa) */
/* You can change the function definition. */
real saturation_vapor_pressure(real temp, int iSpecies) {
	real satVapPress;
	
	switch(iSpecies) {
		case 0 :
			satVapPress = 101325.0*exp(-11.695+10156.0*(1.0/360.15-1.0/temp+0.38/545.0*(1.0+log(360.15/temp)-360.15/temp))); /* h2so4 */
			break;
			
		case 1 :
			satVapPress = exp(77.34491296-7235.424651/temp-8.2*log(temp)+0.0057113*temp); /* h2o */
			break;
			
			
		default :
			Error("Species ID not found in saturation_vapor_pressure function\n");
			satVapPress = 0.0;
	}
	
	return satVapPress;
}

/* Diameter of vapor molecule (m) */
/* You can change the function definition. */
real molecule_diameter(int iSpecies) {
	real d;
	
	switch(iSpecies) {
		case 0 :
			d = 0.527e-9;
			break;
			
		case 1 :
			d = 0.0e-9;
			break;
			
		default :
			Error("Species ID not found in molecule_diameter function\n");
			d = 0.0e-9;
	}
	
	return d;
}

/* Surface tension (N/m) */
/* You can change the function definition. */
real surface_tension(real temp, int ph, const real *moleFractionsInPhase) {
	real surftens;
	
	switch(ph) {
		case 0 :
			surftens = sulfuric_acid_water_surface_tension(temp,moleFractionsInPhase);
			break;


		default :
			Error("Phase ID not found in surface_tension function\n");
			surftens = 0.;
	}
	
	return surftens;
}

/* Surface tension of sulfuric acid-water solution (N/m) */
/* You can change the function definition. */
real sulfuric_acid_water_surface_tension(real temp, const real *moleFractionsInPhase) {
	real surftens;
	real moleFractionSulfuricAcid = 0.0;
	real tempc;
	real surfpara;
	real surfparb;
	real tempone;
	int iSpeciesSulfuricAcid;
	int iFluentSpeciesSulfuricAcid = SV_SpeciesIndex("h2so4"); /* find sulfuric acid ID */
	
	if (iFluentSpeciesSulfuricAcid == -1) {
		Error("h2so4 not found\n");
		return 1.0;
	}
	iSpeciesSulfuricAcid = tutmamSpeciesIdVector[iFluentSpeciesSulfuricAcid];
	
	moleFractionSulfuricAcid = moleFractionsInPhase[iSpeciesSulfuricAcid];

	tempc = 647.15*(1.0-moleFractionSulfuricAcid)*(1.0-moleFractionSulfuricAcid)+900.0*SQR(moleFractionSulfuricAcid)+3156.186*moleFractionSulfuricAcid*(1.0-moleFractionSulfuricAcid);
	
	/*parameters for surface tension*/
	surfpara = 0.2358-0.529*moleFractionSulfuricAcid+4.073*SQR(moleFractionSulfuricAcid)-12.6707*pow(moleFractionSulfuricAcid,3.0)+15.3552*pow(moleFractionSulfuricAcid,4.0)-6.3138*pow(moleFractionSulfuricAcid,5.0);
	surfparb = -0.14738+0.6253*moleFractionSulfuricAcid-5.4808*SQR(moleFractionSulfuricAcid)+17.2366*pow(moleFractionSulfuricAcid,3.0)-21.0487*pow(moleFractionSulfuricAcid,4.0)+8.719*pow(moleFractionSulfuricAcid,5.0);
	tempone = 1.0-temp/tempc;
	
	/*surface tension*/
	surftens = 0.001;
	if (tempone > 0.0) {
		surftens = (surfpara+surfparb*tempone)*pow(tempone,1.256);
	} 	
	if (surftens < 0.001) {
		surftens = 0.001;
	}
	
	return surftens;
}

/* Equilibrium mole fraction of water in particle */
/* Water-sulfuric acid solution. */
/* You can change the function definition. */
real eq_mole_fraction_water(real temp, real rh, real dp) {
	real a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4,A,B,C,xsa,lnrh;
	
	if (temp < 280.) {
		temp = 280.;
	}
	
	if (temp > 400.) {
		temp = 400.;
	}
	
	if (rh < 0.1) {
		rh = 0.1;
	}
	
	lnrh = log(rh);
	
	if (rh > 2.) {
		rh = 2.;
	}
	
	if (dp < 0.4e-9) {
		dp = 0.4e-9;
	}
	
	if (dp > 10.0e-9) {
		dp = 10.0e-9;
	}
	
	a1=-9.942e-23*temp+1.621e-20;
	a2=-4.374e-25*temp*temp-7.751e-23*temp+0.01491e-18;
	a3=-1.025e-24*temp*temp+0.0002947e-18*temp-0.04762e-18;
	a4=-6.501e-25*temp*temp+0.0003273e-18*temp-0.06867e-18;

	b1=6.308e-16*temp*temp-4.887e-14*temp+0.002309e-9;
	b2=3.216e-15*temp*temp-0.0007841e-9*temp+0.1026e-9;
	b3=4.691e-15*temp*temp-0.00169e-9*temp+0.2631e-9;
	b4=2.982e-15*temp*temp-0.001376e-9*temp+0.3312e-9;

	c1=-1.382e-6*temp*temp+0.0005931*temp-0.08432;
	c2=-4.773e-6*temp*temp+0.002017*temp-0.2866;
	c3=-7.658e-6*temp*temp+0.003294*temp-0.5002;
	c4=-1.079e-7*temp*temp+9.526e-5*temp+0.0415;

	
	A=a1*lnrh*lnrh*lnrh+a2*lnrh*lnrh+a3*lnrh+a4;
	B=b1*lnrh*lnrh*lnrh+b2*lnrh*lnrh+b3*lnrh+b4;
	C=c1*lnrh*lnrh*lnrh+c2*lnrh*lnrh+c3*lnrh+c4;

	xsa = A/SQR(dp)+B/dp+C;
	xsa = tutmam_limits(0.0,xsa,1.0);
	
	return 1.0-xsa;
}

/* Fraction of condensing HCs in diesel exhaust */
/* You can change the function definition. */
real diesel_exhaust_hc_fraction(real temp, real partialPressureHC) {

	return tutmam_limits(0.,1.0/(1.0+exp(-5457.0/temp+11.83)*pow(partialPressureHC,-0.7)),1.);
}

