/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	This file contains internal functions for CFD-TUTMAM model.
*	Do not change this.
*/
#include  "udf.h"
#include  "tutmam_constants.h"
#include  "tutmam_settings.h"
#include  "tutmam_func.h"
#include  "tutmam_material_func.h"
#include  "tutmam_condensation.h"
#include  "tutmam_alpha_d_interpolation.h"
#include  "tutmam_alpha_d_iteration.h"

/* Mole fraction of iSpecies in fluid. */
real C_XI(cell_t c, Thread *t, int iSpecies) {
	real massFraction;	/* mass fraction of iSpecies in fluid */
	real molarMass;
	int iFluentSpecies;	/* fluent ID number for iSpecies */
	int i;
	Material *mixMat;
	Material *speMat = NULL;
	real down = 0.0;
	real up = 0.0;

	iFluentSpecies = fluentSpeciesIdVector[iSpecies];
	mixMat = mixture_material(Get_Domain(1));
	
	mixture_species_loop(mixMat, speMat, i)
	{
		massFraction = C_YI(c,t,i);
		molarMass = MATERIAL_PROP(speMat,PROP_mwi);
		
		if (i == iFluentSpecies) {
			up = massFraction/molarMass;
		}
		
		down += massFraction/molarMass;
	}

	return tutmam_limits(0.0,up/down,1.0);
}

/* Mass fraction of iSpecies in particle */
real C_YI_PARTICLE(cell_t c, Thread *t, int iSpecies, int j) {
	real iSpeciesM1;	/* moment 1 of iSpecies (unitless) */
	int i;				/* index for summation */
	real totalM1;		/* total moment 1 (unitless) */

	totalM1 = 0.;
	for (i = 0; i < nTutmamSpecies; ++i) { /* looping over all particle species */
		totalM1 += tutmam_positive(C_UDS_M1(c,t,i,j)); /* summing the mass concentrations */
	}
	
	if (totalM1 < minMassConc) {
		return 1.0/nTutmamSpecies;
	}
	
	iSpeciesM1 = tutmam_positive(C_UDS_M1(c,t,iSpecies,j));
	
	return tutmam_limits(0.0,iSpeciesM1/totalM1,1.0); /* mass fraction */
}

/* Mole fraction of iSpecies in particle */
real C_XI_PARTICLE(cell_t c, Thread *t, int iSpecies, int j) {
	real iSpeciesM1;	/* moment 1 of iSpecies in moles (mol/g) */
	int i;				/* index for summation */
	real totalM1;		/* total moment 1 (mol/g) in moles */

	totalM1 = 0.;
	for (i = 0; i < nTutmamSpecies; ++i) { /* looping over all particle species */
		totalM1 += tutmam_positive(C_UDS_M1(c,t,i,j))/molarMassVector[i]; /* summing the mole concentrations */
	}
	
	if (totalM1 < minMassConc*0.029) {
		return 0.;
	}
	
	iSpeciesM1 = tutmam_positive(C_UDS_M1(c,t,iSpecies,j))/molarMassVector[iSpecies];
	
	return tutmam_limits(0.0,iSpeciesM1/totalM1,1.0); /* mole fraction */
}

/* Mass fraction of ph phase in particle */
real C_YPH_PARTICLE(cell_t c, Thread *t, int ph, int j) {
	real phM1 = 0.0;	/* moment 1 of ph (unitless) */
	int i;				/* index for summation */
	real totalM1 = 0.0;	/* total moment 1 (unitless) */

	for (i = 0; i < nTutmamSpecies; ++i) { /* looping over all particle species */
		totalM1 += tutmam_positive(C_UDS_M1(c,t,i,j)); /* summing the mass concentrations */
	}
	
	if (totalM1 < minMassConc) {
		return 1.0/nTutmamPhases;
	}
	
	for (i = 0; i < nTutmamSpecies; ++i) { /* looping over all particle species */
		if (phaseMatrix[ph][i] == 1) {
			phM1 += tutmam_positive(C_UDS_M1(c,t,i,j));
		}
	}
	
	return tutmam_limits(0.0,phM1/totalM1,1.0); /* mass fraction */
}

/* Volume fraction of ph phase in particle */
real C_VPH_PARTICLE(cell_t c, Thread *t, int ph, int j) {
	real phM1 = 0.0;	/* moment 1 of ph (unitless) */
	int i;				/* index for summation */
	real totalM1 = 0.0;	/* total moment 1 (unitless) */
	real massFractionsInPhase[nTutmamSpecies];
	real temp;
	
	temp = C_T(c,t);

	for (i = 0; i < nTutmamSpecies; ++i) { /* looping over all particle species */
		totalM1 += tutmam_positive(C_UDS_M1(c,t,i,j)); /* summing the mass concentrations */
		
		if (phaseMatrix[ph][i] == 1) {
			phM1 += tutmam_positive(C_UDS_M1(c,t,i,j));
		}
		
		massFractionsInPhase[i] = C_YI_PHASE(c,t,i,ph,j);
	}
	
	if (totalM1 < minMassConc) {
		return 1.0/nTutmamPhases;
	}
			
	return tutmam_limits(0.0,phM1/totalM1*C_R_PARTICLE(c,t,j)/particle_phase_density(temp,ph,massFractionsInPhase),1.0); /* volume fraction */
}

/* Mass fraction of iSpecies in ph phase */
real C_YI_PHASE(cell_t c, Thread *t, int iSpecies, int ph, int j) {
	real iSpeciesM1;	/* moment 1 of iSpecies (unitless) */
	int i;				/* index for summation */
	real totalM1;		/* total moment 1 (unitless) */

	totalM1 = 0.;
	for (i = 0; i < nTutmamSpecies; ++i) { /* looping over all particle species */
		if (phaseMatrix[ph][i] == 1) {
			totalM1 += tutmam_positive(C_UDS_M1(c,t,i,j)); /* summing the mass concentrations */
		}
	}
	
	if (totalM1 < minMassConc) {
		return 1.0/nTutmamSpecies;
	}
	
	iSpeciesM1 = tutmam_positive(C_UDS_M1(c,t,iSpecies,j));
	
	return tutmam_limits(0.0,iSpeciesM1/totalM1,1.0); /* mass fraction */
}

/* Mole fraction of iSpecies in ph phase */
real C_XI_PHASE(cell_t c, Thread *t, int iSpecies, int ph, int j) {
	real iSpeciesM1;	/* moment 1 of iSpecies in moles (mol/g) */
	int i;				/* index for summation */
	real totalM1;		/* total moment 1 (mol/g) in moles */

	totalM1 = 0.;
	for (i = 0; i < nTutmamSpecies; ++i) { /* looping over all particle species */
		if (phaseMatrix[ph][i] == 1) {
			totalM1 += tutmam_positive(C_UDS_M1(c,t,i,j))/molarMassVector[i]; /* summing the mole concentrations */
		}
	}
	
	if (totalM1 < minMassConc) {
		return 1.0/nTutmamSpecies;
	}
	
	iSpeciesM1 = tutmam_positive(C_UDS_M1(c,t,iSpecies,j))/molarMassVector[iSpecies];
	
	return tutmam_limits(0.0,iSpeciesM1/totalM1,1.0); /* mole fraction */
}

/* Condensed mass fraction of diesel exhaust hydrocarbons */
real C_YHC_CONDENSED(cell_t c, Thread *t) {
	real hcInParticles = 0.0;
	real hcInFluid;
	real hcTotal;
	int j;
	
	hcInFluid = C_YI(c,t,fluentSpeciesIdVector[dieselExhaustHC]);
	
	for (j = 0; j < nTutmamModes; ++j) {
		hcInParticles += tutmam_positive(C_UDS_M1(c,t,dieselExhaustHC,j));
	}
	
	hcTotal = hcInFluid+hcInParticles;
	
	if (hcTotal < minMassConc) {
		return 0.0;
	}
	
	return tutmam_limits(0.0,hcInParticles/hcTotal,1.0);
}

/* Particle density function (kg/m^3) */
real particle_density(cell_t c, Thread *t, int j) {
	real denominator;		/* denominator in the last equation (m^3/kg) */
	real massFractionPhase;	/* mass fraction of ph phase in particle */
	real rhoPhase;			/* density of pure ph phase in particle phase (kg/m^3) */
	int i,ph;				/* integers used in summation */
	real temp;				/* temperature (K) */
	real massFractionsInPhase[nTutmamSpecies]; /* mass fractions of species in a phase */
	
	temp = C_T(c,t);
	
	denominator = 0.0;	/* summation starts from zero */
	for (ph = 0; ph < nTutmamPhases; ++ph) { /* looping over all particle species */
		for (i = 0; i < nTutmamSpecies; ++i) {
			massFractionsInPhase[i] = C_YI_PHASE(c,t,i,ph,j);
		}
	
		massFractionPhase = C_YPH_PARTICLE(c,t,ph,j);
		rhoPhase = particle_phase_density(temp,ph,massFractionsInPhase);
		denominator += massFractionPhase/rhoPhase; /* (m^3/kg) */
	}
	
	if (denominator < minMassConc) {
		return 1000.;
	}
	
	return 1.0/denominator;	/* (kg/m^3) */
}

/* Particle temperature function (K) */
real particle_temperature(cell_t c, Thread *t, int j) {
	real tempFluid;					/* fluid temperature (K) */
	real tempParticle;				/* particle temperature (K) */
	real tempDifference;			/* particle temperature minus fluid temperature (K) */
	real latentHeat;				/* latent heat of condensation (J/kg) */
	real sourceTerm;				/* condensation source term (kg/(s m^3)) */
	real latentHeatTimesSourceTerm;	/* latentHeat*sourceTerm (J/(s m^3)) */
	real cmd;						/* particle distribution CMD (m) */
	real thermCondFluid;			/* thermal conductivity of fluid (W/(m K)) */
	real numberConc;				/* particle number concentration (#/m^3) */
	real fuchsSutugin;				/* Fuchs-Sutugin correction factor */
	real pressure;					/* total pressure (Pa) */
	int iSpecies;					/* species ID number */
	
	tempFluid = C_T(c,t);
	numberConc = C_NTOT(c,t,j);
	pressure = C_P_TOT(c,t);
	
	if (numberConc < minNumberConc) {	/* preventing division by zero */
		return tempFluid;
	}
	
	cmd = C_CMD(c,t,j);
	if (j == powerLawDistribution) { /* if power law */
		cmd = calculateCmdOfPowerLaw(cmd,C_LN2S(c,t,j));
	}
	
	thermCondFluid = C_K_L(c,t);
	fuchsSutugin = fuchs_sutugin_for_heat(cmd,tempFluid,pressure);
	
	latentHeatTimesSourceTerm = 0.0;
	for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) {	/* summing over all species */
		latentHeat = latent_heat(iSpecies);
		sourceTerm = C_UDMI(c,t,iUdmCondensationM1First+j*nUdmPerMode+iSpecies);
		
		latentHeatTimesSourceTerm += latentHeat*sourceTerm;
	}
	
	tempDifference = latentHeatTimesSourceTerm/(TUTMAM_2PI*cmd*thermCondFluid*numberConc*fuchsSutugin);
	tempDifference = tutmam_limits(0.0-maxParticleTemperatureDifference,tempDifference,maxParticleTemperatureDifference);
	
	tempParticle = tempFluid + tempDifference;
	
	return tempParticle;
}

/* mean free path of air (m) (Willeke 1976, J. Aerosol Sci. 7, 381-387) */
real gas_mean_free_path(real temp, real pressure) {
	real su = 110.4; 			/* Sutherland constant for air (K) */
	real temp_r = 293.15;		/* Reference temperature (K) */
	real l_r = 66.5e-9;			/* Mean free path in reference temperature (m) */
	real p_0 = 101325.0;		/* atmospheric pressure (Pa) */
	
	return p_0/pressure*l_r*temp/temp_r*(1.0+su/temp_r)/(1.0+su/temp);
}

/* particle Knudsen number */
real knudsen_number(real dp, real temp, real pressure) {
	real lambda = gas_mean_free_path(temp,pressure);	/* gas mean free path (m) */
	
	return 2*lambda/dp;
}

/* particle slip correction coefficent (Allen & Raabe 1985, Aerosol Sci. Technol. 4, 269-286) */
real slip_correction_coefficient(real temp, real dp, real pressure) {
	real lambda; 	/* gas mean free path (m) */
	real lambdaDp;	/* lambda/dp */
	
	if (slipCorrectionLaw == 0) { /* This comes from the settings in CFD-TUTMAM GUI. */
		/* constant slip correction coefficient */
		return constantSlipCorrectionFactor;
	}
	
	lambda = gas_mean_free_path(temp,pressure);
	lambdaDp = lambda/dp;
	
	return 1+lambdaDp*(2.34+1.05*exp(-0.39/lambdaDp));
}

/* limits value to lowerLimit */
real tutmam_lower_limit(real value, real lowerLimit) {
	if (value < lowerLimit) {
		value = lowerLimit;
	}
	return value;
}

/* limits value to upperLimit */
real tutmam_upper_limit(real value, real upperLimit) {
	if (value > upperLimit) {
		value = upperLimit;
	}
	return value;
}

/* limits value to lowerLimit and upperLimit */
real tutmam_limits(real lowerLimit, real value, real upperLimit) {
	if (value < lowerLimit) {
		value = lowerLimit;
	}
	if (value > upperLimit) {
		value = upperLimit;
	}
	return value;
}

/* returns only a positive value */
real tutmam_positive(real value) {
	return tutmam_lower_limit(value,0.);
}

/* returns only a negative value */
real tutmam_negative(real value) {
	return tutmam_upper_limit(value,0.);
}

/* abscissas for Gauss-Hermite quadrature numerical integration (efunda.com) */
real gauss_hermite_abscissas(int i) {
	real abscissa = 0.;
	
	switch(gaussHermiteLevel) {
		case 5 :
			switch(i) {
				case 0 :
				case 4 :
					abscissa = -2.02018287046;
					break;
					
				case 1 :
				case 3 :
					abscissa = -0.958572464614;
					break;
					
				case 2 :
					abscissa = 0.;
					break;
				
				default :
					abscissa = 0.;
			}
			
			if (i > 2) {
				abscissa = 0.-abscissa;
			}
			
			break;
			
		case 7 : 
			switch(i) {
				case 0 :
				case 6 :
					abscissa = -2.65196135684;
					break;
					
				case 1 :
				case 5 :
					abscissa = -1.67355162877;
					break;
					
				case 2 :
				case 4 :
					abscissa = -0.816287882859;
					break;
				
				case 3 :
					abscissa = 0.;
					break;
					
				default :
					abscissa = 0.;
			}
			
			if (i > 3) {
				abscissa = 0.-abscissa;
			}
			
			break;

		case 9 :
			switch(i) {
				case 0 :
				case 8 :
					abscissa = -3.19099320178;
					break;
					
				case 1 :
				case 7 :
					abscissa = -2.26658058453;
					break;
					
				case 2 :
				case 6 :
					abscissa = -1.46855328922;
					break;
				
				case 3 :
				case 5 :
					abscissa = -0.723551018753;
					break;
					
				case 4 :
					abscissa = 0.;
					break;
					
				default :
					abscissa = 0.;
			}
			
			if (i > 4) {
				abscissa = 0.-abscissa;
			}
			
			break;

		case 11 :
			switch(i) {
				case 0 :
				case 10 :
					abscissa = -3.66847084656;
					break;
					
				case 1 :
				case 9 :
					abscissa = -2.78329009978;
					break;
					
				case 2 :
				case 8 :
					abscissa = -2.02594801583;
					break;
				
				case 3 :
				case 7 :
					abscissa = -1.32655708449;
					break;
					
				case 4 :
				case 6 :
					abscissa = -0.656809566882;
					break;
					
				case 5 :
					abscissa = 0.;
					break;
					
				default :
					abscissa = 0.;
			}
			
			if (i > 5) {
				abscissa = 0.-abscissa;
			}
			
			break;

		default :
			Error("Illegal gaussHermiteLevel\n");
	}

	return abscissa;
}

/* weights for Gauss-Hermite quadrature numerical integration (efunda.com) */
real gauss_hermite_weights(int i) {
	real weight = 0.;
	
	switch(gaussHermiteLevel) {
		case 5 :
			switch(i) {
				case 0 :
				case 4 :
					weight = 0.019953242059;
					break;
					
				case 1 :
				case 3 :
					weight = 0.393619323152;
					break;
					
				case 2 :
					weight = 0.945308720483;
					break;
				
				default :
					weight = 1.0;
			}		
			break;
			
		case 7 : 
			switch(i) {
				case 0 :
				case 6 :
					weight = 0.0009717812451;
					break;
					
				case 1 :
				case 5 :
					weight = 0.0545155828191;
					break;
					
				case 2 :
				case 4 :
					weight = 0.42560725261;
					break;
				
				case 3 :
					weight = 0.810264617557;
					break;
					
				default :
					weight = 1.0;
			}
			break;

		case 9 :
			switch(i) {
				case 0 :
				case 8 :
					weight = 3.96069772633e-5;
					break;
					
				case 1 :
				case 7 :
					weight = 0.00494362427554;
					break;
					
				case 2 :
				case 6 :
					weight = 0.0884745273944;
					break;
				
				case 3 :
				case 5 :
					weight = 0.432651559003;
					break;
					
				case 4 :
					weight = 0.720235215606;
					break;
					
				default :
					weight = 1.0;
			}
			break;

		case 11 :
			switch(i) {
				case 0 :
				case 10 :
					weight = 1.43956039371e-6;
					break;
					
				case 1 :
				case 9 :
					weight = 0.000346819466323;
					break;
					
				case 2 :
				case 8 :
					weight = 0.0119113954449;
					break;
				
				case 3 :
				case 7 :
					weight = 0.117227875168;
					break;
					
				case 4 :
				case 6 :
					weight = 0.429359752356;
					break;
					
				case 5 :
					weight = 0.654759286915;
					break;
					
				default :
					weight = 1.0;
			}
			break;

		default :
			Error("Illegal gaussHermiteLevel\n");
	}

	return weight;
}

/* abscissas for Gauss-Olin quadrature numerical integration */
real gauss_olin_abscissas_dp(int i, real D2) {

	switch(i) {
		case 0 :
			return powerLawD1;
						
		case 1 :
			return pow(powerLawD1,TUTMAM_23)*pow(D2,TUTMAM_13);
			
		case 2 :
			return pow(powerLawD1,TUTMAM_13)*pow(D2,TUTMAM_23);
		
		case 3 :
			return D2;
			
		default :
			return powerLawD1;
	}
	
	return powerLawD1;
}

/* weights for Gauss-Olin quadrature numerical integration */
real gauss_olin_weights(int i, real cc) {
	real ecc = exp(cc);
	
	if (fabs(cc) < 0.001) {
		switch(i) {
			case 0 :
			case 3 :
				return 0.125;
				
			case 1 :
			case 2 :
				return 0.375;
								
			default :
				return 0.0;
		}
	}
	
	switch(i) {
		case 0 :
			return 27*(ecc-1)/pow(cc,4) + (-9*ecc-18)/CBC(cc) + (ecc-5.5)/SQR(cc) - 1/cc;
			
		case 1 :
			return -81*(ecc-1)/pow(cc,4) + (36*ecc+45)/CBC(cc) + (-4.5*ecc+9)/SQR(cc);
			
		case 2 :
			return 81*(ecc-1)/pow(cc,4) - (45*ecc+36)/CBC(cc) + (9*ecc-4.5)/SQR(cc);
		
		case 3 :
			return -27*(ecc-1)/pow(cc,4) + (18*ecc+9)/CBC(cc) + (-5.5*ecc+1)/SQR(cc) + ecc/cc;
			
		default :
			return 0.0;
	}		
			
	return 0.0;
}

/* abscissas for Gauss-Olin quadrature numerical integration, D1 replaced with dCoag */
real gauss_olin_abscissas_dp_dCoag(int i, real dCoag, real D2) {

	switch(i) {
		case 0 :
			return dCoag;
						
		case 1 :
			return pow(dCoag,TUTMAM_23)*pow(D2,TUTMAM_13);
			
		case 2 :
			return pow(dCoag,TUTMAM_13)*pow(D2,TUTMAM_23);
		
		case 3 :
			return D2;
			
		default :
			return dCoag;
	}
	
	return dCoag;
}

/* factor for Gauss-Olin quadrature */
real quadrature_factor_olin(real cc) {
	if (fabs(cc) < 0.001) {
		return 1.0;
		
	} else {
		return cc/(exp(cc) - 1);
	}
}

/* density function for power law distribution (dimensionless) */
real powerlaw_density(real alpha, real d2, real dp) {
	real d = d2/powerLawD1;
	
	if (d-1.0 < 1.0e-3 || dp < powerLawD1 || dp > d2) {
		return 0.0;
	}
	
	if (fabs(alpha) < whenAlphaIsZero) {
		return 1.0/log(d);
	}
	
	return alpha/(pow(d2/dp,alpha)-pow(powerLawD1/dp,alpha));
}

/* k-th order of the moment is defined by UDS id */
real k_moment(int iUds) {
	switch(iUds % nUdsPerMode) {
		case iUdsM0 :
			/* k for number concentration */
			return 0.0;
			
		case iUdsM23 :
			/* k for surface area concentration is 2/3 */
			return TUTMAM_23;
			
		default :
			/* k for all mass concentrations */
			return 1.0;
	}
}

/* Calculates ln^2(GSD) from moment concentrations. */
real calculateLn2s(real numberConc,real surfaceConc,real massConc) {
	real ln2s = SQR(log(minGSD));
	
	/* If concentrations are zero, minGSD is used as GSD to prevent division by zero. */
	if (numberConc < minNumberConc) { /* (#/m^3) */
		return ln2s;
	}
	
	if (surfaceConc < minSurfaceConc) { /* (kg^(2/3) / m^3) */
		return ln2s;
	}
	
	if (massConc < minMassConc) { /* (kg/m^3) */
		return ln2s;
	}
	
	ln2s = TUTMAM_23*log(massConc) + TUTMAM_13*log(numberConc) - log(surfaceConc);
	
	/* GSD result limited between minGSD and maxGSD */
	return tutmam_limits(SQR(log(minGSD)),ln2s,SQR(log(maxGSD)));
}

/* Calculates count median diameter (m) from moment concentrations. */
real calculateCmd(real numberConc,real surfaceConc,real massConc,real rhoParticle) {
	real cmd = minCMD;
	
	/* If concentrations are zero, minCMD to prevent division by zero. */
	if (numberConc < minNumberConc) { /* (#/m^3) */
		return cmd;
	}
	
	if (surfaceConc < minSurfaceConc) { /* (kg^(2/3) / m^3) */
		return cmd;
	}
	
	if (massConc < minMassConc) { /* (kg/m^3) */
		return cmd;
	}
	
	cmd = TUTMAM_PI6M13*pow(rhoParticle,-TUTMAM_13)*pow(surfaceConc,1.5)*pow(massConc,-TUTMAM_23)*pow(numberConc,-5./6.);
	
	/* result limited between minCMD and maxCMD */
	return tutmam_limits(minCMD,cmd,maxCMD);
}

/* Calculates A and B of the PL distribution from moment concentrations */
void calculateAAndBMoments(real *AAndB, real numberConc, real surfaceConc, real massConc, real rhoParticle) {
	/* m1^(-1/3) */
	real m1m13 = TUTMAM_PI6M13*pow(rhoParticle,-TUTMAM_13)/powerLawD1;
	
	AAndB[0] = pow(massConc/numberConc,TUTMAM_13)*m1m13;
	AAndB[1] = massConc/surfaceConc*m1m13;
	
	return;
}

/* Calculates A and B of the PL distribution from distribution parameters */
void calculateAAndBParameters(real *AAndB, real *alphaAndD) {
	real a = alphaAndD[0];
	real d = alphaAndD[1];
	real a2;
	real a3;
	real da3;
	
	/* to prevent divisions by zero */
	if (fabs(a) < whenAlphaIsZero || fabs(a+2.0) < whenAlphaIsZero || fabs(a+3.0) < whenAlphaIsZero) {
		a += whenAlphaIsZero;
	}
	
	if (fabs(d-1.0) < whenAlphaIsZero) {
		d += whenAlphaIsZero;
	}
	
	a2 = a+2.0;
	a3 = a+3.0;
	da3 = pow(d,a3);
	
	AAndB[0] = pow(a/a3*(1.0-da3)/(1.0-pow(d,a)),TUTMAM_13);
	AAndB[1] = a2/a3*(1.0-da3)/(1.0-pow(d,a2));
	
	return;
}

/* Calculates JAcobian matrix of A and B of the PL distribution from distribution parameters */
void calculateJAAndB(real *JAAndB, real *alphaAndD) {
	real a = alphaAndD[0];
	real d = alphaAndD[1];
	real a2,a3,da,da1,da2,da3,lnd,dam1,da2m1,da3m1;
	
	/* to prevent divisions by zero */
	if (fabs(a) < whenAlphaIsZero || fabs(a+2.0) < whenAlphaIsZero || fabs(a+3.0) < whenAlphaIsZero) {
		a += whenAlphaIsZero;
	}
	
	if (fabs(d-1.0) < whenAlphaIsZero) {
		d += whenAlphaIsZero;
	}
	
	a2 = a+2.0;
	a3 = a+3.0;
	da = pow(d,a);
	da1 = pow(d,a+1.0);
	da2 = pow(d,a2);
	da3 = pow(d,a3);
	lnd = log(d);
	dam1 = da-1.0;
	da2m1 = da2-1.0;
	da3m1 = da3-1.0;
	
	JAAndB[0] = (((da3*(1.0+a*lnd)-1.0)*(da*a3-a3)-(da*(1.0+a3*lnd)-1.0)*(a*da3m1)))/SQR(da*a3-a3)/(3.0*pow(a/a3*da3m1/dam1,TUTMAM_23));
	JAAndB[1] = (a/a3*(a3*da2*dam1-a*pow(d,a-1.0)*da3m1))/SQR(dam1)/(3.0*pow(a/a3*da3m1/dam1,TUTMAM_23));
	JAAndB[2] = ((da3*(1.0+a2*lnd)-1.0)*(da2*a3-a3)-(da2*(1.0+a3*lnd)-1.0)*a2*da3m1)/SQR(da2*a3-a3);
	JAAndB[3] = a2/a3*(a3*da2*da2m1-a2*da1*da3m1)/SQR(da2m1);
	
	return;
}

/* Calculates alpha from moment concentrations. */
void calculateAlphaAndD2(real *alphaAndD2, real numberConc, real surfaceConc, real massConc, real rhoParticle, real *alphaAndDGuess) {
	real AAndB[2];
	real m1m13;
	real alphaAndD[2];
	
	if (numberConc < minNumberConc) { /* (#/m^3) */
		alphaAndD2[0] = 1.0;
		alphaAndD2[1] = powerLawD1;
		return;
	}
	
	if (surfaceConc < minSurfaceConc) { /* (kg^(2/3) / m^3) */
		alphaAndD2[0] = 1.0;
		alphaAndD2[1] = powerLawD1;
		return;
	}
	
	if (massConc < minMassConc) { /* (kg/m^3) */
		alphaAndD2[0] = 1.0;
		alphaAndD2[1] = powerLawD1;
		return;
	}
	
	calculateAAndBMoments(AAndB,numberConc,surfaceConc,massConc,rhoParticle);
	
	if (powerLawParametersSolvingMethod == 1) {
		interpolateAlphaAndD(alphaAndD,AAndB[0],AAndB[1]);
		
	} else {
		iterateAlphaAndD(alphaAndD,AAndB,alphaAndDGuess);
	}
		
	alphaAndD2[0] = tutmam_limits(minAlpha,alphaAndD[0],maxAlpha);
	alphaAndD2[1] = tutmam_limits(powerLawD1,alphaAndD[1]*powerLawD1,maxD2);
	return;
}

/* Calculates CMD of power law distribution from D1,D2, and alpha */
real calculateCmdOfPowerLaw(real D2, real alpha) {
	
	if (fabs(alpha) < whenAlphaIsZero) {
		return sqrt(powerLawD1*D2);
		
	} else {
		return pow(pow(powerLawD1,alpha)+pow(D2,alpha),1.0/alpha)/pow(2.0,1.0/alpha);
	}
}

/* calculates a multiplier for moment averaged D_p for power law distribution */
real momentAveragingPowerLaw(real alpha, real d, real k, real m) {
	real alpha3k = 3*k+alpha;
	
	if (d < 1.001 || fabs(m) < 0.001) {
		return 1.0;
	} 

	if (fabs(alpha3k+m) < whenAlphaIsZero) {
		return m*log(d)/(1.0-pow(d,0.0-m));
	}
	
	if (fabs(alpha3k) < whenAlphaIsZero) {
		return (pow(d,m)-1.0)/(m*log(d));
	}
	
	return alpha3k/(alpha3k+m)*(1.0-pow(d,alpha3k+m))/(1.0-pow(d,alpha3k));
}

/* Initialization function that sets initial values for UDSs and UDMs, when initialization is done in Fluent. */
void tutmam_initialize(Domain *mixtureDomain) {
 	int i;			/* integer for UDS loop */
	Thread *t;		/* thread variable */
    cell_t c;		/* cell variable */
	real rhoFluid;	/* fluid density */
	int j;			/* mode ID */
	
	/* initial vectors defined in tutmam_settings.h */
	real initialYSpeciesPVar[nTutmamModes][nTutmamSpecies] = { initialYSpeciesP };
	real initialNumberConcVar[nTutmamModes] = { initialNumberConc };
	real initialCMDVar[nTutmamModes] = { initialCMD };
	real initialGSDVar[nTutmamModes] = { initialGSD };
	real initialRhoParticleVar[nTutmamModes] = { initialRhoParticle };
	
	/* looping over all threads of the domain */
	thread_loop_c(t,mixtureDomain)
    {
		 if (THREAD_ID(t) == fluidZoneNumber) { /* if thread is the fluid zone */
		 
			/* looping over all cells of the fluid zone */
			begin_c_loop(c,t)
			
				rhoFluid = C_R(c,t);
				
				for (j = 0; j < nTutmamModes; ++j) { /* looping over all modes */
					/* UDS for number concentration is set by the value in tutmam_settings.h */
					C_UDS_M0(c,t,j) = initialNumberConcVar[j]/rhoFluid;
					
					/* UDS for surface area concentration is calculated from the values in tutmam_settings.h */
					C_UDS_M23(c,t,j) = pow(TUTMAM_PI6*initialRhoParticleVar[j],TUTMAM_23)*initialNumberConcVar[j]*SQR(initialCMDVar[j])*exp(2*SQR(log(initialGSDVar[j])))/rhoFluid;
					
					/* UDSs for species mass concentrations are calculated from the values in tutmam_settings.h */
					for (i = 0; i < nTutmamSpecies; ++i) {
						C_UDS_M1(c,t,i,j) = initialYSpeciesPVar[j][i]*TUTMAM_PI6*initialRhoParticleVar[j]*initialNumberConcVar[j]*pow(initialCMDVar[j],3.)*exp(4.5*pow(log(initialGSDVar[j]),2.0))/rhoFluid;
					}
					
					/* UDMs are set by the values in tutmam_settings.h */
					C_NTOT(c,t,j) = initialNumberConcVar[j];
					
					if (j == powerLawDistribution) {
						C_LN2S(c,t,j) = 1.0;
						C_CMD(c,t,j) = powerLawD1;
						
					} else {
						C_LN2S(c,t,j) = SQR(log(initialGSDVar[j]));
						C_CMD(c,t,j) = initialCMDVar[j];
					}
					
					C_R_PARTICLE(c,t,j) = initialRhoParticleVar[j];
					
					/* particle temperature is set to fluid temperature */
					C_T_PARTICLE(c,t,j) = C_T(c,t);
					
					/* water equilibrium factors are set to 1 */
					C_WATER_KAPPA(c,t,j) = 1.0;
					
					/* Equilibrium rh / Fluid rh is set to -1 that denotes no value */
					C_UDMI(c,t,iRhEqPerRhUdm+j*nUdmPerMode) = -1.0;
					
					/* All source terms are set to 0 */
					C_UDMI(c,t,iUdmCondensationM23+j*nUdmPerMode) = 0.0;
					for (i = 0; i < nTutmamSpecies; ++i) {
						C_UDMI(c,t,iUdmCondensationM1First+i+j*nUdmPerMode) = 0.0;
					}
				}
				
				/* Nucleation source terms are set to 0 */
				C_UDMI(c,t,iUdmNucleationM0) = 0.0;
				C_UDMI(c,t,iUdmNucleationM23) = 0.0;
				for (i = 0; i < nTutmamSpecies; ++i) {
					C_UDMI(c,t,iUdmNucleationM1First+i) = 0.0;
				}
				
				#if powerLawDistribution > -1
					#if nTutmamModes > 1
						/* PL to LN transfer source terms are set to 0 */
						C_UDMI(c,t,iUdmTransferPl2LnM0) = 0.0;
						C_UDMI(c,t,iUdmTransferPl2LnM23) = 0.0;
						C_UDMI(c,t,iUdmTransferPl2LnM1) = 0.0;
					#endif
				#endif
				
			end_c_loop(c,t)
		}
	}
	
	return;
}

/* The following functions are settings transfer functions that sets the values of */
/* C-code varibles to the values of RP-variables. */
/* In parallel solver, C->RP variable transfer is done in host process, */
/* but C->C variable transfer is done from host to node in node processes. */

/* General settings -settings transfer function */
void tutmam_transfer_settings_general_settings() {

#if RP_HOST || !PARALLEL /* serial solver or host process of parallel solver */
	
	/* setting the values of C-code varibles to the values of RP-variables */	
	fluidZoneNumber = RP_Get_Integer("fluid_zone_number_rp");
	gaussHermiteLevel = RP_Get_Integer("gauss_hermite_level_rp");
	condensationIntegrationLevel = RP_Get_Integer("condensation_integration_level_rp");
	condensationIntegrationBins = RP_Get_Integer("condensation_integration_bins_rp");
	coagulationIntegrationLevel = RP_Get_Integer("coagulation_integration_level_rp");
	coagulationIntegrationBins = RP_Get_Integer("coagulation_integration_bins_rp");
	
	Message("\n");
	Message("Fluid zone ID:                %d\n",fluidZoneNumber);
	Message("Gauss-Hermite quad. level:    %d\n",gaussHermiteLevel);
	
	#if powerLawDistribution > -1
		Message("Integration for condensation:\n    ");
		switch(condensationIntegrationLevel) {
			case 1:
				Message("Gauss-Olin");
				break;
				
			case 2:
				Message("Gauss-Olin when alpha > 0.5\n    Numeric (%d bins) when alpha < 0.5",condensationIntegrationBins);
				break;
				
			case 3:
				Message("Numeric (%d bins)",condensationIntegrationBins);
				break;
				
			default:
				Message("Unknown condensationIntegrationLevel");
		}
		Message("\n");
		
		Message("Integration for coagulation:\n    ");
		switch(coagulationIntegrationLevel) {
			case 1:
				Message("Gauss-Olin");
				break;
				
			case 2:
				Message("Gauss-Olin when D2/D1 < 3\n    Numeric (%d bins) when D2/D1 > 3",coagulationIntegrationBins);
				break;
				
			case 3:
				Message("Numeric (%d bins)",coagulationIntegrationBins);
				break;
				
			default:
				Message("Unknown coagulationIntegrationLevel");
		}
		Message("\n");
	#endif
			
	
#endif

#if PARALLEL /* parallel solver */

	/* setting the C-code variables from host process to node processes */
	host_to_node_int_6(fluidZoneNumber,gaussHermiteLevel,condensationIntegrationLevel,condensationIntegrationBins,coagulationIntegrationLevel,coagulationIntegrationBins);
	
#endif

	return;
}

/* Diffusion -settings transfer function */
void tutmam_transfer_settings_diffusion() {

#if RP_HOST || !PARALLEL /* serial solver or host process of parallel solver */
	
	Message("\n");
	
	/* setting the values of C-code varibles to the values of RP-variables */	
	diffusionLaw = RP_Get_Boolean("diffusion_law_rp");
	slipCorrectionLaw = RP_Get_Boolean("slip_correction_law_rp");
	
	particleTurbDiff = RP_Get_Boolean("particle_turb_diff_rp");
	particleTurbSchmidtNumber = RP_Get_Real("particle_turb_schmidt_number_rp");
	if (particleTurbSchmidtNumber < 1.e-9) { /* preventing division by zero */
		Message("Particle turbulent Schmidt number must be higher than 0. Value set to 0.7\n");
		particleTurbSchmidtNumber = 0.7;
	}
	
	constantSlipCorrectionFactor = RP_Get_Real("constant_slip_correction_factor_rp");
	if (constantSlipCorrectionFactor < 1.0) {
		Message("Slip correction coefficent must be higher than 1. Value set to 1\n");
		constantSlipCorrectionFactor = 1.0;
	}
	
	Message("Calculate slip corr. coeff.:        %d\n",slipCorrectionLaw);
	if (slipCorrectionLaw) {
		Message("Diffusion integral parametrization: %d\n",diffusionLaw);
	} else {
		Message("Constant slip corr. coeff.:         %g\n",constantSlipCorrectionFactor);	
	}
	
	Message("Particle turbulent diffusion:       %d\n",particleTurbDiff);
	if (particleTurbDiff) {
		Message("Particle turbulent Schmidt number:  %g\n",particleTurbSchmidtNumber);
	}
	
#endif

#if PARALLEL /* parallel solver */

	/* setting the C-code variables from host process to node processes */
	host_to_node_int_3(diffusionLaw,slipCorrectionLaw,particleTurbDiff);
	host_to_node_real_2(particleTurbSchmidtNumber,constantSlipCorrectionFactor);
	
#endif

	return;
}

/* Aerosol distribution -settings transfer function */
void tutmam_transfer_settings_aerosol_distribution() {

#if RP_HOST || !PARALLEL /* serial solver or host process of parallel solver */
	
	Message("\n");
	
	/* setting the values of C-code varibles to the values of RP-variables */	
	minGSD = RP_Get_Real("min_gsd_rp");
	/* Preventing division by zero. If GSD is equal to 1, ln(GSD) is zero and division by zero occurs. */
	if (minGSD < 1.0000001) {
		Message("Minimum GSD must be higher than 1. Value set to 1.01\n");
		minGSD = 1.01;
	}
	
	maxGSD = RP_Get_Real("max_gsd_rp");
	/* GSD range cannot be zero */
	if (maxGSD-minGSD < 1.0e-6) {
		Message("Maximum GSD must be higher than minimum GSD. Value set to 3\n");
		maxGSD = 3.0;
	}
	
	powerLawParametersSolvingMethod = RP_Get_Boolean("power_law_parameters_solving_method_rp");
	
	Message("Minimum GSD: %g\n",minGSD);
	Message("Maximum GSD: %g\n",maxGSD);
	
	Message("PL distr. param. solving method:\n    ");
	switch(powerLawParametersSolvingMethod) {
		case 0:
			Message("Levenberg-Marquardt iteration algorithm\n");
			break;
			
		case 1:
			Message("Interpolation table\n");
			break;
			
		default:
			Error("Unknown PL distr. param. solving method\n");
	}
	
#endif

#if PARALLEL /* parallel solver */

	/* setting the C-code variables from host process to node processes */
	host_to_node_real_2(minGSD,maxGSD);
	host_to_node_int_1(powerLawParametersSolvingMethod);
	
#endif

	return;
}

/* Nucleation -settings transfer function */
void tutmam_transfer_settings_nucleation() {

#if RP_HOST || !PARALLEL /* serial solver or host process of parallel solver */

	char *nucleatingSpeciesString;
	char *nucleationExponentsString;
	char *nucleationSatVapPresExponentsString;
	char *nMolecClusterVectorString;
			
	Message("\n");
	
	/* setting the values of C-code varibles to the values of RP-variables */	
	uRFNucleation = RP_Get_Real("urf_nucleation_rp");
	if (uRFNucleation < 0.0 || uRFNucleation > 1.0) {
		Message("Under-relaxation factors must be in the range between 0 and 1. Value set to 1\n");
		uRFNucleation = 1.0;
	}

	nucleationLaw = RP_Get_Integer("nucleation_law_rp");
	nucleationCorrectionFactor = RP_Get_Real("nucleation_correction_factor_rp");
	
	nucleatingSpeciesString=RP_Get_String("nucleating_species_rp");
	nucleationExponentsString=RP_Get_String("nucleation_exponents_rp");
	nucleationSatVapPresExponentsString=RP_Get_String("nucleation_sat_vap_pres_exponents_rp");
	nMolecClusterVectorString=RP_Get_String("n_molec_cluster_vector_rp");
	string_to_int_vector(nucleatingSpeciesString,iNucleatingSpecies,nTutmamSpecies);
	string_to_vector(nucleationExponentsString,nucleationExponents);
	string_to_vector(nucleationSatVapPresExponentsString,nucleationSatVapPresExponents);
	string_to_vector(nMolecClusterVectorString,nMolecClusterVector);
	
	Message("URF for nucleation:                          %g\n",uRFNucleation);
	Message("Nucleation law ID number:                    %d\n",nucleationLaw);
	Message("Nucleation correction factor:                %g\n",nucleationCorrectionFactor);
	
	if (nucleationLaw == 2) {
		message_Olin_nucleation_parameters_vector(iNucleatingSpecies,nucleationExponents,nucleationSatVapPresExponents,nMolecClusterVector);
	}
	
#endif

#if PARALLEL /* parallel solver */

	/* setting the C-code variables from host process to node processes */
	host_to_node_real_2(uRFNucleation,nucleationCorrectionFactor);
	host_to_node_int_1(nucleationLaw);
	host_to_node_int(iNucleatingSpecies,nTutmamSpecies);
	host_to_node_real(nucleationExponents,nTutmamSpecies);
	host_to_node_real(nucleationSatVapPresExponents,nTutmamSpecies);
	host_to_node_real(nMolecClusterVector,nTutmamSpecies);
	
#endif

	return;
}

/* Condensation -settings transfer function */
void tutmam_transfer_settings_condensation() {

#if RP_HOST || !PARALLEL /* serial solver or host process of parallel solver */
	
	char *uRFCondensationM1String;
	char *condensingSpeciesString;
	char *condensationDirectionVectorString;
	
	Message("\n");
	
	/* setting the values of C-code varibles to the values of RP-variables */	
	uRFCondensationM23 = RP_Get_Real("urf_condensation_m23_rp");
	if (uRFCondensationM23 < 0.0 || uRFCondensationM23 > 1.0) {
		Message("Under-relaxation factors must be in the range between 0 and 1. Value set to 1\n");
		uRFCondensationM23 = 1.0;
	}
	
	uRFCondensationM1String = RP_Get_String("urf_condensation_m1_rp");
	string_to_vector(uRFCondensationM1String,uRFCondensationM1);
	latentHeatOfCondensation = RP_Get_Boolean("latent_heat_of_condensation_rp");
	condensingSpeciesString = RP_Get_String("condensing_species_rp");
	kelvinEffect = RP_Get_Boolean("kelvin_effect_rp");
	condensationDirectionVectorString = RP_Get_String("condensation_direction_vector_rp");
	phaseActivityModel = RP_Get_Boolean("phase_activity_model_rp");
	waterEq = RP_Get_Boolean("water_eq_rp");
	waterEqSpecies = RP_Get_Integer("water_eq_species_rp");
	minKappa = RP_Get_Real("min_kappa_rp");
	maxKappa = RP_Get_Real("max_kappa_rp");
	dieselExhaustHCFractionModel = RP_Get_Boolean("diesel_exhaust_hc_fraction_model_rp");
	dieselExhaustHC = RP_Get_Integer("diesel_exhaust_hc_rp");
	interModalCondensationFactor = RP_Get_Real("inter_modal_condensation_factor_rp");
	
 	string_to_int_vector(condensingSpeciesString,iCondensingSpecies,nTutmamModes*nTutmamSpecies);
	string_to_int_vector(condensationDirectionVectorString,condensationDirectionVector,nTutmamSpecies);
	message_iCondensingSpecies_vector(iCondensingSpecies);
	message_condensationDirectionVector_vector(condensationDirectionVector);
	
	Message("URF for surface area moment condensation:    %g\n",uRFCondensationM23);
	message_URF_vector(uRFCondensationM1);
	Message("Latent heat effect included in condensation: %d\n",latentHeatOfCondensation);
	Message("Kelvin effect:                               %d\n",kelvinEffect);
	
	if (nTutmamPhases > 1) {
		Message("Phase activity model:                        %d\n",phaseActivityModel);
	}
	
	if (waterEq == 0) {
		Message("Water condensation calculated normally\n");
		
	} else {
		Message("Water condensation calculated through equilibrium\n      Condensation term is connected to species %d\n",waterEqSpecies);
		Message("      Kappa range: %g - %g\n",minKappa,maxKappa);
	}

	Message("Diesel exhaust HC fraction model:            %d\n",dieselExhaustHCFractionModel);
	if (dieselExhaustHCFractionModel == 1) {
		Message("HC species ID:                               %d\n",dieselExhaustHC);
	}
	
	#if powerLawDistribution > -1 && nTutmamModes > 1
		Message("Intermodal condensation factor:              %g\n",interModalCondensationFactor);
	#endif
	
#endif

#if PARALLEL /* parallel solver */

	/* setting the C-code variables from host process to node processes */
	host_to_node_real_1(uRFCondensationM23);
	host_to_node_real(uRFCondensationM1,nTutmamSpecies);
	host_to_node_int_3(latentHeatOfCondensation,kelvinEffect,phaseActivityModel);
	host_to_node_int(iCondensingSpecies,nTutmamSpecies*nTutmamModes);
	host_to_node_int(condensationDirectionVector,nTutmamSpecies);
	host_to_node_int_2(waterEq,waterEqSpecies);
	host_to_node_real_3(minKappa,maxKappa,interModalCondensationFactor);
	host_to_node_int_2(dieselExhaustHCFractionModel,dieselExhaustHC);
	
#endif

	return;
}

/* Aerosol process control -settings transfer function */
void tutmam_transfer_settings_aerosol_process_control() {

#if RP_HOST || !PARALLEL /* serial solver or host process of parallel solver */
	
	Message("\n");
	
	/* setting the values of C-code varibles to the values of RP-variables */	
	diffusionProcess = RP_Get_Boolean("diffusion_process_rp");
	nucleationProcess = RP_Get_Boolean("nucleation_process_rp");
	nucleationComputing = RP_Get_Boolean("nucleation_computing_rp");
	condensationProcess = RP_Get_Boolean("condensation_process_rp");
	condensationComputing = RP_Get_Boolean("condensation_computing_rp");
	coagulationProcess = RP_Get_Boolean("coagulation_process_rp");
	coagulationComputing = RP_Get_Boolean("coagulation_computing_rp");
	selfCoagulationalTransferProcess = RP_Get_Boolean("self_coagulational_transfer_process_rp");
	interModalCondensationProcess = RP_Get_Boolean("inter_modal_condensation_process_rp");
	
	Message("Diffusion    process:    %d\n",diffusionProcess);
	
	Message("Nucleation   process:    %d\n",nucleationProcess);
	
	if (nucleationProcess) {
	Message("             computing:  %d\n",nucleationComputing);
	}
	
	Message("Condensation process:    %d\n",condensationProcess);
	
	if (condensationProcess) {
	Message("             computing:  %d\n",condensationComputing);
	}
	
	Message("Coagulation  process:    %d\n",coagulationProcess);
	
	if (coagulationProcess) {
	Message("             computing:  %d\n",coagulationComputing);
	}

	#if powerLawDistribution > -1
		#if nTutmamModes > 1
			Message("Transfer processes from PL to LN:\n");
			Message("  self-coagulation:        %d\n",selfCoagulationalTransferProcess);
			Message("  intermodal condensation: %d\n",interModalCondensationProcess);
		#endif
	#endif
	
#endif

#if PARALLEL /* parallel solver */

	/* setting the C-code variables from host process to node processes */
	host_to_node_int_7(diffusionProcess,nucleationProcess,nucleationComputing,condensationProcess,condensationComputing,coagulationProcess,coagulationComputing);
	host_to_node_int_2(selfCoagulationalTransferProcess,interModalCondensationProcess);
	
#endif

	return;
}

/* Coagulation -settings transfer function */
void tutmam_transfer_settings_coagulation() {

#if RP_HOST || !PARALLEL /* serial solver or host process of parallel solver */
	int j1,j2;
	
	Message("\n");
	
	/* setting the values of C-code varibles to the values of RP-variables */	
	intraModalCoagulationProcess = RP_Get_Boolean("intra_modal_coagulation_process_rp");
	interModalCoagulationProcess = RP_Get_Boolean("inter_modal_coagulation_process_rp");
	coagulationIterationSkip = RP_Get_Integer("coagulation_iteration_skip_rp");
	
	Message("Intramodal coagulation:                         %d\n",intraModalCoagulationProcess);
	Message("  for modes ");
	for (j1 = 0; j1 < nTutmamModes; ++j1) {
		if (coagulationMatrix[j1][j1] == 1) {
			Message("%d ",j1);
		}
	}
	Message("\n");
	
	
	
	if (coagulationIterationSkip < 1) {
		Message("Iteration skip must be 1 or more. 1 represents every iteration. Value set to 1.\n");
		coagulationIterationSkip = 1;
	}
	Message("Intramodal coagulation iteration skip:          %d\n",coagulationIterationSkip);
	
	uRFCoagulation = RP_Get_Real("urf_coagulation_rp");
	if (uRFCoagulation < 0.0 || uRFCoagulation > 1.0) {
		Message("Under-relaxation factors must be in the range between 0 and 1. Value set to 0.5\n");
		uRFCoagulation = 0.5;
	}
	Message("URF for intramodal coagulation:                 %g\n",uRFCoagulation);
	
	coagulationRobustnessModel = RP_Get_Boolean("coagulation_robustness_model_rp");
	Message("Robustness model:                               %d\n",coagulationRobustnessModel);
	
	Message("Intermodal coagulation:                         %d\n",interModalCoagulationProcess);
	for (j1 = 0; j1 < nTutmamModes; ++j1) {
		Message("  for mode %d:\n",j1);
		Message("    losses to modes ");
		for (j2 = 0; j2 < nTutmamModes; ++j2) {
			if (j1 != j2) {
				if (coagulationMatrix[j1][j2] == 1) {
					Message("%d ",j2);
				}
			}
		}
		Message("\n");
		
		Message("    gains from modes ");
		for (j2 = 0; j2 < nTutmamModes; ++j2) {
			if (j1 > j2) {
				if (coagulationMatrix[j1][j2] == 1) {
					Message("%d ",j2);
				}
			}
		}
		Message("\n");
	}
	
	transitionRegimeCorrectionFactorForCoagulationLaw = RP_Get_Integer("transition_regime_correction_factor_for_coagulation_law_rp");
	Message("Transition regime correction law (coagulation): %d\n",transitionRegimeCorrectionFactorForCoagulationLaw);
	
#endif

#if PARALLEL /* parallel solver */

	/* setting the C-code variables from host process to node processes */
	host_to_node_real_1(uRFCoagulation);
	host_to_node_int_5(intraModalCoagulationProcess,interModalCoagulationProcess,coagulationIterationSkip,transitionRegimeCorrectionFactorForCoagulationLaw,coagulationRobustnessModel);
	
#endif

	return;
}

/* this converts a string to a real vector */
void string_to_vector(const char *string, real *vector) {
	char s[] = " ,;/";						/* string containing delimiters */
	char *token;							/* pointer to parts of the string */
	real value;								/* value of a part of the string */
	int i = 0;								/* index for vector */
	char changingString[nTutmamSpecies*10];	/* a string that is changed in splitting */

	/* copying string to changingString */
	strcpy(changingString,string);

	/* finding the first part of the string */
	token = strtok(changingString, s);

	/* finding other parts of the string until the string or the vector ends */
	while (token != NULL && i < nTutmamSpecies) {
		/* converting string to a real number */
		value = atof(token);
		vector[i] = value;
		++i;

		/* finding the next part of the string */
		token = strtok(NULL, s);
	}
	
	/* if the string is longer than the vector */
	if (token != NULL) {
		Message("There are some extra parts in the string compared to the vector length.\n");
	}
	
	/* If the string was not complete, the last values are set to 0. */
	while (i < nTutmamSpecies) {
		Message("Value for vector[%d] is not found in the string. Value set to 0.\n",i);
		vector[i] = 0.0;
		++i;
	}

	return;
}

/* this converts a string to an integer vector */
void string_to_int_vector(const char *string, int *vector, int length) {
	char s[] = " ,;/";					/* string containing delimiters */
	char *token;						/* pointer to parts of the string */
	int value;							/* value of a part of the string */
	int i = 0;							/* index for vector */
	char changingString[nTutmamModes*nTutmamSpecies*10]; /* a string that is changed in splitting */

	/* copying string to changingString */
	strcpy(changingString,string);

	/* finding the first part of the string */
	token = strtok(changingString, s);

	/* finding other parts of the string until the string or the vector ends */
	while (token != NULL && i < length) {
		/* converting string to a real number */
		value = atoi(token);
		vector[i] = value;
		++i;

		/* finding the next part of the string */
		token = strtok(NULL, s);
	}
	
	/* if the string is longer than the vector */
	if (token != NULL) {
		Message("There are some extra parts in the string compared to the vector length.\n");
	}
	
	/* If the string was not complete, the last values are set to 0. */
	while (i < length) {
		Message("Value for vector[%d] is not found in the string. Value set to 0.\n",i);
		vector[i] = 0;
		++i;
	}

	return;
}

void message_iCondensingSpecies_vector(int *iCondensingSpeciesVector) {
	int i;
	int j;
	char *speciesName;
	Domain *domain = NULL;
	Material *mixMat = NULL;
	domain = Get_Domain(1);
	mixMat = mixture_material(domain);
	
	Message("Condensing species:\n");
	
	for (j = 0; j < nTutmamModes; ++j) { /* looping over all modes */	
		Message("  Mode %d:  ",j);
		
		for (i = 0; i < nTutmamSpecies; ++i) {			/* looping over all species */		
			if (iCondensingSpeciesVector[i+j*nTutmamSpecies] == 1) {
				if (mixMat == NULL) {
					Message("%d  ",i);
				} else {
					speciesName = MIXTURE_SPECIE_NAME(mixMat,fluentSpeciesIdVector[i]);
					Message("%s  ",speciesName);
				}
			}
		}
		
		Message("\n");
	}
	
	Message("\n");
}

void message_condensationDirectionVector_vector(int *condensationDirectionVector) {
	int i;
	int iFluentSpeciesWater = -1;
	char *speciesName;
	Domain *domain = NULL;
	Material *mixMat = NULL;
	domain = Get_Domain(1);
	mixMat = mixture_material(domain);
	
	if (waterEq == 1 && mixMat != NULL) { /* if water-equilibrium calculation is on */

		/* find Fluent species ID for water */
		iFluentSpeciesWater = SV_SpeciesIndex("h2o");
		if (iFluentSpeciesWater == -1) {
			Error("h2o not found! Switch water equilibrium calculation off.\n");
			return;
		}
	}
	
	Message("Condensation directions:\n");
			
	for (i = 0; i < nTutmamSpecies; ++i) {			/* looping over all species */		
		if (waterEq == 0 || fluentSpeciesIdVector[i] != iFluentSpeciesWater) {
			if (condensationDirectionVector[i] == -1) {
				if (mixMat == NULL) {
					Message("%d: evaporation only\n",i);
				} else {
					speciesName = MIXTURE_SPECIE_NAME(mixMat,fluentSpeciesIdVector[i]);
					Message("%s: evaporation only\n",speciesName);
				}
				
			} else if (condensationDirectionVector[i] == 0) {
				if (mixMat == NULL) {
					Message("%d: condensation and evaporation\n",i);
				} else {
					speciesName = MIXTURE_SPECIE_NAME(mixMat,fluentSpeciesIdVector[i]);
					Message("%s: condensation and evaporation\n",speciesName);
				}
				
			} else {
				if (mixMat == NULL) {
					Message("%d: condensation only\n",i);
				} else {
					speciesName = MIXTURE_SPECIE_NAME(mixMat,fluentSpeciesIdVector[i]);
					Message("%s: condensation only\n",speciesName);
				}
			}
		}
	}
	
	if (waterEq == 1) {
		if (mixMat == NULL) {
			Message("h2o: connected to species %d ",waterEqSpecies);
		} else {
			speciesName = MIXTURE_SPECIE_NAME(mixMat,fluentSpeciesIdVector[waterEqSpecies]);
			Message("h2o: connected to %s ",speciesName);
		}
		
		if (condensationDirectionVector[waterEqSpecies] == -1) {
			Message("(evaporation only)\n");
		} else if (condensationDirectionVector[waterEqSpecies] == 0) {
			Message("(condensation and evaporation)\n");
		} else {
			Message("(condensation only)\n");
		}
	}

	Message("\n");
}

/* prints a message of URF vector values and checks that the values are between 0 and 1 */
void message_URF_vector(real *vector) {
	int i;							/* index in for-loop */
	int warningMessagePrinted = 0;	/* flag for preventing warning message to be printed more than once */
	
	Message("Under-relaxation factors for mass moment condensation\n");
	for (i = 0; i < nTutmamSpecies; ++i) {			/* looping over all species */
		if (vector[i] < 0.0 || vector[i] > 1.0) {	/* if illegal value */
			if (warningMessagePrinted == 0) {
				Message("URF must be between 0 and 1. Illegal values set to 1.\n");
				warningMessagePrinted = 1;
			}
			
			vector[i] = 1.0; /* setting illegals to 1 */
		}
		
		Message("    species %d: %g\n",i,vector[i]);
	}

	Message("\n");

	return;
}

/* prints a message of Olin nucleation law parameters */
void message_Olin_nucleation_parameters_vector(int *iNucleatingSpeciesVector, real *nucleationExponentsVector, real *nucleationSatVapPresExponentsVector, real *nMolecClusterVectorVector) {
	int i;							/* index in for-loop */
	char *speciesName;
	Domain *domain = NULL;
	Material *mixMat = NULL;
	domain = Get_Domain(1);
	mixMat = mixture_material(domain);

	Message("Olin nucleation law parameters\n");
	Message("    nucleating species:              ");
	for (i = 0; i < nTutmamSpecies; ++i) {			/* looping over all species */		
		if (iNucleatingSpeciesVector[i] == 1) {
			if (mixMat == NULL) {
				Message("%d  ",i);
			} else {
				speciesName = MIXTURE_SPECIE_NAME(mixMat,fluentSpeciesIdVector[i]);
				Message("%s  ",speciesName);
			}
		}
	}
	Message("\n");
	
	Message("    nucl. exponents:                 ");
	for (i = 0; i < nTutmamSpecies; ++i) {			/* looping over all species */		
		if (iNucleatingSpeciesVector[i] == 1) {
			Message("%g      ",nucleationExponentsVector[i]);
		}
	}
	Message("\n");
	
	Message("    nucl. sat. vap. pres. exponents: ");
	for (i = 0; i < nTutmamSpecies; ++i) {			/* looping over all species */		
		if (iNucleatingSpeciesVector[i] == 1) {
			Message("%g      ",nucleationSatVapPresExponentsVector[i]);
		}
	}
	Message("\n");
	
	Message("    molecules in critical cluster:   ");
	for (i = 0; i < nTutmamSpecies; ++i) {			/* looping over all species */		
		if (iNucleatingSpeciesVector[i] == 1) {
			Message("%g      ",nMolecClusterVectorVector[i]);
		}
	}
	Message("\n");

	Message("\n");

	return;
}

/* This function saves the particle distribution parameters to UDMs. */
/* This is called after each iteration. */
void tutmam_save_udm_variables() {
    Thread *t;				/* thread variable */
    cell_t c;				/* cell variable */
	Domain *mixtureDomain;	/* domain variable */
	real rhoFluid;			/* fluid density (kg/m^3) */
	real rhoParticle;		/* particle density (kg/m^3) */
	real numberConc;		/* number concentration (#/m^3) */
	real surfaceConc;		/* surface area concentration (kg^(2/3) / m^3) */
	real massConc;			/* sum of all mass concentrations (kg/m^3) */
	int iSpecies;			/* integer for particle species loop */
	int j;					/* mode ID */
	real alphaAndD2[2];
	real alphaAndDGuess[2];
	
	/* fluid domain has (always?) id 1 */
    mixtureDomain = Get_Domain(1);

	/* looping over all cell threads of the domain */
	thread_loop_c(t,mixtureDomain) {

		if (THREAD_ID(t) == fluidZoneNumber) { /* if the thread is the fluid zone */
		
			/* looping over all cells of the fluid zone */
			begin_c_loop(c,t)
			{
				rhoFluid = C_R(c,t);
				
				for (j = 0; j < nTutmamModes; ++j) { /* looping over all modes */
					/* negative UDS values ignored */
					numberConc = tutmam_positive(C_UDS_M0(c,t,j)*rhoFluid);  
					surfaceConc = tutmam_positive(C_UDS_M23(c,t,j)*rhoFluid);
					
					rhoParticle = C_R_PARTICLE(c,t,j);
					massConc = 0.; /* summation starts from zero */
					for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) { /* looping over all particle species */
						massConc += tutmam_positive(C_UDS_M1(c,t,iSpecies,j)); /* summing the mass concentrations */
						/* here, the unit is unitless */
					}
					massConc = massConc*rhoFluid; /* the unit is now kg/m^3 */
					
					/* storing the particle distribution parameters to UDMs */
					C_NTOT(c,t,j) = numberConc;
					
					if (j == powerLawDistribution) {
						alphaAndDGuess[0] = C_LN2S(c,t,j);
						alphaAndDGuess[1] = C_CMD(c,t,j)/powerLawD1;
						
						calculateAlphaAndD2(alphaAndD2,numberConc,surfaceConc,massConc,rhoParticle,alphaAndDGuess);
						C_LN2S(c,t,j) = uRFPowerLawParameters*alphaAndD2[0] + (1.0-uRFPowerLawParameters)*C_LN2S(c,t,j);
						C_CMD(c,t,j) = uRFPowerLawParameters*alphaAndD2[1] + (1.0-uRFPowerLawParameters)*C_CMD(c,t,j);
						
					} else {
						C_LN2S(c,t,j) = calculateLn2s(numberConc,surfaceConc,massConc);
						C_CMD(c,t,j) = calculateCmd(numberConc,surfaceConc,massConc,rhoParticle);
					}
					
					/* storing particle density to UDM */
					C_R_PARTICLE(c,t,j) = particle_density(c,t,j);
					
					if (condensationProcess == 1 && condensationComputing == 1 && latentHeatOfCondensation == 1) {
						/* storing particle temperature to UDM */
						C_T_PARTICLE(c,t,j) = particle_temperature(c,t,j);
					}

					if (condensationProcess == 1 && condensationComputing == 1 && waterEq == 1) {
						C_WATER_KAPPA(c,t,j) = water_kappa(c,t,j);
					}
				}				
			}				
			end_c_loop(c,t)
		} 
	}	
		
	return;
}

/* This shows what UDS and UDM slots are mapped. */
void show_uds_udm_mapping() {
	
	/* indices of UDSs */
	int iM0Uds = iUdsM0;
	int iM23Uds = iUdsM23;
	int iFirstM1Uds = iUdsM1First;
	int iLastM1Uds = iUdsM1First+nTutmamSpecies-1;
	
	/* indices of UDMs */
	int iNumberConcUdm = iUdmNumberConc;
	int iCmdUdm = iUdmCmd;
	int iLn2sUdm = iUdmLn2s;
	int iParticleDensityUdm = iUdmParticleDensity;
	int iParticleTemperatureUdm = iUdmParticleTemperature;
	int iNucleationM0Udm = iUdmNucleationM0;
	int iNucleationM23Udm = iUdmNucleationM23;
	int iNucleationM1FirstUdm = iUdmNucleationM1First;
	int iNucleationM1LastUdm = iUdmNucleationM1First+nTutmamSpecies-1;
	int iCondensationM23Udm = iUdmCondensationM23;
	int iCondensationM1FirstUdm = iUdmCondensationM1First;
	int iCondensationM1LastUdm = iUdmCondensationM1First+nTutmamSpecies-1;
	int iCoagulationM0Udm = iUdmCoagulationM0;
	int iCoagulationM23Udm = iUdmCoagulationM23;
	
	#if powerLawDistribution > -1
		#if nTutmamModes > 1
			int iTransferPl2LnM0Udm = iUdmTransferPl2LnM0;
			int iTransferPl2LnM23Udm = iUdmTransferPl2LnM23;
			int iTransferPl2LnM1Udm = iUdmTransferPl2LnM1;
		#endif
	#endif
	
	int j;	/* mode ID */
	
	/* print UDS mapping */
	Message("\nUDS mapping:\n");
	for (j = 0; j < nTutmamModes; ++j) {
		Message("  Mode %d:\n",j);
		Message("    0th moment (1/kg) (number concentration):                             %d \n",iM0Uds+j*nUdsPerMode);
		Message("    0.667th moment (kg^(-1/3)) (surface area concentration):              %d \n",iM23Uds+j*nUdsPerMode);
		
		if (iFirstM1Uds == iLastM1Uds) {
			Message("    1st moment (dimensionless) (mass concentration):                      %d \n",iFirstM1Uds+j*nUdsPerMode);
		} else {
			Message("    1st moments (dimensionless) (mass concentrations) of all species:     %d - %d \n",iFirstM1Uds+j*nUdsPerMode,iLastM1Uds+j*nUdsPerMode);
		}
		
		Message("\n");
	}
	
	/* print UDM mapping */
	Message("\nUDM mapping:\n");
	for (j = 0; j < nTutmamModes; ++j) {
		Message("  Mode %d:\n",j);
		Message("    Number concentration (#/m^3):                                         %d \n",iNumberConcUdm+j*nUdmPerMode);
		
		if (j == powerLawDistribution) {
			Message("    D2 diameter (m):                                                      %d \n",iCmdUdm+j*nUdmPerMode);	
			Message("    Alpha:                                                                %d \n",iLn2sUdm+j*nUdmPerMode);
		} else {
			Message("    Count median diameter (m):                                            %d \n",iCmdUdm+j*nUdmPerMode);
			Message("    Ln^2(Geometric standard deviation):                                   %d \n",iLn2sUdm+j*nUdmPerMode);
		}
		
		Message("    Particle density (kg/m^3):                                            %d \n",iParticleDensityUdm+j*nUdmPerMode);
		Message("    Particle temperature (K):                                             %d \n",iParticleTemperatureUdm+j*nUdmPerMode);
		Message("    Water kappa:                                                          %d \n",iWaterKappaUdm+j*nUdmPerMode);
		Message("    Equilibrium rh / Fluid rh:                                            %d \n",iRhEqPerRhUdm+j*nUdmPerMode);
		
		Message("    Condensation source term for 0.667th moment (kg^(2/3)/(s m^3)):       %d \n",iCondensationM23Udm+j*nUdmPerMode);
		
		if (iCondensationM1FirstUdm == iCondensationM1LastUdm) {
			Message("    Condensation source term for 1st moment (kg/(s m^3)):                 %d \n",iCondensationM1FirstUdm+j*nUdmPerMode);
		} else {
			Message("    Condensation source term for 1st moments of all species (kg/(s m^3)): %d - %d \n",iCondensationM1FirstUdm+j*nUdmPerMode,iCondensationM1LastUdm+j*nUdmPerMode);
		}
		
		Message("    Coagulation source term for 0th moment (1/(s m^3)):                   %d \n",iCoagulationM0Udm+j*nUdmPerMode);
		Message("    Coagulation source term for 0.667th moment (kg^(2/3)/(s m^3)):        %d \n",iCoagulationM23Udm+j*nUdmPerMode);
		
		Message("\n");
	}
	
	Message("  Mode 0\n");
	Message("    Nucleation source term for 0th moment (1/(s m^3)):                    %d \n",iNucleationM0Udm);
	Message("    Nucleation source term for 0.667th moment (kg^(2/3)/(s m^3)):         %d \n",iNucleationM23Udm);
	
	if (iNucleationM1FirstUdm == iNucleationM1LastUdm) {
		Message("    Nucleation source term for 1st moment (kg/(s m^3)):                   %d \n\n",iNucleationM1FirstUdm);
	} else {
		Message("    Nucleation source term for 1st moments of all species (kg/(s m^3)):   %d - %d \n\n",iNucleationM1FirstUdm,iNucleationM1LastUdm);
	}
	
	#if powerLawDistribution > -1
		#if nTutmamModes > 1
			Message("  Transfer from mode 0 (PL) to 1 (LN)\n");
			Message("    For 0th moment (1/(s m^3)):                                           %d \n",iTransferPl2LnM0Udm);
			Message("    For 0.667th moment (kg^(2/3)/(s m^3)):                                %d \n",iTransferPl2LnM23Udm);
			Message("    For 1st moment (kg/(s m^3)):                                          %d \n\n",iTransferPl2LnM1Udm);
		#endif
	#endif
			
	return;	
}

/* This construct matrices containing some settings from settings vectors */
void tutmam_construct_settings_matrices() {
	int i;		/* ID for particle species */
	int iF;		/* ID for Fluent species */
	int ph;		/* ID for particle phase */
	
	for (ph = 0; ph < nTutmamPhases; ++ph) {	/* looping over all particle phases */
		for (i = 0; i < nTutmamSpecies; ++i) {	/* looping over all particle species */
			if (phaseIdVector[i] == ph) {
				phaseMatrix[ph][i] = 1;
			}
		}
	}
	
	for (iF = 0; iF < nFluentSpecies; ++iF) {
		tutmamSpeciesIdVector[iF] = -1;
		
		for (i = 0; i < nTutmamSpecies; ++i) {
			if (fluentSpeciesIdVector[i] == iF) {
				tutmamSpeciesIdVector[iF] = i;
				break;
			}
		}
	}
	
	return;
}
