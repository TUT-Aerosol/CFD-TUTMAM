/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	This file includes particles condensation functions
*	You should not need to change this
*/
#include  "udf.h"
#include  "tutmam_constants.h"
#include  "tutmam_settings.h"
#include  "tutmam_func.h"
#include  "tutmam_condensation.h"
#include  "tutmam_material_func.h"
#include  "tutmam_diffusion.h"

/* Fuchs-Sutugin correction for heat transfer of a particle (Fuchs & Sutugin, 1970) */
real fuchs_sutugin_for_heat(real dp, real temp, real pressure) {
	real kn;	/* knudsen number */
	
	kn = knudsen_number(dp,temp,pressure);
	return (1.0+kn)/(1.333*kn*kn + 1.71*kn + 1.0);
}

/* Fuchs-Sutugin correction for mass transfer of a particle (Fuchs & Sutugin, 1970) */
real fuchs_sutugin_for_mass(real dp, real temp, real pressure, real iSpecies, real diffCoeffGas, real diffCoeffParticle, real dMolecule, real mParticle) {
	real kn = 1.0;	/* knudsen number */

	kn = knudsen_number_for_mass(dp,temp,pressure,iSpecies,diffCoeffGas,diffCoeffParticle,dMolecule,mParticle);
	return (1.0+kn)/(1.333*kn*kn + 1.71*kn + 1.0);
}

/* Knudsen number for mass transfer of species iSpecies */
/* Single-component approximation for mean free path */
real knudsen_number_for_mass(real dp, real temp, real pressure, int iSpecies, real diffCoeffGas, real diffCoeffParticle, real dMolecule, real mParticle) {
	real M;					/* molar mass of species iSpecies (kg/mol) */
	real Mp;				/* molar mass of particle (kg/mol) */
	M = molarMassVector[iSpecies];
	Mp = mParticle*TUTMAM_AVOGADRO;
	
	return 1.303959714852510*(diffCoeffGas+diffCoeffParticle)/((dp+dMolecule)*sqrt(temp*(1.0/M+1.0/Mp)));
}

/* Single particle mass growth rate of species iSpecies (size-independent part) (s^2/m^2) */
real mass_growth_rate_indep(real temp, int iSpecies) {
	real M;				/* molar mass of species iSpecies (kg/mol) */
	M = molarMassVector[iSpecies];
	
	/* size-independent part of mass growth rate (s^2/m^2) */
	/* 2*pi/R = 0.755693750203662 */
	return 0.755693750203662*M/temp;
}

/* Single particle mass growth rate of species iSpecies (size-dependent part) (kg m^2/s^3) */
real mass_growth_rate_dep(real dp, real fuchsSutuginForMass, real dMolecule, real diffCoeffGas, real diffCoeffParticle, real moleFractionInFluid, real moleFractionInPhase, real actCoeff, real saturationVaporPressure, real pressure, real kelvinFactor) {
	return (dp+dMolecule)*fuchsSutuginForMass*(diffCoeffGas+diffCoeffParticle)*(moleFractionInFluid*pressure - moleFractionInPhase*actCoeff*kelvinFactor*saturationVaporPressure);
}

/* Kelvin diameter (m) */
real kelvin_diameter(real temp, int iSpecies, int j, const real *moleFractionsInPhase, const real *massFractionsInPhase) {
	real kelvinDiameter;
	real surfaceTension;
	real rhoLiquid;
	real M;
	int ph = phaseIdVector[iSpecies];
	
	surfaceTension = surface_tension(temp,ph,moleFractionsInPhase);
	M = molarMassVector[iSpecies];
	rhoLiquid = particle_phase_density(temp,ph,massFractionsInPhase);
	
	kelvinDiameter = (4*surfaceTension*M)/(TUTMAM_R*temp*rhoLiquid);
	
	return kelvinDiameter;
}

/* Calculating condensation rate source terms and storing them to sourceTerm vector */
void calculate_condensation_rate(real *sourceTerm, cell_t c, Thread *t, int j) {
	real rhoP;					/* particle density (kg/m^3) */
	real massGrowthRateIndep;	/* size-independent part of mass growth rate (s^2/m^2) */
	real massGrowthRateDep;		/* size-dependent part of mass growth rate (kg m^2/s^3) */
	real tempFluid;				/* fluid temperature (K) */
	real tempParticle;			/* particle temperature (K) */
	real fuchsSutuginForMass;	/* Fuchs-Sutugin correction factor */
	real moleFractionInFluid;	/* iSpecies mole fraction in fluid */
	real actCoeff;				/* iSpecies activity coefficient */
	real phaseActivity = 1.0;	/* activity of phase ph */
	real saturationVaporPressure;/* iSpecies saturation vapor pressure (Pa) */
	real pressure;				/* total pressure (Pa) */
	real dMolecule;				/* molecule diameter (m) */
	
	real integral;				/* value of integrated part of the function for 0.667th moment (kg/(s m)) */
	real integralM1;			/* value of integrated part of the function for 1st moment (kg/s^2) */
	int iIntegral;				/* index used in integration */
	real sumOverSpecies;		/* sum over all particle species (kg/(s m)) */
	int iSpecies;				/* index used in summation over all particle species */
	
	real cmd;					/* count median diameter (m) */
	real ln2s;					/* ln^2(GSD) */
	real dp;					/* current particle diameter used in integration (m) */
	real x;						/* (ln(dp/cmd)) / (sqrt(2)*ln(GSD)) */
	real cc = 0.0;
	real quadratureFactor;
	
	real diffCoeffGas;			/* laminar diffusion coefficient of gas iSpecies (m^2/s) */
	int iFluentSpecies;			/* ID for iSpecies in Fluent */
	real diffCoeffParticleIndep;/* diffusion coefficient of particle: size-independent part (m^3/s) */
	real diffCoeffParticleDep;	/* diffusion coefficient of particle: size-dependent part (1/m) */
	real diffCoeffParticle;		/* diffusion coefficient of particle (m^2/s) */
	real mParticle;				/* mass of a particle (kg) */
	real kelvinFactor = 1.0;	/* Kelvin factor */
	real kelvinDiameter = 0.0;	/* Kelvin diameter (m) */
	
	int iFluentSpeciesWater = -1;
	int iSpeciesWater = -1;
	real kappa = 0.0;
	real condensationMultiplier = 0.0;
	
	real moleFractionsInPhase[nTutmamSpecies];	/* mole fractions of all species in a phase */
	real massFractionsInPhase[nTutmamSpecies];	/* mass fractions of all species in a phase */
	real volumeFractionsInParticle[nTutmamPhases];	/* volume fractions of all phases in a particle */
	int i,ph;	

	cmd = C_CMD(c,t,j);
	ln2s = C_LN2S(c,t,j);
	rhoP = C_R_PARTICLE(c,t,j);
	tempFluid = C_T(c,t);
	pressure = C_P_TOT(c,t);
	diffCoeffParticleIndep = diff_indep(c,t);
	
	if (j == powerLawDistribution) {
		cc = ln2s*log(cmd/powerLawD1);
		quadratureFactor = quadrature_factor_olin(cc);
				
	} else {
		quadratureFactor = 1.0/TUTMAM_SQRTPI;
	}
		
	if (waterEq == 1) { /* if water-equilibrium calculation is on */
	
		/* find Fluent species ID for water */
		iFluentSpeciesWater = SV_SpeciesIndex("h2o");
		if (iFluentSpeciesWater == -1) {
			Error("h2o not found! Switch water equilibrium calculation off.\n");
			return;
		}
		
		iSpeciesWater = tutmamSpeciesIdVector[iFluentSpeciesWater];
		condensationMultiplier = condensation_multiplier(c,t,j);
		kappa = C_WATER_KAPPA(c,t,j);
	}
	
	if (latentHeatOfCondensation == 1) { /* latent heat of condensation taken into account */
		tempParticle = C_T_PARTICLE(c,t,j);
		if (tempParticle < 1.0) { /* if particle temperature is not stored yet */
			tempParticle = tempFluid;
		}

	} else { /* latent heat of condensation neglected */
		tempParticle = tempFluid;
	}
	
	if (phaseActivityModel == 1) {
		for (ph = 0; ph < nTutmamPhases; ++ph) {			
			volumeFractionsInParticle[ph] = C_VPH_PARTICLE(c,t,ph,j);
		}
	}

	
	/* summation starts from zero */
	sumOverSpecies = 0.0;
			
	for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) {
		if (iCondensingSpecies[j*nTutmamSpecies+iSpecies] == 1) { /* if the species is condensing */
			if (waterEq == 0 || fluentSpeciesIdVector[iSpecies] != iFluentSpeciesWater) { /* skip water in water-equilibrium calculation */
			
				ph = phaseIdVector[iSpecies];
		
				for (i = 0; i < nTutmamSpecies; ++i) {
					moleFractionsInPhase[i] = C_XI_PHASE(c,t,i,ph,j);
					massFractionsInPhase[i] = C_YI_PHASE(c,t,i,ph,j);
				}
			
				moleFractionInFluid = C_XI(c,t,iSpecies);
				dMolecule = molecule_diameter(iSpecies);
				iFluentSpecies = fluentSpeciesIdVector[iSpecies];
				diffCoeffGas = diffusion_coefficient_gas(tempFluid,pressure,iFluentSpecies);
				actCoeff = activity_coefficient(tempParticle,iSpecies,moleFractionsInPhase);
				
				if (phaseActivityModel == 1) {
					phaseActivity = phase_activity(volumeFractionsInParticle[ph]);
				}
				
				saturationVaporPressure = saturation_vapor_pressure(tempParticle,iSpecies);
				
				if (kelvinEffect == 1) {
					kelvinDiameter = kelvin_diameter(tempParticle,iSpecies,j,moleFractionsInPhase,massFractionsInPhase);
				}
				
				if (dieselExhaustHCFractionModel == 1 && iSpecies == dieselExhaustHC) { /* if HC in diesel exhaust HC fraction model */
					/* multiplying partial pressure of HC with the fraction of condensing HCs */
					/* subtracting already condesed fraction */
					moleFractionInFluid *= tutmam_limits(0.0,diesel_exhaust_hc_fraction(tempFluid,moleFractionInFluid*pressure/tutmam_lower_limit(phaseActivity,0.001)) - C_YHC_CONDENSED(c,t),1.0);
				}
				
				/* size-independent parts of growth rates */
				massGrowthRateIndep = mass_growth_rate_indep(tempFluid,iSpecies);
							
				/* integration starts from zero */
				integral = 0.0;
				integralM1 = 0.0;
				
				if (j == powerLawDistribution) {
					for (iIntegral = 0; iIntegral < 4; ++iIntegral) {
						dp = gauss_olin_abscissas_dp(iIntegral,cmd);
						mParticle = TUTMAM_PI6*rhoP*CBC(dp);
						
						diffCoeffParticleDep = diff_dep(c,t,dp);
						diffCoeffParticle = diffCoeffParticleDep*diffCoeffParticleIndep;
						
						if (kelvinEffect == 1) {
							kelvinFactor = exp(tutmam_upper_limit(kelvinDiameter/dp,50.0));
						}
						
						fuchsSutuginForMass = fuchs_sutugin_for_mass(dp,tempFluid,pressure,iSpecies,diffCoeffGas,diffCoeffParticle,dMolecule,mParticle);
						
						massGrowthRateDep = mass_growth_rate_dep(dp,fuchsSutuginForMass,dMolecule,diffCoeffGas,diffCoeffParticle,moleFractionInFluid,moleFractionsInPhase[iSpecies],phaseActivity*actCoeff,saturationVaporPressure,pressure,kelvinFactor);
						if (condensationDirectionVector[iSpecies] == 1) {
							massGrowthRateDep = tutmam_positive(massGrowthRateDep);
						} else if (condensationDirectionVector[iSpecies] == -1) {
							massGrowthRateDep = tutmam_negative(massGrowthRateDep);
						}
							
						/* weights are defined by Gauss-Hermite quadrature */
						integral += gauss_olin_weights(iIntegral,cc)/dp*massGrowthRateIndep*massGrowthRateDep; /* (kg m^2/s^3) */
						integralM1 += gauss_olin_weights(iIntegral,cc)*massGrowthRateDep; /* (kg m^2/s^3) */
					}
				
				} else {
					/* integration is done in amount of gaussHermiteLevel steps */
					/* gaussHermiteLevel comes from the settings in CFD-TUTMAM GUI */
					for (iIntegral = 0; iIntegral < gaussHermiteLevel; ++iIntegral) {
						x = gauss_hermite_abscissas(iIntegral); /* x-variable is defined by Gauss-Hermite quadrature */
						dp = cmd*exp(x*sqrt(2.*ln2s)); /* converting x back to dp */
						mParticle = TUTMAM_PI6*rhoP*CBC(dp);
						
						diffCoeffParticleDep = diff_dep(c,t,dp);
						diffCoeffParticle = diffCoeffParticleDep*diffCoeffParticleIndep;
						
						if (kelvinEffect == 1) {
							kelvinFactor = exp(tutmam_upper_limit(kelvinDiameter/dp,50.0));
						}
						
						fuchsSutuginForMass = fuchs_sutugin_for_mass(dp,tempFluid,pressure,iSpecies,diffCoeffGas,diffCoeffParticle,dMolecule,mParticle);
						
						massGrowthRateDep = mass_growth_rate_dep(dp,fuchsSutuginForMass,dMolecule,diffCoeffGas,diffCoeffParticle,moleFractionInFluid,moleFractionsInPhase[iSpecies],phaseActivity*actCoeff,saturationVaporPressure,pressure,kelvinFactor);
						
						if (condensationDirectionVector[iSpecies] == 1) {
							massGrowthRateDep = tutmam_positive(massGrowthRateDep);
						} else if (condensationDirectionVector[iSpecies] == -1) {
							massGrowthRateDep = tutmam_negative(massGrowthRateDep);
						}
							
						/* weights are defined by Gauss-Hermite quadrature */
						integral += gauss_hermite_weights(iIntegral)/dp*massGrowthRateIndep*massGrowthRateDep; /* (kg m^2/s^3) */
						integralM1 += gauss_hermite_weights(iIntegral)*massGrowthRateDep; /* (kg m^2/s^3) */
					}
				}
				
			
				/* summing over species for 0.667th moment */
				sumOverSpecies += integral;
									
				/* storing sourceTerms for 1st moments */
				sourceTerm[1+iSpecies] = C_NTOT(c,t,j)*massGrowthRateIndep*integralM1*quadratureFactor; /* (kg/m^3 s) */
				
				if (waterEq == 1 && iSpecies == waterEqSpecies && iCondensingSpecies[j*nTutmamSpecies+iSpeciesWater] == 1) { /* if species that is connected to water-equilibrium calculation */
					/* storing sourceTerms for water */
					sourceTerm[1+iSpeciesWater] = kappa*condensationMultiplier*sourceTerm[1+iSpecies];
					
					/* addind water part to sumOverSpecies */
					sumOverSpecies += kappa*condensationMultiplier*integral;
				}
			}
		}
	}
	
	/* storing sourceTerm for 0.667th moment */
	sourceTerm[0] = C_NTOT(c,t,j)*TUTMAM_23PI6M13*sumOverSpecies*pow(rhoP,-TUTMAM_13)*quadratureFactor; /* (kg^(2/3)/m^3 s) */

	return;
}

/* Making under-relaxation to sourceTerm vector */
void under_relax_condensation_rate(real *sourceTerm, real uRFCondensationM23, const real *uRFCondensationM1, cell_t c, Thread *t, int j) {
	int iSpecies;	/* particle species ID */
	
	sourceTerm[0] = sourceTerm[0]*uRFCondensationM23 + C_UDMI(c,t,iUdmCondensationM23+j*nUdmPerMode)*(1.0-uRFCondensationM23);
	
	for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) {
		if (iCondensingSpecies[j*nTutmamSpecies+iSpecies] == 1) {
			sourceTerm[iSpecies+1] = sourceTerm[iSpecies+1]*uRFCondensationM1[iSpecies] + C_UDMI(c,t,iUdmCondensationM1First+j*nUdmPerMode+iSpecies)*(1.0-uRFCondensationM1[iSpecies]);
	
		} else {
			sourceTerm[iSpecies+1] = 0.0;
		}
	}

	return;
}

/* Storing condensation rate source terms to UDMs */
void store_condensation_rate(const real *sourceTerm, cell_t c, Thread *t, int j) {
	int iSpecies;	/* particle species ID */
	
	C_UDMI(c,t,iUdmCondensationM23+j*nUdmPerMode) = sourceTerm[0];
	
	for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) {
		C_UDMI(c,t,iUdmCondensationM1First+j*nUdmPerMode+iSpecies) = sourceTerm[iSpecies+1];
	}

	return;
}

/* Source term macro for condensation */
DEFINE_SOURCE(condensation,c,t,dS,eqn)
{	
	real sourceTerm[1+nTutmamSpecies] = { 0.0 };	/* source term (kg^(2/3)/(s m^3)) or (kg/(s m^3)) */
	real sourceTermVapor = 0.0;
	int iSpecies;			/* species ID number */
	int iFluentSpecies;		/* Fluent species ID number */
	int iUds;				/* UDS ID number */
	int iUdsReducedToMode0;	/* UDS ID number reduced to mode 0 */
	int j;					/* mode ID */
	dS[eqn] = 0.0;			/* differential of the source term */

	iUds = eqn-EQ_UDS;
	iUdsReducedToMode0 = iUds % nUdsPerMode;
	j = iUds/nUdsPerMode;
	
	if (iUds < 0 && condensationProcess == 1) { /* species */
		iFluentSpecies = eqn-EQ_SPECIES;
		iSpecies = tutmamSpeciesIdVector[iFluentSpecies];
		
		if (iSpecies > -1) {
			for (j = 0; j < nTutmamModes; ++j) {
				if (iCondensingSpecies[j*nTutmamSpecies+iSpecies] == 1) {
					/* returning opposite number of the source term calculated for particle modes */
					sourceTermVapor -= C_UDMI(c,t,iUdmCondensationM1First+j*nUdmPerMode+iSpecies);
				}
			}
			
			return sourceTermVapor;
			
		} else {
			return 0.0;
		}
	
	} else if (iUdsReducedToMode0 == iUdsM23) { /* UDS for 0.667th moment */
		if (condensationProcess == 1) {
			if (condensationComputing == 1) {
				/* calculate condensation rate and store it to sourceTerm */
				calculate_condensation_rate(sourceTerm,c,t,j);

				/* make under-relaxation to sourceTerm */
				under_relax_condensation_rate(sourceTerm,uRFCondensationM23,uRFCondensationM1,c,t,j);

			} else { /* computing off */
				return C_UDMI(c,t,iUdmCondensationM23+j*nUdmPerMode);
			}
		}
		
		/* store condensation rate to UDMs */
		store_condensation_rate(sourceTerm,c,t,j);

		/* return condensation rate for 0.667th moment */
		return sourceTerm[0];

	} else if (iUdsReducedToMode0 >= iUdsM1First && iUdsReducedToMode0 < nUdsPerMode) { /* UDS for 1st moment */
		iSpecies = iUdsReducedToMode0-iUdsM1First;
		/* obtain condensation rate source term from UDM */
		return C_UDMI(c,t,iUdmCondensationM1First+j*nUdmPerMode+iSpecies);
		
	} else { /* all other equations */
		return 0.0;
	}
}

/* calculates relative humidity */
real tutmam_rh(cell_t c, Thread *t, int iSpeciesWater) {
	real moleFractionInFluid;		/* mole fraction of water in fluid */
	real saturationVaporPressure;	/* saturation vapor pressure of water (Pa) */
	real temp;						/* temperature (K) */
	real pressure;					/* total pressure (Pa) */
	real rh;						/* rh */
	
	pressure = C_P_TOT(c,t);
	temp = C_T(c,t);
	moleFractionInFluid = C_XI(c,t,iSpeciesWater);
	saturationVaporPressure = saturation_vapor_pressure(temp,iSpeciesWater);
	
	rh = pressure*moleFractionInFluid/saturationVaporPressure;
	return rh;
}

/* calculates water equilibrium relative humidity */
real rh_eq(cell_t c, Thread *t, int iSpeciesWater, int j, const real *M1, real ntot, real ln2s) {
	real rhEq;					/* rh in equilibrium */
	real actCoeff;				/* activity coefficient */
	real phaseActivity = 1.0;	/* phase activity */
	real kelvinFactor = 1.0;	/* Kelvin factor */
	real kelvinDiameter;		/* Kelvin diameter (m) */
	real dp;					/* particle diameter (m) */
	real temp;					/* temperature (K) */
	real M1Total = 0.0;			/* total 1st moment (kg/m^3) */
	real M1Phase = 0.0;			/* 1st moment in the phase ph (kg/m^3) */
	real molesPhase = 0.0;		/* moles in the phase ph (mol/m^3) */
	int iSpecies;				/* species ID */
	int ph;						/* phase ID */
	real rhoP;					/* particle density (kg/m^3) */
	real moleFractionsInPhase[nTutmamSpecies] = { 0.0 };	/* mole fractions of all species in a phase */
	real massFractionsInPhase[nTutmamSpecies] = { 0.0 };	/* mass fractions of all species in a phase */
	
	/* find the phase ID */
	ph = phaseIdVector[iSpeciesWater];
	
	for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) { /* looping over all particle species */
		M1Total += M1[iSpecies]; /* adding mass to total mass */
		
		if (phaseIdVector[iSpecies] == ph) { /* does iSpecies belong to phase ph */
			M1Phase += M1[iSpecies]; /* adding mass to phase mass */
			molesPhase += M1[iSpecies]/molarMassVector[iSpecies]; /* adding moles to phase mole amount */
		}
	}
	
	if (M1Total < minMassConc) {
		C_UDMI(c,t,iRhEqPerRhUdm+j*nUdmPerMode) = 1.0;
		return tutmam_rh(c,t,iSpeciesWater);
	}
	
	if (latentHeatOfCondensation == 1) { /* is latent heat on */
		temp = C_T_PARTICLE(c,t,j); /* using particle temperature */
		if (temp < 1.0) { /* if particle temperature is not stored yet */
			temp = C_T(c,t);
		}
		
	} else { /* latent heat off */
		temp = C_T(c,t); /* using fluid temperature */
	}
	
	/* using stored particle density value */
	rhoP = C_R_PARTICLE(c,t,j);
	
	/* calculate dp */
	dp = TUTMAM_PI6M13*pow(M1Total/(rhoP*ntot),TUTMAM_13); /* diameter of average mass */
	if (j == powerLawDistribution) {
		dp = tutmam_limits(powerLawD1,dp,C_CMD(c,t,j)); 
		
	} else {
		dp = tutmam_limits(minCMD,dp,maxCMD);
	}
	
	if (powerLawDistribution < 0) { /* if log-normal distribution only */
		dp = tutmam_limits(minCMD,dp*exp(-1.5*ln2s),maxCMD); /* cmd */
	}
	
	/* calculate mole and mass fractions */
	for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) {
		if (phaseIdVector[iSpecies] == ph) {
			moleFractionsInPhase[iSpecies] = (M1[iSpecies]/molarMassVector[iSpecies])/molesPhase;
			massFractionsInPhase[iSpecies] = M1[iSpecies]/M1Phase;
		}
	}
	
	if (phaseActivityModel == 1) {
		phaseActivity = phase_activity(M1Phase/M1Total*rhoP/particle_phase_density(temp,ph,massFractionsInPhase));
	}
	
	if (kelvinEffect == 1) {
		kelvinDiameter = kelvin_diameter(temp,iSpeciesWater,j,moleFractionsInPhase,massFractionsInPhase);
		kelvinFactor = exp(kelvinDiameter/dp);
	}
	
	actCoeff = activity_coefficient(temp,iSpeciesWater,moleFractionsInPhase);
	
	rhEq = moleFractionsInPhase[iSpeciesWater]*phaseActivity*actCoeff*kelvinFactor;
	C_UDMI(c,t,iRhEqPerRhUdm+j*nUdmPerMode) = rhEq/tutmam_lower_limit(tutmam_rh(c,t,iSpeciesWater),1.0e-9);
	
	return rhEq;
}

/* calculates water equilibrium factor, kappa */
real water_kappa(cell_t c, Thread *t, int j) {
	real kappa_0;				/* kappa in the beginning */
	real kappa;					/* kappa now */
	real cmd;					/* mode cmd (m) */
	real rho;					/* fluid density (kg/m^3) */
	real ntot;					/* mode Ntot (1/m^3) */
	real M1[nTutmamSpecies];	/* species mass vector (kg/m^3) */
	real M1Water_0 = 0.0;		/* water mass in the beginning (kg/m^3) */
	real M1Water;				/* water mass now (kg/m^3) */
	real M1Total = 0.0;			/* total particle mass (kg/m^3) */
	real C1;					/* temp variable */
	int iSpecies;				/* species ID */
	int iFluentSpeciesWater;	/* Fluent species ID for water */
	int iSpeciesWater;			/* species ID for water */
	real rhEq_0;				/* rheq in the beginning */
	real rhEq;					/* rheq now */
	real rh;					/* rh */
	real ln2s;					/* ln^(GSD) */
	int i = 0;					/* integer used in summation */
	
	cmd = C_CMD(c,t,j);
	if (cmd < minCMDForKappaCalculation) { /* returning 1 for too small particles */
		return 1.0;
	}
	
	ntot = C_NTOT(c,t,j);
	if (ntot < minNumberConc) { /* if no particles */
		return 1.0;
	}
	
	ln2s = C_LN2S(c,t,j);
	
	iFluentSpeciesWater = SV_SpeciesIndex("h2o"); /* find water ID */
	if (iFluentSpeciesWater == -1) {
		Error("h2o not found! Switch off water equilibrium calculation.\n");
		return 1.0;
	}
	
	iSpeciesWater = tutmamSpeciesIdVector[iFluentSpeciesWater];
	
	rho = C_R(c,t);
			
	for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) { /* calculating particle species masses */
		M1[iSpecies] = tutmam_positive(C_UDS_M1(c,t,iSpecies,j))*rho;
		M1Total += M1[iSpecies];
		
		if (iSpecies == iSpeciesWater) { /* is it water */
			M1Water_0 = M1[iSpecies];
		}
	}
	
	if (M1Total < minMassConc || M1Water_0 < minMassConc) { /* too low mass */
		return 1.0;
	}
	
	rh = tutmam_rh(c,t,iSpeciesWater); /* fluid rh */
	rhEq_0 = rh_eq(c,t,iSpeciesWater,j,M1,ntot,ln2s); /* rheq */
	
	rhEq = rhEq_0;
	M1Water = M1Water_0;
	
	while (fabs(rh/rhEq - 1.0) > convergenceCriteriumForKappaCalculation && i < maxIterationsInKappaCalculation) { /* while not converged */
		M1Water *= (1.0-uRFKappaCalculation+uRFKappaCalculation*rh/rhEq); /* multiply water amount */
		M1[iSpeciesWater] = M1Water;
		
		rhEq = rh_eq(c,t,iSpeciesWater,j,M1,ntot,ln2s); /* calculate new rheq */

		++i;
	}

	if (fabs(rhEq-rh) < fabs(rhEq_0-rh)) { /* if better result is obtained */
		C1 = M1Water/M1Water_0;
		kappa_0 = C_WATER_KAPPA(c,t,j);
		
		if (C1 >= 1.0 || condensationDirectionVector[iSpeciesWater] == 1) {
			if (kappa_0 < 0.0) {
				kappa = tutmam_limits(0.0,0.01*C1,maxKappa);
				
			} else {
				kappa = ((1.0-uRFKappa) + uRFKappa*tutmam_limits(0.0,kappa_0*C1,maxKappa)); /* under-relaxation */
			}
			
		} else {
			kappa = ((1.0-uRFKappa) + uRFKappa*tutmam_limits(minKappa,kappa_0-1.0/C1,maxKappa)); /* under-relaxation */
		}
				
	} else {
		kappa = C_WATER_KAPPA(c,t,j);
	}

	return kappa;
}

/* calculates water condensation multiplier */
real condensation_multiplier(cell_t c, Thread *t, int j) {
	real eqMoleFractionWater;	/* mole fraction of water in equilibrium particle */
	real eqMassFractionWater;	/* mass fraction of water in equilibrium particle */
	real temp;					/* temperature (K) */
	real rh;					/* relative humidity */
	real dp;					/* dp (m) */
	int iSpeciesWater;			/* water species ID */
	int iSpecies;				/* species ID */
	int iFluentSpeciesWater;	/* water species ID in Fluent */
	real phaseActivity = 1.0;	/* phase activity */
	int ph;						/* phase ID */
	real M1Total = 0.0;			/* total particle mass (kg/m^3) */
	real ntot;					/* mode Ntot (1/m^3) */
	real rhoP;					/* particle density (kg/m^3) */

	iFluentSpeciesWater = SV_SpeciesIndex("h2o"); /* find water ID */
	if (iFluentSpeciesWater == -1) {
		Error("h2o not found! Switch off water equilibrium calculation.\n");
		return 1.0;
	}
	
	iSpeciesWater = tutmamSpeciesIdVector[iFluentSpeciesWater];
	
	if (phaseActivityModel == 1) {
		ph = phaseIdVector[iSpeciesWater];
		phaseActivity = tutmam_lower_limit(phase_activity(C_VPH_PARTICLE(c,t,ph,j)),0.001);
	}
	
	temp = C_T(c,t);
	rh = tutmam_rh(c,t,iSpeciesWater);

	/* calculate dp */
	if (powerLawDistribution < 0) { /* if log-normal distribution only */
		dp = C_CMD(c,t,j); /* cmd */
		
	} else {
		for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) {
			M1Total += tutmam_positive(C_UDS_M1(c,t,iSpecies,j)); /* (unitless) */
		}
		M1Total *= C_R(c,t); /* (kg/m^3) */
		
		ntot = C_NTOT(c,t,j);
		if (ntot < minNumberConc) {
			return 1.0;
		}
		
		/* using stored particle density value */
		rhoP = C_R_PARTICLE(c,t,j);
		
		dp = TUTMAM_PI6M13*pow(M1Total/(rhoP*ntot),TUTMAM_13); /* diameter of average mass */
		if (j == powerLawDistribution) {
			dp = tutmam_limits(powerLawD1,dp,C_CMD(c,t,j)); 
		} 
	}

	eqMoleFractionWater = eq_mole_fraction_water(temp,rh/phaseActivity,dp);
	eqMoleFractionWater = tutmam_lower_limit(eqMoleFractionWater,0.001);
	
	eqMassFractionWater = 1.0/(1.0+(1.0-eqMoleFractionWater)/eqMoleFractionWater*molarMassVector[waterEqSpecies]/molarMassVector[iSpeciesWater]);
	eqMassFractionWater = tutmam_upper_limit(eqMassFractionWater,0.999);
	return eqMassFractionWater/(1.0-eqMassFractionWater);
}
