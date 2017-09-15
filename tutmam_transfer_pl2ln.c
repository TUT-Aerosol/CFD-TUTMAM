/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	This file includes tranfer functions from power law distribution to log-normal distribution
*	You should not need to change this
*/
#include  "udf.h"
#include  "tutmam_constants.h"
#include  "tutmam_settings.h"
#include  "tutmam_func.h"
#include  "tutmam_diffusion.h"
#include  "tutmam_coagulation.h"
#include  "tutmam_transfer_pl2ln.h"
#include  "tutmam_condensation.h"
#include  "tutmam_material_func.h"

/* Calculating self-coagulational rate source terms and storing them to sourceTerm vector */
void calculate_self_coagulational_transfer_rate(real *sourceTerm, cell_t c, Thread *t) {
	real rhoP;					/* particle density (kg/m^3) */
	real tempFluid;				/* fluid temperature (K) */
	
	real integralN1,integralN2;	/* value of integrated part of the function for 0th moment (m^3/s) */
	real integralS1,integralS2;	/* value of integrated part of the function for 0.667th moment (kg^(2/3) m^3/s) */
	real integralM1,integralM2;	/* value of integrated part of the function for 1st moment (kg m^3/s) */
	int iIntegral1,iIntegral2;	/* indices used in integrations */
	
	real ntot;					/* particle number concentration (1/m^3) */
	real D2;					/* D2 (m) */
	real alpha;					/* alpha */
	real dp1,dp2;				/* current particle diameter used in integration (m) */
	real cc1,cc2;
	real dCoag;					/* critical diameter in coagulation (m) */
	real quadratureFactor;
	real factor1;
	real factor2;				/* (1/m^6) */
	
	real diffIndep;/* diffusion coefficient of particle: size-independent part (m^3/s) */
	real diffP1,diffP2;
	real mParticle1,mParticle2;	/* masses of a particle (kg) */
	real correctionFactor;
	real coagCoeff;

	ntot = C_NTOT(c,t,powerLawDistribution);
	D2 = C_CMD(c,t,powerLawDistribution);
	alpha = C_LN2S(c,t,powerLawDistribution);
	rhoP = C_R_PARTICLE(c,t,powerLawDistribution);
	tempFluid = C_T(c,t);
	diffIndep = diff_indep(c,t);
	
	if (D2/powerLawD1-1.0 < 0.001) {
		diffP1 = diffIndep*diff_dep(c,t,powerLawD1);
		
		if (transitionRegimeCorrectionFactorForCoagulationLaw == 1) {
			correctionFactor = fuchs_sutugin_for_coagulation(powerLawD1,powerLawD1,diffP1,diffP1,tempFluid,rhoP,rhoP);
		} else if (transitionRegimeCorrectionFactorForCoagulationLaw == 2) {
			correctionFactor = dahneke_for_coagulation(powerLawD1,powerLawD1,diffP1,diffP1,tempFluid,rhoP,rhoP);
		} else {
			Error("Illegal transitionRegimeCorrectionFactorForCoagulationLaw\n");
			return;
		}
		
		sourceTerm[0] = SQR(ntot)*TUTMAM_PI*coagulation_coefficient(powerLawD1,powerLawD1,diffP1,diffP1,correctionFactor);
		sourceTerm[1] = SQR(ntot)*TUTMAM_PI*coagulation_coefficient(powerLawD1,powerLawD1,diffP1,diffP1,correctionFactor)*TUTMAM_PI623*pow(2*rhoP,TUTMAM_23)*SQR(powerLawD1);
		sourceTerm[2] = SQR(ntot)*TUTMAM_PI*coagulation_coefficient(powerLawD1,powerLawD1,diffP1,diffP1,correctionFactor)*2*TUTMAM_PI6*rhoP*CBC(powerLawD1);
		
		return;
	}

	cc1 = alpha*log(D2/powerLawD1);
	quadratureFactor = quadrature_factor_olin(cc1);
	factor2 = SQR(ntot)*TUTMAM_PI*quadratureFactor;
				
	/* integration starts from zero */
	integralN1 = 0.0;
	integralS1 = 0.0;
	integralM1 = 0.0;
	for (iIntegral1 = 0; iIntegral1 < 4; ++iIntegral1) {
		dp1 = gauss_olin_abscissas_dp(iIntegral1,D2);
		mParticle1 = TUTMAM_PI6*rhoP*CBC(dp1);
		diffP1 = diffIndep*diff_dep(c,t,dp1);
		
		dCoag = tutmam_lower_limit(pow(CBC(D2)-CBC(dp1),TUTMAM_13),powerLawD1);
		cc2 = alpha*log(D2/dCoag);
		
		integralN2 = 0.0;
		integralS2 = 0.0;
		integralM2 = 0.0;
		
		for (iIntegral2 = 0; iIntegral2 < 4; ++iIntegral2) {
			dp2 = gauss_olin_abscissas_dp_dCoag(iIntegral2,dCoag,D2);
			mParticle2 = TUTMAM_PI6*rhoP*CBC(dp2);
			diffP2 = diffIndep*diff_dep(c,t,dp2);
			
			if (transitionRegimeCorrectionFactorForCoagulationLaw == 1) {
				correctionFactor = fuchs_sutugin_for_coagulation(dp1,dp2,diffP1,diffP2,tempFluid,rhoP,rhoP);
			} else if (transitionRegimeCorrectionFactorForCoagulationLaw == 2) {
				correctionFactor = dahneke_for_coagulation(dp1,dp2,diffP1,diffP2,tempFluid,rhoP,rhoP);
			} else {
				Error("Illegal transitionRegimeCorrectionFactorForCoagulationLaw\n");
				return;
			}
			
			coagCoeff = coagulation_coefficient(dp1,dp2,diffP1,diffP2,correctionFactor);
			
			integralN2 += gauss_olin_weights(iIntegral2,cc2)*coagCoeff;
			integralS2 += gauss_olin_weights(iIntegral2,cc2)*coagCoeff*pow(mParticle1+mParticle2,TUTMAM_23);
			integralM2 += gauss_olin_weights(iIntegral2,cc2)*coagCoeff*(mParticle1+mParticle2);
		}
		
		if (fabs(cc2) < 0.001) {
			factor1 = gauss_olin_weights(iIntegral1,cc1)*log(D2/dCoag)/log(D2/powerLawD1);
		
		} else {
			factor1 = gauss_olin_weights(iIntegral1,cc1)*cc2/(exp(cc2)*(1.0-exp(-cc1)));
		}

		integralN1 += integralN2*factor1; 
		integralS1 += integralS2*factor1; 
		integralM1 += integralM2*factor1; 
	}

	sourceTerm[0] = factor2*integralN1; /* (1/(m^3 s)) */ 
	sourceTerm[1] = factor2*integralS1; /* (kg^(2/3)/(m^3 s)) */
	sourceTerm[2] = factor2*integralM1; /* (kg/(m^3 s)) */
	
	return;
}

/* Use robustness coagulation model to coagulational transfer */
void use_robustness_coagulation_model_transfer(real *sourceTerm, cell_t c, Thread *t) {
	real x[ND_ND];
	real deltaT;
	real ntot;
	real stot;
	real mtot = 0.0;
	int iSpecies;
	
	if (coagulationRobustnessModel == 0) {
		return;
	}
	
	ntot = C_UDMI(c,t,powerLawDistribution);
	stot = C_UDS_M23(c,t,powerLawDistribution)*C_R(c,t);
	
	for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) {
		mtot += C_UDS_M1(c,t,iSpecies,powerLawDistribution)*C_R(c,t);
	}
	
	if (ntot < minNumberConc || stot < minSurfaceConc || mtot < minMassConc) {
		return;
	}
	
	#if RP_2D
		C_CENTROID(x,c,t);
		deltaT = sqrt(C_VOLUME(c,t)/x[1]/(SQR(C_U(c,t))+SQR(C_V(c,t))));
	#endif
	
	#if RP_3D
		deltaT = (2.04*pow(C_VOLUME(c,t),TUTMAM_13))/sqrt(SQR(C_U(c,t))+SQR(C_V(c,t))+SQR(C_W(c,t)));
	#endif	
	
	sourceTerm[0] /= (1.0 + sourceTerm[0]*deltaT/ntot);
	sourceTerm[1] /= (1.0 + sourceTerm[1]*deltaT/stot);
	sourceTerm[2] /= (1.0 + sourceTerm[2]*deltaT/mtot);
	
	return;
}

/* Calculating intramodal condensation rate source terms and adding them to sourceTerm vector */
void calculate_inter_modal_condensation_rate(real *sourceTerm, cell_t c, Thread *t) {
	real ntot;										/* number concentration (1/m^3) */
	real rhoP;										/* particle density (kg/m^3) */
	real D2;										/* D2 (m) */
	real massGrowthRate = 0.0;						/* total mass growth rate of a particle (kg/s) */
	real massGrowthRateIndep;						/* size-independent part of mass growth rate of a particle (s^2/m^2) */
	real massGrowthRateDep; 						/* size-dependent part of mass growth rate of a particle (kg m^2/s^3) */
	int i,iSpecies;									/* particle species ID */
	int ph;											/* particle phase ID */
	real alpha;										/* alpha */
	real d;											/* D2/D1 */
	real tempFluid;									/* fluid temperature (K) */
	real tempParticle;								/* particle temperature (K) */
	real pressure;									/* pressure (Pa) */
	real factor;									/* factor used in power law distribution */
	real diffCoeffGas;								/* diffusion coefficient of vapor (m^2/s) */
	real diffCoeffGasDry;							/* diffusion coefficient of dry vapor (m^2/s) */
	real diffCoeffParticle;							/* diffusion coefficient of particle (m^2/s) */
	real dMolecule;									/* vapor molecule diameter (m) */
	real mParticle;									/* particle mass (kg) */
	int iFluentSpecies;								/* Fluent species ID */
	real moleFractionInFluid;						/* mole fraction in fluid */
	real moleFractionInParticle;					/* mole fraction in particle */
	real actCoeff;									/* activity coefficient */
	real fuchsSutuginForMass;						/* Fuchs-Sutugin correction factor */
	real phaseActivity = 1.0;						/* phase activity */
	real moleFractionsInPhase[nTutmamSpecies];		/* mole fractions of all species in a phase */
	real massFractionsInPhase[nTutmamSpecies];		/* mass fractions of all species in a phase */
	real volumeFractionsInParticle[nTutmamPhases];	/* volume fractions of all phases in a particle */
	real saturationVaporPressure;					/* saturation vapor pressure (Pa) */
	real kelvinFactor;								/* Kelvin factor */
	real kelvinDiameter;							/* Kelvin diameter (m) */
	real kappa = 0.0;								/* kappa in water-equilibrium calculation */
	real condensationMultiplier = 0.0;				/* water condensation multiplier */
	int iFluentSpeciesWater = -1;					/* Fluent species ID for water */
	int iFluentSpeciesSulfuricAcid = -1;			/* Fluent species ID for sulfuric acid */
	int iSpeciesWater = -1;							/* species ID for water */
	real relativeHumidity = 0.5;					/* rh */
	
	ntot = C_NTOT(c,t,powerLawDistribution);
	rhoP = C_R_PARTICLE(c,t,powerLawDistribution);
	D2 = C_CMD(c,t,powerLawDistribution);
	alpha = C_LN2S(c,t,powerLawDistribution);
	d = D2/powerLawD1;
	tempFluid = C_T(c,t);
	pressure = C_P_TOT(c,t);
	
	if (ntot < minNumberConc) {
		return;
	}
	
	if (latentHeatOfCondensation == 1) { /* latent heat of condensation taken into account */
		tempParticle = C_T_PARTICLE(c,t,powerLawDistribution);
		if (tempParticle < 1.0) { /* if particle temperature is not stored yet */
			tempParticle = tempFluid;
		}

	} else { /* latent heat of condensation neglected */
		tempParticle = tempFluid;
	}
	
	iFluentSpeciesSulfuricAcid = SV_SpeciesIndex("h2so4"); 
	
	if (waterEq == 1 || iFluentSpeciesSulfuricAcid > -1) { /* if water-equilibrium calculation is on, or if h2so4 exists */
		iFluentSpeciesWater = SV_SpeciesIndex("h2o"); /* find water ID in Fluent */
		if (iFluentSpeciesWater == -1) {
			Error("h2o not found! Switch water equilibrium calculation off.\n");
			return;
		}
		
		iSpeciesWater = tutmamSpeciesIdVector[iFluentSpeciesWater];
	}
		
	if (waterEq == 1) {/* if water-equilibrium calculation is on */
		kappa = C_WATER_KAPPA(c,t,powerLawDistribution);
		condensationMultiplier = condensation_multiplier(c,t,powerLawDistribution);
	}
	
	if (phaseActivityModel == 1) {
		for (ph = 0; ph < nTutmamPhases; ++ph) {			
			volumeFractionsInParticle[ph] = C_VPH_PARTICLE(c,t,ph,powerLawDistribution);
		}
	}
	
	/* calculate rh only when sulfuric acid exists */
	if (iFluentSpeciesSulfuricAcid > -1) {
		relativeHumidity = pressure*C_XI(c,t,iSpeciesWater)/saturation_vapor_pressure(tempFluid,iSpeciesWater);
	}
	
	for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) {
		if (iCondensingSpecies[iSpecies] == 1) { /* is it condensing */
			if (waterEq == 0 || fluentSpeciesIdVector[iSpecies] != iFluentSpeciesWater) { /* skip if water */
				massGrowthRateIndep = mass_growth_rate_indep(tempFluid,iSpecies);
				
				dMolecule = molecule_diameter(iSpecies);
				mParticle = TUTMAM_PI6*rhoP*CBC(D2);
				iFluentSpecies = fluentSpeciesIdVector[iSpecies];
				diffCoeffGas = diffusion_coefficient_gas(tempFluid,pressure,iFluentSpecies,relativeHumidity);
				diffCoeffGasDry = diffusion_coefficient_gas(tempFluid,pressure,iFluentSpecies,0.0);
				diffCoeffParticle = diff_dep(c,t,D2)*diff_indep(c,t);
				fuchsSutuginForMass = fuchs_sutugin_for_mass(D2,tempFluid,pressure,iSpecies,diffCoeffGas,diffCoeffGasDry,diffCoeffParticle,dMolecule,mParticle);
				moleFractionInFluid = C_XI(c,t,iSpecies);
				moleFractionInParticle = C_XI_PARTICLE(c,t,iSpecies,powerLawDistribution);
				
				ph = phaseIdVector[iSpecies];
				for (i = 0; i < nTutmamSpecies; ++i) {
					moleFractionsInPhase[i] = C_XI_PHASE(c,t,i,ph,powerLawDistribution);
					massFractionsInPhase[i] = C_YI_PHASE(c,t,i,ph,powerLawDistribution);
				}
				
				if (phaseActivityModel == 1) {
					phaseActivity = phase_activity(volumeFractionsInParticle[ph]);
				}
				
				actCoeff = activity_coefficient(tempParticle,iSpecies,moleFractionsInPhase);
				saturationVaporPressure = saturation_vapor_pressure(tempParticle,iSpecies);
				
				if (kelvinEffect == 1) {
					kelvinDiameter = kelvin_diameter(tempParticle,iSpecies,powerLawDistribution,moleFractionsInPhase,massFractionsInPhase);
					kelvinFactor = exp(tutmam_upper_limit(kelvinDiameter/D2,10.0));
					
				} else {
					kelvinFactor = 1.0;
				}
				
				if (dieselExhaustHCFractionModel == 1 && iSpecies == dieselExhaustHC) { /* if HC in diesel exhaust HC fraction model */
					/* multiplying partial pressure of HC with the fraction of condensing HCs */
					/* subtracting already condesed fraction */
					moleFractionInFluid *= tutmam_limits(0.0,diesel_exhaust_hc_fraction(tempFluid,moleFractionInFluid*pressure/tutmam_lower_limit(phaseActivity,0.001)) - C_YHC_CONDENSED(c,t),1.0);
				}
				
				massGrowthRateDep = mass_growth_rate_dep(D2,fuchsSutuginForMass,dMolecule,diffCoeffGas,diffCoeffParticle,moleFractionInFluid,moleFractionsInPhase[iSpecies],phaseActivity*actCoeff,saturationVaporPressure,pressure,kelvinFactor);
				if (condensationDirectionVector[iSpecies] == 1) {
					massGrowthRateDep = tutmam_positive(massGrowthRateDep);
				} else if (condensationDirectionVector[iSpecies] == -1) {
					massGrowthRateDep = tutmam_negative(massGrowthRateDep);
				}
				
				massGrowthRate += massGrowthRateIndep*massGrowthRateDep; /* calculate total mass growth rate */
				if (waterEq == 1 && iSpecies == waterEqSpecies && iCondensingSpecies[iSpeciesWater] == 1) { /* if water-connected species in water-equilibrium calculation */
					massGrowthRate += kappa*condensationMultiplier*massGrowthRateIndep*massGrowthRateDep; /* add water part */
				} 
			}
		}
	}
	
	massGrowthRate = tutmam_positive(massGrowthRate); /* only growth is taken into account */
	
	if (d-1.0 < whenAlphaIsZero) {
		factor = 3.0;
		
	} else {
		if (fabs(alpha) < whenAlphaIsZero) {
			factor = 1.5*(1.0-pow(d,-2.0))/log(d);
			
		} else {
			if (fabs(alpha + 2.0) < whenAlphaIsZero) {
				factor = -6.0*log(d)/(1.0-d*d);
				
			} else {
				factor = 3.0*alpha/(alpha+2.0)*(1.0-pow(d,-2.0-alpha))/(1.0-pow(d,0.0-alpha));
			}
		}
	}
	
	factor *= interModalCondensationFactor;

	sourceTerm[0] += 0.636619772367581*ntot/(rhoP*CBC(D2))*massGrowthRate*factor; /* source term for 0th moment */
	sourceTerm[1] += 0.413566993932933*ntot/(pow(rhoP,TUTMAM_13)*D2)*massGrowthRate*factor; /* source term for 0.667th moment */
	sourceTerm[2] += TUTMAM_13*ntot*massGrowthRate*factor;	/* source term for total 1st moment */
	
	return;
}

/* Making under-relaxation to sourceTerm vector */
void under_relax_transfer_pl2ln_rate(real *sourceTerm, real uRFTransferPl2Ln, cell_t c, Thread *t) {
	
	sourceTerm[0] = sourceTerm[0]*uRFTransferPl2Ln + C_UDMI(c,t,iUdmTransferPl2LnM0)*(1.0-uRFTransferPl2Ln);
	sourceTerm[1] = sourceTerm[1]*uRFTransferPl2Ln + C_UDMI(c,t,iUdmTransferPl2LnM23)*(1.0-uRFTransferPl2Ln);
	sourceTerm[2] = sourceTerm[2]*uRFTransferPl2Ln + C_UDMI(c,t,iUdmTransferPl2LnM1)*(1.0-uRFTransferPl2Ln);

	return;
}

/* Storing transfer rate source terms to UDMs */
void store_transfer_rate(real *sourceTerm, cell_t c, Thread *t) {
	C_UDMI(c,t,iUdmTransferPl2LnM0) = sourceTerm[0];
	C_UDMI(c,t,iUdmTransferPl2LnM23) = sourceTerm[1];
	C_UDMI(c,t,iUdmTransferPl2LnM1) = sourceTerm[2];
	
	return;
}

/* Source term macro for transfer */
DEFINE_SOURCE(transfer_pl2ln,c,t,dS,eqn)
{	
	real sourceTerm[3] = { 0.0 };	/* source term (1/m^3), (kg^(2/3) / m^3), or (kg/m^3) */
	int iUds;				/* UDS ID number */
	int iUdsReducedToMode0;	/* UDS ID number reduced to mode 0 */
	int j;					/* mode ID */
	int iSpecies;			/* species ID */
	dS[eqn] = 0.0;			/* differential of the source term */

	iUds = eqn-EQ_UDS;
	iUdsReducedToMode0 = iUds % nUdsPerMode;
	j = iUds/nUdsPerMode;
	
	#if powerLawDistribution == -1 || nTutmamModes == 1
		return 0.0;
	
	#else
		
		if (j > 1 || j < 0) { /* if mode 2 or higher, or not uds */
			return 0.0;
		}
		
		if (j == 1) { /* if LN mode */
			if (iUdsReducedToMode0 == iUdsM0) {
				return C_UDMI(c,t,iUdmTransferPl2LnM0);
				
			} else if (iUdsReducedToMode0 == iUdsM23) {
				return C_UDMI(c,t,iUdmTransferPl2LnM23);
				
			} else {
				iSpecies = iUdsReducedToMode0-iUdsM1First;
				return C_UDMI(c,t,iUdmTransferPl2LnM1)*C_YI_PARTICLE(c,t,iSpecies,0);
			}
			
		}
		
		/* if PL mode */
		
		if (iUds == iUdsM0) { /* UDS for 0th moment */
			if (selfCoagulationalTransferProcess == 1) {
				/* calculate self-coagulation transfer rate and store it to sourceTerm */
				calculate_self_coagulational_transfer_rate(sourceTerm,c,t);
				
				/* Use robustness coagulation model to coagulational transfer */
				use_robustness_coagulation_model_transfer(sourceTerm,c,t);
			}
			
			if (interModalCondensationProcess == 1) {
				/* calculate intermodal condensation rate and add it to sourceTerm */
				calculate_inter_modal_condensation_rate(sourceTerm,c,t);
			}
			
			/* make under-relaxation to sourceTerm */
			under_relax_transfer_pl2ln_rate(sourceTerm,uRFTransferPl2Ln,c,t);
						
			/* store rates to UDMs */
			store_transfer_rate(sourceTerm,c,t);

			/* return transfer rate for 0th moment */
			return 0.0 - sourceTerm[0];
			
		} else if (iUds == iUdsM23) { /* UDS for 23th moment */
			/* return transfer rate for 0.667th moment */
			return 0.0 - C_UDMI(c,t,iUdmTransferPl2LnM23);

		} else if (iUds >= iUdsM1First && iUds < nUdsPerMode) { /* UDS for 1st moment */
			iSpecies = iUds-iUdsM1First;
			/* return transfer rate for 1st moment */
			return 0.0 - C_UDMI(c,t,iUdmTransferPl2LnM1)*C_YI_PARTICLE(c,t,iSpecies,0);
		}
		
		return 0.0;
	
	#endif
}

