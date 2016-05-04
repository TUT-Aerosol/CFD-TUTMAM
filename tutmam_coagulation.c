/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	This file includes particle coagulation functions
*	You should not need to change this
*/
#include  "udf.h"
#include  "tutmam_constants.h"
#include  "tutmam_settings.h"
#include  "tutmam_func.h"
#include  "tutmam_diffusion.h"
#include  "tutmam_coagulation.h"

/* Fuchs-Sutugin correction for coagulation (Fuchs & Sutugin, 1970) */
real fuchs_sutugin_for_coagulation(real dpi, real dpj, real diffPI, real diffPJ, real temp, real rhoPI, real rhoPJ) {
	real kn;	/* knudsen number */
	
	kn = knudsen_number_for_coagulation(dpi,dpj,diffPI,diffPJ,temp,rhoPI,rhoPJ);
	return (1.0+kn)/(1.333*kn*kn + 1.71*kn + 1.0);
}

/* Dahneke correction for coagulation (Dahneke, 1983) */
real dahneke_for_coagulation(real dpi, real dpj, real diffPI, real diffPJ, real temp, real rhoPI, real rhoPJ) {
	real kn;	/* knudsen number */
	
	kn = TUTMAM_23*knudsen_number_for_coagulation(dpi,dpj,diffPI,diffPJ,temp,rhoPI,rhoPJ);
	return (1.0+kn)/(2.0*kn*kn + 2.0*kn + 1.0);
}

/* Knudsen number for coagulation */
real knudsen_number_for_coagulation(real dpi, real dpj, real diffPI, real diffPJ, real temp, real rhoPI, real rhoPJ) {
	return 7.3212235e11/sqrt(temp)*(diffPI+diffPJ)/(dpi+dpj)/sqrt(pow(dpi,-3.0)/rhoPI + pow(dpj,-3.0)/rhoPJ);
}

/* Coagulation coefficient divided by 2pi (m^3/s) */
real coagulation_coefficient(real dpi, real dpj, real diffPI, real diffPJ, real fuchsSutugin) {
	return (dpi+dpj)*(diffPI+diffPJ)*fuchsSutugin;
}

/* Calculating intramodal coagulation rate source terms and storing them to sourceTerm vector */
void calculate_intra_coagulation_rate(real *sourceTerm, cell_t c, Thread *t, int j) {
	real rhoP;					/* particle density (kg/m^3) */
	real tempFluid;				/* fluid temperature (K) */
	
	real integralN1,integralN2;	/* value of integrated part of the function for 0th moment (m^3/s) */
	real integralS1,integralS2;	/* value of integrated part of the function for 0.667th moment (kg^(2/3) m^3/s) */
	int iIntegral1,iIntegral2;	/* indices used in integrations */
	
	real ntot;					/* particle number concentration (1/m^3) */
	real D2;					/* D2 (m) */
	real alpha;					/* alpha */
	real dp1,dp2;				/* current particle diameter used in integration (m) */
	real cc;
	real x1,x2;
	real quadratureFactor;
	real factor1;
	real factor2;				/* (1/m^6) */
	
	real diffIndep;/* diffusion coefficient of particle: size-independent part (m^3/s) */
	real diffP1,diffP2;
	real mParticle1,mParticle2;	/* masses of a particle (kg) */
	real correctionFactor;
	real coagCoeff;

	ntot = C_NTOT(c,t,j);
	if (ntot < minNumberConc || coagulationMatrix[j][j] == 0) {
		return;
	}
	
	D2 = C_CMD(c,t,j);
	alpha = C_LN2S(c,t,j);
	rhoP = C_R_PARTICLE(c,t,j);
	tempFluid = C_T(c,t);
	diffIndep = diff_indep(c,t);
	
	if (j == powerLawDistribution) { /* power law */
		cc = alpha*log(D2/powerLawD1);
		quadratureFactor = quadrature_factor_olin(cc);
		
	} else { /* log-normal */
		quadratureFactor = 1.0/TUTMAM_SQRTPI;
	}
	
	factor2 = TUTMAM_PI*SQR(ntot*quadratureFactor);
	/* integration starts from zero */
	integralN1 = 0.0;
	integralS1 = 0.0;
	
	if (j == powerLawDistribution) { /* power law */
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
			
			sourceTerm[0] = -TUTMAM_PI*SQR(ntot)*coagulation_coefficient(powerLawD1,powerLawD1,diffP1,diffP1,correctionFactor);
			sourceTerm[1] = -0.84206135628927*SQR(ntot)*coagulation_coefficient(powerLawD1,powerLawD1,diffP1,diffP1,correctionFactor)*pow(rhoP,TUTMAM_23)*SQR(powerLawD1);
			
			return;
		}

		for (iIntegral1 = 0; iIntegral1 < 4; ++iIntegral1) {
			dp1 = gauss_olin_abscissas_dp(iIntegral1,D2);
			mParticle1 = TUTMAM_PI6*rhoP*CBC(dp1);
			diffP1 = diffIndep*diff_dep(c,t,dp1);
					
			integralN2 = 0.0;
			integralS2 = 0.0;
			
			for (iIntegral2 = 0; iIntegral2 < 4; ++iIntegral2) {
				dp2 = gauss_olin_abscissas_dp(iIntegral2,D2);
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
				
				integralN2 -= gauss_olin_weights(iIntegral2,cc)*coagCoeff;
				integralS2 += gauss_olin_weights(iIntegral2,cc)*coagCoeff*(pow(mParticle1+mParticle2,TUTMAM_23) - 2*pow(mParticle1,TUTMAM_23));
			}
			
			factor1 = gauss_olin_weights(iIntegral1,cc);

			integralN1 += integralN2*factor1; 
			integralS1 += integralS2*factor1; 
		}
		
	} else { /* log-normal */
		for (iIntegral1 = 0; iIntegral1 < gaussHermiteLevel; ++iIntegral1) {
			x1 = gauss_hermite_abscissas(iIntegral1); /* x-variable is defined by Gauss-Hermite quadrature */
			dp1 = D2*exp(x1*sqrt(2*alpha)); /* converting x back to dp */
			mParticle1 = TUTMAM_PI6*rhoP*CBC(dp1);
			diffP1 = diffIndep*diff_dep(c,t,dp1);
					
			integralN2 = 0.0;
			integralS2 = 0.0;
			
			for (iIntegral2 = 0; iIntegral2 < 4; ++iIntegral2) {
				x2 = gauss_hermite_abscissas(iIntegral2); /* x-variable is defined by Gauss-Hermite quadrature */
				dp2 = D2*exp(x2*sqrt(2*alpha)); /* converting x back to dp */
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
				
				integralN2 -= gauss_hermite_weights(iIntegral2)*coagCoeff;
				integralS2 += gauss_hermite_weights(iIntegral2)*coagCoeff*(pow(mParticle1+mParticle2,TUTMAM_23) - 2*pow(mParticle1,TUTMAM_23));
			}
			
			factor1 = gauss_hermite_weights(iIntegral1);

			integralN1 += integralN2*factor1; 
			integralS1 += integralS2*factor1; 
		}
	}
	
	sourceTerm[0] = factor2*integralN1; /* (1/(m^3 s)) */ 
	sourceTerm[1] = factor2*integralS1; /* (kg^(2/3)/(m^3 s)) */
	
	return;
}

/* Use robustness coagulation model */
void use_robustness_coagulation_model(real *sourceTerm, cell_t c, Thread *t, int j) {
	real x[ND_ND];
	real deltaT;
	real ntot;
	real stot;
	
	if (coagulationRobustnessModel == 0 || powerLawDistribution >= 0) {
		return;
	}
	
	ntot = C_UDMI(c,t,j);
	stot = C_UDS_M23(c,t,j)*C_R(c,t);
	
	if (ntot < minNumberConc || stot < minSurfaceConc) {
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
	
	return;
}

/* Making under-relaxation to sourceTerm vector */
void under_relax_intra_coagulation_rate(real *sourceTerm, real uRFCoagulation, cell_t c, Thread *t, int j) {
	
	sourceTerm[0] = sourceTerm[0]*uRFCoagulation + C_UDMI(c,t,iUdmCoagulationM0+j*nUdmPerMode)*(1.0-uRFCoagulation);
	sourceTerm[1] = sourceTerm[1]*uRFCoagulation + C_UDMI(c,t,iUdmCoagulationM23+j*nUdmPerMode)*(1.0-uRFCoagulation);

	return;
}

/* Storing intramodal coagulation rate source terms to UDMs */
void store_intra_coagulation_rate(const real *sourceTerm, cell_t c, Thread *t, int j) {
	
	C_UDMI(c,t,iUdmCoagulationM0+j*nUdmPerMode) = sourceTerm[0];
	C_UDMI(c,t,iUdmCoagulationM23+j*nUdmPerMode) = sourceTerm[1];

	return;
}

/* Source term macro for intramodal coagulation */
DEFINE_SOURCE(intra_coagulation,c,t,dS,eqn)
{	
	real sourceTerm[2] = {0.0};	/* source term  */
	int iUds;				/* UDS ID number */
	int iUdsReducedToMode0;	/* UDS ID number reduced to mode 0 */
	int j;					/* mode ID */
	dS[eqn] = 0.0;			/* differential of the source term */
	
	iUds = eqn-EQ_UDS;
	iUdsReducedToMode0 = iUds % nUdsPerMode;
	j = iUds/nUdsPerMode;
	
	if (iUdsReducedToMode0 == iUdsM0) {
		if (coagulationProcess == 1 && intraModalCoagulationProcess == 1) {
			if (coagulationComputing == 1 && N_ITER % coagulationIterationSkip == 0) {
				/* calculate coagulation rate and store it to sourceTerm */
				calculate_intra_coagulation_rate(sourceTerm,c,t,j);
				
				/* use robustness coagulation model */
				use_robustness_coagulation_model(sourceTerm,c,t,j);

				/* make under-relaxation to sourceTerm */
				under_relax_intra_coagulation_rate(sourceTerm,uRFCoagulation,c,t,j);
			
			} else {
				return C_UDMI(c,t,iUdmCoagulationM0+j*nUdmPerMode);
			}
		}
		
		/* store coagulation rate to UDMs */
		store_intra_coagulation_rate(sourceTerm,c,t,j);

		/* return coagulation rate for 0th moment */
		return sourceTerm[0];
		
	} else if (iUdsReducedToMode0 == iUdsM23) {
		return C_UDMI(c,t,iUdmCoagulationM23+j*nUdmPerMode);
		
	} else {
		return 0.;
	}
}

/* Returning intermodal coagulation rate source term */
/* For cases with loss terms only */
real inter_coagulation_rate_L(cell_t c, Thread *t, int iUds, int j, int jOther) {
	real rhoPj,rhoPjOther;		/* particle density (kg/m^3) */
	real tempFluid;				/* fluid temperature (K) */
	
	real integral1,integral2;	/* value of integrated part of the function for the moment */
	int iIntegral1,iIntegral2;	/* indices used in integrations */
	
	real cmdj,cmdjOther;		/* cmd or D2 (m) */
	real ln2sj,ln2sjOther;		/* ln2s or alpha */
	real dp1,dp2;				/* current particle diameter used in integration (m) */
	real cc;
	real x1,x2;
	real quadratureFactor;
	real k;						/* k-th moment */
	
	real diffIndep;				/* diffusion coefficient of particle: size-independent part (m^3/s) */
	real diffP1,diffP2;			/* diffusion coefficient of particle (m^2/s) */
	real correctionFactor;		/* transtion regime correction factor */
	real coagCoeff;				/* coagulation coefficient divided by 2pi (m^3/s) */
	real massFraction = 0.0;	/* mass fraction of current species in mode */
	real sourceTerm;

	cmdj = C_CMD(c,t,j);
	cmdjOther = C_CMD(c,t,jOther);
	ln2sj = C_LN2S(c,t,j);
	ln2sjOther = C_LN2S(c,t,jOther);
	rhoPj = C_R_PARTICLE(c,t,j);
	rhoPjOther = C_R_PARTICLE(c,t,jOther);
	tempFluid = C_T(c,t);
	diffIndep = diff_indep(c,t);
	k = k_moment(iUds);
	quadratureFactor = 1.0/TUTMAM_SQRTPI;
	
	if (k > 0.9) { /* k=1 */
		massFraction = C_YI_PARTICLE(c,t,iUds-iUdsM1First,j);
	}
	
	/* integration starts from zero */
	integral1 = 0.0;
	
	if (j == powerLawDistribution) { /* power law */
		cc = ln2sj*log(cmdj/powerLawD1);
		
		for (iIntegral1 = 0; iIntegral1 < 4; ++iIntegral1) {
			dp1 = gauss_olin_abscissas_dp(iIntegral1,cmdj);
			diffP1 = diffIndep*diff_dep(c,t,dp1);
					
			integral2 = 0.0;
			
			for (iIntegral2 = 0; iIntegral2 < gaussHermiteLevel; ++iIntegral2) {
				x2 = gauss_hermite_abscissas(iIntegral2); /* x-variable is defined by Gauss-Hermite quadrature */
				dp2 = cmdjOther*exp(x2*sqrt(2*ln2sjOther)); /* converting x back to dp */
				diffP2 = diffIndep*diff_dep(c,t,dp2);
				
				if (transitionRegimeCorrectionFactorForCoagulationLaw == 1) {
					correctionFactor = fuchs_sutugin_for_coagulation(dp1,dp2,diffP1,diffP2,tempFluid,rhoPj,rhoPjOther);
				} else if (transitionRegimeCorrectionFactorForCoagulationLaw == 2) {
					correctionFactor = dahneke_for_coagulation(dp1,dp2,diffP1,diffP2,tempFluid,rhoPj,rhoPjOther);
				} else {
					Error("Illegal transitionRegimeCorrectionFactorForCoagulationLaw\n");
					return 0.0;
				}
				
				coagCoeff = coagulation_coefficient(dp1,dp2,diffP1,diffP2,correctionFactor);
				integral2 += gauss_hermite_weights(iIntegral2)*coagCoeff;
			}
			
			if (k < 0.1) { /* k=0 */
				integral1 += integral2*gauss_olin_weights(iIntegral1,cc);
				
			} else if (k > 0.9) { /* k=1 */
				integral1 += integral2*gauss_olin_weights(iIntegral1,cc)*CBC(dp1);
				
			} else { /* k=2/3 */
				integral1 += integral2*gauss_olin_weights(iIntegral1,cc)*SQR(dp1);
			}
		}
		
	} else { /* log-normal */
		for (iIntegral1 = 0; iIntegral1 < gaussHermiteLevel; ++iIntegral1) {
			x1 = gauss_hermite_abscissas(iIntegral1); /* x-variable is defined by Gauss-Hermite quadrature */
			dp1 = cmdj*exp(x1*sqrt(2*ln2sj)); /* converting x back to dp */
			diffP1 = diffIndep*diff_dep(c,t,dp1);
					
			integral2 = 0.0;
			
			for (iIntegral2 = 0; iIntegral2 < gaussHermiteLevel; ++iIntegral2) {
				x2 = gauss_hermite_abscissas(iIntegral2); /* x-variable is defined by Gauss-Hermite quadrature */
				dp2 = cmdjOther*exp(x2*sqrt(2*ln2sjOther)); /* converting x back to dp */
				diffP2 = diffIndep*diff_dep(c,t,dp2);
				
				if (transitionRegimeCorrectionFactorForCoagulationLaw == 1) {
					correctionFactor = fuchs_sutugin_for_coagulation(dp1,dp2,diffP1,diffP2,tempFluid,rhoPj,rhoPjOther);
				} else if (transitionRegimeCorrectionFactorForCoagulationLaw == 2) {
					correctionFactor = dahneke_for_coagulation(dp1,dp2,diffP1,diffP2,tempFluid,rhoPj,rhoPjOther);
				} else {
					Error("Illegal transitionRegimeCorrectionFactorForCoagulationLaw\n");
					return 0.0;
				}
				
				coagCoeff = coagulation_coefficient(dp1,dp2,diffP1,diffP2,correctionFactor);
				integral2 += gauss_hermite_weights(iIntegral2)*coagCoeff;
			}
			
			if (k < 0.1) { /* k=0 */
				integral1 += integral2*gauss_hermite_weights(iIntegral1);
				
			} else if (k > 0.9) { /* k=1 */
				integral1 += integral2*gauss_hermite_weights(iIntegral1)*CBC(dp1);
				
			} else { /* k=2/3 */
				integral1 += integral2*gauss_hermite_weights(iIntegral1)*SQR(dp1);
			}
		}
	}

	sourceTerm = quadratureFactor*integral1;
	if (k < 0.1) {
		return sourceTerm;
		
	} else if (k > 0.9) {
		sourceTerm *= TUTMAM_PI6*rhoPj*massFraction;
		
	} else {
		sourceTerm *= TUTMAM_PI623*pow(rhoPj,TUTMAM_23);
	}
	
	return sourceTerm;
}

/* Returning intermodal coagulation rate source term */
/* For cases with loss and gain terms */
real inter_coagulation_rate_LG(cell_t c, Thread *t, int iUds, int j, int jOther) {
	real rhoPj,rhoPjOther;		/* particle density (kg/m^3) */
	real tempFluid;				/* fluid temperature (K) */
	
	real integral1,integral2;	/* value of integrated part of the function for the moment */
	int iIntegral1,iIntegral2;	/* indices used in integrations */
	
	real cmdj,cmdjOther;		/* cmd or D2 (m) */
	real ln2sj,ln2sjOther;		/* ln2s or alpha */
	real dp1,dp2;				/* current particle diameter used in integration (m) */
	real cc = 0.0;
	real x1,x2;
	real quadratureFactor;
	real k;						/* k-th moment */
	
	real diffIndep;				/* diffusion coefficient of particle: size-independent part (m^3/s) */
	real diffP1,diffP2;			/* diffusion coefficient of particle (m^2/s) */
	real correctionFactor;		/* transtion regime correction factor */
	real coagCoeff;				/* coagulation coefficient divided by 2pi (m^3/s) */
	real massFraction = 0.0;	/* mass fraction of current species in mode */
	real sourceTerm;

	cmdj = C_CMD(c,t,j);
	cmdjOther = C_CMD(c,t,jOther);
	ln2sj = C_LN2S(c,t,j);
	ln2sjOther = C_LN2S(c,t,jOther);
	rhoPj = C_R_PARTICLE(c,t,j);
	rhoPjOther = C_R_PARTICLE(c,t,jOther);
	tempFluid = C_T(c,t);
	diffIndep = diff_indep(c,t);
	k = k_moment(iUds);
	
	if (jOther == powerLawDistribution) {
		cc = ln2sjOther*log(cmdjOther/powerLawD1);
		quadratureFactor = quadrature_factor_olin(cc);
		
	} else {
		quadratureFactor = 1.0/TUTMAM_SQRTPI;
	}
	
	if (k < 0.1) { /* k=0 */
		return 0.0;
	}
	
	if (k > 0.9) { /* k=1 */
		massFraction = C_YI_PARTICLE(c,t,iUds-iUdsM1First,jOther);
	}
	
	/* integration starts from zero */
	integral1 = 0.0;

	for (iIntegral1 = 0; iIntegral1 < gaussHermiteLevel; ++iIntegral1) {
		x1 = gauss_hermite_abscissas(iIntegral1); /* x-variable is defined by Gauss-Hermite quadrature */
		dp1 = cmdj*exp(x1*sqrt(2*ln2sj)); /* converting x back to dp */
		diffP1 = diffIndep*diff_dep(c,t,dp1);
				
		integral2 = 0.0;
		
		if (jOther == powerLawDistribution) {
			for (iIntegral2 = 0; iIntegral2 < 4; ++iIntegral2) {
				dp2 = gauss_olin_abscissas_dp(iIntegral2,cmdjOther);
				diffP2 = diffIndep*diff_dep(c,t,dp2);
				
				if (transitionRegimeCorrectionFactorForCoagulationLaw == 1) {
					correctionFactor = fuchs_sutugin_for_coagulation(dp1,dp2,diffP1,diffP2,tempFluid,rhoPj,rhoPjOther);
				} else if (transitionRegimeCorrectionFactorForCoagulationLaw == 2) {
					correctionFactor = dahneke_for_coagulation(dp1,dp2,diffP1,diffP2,tempFluid,rhoPj,rhoPjOther);
				} else {
					Error("Illegal transitionRegimeCorrectionFactorForCoagulationLaw\n");
					return 0.0;
				}
				
				coagCoeff = coagulation_coefficient(dp1,dp2,diffP1,diffP2,correctionFactor);
				
				if (k > 0.9) { /* k=1 */
					integral2 += gauss_olin_weights(iIntegral2,cc)*coagCoeff*CBC(dp1);
				
				} else { /* k=2/3 */
					integral2 += gauss_olin_weights(iIntegral2,cc)*coagCoeff*(pow(rhoPj*CBC(dp1)+rhoPjOther*CBC(dp2),TUTMAM_23) - pow(rhoPj,TUTMAM_23)*SQR(dp1));
				}
			}
			
		} else {
			for (iIntegral2 = 0; iIntegral2 < gaussHermiteLevel; ++iIntegral2) {
				x2 = gauss_hermite_abscissas(iIntegral2); /* x-variable is defined by Gauss-Hermite quadrature */
				dp2 = cmdjOther*exp(x2*sqrt(2.*ln2sjOther)); /* converting x back to dp */
				diffP2 = diffIndep*diff_dep(c,t,dp2);
				
				if (transitionRegimeCorrectionFactorForCoagulationLaw == 1) {
					correctionFactor = fuchs_sutugin_for_coagulation(dp1,dp2,diffP1,diffP2,tempFluid,rhoPj,rhoPjOther);
				} else if (transitionRegimeCorrectionFactorForCoagulationLaw == 2) {
					correctionFactor = dahneke_for_coagulation(dp1,dp2,diffP1,diffP2,tempFluid,rhoPj,rhoPjOther);
				} else {
					Error("Illegal transitionRegimeCorrectionFactorForCoagulationLaw\n");
					return 0.0;
				}
				
				coagCoeff = coagulation_coefficient(dp1,dp2,diffP1,diffP2,correctionFactor);
				
				if (k > 0.9) { /* k=1 */
					integral2 += gauss_hermite_weights(iIntegral2)*coagCoeff*CBC(dp1);
				
				} else { /* k=2/3 */
					integral2 += gauss_hermite_weights(iIntegral2)*coagCoeff*(pow(rhoPj*CBC(dp1)+rhoPjOther*CBC(dp2),TUTMAM_23) - pow(rhoPj,TUTMAM_23)*SQR(dp1));
				}
			}
		}

		integral1 += integral2*gauss_hermite_weights(iIntegral1);
	}

	sourceTerm = quadratureFactor*integral1;

	if (k > 0.9) {
		sourceTerm *= TUTMAM_PI6*rhoPjOther*massFraction;
		
	} else {
		sourceTerm *= TUTMAM_PI623;
	}
	
	return sourceTerm;
}

/* Source term macro for intermodal coagulation */
DEFINE_SOURCE(inter_coagulation,c,t,dS,eqn)
{	
	int iUds;				/* UDS ID number */
	int iUdsReducedToMode0;	/* UDS ID number reduced to mode 0 */
	int j,jOther;			/* mode ID */
	real ntotj,ntotjOther;	/* ntot (1/m^3) */
	real quadratureFactor;
	real sourceTerm = 0.0;
	dS[eqn] = 0.0;			/* differential of the source term */
	
	iUds = eqn-EQ_UDS;
	iUdsReducedToMode0 = iUds % nUdsPerMode;
	j = iUds/nUdsPerMode;
	
	ntotj = C_NTOT(c,t,j);
	if (ntotj < minNumberConc) {
		return 0.0;
	}
	
	if (coagulationProcess == 1 && interModalCoagulationProcess == 1) {
		for (jOther = 0; jOther < nTutmamModes; ++jOther) {
			if (j != jOther) { /* intermodal coagulation */
				if (coagulationMatrix[j][jOther] == 1) {
					ntotjOther = C_NTOT(c,t,jOther);
					if (ntotjOther < minNumberConc) {
						if (j > jOther) { /* loss and gain for mode j */
							sourceTerm += inter_coagulation_rate_LG(c,t,iUdsReducedToMode0,j,jOther)*ntotjOther;
							
						} else { /* loss only for mode j */
							sourceTerm -= inter_coagulation_rate_L(c,t,iUdsReducedToMode0,j,jOther)*ntotjOther;
						}
					}
				}
			}
		}
	}
	
	if (j == powerLawDistribution) {
		quadratureFactor = quadrature_factor_olin(C_LN2S(c,t,j)*log(C_CMD(c,t,j)/powerLawD1));
		
	} else {
		quadratureFactor = 1.0/TUTMAM_SQRTPI;
	}
	
	return TUTMAM_2PI*quadratureFactor*sourceTerm*ntotj;
}


