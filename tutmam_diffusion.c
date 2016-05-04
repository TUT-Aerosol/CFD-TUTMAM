/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	This file defines Brownian diffusion of particles
*	You should not need to change this
*/
#include  "udf.h"
#include  "tutmam_constants.h"
#include  "tutmam_settings.h"
#include  "tutmam_func.h"
#include  "tutmam_diffusion.h"
#include  "tutmam_material_func.h"

/* Diffusion coefficient: size-inpedependent part (m^3/s) */
real diff_indep(cell_t c, Thread *t) {
	real temp; /* temperature (K) */
	real visc; /* laminar dynamic viscosity (Pa*s) */
	
	temp = C_T(c,t);
	visc = C_MU_L(c,t);
	
	/* kT/(3 pi mu) */
	return TUTMAM_BOLTZ*temp/(TUTMAM_3PI*visc);
}

/* Diffusion coefficient: size-dependent part (1/m) */
real diff_dep(cell_t c, Thread *t, real dp) {
	real temp; /* temperature (K) */
	real slipCorrectionCoefficient; /* Cc */
	real pressure; /* total pressure (Pa) */
	
	temp = C_T(c,t);
	pressure = C_P_TOT(c,t);
	slipCorrectionCoefficient = slip_correction_coefficient(temp,dp,pressure);
	
	/* Cc/dp */
	return slipCorrectionCoefficient/dp;
}

/* k-th moment averaged diffusion coefficient */
/* integrated numerically by Gauss-Hermite quadrature method */
real diffusion_coefficient_uds(cell_t c, Thread *t, int i) {
	real rhoFluid;				/* fluid density (kg/m^3) */
	real k;						/* k-th power of the moment */
	real diffIndep;				/* size-independent part of particle diffusion coefficient (m^3/s) */
	real diffDep;				/* size-dependent part of particle diffusion coefficient (1/m) */
	real diffCoeff;				/* diffusion coefficient (kg/(m s)) */
	real integral;				/* value of integrated part of the function (m^3k * 1/m) */
	int iIntegral;				/* index used in integration */
	real cmd;					/* count median diameter (m) */
	real ln2s;					/* ln^2(GSD) */
	real d;						/* D2/D1 for power law distribution */
	real dp;					/* current particle diameter used in integration (m) */
	real x;						/* (ln(dp/cmd)) / (sqrt(2)*ln(GSD)) */
	real cc;
	int j = i/nUdsPerMode;		/* mode ID */
	
	rhoFluid = C_R(c,t);
	cmd = C_CMD(c,t,j);
	ln2s = C_LN2S(c,t,j);
	diffIndep = diff_indep(c,t);
	k = k_moment(i); /* k is defined by UDS number i */
	
	/* integration starts from zero */
	integral = 0.0;
	
	if (j == powerLawDistribution) {
		d = cmd/powerLawD1;
		if (d < 1.001) {
			return rhoFluid*diffIndep*diff_dep(c,t,cmd);
		}
		
		for (iIntegral = 0; iIntegral < 4; ++iIntegral) {
			dp = gauss_olin_abscissas_dp(iIntegral,cmd);
			diffDep = diff_dep(c,t,dp);
			cc = ln2s*log(d);
			
			if (k > 0.1) { /* when k is not 0 */
				/* weights are defined by Gauss-Hermite quadrature */
				integral += gauss_olin_weights(iIntegral,cc)*pow(dp,3.*k)*diffDep; /* (m^3k * 1/m) */
				
			} else { /* when k is 0, shorter expression is enough */
				integral += gauss_olin_weights(iIntegral,cc)*diffDep; /* (1/m) */
			}
		}
				
		if (fabs(ln2s+3.0*k) < whenAlphaIsZero) {
			diffCoeff = rhoFluid*diffIndep*pow(powerLawD1,-3.0*k)*integral;
			
		} else {
			if (k > 0.1) { /* when k is not 0 */
				diffCoeff = rhoFluid*diffIndep*(ln2s+3.0*k)*pow(powerLawD1,-3.0*k)*log(d)*integral/(pow(d,ln2s+3.0*k)-1.0);
				
			} else { /* when k is 0, shorter expression is enough */
				diffCoeff = rhoFluid*diffIndep*ln2s*log(d)*integral/(pow(d,ln2s)-1.0);
			}
		}
		
		return tutmam_lower_limit(diffCoeff,1.0e-20);
	}
	
	/* integration is done in amount of gaussHermiteLevel steps */
	/* gaussHermiteLevel comes from the settings in CFD-TUTMAM GUI */
	for (iIntegral = 0; iIntegral < gaussHermiteLevel; ++iIntegral) {
		x = gauss_hermite_abscissas(iIntegral); /* x-variable is defined by Gauss-Hermite quadrature */
		dp = cmd*exp(x*sqrt(2.*ln2s)); /* converting x back to dp */
		diffDep = diff_dep(c,t,dp);
		
		if (k > 0.1) { /* when k is not 0 */
			/* weights are defined by Gauss-Hermite quadrature */
			integral += gauss_hermite_weights(iIntegral)*pow(dp,3.*k)*diffDep; /* (m^3k * 1/m) */
			
		} else { /* when k is 0, shorter expression is enough */
			integral += gauss_hermite_weights(iIntegral)*diffDep; /* (1/m) */
		}
	}
	
	/* k-th moment diffusion coefficient is divided by k-th moment to get the averaged value. */
	/* k-th moments are obtained from analytical solution of the distribution density integral. */
	/* Multiplied by size-independent variables. */
	if (k > 0.1) { /* when k is not 0 */
		return rhoFluid*diffIndep*pow(cmd,-3.*k)*exp(-4.5*k*k*ln2s)*integral/TUTMAM_SQRTPI; /* (kg/(m s)) */
		
	} else { /* when k is 0, shorter expression is enough */
		return rhoFluid*diffIndep*integral/TUTMAM_SQRTPI; /* (kg/(m s)) */
	}
}

/* k-th moment averaged diffusion coefficient */
/* calculated by analytical solution of integrals */
/* using polynomial-parameterized slip correction coefficient */
real diffusion_coefficient_uds_parametrisation(cell_t c, Thread *t, int i) {
	real rhoFluid;					/* fluid density (kg/m^3) */
	real temp; 						/* temperature (K) */
	real k;							/* k-th power of the moment */
	real diffIndep;					/* size-independent part of particle diffusion coefficient (m^3/s) */
	real sum;						/* value of summation (1/m) */
	real parametrisationFactors[4];	/* factors in parametrisation */
	real cmd;						/* count median diameter (m) */
	real ln2s;						/* ln^2(GSD) */
	int iSum;						/* index used in summation */
	int m;							/* power of dp in parametrisation */
	int j = i/nUdsPerMode;			/* mode ID */
	
	rhoFluid = C_R(c,t);
	
	/* temperature is limited to range where the parametrisation is valid (250-1800 K, error<6%) */
	temp = tutmam_limits(250.0,C_T(c,t),1800.0);

	cmd = C_CMD(c,t,j);
	ln2s = C_LN2S(c,t,j);
	diffIndep = diff_indep(c,t);
	k = k_moment(i); /* k is defined by UDS number i */
	
	parametrisationFactors[0] = -1.15e11;				/* (1/m^2) */
	parametrisationFactors[1] = -25.61*temp+3.447e5;	/* (1/m) */
	parametrisationFactors[2] = -1.451e-4*temp+0.7301;	/* (1) */
	parametrisationFactors[3] = 1.044e-9*temp-8.884e-8;	/* (m) */
	
	/* summation starts from zero */
	sum = 0.0;
	
	/* summed 4 times */
	for (iSum = 0; iSum < 4; ++iSum) {
		m = 1-iSum; /* this comes from the parametrisation */
		
		if (j == powerLawDistribution) {
			sum += parametrisationFactors[iSum]*pow(powerLawD1,m)*momentAveragingPowerLaw(ln2s,cmd/powerLawD1,k,m); /* (1/m) */
		
		} else {
			sum += parametrisationFactors[iSum]*pow(cmd,m)*exp((0.5*m+k)*m*ln2s); /* (1/m) */
		}
	}
	
	/* multiplied by size-independent variables */
	return rhoFluid*diffIndep*sum; /* (kg/(m s)) */
}

/* k-th moment averaged diffusion coefficient */
/* calculated by analytical solution of integrals */
/* using constant slip correction coefficient */
real diffusion_coefficient_uds_constant_slip_correction(cell_t c, Thread *t, int i) {
	real rhoFluid;			/* fluid density (kg/m^3) */
	real k;					/* k-th power of the moment */
	real diffIndep;			/* size-independent part of particle diffusion coefficient (m^3/s) */
	real cmd;				/* count median diameter (m) */
	real ln2s;				/* ln^2(GSD) */
	int j = i/nUdsPerMode;	/* mode ID */
	
	rhoFluid = C_R(c,t);
	cmd = C_CMD(c,t,j);
	ln2s = C_LN2S(c,t,j);
	diffIndep = diff_indep(c,t);
	k = k_moment(i); /* k is defined by UDS number i */
	
	if (j == powerLawDistribution) {
		return rhoFluid*diffIndep*constantSlipCorrectionFactor/powerLawD1*momentAveragingPowerLaw(ln2s,cmd/powerLawD1,k,-1.0);
		
	}
	
	return rhoFluid*diffIndep*constantSlipCorrectionFactor/cmd*exp((0.5-k)*ln2s); /* (kg/(m s)) */
}

/* This is the macro for particle diffusion coefficient, which can be seen in Fluent GUI. */
DEFINE_DIFFUSIVITY(diff_uds,c,t,i)
{	
	real diffCoeff; /* Effective diffusion coefficient (kg/(m s)) */
	
	/* If diffusion process is off */
	if (diffusionProcess == 0) {
		/* Returning very small number (it cannot be zero) */
		return 1.0e-20;
	}
	
	/* Firstly, laminar diffusion coefficient is calculated */
	if (slipCorrectionLaw == 0) { /* this comes from the settings in TUTMAM GUI */
		/* using constant slip correction coefficient */
		diffCoeff = diffusion_coefficient_uds_constant_slip_correction(c,t,i);
		
	} else { 
		/* size-dependent slip correction coefficient requires tougher integration: */
		if (diffusionLaw == 0) { /* this comes from the settings in TUTMAM GUI */
			/* numerical integration using Gauss-Hermite quadrature */
			diffCoeff = diffusion_coefficient_uds(c,t,i);
			
		} else { 
			/* analytical integration using parameterised slip correction coefficient */
			diffCoeff = diffusion_coefficient_uds_parametrisation(c,t,i);
		}
	}
	
	/* Secondly, turbulent diffusion coefficient is added to the laminar one */
	if (particleTurbDiff == 1) { /* this comes from the settings in TUTMAM GUI */
		diffCoeff += C_MU_T(c,t)/particleTurbSchmidtNumber; /* mu_t/Sc_t (kg/(m s)) */
	}

	return diffCoeff; /* (kg/(m s)) */
}

/* This is the macro for gas laminar diffusion coefficient, which can be seen in Fluent GUI. */
DEFINE_DIFFUSIVITY(diff_gas,c,t,iFluentSpecies)
{	
	real temp;	 	/* temperature (K) */
	real pressure;	/* total pressure (Pa) */
	
	temp = C_T(c,t);
	pressure = C_P_TOT(c,t);
	
	return diffusion_coefficient_gas(temp,pressure,iFluentSpecies); /* (m^2/s) */
}

