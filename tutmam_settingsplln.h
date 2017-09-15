/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	Header-file for tutmam_settings.c
*	You need to change the define-variables in this file,
*	if you want to change these settings.
*/
#include  "udf.h"

#ifndef _TUTMAM_SETTINGS_H
#define _TUTMAM_SETTINGS_H

	/* Number of particle modes. */
	/* You must change it to correspond your case. */
	#define nTutmamModes 2
	
	/* Number of particle phases. */
	/* You must change it to correspond your case. */
	#define nTutmamPhases 1

	/* Number of particle species. */
	/* You must change it to correspond your case. */
	#define nTutmamSpecies 2
	
	/* Number of Fluent species. */
	/* You must change it to correspond your case. */
	#define nFluentSpecies 3
	
	/* Molar masses of particle species (kg/mol). */
	/* You must change it to correspond your case. */
	#define molarMasses 0.098079,0.018015
	
	/* Change these to ID values of particle species in Fluent. */
	#define fluentSpeciesIds 0,1
	
	/* Change these to particle phase IDs. */
	#define phaseIds 0,0
	
	/* Particle mass fraction vector used in initialization of the domain. */
	/* You must change it to correspond your case. */
	/* Separate vectors for different modes with , */
	#define initialYSpeciesP {0.3,0.7},{0.3,0.7}
	
	/* ID of power law distribution */
	#define powerLawDistribution 0
	
	/* UDM id number for which number concentration (#/m^3) is saved */
	#define iUdmNumberConc 0
	
	/* UDM id number for which count median diameter (m) is saved, log-normal. */
	/* D2 (m), power law. */
	#define iUdmCmd 1
	
	/* UDM id number for which ln^2(GSD) is saved, log-normal. */
	/* Alpha, power law. */
	#define iUdmLn2s 2
	
	/* UDM id number for which particle density (kg/m^3) is saved */
	#define iUdmParticleDensity 3
	
	/* UDM id number for which particle temperature (K) is saved */
	#define iUdmParticleTemperature 4
	
	/* UDM for water condensation factor, kappa */
	#define iWaterKappaUdm 5
	
	/* UDM for rheq/rh */
	#define iRhEqPerRhUdm 6
	
	/* UDM id numbers for condensation source terms */
	#define iUdmCondensationM23 7
	#define iUdmCondensationM1First 8
	
	/* UDM id numbers for intramodal coagulation source terms */
	#define iUdmCoagulationM0 (iUdmCondensationM1First+nTutmamSpecies)
	#define iUdmCoagulationM23 (iUdmCondensationM1First+nTutmamSpecies+1)
	
	/* The last UDM +1 for mode 0 */
	#define nUdmPerMode (iUdmCoagulationM23+1)
	
	/* UDM id numbers for nucleation source terms */
	#define iUdmNucleationM0 (nTutmamModes*nUdmPerMode)
	#define iUdmNucleationM23 (iUdmNucleationM0+1)
	#define iUdmNucleationM1First (iUdmNucleationM23+1)
	
	/* UDM id numbers for PL to LN transfer source terms */
	#define iUdmTransferPl2LnM0 (iUdmNucleationM1First+nTutmamSpecies)
	#define iUdmTransferPl2LnM23 (iUdmTransferPl2LnM0+1)
	#define iUdmTransferPl2LnM1 (iUdmTransferPl2LnM23+1)
	
	/* UDS id number for number concentration scalar (1/kg) */
	#define iUdsM0 0
	
	/* UDS id number for surface area concentration scalar (kg^(-1/3)) */
	#define iUdsM23 1
	
	/* UDS id number for the first species mass concentration scalar (dimensionless). */
	/* The id numbers for the following species mass concentration scalars follows this number. */
	#define iUdsM1First 2
	
	/* The last UDS +1 for mode 0 */
	#define nUdsPerMode (iUdsM1First+nTutmamSpecies)
	
	/* These minimum values are used to prevent division by zero in some functions. */
	/* The default values are very small for most cases; therefore, their altering is not needed. */
	#define minNumberConc 1.0
	#define minSurfaceConc 1.0e-16
	#define minMassConc 1.0e-25
	#define minCMD 1.15e-9
	#define maxCMD 100.0e-9
	
	/* The distribution values used in initialization of the domain. */
	/* You should change these to correspond your case. */
	/* Number concentration (#/m^3) */
	/* Separate values for different modes with , */
	#define initialNumberConc 0.0e10,0.0e10
	/* Count median diameter (m) */
	#define initialCMD 1.5e-9,2.0e-9
	/* GSD */
	#define initialGSD 1.1,1.2
	
	/* Particle density (kg/m^3) used in initialization of the domain. */
	/* You should change it to correspond your case. */
	/* Separate values for different modes with , */
	#define initialRhoParticle 1200.0,1200.0
	
	/* Maximum particle temperature difference (K) */
	/* You can change it to correspond your case. */
	#define maxParticleTemperatureDifference 100.0
	
	/* Options for water equilibrium calculation */
	#define minCMDForKappaCalculation 1.5e-9
	#define convergenceCriteriumForKappaCalculation 0.01
	#define maxIterationsInKappaCalculation 20
	#define uRFKappaCalculation 0.5
	#define uRFKappa 0.8
		
	/* the value when alpha is considered 0 */
	#define whenAlphaIsZero 0.001
	
	/* min and max values for power law distribution parameters */
	#define minAlpha -5.0
	#define maxAlpha 5.0
	#define powerLawD1 1.15e-9
	#define maxD2 100.0e-9
	
	/* under-relaxation factor for power law distribution parameters */
	#define uRFPowerLawParameters 0.99
	
	/* coagulation processes from mode to another mode */
	/* Separate vectors for different modes with , */
	#define coagulationFromModeToMode {1,1},{1,1}
	
	/* Wall condensation laws for vapors. */
	/* No condensation (0), condensation only if saturation exceeded (1), activity based condensation (2), full condensation(3) */
	#define wallCondensationLaws 2,2
	
	extern real molarMassVector[nTutmamSpecies];			/* a vector storing molar masses of particle species (kg/mol) */
	extern int fluentSpeciesIdVector[nTutmamSpecies];		/* a vector storing IDs of particle species in Fluent */
	extern int tutmamSpeciesIdVector[nFluentSpecies];		/* a vector storing particle IDs of Fluent species */
	extern int phaseIdVector[nTutmamSpecies];				/* a vector storing phase IDs of particle species */
	extern int phaseMatrix[nTutmamPhases][nTutmamSpecies];	/* a matrix storing particle species flags for specific phases */
	extern int coagulationMatrix[nTutmamModes][nTutmamModes];	/* a matrix storing coagulation flags for specific modes */
	extern int wallCondensationLawVector[nTutmamSpecies];	/* a vector storing wall condensation laws */
		
	extern int fluidZoneNumber;					/* cell zone id of the fluid */
	extern int diffusionLaw;					/* Is parametrisation used for slip correction coefficient? */
	extern int slipCorrectionLaw;				/* Is slip correction calculated as size-dependent? */
	extern real constantSlipCorrectionFactor;	/* value for constant slip correction coefficient */
	extern int particleTurbDiff;				/* Is particle turbulent diffusion included? */
	extern real particleTurbSchmidtNumber;		/* turbulent Schmidt number for particles */
	
	extern int gaussHermiteLevel;				/* Gauss-Hermite quadrature level (5, 7, 9, or 11) */
	extern int condensationIntegrationLevel;			/* integration level for condensation */
	extern int condensationIntegrationBins;			/* integration bins for condensation */
	extern int coagulationIntegrationLevel;			/* integration level for coagulation */
	extern int coagulationIntegrationBins;			/* integration bins for coagulation */

	extern real minGSD;							/* Minimum value for GSD. 1 causes division by zero; therefore, over 1 is required. */
	extern real maxGSD;							/* Maximum value for GSD. */
	extern int powerLawParametersSolvingMethod;	/* Method for solving power law distribution parameters. 1: interpolation, 2: iteration */
	
	extern real uRFNucleation;					/* under-relaxation factor for nucleation */
	extern int nucleationLaw;					/* nucleation law ID */
	extern real nucleationCorrectionFactor;		/* correction factor for nucleation */
	extern int iNucleatingSpecies[nTutmamSpecies];				/* vector for nucleating species in Olin nucleation law */
	extern real nucleationExponents[nTutmamSpecies];			/* vector for nucleation exponents in Olin nucleation law */
	extern real nucleationSatVapPresExponents[nTutmamSpecies];	/* vector for nucleation saturation vapor pressure exponents in Olin nucleation law */
	extern real nMolecClusterVector[nTutmamSpecies];			/* vector for number of molecules in a cluster formed by Olin nucleation law */
	extern real clusterGSD;						/* GSD of the distribution formed by Olin nucleation law */
	
	extern real uRFCondensationM23;				/* under-relaxation factor for condensation to moment 2/3 */
	extern real uRFCondensationM1[nTutmamSpecies];/* under-relaxation factor vector for condensation to moment 1 */
	extern int latentHeatOfCondensation;		/* Is latent heat of condensation considered? */
	extern int iCondensingSpecies[nTutmamModes*nTutmamSpecies];	/* vector for condensing species flags */
	extern int kelvinEffect;					/* Is Kelvin effect on ? */
	extern int phaseActivityModel;				/* Is phase activity model on ? */
	extern int condensationDirectionVector[nTutmamSpecies]; /* Condensation directions: evaporation (-1) both (0) condensation (1) */
	extern int waterEq;							/* Is water condensation calculated through equilibrium? */
	extern int waterEqSpecies;					/* Connecting species in water equilibrium calculation */
	extern real minKappa;						/* Minimum kappa value */
	extern real maxKappa;						/* Maximum kappa value */
	extern int dieselExhaustHCFractionModel;	/* Is diesel exhaust hydrocarbon fraction model on ? */
	extern int dieselExhaustHC;					/* Diesel exhaust hydrocarbon fraction model species ID */
	extern real interModalCondensationFactor;	/* Intermodal condensation factor used with combined power law and log-normal distribution */

	extern int diffusionProcess;				/* Is diffusion process on? */
	extern int nucleationProcess;				/* Is nucleation process on? */
	extern int nucleationComputing;				/* Is nucleation computed? */
	extern int condensationProcess;				/* Is condensation process on? */
	extern int condensationComputing;			/* Is condensation computed? */
	extern int coagulationProcess;				/* Is coagulation process on? */
	extern int coagulationComputing;			/* Is coagulation computed? */
	extern int coagulationIterationSkip;		/* How often coagulation is computed? */
	extern int selfCoagulationalTransferProcess;/* Is self-coagulational transfer process on? */
	extern int interModalCondensationProcess;   /* Is intermodal condensation process on? */
	
	extern int intraModalCoagulationProcess;	/* Is intramodal coagulation process on? */
	extern int interModalCoagulationProcess;	/* Is intramodal coagulation process on? */
	extern real uRFCoagulation;					/* under-relaxation factor for coagulation */
	extern int transitionRegimeCorrectionFactorForCoagulationLaw; /* transition regime correction factor law for coagulation */
	extern int coagulationRobustnessModel;		/* Is robustness model used in coagulation? */
	
	extern real uRFTransferPl2Ln;				/* under-relaxation factor for PL to LN transfer */
	
	
	extern int lippu;

#endif
