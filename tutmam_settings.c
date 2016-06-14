/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	This file initializes the variables that are set using CFD-TUTMAM GUI.
*	Do not change anything.
*	These values have no meaning, because CFD-TUTMAM GUI overwrites them with default/saved values. 
*/
#include  "udf.h"
#include  "tutmam_settings.h"

/* However, these are not set by CFD-TUTMAM GUI */
real molarMassVector[] = { molarMasses };				/* a vector storing molar masses of particle species (kg/mol) */
int fluentSpeciesIdVector[] = { fluentSpeciesIds };		/* a vector storing IDs of particle species in Fluent */
int tutmamSpeciesIdVector[] = { 0 };					/* a vector storing particle IDs of Fluent species */
int phaseIdVector[] = { phaseIds };						/* a vector storing phase IDs of particle species */
int phaseMatrix[nTutmamPhases][nTutmamSpecies] = {0};	/* a matrix storing particle species flags for specific phases */
int coagulationMatrix[nTutmamModes][nTutmamModes] = { coagulationFromModeToMode };	/* a matrix storing coagulation flags for specific modes */
int wallCondensationLawVector[nTutmamSpecies] = { wallCondensationLaws };	/* a vector storing wall condensation laws */;	/* a vector storing wall condensation laws */

/* The following values are overwritten by CFD-TUTMAM GUI */
int fluidZoneNumber = 1; 				/* cell zone id of the fluid */
int gaussHermiteLevel = 7;				/* Gauss-Hermite quadrature level (5, 7, 9, or 11) */
int condensationIntegrationLevel = 1;			/* integration level for condensation */
int condensationIntegrationBins = 200;			/* integration bins for condensation */
int coagulationIntegrationLevel = 1;			/* integration level for coagulation */
int coagulationIntegrationBins = 20;			/* integration bins for coagulation */
int diffusionLaw = 1;					/* Is parametrisation used for slip correction coefficient? */
int slipCorrectionLaw = 1;				/* Is slip correction calculated as size-dependent? */
real constantSlipCorrectionFactor = 1.0;/* value for constant slip correction coefficient */
int particleTurbDiff = 1;				/* Is particle turbulent diffusion included? */
real particleTurbSchmidtNumber = 0.7;	/* turbulent Schmidt number for particles */

real minGSD = 1.01;						/* Minimum value for GSD. 1 causes division by zero; therefore, over 1 is required. */
real maxGSD = 3.0;						/* Maximum value for GSD. */

real uRFCondensationM23 = 1.0;			/* under-relaxation factor for condensation to moment 2/3 */
real uRFCondensationM1[] = { 0.0 };		/* under-relaxation factor vector for condensation to moment 1 */

real uRFNucleation = 0.9;				/* under-relaxation factor for nucleation */
int nucleationLaw = 1;					/* nucleation law ID */
real nucleationCorrectionFactor = 1.0;	/* correction factor for nucleation */
int iNucleatingSpecies[] = { 0 };				/* vector for nucleating species in Olin nucleation law */
real nucleationExponents[] = { 0.0 };			/* vector for nucleation exponents in Olin nucleation law */
real nucleationSatVapPresExponents[] = { 0.0 };	/* vector for nucleation saturation vapor pressure exponents in Olin nucleation law */
real nMolecClusterVector[] = { 0.0 };			/* vector for number of molecules in a cluster formed by Olin nucleation law */

int latentHeatOfCondensation = 1;		/* Is latent heat of condensation considered? */
int iCondensingSpecies[] = { 0 };		/* vector for condensing species flags */
int kelvinEffect = 0;					/* Is Kelvin effect on ? */
int phaseActivityModel = 0;				/* Is phase activity model on ? */
int condensationDirectionVector[nTutmamSpecies] = { 0 }; /* Condensation directions: evaporation (-1) both (0) condensation (1) */
int waterEq = 0;						/* Is water condensation calculated through equilibrium? */
int waterEqSpecies = 0;					/* Connecting species in water equilibrium calculation */
real minKappa = 0.01;					/* Minimum kappa value */
real maxKappa = 5.0;					/* Maximum kappa value */
int dieselExhaustHCFractionModel = 0;	/* Is diesel exhaust hydrocarbon fraction model on ? */
int dieselExhaustHC = 0;				/* Diesel exhaust hydrocarbon fraction model species ID */
real interModalCondensationFactor = 0.0;/* Intermodal condensation factor used with combined power law and log-normal distribution */

int diffusionProcess = 1;				/* Is diffusion process on? */
int nucleationProcess = 0;				/* Is nucleation process on? */
int nucleationComputing = 1;			/* Is nucleation computed? */
int condensationProcess = 0;			/* Is condensation process on? */
int condensationComputing = 1;			/* Is condensation computed? */
int coagulationProcess = 0;				/* Is coagulation process on? */
int coagulationComputing = 1;			/* Is coagulation computed? */
int coagulationIterationSkip = 1;		/* How often coagulation is computed? */
int selfCoagulationalTransferProcess = 0;/* Is self-coagulational transfer process on? */
int interModalCondensationProcess = 0;  /* Is intermodal condensation process on? */

int intraModalCoagulationProcess = 1;	/* Is intramodal coagulation process on? */
int interModalCoagulationProcess = 1;	/* Is intermodal coagulation process on? */
real uRFCoagulation = 0.5;				/* under-relaxation factor for coagulation */
int transitionRegimeCorrectionFactorForCoagulationLaw = 2; /* transition regime correction factor law for coagulation */
int coagulationRobustnessModel = 0;		/* Is robustness model used in coagulation? */

int lippu = 2;