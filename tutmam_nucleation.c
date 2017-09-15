/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	This file includes particles nucleation functions
*	You should not need to change this
*/
#include  "udf.h"
#include  "tutmam_constants.h"
#include  "tutmam_settings.h"
#include  "tutmam_func.h"
#include  "tutmam_nucleation.h"
#include  "tutmam_material_func.h"

/* Sulfuric acid-water nucleation rate from Vehkamaki et al. (2002,2003) parametrisation */
void nucleation_rate_vehkamaki(real *sourceTerm, cell_t c, Thread *t) {
	int iFluentSpeciesWater;				/* ID for water vapor in Fluent */
	int iFluentSpeciesSulfuricAcid;			/* ID for sulfuric acid vapor in Fluent */
	int iSpeciesWater;						/* ID for water in TUTMAM */
	int iSpeciesSulfuricAcid;				/* ID for sulfuric acid in TUTMAM */
	real temp;								/* temperature (K) */
	real temp2;								/* temperature^2 */
	real temp3;								/* temperature^3 */
	real pressure;							/* total pressure (Pa) */
	real moleFractionWater;					/* mole fraction of water in fluid */
	real moleFractionSulfuricAcid;			/* mole fraction of sulfuric acid in fluid */
	real satVapPresWater;					/* saturation vapor pressure of water (Pa) */
	real numberConcSulfuricAcid;			/* molecular number concentration of sulfuric acid in fluid (1/cm^3) */
	real lnNumberConcSulfuricAcid;			/* ln numberConcSulfuricAcid */
	real J;									/* nucleation rate (1/(m^3 s)) */
	real rh;								/* relative humidity */
	real lnRh;								/* ln(rh) */
	real lnRh2;								/* ln^2(rh) */
	real lnRh3;								/* ln^3(rh) */
	real moleFractionSulfuricAcidCluster;	/* mole fraction of sulfuric acid in the cluster */
	real moleFractionSulfuricAcidCluster2;	/* moleFractionSulfuricAcidCluster^2 */
	real nMolecCluster;						/* number of molecules in the cluster */
	real mWaterCluster;						/* partial water mass in the cluster (kg) */
	real mSulfuricAcidCluster;				/* partial sulfuric acid mass in the cluster (kg) */
	real sCluster;							/* 0.667th moment of the cluster (kg^(2/3)) */
	int parametrisationRange = 0;			/* Vehkamaki parametrisation: 0 (2002), 1 (2003) */
	real parameters1[10] = { 0.0 };			/* parameters used in Vehkamaki parametrisation */
	real parameters2[10] = { 0.0 };			/* parameters used in Vehkamaki parametrisation */
	
	iFluentSpeciesWater = SV_SpeciesIndex("h2o"); /* find water ID */
	iFluentSpeciesSulfuricAcid = SV_SpeciesIndex("h2so4"); /* find sulfuric acid ID */
	
	if (iFluentSpeciesWater < 0 || iFluentSpeciesSulfuricAcid < 0) { /* if both not found */
		Error("h2o and h2so4 are both needed for nucleationLaw 1\n");
		return;
	}
	
	iSpeciesWater = tutmamSpeciesIdVector[iFluentSpeciesWater];
	iSpeciesSulfuricAcid = tutmamSpeciesIdVector[iFluentSpeciesSulfuricAcid];
	
	temp = C_T(c,t);
	pressure = C_P_TOT(c,t);
	moleFractionSulfuricAcid = C_XI(c,t,iSpeciesSulfuricAcid);
	numberConcSulfuricAcid = 1.0e-6*pressure*moleFractionSulfuricAcid/(TUTMAM_BOLTZ*temp);
	
	/* limits of parametrisations */
	if (temp < 230.15) {
		return;
	} else if (temp > 400.15) {
		return;
	}
	
	/* choosing which parametrisation is to be used, 10 denotes that both are useful */
	if (temp > 305.15) {
		parametrisationRange = 1;
	} else {
		if (temp > 300.15) {
			parametrisationRange = 10;
		}
	}

	if (parametrisationRange == 0) {
		if (numberConcSulfuricAcid < 1.0e4) {
			return;
		} else if (numberConcSulfuricAcid > 1.0e11) {
			numberConcSulfuricAcid = 1.0e11;
		}
	}
	
	if (parametrisationRange == 1) {
		if (numberConcSulfuricAcid < 2.0e9) {
			return;
		} else if (numberConcSulfuricAcid > 2.0e15) {
			numberConcSulfuricAcid = 2.0e15;
		}
	}
	
	if (parametrisationRange == 10) {
		if (numberConcSulfuricAcid < 1.0e4) {
			return;
		} else if (numberConcSulfuricAcid > 2.0e15) {
			numberConcSulfuricAcid = 2.0e15;
		}
		
		if (numberConcSulfuricAcid < 2.0e9) {
			parametrisationRange = 0;
		} else if (numberConcSulfuricAcid > 1.0e11) {
			parametrisationRange = 1;
		}
	}
	
	satVapPresWater = saturation_vapor_pressure(temp,iSpeciesWater);
	moleFractionWater = C_XI(c,t,iSpeciesWater);
	rh = moleFractionWater*pressure/satVapPresWater;
	
	if (parametrisationRange == 0) {
		if (rh < 0.0001) {
			return;
		} else if (rh > 1.0) {
			rh = 1.0;
		}
	}
	
	if (parametrisationRange == 1) {
		if (rh < 0.01) {
			return;
		} else if (rh > 1.0) {
			rh = 1.0;
		}
	}
	
	if (parametrisationRange == 10) {
		if (rh < 0.0001) {
			return;
		} else if (rh > 1.0) {
			rh = 1.0;
		}
		
		if (rh < 0.01) {
			parametrisationRange = 0;
		} else {
			parametrisationRange = 1;
		}
	}
	
	temp2 = SQR(temp);
	temp3 = pow(temp,3.0);
	lnRh = log(rh);
	lnRh2 = SQR(lnRh);
	lnRh3 = pow(lnRh,3.0);
	lnNumberConcSulfuricAcid = log(numberConcSulfuricAcid);
	
	if (parametrisationRange == 0) {
		moleFractionSulfuricAcidCluster = 0.7409967177282139 -0.002663785665140117*temp-0.003499978417957668*lnNumberConcSulfuricAcid+0.0000504021689382576*temp*lnNumberConcSulfuricAcid+0.002010478847383187*lnRh-0.0001832894131464668*temp*lnRh+0.001574072538464286*lnRh2-0.00001790589121766952*temp*lnRh2+0.0001844027436573778*lnRh3-1.503452308794887e-6*temp*lnRh3;

		if (moleFractionSulfuricAcidCluster < 0.17) {
			return;
		}
		if (moleFractionSulfuricAcidCluster > 0.62) {
			return;
		}
		
	} else {
		moleFractionSulfuricAcidCluster = 0.847011529956499-0.002965596681859677*temp-0.006622663242677592*lnNumberConcSulfuricAcid+0.00005878351930618956*temp*lnNumberConcSulfuricAcid+0.05926531222362889*lnRh-0.0003631920982114802*temp*lnRh+0.02300735131008404*lnRh2-0.0000851373757047879*temp*lnRh2+0.002174171592972587*lnRh3-7.923002731097429e-6*temp*lnRh3;
	
		if (moleFractionSulfuricAcidCluster < 0.15) {
			return;
		}
		if (moleFractionSulfuricAcidCluster > 0.54) {
			return;
		}
	}
	
	moleFractionSulfuricAcidCluster2 = SQR(moleFractionSulfuricAcidCluster);
	
	if (parametrisationRange == 0) {
		 parameters1[0] = 0.1430901615568665 +2.219563673425199*temp-0.02739106114964264*temp2+0.00007228107239317088*temp3+5.91822263375044/moleFractionSulfuricAcidCluster;
		 parameters1[1] = 0.1174886643003278+0.4625315047693772*temp-0.01180591129059253*temp2+0.0000404196487152575*temp3+15.79628615047088/moleFractionSulfuricAcidCluster;
		 parameters1[2] = -0.215553951893509-0.0810269192332194*temp+0.001435808434184642*temp2-4.775796947178588e-6*temp3-2.912974063702185/moleFractionSulfuricAcidCluster;
		 parameters1[3] = -3.588557942822751+0.04950795302831703*temp-0.0002138195118737068*temp2+3.108005107949533e-7*temp3-0.02933332747098296/moleFractionSulfuricAcidCluster;
		 parameters1[4] = 1.145983818561277-0.6007956227856778*temp+0.00864244733283759*temp2-0.00002289467254710888*temp3-8.44984513869014/moleFractionSulfuricAcidCluster;
		 parameters1[5] = 2.158548369286559+0.0808121412840917*temp-0.0004073815255395214*temp2-4.019572560156515e-7*temp3+0.7213255852557236/moleFractionSulfuricAcidCluster;
		 parameters1[6] = 1.62409850488771-0.01601062035325362*temp+0.00003771238979714162*temp2+3.217942606371182e-8*temp3-0.01132550810022116/moleFractionSulfuricAcidCluster;
		 parameters1[7] = 9.71681713056504-0.1150478558347306*temp+0.0001570982486038294*temp2+4.009144680125015e-7*temp3+0.7118597859976135/moleFractionSulfuricAcidCluster;
		 parameters1[8] = -1.056105824379897+0.00903377584628419*temp-0.00001984167387090606*temp2+2.460478196482179e-8*temp3-0.05790872906645181/moleFractionSulfuricAcidCluster;
		 parameters1[9] = -0.1487119673397459+0.002835082097822667*temp-9.24618825471694e-6*temp2+5.004267665960894e-9*temp3-0.01270805101481648/moleFractionSulfuricAcidCluster;
		 
	} else {
		 parameters1[0] = -0.001569745387175332-0.1342450856797245*temp+0.1005067182690175*temp2-0.0004601032625944248*temp3+0.1874155345616552/moleFractionSulfuricAcidCluster2+0.0104122021949141/moleFractionSulfuricAcidCluster;
		 parameters1[1] = 0.001950770887486148+0.1680380731937475*temp-0.02257545214945428*temp2+0.0000827148525354003*temp3+0.002502896602782628/moleFractionSulfuricAcidCluster2+0.01552154487373872/moleFractionSulfuricAcidCluster;
		 parameters1[2] = 0.000154083599132041-0.02803005713673646*temp+0.001545869791064423*temp2-4.527012995558728e-6*temp3+0.0915323112970841/moleFractionSulfuricAcidCluster2+0.07116518392455503/moleFractionSulfuricAcidCluster;
		 parameters1[3] = -0.005092666464919691-0.0079684608373626*temp+0.00004468280779999966*temp2-8.79425326170835e-8*temp3+0.1339909859611695/moleFractionSulfuricAcidCluster2+0.831112029665779/moleFractionSulfuricAcidCluster;
		 parameters1[4] = -0.02272226993056786-1.565123155966037*temp+0.003807170798570766*temp2+0.00001641085384332953*temp3+1.294987606109814/moleFractionSulfuricAcidCluster2+0.04748213867867968/moleFractionSulfuricAcidCluster;
		 parameters1[5] = 0.003106456643285619+0.3045177689813154*temp-0.0005640120058411354*temp2-2.032666211033406e-6*temp3-0.3515835365109146/moleFractionSulfuricAcidCluster2+0.1037494194279776/moleFractionSulfuricAcidCluster;
		 parameters1[6] = 0.07754300076253439-0.001963146205019902*temp-0.00001304119073568482*temp2+6.623693300681445e-8*temp3+0.01134695792517722/moleFractionSulfuricAcidCluster2+0.0972804474950142/moleFractionSulfuricAcidCluster;
		 parameters1[7] = -0.1531427628782814+0.05753924492517699*temp-0.0003065107065724253*temp2-2.960965968850208e-8*temp3-0.0982514182025573/moleFractionSulfuricAcidCluster2+0.3362860198323058/moleFractionSulfuricAcidCluster;
		 parameters1[8] = -0.5521732040508044-0.002070431390446516*temp+0.00001440324612102496*temp2+8.83000427843853e-9*temp3+0.01198331714250616/moleFractionSulfuricAcidCluster2-0.07000246307134552/moleFractionSulfuricAcidCluster;
		 parameters1[9] = 0.1265442035303987-0.001360286418822358*temp+5.905978942513288e-6*temp2-4.171501234425879e-9*temp3+0.001708065264911778/moleFractionSulfuricAcidCluster2-0.006432296088358141/moleFractionSulfuricAcidCluster;
	}
	
	J = exp(parameters1[0]+parameters1[1]*lnRh+parameters1[2]*lnRh2+parameters1[3]*lnRh3+parameters1[4]*lnNumberConcSulfuricAcid+parameters1[5]*lnRh*lnNumberConcSulfuricAcid+parameters1[6]*lnRh2*lnNumberConcSulfuricAcid+parameters1[7]*SQR(lnNumberConcSulfuricAcid)+parameters1[8]*lnRh*SQR(lnNumberConcSulfuricAcid)+parameters1[9]*pow(lnNumberConcSulfuricAcid,3.0));
	
	if (parametrisationRange == 0) {
		if (J < 1.0e-7) {
			return;
			
		} else if (J > 1.0e10) {
			J = 1.0e10;
		}
	}
	
    if (parametrisationRange == 1) {
		if (J < 0.1) {
			return;
			
		} else if (J > 1.0e14) {
			J = 1.0e14;
		}
	}
	
	/* from (1/(cm^3 s)) to (1/(m^3 s)) */
	J = J*nucleationCorrectionFactor*1.0e6;
	
	if (parametrisationRange == 0) {
		 parameters2[0] = -0.002954125078716302 -0.0976834264241286*temp+0.001024847927067835*temp2-2.186459697726116e-6*temp3-0.1017165718716887/moleFractionSulfuricAcidCluster;
		 parameters2[1] = -0.002050640345231486-0.007585041382707174*temp+0.0001926539658089536*temp2-6.70429719683894e-7*temp3-0.2557744774673163/moleFractionSulfuricAcidCluster;
		 parameters2[2] = 0.003223076552477191+0.000852636632240633*temp-0.00001547571354871789*temp2+5.666608424980593e-8*temp3+0.03384437400744206/moleFractionSulfuricAcidCluster;
		 parameters2[3] = 0.04743226764572505-0.0006251042204583412*temp+2.650663328519478e-6*temp2-3.674710848763778e-9*temp3-0.0002672510825259393/moleFractionSulfuricAcidCluster;
		 parameters2[4] = -0.01252108546759328+0.005806550506277202*temp-0.0001016735312443444*temp2+2.881946187214505e-7*temp3+0.0942243379396279/moleFractionSulfuricAcidCluster;
		 parameters2[5] = -0.0385459592773097-0.0006723156277391984*temp+2.602884877659698e-6*temp2+1.194163699688297e-8*temp3-0.00851515345806281/moleFractionSulfuricAcidCluster;
		 parameters2[6] = -0.01837488495738111+0.0001720723574407498*temp-3.717657974086814e-7*temp2-5.148746022615196e-10*temp3+0.0002686602132926594/moleFractionSulfuricAcidCluster;
		 parameters2[7] = -0.06199739728812199+0.000906958053583576*temp-9.11727926129757e-7*temp2-5.367963396508457e-9*temp3-0.007742343393937707/moleFractionSulfuricAcidCluster;
		 parameters2[8] = 0.0121827103101659-0.0001066499571188091*temp+2.534598655067518e-7*temp2-3.635186504599571e-10*temp3+0.0006100650851863252/moleFractionSulfuricAcidCluster;
		 parameters2[9] = 0.0003201836700403512-0.0000174761713262546*temp+6.065037668052182e-8*temp2-1.421771723004557e-11*temp3+0.0001357509859501723/moleFractionSulfuricAcidCluster;
		 
    } else {
		 parameters2[0] = 7.510239772212755e-6+0.0005020543717858126*temp-0.00003686018555428167*temp2+1.082558571899109e-6*temp3-0.0002702823459545526/moleFractionSulfuricAcidCluster;
		 parameters2[1] = -4.300475182632171e-6-0.000730132746532732*temp+0.0002520622067147548*temp2-1.016483987242821e-6*temp3-0.001142825737843263/moleFractionSulfuricAcidCluster;
		 parameters2[2] = -4.421559550727262e-6-0.002348602119724249*temp+3.00649923314284e-7*temp2+2.447973651342882e-8*temp3-0.002502258291246978/moleFractionSulfuricAcidCluster;
		 parameters2[3] = -0.0001670569160183319+0.0002075039029766036*temp-1.130128595080389e-6*temp2+1.802682362218632e-9*temp3-0.0168244622970241/moleFractionSulfuricAcidCluster;
		 parameters2[4] = 0.0000985954120966784+0.004512847128248536*temp-0.000051255749468906*temp2+4.607487783359696e-8*temp3-0.002143177270467462/moleFractionSulfuricAcidCluster;
		 parameters2[5] = 0.00006365278270556093-0.002885292129285572*temp+6.51706496257436e-6*temp2+2.326013038006477e-8*temp3-0.01103185560619721/moleFractionSulfuricAcidCluster;
		 parameters2[6] = 0.0004492390899906145+0.00006894159835203836*temp-3.503018074400659e-7*temp2+1.07451436139538e-10*temp3+0.00169646418409145/moleFractionSulfuricAcidCluster;
		 parameters2[7] = 0.000831843754899915-5.351078116383922e-6*temp+1.664320030743163e-6*temp2-3.051079180252276e-9*temp3-0.0003062512887774901/moleFractionSulfuricAcidCluster;
		 parameters2[8] = 0.003553741933599263+0.0000306009050550714*temp-2.110041590842608e-7*temp2-2.114356501960609e-11*temp3+0.0007498898835321218/moleFractionSulfuricAcidCluster;
		 parameters2[9] = -0.00143534432169805+7.855995385362244e-6*temp-3.45127705481677e-8*temp2+5.215472343422993e-11*temp3-0.00002142302928436071/moleFractionSulfuricAcidCluster;
    }
	
	nMolecCluster = exp(parameters2[0]+parameters2[1]*lnRh+parameters2[2]*lnRh2+parameters2[3]*lnRh3+parameters2[4]*lnNumberConcSulfuricAcid+parameters2[5]*lnRh*lnNumberConcSulfuricAcid+parameters2[6]*lnRh2*lnNumberConcSulfuricAcid+parameters2[7]*SQR(lnNumberConcSulfuricAcid)+parameters2[8]*lnRh*SQR(lnNumberConcSulfuricAcid)+parameters2[9]*pow(lnNumberConcSulfuricAcid,3.0));
	
	if (parametrisationRange == 0) {
		nMolecCluster = tutmam_limits(4.0,nMolecCluster,70.0);

	} else {
		nMolecCluster = tutmam_limits(6.0,nMolecCluster,133.0);
	}
	
	mWaterCluster = nMolecCluster*(1.0-moleFractionSulfuricAcidCluster)*2.9915e-26;
	mSulfuricAcidCluster = nMolecCluster*moleFractionSulfuricAcidCluster*1.629e-25;
	sCluster = pow(mWaterCluster+mSulfuricAcidCluster,TUTMAM_23);
	
	sourceTerm[0] = J;
	sourceTerm[1] = J*sCluster;
	sourceTerm[iSpeciesWater+2] = J*mWaterCluster;
	sourceTerm[iSpeciesSulfuricAcid+2] = J*mSulfuricAcidCluster;

	return;
}

/* Olin nucleation rate stored to sourceTerm */
void nucleation_rate_olin(real *sourceTerm, cell_t c, Thread *t) {
	real J;					/* nucleation rate (1/(m^3 s)) */
	real pressure;			/* total pressure (Pa) */
	real temp;				/* temperature (K) */
	int iSpecies;			/* species ID */
	real moleFractionI;		/* mole fraction of vapor i */
	real numberConcI;		/* molecular number concentration of vapor i (1/cm^3) */
	real satVapPresI;		/* saturation vapor pressure of vapor i (Pa) */
	real sCluster;			/* 0.667th moment of the cluster (kg^(2/3)) */
	real mCluster = 0.0;	/* mass of the cluster (kg) */
	real mClusterI;			/* partial mass of the cluster for species i (kg) */
	
	/* nucleation rate converted from (1/(cm^3 s)) to (1/(m^3 s)) */
	J = nucleationCorrectionFactor*1.0e6;
	
	pressure = C_P_TOT(c,t);
	temp = C_T(c,t);
	
	for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) { /* looping over all particle species */
		if (iNucleatingSpecies[iSpecies] == 1) { /* if the species is nucleating */
			
			if (fabs(nucleationExponents[iSpecies]) > 1.0e-3) {	/* if exponent is not 0 */
				moleFractionI = C_XI(c,t,iSpecies);
				numberConcI = 1.0e-6*pressure*moleFractionI/(TUTMAM_BOLTZ*temp);
				J = J*pow(numberConcI,nucleationExponents[iSpecies]); /* multiplying J */
			}
			
			if (fabs(nucleationSatVapPresExponents[iSpecies]) > 1.0e-3) { /* if exponent is not 0 */ 
				satVapPresI = saturation_vapor_pressure(temp,iSpecies);
				J = J/pow(satVapPresI,nucleationSatVapPresExponents[iSpecies]); /* dividing J */
			}
		}
	}
		
	/* source term for 0th moment */
	sourceTerm[0]=J;
	
	/* source terms for 1st moments */
	for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) { /* looping over all particle species */
		if (iNucleatingSpecies[iSpecies] == 1) { /* if the species is nucleating */
			mClusterI = nMolecClusterVector[iSpecies]*molarMassVector[iSpecies]/TUTMAM_AVOGADRO;
			sourceTerm[iSpecies+2] = J*mClusterI;
			mCluster = mCluster + mClusterI; /* total mass of cluster */
		}
	}
	
	/* 0.667th moment calculated from the total mass */
	sCluster = pow(mCluster,TUTMAM_23)/exp(SQR(log(clusterGSD)));
	sourceTerm[1]=J*sCluster;
	
	return;
}

/* Calculating nucleation rate source terms and storing them to sourceTerm vector */
void calculate_nucleation_rate(real *sourceTerm, cell_t c, Thread *t) {
	
	if (nucleationLaw == 1) {
		nucleation_rate_vehkamaki(sourceTerm,c,t);
		
	} else if (nucleationLaw == 2) {
		nucleation_rate_olin(sourceTerm,c,t);
		
	} else {
		Error("Illegal nucleationLaw\n");
	}
	
	return;
}

/* Making under-relaxation to sourceTerm vector */
void under_relax_nucleation_rate(real *sourceTerm, real uRFNucleation, cell_t c, Thread *t) {
	int iSpecies;	/* particle species ID */
	
	sourceTerm[0] = sourceTerm[0]*uRFNucleation + C_UDMI(c,t,iUdmNucleationM0)*(1.0-uRFNucleation);
	sourceTerm[1] = sourceTerm[1]*uRFNucleation + C_UDMI(c,t,iUdmNucleationM23)*(1.0-uRFNucleation);
	
	for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) {
		sourceTerm[iSpecies+2] = sourceTerm[iSpecies+2]*uRFNucleation + C_UDMI(c,t,iUdmNucleationM1First+iSpecies)*(1.0-uRFNucleation);
	}

	return;
}

/* Storing nucleation rate source terms to UDMs */
void store_nucleation_rate(real *sourceTerm, cell_t c, Thread *t) {
	int iSpecies;	/* particle species ID */
	
	C_UDMI(c,t,iUdmNucleationM0) = sourceTerm[0];
	C_UDMI(c,t,iUdmNucleationM23) = sourceTerm[1];
	
	for (iSpecies = 0; iSpecies < nTutmamSpecies; ++iSpecies) {
		C_UDMI(c,t,iUdmNucleationM1First+iSpecies) = sourceTerm[iSpecies+2];
	}

	return;
}

/* Source term macro for nucleation */
DEFINE_SOURCE(nucleation,c,t,dS,eqn)
{	
	real sourceTerm[2+nTutmamSpecies] = { 0.0 };	/* source term (1/(s m^3)), (kg^(2/3)/(s m^3)), or (kg/(s m^3)) */
	int iSpecies;			/* species ID number */
	int iFluentSpecies;		/* Fluent species ID number */
	int iUds;				/* UDS ID number */
	dS[eqn] = 0.0;			/* differential of the source term */
		
	iUds = eqn-EQ_UDS;

	if (iUds == iUdsM0) { /* UDS for 0th moment */
		if (nucleationProcess == 1) {
			if (nucleationComputing == 1) {
				/* calculate nucleation rate and store it to sourceTerm */
				calculate_nucleation_rate(sourceTerm,c,t);
				
				/* make under-relaxation to sourceTerm */
				under_relax_nucleation_rate(sourceTerm,uRFNucleation,c,t);
							
			} else { /* computing off */
				return C_UDMI(c,t,iUdmNucleationM0);
			}
		}
		
		/* store nucleation rate to UDMs */
		store_nucleation_rate(sourceTerm,c,t);
		
		/* return J */
		return sourceTerm[0];
	
	} else if (iUds == iUdsM23) { /* UDS for 0.667th moment */
		/* obtain nucleation rate source term from UDM */
		return C_UDMI(c,t,iUdmNucleationM23);
		
	} else if (iUds >= iUdsM1First && iUds < nUdsPerMode) { /* UDS for 1st moment */
		iSpecies = iUds-iUdsM1First;
		/* obtain nucleation rate source term from UDM */
		return C_UDMI(c,t,iUdmNucleationM1First+iSpecies);
		
	} else if (iUds < 0) { /* species */
		iFluentSpecies = eqn-EQ_SPECIES;
		iSpecies = tutmamSpeciesIdVector[iFluentSpecies];
		
		if (iSpecies > -1) {
			/* returning opposite number of the source term calculated for particle side */
			return 0.0 - C_UDMI(c,t,iUdmNucleationM1First+iSpecies);
			
		} else {
			return 0.0;
		}
		
	} else { /* all other equations */
		return 0.0;
	}
}

