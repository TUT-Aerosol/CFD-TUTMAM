;;CFD-TUTMAM menu generation file
;;Miska Olin (miska.olin@tut.fi), Tampere University of Technology
;;
;;Do not make changes in this file
;;
;;
;;Compilation panel
;;
(define  gui-compilation
(let ((panel #f))
(define (update-cb . args) ())
(define (apply-cb . args) ())
(lambda args  
(if (not panel) 
(let ((table) (form)) ;;Creating the GUI panel.
(set! panel (cx-create-panel "Compilation" apply-cb update-cb))    
(set! table (cx-create-table panel "" 'border #f 'below 0 'right-of 0))
(set! form (cx-create-frame table "" 'border #f))
(cx-create-button table "1. Unload" 'row 1 'col 0 'activate-callback
	(lambda ()
		(ti-menu-load-string "f rj tutmam.unload")
	))
(cx-create-button table "2. Compile" 'row 1 'col 1 'activate-callback
	(lambda ()
		(ti-menu-load-string "f rj tutmam.compile")
	))
(cx-create-button table "3. Load" 'row 1 'col 2 'activate-callback
	(lambda ()
		(ti-menu-load-string "f rj tutmam.load")
	))
)  
) (cx-show-panel panel)    )))
;;
;;
;;General settings panel
;;
;;Setting default values for the RP-variables, if they are not defined.
(if (not (rp-var-object 'fluid_zone_number_rp))(rp-var-define 'fluid_zone_number_rp 1 'int #f))
(if (not (rp-var-object 'gauss_hermite_level_rp))(rp-var-define 'gauss_hermite_level_rp 7 'int #f))
;;
(define  gui-gen-sett
(let ((panel #f)(fluidZoneNumber)(gaussHermiteLevel)) ;;Declaring the scheme-variables
(define (update-cb . args) ;;Update: setting the values of scheme-variables to the values of RP-variables
(cx-set-integer-entry fluidZoneNumber (%rpgetvar 'fluid_zone_number_rp))      
(cx-set-integer-entry gaussHermiteLevel (%rpgetvar 'gauss_hermite_level_rp))
  )    
(define (apply-cb . args) ;;Apply: setting the values of RP-variables to the values of scheme-variables 
(rpsetvar 'fluid_zone_number_rp (cx-show-integer-entry fluidZoneNumber)) 
(rpsetvar 'gauss_hermite_level_rp (cx-show-integer-entry gaussHermiteLevel)) 
(%udf-on-demand "tutmam_transfer_settings_manually_general_settings::tutmamudf") ;;Settings are transferred to C-code
   )    
(lambda args  
(if (not panel) 
(let ((table) (form)) ;;Creating the GUI panel. Scheme variables are connected to panel entries.
(set! panel (cx-create-panel "General settings" apply-cb update-cb))    
(set! table (cx-create-table panel "" 'border #f 'below 0 'right-of 0))
(set! form (cx-create-frame table "" 'border #f))
(cx-create-button table "Show UDS and UDM mapping" 'row 1 'col 0 'activate-callback
	(lambda ()
		(%udf-on-demand "show_uds_udm_mapping_manually::tutmamudf")
	))
(set! fluidZoneNumber (cx-create-integer-entry table "Fluid zone ID" 'width 14 'row 2 'col 0 ))
(set! gaussHermiteLevel (cx-create-integer-entry table "Gauss-Hermite quad. level" 'width 14 'row 3 'col 0 ))
(cx-create-button panel "Apply" 'panel-button #t 'activate-callback apply-cb)
)  
) (cx-show-panel panel)    )))
;;
;;
;;Aerosol distribution panel
;;
(if (not (rp-var-object 'min_gsd_rp))(rp-var-define 'min_gsd_rp 1.01 'real #f))
(if (not (rp-var-object 'max_gsd_rp))(rp-var-define 'max_gsd_rp 3.0 'real #f))
;;
(define  gui-distribution
(let ((panel #f)(minGSD)(maxGSD))
(define (update-cb . args)    
(cx-set-real-entry minGSD (%rpgetvar 'min_gsd_rp))      
(cx-set-real-entry maxGSD (%rpgetvar 'max_gsd_rp))
  )    
(define (apply-cb . args)      
(rpsetvar 'min_gsd_rp (cx-show-real-entry minGSD)) 
(rpsetvar 'max_gsd_rp (cx-show-real-entry maxGSD)) 
(%udf-on-demand "tutmam_transfer_settings_manually_aerosol_distribution::tutmamudf")
   )    
(lambda args  
(if (not panel) 
(let ((table) (form))
(set! panel (cx-create-panel "Log-normal distribution" apply-cb update-cb))    
(set! table (cx-create-table panel "" 'border #f 'below 0 'right-of 0))
(set! form (cx-create-frame table "" 'border #f))
(set! minGSD (cx-create-real-entry table "Minimum GSD" 'width 14 'row 1 'col 0 ))
(set! maxGSD (cx-create-real-entry table "Maximum GSD" 'width 14 'row 2 'col 0 ))
(cx-create-button panel "Apply" 'panel-button #t 'activate-callback apply-cb)
)  
) (cx-show-panel panel)    )))
;;
;;
;;Aerosol process control panel
;;
(if (not (rp-var-object 'diffusion_process_rp))(rp-var-define 'diffusion_process_rp #t 'boolean #f))
(if (not (rp-var-object 'nucleation_process_rp))(rp-var-define 'nucleation_process_rp #f 'boolean #f))
(if (not (rp-var-object 'nucleation_computing_rp))(rp-var-define 'nucleation_computing_rp #t 'boolean #f))
(if (not (rp-var-object 'condensation_process_rp))(rp-var-define 'condensation_process_rp #f 'boolean #f))
(if (not (rp-var-object 'condensation_computing_rp))(rp-var-define 'condensation_computing_rp #t 'boolean #f))
(if (not (rp-var-object 'coagulation_process_rp))(rp-var-define 'coagulation_process_rp #f 'boolean #f))
(if (not (rp-var-object 'coagulation_computing_rp))(rp-var-define 'coagulation_computing_rp #t 'boolean #f))
(if (not (rp-var-object 'self_coagulational_transfer_process_rp))(rp-var-define 'self_coagulational_transfer_process_rp #f 'boolean #f))
(if (not (rp-var-object 'inter_modal_condensation_process_rp))(rp-var-define 'inter_modal_condensation_process_rp #f 'boolean #f))
;;
(define  gui-processcontrol
(let ((panel #f)(diffusionProcess)(nucleationProcess)(nucleationComputing)(condensationProcess)(condensationComputing)(coagulationProcess)(coagulationComputing)(selfCoagulationalTransferProcess)(interModalCondensationProcess))
(define (update-cb . args)    
(cx-set-toggle-button diffusionProcess (%rpgetvar 'diffusion_process_rp))     
(cx-set-toggle-button nucleationProcess (%rpgetvar 'nucleation_process_rp))   
(cx-set-toggle-button nucleationComputing (%rpgetvar 'nucleation_computing_rp))   
(cx-set-toggle-button condensationProcess (%rpgetvar 'condensation_process_rp))   
(cx-set-toggle-button condensationComputing (%rpgetvar 'condensation_computing_rp))   
(cx-set-toggle-button coagulationProcess (%rpgetvar 'coagulation_process_rp))   
(cx-set-toggle-button coagulationComputing (%rpgetvar 'coagulation_computing_rp))   
(cx-set-toggle-button selfCoagulationalTransferProcess (%rpgetvar 'self_coagulational_transfer_process_rp))   
(cx-set-toggle-button interModalCondensationProcess (%rpgetvar 'inter_modal_condensation_process_rp))   
  )    
(define (apply-cb . args)      
(rpsetvar 'diffusion_process_rp (cx-show-toggle-button diffusionProcess)) 
(rpsetvar 'nucleation_process_rp (cx-show-toggle-button nucleationProcess)) 
(rpsetvar 'nucleation_computing_rp (cx-show-toggle-button nucleationComputing)) 
(rpsetvar 'condensation_process_rp (cx-show-toggle-button condensationProcess)) 
(rpsetvar 'condensation_computing_rp (cx-show-toggle-button condensationComputing)) 
(rpsetvar 'coagulation_process_rp (cx-show-toggle-button coagulationProcess)) 
(rpsetvar 'coagulation_computing_rp (cx-show-toggle-button coagulationComputing)) 
(rpsetvar 'self_coagulational_transfer_process_rp (cx-show-toggle-button selfCoagulationalTransferProcess)) 
(rpsetvar 'inter_modal_condensation_process_rp (cx-show-toggle-button interModalCondensationProcess)) 
(%udf-on-demand "tutmam_transfer_settings_manually_aerosol_process_control::tutmamudf")
   )    
(lambda args  
(if (not panel) 
(let ((table) (table2)(form))
(set! panel (cx-create-panel "Aerosol process control" apply-cb update-cb))    
(set! table (cx-create-table panel "" 'border #t 'below 0 'right-of 0))
(set! form (cx-create-frame table "" 'border #f))
(set! diffusionProcess (cx-create-toggle-button table "Diffusion" 'state #t 'row 1 'col 0 ))
(set! nucleationProcess (cx-create-toggle-button table "Nucleation" 'state #f 'row 2 'col 0 'activate-callback
	(lambda ()
		(if (cx-show-toggle-button nucleationProcess)
		(cx-show-item nucleationComputing)
		(cx-hide-item nucleationComputing)
	))
))
(set! nucleationComputing (cx-create-toggle-button table "Compute" 'state #t 'row 2 'col 1 ))
(set! condensationProcess (cx-create-toggle-button table "Condensation" 'state #f 'row 3 'col 0 'activate-callback
	(lambda ()
		(if (cx-show-toggle-button condensationProcess)
		(cx-show-item condensationComputing)
		(cx-hide-item condensationComputing)
	))
))
(set! condensationComputing (cx-create-toggle-button table "Compute" 'state #t 'row 3 'col 1 ))
(set! coagulationProcess (cx-create-toggle-button table "Coagulation" 'state #f 'row 4 'col 0 'activate-callback
	(lambda ()
		(if (cx-show-toggle-button coagulationProcess)
		(cx-show-item coagulationComputing)
		(cx-hide-item coagulationComputing)
	))
))
(set! coagulationComputing (cx-create-toggle-button table "Compute" 'state #t 'row 4 'col 1 ))
(set! table2 (cx-create-table panel "Only if using power-law distribution" 'border #t 'below table 'right-of 0))
(set! selfCoagulationalTransferProcess (cx-create-toggle-button table2 "Self-coagulational transfer" 'state #f 'row 1 'col 0 ))
(set! interModalCondensationProcess (cx-create-toggle-button table2 "Intermodal condensation" 'state #f 'row 2 'col 0 ))
(if (%rpgetvar 'nucleation_process_rp)
	(cx-show-item nucleationComputing)
	(cx-hide-item nucleationComputing)
)
(if (%rpgetvar 'condensation_process_rp)
	(cx-show-item condensationComputing)
	(cx-hide-item condensationComputing)
)
(if (%rpgetvar 'coagulation_process_rp)
	(cx-show-item coagulationComputing)
	(cx-hide-item coagulationComputing)
)
(cx-create-button panel "Apply" 'panel-button #t 'activate-callback apply-cb)
)  
) (cx-show-panel panel)    )))
;;
;;
;;Diffusion panel
;;
(if (not (rp-var-object 'diffusion_law_rp))(rp-var-define 'diffusion_law_rp #t 'boolean #f))
(if (not (rp-var-object 'slip_correction_law_rp))(rp-var-define 'slip_correction_law_rp #t 'boolean #f))
(if (not (rp-var-object 'constant_slip_correction_factor_rp))(rp-var-define 'constant_slip_correction_factor_rp 1.0 'real #f))
(if (not (rp-var-object 'particle_turb_diff_rp))(rp-var-define 'particle_turb_diff_rp #t 'boolean #f))
(if (not (rp-var-object 'particle_turb_schmidt_number_rp))(rp-var-define 'particle_turb_schmidt_number_rp 0.7 'real #f))
;;
(define  gui-diff
(let ((panel #f)(diffusionLaw)(slipCorrectionLaw)(constantSlipCorrectionFactor)(particleTurbDiff)(particleTurbSchmidtNumber))
(define (update-cb . args)    
(cx-set-toggle-button diffusionLaw (%rpgetvar 'diffusion_law_rp))      
(cx-set-toggle-button slipCorrectionLaw (%rpgetvar 'slip_correction_law_rp))   
(cx-set-real-entry constantSlipCorrectionFactor (%rpgetvar 'constant_slip_correction_factor_rp))  
(cx-set-toggle-button particleTurbDiff (%rpgetvar 'particle_turb_diff_rp))   
(cx-set-real-entry particleTurbSchmidtNumber (%rpgetvar 'particle_turb_schmidt_number_rp))  
  )    
(define (apply-cb . args)      
(rpsetvar 'diffusion_law_rp (cx-show-toggle-button diffusionLaw)) 
(rpsetvar 'slip_correction_law_rp (cx-show-toggle-button slipCorrectionLaw))      
(rpsetvar 'constant_slip_correction_factor_rp (cx-show-real-entry constantSlipCorrectionFactor))
(rpsetvar 'particle_turb_diff_rp (cx-show-toggle-button particleTurbDiff))     
(rpsetvar 'particle_turb_schmidt_number_rp (cx-show-real-entry particleTurbSchmidtNumber))
(%udf-on-demand "tutmam_transfer_settings_manually_diffusion::tutmamudf")
   )    
(lambda args  
(if (not panel) 
(let ((table) (form))
(set! panel (cx-create-panel "Diffusion" apply-cb update-cb))    
(set! table (cx-create-table panel "" 'border #f 'below 0 'right-of 0))
(set! form (cx-create-frame table "" 'border #f))
(set! diffusionLaw (cx-create-toggle-button table "Diffusion integral parametrization" 'state #t 'row 3 'col 0 ))
(set! slipCorrectionLaw (cx-create-toggle-button table "Calculate slip corr. coeff." 'state #t 'row 1 'col 0 'activate-callback
	(lambda ()
		(if (cx-show-toggle-button slipCorrectionLaw)
		(cx-hide-item constantSlipCorrectionFactor)
		(cx-show-item constantSlipCorrectionFactor)
	)
		(if (cx-show-toggle-button slipCorrectionLaw)
		(cx-show-item diffusionLaw)
		(cx-hide-item diffusionLaw)
	))
))
(set! constantSlipCorrectionFactor (cx-create-real-entry table "Constant slip corr. coeff." 'width 14 'row 2 'col 0 ))
(set! particleTurbDiff (cx-create-toggle-button table "Particle turbulent diffusion" 'state #t 'row 4 'col 0 'activate-callback
	(lambda ()
		(if (cx-show-toggle-button particleTurbDiff)
		(cx-show-item particleTurbSchmidtNumber)
		(cx-hide-item particleTurbSchmidtNumber)
	))
))
(set! particleTurbSchmidtNumber (cx-create-real-entry table "Particle turbulent Schmidt number" 'width 14 'row 5 'col 0 ))
(if (%rpgetvar 'slip_correction_law_rp)
	(cx-hide-item constantSlipCorrectionFactor)
	(cx-show-item constantSlipCorrectionFactor)
)
(if (%rpgetvar 'slip_correction_law_rp)
	(cx-show-item diffusionLaw)
	(cx-hide-item diffusionLaw)
)
(if (%rpgetvar 'particle_turb_diff_rp)
	(cx-show-item particleTurbSchmidtNumber)
	(cx-hide-item particleTurbSchmidtNumber)
)
(cx-create-button panel "Apply" 'panel-button #t 'activate-callback apply-cb)
)  
) (cx-show-panel panel)    )))
;;
;;
;;Nucleation panel
;;
(if (not (rp-var-object 'urf_nucleation_rp))(rp-var-define 'urf_nucleation_rp 0.9 'real #f))
(if (not (rp-var-object 'nucleation_correction_factor_rp))(rp-var-define 'nucleation_correction_factor_rp 1.0 'real #f))
(if (not (rp-var-object 'nucleation_law_rp))(rp-var-define 'nucleation_law_rp 1 'int #f))
(if (not (rp-var-object 'nucleating_species_rp))(rp-var-define 'nucleating_species_rp "1,1" 'string #f))
(if (not (rp-var-object 'nucleation_exponents_rp))(rp-var-define 'nucleation_exponents_rp "1,1" 'string #f))
(if (not (rp-var-object 'nucleation_sat_vap_pres_exponents_rp))(rp-var-define 'nucleation_sat_vap_pres_exponents_rp "1,0" 'string #f))
(if (not (rp-var-object 'n_molec_cluster_vector_rp))(rp-var-define 'n_molec_cluster_vector_rp "15,20" 'string #f))
;;
(define  gui-nucleation
(let ((panel #f)(uRFNucleation)(nucleationCorrectionFactor)(nucleationLaw)(nucleatingSpecies)(nucleationExponents)(nucleationSatVapPresExponents)(nMolecClusterVector))
(define (update-cb . args)    
(cx-set-real-entry uRFNucleation (%rpgetvar 'urf_nucleation_rp))        
(cx-set-real-entry nucleationCorrectionFactor (%rpgetvar 'nucleation_correction_factor_rp))        
(cx-set-integer-entry nucleationLaw (%rpgetvar 'nucleation_law_rp))     
(cx-set-text-entry nucleatingSpecies (%rpgetvar 'nucleating_species_rp))     
(cx-set-text-entry nucleationExponents (%rpgetvar 'nucleation_exponents_rp))     
(cx-set-text-entry nucleationSatVapPresExponents (%rpgetvar 'nucleation_sat_vap_pres_exponents_rp))     
(cx-set-text-entry nMolecClusterVector (%rpgetvar 'n_molec_cluster_vector_rp)) 
  )    
(define (apply-cb . args)      
(rpsetvar 'urf_nucleation_rp (cx-show-real-entry uRFNucleation)) 
(rpsetvar 'nucleation_correction_factor_rp (cx-show-real-entry nucleationCorrectionFactor)) 
(rpsetvar 'nucleation_law_rp (cx-show-integer-entry nucleationLaw)) 
(rpsetvar 'nucleating_species_rp (cx-show-text-entry nucleatingSpecies)) 
(rpsetvar 'nucleation_exponents_rp (cx-show-text-entry nucleationExponents)) 
(rpsetvar 'nucleation_sat_vap_pres_exponents_rp (cx-show-text-entry nucleationSatVapPresExponents)) 
(rpsetvar 'n_molec_cluster_vector_rp (cx-show-text-entry nMolecClusterVector)) 
(%udf-on-demand "tutmam_transfer_settings_manually_nucleation::tutmamudf")
   )    
(lambda args  
(if (not panel) 
(let ((table)(table2) (form))
(set! panel (cx-create-panel "Nucleation" apply-cb update-cb))    
(set! table (cx-create-table panel "" 'border #f 'below 0 'right-of 0))
(set! form (cx-create-frame table "" 'border #f))
(set! uRFNucleation (cx-create-real-entry table "URF for nucleation" 'width 14 'row 1 'col 0 ))
(set! nucleationCorrectionFactor (cx-create-real-entry table "Correction factor" 'width 14 'row 2 'col 0 ))
(set! nucleationLaw (cx-create-integer-entry table "Nucleation law ID number" 'width 14 'row 3 'col 0 'activate-callback
	(lambda ()
		(if (= (cx-show-integer-entry nucleationLaw) 2)
		(cx-show-item table2)
		(cx-hide-item table2)
	))
))
(set! table2 (cx-create-table panel "Olin nucleation law" 'border #t 'below table 'right-of 0))
(set! nucleatingSpecies (cx-create-text-entry table2 "Nucleating species flags" 'width 14 'row 4 'col 0 ))
(set! nucleationExponents (cx-create-text-entry table2 "Nucleation exponents" 'width 14 'row 5 'col 0 ))
(set! nucleationSatVapPresExponents (cx-create-text-entry table2 "Nucleation sat. vap. pres. exponents" 'width 14 'row 6 'col 0 ))
(set! nMolecClusterVector (cx-create-text-entry table2 "Number of molecules in a cluster" 'width 14 'row 7 'col 0 ))
(if (= (%rpgetvar 'nucleation_law_rp) 2)
	(cx-show-item table2)
	(cx-hide-item table2)
)
(cx-create-button panel "Apply" 'panel-button #t 'activate-callback apply-cb)
)  
) (cx-show-panel panel)    )))
;;
;;
;;Condensation panel
;;
(if (not (rp-var-object 'urf_condensation_m23_rp))(rp-var-define 'urf_condensation_m23_rp 1.0 'real #f))
(if (not (rp-var-object 'urf_condensation_m1_rp))(rp-var-define 'urf_condensation_m1_rp "1,1" 'string #f))
(if (not (rp-var-object 'latent_heat_of_condensation_rp))(rp-var-define 'latent_heat_of_condensation_rp #t 'boolean #f))
(if (not (rp-var-object 'condensing_species_rp))(rp-var-define 'condensing_species_rp "1,1" 'string #f))
(if (not (rp-var-object 'kelvin_effect_rp))(rp-var-define 'kelvin_effect_rp #f 'boolean #f))
(if (not (rp-var-object 'phase_activity_model_rp))(rp-var-define 'phase_activity_model_rp #f 'boolean #f))
(if (not (rp-var-object 'condensation_direction_vector_rp))(rp-var-define 'condensation_direction_vector_rp "0,0" 'string #f))
(if (not (rp-var-object 'water_eq_rp))(rp-var-define 'water_eq_rp #f 'boolean #f))
(if (not (rp-var-object 'water_eq_species_rp))(rp-var-define 'water_eq_species_rp 0 'int #f))
(if (not (rp-var-object 'min_kappa_rp))(rp-var-define 'min_kappa_rp 0.01 'real #f))
(if (not (rp-var-object 'max_kappa_rp))(rp-var-define 'max_kappa_rp 5.0 'real #f))
(if (not (rp-var-object 'diesel_exhaust_hc_fraction_model_rp))(rp-var-define 'diesel_exhaust_hc_fraction_model_rp #f 'boolean #f))
(if (not (rp-var-object 'diesel_exhaust_hc_rp))(rp-var-define 'diesel_exhaust_hc_rp 0 'int #f))
(if (not (rp-var-object 'inter_modal_condensation_factor_rp))(rp-var-define 'inter_modal_condensation_factor_rp 0.0 'real #f))
;;
(define  gui-condensation
(let ((panel #f)(uRFCondensationM23)(uRFCondensationM1)(latentHeatOfCondensation)(condensingSpecies)(kelvinEffect)(phaseActivityModel)(condensationDirectionVector)(waterEq)(waterEqSpecies)(minKappa)(maxKappa)(dieselExhaustHCFractionModel)(dieselExhaustHC)(interModalCondensationFactor))
(define (update-cb . args)    
(cx-set-real-entry uRFCondensationM23 (%rpgetvar 'urf_condensation_m23_rp))      
(cx-set-text-entry uRFCondensationM1 (%rpgetvar 'urf_condensation_m1_rp))         
(cx-set-toggle-button latentHeatOfCondensation (%rpgetvar 'latent_heat_of_condensation_rp))  
(cx-set-text-entry condensingSpecies (%rpgetvar 'condensing_species_rp))    
(cx-set-toggle-button kelvinEffect (%rpgetvar 'kelvin_effect_rp))   
(cx-set-toggle-button phaseActivityModel (%rpgetvar 'phase_activity_model_rp))   
(cx-set-text-entry condensationDirectionVector (%rpgetvar 'condensation_direction_vector_rp))    
(cx-set-toggle-button waterEq (%rpgetvar 'water_eq_rp))   
(cx-set-integer-entry waterEqSpecies (%rpgetvar 'water_eq_species_rp))   
(cx-set-real-entry minKappa (%rpgetvar 'min_kappa_rp))   
(cx-set-real-entry maxKappa (%rpgetvar 'max_kappa_rp))   
(cx-set-toggle-button dieselExhaustHCFractionModel (%rpgetvar 'diesel_exhaust_hc_fraction_model_rp))   
(cx-set-integer-entry dieselExhaustHC (%rpgetvar 'diesel_exhaust_hc_rp))   
(cx-set-real-entry interModalCondensationFactor (%rpgetvar 'inter_modal_condensation_factor_rp))   
  )    
(define (apply-cb . args)      
(rpsetvar 'urf_condensation_m23_rp (cx-show-real-entry uRFCondensationM23)) 
(rpsetvar 'urf_condensation_m1_rp (cx-show-text-entry uRFCondensationM1)) 
(rpsetvar 'latent_heat_of_condensation_rp (cx-show-toggle-button latentHeatOfCondensation)) 
(rpsetvar 'condensing_species_rp (cx-show-text-entry condensingSpecies)) 
(rpsetvar 'kelvin_effect_rp (cx-show-toggle-button kelvinEffect)) 
(rpsetvar 'phase_activity_model_rp (cx-show-toggle-button phaseActivityModel)) 
(rpsetvar 'condensation_direction_vector_rp (cx-show-text-entry condensationDirectionVector)) 
(rpsetvar 'water_eq_rp (cx-show-toggle-button waterEq)) 
(rpsetvar 'water_eq_species_rp (cx-show-integer-entry waterEqSpecies)) 
(rpsetvar 'min_kappa_rp (cx-show-real-entry minKappa)) 
(rpsetvar 'max_kappa_rp (cx-show-real-entry maxKappa)) 
(rpsetvar 'diesel_exhaust_hc_fraction_model_rp (cx-show-toggle-button dieselExhaustHCFractionModel)) 
(rpsetvar 'diesel_exhaust_hc_rp (cx-show-integer-entry dieselExhaustHC)) 
(rpsetvar 'inter_modal_condensation_factor_rp (cx-show-real-entry interModalCondensationFactor)) 
(%udf-on-demand "tutmam_transfer_settings_manually_condensation::tutmamudf")
   )    
(lambda args  
(if (not panel) 
(let ((table)(table2)(table2s)(table3)(table3s) (form))
(set! panel (cx-create-panel "Condensation" apply-cb update-cb))    
(set! table (cx-create-table panel "" 'border #f 'below 0 'right-of 0))
(set! form (cx-create-frame table "" 'border #f))
(set! uRFCondensationM23 (cx-create-real-entry table "URF for surface area moment condensation" 'width 14 'row 1 'col 0 ))
(set! uRFCondensationM1 (cx-create-text-entry table "URF vector for mass moment condensation" 'width 14 'row 2 'col 0 ))
(set! latentHeatOfCondensation (cx-create-toggle-button table "Latent heat effect included in condensation" 'state #t 'row 3 'col 0 ))
(set! condensingSpecies (cx-create-text-entry table "Condensing species flags" 'width 14 'row 4 'col 0 ))
(set! kelvinEffect (cx-create-toggle-button table "Kelvin effect" 'state #f 'row 5 'col 0 ))
(set! condensationDirectionVector (cx-create-text-entry table "Condensing directions" 'width 14 'row 6 'col 0 ))
(set! phaseActivityModel (cx-create-toggle-button table "Phase activity model" 'state #f 'row 7 'col 0 ))
(set! interModalCondensationFactor (cx-create-real-entry table "Intermodal condensation factor (PL+LN)" 'width 14 'row 8 'col 0 ))
(set! table2 (cx-create-table panel "" 'border #t 'below 0 'right-of table))
(set! waterEq (cx-create-toggle-button table2 "Water equilibrium model" 'state #f 'row 1 'col 0 'activate-callback
	(lambda ()
		(if (cx-show-toggle-button waterEq)
		(cx-show-item table2s)
		(cx-hide-item table2s)
	))
))
(set! table2s (cx-create-table table2 "" 'border #f 'row 2 'col 0))
(set! waterEqSpecies (cx-create-integer-entry table2s "Connecting species" 'width 14 'row 1 'col 0 ))
(set! minKappa (cx-create-real-entry table2s "Min. Kappa value" 'width 14 'row 2 'col 0 ))
(set! maxKappa (cx-create-real-entry table2s "Max. Kappa value" 'width 14 'row 3 'col 0 ))
(set! table3 (cx-create-table panel "" 'border #t 'below 0 'right-of table2))
(set! dieselExhaustHCFractionModel (cx-create-toggle-button table3 "Diesel exhaust HC fraction model" 'state #f 'row 1 'col 0 'activate-callback
	(lambda ()
		(if (cx-show-toggle-button dieselExhaustHCFractionModel)
		(cx-show-item table3s)
		(cx-hide-item table3s)
	))
))
(set! table3s (cx-create-table table3 "" 'border #f 'row 2 'col 0))
(set! dieselExhaustHC (cx-create-integer-entry table3s "HC species ID" 'width 14 'row 1 'col 0 ))
(if (%rpgetvar 'water_eq_rp)
	(cx-show-item table2s)
	(cx-hide-item table2s)
)
(if (%rpgetvar 'diesel_exhaust_hc_fraction_model_rp)
	(cx-show-item table3s)
	(cx-hide-item table3s)
)
(cx-create-button panel "Apply" 'panel-button #t 'activate-callback apply-cb)
)  
) (cx-show-panel panel)    )))
;;
;;
;;Coagulation panel
;;
(if (not (rp-var-object 'intra_modal_coagulation_process_rp))(rp-var-define 'intra_modal_coagulation_process_rp #t 'boolean #f))
(if (not (rp-var-object 'inter_modal_coagulation_process_rp))(rp-var-define 'inter_modal_coagulation_process_rp #t 'boolean #f))
(if (not (rp-var-object 'urf_coagulation_rp))(rp-var-define 'urf_coagulation_rp 0.5 'real #f))
(if (not (rp-var-object 'coagulation_iteration_skip_rp))(rp-var-define 'coagulation_iteration_skip_rp 10 'int #f))
(if (not (rp-var-object 'transition_regime_correction_factor_for_coagulation_law_rp))(rp-var-define 'transition_regime_correction_factor_for_coagulation_law_rp 2 'int #f))
(if (not (rp-var-object 'coagulation_robustness_model_rp))(rp-var-define 'coagulation_robustness_model_rp #f 'boolean #f))
;;
(define  gui-coagulation
(let ((panel #f)(intraModalCoagulationProcess)(interModalCoagulationProcess)(uRFCoagulation)(coagulationIterationSkip)(transitionRegimeCorrectionFactorForCoagulationLaw)(coagulationRobustnessModel))
(define (update-cb . args)    
(cx-set-toggle-button intraModalCoagulationProcess (%rpgetvar 'intra_modal_coagulation_process_rp))   
(cx-set-toggle-button interModalCoagulationProcess (%rpgetvar 'inter_modal_coagulation_process_rp))   
(cx-set-real-entry uRFCoagulation (%rpgetvar 'urf_coagulation_rp))      
(cx-set-integer-entry coagulationIterationSkip (%rpgetvar 'coagulation_iteration_skip_rp))   
(cx-set-integer-entry transitionRegimeCorrectionFactorForCoagulationLaw (%rpgetvar 'transition_regime_correction_factor_for_coagulation_law_rp)) 
(cx-set-toggle-button coagulationRobustnessModel (%rpgetvar 'coagulation_robustness_model_rp))   
  )    
(define (apply-cb . args)      
(rpsetvar 'intra_modal_coagulation_process_rp (cx-show-toggle-button intraModalCoagulationProcess)) 
(rpsetvar 'inter_modal_coagulation_process_rp (cx-show-toggle-button interModalCoagulationProcess)) 
(rpsetvar 'urf_coagulation_rp (cx-show-real-entry uRFCoagulation)) 
(rpsetvar 'coagulation_iteration_skip_rp (cx-show-integer-entry coagulationIterationSkip)) 
(rpsetvar 'transition_regime_correction_factor_for_coagulation_law_rp (cx-show-integer-entry transitionRegimeCorrectionFactorForCoagulationLaw)) 
(rpsetvar 'coagulation_robustness_model_rp (cx-show-toggle-button coagulationRobustnessModel)) 
(%udf-on-demand "tutmam_transfer_settings_manually_coagulation::tutmamudf")
   )    
(lambda args  
(if (not panel) 
(let ((table)(table2)(table2s)(table3) (form))
(set! panel (cx-create-panel "Coagulation" apply-cb update-cb))    
(set! table (cx-create-table panel "" 'border #f 'below 0 'right-of 0))
(set! form (cx-create-frame table "" 'border #f))
(set! transitionRegimeCorrectionFactorForCoagulationLaw (cx-create-integer-entry table "Transition regime correction law" 'width 14 'row 1 'col 0 ))
(set! table2 (cx-create-table panel "" 'border #t 'below 0 'right-of table))
(set! intraModalCoagulationProcess (cx-create-toggle-button table2 "Intramodal" 'state #t 'row 1 'col 0 'activate-callback
	(lambda ()
		(if (cx-show-toggle-button intraModalCoagulationProcess)
		(cx-show-item table2s)
		(cx-hide-item table2s)
	))
))
(set! table2s (cx-create-table table2 "" 'border #f 'row 2 'col 0))
(set! coagulationIterationSkip (cx-create-integer-entry table2s "Iteration skip" 'width 14 'row 1 'col 0 ))
(set! uRFCoagulation (cx-create-real-entry table2s "URF" 'width 14 'row 2 'col 0 ))
(set! coagulationRobustnessModel (cx-create-toggle-button table2s "Robustness model" 'state #f 'row 3 'col 0 ))
(set! table3 (cx-create-table panel "" 'border #t 'below 0 'right-of table2))
(set! interModalCoagulationProcess (cx-create-toggle-button table3 "Intermodal" 'state #t 'row 1 'col 0 ))

(cx-create-button panel "Apply" 'panel-button #t 'activate-callback apply-cb)
)  
) (cx-show-panel panel)    )))
;;
;;
;;CFD-TUTMAM menu is created to Fluent's GUI
;;
(let ((menu (cx-add-menu "CFD-TUTMAM" #\U )))
(cx-add-item menu "Compilation" #\U #f cx-client? gui-compilation)
(cx-add-item menu "General settings" #\U #f cx-client? gui-gen-sett)
(cx-add-item menu "Log-normal distribution" #\U #f cx-client? gui-distribution)
(cx-add-item menu "Aerosol process control" #\U #f cx-client? gui-processcontrol)
(cx-add-item menu "Diffusion" #\U #f cx-client? gui-diff)
(cx-add-item menu "Nucleation" #\U #f cx-client? gui-nucleation)
(cx-add-item menu "Condensation" #\U #f cx-client? gui-condensation)
(cx-add-item menu "Coagulation" #\U #f cx-client? gui-coagulation)
(set! *cx-exit-on-error* #f)
) 
;;Do not care about the error text "cx-exit-on-error" if you saw it here, it is not an error message.
;;
(display "************************************************************")
(newline)
(display "      TUT Modal Aerosol Model for CFD 1.0 menus loaded      ")
(newline)
(display "************************************************************")
(newline)
