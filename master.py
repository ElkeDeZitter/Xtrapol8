# -*- coding: utf-8 -*-
"""
Use a phil file with all input parameters

-------

authors and contact information
-------
Elke De Zitter - elke.de-zitter@ibs.fr
Nicolas Coquelle - nicolas.coquelle@esrf.fr
Thomas Barends - Thomas.Barends@mpimf-heidelberg.mpg.de
Jacques Philippe Colletier - jacques-Philippe.colletier@ibs.fr

-------

license information
-------
Copyright (c) 2021 Elke De Zitter, Nicolas Coquelle, Thomas Barends and Jacques-Philippe Colletier
see https://github.com/ElkeDeZitter/Xtrapol8/blob/main/LICENSE

-------
"""
import iotbx.phil

master_phil = iotbx.phil.parse("""
input{
    reference_mtz = None
        .type = path
        .help = Reference data in mtz or mmcif format (merged).
        .expert_level = 0
    triggered_mtz = None
        .type = path
        .help = Triggered data in mtz or mmcif format (merged).
        .expert_level = 0
    reference_pdb = None
        .type = path
        .help = Reference coordinates in pdb or mmcif format.
        .expert_level = 0
    additional_files = None
        .type = path
        .multiple = True 
        .help = Additional files required for refinement, e.g ligand cif file, restraints file.
        .expert_level = 0
    high_resolution = None
        .type = float
        .help = High resolution cutoff (Angstrom). Will only be used if high resolution of the input data files extends to this value.
        .expert_level = 0
    low_resolution = None
        .type = float
        .help = Low resolution cutoff (Angstrom).
        .expert_level = 0
    }
occupancies{
    low_occ = 0.1
        .type = float(value_min=0, value_max=1)
        .help = Lowest occupancy to test (fractional).
        .expert_level = 0
    high_occ = 0.5
        .type = float(value_min=0, value_max=1)
        .help = Highest occupancy to test (fractional).
        .expert_level = 0
    steps = 5
        .type = int
        .help = Amount of equaly spaced occupancies to be tested.
        .expert_level = 0
    list_occ = None
        .type = floats(size_min=1, value_min=0, value_max=1)
        .help = List of occupancies to test (fractional). Will overwrite low_occ, high_occ and steps if defined.
        .expert_level = 0
    }
scaling{
    b_scaling = no isotropic *anisotropic 
        .type = choice(multi=False)
        .help = B-factor scaling for scaling triggered data vs reference data using scaleit.
        .expert_level = 0
    high_resolution = None
        .type = float
        .help = High resolution for scaling triggered data vs reference data using scaleit (Angstrom). Will only be used if high resolution of the input data files extends to this value. This only implies scaling, the data will not be cut. If not specified, then input.high_resolution will be used.
        .expert_level = 3
    low_resolution = None
        .type = float
        .help = Low resolution for scaling triggered data vs reference data using scaleit (Angstrom). Will only be used if low resolution of the input data files extends to this value. This only implies scaling, the data will not be cut. If not specified, then input.low_resolution will be used.
        .expert_level = 3
    }
f_and_maps{
    fofo_type = *qfofo fofo kfofo
        .type = choice(multi=False)
        .help = Calculate q-weighted, non-weighted or k-weighted Fourier difference (Fo-Fo) map.
        .expert_level = 1
    kweight_scale = 0.05
        .type = float(value_min=0, value_max=1)
        .help = scale factor for structure factor difference outlier rejection in k-weigting scheme.
        .expert_level = 3
    f_extrapolated_and_maps = *qfextr fextr kfextr qfgenick fgenick kfgenick qfextr_calc fextr_calc kfextr_calc
        .type = choice(multi=True)
        .help = Extrapolated structure factors and map types: qFextr, kFextr, Fextr: (q/k-weighted)-Fextr structure factors and maps by Coquelle method (Fextr= alpha*(Fobs,triggered-Fobs,reference)+Fobs,triggered, map 2mFextr|-D|Fcalc|, phi_model). qFgenick, kFgenick, Fgenick: (q/k-weighted)-Fextr structure factors and maps by Genick method (|Fextr|= alpha*(|Fobs,triggered|-|Fobs,reference|)+|Fobs,triggered|, map: m|Fextr|, phi_model). qFextr_calc, kFextr_calc, Fextr_calc: (q-weighted)-Fextr structure factors and maps by Fcalc method (|Fextr|= alpha*(|Fobs,triggered|-|Fobs,reference|)+|Fcalc|, map" 2m|Fextr|-D|Fcalc|, phi_model).
        .expert_level = 0
    all_maps = False
        .type = bool
        .help = Calculate all extrapolated structure factors and maps.
        .expert_level = 0
    only_qweight = False
        .type = bool
        .help = Calculate all extrapolated structure factors and maps with q-weighting.
        .expert_level = 0
    only_kweight = False
        .type = bool
        .help = Calculate all extrapolated structure factors and maps with k-weighting.
        .expert_level = 0
    only_no_weight = False
        .type = bool
        .help = Calculate all extrapolated structure factors and maps without q/k-weighting.
        .expert_level = 0
    fast_and_furious = False
        .type = bool
        .help = Run fast and furious (aka without supervision). Will only calculate qFextr and associated maps and run refinement with finally with derived alpha/occupancy. Default parameters will be used for fofo_type and negative_and_missing. Usefull for a first quick evaluation.
        .expert_level = 0
    negative_and_missing = *truncate_and_fill truncate_no_fill fref_and_fill fref_no_fill fcalc_and_fill fcalc_no_fill keep_and_fill keep_no_fill reject_and_fill reject_no_fill zero_and_fill zero_no_fill fill_missing no_fill absolute_and_fill absolute_no_fill
        .type = choice(multi=False)
        .help = Handling of negative and missing extrapolated structure factor amplitudes (ESFAs) Note that this will not be applied on the Fourier difference map. If selected, filling of missing reflections is only carried on maps of the 2mFextr-DFcalc type. This parameters is NOT applicable for (q/k)Fgenick because negative reflections are rejected anyway. For refinement, default phenix.refine or refmac handling of negative/missing reflections is applied. keep_no_fill maps will be calculated in addition in all cases. keep_and_fill and keep_no_fill replace the old fill_missing and no_fill arguments which will become invalid keywords in future Xtrapol8 versions. Please check the manual for more information.
        .expert_level = 2
    }
map_explorer{
    peak_integration_floor = 3.5
        .type = float
        .help = Floor value for peak integration (sigma). Peaks will be integrated from their maximum value towards this lower bound to avoid integration of noise.
        .expert_level = 0
    peak_detection_threshold = 4.0
        .type = float
        .help = Peak detection threshold (sigma). Only peaks with an absolute value equal or above this value will be integrated.
    radius = None
        .type = float
        .help = Maximum radius (A) to allocate a density blob to a protein atom in map explorer. Resolution will be used if not specified.
        .expert_level = 0
    z_score = 2.0
        .type = float
        .help = Z-score to determine residue list with only highest peaks.
        .expert_level = 0
    use_occupancy_from_distance_analysis = False
        .type = bool
        .help = Use occupancy as estimated by the distance analysis method (only in calm_and_curious mode) instead of the differrence map analysis. This keyword will become obselete in future Xtrapol8 versions, use occupancy_estimation choice instead.
        .expert_level = 1
    occupancy_estimation = *difference_map_maximization difference_map_PearsonCC distance_analysis
        .type = choice(multi=False)
        .help = Select a main method for the occupancy estimation in Xtrapol8. Take care that the distance_analysis method can only be used in calm_and_curious mode. This keyword replaces the use_occupancy_from_distance_analysis keyword which will become obsolete in future Xtrapol8 versions.
    }
refinement{
    run_refinement = True
    .type = bool
    .help = Run the automatic refinements. Setting this parameter to False can be useful when a manual intervention is required before running the refinements. The Refiner.py script can be used to run the refinements and subsequent analysis afterwards.
    .expert_level = 1
    reciprocal_space = *phenix.refine refmac5
        .type = choice(multi=False)
        .help = Program for reciprocal space refinement: phenix.refine or refmac5. This keyword replaces the use_refmac_instead_of_phenix to make a choice of refinement program. 
    real_space = *phenix.real_space_refine coot
        .type = choice(multi=False)
        .help = Program for real space refinement: phenix.real_space_refine or COOT. This keyword replaces the use_refmac_instead_of_phenix to make a choice of refinement program.
    use_refmac_instead_of_phenix = False
        .type = bool
        .help = use Refmac for reciprocal space refinement and COOT for real-space refinement instead of phenix.refine and phenix.real_space_refine. This keyword will becone obsolete in furure Xtrapol8 versions, use reciprocal_space and real_space instead.
        .expert_level = 0
    phenix_keywords{
            real_space_refine{
            cycles = 5
            .type = int
            .help = Number of refinement cycles for real space refinement.
            .expert_level = 0
            }
        target_weights{
            wxc_scale = 0.5
            .type = float
            .help = phenix.refine refinement.target_weights.wxc_scale.
            .expert_level = 2
            wxu_scale = 1.0
            .type = float
            .help = phenix.refine refinement.target_weights.wxu_scale.
            .expert_level = 2
            weight_selection_criteria{
                bonds_rmsd = None
                .type = float
                .help = phenix.refine refinement.target_weights.weight_selection_criteria.bonds_rmsd.
                .expert_level = 3
                angles_rmsd = None
                .type = float
                .help = phenix.refine refinement.target_weights.weight_selection_criteria.angles_rmsd.
                .expert_level = 3
                r_free_minus_r_work = None
                .type = float
                .help = phenix.refine refinement.target_weights.weight_selection_criteria.r_free_minus_r_work.
                .expert_level = 3
                }
            }
        refine{
            strategy = *individual_sites individual_sites_real_space rigid_body *individual_adp group_adp tls occupancies group_anomalous
            .type = choice(multi=True)
            .help = phenix.refine refinement.refine.strategy.
            .expert_level = 1
            }
        main{
            cycles = 5
            .type = int
            .help = Number of refinement macro cycles for reciprocal space refinement.
            .expert_level = 0
            ordered_solvent = False
            .type = bool
            .help = Add and remove ordered solvent during reciprocal space refinement (refinement.refine.main.ordered_solvent).
            .expert_level = 0
            simulated_annealing = False
            .type = bool
            .help = Simulated annealing during refinement.
            .expert_level = 1
            }
        simulated_annealing{
            start_temperature = 5000
            .type = float
            .help = start temperature for simulated annealing.
            .expert_level = 2
            final_temperature = 300
            .type = float
            .help = final temperature for simulated annealing.
            .expert_level = 2
            cool_rate = 100
            .type = float
            .help = cool rate for simulated annealing.
            .expert_level = 2
            mode = every_macro_cycle *second_and_before_last once first first_half
            .type = choice(multi=False)
            .help = simulated annealing mode.
            .expert_level = 2
            }
        map_sharpening{
            map_sharpening = False
            .type = bool
            .help = phenix map sharpening.
            .expert_level = 1
            }
        additional_reciprocal_space_keywords = None
            .type = str
            .multiple = True 
            .help = Additional phenix.refine keywords which cannot be altered via included options (e.g. ncs_search.enabled=True).
            .expert_level = 2
        additional_real_space_keywords = None
            .type = str
            .multiple = True 
            .help = Additional phenix_real_space.refine keywords which cannot be altered via included options (e.g. ncs_constraints=False).
            .expert_level = 2
        density_modification{
            density_modification = False
            .type = bool
            .help = use dm (ccp4) for density modification.
            .expert_level = 2
            combine = PERT *OMIT
            .type = choice(multi=False)
            .help = dm combine mode.
            .expert_level = 2
            cycles = 3
            .type = int
            .help = number of dm cycles (ncycle keyword). Use a lot of cycles when combine=PERT and only few cycles when combine=OMIT.
            .expert_level = 2
            }
        }
    refmac_keywords{
        target_weights{
            weight = *AUTO MATRIx
            .type = choice(multi=False)
            .help = refmac WEIGHT.
            .expert_level = 1
            weighting_term = 0.2
            .type = float
            .help = refmac weighting term in case of weight matrix.
            .expert_level = 2
            experimental_sigmas = *NOEX EXPE
            .type = choice(multi=False)
            .help = refmac use experimental sigmas to weight Xray terms.
            .expert_level = 2
            }
        restraints{
            jelly_body_refinement = False
            .type = bool
            .help = run refmac ridge regression, also known as jelly body jelly body refinement. Slow refinement convergence, so take at least 50 refinement cycles.
            .expert_level = 1
            jelly_body_sigma = 0.03
            .type = float
            .help = sigma parameter in case of jelly body refinement ('RIDG DIST SIGM' parameter).
            .expert_level = 2
            jelly_body_additional_restraints = None
            .type = str
            .multiple = True
            .help = additional jelly body parameters (will be added to keyword 'RIDG').
            .expert_level = 2
            external_restraints = None
            .type = str
            .multiple = True
            .help = refmac external restraints (will be added to keyword 'external', e.g. 'harmonic residues from 225 A to 250 A atom CA sigma 0.02').
            .expert_level = 2
            }
        refine{
            type = *RESTrained UNREstrained RIGId
            .type = choice(multi=False)
            .help = refmac refinement type refinement.
            .expert_level = 1
            TLS = False
            .type = bool
            .help = tls refinement before coordinate and B-factor refinement.
            .expert_level = 1
            TLS_cycles = 20
            .type = int
            .help = number of TLS cycles in case of TLS refinement.
            .expert_level = 2
            bfac_set = 30
            .type = float
            .help = reset individual B-factors to constant value before running TLS. Will only be applied in case TLS is run.
            .expert_level = 2
            twinning = False
            .type = bool
            .help = do refmac twin refinement.
            .expert_level = 1
            Brefinement = OVERall *ISOTropic
            .type = choice(multi=False)
            .help = refmac B-factor refinement.
            .expert_level = 1
            cycles = 20
            .type = int
            .help = Number of refinement cycles for reciprocal space refinement.
            .expert_level = 0
            }
        map_sharpening{
            map_sharpening = False
            .type = bool
            .help = refmac map sharpening.
            .expert_level = 1
            }
        additional_refmac_keywords = None
            .type = str
            .multiple = True 
            .help = Additional refmac keywords which cannot be altered via included options (e.g. ncsr local).
            .expert_level = 2        
        density_modification{
            density_modification = False
            .type = bool
            .help = use dm for density modification.
            .expert_level = 2
            combine = PERT *OMIT
            .type = choice(multi=False)
            .help = dm combine mode.
            .expert_level = 2
            cycles = 3
            .type = int
            .help = number of dm cycles (ncycle keyword). Use a lot of cycles when combine=PERT and only few cycles when combine=OMIT
            .expert_level = 2
            }
        }
    }
output{
    outdir = None
        .type = str
        .help = Output directory. 'Xtrapol8' be used if not specified.
        .expert_level = 0
    outname = None
        .type = str
        .help = Prefix or suffix for output files. The prefix of triggered_mtz will be used if not specified.
        .expert_level = 0
    generate_phil_only = False
        .type = bool
        .help = Generate input phil-file and quit.
        .expert_level = 0
    generate_fofo_only = False
        .type = bool
        .help = Stop Xtrapol8 after generation of Fourier Difference map.
        .expert_level = 0
    open_coot = True
        .type = bool
        .help = Automatically open COOT at the end.
        .expert_level = 0
    ddm_scale = None
        .type = float
        .help = The ddm colors will range from -scale to +scale.
        .expert_level = 2
    GUI = False
        .type = bool
        .help = Xtrapol8 launched from GUI. In order to work correctly, this should never be manually changed.
        .expert_level = 3
    }
""", process_includes=True)
