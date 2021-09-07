"""
Main script to run Xtrapol8

authors and contact information
-------
Elke De Zitter - elke.de-zitter@ibs.fr
Nicolac Coquelle - nicolas.coquelle@esrf.fr
Thomas Barends - Thomas.Barends@mpimf-heidelberg.mpg.de
Jacques Philippe Colletier - jacques-Philippe.colletier@ibs.fr
Please start your mail with [X8]

licence
-------
Some free licence
Xtrapol8 requires a Phenix and CCP4 installation with proper licence.

usage
-------
Fextr.py is the main script to be called with phenix.python

Run without any argument to see all options and explanation:
>>> phenix.python <wherever>/Fextr.py
Parameters can be added using an input file or via command line

example using input file (preferable)
-------
1) Change the nano Xtrapol8.phil using your favourite editor, e.g.
>>> nano trapol8.phil
2) Run Xtrapol8
>>> phenix.python <wherever>/Fextr.py trapol8.phil

example using command line only
-------
1) Run Xtrapol8 with all your arguments
>>> phenix.python <wherever>/Fextr.py input.reference_mtz=hiephiep.mtz input.triggered_mtz=hieperdepiep.mtz input.reference_pdb=hoera.pdb input.additional_files=jeej.cif input.additional_files=another.cif occupancies.list_occ=0.1,0.3,0.5 f_and_maps.f_extrapolated_and_maps=qfextr,qfgenick map_explorer.threshold=3.5 map_explorer.peak=4 output.outdir=fancy_party

example using input file and command line
-------
1) Change the nano Xtrapol8.phil using your favourite editor, e.g.
>>> nano trapol8.phil
2) Xtrapol8 with additional arguments. The order of arguments determines how paramters will be overwriten:
>>> phenix.python <wherever>/Fextr.py trapol8.phil refinement.phenix_keywords.refine.cycles=3

-------

TODO:
- clean up Fextr_utils: remove unnecessary functions
- Resolve problem with 90% multiplicity estimation
- option to automatically create cif with phenix.elbow
- Map explorer analysis: should be working for DNA.
- GUI: warning message if people want to run crazy stuff that they better run on a powerfull machine in command line mode
- update example phil and manual
- Suggest resolution cutoff based on Riso and CCiso, write resolution suggestions based on true completeness and <F/SIG(F)> to logfile
- add comments to various steps of script (especially the weird steps)
- nucleic acid compatibility
- keep anomalous data
- If model hasen't been refined with the mtz_ref,we should first run a rigid body refinement and use the output pdb from that
- make an object to store the results in a clean and transparant way and get rid of lists which may mess up the analysis if refinements have failed
- check if refinement has finished and if not, check for .eff file to retry.
- ss-plot, share ylabel for graphs in orderto avoid overlap of the ylabels (which are identical anyway)
- Refmac BREF no
- option to use maps in absolute electron density instead of sigma-scaled
- absolute q-weighting (no normalization by devision of the average)
- usefullness of the sig(FoFo)_vs_sig(Fextr) plot
- check scaling no with data from XSCALE
"""
from __future__ import division, print_function
import re
import os
import sys
import random
import subprocess
import shutil
from select import select
from datetime import datetime
import numpy as np
from iotbx.file_reader import any_file
from iotbx import symmetry
from iotbx.pdb import hierarchy
from cctbx.array_family import flex
from cctbx import maptbx, miller, crystal, xray
import iotbx.phil
import iotbx.map_tools
from libtbx.utils import Usage
import mmtbx.f_model
import mmtbx.map_tools
import mmtbx.maps.utils
from mmtbx import utils
from cctbx import sgtbx
from iotbx import pdb
from mmtbx.scaling.matthews import p_vm_calculator
import multiprocessing
from functools import partial

#sys.path.append("/Users/edezitter/Scripts/Fextrapolation")

from calculate_q import calculate_q, outlier_rejection_only
from calculate_k import calculate_k
from map_explorer import map_explorer
from map_explorer_analysis import Map_explorer_analysis
from ccp4_scaleit import run_scaleit
from plotalpha import plotalpha
import map_tools_fomsource
import map_tools_Millerset
import ccp4_refmac
import phenix_refinements
from column_extraction import Column_extraction, Extrapolated_column_extraction
from pymol_visualization import Pymol_visualization, Pymol_movie
from ddm import Difference_distance_analysis
from distance_analysis import *
from Fextr_utils import *
from Fextr import SymManager, DataHandler, FobsFobs, Fextrapolate, Filesandmaps
import main_X8
import main_JK

import JK_utils
from JK_utils import print_terminal_and_log as print_in_T_and_log
import Parameters as Para
from Parameters import Parameters
import JK_results
import JK_simple_refinement

#TODO: Check if the input arguments for Crystfel are not overwritten if for example an other hkl file is given in 'other_stats_compare_hkl'. This is linked to the problem where we can not check if the additional parameters for the Crystfel programs are correct.
#TODO: If UC and SG are different for reference and triggered, what happens?
#TODO: Create graphs for the normal distribution of the UC parameters if the means are taken
#TODO: Check/Improve outputs printed in log files and terminal
#TODO: Add the parameters changed from the phil to the output phil: done in Parameters, needs to be checked, commented as 'return in output phil'
#TODO: Be carefull with the unique axis given, if it is wrong, f2mtz will not work!
#TODO: Add possibility for many triggered mtz files or stream files: done

master_phil = iotbx.phil.parse("""
options{
    programs = run_JackKnife run_Xtrapol8
            .type = choice(multi=True)
            .help = Program to run
            .expert_level = 0
    processors = 10
        .type = int
        .help = Number of processors to be used for the program
        .expert_level = 1
    }
JackKnife{
    input{
        repeats = None
            .type = int
            .help = Number of times JackKnife will be repeated
            .expert_level = 0
        fraction = 0.9
            .type = float(value_min=0, value_max=1)
            .help = Percentage of images used for JackKnife (fractional)
            .expert_level = 0
        directory_crystfel_programs = None
            .type = path
            .help = path where the CrystFEL programs like process_hkl are found (example: /usr/local/bin)
            .expert_level = 0
        }
    Reference_state{
        reference_stream_file = None
            .type = path
            .help = File with images to process for the off state or only state in stream format
            .expert_level = 0            
        use_UC_and_SG_from_pdb = False
            .type = bool
            .help = use the unit cell and space group from pdb file given
            .expert_level = 0
        unit_cell{
            a = None
                .type = float
                .help = lattice a axis length (A)
                .expert_level = 0
            b = None
                .type = float
                .help = lattice b axis length (A)
                .expert_level = 0
            c = None
                .type = float
                .help = lattice c axis length (A)
                .expert_level = 0
            alpha = None
                .type = float
                .help = lattice alpha angle (degree)
                .expert_level = 0
            beta = None
                .type = float
                .help = lattice beta angle (degree)
                .expert_level = 0
            gamma = None
                .type = float
                .help = lattice gamma angle (degree)
                .expert_level = 0
            }    
            spacegroup = None
                .type = str
                .help = space group of the processed crystal
                .expert_level = 0
            unique_axis = *default a b c
                .type = choice(multi=False)
                .help = symmetry axis
                .expert_level = 0
        }
    Triggered_state{
        triggered_stream_file = None
            .type = path
            .multiple = True
            .help = File with images to process for the on state in stream format
            .expert_level = 0            
        use_UC_and_SG_from_pdb = False
            .type = bool
            .help = use the unit cell and space group from reference pdb file given
            .expert_level = 0
        unit_cell{
            a = None
                .type = float
                .help = lattice a axis length (A)
                .expert_level = 0
            b = None
                .type = float
                .help = lattice b axis length (A)
                .expert_level = 0
            c = None
                .type = float
                .help = lattice c axis length (A)
                .expert_level = 0
            alpha = None
                .type = float
                .help = lattice alpha angle (degree)
                .expert_level = 0
            beta = None
                .type = float
                .help = lattice beta angle (degree)
                .expert_level = 0
            gamma = None
                .type = float
                .help = lattice gamma angle (degree)
                .expert_level = 0
            }
            spacegroup = None
                .type = str
                .help = space group of the processed crystal
                .expert_level = 0
            unique_axis = *default a b c
                .type = choice(multi=False)
                .help = symmetry axis
                .expert_level = 0
        }
    Scaling_and_merging{
        algorithm = *process_hkl partialator
            .type = choice(multi=False)
            .help = Method for intensity merging
            .expert_level = 0
        process_hkl{
            other_process_hkl = --max-adu=17000 --scale
                .type = str
                .help = other inputs for the process_hlk merging of images. For more infromation go to https://www.desy.de/~twhite/crystfel/manual.html (written like inputs for Crystfel, example: --min-snr=1)
                .expert_level = 0
            }
        partialator{
            other_partialator = --iterations=1 --model=unity
                .type = str
                .help = other inputs for the process_hlk merging of images.
                .expert_level = 0
            }
        }
    Statistics{
        other_stats_compare_hkl = None
            .type = str
            .help = other inputs for the compare_hkl getting statistics.
            .expert_level = 0
        }
    }
Xtrapol8{
    input{
        reference_mtz = None
            .type = path
            .help = Reference data in mtz or mmcif format (merged).
            .expert_level = 0
        triggered_mtz = None
            .type = path
            .multiple = True
            .help = Triggered data in mtz or mmcif format (merged).
            .expert_level = 0
        reference_pdb = None
            .type = path
            .help = Reference coordinates in pdb or mmcif format. (in former versions this was called model_pdb)
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
            .help = Lowest occupancy to test (fractional)
            .expert_level = 0
        high_occ = 0.5
            .type = float(value_min=0, value_max=1)
            .help = Highest occupancy to test (fractional)
            .expert_level = 0
        steps = 3
            .type = int
            .help = Amount of equaly spaced occupancies to be tested
            .expert_level = 0
        list_occ = None
            .type = floats(size_min=1, value_min=0, value_max=1)
            .help = List of occupancies to test (fractional). Will overwrite low_occ, high_occ and steps if defined
            .expert_level = 0
        }
    scaling{
        b_scaling = no isotropic *anisotropic 
            .type = choice(multi=False)
            .help = B-factor scaling for scaling triggered data vs reference data. Cannot be used for reference data with fcalc when using mmtbx.fmodel.manager.
            .expert_level = 0
        }
    f_and_maps{
        fofo_type = *qfofo fofo kfofo
            .type = choice(multi=False)
            .help = Calculate q-weighted or non-q-weighted Fo-Fo difference map. Q-weighted is highly recommended. K-weighting is under development
            .expert_level = 1
        kweight_scale = 0.05
            .type = float(value_min=0, value_max=1)
            .help = scale factor for structure factor difference in k-weigting scheme (for calculation of kfofo)
            .expert_level = 3
        f_extrapolated_and_maps = *qfextr fextr kfextr qfgenick fgenick kfgenick qfextr_calc fextr_calc kfextr_calc
            .type = choice(multi=True)
            .help = Extrapolated structure factors and map types: qFextr, kFextr, Fextr: (q/k-weighted)-Fextr structure factors and maps by Coquelle method (Fextr= alpha*(Fobs,triggered-Fobs,reference)+Fobs,triggered, map 2mFextr|-D|Fcalc|, phi_model). qFgenick, kFgenick, Fgenick: (q/k-weighted)-Fextr structure factors and maps by Genick method (|Fextr|= alpha*(|Fobs,triggered|-|Fobs,reference|)+|Fobs,triggered|, map: m|Fextr|, phi_model). qFextr_calc, kFextr_calc, Fextr_calc: (q-weighted)-Fextr structure factors and maps by Fcalc method (|Fextr|= alpha*(|Fobs,triggered|-|Fobs,reference|)+|Fcalc|, map" 2m|Fextr|-D|Fcalc|, phi_model).
            .expert_level = 0
        all_maps = False
            .type = bool
            .help = Calculate all extrapolated structure factors and maps
            .expert_level = 0
        only_qweight = False
            .type = bool
            .help = Calculate all extrapolated structure factors and maps with q-weighting
            .expert_level = 0
        only_kweight = False
            .type = bool
            .help = Calculate all extrapolated structure factors and maps with q-weighting
            .expert_level = 0
        only_no_weight = False
            .type = bool
            .help = Calculate all extrapolated structure factors and maps without q/k-weighting
            .expert_level = 0
        fast_and_furious = False
            .type = bool
            .help = Run fast and furious (aka without supervision). Will only calculate qFextr and associated maps, use highest peaks for alpha/occupancy determination (alpha/occupancy will be nonsense if map_explorer parameters being bad), run refinement with finally with derived alpha/occupancy, use truncate_and_fill for negative and missing handling. Usefull for a first quick evaluation.
            .expert_level = 0
        negative_and_missing = *truncate_and_fill truncate_no_fill fref_and_fill fref_no_fill fcalc_and_fill fcalc_no_fill fill_missing no_fill reject_and_fill reject_no_fill zero_and_fill zero_no_fill
            .type = choice(multi=False)
            .help = Handling of negative and missing extrapolated reflections (note that this will not be applied on FoFo difference maps). Please check the manual for more information. This parameters is NOT applicable for (q)Fgenick because negative reflections are rejected anyway. For refinement, default phenix.refine or refmac handling of negative/missing reflections is applied.
            .expert_level = 2
        }
    map_explorer{
        threshold = 3.5
            .type = float
            .help = Integration threshold (in sigma) 
            .expert_level = 0
        peak = 4.0
            .type = float
            .help = Peak detection threshold (sigma)
        radius = None
            .type = float
            .help = Maximum radius (A) to allocate a density blob to a protein atom in map explorer. Resolution will be used if not specified.
            .expert_level = 0
        z_score = 2.0
            .type = float
            .help = Z-score to determine residue list with only highest peaks
            .expert_level = 0
        use_occupancy_from_distance_analysis = False
            .type = bool
            .help = Use occupancy from determination based on the differences between reference_pdb and real-space refined model (only in calm_and_curious mode) instead of map explorer
            .expert_level = 1
        }
    refinement{
        run_refinement = True
        .type = bool
        .help = Run the automatic refinements. Setting this parameter to False can be useful when a manual intervention is required before running the refinements. The Refiner.py script can be used to run the refinements and subsequent analysis afterwards.
        .expert_level = 1
        use_refmac_instead_of_phenix = False
            .type = bool
            .help = use Refmac for reciprocal space refinement and COOT for real-space refinement instead of phenix.refine and phenix.real_space_refine
            .expert_level = 0
        phenix_keywords{
            target_weights{
                wxc_scale = 0.5
                .type = float
                .help = see phenix.refine refinement.target_weights.wxc_scale
                .expert_level = 2
                wxu_scale = 1.0
                .type = float
                .help = phenix.refine refinement.target_weights.wxu_scale
                .expert_level = 2
                weight_selection_criteria{
                    bonds_rmsd = None
                    .type = float
                    .help = phenix.refine refinement.target_weights.weight_selection_criteria.bonds_rmsd
                    .expert_level = 3
                    angles_rmsd = None
                    .type = float
                    .help = phenix.refine refinement.target_weights.weight_selection_criteria.angles_rmsd
                    .expert_level = 3
                    r_free_minus_r_work = None
                    .type = float
                    .help = phenix.refine refinement.target_weights.weight_selection_criteria.r_free_minus_r_work
                    .expert_level = 3
                    }
                }
            refine{
                strategy = *individual_sites individual_sites_real_space rigid_body *individual_adp group_adp tls occupancies group_anomalous
                .type = choice(multi=True)
                .help = see phenix.refine refinement.refine.strategy
                .expert_level = 1
                }
            main{
                cycles = 5
                .type = int
                .help = Number of refinement macro cycles for reciprocal space refinement
                .expert_level = 0
                ordered_solvent = False
                .type = bool
                .help = Add and remove ordered solvent during reciprocal space refinement
                .expert_level = 0
                simulated_annealing = False
                .type = bool
                .help = Simulated annealing during refinement
                .expert_level = 1
                }
            simulated_annealing{
                start_temperature = 5000
                .type = float
                .help = start temperature for simulated annealing
                .expert_level = 2
                final_temperature = 300
                .type = float
                .help = final temperature for simulated annealing
                .expert_level = 2
                cool_rate = 100
                .type = float
                .help = cool rate for simulated annealing
                .expert_level = 2
                mode = every_macro_cycle *second_and_before_last once first first_half
                .type = choice(multi=False)
                .help = simulated annealing mode
                .expert_level = 2
                }
            map_sharpening{
                map_sharpening = False
                .type = bool
                .help = phenix map sharpening
                .expert_level = 1
                }
            real_space_refine{
                cycles = 5
                .type = int
                .help = Number of refinement cycles for real space refinement
                .expert_level = 0
                }
            density_modification{
                density_modification = False
                .type = bool
                .help = use dm (ccp4) for density modification
                .expert_level = 2
                combine = *PERT OMIT
                .type = choice(multi=False)
                .help = dm combine mode
                .expert_level = 2
                cycles = 10
                .type = int
                .help = number of dm cycles (ncycle keyword). Use only few cycles in case of combine=OMIT
                .expert_level = 2
                }
            }
        refmac_keywords{
            target_weights{
                weight = *AUTO MATRIx
                .type = choice(multi=False)
                .help = refmac WEIGHT
                .expert_level = 1
                weighting_term = 0.2
                .type = float
                .help = refmac weighting term in case of weight matrix
                .expert_level = 2
                experimental_sigmas = *NOEX EXPE
                .type = choice(multi=False)
                .help = refmac use experimental sigmas to weight Xray terms
                .expert_level = 2
                }
            restraints{
                jelly_body_refinement = False
                .type = bool
                .help = run refmac ridge regression, also known as jelly body jelly body refinement. Slow refinement convergence, so take at least 50 refinement cycles.
                .expert_level = 1
                jelly_body_sigma = 0.03
                .type = float
                .help = sigma parameter in case of jelly body refinement ('RIDG DIST SIGM' parameter)
                .expert_level = 2
                jelly_body_additional_restraints = None
                .type = str
                .multiple = True
                .help = additional jelly body parameters (will be added to keyword 'RIDG')
                .expert_level = 2
                external_restraints = None
                .type = str
                .multiple = True
                .help = refmac external restraints (will be added to keyword 'external', e.g. 'harmonic residues from 225 A to 250 A atom CA sigma 0.02')
                .expert_level = 2
                }
            refine{
                type = *RESTrained UNREstrained RIGId
                .type = choice(multi=False)
                .help = refmac refinement type refinement
                .expert_level = 1
                TLS = False
                .type = bool
                .help = tls refinement before coordinate and B-factor refinement
                .expert_level = 1
                TLS_cycles = 20
                .type = int
                .help = number of TLS cycles in case of TLS refinement
                .expert_level = 2
                bfac_set = 30
                .type = float
                .help = reset individual B-factors to constant value before running TLS. Will only be applied in case TLS is run
                .expert_level = 2
                twinning = False
                .type = bool
                .help = do refmac twin refinement
                .expert_level = 1
                Brefinement = OVERall *ISOTropic
                .type = choice(multi=False)
                .help = refmac B-factor refinement
                .expert_level = 1
                cycles = 20
                .type = int
                .help = Number of refinement cycles for reciprocal space refinement
                .expert_level = 0
                }
            map_sharpening{
                map_sharpening = False
                .type = bool
                .help = refmac map sharpening
                .expert_level = 1
                }
            density_modification{
                density_modification = False
                .type = bool
                .help = use dm for density modification
                .expert_level = 2
                combine = *PERT OMIT
                .type = choice(multi=False)
                .help = dm combine mode
                .expert_level = 2
                cycles = 10
                .type = int
                .help = number of dm cycles (ncycle keyword). Use only few cycles in case of combine=OMIT
                .expert_level = 2
                }
            }
        }
    output{
        generate_phil_only = False
            .type = bool
            .help = Generate input phil-file and quit.
            .expert_level = 0
        generate_fofo_only = False
            .type = bool
            .help = Stop Xtrapol8 after generation of Fourier Difference map
            .expert_level = 0
        open_coot = True
            .type = bool
            .help = Automatically open COOT at the end.
            .expert_level = 0
        ddm_scale = 1.5
            .type = float
            .help = The ddm colors will range from -scale to +scale.
            .expert_level = 2
        }
    }
output{
    outdir = None
        .type = str
        .help = Output directory. Current directory directory will be used if not specified.
        .expert_level = 0
    outname = None
        .type = str
        .multiple = True
        .help = Output prefix. Prefix of triggered_mtz will be used if not specified.
        .expert_level = 0
    GUI = False
        .type = bool
        .help = Xtrapol8 launched from GUI.
        .expert_level = 3
    }
""", process_includes=True)

def run(args):

    print('################################################################################\nLAUNCH XTRAPOL8\n################################################################################')
    print('LAUNCH XTRAPOL8\n==============================================================')
    # If no input, show complete help, should be changed in order to give help depending on the attribute level
    if len(args) == 0 :
        master_phil.show(attributes_level=1)
        raise Usage("phenix.python main_JK_X8.py + [.phil] + [arguments]\n arguments only overwrite .phil if provided last")

    print('CREATING LOG FILE\n==============================================================')
    # Generate log-file. Needs to be created before the output directory is created and to be a global parameter in order to be easily used in all classes and functions
    #create log file and get global log = open(logname, "w")
    log_dir = os.getcwd() #directory where the log file will be created
    logdir, log = JK_utils.create_log(log_dir, global_log=True) #full directory of the log file

    print_in_T_and_log('EXTRACTING INPUT PHIL FILE\n==============================================================')
    # Extract input from inputfile and command line
    argument_interpreter = master_phil.command_line_argument_interpreter(home_scope="options")
    input_objects = iotbx.phil.process_command_line_with_files(
        args=args,
        master_phil=master_phil
    )
    params = input_objects.work.extract()
#Elke    # modified_phil = master_phil.format(python_object=params)

    print_in_T_and_log('GETTING PROGRAM TO EXECUTE\n==============================================================')
    P = Parameters(params) #init Class to get parameters
    P.get_parameters_JK_X8() #get parameters common to JK and Xtrapol8

    print_in_T_and_log('CHECKING PROGRAMS\n==============================================================')
    #Check existance and execution of necessary programs
    if P.run_JackKnife:
        Para.check_programs_JK()
    if P.run_Xtrapol8:
        Para.check_programs_X8(params)

    JK_utils.print_terminal_and_log('GETTING PARAMETERS\n==============================================================')
    #get all input parameters from phil and stream file

    #Order of getting parameters is IMPORTANT!!
    if P.run_Xtrapol8:
         P.get_parameters_X8()

    if P.run_JackKnife:
        P.get_parameters_JK()

    ################################################################

    # Add all arguments to log-file
    print_in_T_and_log('-----------------------------------------\nNON DEFAULT ARGUMENTS\n-----------------------------------------')

    modified_phil = master_phil.format(python_object=params)
#Elke    # modified_phil.show(out=log)
#Elke    # get the differences with the default values and only show these in the log-file
    diff_phil = master_phil.fetch_diff(source=modified_phil)
    diff_phil.show()
    diff_phil.show(out=log)

    ################################################################

    #Run JackKnife and Xtrapol8
    #Run JackKnife
    if P.run_JackKnife:
    #get parameters and move to directory
        index=0 #init the index
        for stream_file_on in P.stream_file_on: #for each stream file (multiple stream files possible)
            P.get_parameters_multiple_JK_X8_output(params, index)  # getting output parameters common to JK and X8: output and outname
            # Change directory to output directory
            outdir = P.outdir
            outname = P.outname
            print('\nthe output directory is %s' % (outdir))
            # Move log file to output directory
            if os.path.isfile(logdir):
                shutil.move(logdir, logdir.replace(log_dir, outdir))
            os.chdir(outdir)

            P.get_parameters_JK_triggered_stream_file(stream_file_on) #getting the parameters of the triggered stream file

    #run Jackknife only
            print_in_T_and_log('\n################################################################################\nLAUNCH JACK KNIFE\n################################################################################')

    #run Jackknife

            print_in_T_and_log('JACK KNIFE FOR OFF STATE\n==============================================================')
            #run JackKnife and get list of [directory where to find mtz file, mtz file complete directory, nb_JK] for the off state
            if P.fraction == 1: #if fraction wanted is 1, no total needed (same thing)#TODO: if the fraction =1, JK is ran as many times as 'repeats', change that into 1? If so, go to main_JK>Run_JK and change 'P.repeats'=1
                mtzoutdirs_off,_ = main_JK.run_JK(P, outdir, P.stream_file_off, P.stream_file_name_off,
                                                P.n_frames_to_keep_off, P.system_off, P.pointgroup_off,
                                                P.unique_axis_off, P.a_off, P.b_off, P.c_off, P.alpha_off, P.beta_off,
                                                P.gamma_off, P.spacegroup_off, log, state='off', total=False)
            else:
                mtzoutdirs_off, mtzoutdirs_off_total = main_JK.run_JK(P, outdir, P.stream_file_off, P.stream_file_name_off, P.n_frames_to_keep_off, P.system_off, P.pointgroup_off, P.unique_axis_off, P.a_off, P.b_off, P.c_off, P.alpha_off, P.beta_off, P.gamma_off, P.spacegroup_off, log, state='off', total=True)
            print_in_T_and_log('JACK KNIFE FOR ON STATE\n==============================================================')

            #run JackKnife and get list of [directory where to find mtz file, mtz file complete directory, nb_JK] for the on state
            if P.fraction == 1:
                mtzoutdirs_on, _ = main_JK.run_JK(P, outdir, stream_file_on, P.stream_file_name_on,
                                               P.n_frames_to_keep_on, P.system_on, P.pointgroup_on, P.unique_axis_on,
                                               P.a_on, P.b_on, P.c_on, P.alpha_on, P.beta_on, P.gamma_on,
                                               P.spacegroup_on, log, state='on', total=False)
            else:
                mtzoutdirs_on, mtzoutdirs_on_total = main_JK.run_JK(P, outdir, stream_file_on, P.stream_file_name_on, P.n_frames_to_keep_on, P.system_on, P.pointgroup_on, P.unique_axis_on, P.a_on, P.b_on, P.c_on, P.alpha_on, P.beta_on, P.gamma_on, P.spacegroup_on, log, state='on', total=True)

            ################################################################

        #create the list mtzoutdirs_dir_off_on with [outdirJKX8, mtz_off, mtz_on, newlog, nb_JK, outname]
            def create_multiprocessing_list(mtzoutdirs_off, mtzoutdirs_on, outname):
                '''
                Create the list of list to run Xtrapol8 in multiprocessing for all JK files
                Args:
                    mtzoutdirs_off: list of [directory where to find mtz file, mtz file complete directory, True/False:the file is for total data] for off state
                    mtzoutdirs_on: list of [directory where to find mtz file, mtz file complete directory, True/False:the file is for total data] for on state

                Returns:
                    mtzoutdirs_dir_off_on: list of [outdirJKX8, mtz_off, mtz_on, newlog, nb_JK, outname]
                '''
                mtzoutdirs_dir_off_on=[]

                for i in range(0,len(mtzoutdirs_off)): #for all the datasets
                    i_mtzoutdirs_off = mtzoutdirs_off[i] #get the list with [directory where to find mtz file, mtz file complete directory, True/False:the file is for total data] for each dataset of the off state
                    i_mtzoutdirs_on = mtzoutdirs_on[i] #get the list with [directory where to find mtz file, mtz file complete directory, True/False:the file is for total data] for each dataset of the on state

                    outdirJKX8 = i_mtzoutdirs_on[0] + '/Xtrapol8' #take output directory from the mtz_on file to put results of Xtrapol8
                    os.mkdir(outdirJKX8)
                    mtz_off = i_mtzoutdirs_off[1]
                    JK_i_off = i_mtzoutdirs_off[2] #nb_JK
                    mtz_on = i_mtzoutdirs_on[1]
                    JK_i_on = i_mtzoutdirs_on[2]

                    # Check if the two total mtz files are correctly paired else get the two numbers of JK
                    if JK_i_off == JK_i_on :  # if the mtz files are correctly paired
                        # get log file: create new log files
                        newlogname, newlog = JK_utils.create_log(outdirJKX8)
                        print('The log file of Xtrapol8 launched with reference mtz = %s and triggered mtz = %s is: %s' % (
                            mtz_off, mtz_on, newlogname), file=log)
                        # add files to list
                        mtzoutdirs_dir_off_on.append([outdirJKX8, mtz_off, mtz_on, newlogname, JK_i_off, outname])

                    else:
                        # get log file: create new log files
                        newlogname, newlog = JK_utils.create_log(outdirJKX8)
                        print('The log file of Xtrapol8 launched with reference mtz = %s and triggered mtz = %s is: %s' % (
                            mtz_off, mtz_on, newlogname), file=log)
                        # add files to list
                        mtzoutdirs_dir_off_on.append([outdirJKX8, mtz_off, mtz_on, newlogname, str(JK_i_off)+','+str(JK_i_on), outname])

                return(mtzoutdirs_dir_off_on)

            mtzoutdirs_dir_off_on = create_multiprocessing_list(mtzoutdirs_off, mtzoutdirs_on, outname)

        # create the list mtzoutdirs_dir_off_on_total with [outdirJKX8_total, mtz_off_total, mtz_on_total, newlogname, 'total', outname]
            if P.fraction != 1: #if fraction==1, no total created
                # create the list mtzoutdirs_dir_off_on_total with [outdirJKX8, mtz_off, mtz_on, newlog, 'total', outname]
                def create_total_multiprocessing_list(mtzoutdirs_off_total, mtzoutdirs_on_total, outname):
                    '''
                    Create the list to run Xtrapol8 with the total files
                    Args:
                        mtzoutdirs_off_total: [directory where to find total mtz file, total mtz file complete directory, 'total'] for total off state
                        mtzoutdirs_on_total: [directory where to find total mtz file, total mtz file complete directory, 'total'] for total on state

                    Returns:
                        mtzoutdirs_dir_off_on_total: [outdirJKX8, mtz_off, mtz_on, newlog, 'total', outname]
                    '''
                    outdirJKX8_total = mtzoutdirs_on_total[0] + '/Xtrapol8'  # take output directory from the mtz_on file to put results of Xtrapol8
                    os.mkdir(outdirJKX8_total)
                    # get log file: create new log files
                    newlogname, newlog = JK_utils.create_log(outdirJKX8_total)

                    mtz_off_total = mtzoutdirs_off_total[1]
                    mtz_on_total = mtzoutdirs_on_total[1]
                    print('The log file of Xtrapol8 launched with reference mtz = %s and triggered mtz = %s is: %s' % (
                                mtz_off_total, mtz_on_total, newlogname), file=log)
                    # add files to list
                    mtzoutdirs_dir_off_on_total = [outdirJKX8_total, mtz_off_total, mtz_on_total, newlogname, 'total', outname]
                    return (mtzoutdirs_dir_off_on_total)

                mtzoutdirs_dir_off_on_total = create_total_multiprocessing_list(mtzoutdirs_off_total, mtzoutdirs_on_total, outname)

                print(mtzoutdirs_dir_off_on_total)
            print(mtzoutdirs_dir_off_on)


    #Run Xtrapol8
            if P.run_Xtrapol8:
                print_in_T_and_log('\n################################################################################\nLAUNCH XTRAPOL8\n################################################################################')

                # #Loop to run Xtrapol8 one by one for debugging
                # tab_list=[]
                # for outdir_and_mtz_file_off_on in mtzoutdirs_dir_off_on:
                #     # init
                #     # get output directory for Xtrapol8
                #     outdir1 = outdir_and_mtz_file_off_on[0]
                #     os.chdir(outdir1)
                #     # get mtz files
                #     mtz_off = outdir_and_mtz_file_off_on[1]
                #     mtz_on = outdir_and_mtz_file_off_on[2]
                #     newlogname=outdir_and_mtz_file_off_on[3]
                #     total = outdir_and_mtz_file_off_on[4]
                #
                #     # JK_utils.run_in_terminal('sftools << eof\n'
                #     #                  'READ %s\n'
                #     #                  'REDUCE CCP4\n'
                #     #                  'WRITE %s_sftools.mtz\n'
                #     #                  'STOP\n'
                #     #                  'eof'% (mtz_off,mtz_off),
                #     #                  existing_files=[mtz_off +'_sftools.mtz'])
                #     # mtz_off=mtz_off +'_sftools.mtz'
                #     # JK_utils.run_in_terminal('sftools << eof\n'
                #     #                  'READ %s\n'
                #     #                  'REDUCE CCP4\n'
                #     #                  'WRITE %s_sftools.mtz\n'
                #     #                  'STOP\n'
                #     #                  'eof'% (mtz_on,mtz_on),
                #     #                  existing_files=[mtz_on +'_sftools.mtz'])
                #     # mtz_on=mtz_on +'_sftools.mtz'
                #     # check mtz files
                #     check_all_files([mtz_off, mtz_on])
                #
                #     tab_listi=main_X8.run_X8([outdir1, mtz_off, mtz_on, newlogname, total], params, P, master_phil, P.startdir)
                #     tab_list.append(tab_listi)

                pool_process = multiprocessing.Pool(P.processors)  # Create a multiprocessing Pool with the number of processors defined
                tab_list = pool_process.map(partial(main_X8.run_X8, params=params, P=P, master_phil=master_phil, startdir=P.startdir), mtzoutdirs_dir_off_on)  # process the Xtrapol8 with mtz files iterable with pool in mtzoutdirs_dir_off_on


                if P.fraction != 1:
                    tab_total = main_X8.run_X8(mtzoutdirs_dir_off_on_total, params, P, master_phil, P.startdir)  # run Xtrapol8 for the total files

                    print(tab_total)
                print(tab_list)

    #Run simple refinements
            else: #only JK is run so a simple refinement is done
                print_in_T_and_log('\n################################################################################\nLAUNCH SIMPLE REFINEMENTS\n################################################################################')
                for mtzoutdirs_dir_off_on_outname in mtzoutdirs_dir_off_on:
                    tab_list = JK_simple_refinement.run_simple_refinement(params, P, mtzoutdirs_dir_off_on_outname) #run simple refinement for the JK files

                if P.fraction != 1:
                    tab_total = JK_simple_refinement.run_simple_refinement(params, P, mtzoutdirs_dir_off_on_outname)  # run simple refinement for the total files


    #Get JK results
            if P.fraction != 1:
                #compare the files from JK and total and get the results (CC, RMSD, Plot of differences between models)
                JK_results_outdir = outdir + '/JK_average_and_comparison_results' #name of new output directory
                print(JK_results_outdir)
                os.mkdir(JK_results_outdir) #create new output directory
                JK_utils.print_terminal_and_log('CALCULATING JK RESULTS\n==============================================================')
                JK_results.get_JK_results(tab_total, tab_list, JK_results_outdir, params.Xtrapol8.refinement.phenix_keywords.main.ordered_solvent) #calculate all the results and print them in a table and create plot
    ################################################################

    # run only Xtrapol8
    elif P.run_Xtrapol8 == True:
        index=0 #index of the triggered mtz
        list_for_X8 = [] #create the list for multiprocessing
        for triggered_mtz in params.Xtrapol8.input.triggered_mtz:#in case of multiple triggered mtz given
            P.get_parameters_multiple_JK_X8_output(params, index) #get the outdir and outname
            if len(P.triggered_mtz) > 1: logdir, log = JK_utils.create_log(P.output, global_log=True) #if there are many triggered mtz files create a new log file and get full directory of the log file
            list_for_X8.append([P.output, params.Xtrapol8.input.reference_mtz, triggered_mtz, logdir, 'total', P.outname])
            index += 1 #get the next index


            print_in_T_and_log(
                '################################################################################\nLAUNCH XTRAPOL8 ONLY\n################################################################################')
            pool_process = multiprocessing.Pool(P.processors)  # Create a multiprocessing Pool with the number of processors defined
            tab_list = pool_process.map(partial(main_X8.run_X8, params=params, P=P, master_phil=master_phil, startdir=P.startdir), list_for_X8)  # process the Xtrapol8 with mtz files iterable with pool in list_for_X8


    ################################################################

    # Write all input parameters to a phil file.
    modified_phil.show(out=open("Xtrapol8_in.phil", "w"))
    if params.Xtrapol8.output.generate_phil_only:
        log.close()
        sys.exit()
########################################################################################################################
if __name__ == "__main__":
    run(sys.argv[1:])
    print('End')
