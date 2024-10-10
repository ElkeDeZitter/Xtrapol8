# -*- coding: utf-8 -*-
"""
Main script to run Xtrapol8

-------

usage
-------
Fextr.py is the main script to be called with phenix.python

Run without any argument to see all options and explanation:
>>> phenix.python <wherever>/Fextr.py
Parameters can be added using an input file or via command line

example using input file (preferable)
-------
1) Change the nano Xtrapol8.phil using your favourite editor, e.g.
>>> nano Xtrapol8.phil
2) Run Xtrapol8
>>> phenix.python <wherever>/Fextr.py Xtrapol8.phil

example using command line only
-------
1) Run Xtrapol8 with all your arguments
>>> phenix.python <wherever>/Fextr.py input.reference_mtz=hiephiep.mtz input.triggered_mtz=hieperdepiep.mtz input.reference_pdb=hoera.pdb input.additional_files=jeej.cif input.additional_files=another.cif occupancies.list_occ=0.1,0.3,0.5 f_and_maps.f_extrapolated_and_maps=qfextr,qfgenick map_explorer.peak_integration_floor=3.5 map_explorer.peak_detection_threshold=4 output.outdir=fancy_party

example using input file and command line
-------
1) Change the nano Xtrapol8.phil using your favourite editor, e.g.
>>> nano trapol8.phil
2) Xtrapol8 with additional arguments. The order of arguments determines how paramters will be overwriten:
>>> phenix.python <wherever>/Fextr.py trapol8.phil refinement.phenix_keywords.refine.cycles=3

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

TODO:
- clean up Fextr_utils: remove unnecessary functions
- Resolve problem with 90% multiplicity estimation
- option to automatically create cif with phenix.elbow
- GUI: warning message if people want to run crazy stuff that they better run on a powerfull machine in command line mode
- update example phil and manual
- Suggest resolution cutoff based on Riso and CCiso, write resolution suggestions based on true completeness and <F/SIG(F)> to logfile
- add comments to various steps of script (especially the weird steps)
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
- plot the F and sig(F) after negative handling
- ddm with linear scale / sum of differences instead of plotting for each atom
"""
from __future__ import division, print_function
import re
import os
import sys
import random
import subprocess
import shutil
import pickle
from select import select
from datetime import datetime
import numpy as np
import uuid

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
from iotbx import ccp4_map
from scipy.stats import pearsonr
from wx.lib.pubsub import pub

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
import version

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
scattering_table = *n_gaussian wk1995 it1992 electron neutron
    .type = choice(multi=False)
    .help = Scattering table. For testing only. It is unclear if Xtrapol8 will give proper results with data other than originating from X-ray diffraction and a scattering table different from n_gaussian. If using anything else than n_gaussian, use only phenix for refinement.
    .expert_level = 3
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
    negative_and_missing = *truncate_and_fill truncate_no_fill fref_and_fill fref_no_fill fcalc_and_fill fcalc_no_fill keep_and_fill keep_no_fill reject_and_fill reject_no_fill zero_and_fill zero_no_fill fill_missing no_fill
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
    use_refmac_instead_of_phenix = False
        .type = bool
        .help = use Refmac for reciprocal space refinement and COOT for real-space refinement instead of phenix.refine and phenix.real_space_refine.
        .expert_level = 0
    phenix_keywords{
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
        real_space_refine{
            cycles = 5
            .type = int
            .help = Number of refinement cycles for real space refinement.
            .expert_level = 0
            }
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

class SymManager(symmetry.manager):
    """
    Use symmetry manager to compare unit cell and space group with that of the model.
    """

    def __init__(self, model_in):
        super(SymManager, self).__init__()

        self.model = model_in
        self.process_pdb_file(self.model)

    def check_symm(self,*hkls):

        for reflections in hkls:
            #the first check can only be applied for the unit cell as space groups belonging to the same point group can be allowed, while we need to have exactly the same space groups
            _, UC_err =self.process_reflections_file(reflections)  # Not sure which tolerance range is applied
            SG_ok = (str(self.current_space_group) == str(reflections.crystal_symmetry().space_group_info()))
            err = 0
            err_m = 'Space group and unit cell compatible'

            if UC_err:
                err += 1
                err_m = '!!!Unit cell incompatibility!!! Do you want to continue? This might just be nonsense.'
            if SG_ok == False:
                err += 1
                err_m = '!!!Space group incompatibility: %s vs %s!!! Do you want to continue? This might just be nonsense.'                %(str(self.current_space_group), str(reflections.crystal_symmetry().space_group_info()))
        return err_m, err
    
        

class DataHandler(object):
    """
    Handle all input file and generate objects to be used in map calculations and analyses.
    """

    def __init__(self, pdb_in, mtz_off, additional, outdir, mtz_on):

        self.pdb_in            = pdb_in
        self.mtz_off           = mtz_off
        self.mtz_on            = mtz_on
        self.additional        = additional
        self.outdir            = outdir
        
    def check_outdir(self):
        """
        Make output directory:
        - if no name specified, then it will be called Xtrapol8
        - if the output directory already exists, then a number will be added
         This way creates a maximum of 1000 Xtrapol8 output directories
        """
        if self.outdir == None:
            #self.outdir = os.getcwd()
            self.outdir = "Xtrapol8"
            
        #else:
            #if os.path.exists(self.outdir) == False:
                #try:
                    #os.mkdir(self.outdir)
                    #print('Output directory not present thus being created: %s'%(self.outdir))
                #except OSError:
                    #os.makedirs(self.outdir)
            #self.outdir = os.path.abspath(self.outdir)
            
        outdir = self.outdir
        i = 1
        while os.path.exists(outdir):
            if os.path.isdir(outdir):
                if len(os.listdir(outdir)) ==0:
                    #outdir = self.outdir
                    break
            ##Keep outdir given by user if it only contains Xtrapol8 log-files:
            #if len([fle for fle in os.listdir(self.outdir) if fle.endswith("Xtrapol8.log")]) == len(os.listdir(self.outdir)):
                #outdir = self.outdir
                #break
            outdir = "%s_%d" %(self.outdir, i)
            i += 1
            if i == 1000: #to avoid endless loop, but this leads to a max of 1000 Xtrapol8 runs
                break
            
        try:
            os.mkdir(outdir)
            print('Output directory being created: %s'%(outdir))
        except OSError:
            try:
                os.makedirs(outdir)
                print('Output directory being created: %s'%(outdir))
            except OSError:
                print("Output directory: %s" %(outdir))
            
        self.outdir = os.path.abspath(outdir)
                

    def open_files(self):
        """
        check and read input files and output directory
        """
        self.check_outdir()
        self.reflections_off = any_file(self.mtz_off, force_type="hkl", raise_sorry_if_errors=True)
        self.reflections_on = any_file(self.mtz_on, force_type="hkl", raise_sorry_if_errors=True)
        self.from_cif_create_pdb_file()
        self.model_in = any_file(self.check_and_delete_hydrogen(), force_type="pdb", raise_sorry_if_errors=True)
        # No H deletion for electron diffraction, however leads to crash of fmodel update all scales
        # self.model_in = any_file(self.pdb_in, force_type="pdb", raise_sorry_if_errors=True)

    def check_all_files(self):

        err = 0
        err_m = ''
        for fle in [self.pdb_in, self.mtz_off,self.mtz_on]:
            if not self.check_single_file(fle):
                err = 1
                err_m += '\nFile not found: %s'%fle
        if err == 0: self.open_files()
        return err, err_m

    def check_single_file(self, fle):
        if fle == None:
            return False
        else:
            return os.path.isfile(fle)

    def from_cif_create_pdb_file(self):
        """
        If the input model file (pdb_in) has a format cif, a new pdb file of the model is created and replaces the cif file (pdb_in)

        Parameters:
        -----------
        self.pdb_in (file)
            model file input by the user
        self.outdir (string)
            output directory
        file_format (format, example: .cif)
            format of the model file input by the user

        Returns:
        -----------
        self.pdb_in (pdb file)
            Created model pdb file from the model cif input file
        """
        
        try:  # iotbx.file_reader.splitext from phenix 1.19
            _, file_format, _ = iotbx.file_reader.splitext(self.pdb_in)
        except ValueError:  # iotbx.file_reader.splitext from phenix 1.18
            _, file_format = os.path.splitext(self.pdb_in)
            
        # get the file path and type
        if file_format == '.cif':
        # if the file type is cif
            pdb_hier = hierarchy.input(file_name=self.pdb_in)
            p = open(self.outdir+'/'+get_name(self.pdb_in)+'.pdb', 'w')
            # create and open a new pdb file to write in, the file will be in the new outdir and have the same name as the input file with a pdb format
            p.write(pdb_hier.hierarchy.as_pdb_string(crystal_symmetry=pdb_hier.input.crystal_symmetry()))
            # write the columns of the cif file into the pdb file
            p.close()
            if self.check_single_file(self.outdir+'/'+get_name(self.pdb_in)+'.pdb'):
            # check if the created pdb file exists
                self.pdb_in = os.path.abspath(self.outdir+'/'+get_name(self.pdb_in)+'.pdb')
                # the input model used is the created pdb file with the info from cif file

    def check_and_delete_hydrogen(self):
        """
        If hydrogen exists, make new pdb-file with removed hydrogens.
        Possibility to achieve extrapolated structure factors with such high precision as to see hydrogens is very low.
        """
        pdb_hier = hierarchy.input(file_name=self.pdb_in)
        outname = 'model_edit.pdb'
        if pdb_hier.hierarchy.remove_hd() != 0:
            print("Model contains hydrogen atoms. Create a new model without these atoms in the output directory: %s" %(outname))
            print("Model contains hydrogen atoms. Create a new model without these atoms in the output directory: %s" % (outname), file=log)
            pdb_hier.hierarchy.remove_hd()
            p = open('%s/%s' %(self.outdir, outname), 'w')
            p.write(pdb_hier.hierarchy.as_pdb_string(crystal_symmetry=pdb_hier.input.crystal_symmetry()))
            p.close()
            if self.check_single_file('%s/%s' %(self.outdir, outname)):
                self.pdb_in = os.path.abspath('%s/%s' %(self.outdir, outname))
            else:
                self.pdb_in = os.path.abspath(self.pdb_in)
        else:
            self.pdb_in = os.path.abspath(self.pdb_in) #We will need absolute path of pdb file for later refinements in subdirectories
        return self.pdb_in

    def extract_fobs(self, low_res, high_res):
        """
        Extract the actual reflections from the data files and cut at resolution limits (if set)
        For now Friedel pairs will have to be merged.
        """
        self.fobs_off, self.fobs_on = Column_extraction(self.reflections_off,
                                                        self.reflections_on,
                                                        low_res,
                                                        high_res,
                                                        log = log).extract_columns()
        
        if self.fobs_off.anomalous_flag():
            print("I promised to keep the anomalous flags, but that was a lie. Xtrapol8 is not yet ready to handle anomalous data. For now, your Friedel pairs will be merged.", file=log)
            print("I promised to keep the anomalous flags, but that was a lie. Xtrapol8 is not yet ready to handle anomalous data. For now, your Friedel pairs will be merged.")
            self.fobs_off = self.fobs_off.average_bijvoet_mates()
            self.fobs_on  = self.fobs_on.average_bijvoet_mates()
        
        self.fobs_off = self.fobs_off.map_to_asu()
        #self.fobs_off = self.resolution_cutoff(self.fobs_off, low_res, high_res)
        self.fobs_on  = self.fobs_on.map_to_asu()
        #self.fobs_on  = self.resolution_cutoff(self.fobs_on, low_res, high_res)

        #self.fobs_off = self.extract_colums(self.reflections_off, low_res, high_res)
        #self.fobs_on  = self.extract_colums(self.reflections_on, low_res, high_res)
        #self.fobs_on = []
        #for on in self.reflections_on:
            #f_on = self.extract_colums(on, res)
            #self.fobs_on.append(f_on)
            
    def get_UC_and_SG(self):
        """
        Extract unit cell and space group from the model.
        """
        #self.SG = re.search(r"(.+?)\(No",self.fobs_off.space_group_info().symbol_and_number()).group(1)
        #self.UC = self.fobs_off.unit_cell()
        self.SG = str(self.model_in.crystal_symmetry().space_group_info())
        self.UC = self.model_in.crystal_symmetry().unit_cell()
                    
    def resolution_cutoff(self, f_obs, low_res, high_res):
        """
        Cut data at low and high resolution, only if the data extend beyond the limit.
        """
        dmax, dmin = f_obs.d_max_min()
        if (high_res != None and high_res > dmin):
            dmin = high_res
        if (low_res != None and low_res < dmax):
            dmax = low_res
        return f_obs.resolution_filter(dmax, dmin)
        
    def check_additional_files(self):
        """
        Check for all additional files if they exist, if the minimum of additional files is present for phenix.refine and convert to string of additional files to add as an argument to phenix.refine later.
        """
        self.cif_objects = []
        additional = " "
        if self.additional!= None:
            for fle in self.additional:
                if 'cif' in fle:
                    try:
                        cif_object=mmtbx.monomer_library.server.read_cif(file_name=fle)
                    except Exception:
                        raise  AssertionError,"Unable to read the cif file "+fle
                    else:
                        if (len(cif_object) > 0):
                            self.cif_objects.append((fle,cif_object))
                if len(fle) > 0:
                    if self.check_single_file(fle):
                        #new_add = os.path.abspath(fle)
                        new_add = os.path.realpath(fle)
                        additional = additional + "%s " %(new_add)
        SpaceGroup=sgtbx.space_group_info(symbol=self.SG)
        crystal_symmetry=crystal.symmetry(unit_cell=self.UC,space_group_info=SpaceGroup)
        processed_pdb_files_srv = utils.process_pdb_file_srv(
            crystal_symmetry          = crystal_symmetry,
            pdb_parameters            = pdb.input(self.pdb_in),
            cif_objects               = self.cif_objects)
            #log                       = log)
        processed_pdb_files_srv.process_pdb_files(pdb_file_names = [self.pdb_in])
        #processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(pdb_file_names = [self.pdb_in])
        #processed_pdb_file and pdb_inp can also be used to get x_ray_structure etc
        self.additional = additional
        
    def extract_ligand_codes(self):
        """
        Get list with three-letter codes for ligands. Needed in further analysis.
        """
        cif_list = []
        for cif in self.cif_objects:
            normal_cif = False
            for comp in cif[1]:
                try: #This should work for a proper cif file
                    ligand = re.search(r'comp_(.+?)$', comp).group(1)
                    if len(ligand) == 3 and ligand not in cif_list:
                        cif_list.append(ligand)
                    normal_cif == True
                except AttributeError:
                    continue
            if normal_cif: #when multiple ligands are merged into one cif-file
                for ligand in cif[1]['comp_list']['_chem_comp.id']:
                    if ligand not in cif_list:
                        cif_list.append(ligand)
            else:
                try: #This is for a raw downloaded cif file (bugs may appear later in Xtrapol8). Requires a single ligand per file for now
                    ligand = cif[1].keys()[0]
                    if len(ligand) == 3 and ligand not in cif_list:
                        cif_list.append(ligand)
                except AttributeError:
                    continue
        return cif_list

    def generate_Rfree(self, array, fraction):
        """
        Generate random Rfree reflections.
        Not so clever because only created at random without taking symmetry into account.
        """
        #tot = self.fobs_off.data().size()
        #num_free = int(tot * fraction)
        #free_ind = random.sample(range(tot), num_free)
        #free_col = np.zeros(tot, dtype=np.int32)
        #for i in range(tot):
            #if i in free_ind:
                #free_col[i] = 1
        #self.rfree = miller.array(miller_set=self.fobs_off, data=flex.int(free_col))
        
        self.rfree = array.generate_r_free_flags(fraction=fraction)

    def generate_f_model(self, fobs, scattering_table = "n_gaussian"):
        """
        Make an fmodel from data and input model
        """
        #Not sure if neg reflections are already removed from start. Next might always be 0 neg regflections instead of real number
        #print("Generate fmodel from off-state dataset. %d negative reflections will be converted to positive reflections. This is %.2f %% of the total amount of reflections and might impact calculation of Fo-Fo and extrapolated structure factors and maps." %(fobs.select(~(fobs.data() >=0)).data().size(), fobs.select(~(fobs.data() >=0)).data().size()/fobs.data().size() *100))
        #print("Generate fmodel from off-state dataset. %d negative reflections. This is %.2f %% of the total amount of reflections and might impact calculation of Fo-Fo and extrapolated structure factors and maps." %(fobs.select(~(fobs.data() >=0)).data().size(), fobs.select(~(fobs.data() >=0)).data().size()/fobs.data().size() *100), file=log)
        
        #r_free_flags = miller.array(miller_set=self.rfree, data=self.rfree.data().as_bool())
        
        pdb_ini = iotbx.pdb.input(self.pdb_in)
        xray_structure = pdb_ini.xray_structure_simple()
        
        xray_structure.scattering_type_registry(table=scattering_table)
        
        #to change of basis to reference settings is a bad idea because it will mess up things
        #phil_xs = crystal.symmetry(unit_cell=fobs.unit_cell(),
                                    #space_group_info=fobs.space_group_info())
        #to_reference = phil_xs.change_of_basis_op_to_reference_setting()
        #xray_structure = xray_structure.change_basis(to_reference)

        self.fmodel = mmtbx.f_model.manager(
            f_obs          = fobs,
            r_free_flags   = self.rfree,
            xray_structure = xray_structure)
                
    def update_fmodel(self, fobs):
        """
        to change f_obs or remove reflections.
        if the number of Rfree-flagged reflections is too low (consequennce of outlier rejection during scaling), re-assign
        """
        if self.fmodel.f_obs().data().size() > fobs.data().size():
            #print('updating fmodel, this can take a while...')
            #selection = flex.bool()
            #for ind in DH.fmodel.f_obs().indices():
                #if ind in fobs.indices():
                    #selection.append(True)
                #else:
                    #selection.append(False)
            #self.fmodel.select(selection)
            dt = np.dtype([('h', np.int) ,('k',np.int),('l',np.int)])
            f1 = np.array(self.fmodel.f_obs().indices(), dtype=dt)
            f2 = np.array(fobs.indices(), dtype=dt)
            selection = np.nonzero(np.in1d(f1, f2))[0]
            assert selection.shape[0] == fobs.indices().size()
            
            sel = np.zeros(f1.shape[0])
            np.put(sel, selection, 1)
            sel_bool = map(lambda x: bool(int(x)), list(sel))
            selection_flags = flex.bool(sel_bool)
            self.fmodel = self.fmodel.select(selection_flags)
            
            self.fmodel.update(f_obs=fobs)
            
            
            tot  = self.fmodel.r_free_flags().size()
            free = self.fmodel.r_free_flags().select(self.fmodel.r_free_flags().data()==True).size()
            frac = free/tot
            if frac < 0.048:
                print("Rfree fraction too low, re-assign Rfree flags")
                self.rfree = self.fmodel.f_obs().generate_r_free_flags(fraction=0.05) #during the two scaling steps some structure factors might be removed and thus new Rfree reflections have to be chosen
                self.fmodel.update(r_free_flags = self.rfree)
        
            
    def get_common_indices_and_Fobs_off(self, f_obs_on):
        """
        Ugly function to compare all reflections and only keep those that are common between the datasets
        """
        
        f_model_scaled   = self.fmodel.f_model()
        f_obs_off        = self.fobs_off
        f_obs_off_scaled = self.fmodel.f_obs()
        rfree            = self.rfree
        
        off_ini = f_obs_off.data().size()
        on_ini  = f_obs_on.data().size()
        
        f_model_scaled, f_obs_on = f_model_scaled.common_sets(f_obs_on)
        f_model_scaled, f_obs_off = f_model_scaled.common_sets(f_obs_off)
        f_model_scaled, f_obs_off_scaled = f_model_scaled.common_sets(f_obs_off_scaled)
        f_model_scaled, rfree = f_model_scaled.common_sets(rfree)
        if check_common_indices(
                [f_model_scaled, f_obs_on, f_obs_off, f_obs_off_scaled, rfree]) == False:
            f_model_scaled, f_obs_on = f_model_scaled.common_sets(f_obs_on)
            f_model_scaled, f_obs_off = f_model_scaled.common_sets(f_obs_off)
            f_model_scaled, f_obs_off_scaled = f_model_scaled.common_sets(f_obs_off_scaled)
            f_model_scaled, rfree = f_model_scaled.common_sets(rfree)
            if check_common_indices(
                    [f_model_scaled, f_obs_on, f_obs_off, f_obs_off_scaled, rfree]) == False:
                print (
                "I tried to maintain only those indices which all three data sets have in common. Nevertheless, the two data sets have not the same indices. The program will probably stop with and error and/or output will be nonsense.")
                
        off_fin = f_obs_off.data().size()
        on_fin  = f_obs_on.data().size()
        print("%d reflections rejected from off-state dataset, %d reflections left" %(off_ini-off_fin, off_fin))
        print("%d reflections rejected from on-state dataset, %d reflections left" %(on_ini-on_fin, on_fin))
        print("%d reflections rejected from off-state dataset, %d reflections left"
              %(off_ini-off_fin, off_fin), file=log)
        print("%d reflections rejected from on-state dataset, %d reflections left"
              %(on_ini-on_fin, on_fin), file=log)
                
        self.f_model_scaled  = f_model_scaled
        self.fobs_off        = f_obs_off
        self.fobs_off_scaled = f_obs_off_scaled
        self.rfree           = rfree
        self.indices         = self.fobs_off.indices()
                
        self.scale_sigmas_from_fmodel()        
        
        return f_obs_on

    def check_common_indices(self, refl_lst):
        for a in refl_lst:
            for b in refl_lst:
                if b != a:
                    if a.indices().all_eq(b.indices()) == False or a.data().size() != b.data().size():
                        compatible = False
                        break
                    else:
                        compatible = True
        return compatible
    
    def scale_sigmas_from_fmodel(self):
        """
        Scale sigmas. As performed in phenix.FoFo
        """
        sc = flex.sum(self.fobs_off_scaled.data()*self.f_model_scaled.amplitudes().data())/flex.sum(self.f_model_scaled.amplitudes().data()*self.f_model_scaled.amplitudes().data())
        #self.fobs_off_scaled = make_miller_array(self.fobs_off_scaled.data(),self.fobs_off.sigmas()/sc, self.SG, self.UC, self.indices)
        self.fobs_off_scaled = miller.array(miller_set=self.fobs_off,
                                            data=self.fobs_off_scaled.data(),
                                            sigmas=self.fobs_off.sigmas()/sc)
        
    def scale_fobss(self, b_scaling, low_res=None, high_res=None):
        """
        Scale triggered mtz with internally scaled reference mtz using scaleit.
        """
        if b_scaling == 'no':
            print("Fobs,reference and Fobs,triggered not being scaled", file=log)
            print("Fobs,reference and Fobs,triggered not being scaled")
            self.fobs_on_scaled = self.fobs_on #don't scale
            #define the scaling resolution boundaries so that variables exist.
            dmax_off, dmin_off = self.fobs_off_scaled.d_max_min()
            dmax_on, dmin_on = self.fobs_on.d_max_min()
            self.scaling_dmin = np.max([dmin_off, dmin_on])
            self.scaling_dmax = np.min([dmax_off, dmax_on])
        else:
            dmax_off, dmin_off = self.fobs_off_scaled.d_max_min()
            dmax_on, dmin_on = self.fobs_on.d_max_min()
            #get the high resolution edge for scaling (no data truncation)
            if high_res != None:
                self.scaling_dmin = np.max([high_res, dmin_off, dmin_on])
            else:
                self.scaling_dmin = np.max([dmin_off, dmin_on])
            #get the low resolution edge for scaling (no data truncation)
            if low_res != None:
                self.scaling_dmax = np.min([low_res, dmax_off, dmax_on])
            else:
                self.scaling_dmax = np.min([dmax_off, dmax_on])
                
            print("Fobs,reference and Fobs,triggered scaled using scaleit", file=log)
            print("Fobs,reference and Fobs,triggered scaled using scaleit")
            #self.fobs_on_scaled = scalef_cnslike(self.fobs_off_scaled, self.fobs_on, self.SG, self.rfree, bscale=b_scaling) #run CNS-like scaling
            self.fobs_on_scaled = run_scaleit(self.fobs_off_scaled, self.fobs_on, b_scaling, low_res=self.scaling_dmax, high_res=self.scaling_dmin) #prepare mtz-file and run scaleit

class FobsFobs(object):
    """
    Class for the calculation of weighting difference structure factors and maps.
    """
    def __init__(self, fobs_on, fobs_off):
        self.fobs_on    = fobs_on
        self.fobs_off   = fobs_off
        self.indices    = fobs_off.indices()
        self.get_UC_and_SG()
        
    def get_UC_and_SG(self):
        """
        Extract unit cell and space group from the reference data set
        """
        self.SG = re.search(r"(.+?)\(No",self.fobs_off.space_group_info().symbol_and_number()).group(1)
        self.UC = self.fobs_off.unit_cell()
        
    def get_fdif(self):
        """
        Unweighted difference structure factors
        """
        self.fdif = (self.fobs_on.data()-self.fobs_off.data())
        
    def get_sigf(self):
        """
        Unweighted difference_sigmas
        """
        self.sigf_diff = (self.fobs_on.sigmas()**2+self.fobs_off.sigmas()**2)**(0.5)
                
    def q_weighting(self):
        """
        Calculate q-weight
        """
        #q = calculate_q(self.fobs_off, self.fobs_on)
        #q_ms   = make_miller_array(self.fobs_off.data(), q, self.SG, self.UC, self.indices)
        q_ms, self.q_av = calculate_q(self.fobs_off, self.fobs_on, log=log)
        self.q = q_ms.sigmas()
        
    def k_weighting(self, kweight_scale):
        """
        Calculate k-weight, still under development
        """
        k_ms, self.k_av = calculate_k(self.fobs_off, self.fobs_on, kweight_scale = kweight_scale, log=log)
        self.k = k_ms.sigmas()
        
    def outlier_rejection(self):
        """
        Outlier rejection if no weighting is performed.
        """
        c_ms = outlier_rejection_only(self.fobs_off, self.fobs_on, log=log)
        self.c = c_ms.sigmas()
    
    def calculate_fdiff(self, kweight_scale = 0.05):
        """
        Actual calculation of the differences
        """
        self.get_fdif()
        self.get_sigf()
        
        print("----Calculating weight factors----", file=log)
        #calculate with qweight
        self.q_weighting()
        weight      = self.q/self.q_av
        self.fdif_q = weight*self.fdif
        sigf_diff_q = weight*self.sigf_diff
        self.fdif_q_ms= make_miller_array(self.fdif_q, sigf_diff_q, self.SG, self.UC, self.indices)
        
        #calculate with kweight
        self.k_weighting(kweight_scale)
        weight      = self.k/self.k_av
        self.fdif_k = weight*self.fdif
        sigf_diff_k = weight*self.sigf_diff
        self.fdif_k_ms= make_miller_array(self.fdif_k, sigf_diff_k, self.SG, self.UC, self.indices)
        
        #calculate without weight
        self.outlier_rejection()
        self.fdif_c = self.c*self.fdif
        self.fdif_c_ms= make_miller_array(self.fdif_c,self.sigf_diff, self.SG, self.UC, self.indices)
        
    def write_maps(self, fmodel, rfree, outname, qweighting = True, kweighting=False):
        """
        Write FoFo maps
        """
        print("\n************Write reflection files and calculate maps************")
        print("\n************Write reflection files and calculate maps************", file=log)
        if qweighting:
            maptype = 'qFoFo'
            F = Filesandmaps(self.fdif_q_ms, rfree, maptype, outname, fmodel)
            #self.mtz_name, self.ccp4_name, self.xplor_name = Filesandmaps(self.fdif_q_ms, rfree, maptype, outname, fmodel).write_FoFo_output()
        elif kweighting:
            maptype = 'kFoFo'
            F = Filesandmaps(self.fdif_k_ms, rfree, maptype, outname, fmodel)
            #self.mtz_name, self.ccp4_name, self.xplor_name = Filesandmaps(self.fdif_k_ms, rfree, maptype, outname, fmodel).write_FoFo_output()
        else:
            maptype = 'FoFo'
            F = Filesandmaps(self.fdif_c_ms, rfree, maptype, outname, fmodel)
            #self.mtz_name, self.ccp4_name, self.xplor_name = Filesandmaps(self.fdif_c_ms, rfree, maptype, outname, fmodel).write_FoFo_output()
        #self.mtz_name, self.ccp4_name, self.xplor_name = F.write_FoFo_output()
        self.mtz_name, self.ccp4_name = F.write_FoFo_output()
        print("{:s} maps:".format(maptype))
        print("  mtz-format: {:s}".format(self.mtz_name))
        print("  ccp4-format: {:s}".format(self.ccp4_name))
        #print("  xplor-format: {:s}".format(self.xplor_name))
        print("{:s} maps:".format(maptype), file=log)
        print("  mtz-format: {:s}".format(self.mtz_name), file=log)
        print("  ccp4-format: {:s}".format(self.ccp4_name), file=log)
        #print("  xplor-format: {:s}".format(self.xplor_name), file=log)


        self.crystal_gridding = F.crystal_gridding

        
    def get_crystal_gridding(self):
        """
        Return the crystal_gridding of the FoFo map so that the same gridding can be used to calculate other maps on the same grid.
        """
        if self.crystal_gridding:
            return self.crystal_gridding
        else:
            print("Crystal gridding not defined yet. Have you calculated maps before?")
            return 1
       
       
class Fextrapolate(object):
    """
    Class for the calculation, analysis and usage of extrapolated structure factors.
    """
    def __init__(self,
                 fdif,
                 fdif_q,
                 fdif_k,
                 sigf_diff,
                 q_values,
                 k_values,
                 fobs_off,
                 fobs_on,
                 fmodel_fobs_off,
                 rfree,
                 occ              = 1,
                 name_out         = 'Fextrapolate',
                 neg_refl_handle  = 'fill_missing',
                 crystal_gridding = None):
        self.fdif            = fdif
        self.fdif_q          = fdif_q
        self.fdif_k          = fdif_k
        self.sigf_diff       = sigf_diff
        self.q_values        = q_values
        self.k_values        = k_values
        self.fobs_on         = fobs_on
        self.fobs_off        = fobs_off
        self.fmodel_fobs_off = fmodel_fobs_off
        self.rfree           = rfree
        self.occ             = occ
        self.alf             = 1/(self.occ)
        self.neg_refl_handle = neg_refl_handle
        self.indices         = fobs_off.indices()
        self.crystal_gridding= crystal_gridding
        self.get_UC_and_SG()
        #keep_no_fill replaces the old "no_fill" since it becomes more clear what we do with negatives (we dont' do anything, we keep them negative)
        #keep_and_fill replaces the old "fill_missing" since it becomes more clear what we do with negatives (we dont' do anything, we keep them negative)
        #the argument "no_fill" and "fill_missing" remain valid to maintain backwards compatibility
        if neg_refl_handle in ['no_fill', 'keep_no_fill', 'reject_no_fill', 'zero_no_fill', 'fcalc_no_fill', 'fref_no_fill', 'truncate_no_fill']: #'addconstant_no_fill', 'massage_no_fill'
            self.fill_missing = False
        else:
            self.fill_missing = True

        self.name_out = "%s_occ%.3f" %(name_out, occ)
    
    def get_UC_and_SG(self):
        """
        Extract unit cell and space group from the reference data set
        """
        self.SG = re.search(r"(.+?)\(No",self.fobs_off.space_group_info().symbol_and_number()).group(1)
        self.UC = self.fobs_off.unit_cell()
    
    def create_output_dirs(self, outdir):
        """
        Create directories - if not existing yet - for the occupancy and whether or not q-weighting is applied.
        """
        new_dir_q = "qweight_occupancy_%.3f" %(self.occ)
        new_dir_k = "kweight_occupancy_%.3f" %(self.occ)
        new_dir   = "occupancy_%.3f" %(self.occ)
        new_dirpath_q = "%s/%s" %(outdir, new_dir_q)
        new_dirpath_k = "%s/%s" %(outdir, new_dir_k)
        new_dirpath   = "%s/%s" %(outdir, new_dir)
        if os.path.exists(new_dirpath_q) == False:
            os.mkdir(new_dirpath_q)
        if os.path.exists(new_dirpath_k) == False:
            os.mkdir(new_dirpath_k)
        if os.path.exists(new_dirpath) == False:
            os.mkdir(new_dirpath)
        return new_dirpath_q, new_dirpath_k, new_dirpath
    
    def fextr_q(self):
        """
        Calculate q-weighted extrapolated structure factors of type fextr or fgenick
        """
        fextr_q = self.alf*self.fdif_q+self.fobs_off.data()
        self.fextr_ms = make_miller_array(fextr_q, self.fextr_sigmas_q(), self.SG, self.UC, self.indices)

    def fextr_k(self):
        """
        Calculate k-weighted extrapolated structure factors of type fextr or fgenick
        """
        fextr_k = self.alf*self.fdif_k+self.fobs_off.data()
        self.fextr_ms = make_miller_array(fextr_k, self.fextr_sigmas_k(), self.SG, self.UC, self.indices)

    def fextr_noq(self):
        """
        Calculate non-q/k-weighted extrapolated structure factors of type fextr or fgenick
        """
        fextr = self.alf*self.fdif+self.fobs_off.data()
        self.fextr_ms = make_miller_array(fextr, self.fextr_sigmas_noq(), self.SG, self.UC, self.indices)
        
    def fextr_calc_q(self):
        """
        Calculate q-weighted extrapolated structure factors of type fextr_calc
        """
        fextr_q = self.alf*self.fdif_q+self.fmodel_fobs_off.f_model().amplitudes().data()
        self.fextr_calc_ms = make_miller_array(fextr_q, self.fextr_sigmas_q(), self.SG, self.UC, self.indices)

    def fextr_calc_k(self):
        """
        Calculate k-weighted extrapolated structure factors of type fextr_calc
        """
        fextr_k = self.alf*self.fdif_k+self.fmodel_fobs_off.f_model().amplitudes().data()
        self.fextr_calc_ms = make_miller_array(fextr_k, self.fextr_sigmas_k(), self.SG, self.UC, self.indices)

    def fextr_calc_noq(self):
        """
        Calculate non-q/k-weighted extrapolated structure factors of type fextr_calc
        """
        fextr = self.alf*self.fdif+self.fmodel_fobs_off.f_model().amplitudes().data()
        self.fextr_calc_ms = make_miller_array(fextr, self.fextr_sigmas_noq(), self.SG, self.UC, self.indices)
                
    def fextr_sigmas_q(self):
        """
        Calculate q-weighted extrapolated sigmas
        """
        q_factor = self.q_values/flex.mean(self.q_values)
        return ((self.fobs_on.sigmas()*self.alf*q_factor)**2+(self.fobs_off.sigmas()*(1-(q_factor*self.alf)))**2)**(0.5)
    
    def fextr_sigmas_k(self):
        """
        Calculate k-weighted extrapolated sigmas
        """
        k_factor = self.k_values/flex.mean(self.k_values)
        return ((self.fobs_on.sigmas()*self.alf*k_factor)**2+(self.fobs_off.sigmas()*(1-(k_factor*self.alf)))**2)**(0.5)

    def fextr_sigmas_noq(self):
        """
        Calculate non-q-weighted extrapolated sigmas
        """
        return ((self.fobs_on.sigmas()*self.alf)**2+(self.fobs_off.sigmas()*(1-self.alf))**2)**(0.5)
    
    def message(self):
        print("\n---CALCULATING %s TYPE OF ESFAS AND MAPS FOR OCCUPANCY %.3f (ALPHA = %.3f)---" %(self.maptype.upper(), self.occ, self.alf), file=log)
        print("\n---CALCULATING %s TYPE OF ESFAS AND MAPS FOR OCCUPANCY %.3f (ALPHA = %.3f)---" %(self.maptype.upper(), self.occ, self.alf))
        
    #Next come several functions for the handling of negative reflections
    def negatives_reject(self, ms):
        """
        Remove negative reflections. Should be called in case of "reject_and_fill" or "reject_no_fill"
        """
        print("Negative ESFAs will be removed")
        print("Negative ESFAs will be removed",file=log)
        return ms.select(ms.data()>=0)
        
    def negatives_zero(self, ms):
        """
        Set negative reflections to zero. Should be called in case of "zero_and_fill" or "zero_no_fill"
        """
        print("Negative ESFAs will be set to 0")
        print("Negative ESFAs will be set to 0",file=log)
        data = ms.data().deep_copy()
        data.set_selected(~(ms.data()>=0),0)
        sigmas = ms.sigmas().deep_copy()
        sigmas.set_selected(~(ms.data()>=0),0)
        return miller.array(miller_set = ms,
                            data       = data,
                            sigmas     = sigmas)

    #def add_minimum_to_all(self, ms):
        """
        Add the minimum value to all refletcions. This makes no sense
        """
        #print("Add constant to all relfections to avoid negative reflections")
        #print("Add constant to all relfections to avoid negative reflections",file=log)        
        #return miller.array(miller_set = ms,
                            #data       = ms.data()-flex.min(ms.data()),
                            #sigmas     = ms.sigmas())

    def negatives_fcalc(self, ms):
        """
        Replace the negative reflections by fcalc. Should be called in case of "fcalc_and_fill" or "fcalc_no_fill"
        """
        print("Negative ESFAs replaced by Fcalc")
        print("Negative ESFAs replaced by Fcalc",file=log)        
        data = ms.data().deep_copy()
        sel = ms.data()<0
        data.set_selected(sel, self.fmodel_fobs_off.f_model().amplitudes().data())
        return miller.array(miller_set = ms, 
                            data       = data,
                            sigmas     = ms.sigmas())
    
    def negatives_foff(self, ms):
        """
        Replace the negative reflections by Fobs_reference. Should be called in case of "fref_and_fill" or "fref_no_fill"
        """
        print("Negative ESFAs replaced by Foff")
        print("Negative ESFAs replaced by Foff",file=log)        
        data = ms.data().deep_copy()
        sel = ms.data()<0
        data.set_selected(sel, self.fmodel_fobs_off.f_obs().data())
        return miller.array(miller_set = ms, 
                            data       = data,
                            sigmas     = ms.sigmas())
    
    def convert_to_I_then_to_F(self, ms, maptype, algorithm='reflection_file_converter'):
        """
        Use truncate to estimate the negative reflections. This strategy to get rid of negative reflections takes the following steps:
        1) square the Fss to get Is, but keep sign of F so that that negative Fs end up at negative Is
        2) run truncate to convert from I back to F and thus estimate Fs based on French-Wilson distrubtion and relations
        This means that also the values of the possitive reflecions are altered in order to fulfill the Wilson distrubtions
        """
        print("Square ESFAs to estimate Is, but keep sign of Fs")
        print("Square ESFAs to estimate Is, but keep sign of Fs", file=log)
        I_data = ms.data()**2
        I_data.set_selected(ms.data()<0, I_data*(-1))
        I_sigmas = ms.sigmas()*ms.data()
        I_sigmas.set_selected(ms.data()<0, I_sigmas*(-1))
        Iextr = miller.array(miller_set = ms, 
                             data       = I_data,
                             sigmas     = I_sigmas)
        Iextr.set_observation_type_xray_intensity()
        
        labels = re.sub("f","i", self.maptype.lower()).upper()
        suffix = re.sub("f","i", self.maptype.lower())
        if algorithm == 'truncate':
            outname = '%s_fortruncate.mtz' %(suffix)
        else:
            outname = '%s_forconverter.mtz' %(suffix)
        mtz_dataset = Iextr.as_mtz_dataset(column_root_label=labels)
        mtz_dataset.mtz_object().write(outname)
        
        if algorithm == 'truncate':
            Corrected_Fs = Extrapolated_column_extraction(outname, labels, log).get_Fs_from_truncate()
        else:
            Corrected_Fs = Extrapolated_column_extraction(outname, labels, log).get_Fs_from_reflection_file_converter()
        
        if (Corrected_Fs==None or Corrected_Fs.size()==0):
            if "no_fill" in self.neg_refl_handle:
                new_neg_refl_handle = "reject_no_fill"
                #new_neg_refl_handle = "keep_no_fill"
            else:
                new_neg_refl_handle = "reject_and_fill"
                #new_neg_refl_handle = "keep_and_fill"
            print("   Cannot successfully run truncate/phenix.reflection_file_converter. The reason might be the high number of negative ESFAs and their very high absolute value. The negative ESFAs will be removed and we try again. This means that %s will be used for this dataset instead of %s. This might impact further analysis and comparison of the electron density maps." %(new_neg_refl_handle, self.neg_refl_handle))
            print("   Cannot successfully run truncate/phenix.reflection_file_converter. The reason might be the high number of negative ESFAs and their very high absolute value. The negative ESFAs will be removed and we try again. This means that %s will be used for this dataset instead of %s. This might impact further analysis and comparison of the electron density maps." %(new_neg_refl_handle, self.neg_refl_handle), file=log)
            Corrected_Fs = self.negatives_reject(ms)
        else:
            Corrected_Fs = Corrected_Fs.map_to_asu()
        return Corrected_Fs
    
    def initial_maps(self, ms):
        """
        Calculate the initial maps in which the negative ESFAs are maintained (only true if this is run BEFORE negatives handling).
        Maps are calculated in mtz and ccp4 format.
        If the map cannot be calculated because the number of negatives is too high, then the map will not be calculated. This is
        different from the behavior for the negatives in which the negatives will be rejected if maps cannot be calculated with the chosen strategy.
        Missing reflections will not be filled.
        Maps are stored in seperate directory
        """
        
        indir  = os.getcwd()
        outdir = "maps-keep_no_fill"
        if os.path.isdir(outdir) == False:
            os.mkdir(outdir)
        os.chdir(outdir)
            
        print("\n************Calculate maps************")
        print("\n************Calculate maps************", file=log)
        name_out = "{:s}_keep_no_fill".format(self.name_out)
        FM = Filesandmaps(ms, self.rfree, self.maptype, name_out, self.fmodel_fobs_off, crystal_gridding=self.crystal_gridding)
        print("Electron density map in which the negative calculated ESFAs are preserved, and without filling of missing reflections:")
        print("Electron density map in which the negative calculated ESFAs are preserved, and without filling of missing reflections:", file=log)
        try:
            mtz_name, ccp4_name_2FoFc, ccp4_name_FoFc = FM.write_Fextr_maps(fill_missing=False)
            print("  mtz-format: {:s}".format(mtz_name))
            print("  mtz-format: {:s}".format(mtz_name), file=log)
            print("  ccp4-format: {:s} and {:s}".format(ccp4_name_2FoFc, ccp4_name_FoFc))
            print("  ccp4-format: {:s} and {:s}".format(ccp4_name_2FoFc, ccp4_name_FoFc), file=log)
        except AssertionError:
            print("  Cannot update and calculate scales for electron density maps. The reason might be the high number of negative ESFAs. The maps cannot be calculated.")
            print("  Cannot update and calculate scales for electron density maps. The reason might be the high number of negative ESFAs. The maps cannot be calculated.", file=log)
            
        os.chdir(indir)
            
    def negatives_and_maps_Fextr_Fextr_calc(self, ms, outdir_for_negstats = os.getcwd()):
        """
        Calculate ESFSs and extrapolated maps in which the negatives are treated as indicated by the negative_and_missing strategy.
        Whether missing reflections are automatically filled depends on the negative_and_missing strategy.
        If the map cannot be calculated because the number of negatives is too high, then negatives will be rejected and a new attempt for map calculation will be carried out. A warning will be given.
        ESFAs are written in mtz format.
        Maps are calculated in mtz, ccp4 format. Xplor format are removed from Xtrapol8 and all lines are commented
        """
        neg_neflecions_binning(ms, self.maptype, log=log)
        neg_reflections = ms.select(~(ms.data()>=0))
        dump_negative_stats(self.occ, self.maptype, ms.data().size(), neg_reflections.data().size(), outdir_for_negstats)
        print("Total number of negative ESFAs: %d (%.2f %% of the data)"
                %(neg_reflections.data().size(), neg_reflections.data().size()/ms.data().size() *100))
        print("Total number of negative ESFAs: %d (%.2f %% of the data)"
                %(neg_reflections.data().size(), neg_reflections.data().size()/ms.data().size() *100),file=log)
        
        print("\n************Write reflection files and calculate maps************")
        print("\n************Write reflection files and calculate maps************", file=log)
        if neg_reflections.data().size() > 0:
            print("Negative ESFAs handling:", file=log)
            print("Negative ESFAs handling:")
            if self.neg_refl_handle in ['reject_no_fill', 'reject_and_fill']:
                fextr_ms = self.negatives_reject(ms)
                rfree, fmodel_fobs_off = self.get_updated_fmodel_fobs_off(fextr_ms)
                self.FM = Filesandmaps(fextr_ms, rfree, self.maptype, self.name_out, fmodel_fobs_off, crystal_gridding=self.crystal_gridding)
            elif self.neg_refl_handle in ['zero_and_fill', 'zero_no_fill']:
                fextr_ms = self.negatives_zero(ms)
                self.FM = Filesandmaps(fextr_ms, self.rfree, self.maptype, self.name_out, self.fmodel_fobs_off, crystal_gridding=self.crystal_gridding)
            #elif self.neg_refl_handle in ['addconstant_and_fill', 'addconstant_no_fill']:
                #fextr_ms = self.add_minimum_to_all(ms)
                ##No need to update fmodel because for fextr and fextr_calc maps we only use the xray-model from it
                #self.FM = Filesandmaps(fextr_ms, self.rfree, self.maptype, self.name_out, self.fmodel_fobs_off, crystal_gridding=self.crystal_gridding)
            elif self.neg_refl_handle in ['fcalc_and_fill', 'fcalc_no_fill']:
                fextr_ms = self.negatives_fcalc(ms)
                self.FM = Filesandmaps(fextr_ms, self.rfree, self.maptype, self.name_out, self.fmodel_fobs_off, crystal_gridding=self.crystal_gridding)
            elif self.neg_refl_handle in ['fref_and_fill', 'fref_no_fill']:
                fextr_ms = self.negatives_foff(ms)
                self.FM = Filesandmaps(fextr_ms, self.rfree, self.maptype, self.name_out, self.fmodel_fobs_off, crystal_gridding=self.crystal_gridding)
            elif self.neg_refl_handle in ['truncate_and_fill', 'truncate_no_fill']:
                fextr_ms = self.convert_to_I_then_to_F(ms, self.maptype, algorithm='truncate')
                rfree, fmodel_fobs_off = self.get_updated_fmodel_fobs_off(fextr_ms)
                self.FM = Filesandmaps(fextr_ms, rfree, self.maptype, self.name_out, fmodel_fobs_off, crystal_gridding=self.crystal_gridding)
            # elif self.neg_refl_handle in ['massage_and_fill', 'massage_no_fill']:
            #     fextr_ms = self.convert_to_I_then_to_F(ms, self.maptype, algorithm='reflection_file_converter')
            #     rfree, fmodel_fobs_off = self.get_updated_fmodel_fobs_off(fextr_ms)
            #     self.FM = Filesandmaps(fextr_ms, rfree, self.maptype, self.name_out, fmodel_fobs_off)
            #The following else should be not required if self.initial_maps() has been run before
            else: #in case of keep_no_fill (old "no_fill") and keep_and_fill (old "fill_missing") 
                self.FM = Filesandmaps(ms, self.rfree, self.maptype, self.name_out, self.fmodel_fobs_off, crystal_gridding=self.crystal_gridding)
        #he following else should be not required if self.initial_maps() has been run before
        else: #in case of zero negatives 
            self.FM = Filesandmaps(ms, self.rfree, self.maptype, self.name_out, self.fmodel_fobs_off, crystal_gridding=self.crystal_gridding)
            
        #try to calculate the map with the chosen negative handling
        try:
            fm = self.FM.write_Fextr_Fextr_calc_output(self.fill_missing)
        except AssertionError:
            if "no_fill" in self.neg_refl_handle:
                new_neg_refl_handle = "reject_no_fill"
            else:
                new_neg_refl_handle = "reject_and_fill"
            print("  Cannot update and calculate scales for electron density maps. The reason might be the high number of negative ESFAs. The negative ESFAs will be removed and we try again. This means that %s will be used for this dataset instead of %s. This might impact further analysis and comparison of the electron density maps." %(new_neg_refl_handle, self.neg_refl_handle))
            print("  Cannot update and calculate scales for electron density maps. The reason might be the high number of negative ESFAs. The negative ESFAs will be removed and we try again. This means that %s will be used for this dataset instead of %s. This might impact further analysis and comparison of the electron density maps." %(new_neg_refl_handle, self.neg_refl_handle), file=log)
            #fextr_ms = self.negatives_reject(self.fextr_ms)
            fextr_ms = self.negatives_reject(ms)
            rfree, fmodel_fobs_off = self.get_updated_fmodel_fobs_off(fextr_ms)
            self.FM = Filesandmaps(fextr_ms, rfree, self.maptype, self.name_out, fmodel_fobs_off, crystal_gridding=self.crystal_gridding)
            fm = self.FM.write_Fextr_Fextr_calc_output(self.fill_missing)
                  
        #self.F_name, self.mtz_name, self.ccp4_name_2FoFc, self.ccp4_name_FoFc, self.xplor_name_2FoFc, self.xplor_name_FoFc = fm
        self.F_name, self.mtz_name, self.ccp4_name_2FoFc, self.ccp4_name_FoFc = fm
        print("\nESFAs and extraplated electron density maps:")
        print("\nESFAs and extraplated electron density maps:", file=log)
        print("  ESFAs in mtz-format: {:s}".format(self.F_name))
        print("  ESFAs in mtz-format: {:s}".format(self.F_name), file=log)
        print("  electron density in mtz-format: {:s}".format(self.mtz_name))
        print("  electron density in mtz-format: {:s}".format(self.mtz_name), file=log)
        print("  electron density in ccp4-format: {:s} and {:s}".format(self.ccp4_name_2FoFc, self.ccp4_name_FoFc))
        print("  electron density in ccp4-format: {:s} and {:s}".format(self.ccp4_name_2FoFc, self.ccp4_name_FoFc), file=log)
        #print("  electron density in xplor-format: {:s} and {:s}".format(self.xplor_name_2FoFc, self.xplor_name_FoFc))
        #print("  electron density in xplor-format: {:s} and {:s}".format(self.xplor_name_2FoFc, self.xplor_name_FoFc), file=log)
        
    #Next come several functions to calculate the three types of extrapolated structure factors 
    def fextr(self, qweight=False, kweight=False, outdir_for_negstats = os.getcwd()):
        """
        Calculation of (q-weighted) Fextr. Should be called in case extrapolated Fs of type "qFextr" or "Fextr"
        Different steps:
        1) calculate extrapolated structure factors
        2) handle the negative reflections
        3) write mtz-files and maps
        """
        if qweight:
            self.fextr_q()
            self.maptype = 'qFextr'
        elif kweight:
            self.fextr_k()
            self.maptype = 'kFextr'            
        else:
            self.fextr_noq()
            self.maptype = 'Fextr'
        self.message()
        
        #write the map with all reflections independent of the negative handling choses
        self.initial_maps(self.fextr_ms)
        
        #Carry out handling of negative ESFAs and write files and maps.
        self.negatives_and_maps_Fextr_Fextr_calc(self.fextr_ms, outdir_for_negstats = outdir_for_negstats)
        

    def fgenick(self, qweight=False, kweight=False, outdir_for_negstats = os.getcwd()):
        """"
        Calculation of (q-weighted) Fgenick. Should be called in case extrapolated Fs of type "qFgenick" or "Fgenick" (Genick, 2007)
        As this method follows the description by Genick, 2007, negtive reflections are automatically rejected
        Different steps:
        1) calculate extrapolated structure factors
        2) reject negative reflections.
        3) write mtz-files and maps
        """
        if qweight:
            self.fextr_q()
            self.fgenick_ms = self.fextr_ms.deep_copy()
            self.maptype = 'qFgenick'
        elif kweight:
            self.fextr_k()
            self.fgenick_ms = self.fextr_ms.deep_copy()
            self.maptype = 'kFgenick'
        else:
            self.fextr_noq()
            self.fgenick_ms = self.fextr_ms.deep_copy()
            self.maptype = 'Fgenick'
        self.message()
        
        #write the map with all reflections independent of the negative handling choses
        self.initial_maps(self.fgenick_ms)

        #Carry out handling of negative ESFAs and write files and maps.
        #In case of Fgenick the negatives are always rejected
        neg_neflecions_binning(self.fgenick_ms, self.maptype, log=log)
        neg_reflections = self.fgenick_ms.select(~(self.fgenick_ms.data()>=0))
        dump_negative_stats(self.occ, self.maptype, self.fgenick_ms.data().size(), neg_reflections.data().size(), outdir_for_negstats)
        print("Total number of negative ESFAs: %d (%.2f %% of the data)" %(neg_reflections.data().size(), neg_reflections.data().size()/self.fgenick_ms.data().size() *100))
        print("Total number of negative ESFAs: %d (%.2f %% of the data)"
                %(neg_reflections.data().size(), neg_reflections.data().size()/self.fgenick_ms.data().size() *100), file=log)
        
        print("\n************Write reflection files and calculate maps************")
        print("\n************Write reflection files and calculate maps************", file=log)
        if neg_reflections.data().size() > 0:
            print("Negative ESFAs handling:", file=log)
            print("Negative ESFAs handling:")
            self.fgenick_ms = self.negatives_reject(self.fgenick_ms)
            rfree, fmodel_fobs_off = self.get_updated_fmodel_fobs_off(self.fgenick_ms)
            self.FM = Filesandmaps(self.fgenick_ms, rfree, self.maptype, self.name_out, fmodel_fobs_off, crystal_gridding=self.crystal_gridding)
        else:
            self.FM = Filesandmaps(self.fgenick_ms, self.rfree, self.maptype, self.name_out, self.fmodel_fobs_off, crystal_gridding=self.crystal_gridding)
                
        fm = self.FM.write_Fgenick_output()
        #self.F_name, self.mtz_name, self.ccp4_name_2FoFc, self.ccp4_name_FoFc, self.xplor_name_2FoFc, self.xplor_name_FoFc = fm
        self.F_name, self.mtz_name, self.ccp4_name_2FoFc, self.ccp4_name_FoFc = fm

        print("\nESFAs and extraplated electron density maps:")
        print("\nESFAs and extraplated electron density maps:", file=log)
        print("  ESFAs in mtz-format: {:s}".format(self.F_name))
        print("  ESFAs in mtz-format: {:s}".format(self.F_name), file=log)
        print("  electron density in mtz-format: {:s}".format(self.mtz_name))
        print("  electron density in mtz-format: {:s}".format(self.mtz_name), file=log)
        print("  electron density in ccp4-format: {:s} and {:s}".format(self.ccp4_name_2FoFc, self.ccp4_name_FoFc))
        print("  electron density in ccp4-format: {:s} and {:s}".format(self.ccp4_name_2FoFc, self.ccp4_name_FoFc), file=log)
        #print("  electron density in xplor-format: {:s} and {:s}".format(self.xplor_name_2FoFc, self.xplor_name_FoFc))
        #print("  electron density in xplor-format: {:s} and {:s}".format(self.xplor_name_2FoFc, self.xplor_name_FoFc), file=log)

        
    def fextr_calc(self, qweight=False, kweight=False, outdir_for_negstats = os.getcwd()):
        """"
        Calculation of (q-weighted) Fextr_calc. Should be called in case extrapolated Fs of type "qFextr_calc" or "Fextr_calc"
        Different steps:
        1) calculate extrapolated structure factors
        2) handle the negative reflections
        3) write mtz-files and maps
        """
        if qweight:
            self.fextr_calc_q()
            self.maptype = 'qFextr_calc'
        elif kweight:
            self.fextr_calc_k()
            self.maptype = 'kFextr_calc'
        else:
            self.fextr_calc_noq()
            self.maptype = 'Fextr_calc'
        self.message()
            
        #write the map with all reflections independent of the negative handling choses
        self.initial_maps(self.fextr_calc_ms)
        
        #Carry out handling of negative ESFAs and write files and maps.
        self.negatives_and_maps_Fextr_Fextr_calc(self.fextr_calc_ms, outdir_for_negstats = outdir_for_negstats)
                
    def get_updated_fmodel_fobs_off(self, miller_array):
        """
        Update the Fmodel that is associated with the off-state reflections and model in order to have the same reflections as the extrapolated structure factors.
        Generally, this is necessary when the negative reflections are rejected and hence the indices of the fmodel and extrapolated strcucture factors do not correspond anymore.
        """
        fmodel_f_obs, miller_array = self.fmodel_fobs_off.f_obs().common_sets(miller_array)
        if fmodel_f_obs.data().size() < self.fmodel_fobs_off.f_obs().data().size():
            rfree, miller_array = self.rfree.common_sets(miller_array)
            
            dt = np.dtype([('h', np.int) ,('k',np.int),('l',np.int)])
            f1 = np.array(self.fmodel_fobs_off.f_obs().indices(), dtype=dt)
            f2 = np.array(fmodel_f_obs.indices(), dtype=dt)
            selection = np.nonzero(np.in1d(f1, f2))[0]
            assert selection.shape[0] == fmodel_f_obs.indices().size()
            
            sel = np.zeros(f1.shape[0])
            np.put(sel, selection, 1)
            sel_bool = map(lambda x: bool(int(x)), list(sel))
            selection_flags = flex.bool(sel_bool)
            fmodel_fobs_off = self.fmodel_fobs_off.select(selection_flags)
            
            fmodel_fobs_off.update(f_obs=fmodel_f_obs)
            return rfree, fmodel_fobs_off
        else:
            return self.rfree, self.fmodel_fobs_off
        
        
    def phenix_phenix_refinements(self,
                           mtz_F=None,
                           mtz_map=None,
                           pdb_in=None,
                           additional=None,
                           F_column_labels='QFEXTR',
                           column_labels='2FOFCWT,PH2FOFCWT',
                           scattering_table= "n_gaussian",
                           keywords = {}):
        
        """
        Reciprocal space and real space refinement in the extrapolated structure factors and map coefficients, respectively, using Phenix.
        1) reciprocal space refinement in extrapolated structure factors (mtz_F and pdb_in)
        2) real-space refinement in Fextr map coefficients (mtz_map and pdb_in)
        3) real space refinement with results reciprocal space refinement (mtz_out and pd_out = output of step 1)
        !!! take care: different definition of column labels as compared to refmac_coot_refinements!!! here: values and phases for 2FoFc kind of map only
        """
        
        if mtz_F == None:
            mtz_F = self.F_name
        if mtz_map == None:
            mtz_map = self.mtz_name
        assert pdb_in != None, 'Specify pdb for refinement'
        
        ref = phenix_refinements.Phenix_refinements(mtz_F,
                                                    pdb_in, 
                                                    additional         = additional,
                                                    F_column_labels    = F_column_labels,
                                                    strategy           = keywords.refine.strategy,
                                                    rec_cycles         = keywords.main.cycles,
                                                    real_cycles        = keywords.real_space_refine.cycles,
                                                    wxc_scale          = keywords.target_weights.wxc_scale,
                                                    wxu_scale          = keywords.target_weights.wxu_scale,
                                                    solvent            = keywords.main.ordered_solvent,
                                                    sim_annealing      = keywords.main.simulated_annealing,
                                                    sim_annealing_pars = keywords.simulated_annealing,
                                                    map_sharpening     = keywords.map_sharpening.map_sharpening,
                                                    scattering_table   = scattering_table,
                                                    weight_sel_crit    = keywords.target_weights.weight_selection_criteria,
                                                    additional_reciprocal_keywords = keywords.additional_reciprocal_space_keywords,
                                                    additional_real_keywords       = keywords.additional_real_space_keywords,
                                                    log                = log)
        
        #print("Refinements:", file=log)
        #print("Refinements:")
        
        print("RECIPROCAL SPACE REFINEMENT WITH %s AND %s" %(mtz_F, pdb_in))
        mtz_out_rec, pdb_out_rec = ref.phenix_reciprocal_space_refinement()
        print("Output reciprocal space refinement:", file=log)
        print("----------------")
        print("Output reciprocal space refinement:")
        if os.path.isfile(pdb_out_rec):
            print("    pdb-file: %s"%(pdb_out_rec), file=log)
            print("    pdb-file: %s"%(pdb_out_rec))
        else:
            print("    pdb-file not found, %s incorrectly returned" %(pdb_in), file=log)
            print("    pdb-file not found, %s incorrectly returned" %(pdb_in))
            pdb_out_rec = pdb_in
        if os.path.isfile(mtz_out_rec):
            print("    mtz-file: %s"%(mtz_out_rec), file=log)
            print("    mtz-file: %s"%(mtz_out_rec))
        else:
            print("    mtz-file not found. Refinement failed.", file=log)
            print("    mtz-file not found. Refinement failed.")
        print("----------------")
        
        if keywords.density_modification.density_modification:
            print("DENSITY MODIFICATION WITH %s AND %s" %(mtz_F, pdb_out_rec))
            # mtz_dm = ref.phenix_density_modification(mtz_out_rec, pdb_out_rec)
            mtz_dm = ref.ccp4_dm(pdb_out_rec, keywords.density_modification.combine, keywords.density_modification.cycles)
            print("Output density modification:", file=log)
            print("Output density modification:")
            if os.path.isfile(mtz_dm):
                print("    mtz-file: %s"%(mtz_dm), file=log)
                print("    mtz-file: %s" % (mtz_dm))
            else:
                print("    mtz-file not found. Density modification failed.", file=log)
                print("    mtz-file not found. Density modification failed.")

        print("REAL SPACE REFINEMENT WITH %s AND %s" %(mtz_map, pdb_in))
        pdb_out_real = ref.phenix_real_space_refinement_mtz(mtz_map, pdb_in, column_labels)
        print("Output real space refinement:", file=log)
        print("----------------")
        print("Output real space refinement:")
        if os.path.isfile(pdb_out_real):
            print("    pdb-file: %s"%(pdb_out_real), file=log)
            print("    pdb-file: %s"%(pdb_out_real))
        else:
            print("    pdb-file not found, %s incorrectly returned" %(pdb_in), file=log)
            print("    pdb-file not found, %s incorrectly returned" %(pdb_in))
            pdb_out_real = pdb_in
        print("----------------")

        if (keywords.density_modification.density_modification and os.path.isfile(mtz_dm)):
            ccp4_dm = re.sub(r".mtz$", ".ccp4", mtz_dm)
            _, high_res = ref.get_mtz_resolution(mtz_dm)
            print("REAL SPACE REFINEMENT WITH %s AND %s" %(ccp4_dm, pdb_out_rec))
            pdb_out_rec_real = ref.phenix_real_space_refinement_ccp4(ccp4_dm, pdb_out_rec, high_res)
            print("Output real space refinement after reciprocal space refinement:", file=log)
            print("----------------")
            print("Output real space refinement after reciprocal space refinement:")
            if os.path.isfile(pdb_out_rec_real):
                print("    pdb-file: %s"%(pdb_out_rec_real), file=log)
                print("    pdb-file: %s"%(pdb_out_rec_real))
            else:
                print("    pdb-file not found, %s incorrectly returned" %(pdb_out_rec), file=log)
                print("    pdb-file not found, %s incorrectly returned" %(pdb_out_rec))
                pdb_out_rec_real = pdb_out_rec
            print("----------------")
        else:
            print("REAL SPACE REFINEMENT WITH %s AND %s" %(mtz_out_rec, pdb_out_rec))
            pdb_out_rec_real = ref.phenix_real_space_refinement_mtz(mtz_out_rec, pdb_out_rec, '2FOFCWT,PH2FOFCWT')
            print("Output real space refinement after reciprocal space refinement:", file=log)
            print("----------------")
            print("Output real space refinement after reciprocal space refinement:")
            if os.path.isfile(pdb_out_rec_real):
                print("    pdb-file: %s"%(pdb_out_rec_real), file=log)
                print("    pdb-file: %s"%(pdb_out_rec_real))
            else:
                print("    pdb-file not found, %s incorrectly returned" %(pdb_out_rec), file=log)
                print("    pdb-file not found, %s incorrectly returned" %(pdb_out_rec))
                pdb_out_rec_real = pdb_out_rec
            print("----------------")
               
        return mtz_out_rec, pdb_out_rec, pdb_out_real, pdb_out_rec_real
    
    def refmac_coot_refinements(self,
                                mtz_F=None,
                                mtz_map=None,
                                pdb_in=None,
                                additional=None,
                                ligands_list=None,
                                F_column_labels='QFEXTR',
                                rfree_col = 'FreeR_flag',
                                map_column_labels='2FOFCWT, PH2FOFCWT, FOFCWT, PHFOFCWT',
                                keywords = {}):
        """
        Reciprocal space and real space refinement in the extrapolated structure factors using Refmac and map coefficients using COOT, respectively.
        1) reciprocal space refinement in extrapolated structure factors (mtz_F and pdb_in)
        2) real-space refinement in Fextr map coefficients (mtz_map and pdb_in)
        3) real space refinement with results reciprocal space refinement (mtz_out and pd_out = output of step 1)
        !!! take care: different definition of column labels as compared to phenix_refinements!!! here: values and phases for 2FoFc and FoFc kind of map
        """
        if mtz_F == None:
            mtz_F = self.F_name
        if mtz_map == None:
            mtz_map = self.mtz_name
        assert pdb_in != None, 'Specify pdb for refinement'
           
        ref = ccp4_refmac.refmac_refinements(mtz_F,
                 pdb_in,
                 additional,
                 ligands_list,
                 F_column_labels = F_column_labels,
                 rfree_col       = rfree_col,
                 fill_missing    = False,
                 add_rfree       = False,
                 refinement_weight         = keywords.target_weights.weight,
                 refinement_weight_sigmas  = keywords.target_weights.experimental_sigmas,
                 refinement_weighting_term = keywords.target_weights.weighting_term,
                 refinement_type       = keywords.refine.type,
                 TLS                   = keywords.refine.TLS,
                 TLS_cycles            = keywords.refine.TLS_cycles,
                 bfac_set              = keywords.refine.bfac_set,
                 twinning              = keywords.refine.twinning,
                 Brefinement           = keywords.refine.Brefinement,
                 cycles                = keywords.refine.cycles,
                 external_restraints   = keywords.restraints.external_restraints,
                 jelly_body_refinement = keywords.restraints.jelly_body_refinement,
                 jelly_body_sigma      = keywords.restraints.jelly_body_sigma,
                 jelly_body_additional_restraints = keywords.restraints.jelly_body_additional_restraints,
                 map_sharpening        = keywords.map_sharpening.map_sharpening,
                 density_modification  = keywords.density_modification.density_modification,
                 dm_combine            = keywords.density_modification.combine,
                 dm_ncycle             = keywords.density_modification.cycles,
                 additional_reciprocal_keywords = keywords.additional_refmac_keywords)

        print("Refinements:", file=log)
        #print("Refinements:")
            
        print("RECIPROCAL SPACE REFINEMENT WITH %s AND %s" %(mtz_F, pdb_in))
        if keywords.density_modification.density_modification:
            mtz_out_rec, pdb_out_rec, mtz_dm = ref.refmac_reciprocal_space_refinement()
        else:
            mtz_out_rec, pdb_out_rec,_ = ref.refmac_reciprocal_space_refinement()
        print("output reciprocal space refinement:", file=log)
        print("----------------")
        print("output reciprocal space refinement:")
        if os.path.isfile(pdb_out_rec):
            print("    pdb-file: %s"%(pdb_out_rec), file=log)
            print("    pdb-file: %s"%(pdb_out_rec))
        else:
            print("    pdb-file not found, %s incorrectly returned" %(pdb_in), file=log)
            print("    pdb-file not found, %s incorrectly returned" %(pdb_in))
            pdb_out_rec = pdb_in
        if os.path.isfile(mtz_out_rec):
            print("    mtz-file: %s"%(mtz_out_rec), file=log)
            print("    mtz-file: %s"%(mtz_out_rec))
        elif mtz_out_rec == mtz_F:
            print("    mtz-file not found. Refinement failed.", file=log)
            print("    mtz-file not found. Refinement failed.")
        else:
            print("    mtz-file not found. Refinement failed.", file=log)
            print("    mtz-file not found. Refinement failed.")
        if keywords.density_modification.density_modification:
            print("Output density modification:", file=log)
            print("Output density modification:")
            if mtz_dm == mtz_out_rec:
                print("    mtz-file not found. Density modification failed.", file=log)
                print("    mtz-file not found. Density modification failed.")                
            elif os.path.isfile(mtz_dm):
                print("    mtz-file: %s" % (mtz_dm), file=log)
                print("    mtz-file: %s" % (mtz_dm))
            else:
                print("    mtz-file not found. Density modification failed.", file=log)
                print("    mtz-file not found. Density modification failed.")
        print("----------------")

        print("REAL SPACE REFINEMENT WITH %s AND %s" %(mtz_map, pdb_in))
        pdb_out_real = ref.coot_real_space_refinement_mtz(pdb_in, mtz_map, map_column_labels)
        print("output real space refinement:", file=log)
        print("----------------")
        print("output real space refinement:")
        if os.path.isfile(pdb_out_real):
            print("    pdb-file: %s"%(pdb_out_real), file=log)
            print("    pdb-file: %s"%(pdb_out_real))
        else:
            print("    pdb-file not found, %s incorrectly returned" %(pdb_in), file=log)
            print("    pdb-file not found, %s incorrectly returned" %(pdb_in))
            pdb_out_real = pdb_in
        print("----------------")

        if keywords.density_modification.density_modification:
            ccp4_dm = re.sub(r".mtz$", ".ccp4", mtz_dm)
            print("REAL SPACE REFINEMENT WITH %s AND %s" %(ccp4_dm, pdb_out_rec))
            pdb_out_rec_real = ref.coot_real_space_refinement_ccp4(pdb_out_rec, ccp4_dm)
            print("output real space refinement after reciprocal space refinement:", file=log)
            print("----------------")
            print("output real space refinement after reciprocal space refinement:")
            if os.path.isfile(pdb_out_rec_real):
                print("    pdb-file: %s"%(pdb_out_rec_real), file=log)
                print("    pdb-file: %s"%(pdb_out_rec_real))
            else:
                print("    pdb-file not found, %s incorrectly returned" %(pdb_in), file=log)
                print("    pdb-file not found, %s incorrectly returned" %(pdb_in))
                pdb_out_rec_real = pdb_out_rec
            print("----------------")
        else:
            print("REAL SPACE REFINEMENT WITH %s AND %s" %(mtz_out_rec, pdb_out_rec))
            pdb_out_rec_real = ref.coot_real_space_refinement_mtz(pdb_out_rec, mtz_out_rec, '2FOFCWT, PH2FOFCWT, FOFCWT, PHFOFCWT')
            print("output real space refinement after reciprocal space refinement:", file=log)
            print("----------------")
            print("output real space refinement after reciprocal space refinement:")
            if os.path.isfile(pdb_out_rec_real):
                print("    pdb-file: %s"%(pdb_out_rec_real), file=log)
                print("    pdb-file: %s"%(pdb_out_rec_real))
            else:
                print("    pdb-file not found, %s incorrectly returned" %(pdb_in), file=log)
                print("    pdb-file not found, %s incorrectly returned" %(pdb_in))
                pdb_out_rec_real = pdb_out_rec     
            print("----------------")
            
        return mtz_out_rec, pdb_out_rec, pdb_out_real, pdb_out_rec_real

               
class Filesandmaps(object):
    """
    Class to generate mtz files with structure factors; and mtz,ccp4 and xplor files with map coefficients
    """
    def __init__(self, miller_array, rfree, maptype, prefix, fmodel_ref, crystal_gridding = None):
        self.ms         = miller_array
        self.rfree      = rfree
        self.maptype    = maptype
        self.prefix     = prefix
        self.fmodel_ref = fmodel_ref
                 
        #self.r_free_flags = miller.array(miller_set=self.rfree, data=self.rfree.data().as_bool())

        self.fmodel = mmtbx.f_model.manager(
            f_obs          = miller_array,
            r_free_flags   = self.rfree,
            xray_structure = self.fmodel_ref.xray_structure)
        
        self.sites_cart = self.fmodel_ref.xray_structure.sites_cart()
        
        if crystal_gridding == None:
            self.crystal_gridding = self.ms.crystal_gridding(
                d_min             = self.ms.d_min(),
                symmetry_flags    = maptbx.use_space_group_symmetry,
                resolution_factor = 0.25)
        else:
            self.crystal_gridding = crystal_gridding
            
        #print("Crystal_gridding gridpoints:", self.crystal_gridding.n_grid_points())
        
        if self.maptype.lower() == 'qfofo':
            self.labels = {'data':"QFDIFF", 'map_coefs_diff': "QFOFOWT"}
        elif self.maptype.lower() == 'kfofo':
            self.labels = {'data':"KFDIFF", 'map_coefs_diff': "KFOFOWT"}
        elif self.maptype.lower() == 'fofo':
            self.labels = {'data':"FDIFF", 'map_coefs_diff': "FOFOWT"}
        elif self.maptype.lower() == 'qfextr':
            self.labels = {'data':"QFEXTR", 'map_coefs_map': "2QFEXTRFCWT", 'map_coefs_diff': "QFEXTRFCWT"}
        elif self.maptype.lower() == 'kfextr':
            self.labels = {'data':"KFEXTR", 'map_coefs_map': "2KFEXTRFCWT", 'map_coefs_diff': "KFEXTRFCWT"}
        elif self.maptype.lower() == 'fextr':
            self.labels = {'data':"FEXTR", 'map_coefs_map': "2FEXTRFCWT", 'map_coefs_diff': "FEXTRFCWT"}
        elif self.maptype.lower() == 'qfgenick':
            self.labels = {'data': "QFGENICK", 'map_coefs_map': "QFGENICKWT", 'map_coefs_diff': "QFGENICKFCWT"}
        elif self.maptype.lower() == 'kfgenick':
            self.labels = {'data': "KFGENICK", 'map_coefs_map': "KFGENICKWT", 'map_coefs_diff': "KFGENICKFCWT"}
        elif self.maptype.lower() == 'fgenick':
            self.labels = {'data': "FGENICK", 'map_coefs_map': "FGENICKWT", 'map_coefs_diff': "FGENICKFCWT"}
        elif self.maptype.lower() == 'qfextr_calc':
            self.labels = {'data':"QFEXTR_CALC", 'map_coefs_map': "2QFEXTR_CALCFCWT", 'map_coefs_diff': "QFEXTR_CALCFCWT"}
        elif self.maptype.lower() == 'kfextr_calc':
            self.labels = {'data':"KFEXTR_CALC", 'map_coefs_map': "2KFEXTR_CALCFCWT", 'map_coefs_diff': "KFEXTR_CALCFCWT"}
        elif self.maptype.lower() == 'fextr_calc': 
            self.labels = {'data':"FEXTR_CALC", 'map_coefs_map': "2FEXTR_CALCFCWT", 'map_coefs_diff': "FEXTR_CALCFCWT"}
        else:
            print('map type %s not allowed. Using default dummy labels for the outpur files.' %(self.maptype), file=log)
            print('map type %s not allowed. Using default dummy labels for the outpur files.' %(self.maptype))
            self.labels = {'data':"EXTR", 'map_coefs_map': "2EXTRFCWT", 'map_coefs_diff': "EXTRFCWT"}
            
    def mtz_file_Fs(self):
        """
        Generate mtz file with structure factors
        """
        mtz_dataset = self.ms.as_mtz_dataset(column_root_label=self.labels['data'])
        mtz_dataset.add_miller_array(miller_array = self.rfree, column_root_label="FreeR_flag")
        mtz_dataset.mtz_object().write(self.F_name)
        
    def mtz_mapcoefs(self, mc_map, mc_diff):
        """
        Generate mtz file with map coefficients. First column is 2mFo-DFc type, second column is mFo-DFc type
        """
        mtz_dataset = mc_map.as_mtz_dataset(column_root_label=self.labels['map_coefs_map'])
        mtz_dataset.add_miller_array(miller_array = mc_diff, column_root_label = self.labels['map_coefs_diff'])
        mtz_object = mtz_dataset.mtz_object()
        mtz_object.write(file_name = self.mtz_name)

    def ccp4_xplor_mapcoefs(self, mc_map, mc_diff):
        """
        Generate ccp4 and xplor files with map coefficients. First file is 2mFo-DFc type, second file is mFo-DFc type
        resolution factor =0.25 and buffer = 5.0 default values in mtz2map
        """
        fft_map_2mfodfc = mc_map.fft_map(
            crystal_gridding =self.crystal_gridding, resolution_factor=0.25).apply_sigma_scaling()
        #fft_map_2mfodfc.as_ccp4_map(file_name = self.ccp4_name_2FoFc) #Origin not correct, causes problems with pymol
        iotbx.map_tools.write_ccp4_map(
            sites_cart = self.sites_cart,
            unit_cell  = fft_map_2mfodfc.unit_cell(),
            map_data   = fft_map_2mfodfc.real_map(),
            n_real     = fft_map_2mfodfc.n_real(),
            buffer     = 5.0,
            file_name  = self.ccp4_name_2FoFc)
        #generation of xplor map for 2mfodfc type outcommented as it is not necessary for analysis but takes a lot of space
        #mmtbx.maps.utils.write_xplor_map(
            #sites_cart = self.sites_cart,
            #unit_cell  = fft_map_2mfodfc.unit_cell(),
            #map_data   = fft_map_2mfodfc.real_map(),
            #n_real     = fft_map_2mfodfc.n_real(),
            #buffer     = 5.0,
            #file_name  = self.xplor_name_2FoFc)
        
        fft_map_mfodfc = mc_diff.fft_map(
            crystal_gridding = self.crystal_gridding, resolution_factor=0.25).apply_sigma_scaling()
        #fft_map_mfodfc = mc_diff.fft_map(resolution_factor=0.25).apply_sigma_scaling()
        #fft_map_mfodfc.as_ccp4_map(file_name = self.ccp4_name_FoFc)
        iotbx.map_tools.write_ccp4_map(
            sites_cart = self.sites_cart,
            unit_cell  = fft_map_mfodfc.unit_cell(),
            map_data   = fft_map_mfodfc.real_map(),
            n_real     = fft_map_mfodfc.n_real(),
            buffer     = 5.0,
            file_name  = self.ccp4_name_FoFc)
 
        #fft_map_mfodfc.as_xplor_map(file_name = self.xplor_name_FoFc) 
        #Next is how xplor map are calculated in mtz2map although sites_cart are extracted from a pdb instead of fmodel. In our case coordinates are fractional and hence it is not working
        mmtbx.maps.utils.write_xplor_map( 
            sites_cart = self.sites_cart,
            unit_cell  = fft_map_mfodfc.unit_cell(),
            map_data   = fft_map_mfodfc.real_map(),
            n_real     = fft_map_mfodfc.n_real(),
            buffer     = 5.0,
            file_name  = self.xplor_name_FoFc)
        
    def ccp4_mapcoefs(self, mc_map, mc_diff):
        """
        Generate ccp4 files with map coefficients. First file is 2mFo-DFc type, second file is mFo-DFc type
        resolution factor =0.25 and buffer = 5.0 default values in mtz2map
        Same as ccp4_xplor_mapcoefs but no calculation of xplor maps
        """
        #Calculate 2mFo-DFc type of map
        fft_map_2mfodfc = mc_map.fft_map(
            crystal_gridding =self.crystal_gridding, resolution_factor=0.25).apply_sigma_scaling()
        #fft_map_2mfodfc.as_ccp4_map(file_name = self.ccp4_name_2FoFc) #Origin not correct, causes problems with pymol
        iotbx.map_tools.write_ccp4_map(
            sites_cart = self.sites_cart,
            unit_cell  = fft_map_2mfodfc.unit_cell(),
            map_data   = fft_map_2mfodfc.real_map(),
            n_real     = fft_map_2mfodfc.n_real(),
            buffer     = 5.0,
            file_name  = self.ccp4_name_2FoFc)
        
        #Calculate mFo-DFc type of map
        fft_map_mfodfc = mc_diff.fft_map(
            crystal_gridding = self.crystal_gridding, resolution_factor=0.25).apply_sigma_scaling()
        #fft_map_mfodfc = mc_diff.fft_map(resolution_factor=0.25).apply_sigma_scaling()
        #fft_map_mfodfc.as_ccp4_map(file_name = self.ccp4_name_FoFc)
        iotbx.map_tools.write_ccp4_map(
            sites_cart = self.sites_cart,
            unit_cell  = fft_map_mfodfc.unit_cell(),
            map_data   = fft_map_mfodfc.real_map(),
            n_real     = fft_map_mfodfc.n_real(),
            buffer     = 5.0,
            file_name  = self.ccp4_name_FoFc)

    def write_FoFo_output(self):
        """
        Write map coefficients for Fo-Fo in mtz, ccp4 and xplor format.
        Use of adapted verion of map_tools to combine fom from reference data with difference amplitudes
        """
        mtz_name   = '%s_m%s.mtz' %(self.prefix, self.maptype)
        ccp4_name  = '%s_m%s.ccp4' %(self.prefix, self.maptype)
        #xplor_name = '%s_m%s.map' %(self.prefix, self.maptype)
        
        edm = map_tools_Millerset.electron_density_map(fobs_in=self.ms, fmodel_2=self.fmodel_ref)
        mc_mfofo = edm.map_coefficients(map_type='mfo')
        mtz_dataset = mc_mfofo.as_mtz_dataset(column_root_label=self.labels['map_coefs_diff'])
        mtz_object = mtz_dataset.mtz_object()
        mtz_object.write(file_name = mtz_name)
        fft_map_mfofo = mc_mfofo.fft_map(
            crystal_gridding = self.crystal_gridding, resolution_factor=0.25).apply_sigma_scaling()
        #fft_map_mfofo.as_ccp4_map(file_name  = ccp4_name) #works too but below is how xplor-files arewritten with mtz2map
        iotbx.map_tools.write_ccp4_map(
          sites_cart = self.sites_cart,
          unit_cell  = fft_map_mfofo.unit_cell(),
          map_data   = fft_map_mfofo.real_map(),
          n_real     = fft_map_mfofo.n_real(),
          buffer     = 5.0,
          file_name  = ccp4_name)

        ##fft_map_mfofo.as_xplor_map(file_name = xplor_name) #works too but below is how xplor-files arewritten with mtz2map
        #mmtbx.maps.utils.write_xplor_map(
            #sites_cart = self.sites_cart,
            #unit_cell  = fft_map_mfofo.unit_cell(),
            #map_data   = fft_map_mfofo.real_map(),
            #n_real     = fft_map_mfofo.n_real(),
            #buffer     = 5.0,
            #file_name  = xplor_name)
        
        return mtz_name, ccp4_name #, xplor_name
        
    def write_Fextr_maps(self, fill_missing=True):
        """
        Write map coefficients to to mtz, ccp4 format.
        fom and D are calculated as usual.
        Same as write_Fextr_Fextr_calc_output but without writing the structure factor amplitudes.
        Xplor maps are not written since they have become obsolete in Xtrapol8.
        """
        self.mtz_name         = '%s_2m%s-DFc_m%s-DFc.mtz' %(self.prefix, self.maptype, self.maptype)
        self.ccp4_name_2FoFc  = '%s_2m%s-DFc.ccp4' %(self.prefix, self.maptype)
        self.ccp4_name_FoFc   = '%s_m%s-DFc.ccp4' %(self.prefix, self.maptype)
        
        fmodel_update = self.fmodel.deep_copy()
        print("\nUpdating scales for map calculations, can take a few minutes")
        try:
            fmodel_update.update_all_scales(remove_outliers=False, log=sys.stdout) #need update_all in order to update alpha
            fmodel_update.show()
        except RuntimeError:
            print("The model-based structure factors could not be scaled using the fast method. Try again with slow method")
            fmodel_update.update_all_scales(remove_outliers=False, log=sys.stdout, fast=False) #need update_all in order to update alpha
            fmodel_update.show()
        edm = mmtbx.map_tools.electron_density_map(fmodel=fmodel_update) #electron density map object created
        mc_mfodfc  = edm.map_coefficients(map_type='mfo-dfc', isotropize=True, fill_missing = False) #mfodfc mapcoefs defined
        mc_2mfodfc = edm.map_coefficients(map_type='2mfo-dfc', isotropize=True, fill_missing = fill_missing) #2mfodfc mapcoefs defined
        self.mtz_mapcoefs(mc_2mfodfc, mc_mfodfc)
        self.ccp4_mapcoefs(mc_2mfodfc, mc_mfodfc)
        
        return self.mtz_name, self.ccp4_name_2FoFc, self.ccp4_name_FoFc
    
    def write_Fextr_Fextr_calc_output(self, fill_missing=True):
        """
        Write structure factors of (q-weighted) Fetxr and Fextr_calc to mtz file, and calculate map coefficients to mtz, ccp4 and xplot format.
        fom and D are calculated as usual.
        """
        self.F_name           = '%s_%s.mtz' %(self.prefix, self.maptype)
        self.mtz_name         = '%s_2m%s-DFc_m%s-DFc.mtz' %(self.prefix, self.maptype, self.maptype)
        self.ccp4_name_2FoFc  = '%s_2m%s-DFc.ccp4' %(self.prefix, self.maptype)
        self.ccp4_name_FoFc   = '%s_m%s-DFc.ccp4' %(self.prefix, self.maptype)
        #self.xplor_name_2FoFc = '%s_2m%s-DFc.map' %(self.prefix, self.maptype)
        #self.xplor_name_FoFc  = '%s_m%s-DFc.map' %(self.prefix, self.maptype)
        
        self.mtz_file_Fs()
                
        fmodel_update = self.fmodel.deep_copy()
        print("\nUpdating scales for map calculations, can take a few minutes")
        try:
            fmodel_update.update_all_scales(remove_outliers=False, log=sys.stdout) #need update_all in order to update alpha
            fmodel_update.show()
        except RuntimeError:
            print("The model-based structure factors could not be scaled using the fast method. Try again with slow method")
            fmodel_update.update_all_scales(remove_outliers=False, log=sys.stdout, fast=False) #need update_all in order to update alpha
            fmodel_update.show()
        #fmodel_update.update(f_obs=self.fmodel.f_obs()) #to avoid Fextr being altered
        edm = mmtbx.map_tools.electron_density_map(fmodel=fmodel_update) #electron density map object created
        mc_mfodfc  = edm.map_coefficients(map_type='mfo-dfc', isotropize=True, fill_missing = False) #mfodfc mapcoefs defined
        mc_2mfodfc = edm.map_coefficients(map_type='2mfo-dfc', isotropize=True, fill_missing = fill_missing) #2mfodfc mapcoefs defined
        self.mtz_mapcoefs(mc_2mfodfc, mc_mfodfc)
        #self.ccp4_xplor_mapcoefs(mc_2mfodfc, mc_mfodfc)
        self.ccp4_mapcoefs(mc_2mfodfc, mc_mfodfc)
        
        return self.F_name, self.mtz_name, self.ccp4_name_2FoFc, self.ccp4_name_FoFc #, self.xplor_name_2FoFc, self.xplor_name_FoFc
        
    def write_Fgenick_output(self):
        """
        Write structure factors of (q-weighted) Fetxr and Fextr_calc to mtz file, and calculate map coefficients to mtz, ccp4 and xplot format.
        fom and D are calculated with adapted map_tools in order to combine fom from reference data with extrapolated structure factors, as described in Genick, 2007.
        """
        self.F_name           = '%s_%s.mtz' %(self.prefix, self.maptype)
        self.mtz_name         = '%s_m%s_m%s-DFc.mtz' %(self.prefix, self.maptype, self.maptype)
        self.ccp4_name_2FoFc  = '%s_m%s.ccp4' %(self.prefix, self.maptype)
        self.ccp4_name_FoFc   = '%s_m%s-DFc.ccp4' %(self.prefix, self.maptype)
        #self.xplor_name_2FoFc = '%s_m%s.map' %(self.prefix, self.maptype)
        #self.xplor_name_FoFc  = '%s_m%s-DFc.map' %(self.prefix, self.maptype)
        
        self.mtz_file_Fs()
        
        fmodel_update = self.fmodel.deep_copy()
        print("\nUpdating scales for map calculations, can take a few minutes")
        try:
            fmodel_update.update_all_scales(remove_outliers=False, log=sys.stdout) #need update all in order to update alpha
            fmodel_update.show()
        except RuntimeError:
            print("The model-based structure factors could not be scaled using the fast method. Try again with slow method")
            fmodel_update.update_all_scales(remove_outliers=False, log=sys.stdout, fast=False) #need update all in order to update alpha
            fmodel_update.show()
        fmodel_update.update(f_obs=self.fmodel.f_obs()) #to avoid Fextr being altered
        edm = map_tools_fomsource.electron_density_map(fmodel_1=fmodel_update, fmodel_2=self.fmodel_ref)
        mc_mfodfc  = edm.map_coefficients(map_type='mfo-dfc', isotropize=True, fill_missing= False)
        mc_mfo = edm.map_coefficients(map_type='mfo', isotropize=True, fill_missing= False)
        self.mtz_mapcoefs(mc_mfo, mc_mfodfc)
        #self.ccp4_xplor_mapcoefs(mc_mfo, mc_mfodfc)
        self.ccp4_mapcoefs(mc_mfo, mc_mfodfc)
        
        return self.F_name, self.mtz_name, self.ccp4_name_2FoFc, self.ccp4_name_FoFc #, self.xplor_name_2FoFc, self.xplor_name_FoFc

#def get_unique_id(id_length=20):
    #"""
    #Function to get a unique id based on UUID with length id_length
    #"""
    #if id_length > 36:
        #id_length == 36
    #return str(uuid.uuid4())[:id_length]

#def generate_log_name(time_stamp):
    #"""
    #Generate a unique name for the Xtrapol8 logfile.
    #A short uuid of 20 characters is added to the logfile name.
    #"""
    #uuid = get_unique_id(36)
    #logname = "%s_Xtrapol8_%s.log" %(time_stamp, uuid)
    
    #return logname

#def remove_unique_id_from_log():
    #"""
    #Remove the unqiue sequence from the log file
    #"""
    #index = log.name.find("Xtrapol8")+len("Xtrapol8")
    #new_name = log.name[:index]+".log"
    #os.rename(log.name, new_name)
            
def run(args):
    
    version.VERSION
    now = datetime.now().strftime('%Y-%m-%d_%Hh%M')
    print('-----------------------------------------')
    print("Xtrapol8 -- version %s -- run date: %s" %(version.VERSION, now))
    print("Phenix version: {:s}".format(get_phenix_version()))
    print("CCP4 version: {:s}".format(get_ccp4_version()))

    #If no input, show complete help, should be changed in order to give help depending on the attribute level
    if len(args) == 0 :
        print('-----------------------------------------')
        master_phil.show(attributes_level=1)
        raise Usage("phenix.python Fextr.py + [.phil] + [arguments]\n arguments only overwrite .phil if provided last")
    
    #Generate log-file. Needs to be created before the output directory is created and to be a global parameter in order to be easily used in all classes and functions
    logname = generate_log_name(now)
    global log
    log = open(logname, "w")
    print("Xtrapol8 -- version %s -- run date: %s" %(version.VERSION, now), file=log)
    print("Phenix version: {:s}".format(get_phenix_version()), file=log)
    print("CCP4 version: {:s}".format(get_ccp4_version()), file=log)

    log_dir = os.getcwd()
    
    #Extract input from inputfile and command line
    argument_interpreter = master_phil.command_line_argument_interpreter(home_scope="input")
    input_objects = iotbx.phil.process_command_line_with_files(
        args=args,
        master_phil=master_phil
        )
    params = input_objects.work.extract()
    #modified_phil = master_phil.format(python_object=params)
    
    remarks = []
    #Check if non-phenix programs can be found:
    if check_program_path('coot')[1] == False:
        remark = "COOT not found."
        remarks.append(remark)
    if check_program_path('scaleit')[1] == False:
        remark = "scaleit not found. Data will not be scaled."
        remarks.append(remark)
        params.scaling.b_scaling = 'no'
    if params.refinement.use_refmac_instead_of_phenix:
        if (check_program_path('refmac5')[1] and check_program_path('coot')[1]) == False:
            remark = "refmac and/or COOT not found. Phenix will be used for refinement."
            remarks.append(remark)
            params.refinement.use_refmac_instead_of_phenix = False
            
    #specify extrapolated structure factors and map types
    qFextr_map = qFgenick_map = qFextr_calc_map = Fextr_map = Fgenick_map = Fextr_calc_map = kFextr_map = kFgenick_map = kFextr_calc_map = False
    #calculate only the specified map types
    if 'qfextr' in params.f_and_maps.f_extrapolated_and_maps:
        qFextr_map      = True
    if 'fextr' in params.f_and_maps.f_extrapolated_and_maps:
        Fextr_map       = True
    if 'kfextr' in params.f_and_maps.f_extrapolated_and_maps:
        kFextr_map      = True
    if 'qfgenick' in params.f_and_maps.f_extrapolated_and_maps:
        qFgenick_map       = True
    if 'kfgenick' in params.f_and_maps.f_extrapolated_and_maps:
        kFgenick_map       = True
    if 'fgenick' in params.f_and_maps.f_extrapolated_and_maps:
        Fgenick_map        = True
    if ('qfextr_calc') in params.f_and_maps.f_extrapolated_and_maps:
        qFextr_calc_map = True
    if ('kfextr_calc') in params.f_and_maps.f_extrapolated_and_maps:
        kFextr_calc_map = True
    if ('fextr_calc') in params.f_and_maps.f_extrapolated_and_maps:
        Fextr_calc_map  = True
            
    if params.f_and_maps.fofo_type == 'fofo':
        qFoFo_weight = False
        kFoFo_weight = False
    elif params.f_and_maps.fofo_type == 'qfofo':
        qFoFo_weight = True
        kFoFo_weight = False
    elif params.f_and_maps.fofo_type == 'kfofo':
        qFoFo_weight = False
        kFoFo_weight = True
    else:
        remark = "%s not defined. Q-weighting will be applied." %(params.f_and_maps.fofo_type.__phil_path__())
        remarks.append(remark)
        qFoFo_weight = True
        kFoFo_weight = False

    #Only run the non-weighted
    if params.f_and_maps.only_no_weight:
        qFextr_map = qFgenick_map = qFextr_calc_map = kFextr_map = kFgenick_map = kFextr_calc_map = False
        Fextr_map  = Fgenick_map  = Fextr_calc_map = True
    #Only run the k-weighted
    if params.f_and_maps.only_kweight:
        qFextr_map = qFgenick_map = qFextr_calc_map = Fextr_map = Fgenick_map = Fextr_calc_map = False
        kFextr_map = kFgenick_map = kFextr_calc_map = True
    #Only run the q-weighted
    if params.f_and_maps.only_qweight: 
        Fextr_map = Fgenick_map = Fextr_calc_map = kFextr_map = kFgenick_map = kFextr_calc_map = False
        qFextr_map = qFgenick_map = qFextr_calc_map = True
    #Run all maps
    if params.f_and_maps.all_maps: #calculate all Fextr map types
        qFextr_map = qFgenick_map = qFextr_calc_map = Fextr_map = Fgenick_map = Fextr_calc_map = kFextr_map = kFgenick_map = kFextr_calc_map = True
        
    #convert the old use_occupancy_from_distance_analysis to the new occupancy_estimation keyword
    if params.map_explorer.use_occupancy_from_distance_analysis == True:
        remark = "map_explorer.use_occupancy_from_distance_analysis will become obsolete. Use map_explorer.occupancy_estimation=distance_analysis in the future."
        remarks.append(remark)
        params.map_explorer.occupancy_estimation = "distance_analysis"
        
    #change occupancy estimation method to the difference_map_maximization if distance_analysis is set but refinement not run
    if (params.refinement.run_refinement == False and params.map_explorer.occupancy_estimation == 'distance_analysis'):
        remark = "Distance_analysis cannot be carried out because refinement.run_refinement = False. The difference_map_maximization method will be used."
        remarks.append(remark)
        params.map_explorer.occupancy_estimation = 'difference_map_maximization'
        
    #if all map types being false:
    if qFextr_map == qFgenick_map == qFextr_calc_map == Fextr_map == Fgenick_map == Fextr_calc_map == kFextr_map == kFgenick_map == kFextr_calc_map == False and params.f_and_maps.fast_and_furious == False:
        remark = 'The combination of arguments used to define extrapolated structure factors and maps leads to no calculations at all. The default will be applied: qFextr.'
        remarks.append(remark)
        qFextr_map = True
    #if fast_and_furious mode: overwrite all F and map setting to default:
    if params.f_and_maps.fast_and_furious:
        remark = "fast_and_furious mode: some parameters will be reset to their default values."
        remarks.append(remark)
        #change parameters for Xtrapol8_out.phil
        params.f_and_maps.fofo_type = 'qfofo'
        params.f_and_maps.f_extrapolated_and_maps = ['qfextr']
        params.f_and_maps.only_no_weight = params.f_and_maps.all_maps = params.f_and_maps.only_kweight = params.f_and_maps.only_qweight = False
        params.map_explorer.use_occupancy_from_distance_analysis = False
        if params.map_explorer.occupancy_estimation == 'distance_analysis':
            params.map_explorer.occupancy_estimation = 'difference_map_maximization'
        params.f_and_maps.negative_and_missing='truncate_and_fill'
        #change other working parameters
        qFoFo_weight = qFextr_map = True
        kFoFo_weight = qFgenick_map = qFextr_calc_map = Fextr_map = Fgenick_map = Fextr_calc_map = kFextr_map = kFgenick_map = kFextr_calc_map = False

                
    #Bring all maptypes to be calculated together in list instead of using loose varaibles:
    all_maptypes   = ['qFextr_map','Fextr_map', 'qFgenick_map', 'Fgenick_map', 'qFextr_calc_map', 'Fextr_calc_map', 'kFextr_map', 'kFgenick_map', 'kFextr_calc_map']
    all_maps       = [qFextr_map, Fextr_map, qFgenick_map, Fgenick_map, qFextr_calc_map, Fextr_calc_map, kFextr_map, kFgenick_map, kFextr_calc_map]
    maptypes_zip   = zip(all_maptypes, all_maps)
    final_maptypes = [mp[0] for mp in maptypes_zip if mp[1] == True]
    
    #print("final_maptypes", final_maptypes)
    
    if qFoFo_weight:
        params.f_and_maps.fofo_type = 'qfofo'
    elif kFoFo_weight:
        params.f_and_maps.fofo_type = 'kfofo'
    else:
        params.f_and_maps.fofo_type = 'fofo'
    
    #get list with occupancies from start-end-steps or list
    occ_step  = (params.occupancies.high_occ - params.occupancies.low_occ)/params.occupancies.steps
    if params.occupancies.list_occ == None:
        occ = params.occupancies.low_occ
        occ_lst = []
        while occ <= params.occupancies.high_occ:
            occ_lst.append(occ)
            occ+=occ_step
    else:
        occ_lst = params.occupancies.list_occ
    occ_lst.sort()
    if len(occ_lst) == 0:
        remark = "No input occupancies found. Xtrapol8 will stop after the FoFo calculation."
        remarks.append(remark)
        params.output.generate_fofo_only = True
        params.occupancies.list_occ = None
    else:
        params.occupancies.list_occ = occ_lst
    ################################################################
    if len(remarks) > 0:
        print('-----------------------------------------', file=log)
        print("REMARKS", file=log)
        print('-----------------------------------------', file=log)
        print("\n".join(remarks), file=log)
        
        print('-----------------------------------------')
        print("REMARKS")
        print('-----------------------------------------')
        print("\n".join(remarks))
        #############################################################
    
    #Add all arguments to log-file
    print('-----------------------------------------', file=log)
    print("NON DEFAULT ARGUMENTS", file=log)
    print('-----------------------------------------', file=log)
    
    print('-----------------------------------------')
    print("NON DEFAULT ARGUMENTS")
    print('-----------------------------------------')

    modified_phil = master_phil.format(python_object=params)
    #modified_phil.show(out=log)
    #get the differences with the default values and only show these in the log-file
    diff_phil = master_phil.fetch_diff(source=modified_phil)
    diff_phil.show()
    diff_phil.show(out=log)

    ################################################################
    print('-----------------------------------------')
    print('DATA PREPARATION')
    print('-----------------------------------------')
    
    print('-----------------------------------------', file=log)
    print('DATA PREPARATION', file=log)
    print('-----------------------------------------', file=log)

    DH = DataHandler(params.input.reference_pdb, params.input.reference_mtz, params.input.additional_files, params.output.outdir, params.input.triggered_mtz)

    
    #Check if all input files exists and are of correct type
    err , err_m = DH.check_all_files()
    if err == 1:
        print(err_m, file=log)
        print(err_m)
        print("Check input files and rerun")
        sys.exit()
            
    #get output name, or take prefix of triggered mtz or take dummy name if name is too long.
    if params.output.outname == None:
        params.output.outname = get_name(DH.mtz_on)
    if len(params.output.outname) > 50:
        outname = 'triggered'
        print("output.outname='%s' is too long. It will be substituted by '%s' during the excecution of Xtrapol8 and at the end the file names will be converted to the real output names. Please take this into account when inspecting log files." %(params.output.outname, outname), file=log)
        print('---------------------------', file=log)
        print("output.outname='%s' is too long. It will be substituted by '%s' during the excecution of Xtrapol8 and at the end the file names will be converted to the real output names. Please take this into account when inspecting log files." %(params.output.outname, outname))
        print('---------------------------')
    else:
        outname = params.output.outname

    #Extract space group and unit cell from model pdb file
    DH.get_UC_and_SG()
    
    #Check if all cif files are given and extract the three-letter codes
    print("----Check additional files----")
    DH.check_additional_files()
    print('---------------------------')

    #change the names so that they appear as absolute path in the Xtrapol8_out:
    params.input.reference_pdb = DH.pdb_in #This should already be the absolute path
    params.input.reference_mtz = os.path.abspath(DH.mtz_off)
    params.input.triggered_mtz = os.path.abspath(DH.mtz_on)
    params.output.outdir = DH.outdir #This should already be the absolute path
    params.input.additional_files = list(map(lambda x: os.path.abspath(x), params.input.additional_files))

    #change to output directory
    startdir = os.getcwd()
    outdir = DH.outdir
    os.chdir(outdir)
    
    #Move log file to output directory
    full_log = "%s/%s" %(log_dir, log.name)
    if os.path.isfile(full_log):
        shutil.move(full_log, full_log.replace(log_dir,outdir))
    log_name = remove_unique_id_from_log(log.name)
    full_log = "%s/%s" %(outdir, log_name)
    
    #extract columns from mtz files that needs to be substracted
    print("----Column extraction from reflection files----")
    DH.extract_fobs(params.input.low_resolution,params.input.high_resolution)
    print('---------------------------')
    
    # compatibilty test between mtz-files and model
    print('Extracting the unit cell and space group from input files. Check if these are similar for all input files:')
    symm_manager = SymManager(DH.model_in)
    err_m_off, err_off = symm_manager.check_symm(DH.reflections_off)
    err_m_on, err_on = symm_manager.check_symm(DH.reflections_on)
    
    if err_off + err_on == 0:
        print('Space group and unit cell compatible', file=log)
        print('Space group and unit cell compatible')
    else:
        print('Reference:', err_m_off, file=log)
        print('Triggered:', err_m_on, file=log)
        print('Reference:', err_m_off)
        print('Triggered:', err_m_on)
        symm_manager.show(out=log)
        if err_on == 1:
            #err_on == 1 while err_off == 0 can only happen if pointless failed the reindexing
            print("Changing space group and/or unit cell parameters of the triggered data in order to match with the model (if not already done by Pointless).", file=log)
            print("Changing space group and/or unit cell parameters of the triggered data in order to match with the model (if not already done by Pointless).")
            DH.fobs_on = make_miller_array(DH.fobs_on.data(),
                                           DH.fobs_on.sigmas(),
                                           DH.SG,
                                           DH.UC,
                                           DH.fobs_on.indices())
        if err_off == 1: #very weird if this would happen
            print("Changing space group and/or unit cell parameters of the reference data in order to match with the model. This is weird! Was the model generated after refinement in the reference mtz??", file=log)
            print("Changing space group and/or unit cell parameters of the reference data in order to match with the model. This is weird! Was the model generated after refinement in the reference mtz??")
            DH.fobs_off = make_miller_array(DH.fobs_off.data(),
                                            DH.fobs_off.sigmas(),
                                            DH.SG,
                                            DH.UC,
                                            DH.fobs_off.indices())
            
            # As a consequence of the running pointless, the triggered data has been reindexed in the reference space group and thus needs to be altered too.
            print("Changing space group and/or unit cell parameters of the triggered data in order to match with the model as these might have been changed by Pointless.", file=log)
            print("Changing space group and/or unit cell parameters of the triggered data in order to match with the model as these might have been changed by Pointless.")
            DH.fobs_on = make_miller_array(DH.fobs_on.data(),
                                           DH.fobs_on.sigmas(),
                                           DH.SG,
                                           DH.UC,
                                           DH.fobs_on.indices())
    
    
    #Write all input paramters to a phil file.
    modified_phil.show(out=open("Xtrapol8_in.phil", "w"))
    if params.output.generate_phil_only:
        params.output.GUI = False
        log.close()
        sys.exit()
    
    #Print reflection statistics to screen and log-file
    print("----Summary reference data:----")
    print(DH.fobs_on.show_comprehensive_summary())
    print("----Summary reference data:----", file=log)
    print(DH.fobs_on.show_comprehensive_summary(f=log), file=log)
    print("----Summarry triggered data:----")
    print(DH.fobs_off.show_comprehensive_summary())
    print("------------------------------------")
    print("----Summarry triggered data:----", file=log)
    print(DH.fobs_off.show_comprehensive_summary(f=log), file=log)
    print("------------------------------------", file=log)

    #Scale the data sets. Different steps taken:
    #1) generate a set of Rfree reflections
    #2) generate fmodel with reference reflections, model_in and Rfree flag
    #3) update fmodel with bulk-solvent scaling, i.e. scale fobs_off with fcalc
    #4) print R-factors and correlation coefficients between fobs_off and fcalc (this are the traditional Rwork/Rfree values)
    #5) check the reflections agreement between the reference and triggered reflections (and all other generated arrays)
    #6) scale triggered reflections with reference reflections
    #7) check the reflections agreement between the reference and triggered reflections (and all other generated arrays)
    #8) update fmodel with scaled reference and remove reflections that were rejected during the triggered scaling
    #9) print R-factors and correlation coefficients between fobs_off and fcalc which an indicator for isomorfism (this is the famous Riso)
    print("----Scaling Fmodel with Freference----", file=log)
    print("----Scaling Fmodel with Freference----")
    DH.generate_Rfree(DH.fobs_off, 0.05)
    print("fmodel generation, can take a few minutes")
    #Add scattering_table argument to make f_model with adapted scattering table
    DH.generate_f_model(DH.fobs_off, scattering_table=params.scattering_table)
    #DH.fmodel.show()
    if params.scaling.b_scaling != "no":
        print("Updating all fmodel scales.")
        try:
            DH.fmodel.update_all_scales(show=True)#, log=log)
        except RuntimeError:
            print("Fast method failed. Try again with slow method. This may lead to wrong scaling.")
            DH.fmodel.update_all_scales(show=True,
                                        fast=False)
    else: #only update the Fmodel part but keep the Fobs as they were
        print("Updating only the Fmodel part of the fmodel scales.")
        try:
            DH.fmodel.update_all_scales(update_f_part1=True,
                                        remove_outliers=False,
                                        bulk_solvent_and_scaling=True,
                                        apply_scale_k1_to_f_obs=False,
                                        show=True)
        except RuntimeError:
            print("Fast method failed. Try again with slow method. This may lead to wrong scaling.")
            DH.fmodel.update_all_scales(update_f_part1=True,
                                        remove_outliers=False,
                                        bulk_solvent_and_scaling=True,
                                        apply_scale_k1_to_f_obs=False,
                                        fast=False,
                                        show=True)
    #DH.fmodel.show()
    DH.fmodel.info().show_rfactors_targets_scales_overall(out=sys.stdout)
    print("Fobs,reference and Fcalc,reference scaled using mmtbx f_model")
    print('R_work:', DH.fmodel.r_work())
    print('R_free:', DH.fmodel.r_free())
    print("Fobs,reference and Fcalc,reference scaled using mmtbx f_model", file=log)
    print('R_work:', DH.fmodel.r_work(), file=log)
    print('R_free:', DH.fmodel.r_free(), file=log)
    DH.fobs_on = DH.get_common_indices_and_Fobs_off(DH.fobs_on) #compare reflections and reassemble off_state data set after scaling of f_obs_off 
    print("----Scaling Ftriggered with Freference----", file=log)
    print("----Scaling Ftriggered with Freference----")
    DH.scale_fobss(params.scaling.b_scaling,
                params.scaling.low_resolution,
                params.scaling.high_resolution)
    #update the parameters so that they appear correct in the output phil files
    params.scaling.high_resolution = DH.scaling_dmin
    params.scaling.low_resolution  = DH.scaling_dmax
    DH.fobs_on_scaled = DH.get_common_indices_and_Fobs_off(DH.fobs_on_scaled) #compare reflections and reassemble off_state data set after scaling of f_obs_on
    DH.update_fmodel(DH.fobs_off_scaled) #alter fmodel to remove reflections that were removed during fobs_on scaling
    riso, cciso = compute_r_factors(DH.fobs_off_scaled, DH.fobs_on_scaled, DH.rfree, log=log)
    print("Overall R_iso = %.4f " %(riso))
    print("Overall cc_iso = %.4f " %(cciso))
    print("Overall R_iso = %.4f " %(riso), file=log)
    print("Overall cc_iso = %.4f " %(cciso), file=log)


    #Extract the minimum and maximum resolution, as defined by the input parameters and/or the common reflections between the two scaled datasets
    dmax, dmin = DH.fobs_off_scaled.d_max_min()
    
    #set the params.resolution boundaries so that they appear in the output.phil
    params.input.high_resolution = dmin
    params.input.low_resolution = dmax

    #Print scaled reflections statistics to screen and log-file. Stupif that I'm printing everything to the screen and to the log-file. There must be something more elegant
    print("---------------------------", file=log)
    print("Resolution range= %2f - %2f A" %(dmax, dmin), file=log)
    print("---------------------------", file=log)
    print('---------------------------')
    print("Resolution range= %2f - %2f A" %(dmax, dmin))
    print('---------------------------')

    
    print("----Summary triggered scaled and sorted common reflections with reference data:----")
    print(DH.fobs_on_scaled.show_comprehensive_summary())
    print("----Summary triggered scaled and sorted common reflections with reference data:----", file=log)
    print(DH.fobs_on_scaled.show_comprehensive_summary(f=log), file=log)
    print("----Summary reference scaled and sorted common reflections with triggered data:----")
    print(DH.fobs_off_scaled.show_comprehensive_summary())
    print("----Summary reference scaled and sorted common reflections with triggered data:----", file=log)
    print(DH.fobs_off_scaled.show_comprehensive_summary(f=log), file=log)
    #time.sleep(3)

    print('-----------------------------------------')
    print('DATA PREPARATION DONE')
    print('-----------------------------------------')
    
    ################################################################
    print("CALCULATE (WEIGHTED) FO-FO FOURIER DIFFERENCE MAP")
    print('-----------------------------------------')
    
    print('-----------------------------------------', file=log)
    print("CALCULATE (WEIGHTED) FO-FO FOURIER DIFFERENCE MAP", file=log)
    print('-----------------------------------------', file=log)
        
    #calculate Fo-Fo and write map coefficients to mtz, ccp4
    FoFo = FobsFobs(DH.fobs_on_scaled, DH.fobs_off_scaled)
    FoFo.calculate_fdiff(kweight_scale = params.f_and_maps.kweight_scale)
    FoFo.write_maps(DH.fmodel, DH.rfree, outname, qweighting=qFoFo_weight, kweighting=kFoFo_weight)
    #print("%s maps generated in mtz,ccp4 and xplor format"%(params.f_and_maps.fofo_type), file=log)
    ##print("------------------------------------", file=log)
    #print("%s maps generated in mtz,ccp4 and xplor format"%(params.f_and_maps.fofo_type))
    ##print("------------------------------------")
    
    ################################################################

    #Use ccp4 map to find and integrate the peaks, annotate the peaks to residues and keep a list with the most important residues only based on a user-defined Z-score level
    print("\n************Map explorer************", file=log)
    print("\n************Map explorer************")
    if params.map_explorer.radius == None:
        params.map_explorer.radius = dmin
    map_expl_out_FoFo, residlist_zscore, mask, FoFo_ref = map_explorer(FoFo.ccp4_name, DH.pdb_in, params.map_explorer.radius, params.map_explorer.peak_integration_floor, params.map_explorer.peak_detection_threshold, params.map_explorer.z_score)
    #residlist_zscore  = Map_explorer_analysis(peakintegration_file = map_expl_out_FoFo,log=log).residlist_top(Z=params.map_explorer.z_score)
    print("FoFo map explored. Results in %s, residue list in residlist.txt and residues associated to highestpeaks in %s\n"
          %(map_expl_out_FoFo, residlist_zscore), file=log)
    print("FoFo map explored. Results in %s, residue list in residlist.txt and residues associated to highestpeaks in %s\n"
          %(map_expl_out_FoFo, residlist_zscore))
    #print('---------------------------')
    map_expl_out_FoFo = os.path.abspath(check_file_existance(map_expl_out_FoFo))
    residlist_zscore  = os.path.abspath(check_file_existance(residlist_zscore))
    
    #Generate plot which indicates the secondary structure and integrated peak volume
    print("----Generate difference map plot----")
    Map_explorer_analysis(peakintegration_file = map_expl_out_FoFo, ligands = DH.extract_ligand_codes(),log=log).get_ss(DH.pdb_in)
    #print('---------------------------')
    
    #Set parameter for whether or not q_weighting is applied on the FoFo calculation. This will be used later in extrapolated map calculation.
    #k-weighting is not yet set, this might pose problems
    if qFoFo_weight == True:
        FoFo_type = 'qFo-Fo'
    elif kFoFo_weight == True:
        FoFo_type = 'kFo-Fo'
    else:
        FoFo_type = 'Fo-Fo'
    
    print('-----------------------------------------')
    print("CALCULATE (WEIGHTED) FO-FO FOURIER DIFFERENCE MAP DONE")
    print('-----------------------------------------')
    
    if params.output.generate_fofo_only:
        print('-----------------------------------------')
        print("DONE! FOURIER DIFFERENCE MAP CALCULATED AND ANALYSIS PERFORMED.")
        print('-----------------------------------------')
        
        print('-----------------------------------------', file=log)
        print("DONE! FOURIER DIFFERENCE MAP CALCULATED AND ANALYSIS PERFORMED.", file=log)
        print('-----------------------------------------',file=log)
        #change names to real output name in case the dummy name was used
        if outname == 'triggered':
            #print("Replacing the dummpy outname ('%s') with true outname('%s')" %(outname, params.output.outname), file=log)
            print("Replacing the dummpy outname ('%s') with true outname('%s')" %(outname, params.output.outname))
            tr = [os.path.join(root, fle) for root, dirs, files in os.walk(outdir) for fle in files if outname in fle]
            _ = [os.rename(fle, fle.replace(outname, params.output.outname)) for fle in tr]
            FoFo.mtz_name = FoFo.mtz_name.replace(outname, params.output.outname)
        
        script_coot = open_all_in_coot("%s/%s"%(outdir,FoFo.mtz_name), [DH.pdb_in], [], DH.additional, outdir, FoFo_type)
        
        modified_phil = master_phil.format(python_object=params)
        modified_phil.show(out=open("Xtrapol8_out.phil", "w"))
        
        log.close()
        if (params.output.open_coot and check_program_path('coot')[1]):
            os.system("coot --script %s" %(script_coot))
            
        sys.exit()
    
    ################################################################
    print('-----------------------------------------')
    print("CALCULATE (Q/K-WEIGHTED) FEXTRAPOLATED MAPS")
    print('-----------------------------------------')
    
    print('-----------------------------------------', file=log)
    print("CALCULATE (Q/K-WEIGHTED) FEXTRAPOLATED MAPS", file=log)
    print('-----------------------------------------', file=log)

    #For each map type, keep track of the map_explorer files. Therefore generate lists with output file names. Not elegant, but works.
    if qFextr_map:
        qFextr_map_expl_fles  = [FoFo_ref]
        qFextr_recref_mtz_lst = []
        qFextr_recref_pdb_lst = [DH.pdb_in]
        qFextr_realref_lst    = [DH.pdb_in]
        qFextr_recrealref_lst = [DH.pdb_in]
    if kFextr_map:
        kFextr_map_expl_fles  = [FoFo_ref]
        kFextr_recref_mtz_lst = []
        kFextr_recref_pdb_lst = [DH.pdb_in]
        kFextr_realref_lst    = [DH.pdb_in]
        kFextr_recrealref_lst = [DH.pdb_in]
    if Fextr_map:
        Fextr_map_expl_fles  = [FoFo_ref]
        Fextr_recref_mtz_lst = []
        Fextr_recref_pdb_lst = [DH.pdb_in]
        Fextr_recrealref_lst = [DH.pdb_in]
        Fextr_realref_lst    = [DH.pdb_in]
    if qFgenick_map:
        qFgenick_map_expl_fles  = [FoFo_ref]
        qFgenick_recref_mtz_lst = []
        qFgenick_recref_pdb_lst = [DH.pdb_in]
        qFgenick_realref_lst    = [DH.pdb_in]
        qFgenick_recrealref_lst = [DH.pdb_in]
    if kFgenick_map:
        kFgenick_map_expl_fles  = [FoFo_ref]
        kFgenick_recref_mtz_lst = []
        kFgenick_recref_pdb_lst = [DH.pdb_in]
        kFgenick_realref_lst    = [DH.pdb_in]
        kFgenick_recrealref_lst = [DH.pdb_in]
    if Fgenick_map:
        Fgenick_map_expl_fles  = [FoFo_ref]
        Fgenick_recref_mtz_lst = []
        Fgenick_recref_pdb_lst = [DH.pdb_in]
        Fgenick_realref_lst    = [DH.pdb_in]
        Fgenick_recrealref_lst = [DH.pdb_in]
    if qFextr_calc_map:
        qFextr_calc_map_expl_fles  = [FoFo_ref]
        qFextr_calc_recref_mtz_lst = []
        qFextr_calc_recref_pdb_lst = [DH.pdb_in]
        qFextr_calc_realref_lst    = [DH.pdb_in]
        qFextr_calc_recrealref_lst = [DH.pdb_in]
    if kFextr_calc_map:
        kFextr_calc_map_expl_fles  = [FoFo_ref]
        kFextr_calc_recref_mtz_lst = []
        kFextr_calc_recref_pdb_lst = [DH.pdb_in]
        kFextr_calc_realref_lst    = [DH.pdb_in]
        kFextr_calc_recrealref_lst = [DH.pdb_in]
    if Fextr_calc_map:
        Fextr_calc_map_expl_fles  = [FoFo_ref]
        Fextr_calc_recref_mtz_lst = []
        Fextr_calc_recref_pdb_lst = [DH.pdb_in]
        Fextr_calc_realref_lst    = [DH.pdb_in]
        Fextr_calc_recrealref_lst = [DH.pdb_in]
    #fast_and_furious mode: no refinement, but need to keep track of structure factor files
    if (params.f_and_maps.fast_and_furious): Fextr_mtz_lst = []
    
    #Remove pickle file with Festr stats because otherwise plot will consist results from previous runs
    if os.path.isfile('%s/Fextr_binstats.pickle' %(outdir)):
        os.remove('%s/Fextr_binstats.pickle' %(outdir))
    #Remove pickle file with Number of negatove reflections because otherwise plot will consist results from previous runs
    if os.path.isfile('%s/Fextr_negative.pickle' %(outdir)):
        os.remove('%s/Fextr_negative.pickle' %(outdir))
    #Remove python movie movie because otherwise old files will be shown in the pymol session
    if os.path.isfile('%s/pymol_movie.py' %(outdir)):
        os.remove('%s/pymol_movie.py' %(outdir))

    fofo_data = ccp4_map.map_reader(file_name=FoFo.ccp4_name).data.as_numpy_array()
    ################################################################
    #Loop over occupancies and calculate extrapolated structure factors
    for occ in params.occupancies.list_occ:
        Fextr = Fextrapolate(FoFo.fdif,
                             FoFo.fdif_q,
                             FoFo.fdif_k,
                             FoFo.sigf_diff,
                             FoFo.q,
                             FoFo.k,
                             DH.fobs_off_scaled,
                             DH.fobs_on_scaled,
                             DH.fmodel,
                             DH.rfree,
                             occ,
                             name_out         = outname,
                             neg_refl_handle  = params.f_and_maps.negative_and_missing,
                             crystal_gridding = FoFo.get_crystal_gridding())
        
        #Results stored in folder depending on occupancy and whether q-weighting is applied.
        new_dirpath_q, new_dirpath_k, new_dirpath = Fextr.create_output_dirs(outdir)
        
        #for each F_and_map type calculate structure factors and maps
        for mp in final_maptypes:
            print(mp)
            if mp in ('qFextr_map','qFgenick_map','qFextr_calc_map'):
                os.chdir(new_dirpath_q)
            elif mp in ('kFextr_map','kFgenick_map','kFextr_calc_map'):
                os.chdir(new_dirpath_k)
            else:
                os.chdir(new_dirpath)
            
            #Depending on maptype, do following steps:
            #1) calculate the structure factors, write out to mtz file, handle negative reflections and generate associated plots or write to pickle file for later usuage
            #   calculate map coefficients and write to mtz and ccp4 files (latter only for mFo-DFc type)
            #2) get stats and write to pickle file in order to generate plot afterwards with all map types and occupancies
            #3) compute signal to noise and plot
            #As same steps are repeated, but with a sifferent Fextr-function, I should think of a more clever way to reduce the redundancy here (and in following if clauses)
            if mp == 'qFextr_map':
                Fextr.fextr(qweight=True, kweight=False, outdir_for_negstats = outdir)
                get_Fextr_stats(occ, Fextr.fextr_ms, Fextr.maptype, FoFo.fdif_q_ms, FoFo_type, outdir)
                compute_f_sigf(Fextr.fextr_ms, '%s' %(Fextr.maptype), log=log)
                #cc_list.append(plot_F1_F2(DH.fobs_off_scaled,Fextr.fextr_ms, F1_name = "Freference",F2_name = "Fextr"))
            elif mp == 'qFgenick_map':
                Fextr.fgenick(qweight=True, kweight=False,outdir_for_negstats = outdir)
                get_Fextr_stats(occ, Fextr.fgenick_ms, Fextr.maptype, FoFo.fdif_q_ms, FoFo_type, outdir)
                compute_f_sigf(Fextr.fgenick_ms, '%s' %(Fextr.maptype), log=log)
            elif mp == 'qFextr_calc_map':
                Fextr.fextr_calc(qweight=True, kweight=False,outdir_for_negstats = outdir)
                get_Fextr_stats(occ, Fextr.fextr_calc_ms, Fextr.maptype, FoFo.fdif_q_ms, FoFo_type, outdir)
                compute_f_sigf(Fextr.fextr_calc_ms, '%s' %(Fextr.maptype), log=log)
            elif mp == 'kFextr_map':
                Fextr.fextr(qweight=False, kweight=True, outdir_for_negstats = outdir)
                get_Fextr_stats(occ, Fextr.fextr_ms, Fextr.maptype, FoFo.fdif_k_ms, FoFo_type, outdir)
                compute_f_sigf(Fextr.fextr_ms, '%s' %(Fextr.maptype), log=log)
            elif mp == 'kFgenick_map':
                Fextr.fgenick(qweight=False, kweight=True, outdir_for_negstats = outdir)
                get_Fextr_stats(occ, Fextr.fgenick_ms, Fextr.maptype, FoFo.fdif_k_ms, FoFo_type, outdir)
                compute_f_sigf(Fextr.fgenick_ms, '%s' %(Fextr.maptype), log=log)
            elif mp == 'kFextr_calc_map':
                Fextr.fextr_calc(qweight=False, kweight=True, outdir_for_negstats = outdir)
                get_Fextr_stats(occ, Fextr.fextr_calc_ms, Fextr.maptype, FoFo.fdif_k_ms, FoFo_type, outdir)
                compute_f_sigf(Fextr.fextr_calc_ms, '%s' %(Fextr.maptype), log=log)
            elif mp == 'Fextr_map':
                Fextr.fextr(qweight=False, kweight=False,outdir_for_negstats = outdir)
                get_Fextr_stats(occ, Fextr.fextr_ms, Fextr.maptype, FoFo.fdif_c_ms, FoFo_type, outdir)
                compute_f_sigf(Fextr.fextr_ms, '%s' %(Fextr.maptype), log=log)
            elif mp == 'Fgenick_map':
                Fextr.fgenick(qweight=False, kweight=False,outdir_for_negstats = outdir)
                get_Fextr_stats(occ, Fextr.fgenick_ms, Fextr.maptype, FoFo.fdif_c_ms, FoFo_type, outdir)
                compute_f_sigf(Fextr.fgenick_ms, '%s' %(Fextr.maptype), log=log)
            elif mp == 'Fextr_calc_map':
                Fextr.fextr_calc(qweight=False, kweight=False,outdir_for_negstats = outdir)
                get_Fextr_stats(occ, Fextr.fextr_calc_ms, Fextr.maptype, FoFo.fdif_c_ms, FoFo_type, outdir)
                compute_f_sigf(Fextr.fextr_calc_ms, '%s' %(Fextr.maptype), log=log)
            else:
                print("%s not recognised as extrapolated map type" % mp)
                        
            #Use ccp4 map of type mFo-DFc to integrate the masked map
            print("\n************Map explorer************", file=log)
            print("\n************Map explorer************")
            #map_expl_out = map_explorer(Fextr.ccp4_name_FoFc, DH.pdb_in, params.map_explorer.radius, params.map_explorer.peak_integration_floor, params.map_explorer.peak_detection_threshold, maptype=Fextr.maptype)
            data = ccp4_map.map_reader(file_name=Fextr.ccp4_name_FoFc).data.as_numpy_array()
            pos = 0
            neg = 0
            #print(mask[0,0])
            
            for i in range(mask.shape[1]):
                tmp = data[mask[0, i]].sum()
                if tmp > 0: pos+= tmp
                else: neg -= tmp
                
            #print("data",data)
            #print("data.shape", data.shape)
            #integrated_values.append([pos, neg, pos+neg])
            try:
                CC = pearsonr(fofo_data.flatten(), data.flatten())[0]
            except ValueError:
                #the fft_map function might change the cystal gridding. In that case the maps do not have the same shape
                #and thus the CC cannot be calculated. Not nice, but at least Xtrapol8 can continue.
                #in this case, also the output from plotalpha is wrong since the mask will be incorrectly projected!!!
                print("Pearson correlation factor could not be calculated. The CC will be set to zero.")
                CC = 0
            #pearsonCC.append(scipy.stats.pearsonr(fofo_data.flatten(), data.flatten())[0])

            #depending on the map-type, append the output-file of mapexplorer to the correct list
            if mp == 'qFextr_map':
                #append_if_file_exist(qFextr_map_expl_fles, os.path.abspath(map_expl_out))
                qFextr_map_expl_fles.append([CC, pos, neg, pos+neg])
            elif mp == 'qFgenick_map':
                #append_if_file_exist(qFgenick_map_expl_fles, os.path.abspath(map_expl_out))
                qFgenick_map_expl_fles.append([CC, pos, neg, pos+neg])
            elif mp == 'qFextr_calc_map':
                #append_if_file_exist(qFextr_calc_map_expl_fles, os.path.abspath(map_expl_out))
                qFextr_calc_map_expl_fles.append([CC, pos, neg, pos+neg])
            elif mp == 'kFextr_map':
                #append_if_file_exist(kFextr_map_expl_fles, os.path.abspath(map_expl_out))
                kFextr_map_expl_fles.append([CC, pos, neg, pos + neg])
            elif mp == 'kFgenick_map':
                #append_if_file_exist(kFgenick_map_expl_fles, os.path.abspath(map_expl_out))
                kFgenick_map_expl_fles.append([CC, pos, neg, pos + neg])
            elif mp == 'kFextr_calc_map':
                #append_if_file_exist(kFextr_calc_map_expl_fles, os.path.abspath(map_expl_out))
                kFextr_calc_map_expl_fles.append([CC, pos, neg, pos + neg])
            elif mp == 'Fextr_map':
                #append_if_file_exist(Fextr_map_expl_fles, os.path.abspath(map_expl_out))
                Fextr_map_expl_fles.append([CC, pos, neg, pos + neg])
            elif mp == 'Fgenick_map':
                #append_if_file_exist(Fgenick_map_expl_fles, os.path.abspath(map_expl_out))
                Fgenick_map_expl_fles.append([CC, pos, neg, pos + neg])
            elif mp == 'Fextr_calc_map':
                #append_if_file_exist(Fextr_calc_map_expl_fles, os.path.abspath(map_expl_out))
                Fextr_calc_map_expl_fles.append([CC, pos, neg, pos + neg])
            print("m%s-DFcalc map explored" % Fextr.maptype,file=log)
            print("m%s-DFcalc map explored" % Fextr.maptype)

            #In case of running in slow_and_rigorous / slow_and_curious (whatever you like to call the full way mode):
            #run refinement with phenix or refmac/coot
            if (params.f_and_maps.fast_and_furious == False and params.refinement.run_refinement):
                print("\n************Refinements************")
                print("\n************Refinements************", file=log)
                if params.refinement.use_refmac_instead_of_phenix:
                    mtz_out, pdb_rec, pdb_real, pdb_rec_real = Fextr.refmac_coot_refinements(pdb_in = DH.pdb_in,
                                                               additional        = DH.additional,
                                                               ligands_list      = DH.extract_ligand_codes(),
                                                               F_column_labels   = Fextr.FM.labels['data'],
                                                               map_column_labels = '%s, PHI%s, %s, PHI%s'
                                                               %(Fextr.FM.labels['map_coefs_map'],Fextr.FM.labels['map_coefs_map'],
                                                                 Fextr.FM.labels['map_coefs_diff'], Fextr.FM.labels['map_coefs_diff']),
                                                               keywords          = params.refinement.refmac_keywords)
                else:
                    mtz_out, pdb_rec, pdb_real, pdb_rec_real = Fextr.phenix_phenix_refinements(pdb_in = DH.pdb_in,
                                                    additional      = DH.additional,
                                                    F_column_labels = Fextr.FM.labels['data'],
                                                    column_labels   = '%s,PHI%s'%(Fextr.FM.labels['map_coefs_map'],Fextr.FM.labels['map_coefs_map']),
                                                    scattering_table= params.scattering_table,
                                                    keywords        = params.refinement.phenix_keywords)
                print("--------------", file=log)
                #depending on the map-type, append the refinement output to the correct list
                #this is ugly, TODO: make an object to store the results in a clean and transparant way
                if mp == 'qFextr_map':
                    append_if_file_exist(qFextr_recref_mtz_lst, os.path.abspath(mtz_out))
                    append_if_file_exist(qFextr_recref_pdb_lst, os.path.abspath(pdb_rec))
                    append_if_file_exist(qFextr_recrealref_lst, os.path.abspath(pdb_rec_real))
                    append_if_file_exist(qFextr_realref_lst, os.path.abspath(pdb_real))
                elif mp == 'qFgenick_map':
                    append_if_file_exist(qFgenick_recref_mtz_lst, os.path.abspath(mtz_out))
                    append_if_file_exist(qFgenick_recref_pdb_lst, os.path.abspath(pdb_rec))
                    append_if_file_exist(qFgenick_recrealref_lst, os.path.abspath(pdb_rec_real))
                    append_if_file_exist(qFgenick_realref_lst, os.path.abspath(pdb_real))
                elif mp == 'qFextr_calc_map':
                    append_if_file_exist(qFextr_calc_recref_mtz_lst, os.path.abspath(mtz_out))
                    append_if_file_exist(qFextr_calc_recref_pdb_lst, os.path.abspath(pdb_rec))
                    append_if_file_exist(qFextr_calc_recrealref_lst, os.path.abspath(pdb_rec_real))
                    append_if_file_exist(qFextr_calc_realref_lst, os.path.abspath(pdb_real))
                elif mp == 'kFextr_map':
                    append_if_file_exist(kFextr_recref_mtz_lst, os.path.abspath(mtz_out))
                    append_if_file_exist(kFextr_recref_pdb_lst, os.path.abspath(pdb_rec))
                    append_if_file_exist(kFextr_recrealref_lst, os.path.abspath(pdb_rec_real))
                    append_if_file_exist(kFextr_realref_lst, os.path.abspath(pdb_real))
                elif mp == 'kFgenick_map':
                    append_if_file_exist(kFgenick_recref_mtz_lst, os.path.abspath(mtz_out))
                    append_if_file_exist(kFgenick_recref_pdb_lst, os.path.abspath(pdb_rec))
                    append_if_file_exist(kFgenick_recrealref_lst, os.path.abspath(pdb_rec_real))
                    append_if_file_exist(kFgenick_realref_lst, os.path.abspath(pdb_real))
                elif mp == 'kFextr_calc_map':
                    append_if_file_exist(kFextr_calc_recref_mtz_lst, os.path.abspath(mtz_out))
                    append_if_file_exist(kFextr_calc_recref_pdb_lst, os.path.abspath(pdb_rec))
                    append_if_file_exist(kFextr_calc_recrealref_lst, os.path.abspath(pdb_rec_real))
                    append_if_file_exist(kFextr_calc_realref_lst, os.path.abspath(pdb_real))
                elif mp == 'Fextr_map':
                    append_if_file_exist(Fextr_recref_mtz_lst, os.path.abspath(mtz_out))
                    append_if_file_exist(Fextr_recref_pdb_lst, os.path.abspath(pdb_rec))
                    append_if_file_exist(Fextr_recrealref_lst, os.path.abspath(pdb_rec_real))
                    append_if_file_exist(Fextr_realref_lst, os.path.abspath(pdb_real))
                elif mp == 'Fgenick_map':
                    append_if_file_exist(Fgenick_recref_mtz_lst, os.path.abspath(mtz_out))
                    append_if_file_exist(Fgenick_recref_pdb_lst, os.path.abspath(pdb_rec))
                    append_if_file_exist(Fgenick_recrealref_lst, os.path.abspath(pdb_rec_real))
                    append_if_file_exist(Fgenick_realref_lst, os.path.abspath(pdb_real))
                elif mp == 'Fextr_calc_map':
                    append_if_file_exist(Fextr_calc_recref_mtz_lst, os.path.abspath(mtz_out))
                    append_if_file_exist(Fextr_calc_recref_pdb_lst, os.path.abspath(pdb_rec))
                    append_if_file_exist(Fextr_calc_recrealref_lst, os.path.abspath(pdb_rec_real))
                    append_if_file_exist(Fextr_calc_realref_lst, os.path.abspath(pdb_real))
                    
            #fast_and_furious mode: no refinement, but need to keep track of structure factor files
            elif params.f_and_maps.fast_and_furious:
                append_if_file_exist(Fextr_mtz_lst, os.path.abspath(Fextr.F_name))

            print("\n---> Results in %s" %(os.getcwd()), file=log)
            print("------------------------------------", file=log)
            print("\n---> Results in %s" %(os.getcwd()))
            print("------------------------------------")
            
        ################################################################
        #remove empty directories to avoid any confusion. secure because os.rmdir can only remove empty directories.
        if len(os.listdir(new_dirpath_q)) == 0:
            os.rmdir(new_dirpath_q)
        if len(os.listdir(new_dirpath)) == 0:
            os.rmdir(new_dirpath)
        if len(os.listdir(new_dirpath_k)) == 0:
            os.rmdir(new_dirpath_k)

        #Go back to output directory and generate the plots from the pickle files
        os.chdir(outdir)
        plot_Fextr_sigmas()
        #plot_sigmas(maptype_lst=list(map(lambda x: re.sub(r"\_map$", "", x), final_maptypes)))
        plot_negative_reflections()

    #free some memory by deleting the fofo_map
    del fofo_data
    ################################################################
    #Repeat the last step in order to make sure we have all plots correctly
    #Go back to output directory and generate the plots from the pickle files
    os.chdir(outdir)
    plot_Fextr_sigmas()
    #plot_sigmas(maptype_lst=list(map(lambda x: re.sub(r"\_map$", "", x), final_maptypes)))
    plot_negative_reflections()
    
    #plot_correlations(params.occupancies.list_occ, cc_list)
    
    print('-----------------------------------------')
    print("CALCULATE (Q/K-WEIGHTED) FEXTRAPOLATED MAPS DONE")
    print('-----------------------------------------')
    
    ################################################################
    print("ESTIMATE OPTIMAL OCCUPANCY")
    print('-----------------------------------------')
    
    print('-----------------------------------------', file=log)
    print("ESTIMATE OPTIMAL OCCUPANCY", file=log)
    print('-----------------------------------------', file=log)

    #Use Z-score list to select the residues. The standalone version can be used if a specific residue list needs to be used.
    residlst = residlist_zscore
    print("Residue list used for estimation of occupancy of triggered state: %s\n" %(residlst), file=log)
    print("Residue list used for estimation of occupancy of triggered state: %s\n" %(residlst))
    
    ################################################################
    #occupancy estimation will be performed for all map types, but the final automatic descision will be based on a priority list.
    #This means that if maptypes Fextr and qFextr are selected, the output of qFextr will have priority.
    #Therefore, here we loop over the maptypes in a reversed order to
    #1) transfer the stored output files to a generic list
    #2) estimate occupacy based on the peakintegration volumes
    #If slow_and_rigorous:
    #   3) plot the refinement R-factors
    #   4) estimate the occupancy based on the distances of the real space refined models with the input model
    #   5) append the refinement models and maps to the Pymol_movie script
    #   6) make ddm plot
    #   7) write coot script
    ddm_out = None
    script_coot = None
    reversed_maptypes = final_maptypes[:]
    reversed_maptypes.reverse()
    occ_overview = {}
    for mp in reversed_maptypes:
        mp_type = mp.split("_map")[0]
        print("---%s---\n"%(mp_type.upper()))
        print("---%s---\n"%(mp_type.upper()), file= log)
        
        if mp in ('qFextr_map','qFgenick_map','qFextr_calc_map'):
            dir_prefix = 'qweight_occupancy'
        elif mp in ('kFextr_map','kFgenick_map','kFextr_calc_map'):
            dir_prefix = 'kweight_occupancy'
        else:
            dir_prefix = 'occupancy' 
            
        if mp == 'qFextr_map':
            recref_mtz_lst = qFextr_recref_mtz_lst
            recref_pdb_lst = qFextr_recref_pdb_lst
            recrealref_lst = qFextr_recrealref_lst
            realref_lst    = qFextr_realref_lst
            map_expl_lst   = qFextr_map_expl_fles
        elif mp == 'qFgenick_map':
            recref_mtz_lst = qFgenick_recref_mtz_lst
            recref_pdb_lst = qFgenick_recref_pdb_lst
            recrealref_lst = qFgenick_recrealref_lst
            realref_lst    = qFgenick_realref_lst   
            map_expl_lst   = qFgenick_map_expl_fles 
        elif mp == 'qFextr_calc_map':
            recref_mtz_lst = qFextr_calc_recref_mtz_lst
            recref_pdb_lst = qFextr_calc_recref_pdb_lst
            recrealref_lst = qFextr_calc_recrealref_lst
            realref_lst    = qFextr_calc_realref_lst   
            map_expl_lst   = qFextr_calc_map_expl_fles
        elif mp == 'kFextr_map':
            recref_mtz_lst = kFextr_recref_mtz_lst
            recref_pdb_lst = kFextr_recref_pdb_lst
            recrealref_lst = kFextr_recrealref_lst
            realref_lst    = kFextr_realref_lst
            map_expl_lst   = kFextr_map_expl_fles
        elif mp == 'kFgenick_map':
            recref_mtz_lst = kFgenick_recref_mtz_lst
            recref_pdb_lst = kFgenick_recref_pdb_lst
            recrealref_lst = kFgenick_recrealref_lst
            realref_lst    = kFgenick_realref_lst   
            map_expl_lst   = kFgenick_map_expl_fles 
        elif mp == 'kFextr_calc_map':
            recref_mtz_lst = kFextr_calc_recref_mtz_lst
            recref_pdb_lst = kFextr_calc_recref_pdb_lst
            recrealref_lst = kFextr_calc_recrealref_lst
            realref_lst    = kFextr_calc_realref_lst   
            map_expl_lst   = kFextr_calc_map_expl_fles 
        elif mp == 'Fextr_map':
            recref_mtz_lst = Fextr_recref_mtz_lst
            recref_pdb_lst = Fextr_recref_pdb_lst
            recrealref_lst = Fextr_recrealref_lst
            realref_lst    = Fextr_realref_lst   
            map_expl_lst   = Fextr_map_expl_fles 
        elif mp == 'Fgenick_map':
            recref_mtz_lst = Fgenick_recref_mtz_lst
            recref_pdb_lst = Fgenick_recref_pdb_lst
            recrealref_lst = Fgenick_recrealref_lst
            realref_lst    = Fgenick_realref_lst   
            map_expl_lst   = Fgenick_map_expl_fles 
        elif mp == 'Fextr_calc_map':
            recref_mtz_lst = Fextr_calc_recref_mtz_lst
            recref_pdb_lst = Fextr_calc_recref_pdb_lst
            recrealref_lst = Fextr_calc_recrealref_lst
            realref_lst    = Fextr_calc_realref_lst   
            map_expl_lst   = Fextr_calc_map_expl_fles 

        if params.map_explorer.occupancy_estimation in ("difference_map_maximization", "distance_analysis"):
            #in case of distance_analysis alpha and occ will be overwritten if the requirements for distance_analysis are met (calm-and-curious, run_refinement)
            alpha, occ, _, _ = plotalpha(params.occupancies.list_occ, map_expl_lst[1:], map_expl_lst[0], mp_type, log=log).estimate_alpha()
        elif params.map_explorer.occupancy_estimation == "difference_map_PearsonCC":
            _, _, alpha, occ = plotalpha(params.occupancies.list_occ, map_expl_lst[1:], map_expl_lst[0], mp_type, log=log).estimate_alpha()
        
        if (params.f_and_maps.fast_and_furious == False and params.refinement.run_refinement):
            # If water molecules are updated during refinement, waters will be added and removed and their numbers are not in relation to the original waters in the input model
            # therefore keeping them during the distance analysis would make no sense
            if params.refinement.phenix_keywords.main.ordered_solvent:
                distance_use_waters = False
            else:
                distance_use_waters = True
                
            #Make a plot of the refinement R-factors, related to the specific maptype. The log-files should have the same prefix as the mtz-files.
            #This assumption is made in order to avoid storing the log-files in even another list
            plot_Rfactors_per_alpha(map(lambda fle: re.sub(r'mtz$','log', fle), recref_mtz_lst), mp_type)
            print("", file=log)
            print("")
            
            #Check if the list with PDB files from real space refinement after reciprocal space refinement is complete
            #   If this is not the case (when reciprocal space or real space refinement had failed), 
            #   the PDB files from real space refinment (direct real space refinement in the extrapolated maps) are used
            if ((len(recrealref_lst)==len(set(recrealref_lst)))):
                pdb_list = recrealref_lst
            else:
                pdb_list = realref_lst
                
            #Estimate alpha and occupancy based on the distances of the list with PDB files with the input model.
            #Overwrite the earlier found alpha and occupancy if params.map_explorer.occupancy_estimation = "distance_analysis"
            #if params.map_explorer.use_occupancy_from_distance_analysis:
            if params.map_explorer.occupancy_estimation == "distance_analysis":
                alpha, occ = Distance_analysis(pdb_list, params.occupancies.list_occ, resids_lst= residlst, use_waters = distance_use_waters, outsuffix = mp_type, log = log).extract_alpha()
            else:
                _,_ = Distance_analysis(pdb_list, params.occupancies.list_occ, resids_lst= residlst, use_waters = distance_use_waters, outsuffix = mp_type, log = log).extract_alpha()
            #print("---------", file=log)
            print("", file=log)
            print("")            
            
            #add models and maps to pymol movie.
            #Not elegant with all the with hard coded regular expressions but avoids storing other lists
            pymol_pdb_list = pdb_list[:]
            pymol_pdb_list.remove(DH.pdb_in)
            pymol_mtz_list = recref_mtz_lst[:]
            if outname == 'triggered': #if dummy name applied, the files still contain the dummy name
                pymol_mtz_list = map(lambda fle: re.sub(r"triggered",params.output.outname, fle), pymol_mtz_list)
                pymol_pdb_list = map(lambda fle: re.sub(r"triggered",params.output.outname, fle), pymol_pdb_list)
            #Make Pymol movie with the reciprocal space refined maps if recrealref_lst is complete
            #Otherwise use the real space refined models + direct maps
            if pdb_list == recrealref_lst:
                ccp4_list = map(lambda fle: re.sub(r".mtz$", "_2mFo-DFc_filled.ccp4", fle), pymol_mtz_list)
                model_label = '%s_reciprocal_real_space'%(mp_type)
                ccp4_map_label = '%s_reciprocal_space'%(mp)
                #Pymol_movie(params.occupancies.list_occ, pdblst=pymol_pdb_list, ccp4_maps = ccp4_list, resids_lst = residlst, model_label='%s_reciprocal_real_space'%(mp_type), ccp4_map_label='%s_reciprocal_space'%(mp)).write_pymol_script()
            else: 
                if mp == 'qFgenick_map':
                    ccp4_list = map(lambda fle: re.search("(.+?)2mqFgenick-DFc_reciprocal", fle).group(1)+"mqFgenick-DFc.ccp4", pymol_mtz_list)
                elif mp == 'kFgenick_map':
                    ccp4_list = map(lambda fle: re.search("(.+?)2mkFgenick-DFc_reciprocal", fle).group(1)+"mkFgenick-DFc.ccp4", pymol_mtz_list)
                elif mp == 'Fgenick_map':
                    ccp4_list = map(lambda fle: re.search("(.+?)2mFgenick-DFc_reciprocal", fle).group(1)+"mFgenick-DFc.ccp4", pymol_mtz_list)
                else:
                    ccp4_list = map(lambda fle: re.search("(.+?)\_reciprocal", fle).group(1)+".ccp4", pymol_mtz_list)
                model_label='%s_real_space'%(mp_type)
                ccp4_map_label='%s'%(mp)
            if len(ccp4_list) == len(pymol_pdb_list) == len(params.occupancies.list_occ):
                Pymol_movie(params.occupancies.list_occ, pdblst=pymol_pdb_list, ccp4_maps = ccp4_list, resids_lst = residlst, model_label=model_label, ccp4_map_label=ccp4_map_label).write_pymol_script()
            else:
                Pymol_movie(params.occupancies.list_occ, pdblst=pymol_pdb_list, resids_lst = residlst, model_label=model_label).write_pymol_script()
                
            #If estimated occupancy if not in list (will be case when using distance analysis or when plotalpha fails), take the closest occupancy from the list
            if occ not in params.occupancies.list_occ:
                occ = min(params.occupancies.list_occ, key = lambda x: abs(x-occ))
                alpha = 1/occ
            occ_dir = "%s/%s_%.3f" %(outdir, dir_prefix, occ)
            
            #ddm calculation
            pdb_for_ddm = pdb_list[params.occupancies.list_occ.index(occ)+1]
            print("----Generate distance difference plot----")
            ddm_out = Difference_distance_analysis(DH.pdb_in, pdb_for_ddm, ligands = DH.extract_ligand_codes(), outdir=occ_dir, scale=params.output.ddm_scale).ddms()
            print('---------------------------')
            
            #Coot script
            mtzs_for_coot  = []
            if len(recref_mtz_lst) != len(params.occupancies.list_occ):
                occ_list_mtz = [(re.search(r'occupancy\_(.+?)\/',mtz).group(1)) for mtz in recref_mtz_lst]
                try:
                    mtz_rec = recref_mtz_lst[occ_list_mtz.index(occ)]
                except ValueError:
                    mtz_rec = ""
            else:
                mtz_rec = recref_mtz_lst[params.occupancies.list_occ.index(occ)]
            append_if_file_exist(mtzs_for_coot, mtz_rec)
            if ( params.refinement.phenix_keywords.density_modification.density_modification or params.refinement.refmac_keywords.density_modification.density_modification):
                #mtz_dm = re.sub(".mtz$","_densitymod.mtz", mtz_rec)
                mtz_dm = re.sub(".mtz$","_dm.mtz", mtz_rec)
                append_if_file_exist(mtzs_for_coot, mtz_dm)

            #if params.refinement.refmac_keywords.density_modification.density_modification:
                #mtz_dm = re.sub(".mtz$","_dm.mtz", mtz_rec)
                #append_if_file_exist(mtzs_for_coot, mtz_dm)
            mtz_extr = ["%s/%s"%(occ_dir,fle) for fle in os.listdir(occ_dir) if outname in fle and fle.endswith('m%s-DFc.mtz'%(mp_type))][0]
            append_if_file_exist(mtzs_for_coot,os.path.abspath(mtz_extr))

            pdbs_for_coot = [DH.pdb_in,
                    recref_pdb_lst[params.occupancies.list_occ.index(occ)+1],
                    realref_lst[params.occupancies.list_occ.index(occ)+1],
                    recrealref_lst[params.occupancies.list_occ.index(occ)+1]]
            #if outname == 'triggered': #if dummy name applied, the files still contain the dummy name
                #mtzs_for_coot = map(lambda fle: re.sub(r"triggered",params.output.outname, fle), mtzs_for_coot)
                #pdbs_for_coot = map(lambda fle: re.sub(r"triggered",params.output.outname, fle), pdbs_for_coot)
            script_coot = open_all_in_coot(outdir+"/"+FoFo.mtz_name, pdbs_for_coot, mtzs_for_coot, DH.additional, occ_dir, mp_type)
        
        elif (params.f_and_maps.fast_and_furious == False and params.refinement.run_refinement == False):
            #If estimated occupancy if not in list (will be case when using distance analysis or when plotalpha fails), take the closest occupancy from the list
            if occ not in params.occupancies.list_occ:
                occ = min(params.occupancies.list_occ, key = lambda x: abs(x-occ))
                alpha = 1/occ
            occ_dir = "%s/%s_%.3f" %(outdir, dir_prefix, occ)
            mtzs_for_coot = []
            mtz_extr = ["%s/%s"%(occ_dir,fle) for fle in os.listdir(occ_dir) if outname in fle and fle.endswith('m%s-DFc.mtz'%(mp_type))][0]
            append_if_file_exist(mtzs_for_coot,os.path.abspath(mtz_extr))
            pdbs_for_coot = [DH.pdb_in]
            script_coot = open_all_in_coot(outdir+"/"+FoFo.mtz_name, pdbs_for_coot, mtzs_for_coot, DH.additional, occ_dir, mp_type)
         
        else:
            #If estimated occupancy not in list (will be case when using distance analysis or when plotalpha fails), take the closest occupancy from the list
            if occ not in params.occupancies.list_occ:
                occ = min(params.occupancies.list_occ, key = lambda x: abs(x-occ))
                alpha = 1/occ
            occ_dir = "%s/%s_%.3f" %(outdir, dir_prefix, occ)
            
        occ_overview[mp_type] = [float("%.3f"%(occ)), script_coot, ddm_out]
            
        print("------------------------------------")
        print("------------------------------------", file=log)
    
    #Add final lines to Pymol_script
    if os.path.isfile('%s/pymol_movie.py' %(outdir)):
        Pymol_movie(params.occupancies.list_occ, resids_lst = residlst).write_pymol_appearance('%s/pymol_movie.py' %(outdir))
        
    #Send the dictonary to the GUI in order to recuperate the occupancies found for each maptype -> write as pickle file to be opened by the GUI
    #if params.output.GUI:
        #pub.sendMessage("Best occ", occ_overview=occ_overview)
        #print("Message sent to GUI")
    #Write the pickle file always as this will also be used when the GUI is launched as resultsloader
    occ_pickle = open("occupancy_recap.pickle", "wb")
    pickle.dump(occ_overview, occ_pickle)
    occ_pickle.close()
        
      
    print("Summary of occupancy estimation:", file=log)
    print("Method:  {:s}".format(params.map_explorer.occupancy_estimation),file=log)
    print("Map type       Occupancy", file=log)
    for k in occ_overview:
        print("{:<15} {:>5.3f}".format(k, occ_overview[k][0]), file=log)
    print("-> Optimal occupancy of triggered state %.3f." %(occ), file=log)
    
    print("OCCUPANCY ESTIMATION SUMMARY")
    print("Method:  {:s}".format(params.map_explorer.occupancy_estimation))
    print("Map type       Occupancy")
    for k in occ_overview:
        print("{:<15} {:>5.3f}".format(k, occ_overview[k][0]))
    print("-> Optimal occupancy of triggered state %.3f." %(occ))

    
    print('-----------------------------------------')
    print("ESTIMATE OPTIMAL OCCUPANCY DONE")
    print('-----------------------------------------')
    
    ################################################################
    #Run refinements for chosen occupancy when in faf mode (as no refinement has yet run):
    #1) move to directory where refinements have to be run
    #2) find correct input files
    #3) Run refinements
    #4) make ddm plot
    #5) write coot script
    if (params.f_and_maps.fast_and_furious and params.refinement.run_refinement):
        print('-----------------------------------------', file=log)
        print("FAST AND FURIOUS REFINEMENT", file=log)
        print('-----------------------------------------', file=log)
        
        print('-----------------------------------------')
        print("FAST AND FURIOUS REFINEMENT")
        print('-----------------------------------------')


        print("----Refinements with automatically found occupancy values----", file=log)
        print("----Refinements with automatically found occupancy values----")
        os.chdir(occ_dir)

        qFext_mtz_F = [fle for fle in Fextr_mtz_lst if 'occupancy_%.3f' %(occ) in fle and fle.lower().endswith('_qfextr.mtz')][0]
        assert os.path.isfile(qFext_mtz_F)
        print("Extrapolated structure factors for reciprocal space refinement :%s" %(qFext_mtz_F), file=log)
        print("Extrapolated structure factors for reciprocal space refinement :%s" %(qFext_mtz_F))
        
        qFext_mtz_map = re.sub(r'qFextr.mtz$','2mqFextr-DFc_mqFextr-DFc.mtz', qFext_mtz_F)
        assert os.path.isfile(qFext_mtz_map)
        print("Extrapolated map coefficients for real space refinement: %s" %(qFext_mtz_map), file=log)
        print("Extrapolated map coefficients for real space refinement: %s\n" %(qFext_mtz_map))
        
        Fextr.name_out = "%s_occ%.3f" %(outname, occ)
        
        if params.refinement.use_refmac_instead_of_phenix:
            mtz_out, pdb_rec, pdb_real, pdb_rec_real = Fextr.refmac_coot_refinements(mtz_F = qFext_mtz_F,
                                                       mtz_map           = qFext_mtz_map,
                                                       pdb_in            = DH.pdb_in,
                                                       additional        = DH.additional,
                                                       ligands_list      = DH.extract_ligand_codes(),
                                                       F_column_labels   = Fextr.FM.labels['data'],
                                                       map_column_labels = '%s, PHI%s, %s, PHI%s'
                                                               %(Fextr.FM.labels['map_coefs_map'],Fextr.FM.labels['map_coefs_map'],
                                                                 Fextr.FM.labels['map_coefs_diff'], Fextr.FM.labels['map_coefs_diff']),
                                                       keywords          = params.refinement.refmac_keywords)
        else:
            mtz_out, pdb_rec, pdb_real, pdb_rec_real = Fextr.phenix_phenix_refinements(mtz_F = qFext_mtz_F,
                                                  mtz_map         = qFext_mtz_map,
                                                  pdb_in          = DH.pdb_in,
                                                  additional      = DH.additional,
                                                  F_column_labels = Fextr.FM.labels['data'],
                                                  column_labels   = '%s,PHI%s' %(Fextr.FM.labels['map_coefs_map'],Fextr.FM.labels['map_coefs_map']),
                                                  scattering_table= params.scattering_table,
                                                  keywords        = params.refinement.phenix_keywords)
              
        print("---> Results in %s_%.3f"%(dir_prefix, occ), file=log)
        print("---> Results in %s_%.3f"%(dir_prefix, occ))
        
        
        #Coot script
        mtzs_for_coot = [os.path.abspath(qFext_mtz_map), os.path.abspath(mtz_out)]
        if params.refinement.phenix_keywords.density_modification.density_modification:
            mtz_dm = re.sub(".mtz$","_densitymod.mtz", mtz_out)
            if os.path.isfile(mtz_dm):
                mtzs_for_coot.append(os.path.abspath(mtz_dm))
        if params.refinement.phenix_keywords.density_modification.density_modification:
            mtz_dm = re.sub(".mtz$","_dm.mtz", mtz_out)
            if os.path.isfile(mtz_dm):
                mtzs_for_coot.append(os.path.abspath(mtz_dm))
        append_if_file_exist(mtzs_for_coot, os.path.abspath(qFext_mtz_F))
        pdbs_for_coot = [DH.pdb_in, os.path.abspath(check_file_existance(pdb_rec)), os.path.abspath(check_file_existance(pdb_real)), os.path.abspath(check_file_existance(pdb_rec_real))]
        #if outname == 'triggered': #if dummy name applied, the files still contain the dummy name
            #mtzs_for_coot = map(lambda fle: re.sub(r"triggered",params.output.outname, fle), mtzs_for_coot)
            #pdbs_for_coot = map(lambda fle: re.sub(r"triggered",params.output.outname, fle), pdbs_for_coot)
        script_coot = open_all_in_coot(outdir+"/"+FoFo.mtz_name, pdbs_for_coot, mtzs_for_coot, DH.additional, occ_dir, "qFextr")

        print('-----------------------------------------')
        print("FAST AND FURIOUS REFINEMENT DONE")
        print('-----------------------------------------')
        print("----Generate distance difference plot----")
        ddm_out = Difference_distance_analysis(DH.pdb_in, pdb_rec_real, ligands = DH.extract_ligand_codes(), outdir=occ_dir, scale=params.output.ddm_scale).ddms()
        print('---------------------------')
        
        #Rewrite the pickle file as to include the ddm_out path
        os.chdir(outdir)
        occ_overview[mp_type] = [float("%.3f"%(occ)), script_coot, ddm_out]
        #if os.path.isfile("occupancy_recap.pickle"):
            #os.remove("occupancy_recap.pickle")
        occ_pickle = open("occupancy_recap.pickle", "wb")
        pickle.dump(occ_overview, occ_pickle)
        occ_pickle.close()

        
    ################################################################
    #Make sure we are in the output directory
    if os.getcwd() != outdir:
        os.chdir(outdir)
    
    #change names to real output name in case the dummy name was used
    if outname == 'triggered':
        print('---------------------------', file=log)
        #print("Replacing the dummpy outname ('%s') with true outname('%s')" %(outname, params.output.outname), file=log)
        print("Replacing the dummpy outname ('%s') with true outname('%s')" %(outname, params.output.outname))
        tr = [os.path.join(root, fle) for root, dirs, files in os.walk(outdir) for fle in files if outname in fle]
        _ = [os.rename(fle, fle.replace(outname, params.output.outname)) for fle in tr]
        FoFo.mtz_name = FoFo.mtz_name.replace(outname, params.output.outname)
        
        coot_scripts = [os.path.join(root, fle) for root, dirs, files in os.walk(outdir) for fle in files if fle.startswith("coot_all_")]
        for coot_script in coot_scripts:
            with open(coot_script, "r") as f:
                fle = f.read().split("\n")
            o = open(coot_script,"w")
            for lne in fle:
                if outname in lne:
                    lne = re.sub(outname, params.output.outname, lne)
                    o.write("%s\n"%lne)
                else:
                    o.write("%s\n"%lne)
            o.close()

    #Make sure we are in the output directory
    if os.getcwd() != outdir:
        os.chdir(outdir)
    #Make a list of all temporal files. Upon their removal disk space can be saved
    list_redundant_files(outdir)
    
    #################################################################
    #print('-----------------------------------------', file=log)
    #print("WRITE PYMOL AND COOT SCRIPTS", file=log)
    #print('-----------------------------------------', file=log)
    
    ##The pymol movie script is already made
    #print("pymol movie with models and maps in %s/pymol_movie.py. Run from terminal or pymol command line"%(outdir), file=log)
    
    ##The pymol script with just all models and maps is still missing
    ##This script is not working properly and models and maps are in wrong session subfolers or even missing
    #pymol_script = Pymol_visualization(DH.pdb_in, outdir).open_all_in_pymol()
    #print("pymol script with all models written in %s/%s. Can be run from a terminal with 'Pymol %s/%s' if Pymol is in your PATH or 'run %s/%s'from a pymol command line"%(outdir, pymol_script, outdir, pymol_script, outdir, pymol_script), file=log)
        
    #print('-----------------------------------------', file=log)
    #print("WRITE PYMOL AND COOT SCRIPTS DONE", file=log)
    #print('-----------------------------------------', file=log)

    ################################################################
    #Write a phil file. Parameters changed during the excecution of Xtrapol8 are possible when problems appeared
    if params.f_and_maps.negative_and_missing == "fill_missing":
        params.f_and_maps.negative_and_missing = "keep_and_fill"
    if params.f_and_maps.negative_and_missing == "no_fill":
        params.f_and_maps.negative_and_missing = "keep_no_fill"
    
    modified_phil = master_phil.format(python_object=params)
    modified_phil.show(out=open("Xtrapol8_out.phil", "w"))
    
    print('-----------------------------------------')
    print("XTRAPOL8 DONE!")
    print('-----------------------------------------')
    print("Please inspect the output and rerun with different parameters if necessary")
    
    print('-----------------------------------------', file=log)
    print("XTRAPOL8 DONE!", file=log)
    print('-----------------------------------------',file=log)
    print("Please inspect the output and rerun with different parameters if necessary", file=log)
    
    if (params.f_and_maps.fast_and_furious == False and params.refinement.run_refinement):
        print("Pymol script with models and maps: %s/pymol_movie.py."%(outdir), file=log)
        print("Pymol script with models and maps: %s/pymol_movie.py."%(outdir))
    print("Coot script with models and maps in each output directory associated with estimated occupancy (coot_all_<maptype>.py).", file=log)
    print("Coot script with models and maps in each output directory associated with estimated occupancy (coot_all_<maptype>.py).")

    log.close()
    
    if (params.output.open_coot and check_program_path('coot')[1]):
        os.system("coot --script %s" %(script_coot))

    ################################################################
if __name__ == "__main__":
    run(sys.argv[1:])
    print('Hoera')
