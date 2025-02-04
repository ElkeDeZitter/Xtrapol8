# -*- coding: utf-8 -*-
"""
Script to perform the last steps of Xtrapol8 independently:
- refinement (using phenix.refine and phenix.real_space_refine -- Refmac and Coot are not implemented yet)
- occupapancy determination using the distance_analysis method

This can be very usefull if refinement with the input model as such does makes sense, e.g.
- if large conformation changes take place that cannot be fitted by the automatic real space refinement
- if bond are broken or formed
- if addiational molecules needs to be added, e.g. in case of ligand binding
- if you want to add specific refinement parameters that are not covered in Xtrapol8 (yet -> please contact us if you want new parameters being added in Xtrapol8)

This script requires that Xtrapol8 has been completely run and files not renamed.

What do you need?
- The Xtrapol8 output phil file which can be found in the Xtrapol8 output folder and is called Xtrapol8_out.phil
- independent input files for phenix.refine and phenix.real_space_refine (the refinement parameters in Xtrapol8_out.phil will be ignored)
- a model to start the refinements with
- a list with the residues on the occupancy estimation will be based (optional)

usage
-----
Run without any argument to see all options and explanation:
>>> phenix.python <wherever>/refiner.py
Parameters can be added using an input file or via command line

example using input file (preferable)
-------
1) Change the nano refiner.phil using your favourite editor, e.g.
>>> nano refiner.phil
2) Run refiner
>>> phenix.python <wherever>/refiner.py refiner.phil

example using command line only
-------
1) Run refiner with all your arguments
>>> phenix.python <wherever>/refiner.py input.model_pdb=hoera_manually_changed.pdb input.Xtrapol8_out=fancy_party/Xtrapol8_out.phil map_explorer.residue_list=fancy_party/residlist_Zscore2.00.txt refinement.reciprocal_space_phil=reciprocal.phil refinement.real_space_phil=real.phil

example using input file and command line
-------
1) Change the nano refiner.phil using your favourite editor, e.g.
>>> nano refiner.phil
2) refiner with additional arguments. The order of arguments determines how paramters will be overwriten:
>>> phenix.python <wherever>/refiner.py refiner.phil map_explorer.do_distance_analysis=False

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


from __future__ import print_function
import os, re, sys
import glob
from datetime import datetime
import shutil
import iotbx.phil
from libtbx.utils import Usage
from iotbx.file_reader import any_file
from mmtbx.scaling.matthews import p_vm_calculator
import uuid

from Fextr_utils import *
from distance_analysis import Distance_analysis
from pymol_visualization import Pymol_movie
from ddm import Difference_distance_analysis
import version

from master import master_phil
Xtrapol8_master_phil = master_phil

master_phil = iotbx.phil.parse("""
input{
    model_pdb = None
        .type = path
        .help = Model to start all refinements from.
        .expert_level = 0
    Xtrapol8_out = None
        .type = path
        .help = Xtrapol8_out.phil which can be found in the Xtrapol8 output directory
        .expert_level = 0
    }
map_explorer{
    do_distance_analysis = True
        .type = bool
        .help = Estimate occupancy based on the differences between model_pdb and real-space refined model (only in slow_and_curious mode).
        .expert_level = 0
    residue_list = None
        .type = path
        .help = Residue list to use in the distance analysis (e.g. residue list coming from map-explorer Z-score analysis).
        .expert_level = 0
}
refinement{
    reciprocal_space_phil = None
        .type = path
        .help = input file with (non-default) parameters for reciprocal space refinement. Replaces the refinement parameters from the Xtrapol8 run. Input files (mtz, pdb and additional files) will automatically be filled but extra additional cif files can be added but should be provided with their full path. Do not specify output parameters to avoid problems in the downstream analysis.
        .expert_level = 0
    real_space_phil = None
        .type = path
        .help = input file with (non-default) parameters for real space refinements. Replaces the refinement parameters from the Xtrapol8 run. Input files (mtz, pdb and additional files) will automatically be filled but extra additional cif files can be added but should be provided with their full path. Do not specify output parameters to avoid problems in the downstream analysis.
        .expert_level = 0
}
""", process_includes=True)

class DataHandler(object):
    def __init__(self,
                 refine_pdb_in,
                 additional,
                 outdir,
                 X8_pdb_in,
                 reciprocal_space_phil = None,
                 real_space_phil= None,
                 residue_list = None):
        self.refine_pdb_in         = refine_pdb_in
        self.X8_pdb_in             = X8_pdb_in
        self.additional            = additional
        self.outdir                = outdir
        self.reciprocal_space_phil = reciprocal_space_phil
        self.real_space_phil       = real_space_phil
        self.residue_list          = residue_list

    def check_all_files_and_dirs(self):
        """
        Check if files exist and return their absolute path
        """
        err = 0
        err_m = ''
        warning = 0
        warning_m = ''
        # Check the pdb file for refinement
        if self.refine_pdb_in == None:
            err = 1
            err_m += '\nPdb file should be supplied'
        else:
            if self.check_single_file(self.refine_pdb_in):
                self.refine_pdb_in = os.path.abspath(self.refine_pdb_in)
            else:
                err = 1
                err_m += '\nFile not found: %s' %(self.refine_pdb_in)

        # Check the pdb file for distance analysis
        if self.check_single_file(self.X8_pdb_in):
            self.X8_pdb_in = os.path.abspath(self.X8_pdb_in)
        else:
            self.X8_pdb_in != None
            warning = 1
            warning_m += '\nXtrapol8 pdb_in not found. No additional analysis will be applied'

        # Check additional files and append them to a string
        additional = ""
        for fle in self.additional:
            if len(fle)>0:
                if self.check_single_file(fle):
                    new_add = os.path.abspath(fle)
                    additional = additional + "%s " % (new_add)
                else:
                    err = 1
                    err_m += '\nFile not found: %s' %(fle)
        self.additional = additional

        #Check the output directory
        if os.path.isdir(self.outdir):
            self.outdir = os.path.abspath(self.outdir)
        else:
            err = 1
            err_m += "\nXtrapol8 output directory cannot be found." \
                    "Please run this from the same directory from which you ran Xtrapol8."

        #Check the phil file for reciprocal space refinement
        if self.check_single_file(self.reciprocal_space_phil):
            self.reciprocal_space_phil = os.path.abspath(self.reciprocal_space_phil)
        else:
            self.reciprocal_space_phil = ''
            warning = 1
            warning_m += '\nPhil for reciprocal space refinement not found. Refinement will use default parameters.'


        #Check the phil file for real space refinement
        if self.check_single_file(self.real_space_phil):
            self.real_space_phil = os.path.abspath(self.real_space_phil)
        else:
            self.real_space_phil = ''
            warning = 1
            warning_m += '\nPhil for real space refinement not found. Refinement will use default parameters.'

        #Check the residue list for distance analysis
        if self.check_single_file(self.residue_list):
            self.residue_list = os.path.abspath(self.residue_list)
        else:
            self.residue_list = None
            warning = 1
            warning_m += '\nResidue list not found. Distance analysis (if required) will be performed without residue list.'

        return err, err_m, warning, warning_m

    def check_single_file(self, fle):
        if fle == None:
            return False
        else:
            return os.path.isfile(fle)


class Refiner(object):
    def __init__ (self,
                  pdb_in,
                  maptype,
                  outname,
                  occ,
                  additional = '',
                  reciprocal_phil = '',
                  real_phil = '',
                  scattering_table = "n_gaussian",
                  density_modification = {}):
        self.pdb_in               = pdb_in
        self.maptype              = maptype.split("_map")[0]
        self.outname              = outname
        self.occ                  = occ
        self.additional           = additional
        self.reciprocal_phil      = reciprocal_phil
        self.real_phil            = real_phil
        self.scattering_table     = scattering_table
        self.density_modification = density_modification
        
    def check_single_file(self, fle):
        if fle == None:
            return False
        else:
            return os.path.isfile(fle)

    def search_mtz_file_reciprocal_space_refinement(self):
        """
        search the mtz-file in order
        """
        mtz_file = "%s_occ%.3f_%s.mtz" %(self.outname, self.occ, self.maptype)
        if self.check_single_file(mtz_file):
            self.mtz_file = mtz_file
        else:
            print("File not for reciprocal space refinement not found: %s" %(mtz_file))
            print("File not for reciprocal space refinement not found: %s" %(mtz_file), file=log)
            self.mtz_file = None
            
    def search_mtz_map_real_space_refinement(self):
        mtz_map = [fle for fle in os.listdir(os.getcwd()) if fle.startswith("%s_occ%.3f" %(self.outname, self.occ)) and fle.endswith("%s-DFc.mtz" %(self.maptype))][0]
        if self.check_single_file(mtz_map):
            self.mtz_map = mtz_map
        else:
            print("File not for real space refinement not found: %s" %(mtz_map))
            print("File not for real space refinement not found: %s" %(mtz_map), file=log)
            self.mtz_map = None
            
        
    def get_mtz_map_columns(self):
        """
        Get the labels of the columns in an mtz file containing map coefficients from Xtrapol8
        """
        coefs = any_file(self.mtz_map, force_type="hkl", raise_sorry_if_errors=False)
        #the first array should contain the map coefs for the 2Fextr-DFc map
        maplabels = coefs.file_object.as_miller_arrays()[0].info().labels
        #the second array should contain the map coefs for the Fextr-DFc map
        #maplabels_diff = coefs.file_object.as_miller_arrays()[1].info().labels
        
        return ",".join(maplabels)

    def phenix_reciprocal_space_refinement(self):
        """
        Run phenix refinements
        """
        mtz_name = get_name(self.mtz_file)
        try:
            outprefix = re.sub(r"%s"%(self.maptype), "2m%s-DFc_independent_reciprocal_space"%(self.maptype), mtz_name)
            if outprefix == mtz_name:
                raise AttributeError
        except AttributeError:
            outprefix = "%s_independent_reciprocal_space"%(mtz_name)
            
        #Specify phenix version dependent parameters
        phenix_version, phenix_subversion = self.get_phenix_version()
        print("Phenix version: %s" %(phenix_version))

        if phenix_subversion <= 20:
            r_free_flag_parameters = "refinement.input.xray_data.r_free_flags.disable_suitability_test=True refinement.input.xray_data.r_free_flags.ignore_pdb_hexdigest=True refinement.input.xray_data.r_free_flags.label='FreeR_flag' refinement.input.xray_data.r_free_flags.test_flag_value=1 "
        else: #phenix version 1.21
            r_free_flag_parameters = "data_manager.fmodel.xray_data.r_free_flags.ignore_pdb_hexdigest=True data_manager.fmodel.xray_data.r_free_flags.test_flag_value=1 "
            #data_manager.fmodel.xray_data.r_free_flags.disable_suitability_test=True
            #Disable_suitability_test cannot be done. It keeps on giving an error message about the label and value. All combination have been tested, it seems that this does not work


        reciprocal = os.system("phenix.refine --overwrite %s %s %s  %s output.prefix=%s refinement.output.write_model_cif_file=False %s refinement.main.nproc=4 write_maps=true refinement.main.scattering_table=%s " %(self.reciprocal_phil, self.mtz_file, self.additional, self.pdb_in, outprefix, r_free_flag_parameters, self.scattering_table)) # wxc_scale=0.021 #target_weights.optimize_xyz_weight=True
         
        #Find output files
        if reciprocal == 0: #os.system has correctly finished, then search for the last refined structure
            try:
                mtz_fles = glob.glob("%s_???.mtz" %(outprefix))
                # [fle for fle in os.listdir(os.getcwd()) if outprefix+"_0" in fle and fle.endswith('mtz') and not
                # fle.endswith("_densitymod.mtz")]
                mtz_fles.sort()
                mtz_out = mtz_fles[-1]
                pdb_fles = glob.glob("%s_???.pdb" %(outprefix))
                # [fle for fle in os.listdir(os.getcwd()) if outprefix+"_0" in fle and fle.endswith('pdb')]
                pdb_fles.sort()
                pdb_out = pdb_fles[-1]
            except IndexError:
                mtz_out = "%s_001.mtz"%(outprefix)
                pdb_out = "%s_001.pdb"%(outprefix)
        else: #os.system has not correctly finished
            mtz_out = "not_a_file"
            pdb_out = "refinement_did_not_finish_correcty"

        return mtz_out, pdb_out
    
    def get_solvent_content(self, pdb_in):
        """
        Extract solvent content
        """
        pdb_hier = iotbx.pdb.hierarchy.input(file_name=pdb_in)
        hier = pdb_hier.hierarchy
        overall_counts = hier.overall_counts()
        
        pdb_ini = iotbx.pdb.input(pdb_in)
        xray_structure = pdb_ini.xray_structure_simple()
        
        #scattering table should not be specified here because the structure is onluy used
        #to calculate the solvent content. No usage to calculate f_model
        
        vm_calc = p_vm_calculator(xray_structure.crystal_symmetry(),
            n_residues=overall_counts.resname_classes.get("common_amino_acid", 0))
            #can be estended with n_bases=.overall_countsresname_classes.get("common_rna_dna", 0)
        return vm_calc.solc(vm=vm_calc.vm(copies=1))


    def phenix_density_modification(self, mtz_in, pdb_in):
        """
        Run phenix.density_modification
        """
        
        solc = self.get_solvent_content(pdb_in)
        outname = re.sub(r".mtz$", "_densitymod.mtz", mtz_in)
        log_file = re.sub(r".mtz$", ".log", outname)
        
        print('Running density modification, output written to %s. Please wait...'%(log_file))
        os.system("phenix.density_modification %s %s input_files.map_coeffs_file=%s solvent_content=%.3f denmod.mask_type=histograms  output_files.output_mtz=%s > %s" %(mtz_in, pdb_in, mtz_in, solc, outname, log_file))
        
        return outname
    
    def get_F_column_labels(self, mtz):
        """
        get the labels of the first array in an mtz file
        """
        
        hkl = any_file(mtz,force_type="hkl", raise_sorry_if_errors=False)
        
        try:
            label = hkl.file_object.as_miller_arrays()[0].info().labels[0]
        except IndexError:
            print("Problem with reading mtz file: %s" %(mtz))
            print("Problem with reading mtz file: %s" %(mtz), file = log)
            label = 'QFEXTR'
        
        return label
        
    
    def write_refmac_for_dm(self, pdb_in):
        """
        Write and excecute a bash script to run Refmac in order to get an mtz-file that can be used by dm.
        """

        additional_lines = ''
        for cif in self.additional.split():
            if cif.endswith(".cif"):
                additional_lines += 'LIB_IN %s ' % (cif)

        mtz_out = "%s_for_dm.mtz" % (get_name(self.mtz_file))
        pdb_out = "%s_for_dm.pdb" % (get_name(self.mtz_file))
        log_file = "%s_refmac_for_dm.log" % (get_name(self.mtz_file))
        
        F_column_labels = self.get_F_column_labels(self.mtz_file)

        script_out = 'launch_refmac_for_dm.sh'
        i = open(script_out, 'w')
        i.write("#!/bin/sh\n\
refmac5 XYZIN %s HKLIN %s XYZOUT %s HKLOUT %s %s<<eof > %s \n\
LABIN  FP=%s SIGFP=SIG%s FREE=FreeR_flag\n\
REFI TYPE REST RESI MLKF BREF ISOT METH CGMAT \n\
nfree include 1 \n\
make check NONE \n\
ncyc 0 \n\
MAPC SHAR \n\
NOHARVEST \n\
END \n\
eof" % (pdb_in, self.mtz_file, pdb_out, mtz_out, additional_lines, log_file, F_column_labels, F_column_labels))

        i.close()
        os.system("chmod +x %s" % (script_out))
        print('Running refmac with zero cycles to prepare files suitable for running dm afterwards, output written to '
              '%s. Please wait...' % (log_file))
        os.system("./%s" % (script_out))

        return mtz_out, pdb_out

    def write_density_modification_script(self, mtz_in, pdb_in, mtz_out, log_file):
        """
        Write script to perform density modification with dm. In order to have the correct columns, the mtz-file
        should original from refmac.
        """
        solc = self.get_solvent_content(pdb_in)
        
        F_column_labels = self.get_F_column_labels(self.mtz_file) #column labels are those inherited from the extrapolated structure factors

        script_out = 'launch_dm.sh'
        i = open(script_out, 'w')
        i.write('#!/bin/sh \n\
\n\
#dm:\n\
dm hklin %s hklout %s <<eor > %s \n\
SOLC %.3f\n\
MODE SOLV HIST MULTI SAYR\n\
COMBINE %s\n\
NCYC %d\n\
LABI FP=%s SIGFP=SIG%s PHIO=PHIC_ALL FOMO=FOM\n\
LABO FDM=FDM PHIDM=PHIDM\n\
eor\n' % (mtz_in, mtz_out, log_file, solc, self.density_modification.combine, self.density_modification.cycles, F_column_labels, F_column_labels))

        ccp4_map_name = re.sub(r".mtz$", ".ccp4", mtz_out)

        i.write('#generate map in ccp4 format\n\
fft hklin %s mapout %s <<eof > fft.log\n\
LABI F1=FDM PHI=PHIDM\n\
eof' % (mtz_out, ccp4_map_name))

        i.close()
        os.system("chmod +x %s" % (script_out))
        return script_out

    def ccp4_dm(self, pdb_in):
        """
        Use dm to perform density modification. In order to get the correct columns, we first run refmac with 0
        cycles with the refined model from phenix.refine.
        """
        mtz_for_dm, pdb_for_dm = self.write_refmac_for_dm(pdb_in)
        #mtz_out_dm = re.sub(r"for_dm.mtz$", "dm.mtz", mtz_for_dm)
        mtz_out_dm = re.sub(r".pdb$", "_dm.mtz", pdb_in)
        if (os.path.isfile(mtz_for_dm) and os.path.isfile(pdb_for_dm)):
            log_file = re.sub(r".mtz$", ".log", mtz_out_dm)
            script_dm = self.write_density_modification_script(mtz_for_dm, pdb_for_dm, mtz_out_dm, log_file)
            print('Running density modification, output written to %s. Please wait...' % (log_file))
            os.system("./%s" % (script_dm))

        return mtz_out_dm
    
    def get_phenix_version(self):
        """
        Weird construction to get the phenix version. This is required since some parameter names change between versions
        """
        # try:
        #     phenix_version = int(re.search(r"phenix-1\.(.+?)\.", miller.__file__).group(1)) #This is not so robust. relies on the format being 'phenix.1.18.something' or 'phenix.1.18-something'
        # except ValueError:
        #         phenix_version = int(re.search(r"phenix-1\.(.+?)\-", miller.__file__).group(1))
        # except AttributeError:
        #     print('Update phenix! Verify that you are using at least Phenix.1.19.')
        #     phenix_version = 20 #let's assume then that the latest phenix is installed in case this fails for other reasons than a very old phenix version
        
        #get phenix subversion, important since syntax can differ between versions
        phenix_version = get_phenix_version()
        try:
            phenix_subversion = int(phenix_version[2:4])
        except ValueError: #nightly versions can have no number (phenix-dev versions). Assume latest version so put number high
            phenix_subversion = 100

        return phenix_version, phenix_subversion
    
    def get_mtz_resolution(self, mtz_in):
        """
        Easy extraction of the resolution boundaries of an mtz file
        """
        reflections = any_file(mtz_in, force_type="hkl", raise_sorry_if_errors=True)
        low_res, high_res = reflections.file_content.file_content().max_min_resolution()
         
        return low_res, high_res

    
    def phenix_real_space_refinement_mtz(self, mtz_in, pdb_in, column_labels):
        """
        Real space refinement based on mtz file and specified column labels
        use Bash line to run phenix.real_space_refine as usual (use of os.system is bad practice).
        Some parameters have changed between version 1.7, 1.8 and 1.9 hence the weird construction to grep the version
        """       
        
        mtz_name = get_name(mtz_in)
            
        #Specify phenix version dependent parameters
        phenix_version, phenix_subversion = self.get_phenix_version()
        print("Phenix version: %s" %(phenix_version))

        if phenix_subversion >= 19:
            output_prefix = 'output.prefix=%s_independent'%(mtz_name)
            model_format  = 'model_format=pdb'
            # outpdb        = "%s_independent_real_space_refined_000.pdb"%(mtz_name)
        else:
            output_prefix = 'output.file_name_prefix=%s_independent'%(mtz_name)
            model_format  = 'output.model_format=pdb'
            # outpdb        = "%s_independent_real_space_refined.pdb"%(mtz_name)

        # print("phenix.real_space_refine %s %s %s %s %s %s label='%s'" % (self.real_phil,
        # mtz_in, self.additional, pdb_in, output_prefix, model_format, column_labels))
        real = os.system("phenix.real_space_refine %s %s %s %s %s %s label='%s' scattering_table=%s" %(
            self.real_phil, mtz_in, self.additional, pdb_in, output_prefix, model_format, column_labels, self.scattering_table))

        #Find output file
        if real == 0 : #os.system has correctly finished. Then search for the last refined structure
            if phenix_subversion >= 19:
                try:
                    pdb_fles = glob.glob("%s_independent_real_space_refined_???.pdb"%(mtz_name))
                    # [fle for fle in os.listdir(os.getcwd()) if "%s_independent_real_space_refined_0"%(mtz_name) in
                    #            fle and fle.endswith('pdb')]
                    pdb_fles.sort()
                    outpdb = pdb_fles[-1]
                except IndexError:
                    outpdb = "%s_independent_real_space_refined_000.pdb" % (mtz_name)
            else:
                try:
                    pdb_fles = glob.glob("%s_independent_real_space_refined.pdb"%(mtz_name))
                    # [fle for fle in os.listdir(os.getcwd()) if
                    #            "%s_independent_real_space_refined" % (mtz_name) in fle and fle.endswith('pdb')]
                    pdb_fles.sort()
                    outpdb = pdb_fles[-1]
                except IndexError:
                    outpdb = "%s_independent_real_space_refined.pdb" % (mtz_name)
        else: # os.systenm has nog correctly finished
            outpdb = "not_a_file"

        return outpdb
    
    
    def phenix_real_space_refinement_ccp4(self, ccp4_in, pdb_in, resolution):
        """
        Real space refinement based on ccp4 file
        use Bash line to run phenix.real_space_refine as usual (use of os.system is bad practice).
        Some parameters have changed between version 1.7, 1.8 and 1.9 hence the weird construction to grap the version
        """       
        
        ccp4_name = get_name(ccp4_in)
            
        #Specify phenix version dependent parameters
        phenix_version, phenix_subversion = self.get_phenix_version()
        print("Phenix version: %s" %(phenix_version))

        if phenix_subversion >= 19:
            output_prefix = 'output.prefix=%s_independent'%(ccp4_name)
            model_format  = 'model_format=pdb'
            # outpdb        = "%s_independent_real_space_refined_000.pdb"%(ccp4_name)
        else:
            output_prefix = 'output.file_name_prefix=%s_independent'%(ccp4_name)
            model_format  = 'output.model_format=pdb'
            # outpdb        = "%s_independent_real_space_refined.pdb"%(ccp4_name)

        #launch phenix.real_space_refine
        real = os.system("phenix.real_space_refine %s %s %s %s %s %s resolution=%.2f scattering_table=%s" %(
            self.real_phil, ccp4_in, self.additional, pdb_in, output_prefix, model_format, resolution, self.scattering_table))

        #Find output file
        if real == 0 : #os.system has correctly finished. Then search for the last refined structure
            if phenix_subversion >= 19:
                try:
                    pdb_fles = glob.glob("%s_independent_real_space_refined_???.pdb"%(ccp4_name))
                    # [fle for fle in os.listdir(os.getcwd()) if "%s_independent_real_space_refined_0"%(ccp4_name) in
                    #            fle and fle.endswith('pdb')]
                    pdb_fles.sort()
                    outpdb = pdb_fles[-1]
                except IndexError:
                    outpdb = "%s_independent_real_space_refined_000.pdb" % (ccp4_name)
            else:
                try:
                    pdb_fles = glob.glob("%s_independent_real_space_refined.pdb"%(ccp4_name))
                    # [fle for fle in os.listdir(os.getcwd()) if
                    #            "%s_independent_real_space_refined" % (ccp4_name) in fle and fle.endswith('pdb')]
                    pdb_fles.sort()
                    outpdb = pdb_fles[-1]
                except IndexError:
                    outpdb = "%s_independent_real_space_refined.pdb" % (ccp4_name)
        else: # os.systenm has nog correctly finished
            outpdb = "not_a_file"

        return outpdb

        
    def run_refinements(self):

        print("---Refinement %s with occupancy %.3f---" %(self.maptype, self.occ))
        print("---Refinement %s with occupancy %.3f---" %(self.maptype, self.occ), file=log)
        self.search_mtz_file_reciprocal_space_refinement()
        print("RECIPROCAL SPACE REFINEMENT WITH %s AND %s" %(self.mtz_file, self.pdb_in))
        mtz_out_rec, pdb_out_rec = self.phenix_reciprocal_space_refinement()
        print("Output reciprocal space refinement:", file=log)
        print("----------------")
        print("Output reciprocal space refinement:")
        if os.path.isfile(pdb_out_rec):
            print("    pdb-file: %s"%(pdb_out_rec), file=log)
            print("    pdb-file: %s"%(pdb_out_rec))
        else:
            print("    pdb-file not found, %s incorrectly returned" %(self.pdb_in), file=log)
            print("    pdb-file not found, %s incorrectly returned" %(self.pdb_in))
            pdb_out_rec = self.pdb_in
        if os.path.isfile(mtz_out_rec):
            print("    mtz-file: %s"%(mtz_out_rec), file=log)
            print("    mtz-file: %s"%(mtz_out_rec))
        else:
            print("    mtz-file not found. Refinement failed.", file=log)
            print("    mtz-file not found. Refinement failed.")
        print("----------------")
        
        if self.density_modification.density_modification:
            print("DENSITY MODIFICATION WITH %s AND %s" %(mtz_out_rec, pdb_out_rec))
            #mtz_dm = self.phenix_density_modification(mtz_out_rec, pdb_out_rec)
            mtz_dm = self.ccp4_dm(pdb_out_rec)
            print("Output density modification:", file=log)
            print("Output density modification:")
            if os.path.isfile(mtz_dm):
                print("    mtz-file: %s"%(mtz_dm), file=log)
                print("    mtz-file: %s" % (mtz_dm))
            else:
                print("    mtz-file not found. Density modification failed.", file=log)
                print("    mtz-file not found. Density modification failed.")

            
        self.search_mtz_map_real_space_refinement()
        maplabels = self.get_mtz_map_columns()
        print("REAL SPACE REFINEMENT WITH %s AND %s" %(self.mtz_map, self.pdb_in))
        pdb_out_real = self.phenix_real_space_refinement_mtz(self.mtz_map, self.pdb_in, maplabels)
        print("Output real space refinement:", file=log)
        print("----------------")
        print("Output real space refinement:")
        if os.path.isfile(pdb_out_real):
            print("    pdb-file: %s"%(pdb_out_real), file=log)
            print("    pdb-file: %s"%(pdb_out_real))
        else:
            print("    pdb-file not found, %s incorrectly returned" %(self.pdb_in), file=log)
            print("    pdb-file not found, %s incorrectly returned" %(self.pdb_in))
            pdb_out_real = self.pdb_in
        print("pdb_out_real", pdb_out_real)
        print("----------------")

        if self.density_modification.density_modification:
            ccp4_dm = re.sub(r".mtz$", ".ccp4", mtz_dm)
            _, high_res = self.get_mtz_resolution(mtz_dm)
            print("REAL SPACE REFINEMENT WITH %s AND %s" %(ccp4_dm, pdb_out_rec))
            pdb_out_rec_real = self.phenix_real_space_refinement_ccp4(ccp4_dm, pdb_out_rec, high_res)
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
            pdb_out_rec_real = self.phenix_real_space_refinement_mtz(mtz_out_rec, pdb_out_rec, '2FOFCWT,PH2FOFCWT')
            print("pdb_out_rec_real", pdb_out_rec_real)
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

# def get_unique_id(id_length=20):
#     """
#     Function to get a unique id based on UUID with length id_length
#     """
#     if id_length > 36:
#         id_length == 36
#     return str(uuid.uuid4())[:id_length]
# 
def generate_log_name(time_stamp):
    """
    Generate a unique name for the Xtrapol8 logfile.
    A short uuid of 20 characters is added to the logfile name.
    """
    uuid = get_unique_id(36)
    logname = "%s_Xtrapol8_refinements_%s.log" %(time_stamp, uuid)
    
    return logname

def remove_unique_id_from_log(log_name):
    """
    Remove the unqiue sequence from the log file
    """
    index = log_name.find("Xtrapol8_refinements")+len("Xtrapol8_refinements")
    new_name = log_name[:index]+".log"
    os.rename(log_name, new_name)
 
    return new_name

def run(args):
    version.VERSION
    now = datetime.now().strftime('%Y-%m-%d_%Hh%M')
    #Generate log-file. Needs to be created before the output directory is created and to be a global parameter in order to be easily used in all classes and functions
    logname = generate_log_name(now)
    global log
    log = open(logname, "w")

    print("Xtrapol8 refiner -- Xtrapol8 version %s -- run date: %s" %(version.VERSION, now), file=log)
    print("Phenix version: {:s}".format(get_phenix_version()), file=log)
    print("CCP4 version: {:s}".format(get_ccp4_version()), file=log)
    print("Xtrapol8 refiner -- Xtrapol8 version %s -- run date: %s" %(version.VERSION, now))
    print("Phenix version: {:s}".format(get_phenix_version()))
    print("CCP4 version: {:s}".format(get_ccp4_version()))

    print('-----------------------------------------')
    #parse arguments
    if len(args) == 0 :
        master_phil.show(attributes_level=1)
        raise Usage("phenix.python Xtrapol8_refinements.py + [.phil] + [arguments]\n arguments only overwrite .phil if provided last")
    
    log_dir = os.getcwd()
    
    #Extract input from inputfile and command line
    #argument_interpreter = master_phil.command_line_argument_interpreter(home_scope="input")
    input_objects = iotbx.phil.process_command_line_with_files(
        args=args,
        master_phil=master_phil
        )
    params = input_objects.work.extract()
    #modified_phil = master_phil.format(python_object=params)

    #Extract info form Xtrapol8 run
    Xtrapol8_input_objects = iotbx.phil.process_command_line_with_files(
        args = [params.input.Xtrapol8_out],
        master_phil = Xtrapol8_master_phil
        )
    Xtrapol8_params = Xtrapol8_input_objects.work.extract()
    
    
    #Extract info on extrapolated structure factor types:
    #specify extrapolated structure factors and map types
    qFextr_map = qFgenick_map = qFextr_calc_map = Fextr_map = Fgenick_map = Fextr_calc_map = kFextr_map = kFgenick_map = kFextr_calc_map = False
    if Xtrapol8_params.f_and_maps.all_maps: #calculate all Fextr map types
        qFextr_map = qFgenick_map = qFextr_calc_map = Fextr_map = Fgenick_map = Fextr_calc_map = kFextr_map = kFgenick_map = kFextr_calc_map = True
    else: #calculate only the specified map types
        if 'qfextr' in Xtrapol8_params.f_and_maps.f_extrapolated_and_maps:
            qFextr_map      = True
        if 'fextr' in Xtrapol8_params.f_and_maps.f_extrapolated_and_maps:
            Fextr_map       = True
        if 'kfextr' in Xtrapol8_params.f_and_maps.f_extrapolated_and_maps:
            kFextr_map      = True
        if 'qfgenick' in Xtrapol8_params.f_and_maps.f_extrapolated_and_maps:
            qFgenick_map       = True
        if 'kfgenick' in Xtrapol8_params.f_and_maps.f_extrapolated_and_maps:
            kFgenick_map       = True
        if 'fgenick' in Xtrapol8_params.f_and_maps.f_extrapolated_and_maps:
            Fgenick_map        = True
        if ('qfextr_calc') in Xtrapol8_params.f_and_maps.f_extrapolated_and_maps:
            qFextr_calc_map = True
        if ('fextr_calc') in Xtrapol8_params.f_and_maps.f_extrapolated_and_maps:
            Fextr_calc_map  = True
        if ('kfextr_calc') in Xtrapol8_params.f_and_maps.f_extrapolated_and_maps:
            kFextr_calc_map = True
            
    #Only run the non-weighted
    if Xtrapol8_params.f_and_maps.only_no_weight:
        qFextr_map = qFgenick_map = qFextr_calc_map = kFextr_map = kFgenick_map = kFextr_calc_map = False
        Fextr_map  = Fgenick_map  = Fextr_calc_map = True
    #Only run the k-weighted
    if Xtrapol8_params.f_and_maps.only_kweight:
        qFextr_map = qFgenick_map = qFextr_calc_map = Fextr_map = Fgenick_map = Fextr_calc_map = False
        kFextr_map = kFgenick_map = kFextr_calc_map = True
    #Only run the q-weighted
    if Xtrapol8_params.f_and_maps.only_qweight: 
        Fextr_map = Fgenick_map = Fextr_calc_map = kFextr_map = kFgenick_map = kFextr_calc_map = False
        qFextr_map = qFgenick_map = qFextr_calc_map = True
    #Run all maps
    if Xtrapol8_params.f_and_maps.all_maps: #calculate all Fextr map types
        qFextr_map = qFgenick_map = qFextr_calc_map = Fextr_map = Fgenick_map = Fextr_calc_map = kFextr_map = kFgenick_map = kFextr_calc_map = True

    #convert the old use_occupancy_from_distance_analysis to the new occupancy_estimation keyword
    if Xtrapol8_params.map_explorer.use_occupancy_from_distance_analysis == True:
        Xtrapol8_params.map_explorer.occupancy_estimation = "distance_analysis"

    #if all map types being false:
    if qFextr_map == qFgenick_map == qFextr_calc_map == Fextr_map == Fgenick_map == Fextr_calc_map == kFextr_map == kFgenick_map == kFextr_calc_map == False and Xtrapol8_params.f_and_maps.fast_and_furious == False: 
        print('The combination of arguments used to define extrapolated structure factors and maps leads to no calculations at all. The default will be applied: qFextr', file=log)
        print('The combination of arguments used to define extrapolated structure factors and maps leads to no calculations at all. The default will be applied: qFextr')
        qFextr_map = True
    #if fast_and_furious mode: overwrite all F and map setting to default:
    if Xtrapol8_params.f_and_maps.fast_and_furious:
        #change parameters for Xtrapol8_out.phil
        Xtrapol8_params.f_and_maps.fofo_type = 'qfofo'
        Xtrapol8_params.f_and_maps.f_extrapolated_and_maps = ['qfextr']
        Xtrapol8_params.f_and_maps.only_no_weight = Xtrapol8_params.f_and_maps.all_maps = Xtrapol8_params.f_and_maps.only_kweight = Xtrapol8_params.f_and_maps.only_qweight = False
        #change working parameters
        qFoFo_weight = qFextr_map = True
        kFoFo_weight = qFgenick_map = qFextr_calc_map = Fextr_map = Fgenick_map = Fextr_calc_map = kFextr_map = kFgenick_map = kFextr_calc_map = False
        Xtrapol8_params.map_explorer.use_occupancy_from_distance_analysis = False
        if Xtrapol8_params.map_explorer.occupancy_estimation == 'distance_analysis':
            Xtrapol8_params.map_explorer.occupancy_estimation = 'difference_map_maximization'
        Xtrapol8_params.f_and_maps.negative_and_missing='truncate_and_fill'
                
    #Bring all maptypes to be calculated together in list instead of using loose varaibles:
    all_maptypes   = ['qFextr_map','Fextr_map', 'qFgenick_map', 'Fgenick_map', 'qFextr_calc_map', 'Fextr_calc_map', 'kFextr_map', 'kFgenick_map', 'kFextr_calc_map']
    all_maps       = [qFextr_map, Fextr_map, qFgenick_map, Fgenick_map, qFextr_calc_map, Fextr_calc_map, kFextr_map, kFgenick_map, kFextr_calc_map]
    maptypes_zip   = zip(all_maptypes, all_maps)
    final_maptypes = [mp[0] for mp in maptypes_zip if mp[1] == True]
    
    
    #Add all arguments to log-file
    print('-----------------------------------------', file=log)
    print("ARGUMENTS", file=log)
    print('-----------------------------------------', file=log)
    
    print('-----------------------------------------')
    print("ARGUMENTS")
    print('-----------------------------------------')

    modified_phil = master_phil.format(python_object=params)
    modified_phil.show()
    modified_phil.show(out=log)
    ##get the differences with the default values and only show these in the log-file
    #diff_phil = master_phil.fetch_diff(source=modified_phil)
    #diff_phil.show()
    #diff_phil.show(out=log)


    print('-----------------------------------------')
    print('DATA PREPARATION')
    print('-----------------------------------------')

    print('-----------------------------------------', file=log)
    print('DATA PREPARATION', file=log)
    print('-----------------------------------------', file=log)

    #remember starting directory:
    startdir = os.getcwd()
    DH = DataHandler(params.input.model_pdb,
                     Xtrapol8_params.input.additional_files,
                     Xtrapol8_params.output.outdir,
                     Xtrapol8_params.input.reference_pdb,
                     params.refinement.reciprocal_space_phil,
                     params.refinement.real_space_phil,
                     params.map_explorer.residue_list)
    err, err_m, warning, warning_m = DH.check_all_files_and_dirs()
    if err == 1:
        print(err_m)
        sys.exit()
    elif warning == 1:
        print(warning_m)
        print(warning_m, file=log)
    else:
        print("All files found.")
        print("All files found.", file=log)

    if DH.X8_pdb_in == None:
        params.map_explorer.do_distance_analysis = False
        
    #Move to the output directory
    os.chdir(DH.outdir)

    #Move log file to output directory
    # print(log.name)
    full_log = "%s/%s" %(log_dir, log.name)
    #print(os.path.isfile(full_log))
    if os.path.isfile(full_log):
        shutil.move(full_log, full_log.replace(log_dir,DH.outdir))
    log_name = remove_unique_id_from_log(log.name)
    full_log = "%s/%s" %(DH.outdir, log_name)
        
    pymol_script_out = "%s/pymol_movie_refiner.py" %(DH.outdir)



    #get the output filename from the Xtrapol8_out.phil
    outname = Xtrapol8_params.output.outname

    #Initiate the lists to store the output
    if params.map_explorer.do_distance_analysis:
        if qFextr_map:
            qFextr_recref_mtz_lst = []
            qFextr_recref_pdb_lst = [DH.X8_pdb_in]
            qFextr_realref_lst    = [DH.X8_pdb_in]
            qFextr_recrealref_lst = [DH.X8_pdb_in]
        if Fextr_map:
            Fextr_recref_mtz_lst = []
            Fextr_recref_pdb_lst = [DH.X8_pdb_in]
            Fextr_recrealref_lst = [DH.X8_pdb_in]
            Fextr_realref_lst    = [DH.X8_pdb_in]
        if kFextr_map:
            kFextr_recref_mtz_lst = []
            kFextr_recref_pdb_lst = [DH.X8_pdb_in]
            kFextr_realref_lst    = [DH.X8_pdb_in]
            kFextr_recrealref_lst = [DH.X8_pdb_in]
        if qFgenick_map:
            qFgenick_recref_mtz_lst = []
            qFgenick_recref_pdb_lst = [DH.X8_pdb_in]
            qFgenick_realref_lst    = [DH.X8_pdb_in]
            qFgenick_recrealref_lst = [DH.X8_pdb_in]
        if Fgenick_map:
            Fgenick_recref_mtz_lst = []
            Fgenick_recref_pdb_lst = [DH.X8_pdb_in]
            Fgenick_realref_lst    = [DH.X8_pdb_in]
            Fgenick_recrealref_lst = [DH.X8_pdb_in]
        if kFgenick_map:
            kFgenick_recref_mtz_lst = []
            kFgenick_recref_pdb_lst = [DH.X8_pdb_in]
            kFgenick_realref_lst    = [DH.X8_pdb_in]
            kFgenick_recrealref_lst = [DH.X8_pdb_in]
        if qFextr_calc_map:
            qFextr_calc_recref_mtz_lst = []
            qFextr_calc_recref_pdb_lst = [DH.X8_pdb_in]
            qFextr_calc_realref_lst    = [DH.X8_pdb_in]
            qFextr_calc_recrealref_lst = [DH.X8_pdb_in]
        if Fextr_calc_map:
            Fextr_calc_recref_mtz_lst = []
            Fextr_calc_recref_pdb_lst = [DH.X8_pdb_in]
            Fextr_calc_realref_lst    = [DH.X8_pdb_in]
            Fextr_calc_recrealref_lst = [DH.X8_pdb_in]
        if kFextr_calc_map:
            kFextr_calc_recref_mtz_lst = []
            kFextr_calc_recref_pdb_lst = [DH.X8_pdb_in]
            kFextr_calc_realref_lst    = [DH.X8_pdb_in]
            kFextr_calc_recrealref_lst = [DH.X8_pdb_in]
    else:
        if qFextr_map:
            qFextr_recref_mtz_lst = []
            qFextr_recref_pdb_lst = []
            qFextr_realref_lst    = []
            qFextr_recrealref_lst = []
        if Fextr_map:
            Fextr_recref_mtz_lst = []
            Fextr_recref_pdb_lst = []
            Fextr_recrealref_lst = []
            Fextr_realref_lst    = []
        if kFextr_map:
            kFextr_recref_mtz_lst = []
            kFextr_recref_pdb_lst = []
            kFextr_realref_lst    = []
            kFextr_recrealref_lst = []
        if qFgenick_map:
            qFgenick_recref_mtz_lst = []
            qFgenick_recref_pdb_lst = []
            qFgenick_realref_lst    = []
            qFgenick_recrealref_lst = []
        if Fgenick_map:
            Fgenick_recref_mtz_lst = []
            Fgenick_recref_pdb_lst = []
            Fgenick_realref_lst    = []
            Fgenick_recrealref_lst = []
        if kFgenick_map:
            kFgenick_recref_mtz_lst = []
            kFgenick_recref_pdb_lst = []
            kFgenick_realref_lst    = []
            kFgenick_recrealref_lst = []
        if qFextr_calc_map:
            qFextr_calc_recref_mtz_lst = []
            qFextr_calc_recref_pdb_lst = []
            qFextr_calc_realref_lst    = []
            qFextr_calc_recrealref_lst = []
        if Fextr_calc_map:
            Fextr_calc_recref_mtz_lst = []
            Fextr_calc_recref_pdb_lst = []
            Fextr_calc_realref_lst    = []
            Fextr_calc_recrealref_lst = []
        if kFextr_calc_map:
            kFextr_calc_recref_mtz_lst = []
            kFextr_calc_recref_pdb_lst = []
            kFextr_calc_realref_lst    = []
            kFextr_calc_recrealref_lst = []
        

    print('-----------------------------------------')
    print('DATA PREPARATION DONE')
    print('-----------------------------------------')
    print("REFINE THE EXTRAPOLATED MODELS")
    print('-----------------------------------------')

    print('-----------------------------------------', file=log)
    print("REFINE THE EXTRAPOLATED MODELS", file=log)
    print('-----------------------------------------', file=log)

    #Run the refinements:
    for occ in Xtrapol8_params.occupancies.list_occ:
        for maptype in final_maptypes:
            if maptype in ('qFextr_map','qFgenick_map','qFextr_calc_map'):
                os.chdir("qweight_occupancy_%.3f" %(occ))
            elif maptype in ('kFextr_map','kFgenick_map','kFextr_calc_map'):
                os.chdir("kweight_occupancy_%.3f" %(occ))
            else:
                os.chdir("occupancy_%.3f" %(occ))
            mtz_out, pdb_rec, pdb_real, pdb_rec_real = Refiner(DH.refine_pdb_in,
                                                        maptype,
                                                        outname,
                                                        occ,
                                                        DH.additional,
                                                        reciprocal_phil = DH.reciprocal_space_phil,
                                                        real_phil = DH.real_space_phil,
                                                        scattering_table = Xtrapol8_params.scattering_table,
                                                        density_modification = Xtrapol8_params.refinement.phenix_keywords.density_modification).run_refinements()
            
            print("--------------", file=log)
            #depending on the map-type, append the refinement output to the correct list
            #this is ugly, TODO: make an object to store the results in a clean and transparant way
            #if params.map_explorer.do_distance_analysis:
            if maptype == 'qFextr_map':
                append_if_file_exist(qFextr_recref_mtz_lst, os.path.abspath(mtz_out))
                append_if_file_exist(qFextr_recref_pdb_lst, os.path.abspath(pdb_rec))
                append_if_file_exist(qFextr_recrealref_lst, os.path.abspath(pdb_rec_real))
                append_if_file_exist(qFextr_realref_lst, os.path.abspath(pdb_real))
            elif maptype == 'qFgenick_map':
                append_if_file_exist(qFgenick_recref_mtz_lst, os.path.abspath(mtz_out))
                append_if_file_exist(qFgenick_recref_pdb_lst, os.path.abspath(pdb_rec))
                append_if_file_exist(qFgenick_recrealref_lst, os.path.abspath(pdb_rec_real))
                append_if_file_exist(qFgenick_realref_lst, os.path.abspath(pdb_real))
            elif maptype == 'qFextr_calc_map':
                append_if_file_exist(qFextr_calc_recref_mtz_lst, os.path.abspath(mtz_out))
                append_if_file_exist(qFextr_calc_recref_pdb_lst, os.path.abspath(pdb_rec))
                append_if_file_exist(qFextr_calc_recrealref_lst, os.path.abspath(pdb_rec_real))
                append_if_file_exist(qFextr_calc_realref_lst, os.path.abspath(pdb_real))
            elif maptype == 'Fextr_map':
                append_if_file_exist(Fextr_recref_mtz_lst, os.path.abspath(mtz_out))
                append_if_file_exist(Fextr_recref_pdb_lst, os.path.abspath(pdb_rec))
                append_if_file_exist(Fextr_recrealref_lst, os.path.abspath(pdb_rec_real))
                append_if_file_exist(Fextr_realref_lst, os.path.abspath(pdb_real))
            elif maptype == 'Fgenick_map':
                append_if_file_exist(Fgenick_recref_mtz_lst, os.path.abspath(mtz_out))
                append_if_file_exist(Fgenick_recref_pdb_lst, os.path.abspath(pdb_rec))
                append_if_file_exist(Fgenick_recrealref_lst, os.path.abspath(pdb_rec_real))
                append_if_file_exist(Fgenick_realref_lst, os.path.abspath(pdb_real))
            elif maptype == 'Fextr_calc_map':
                append_if_file_exist(Fextr_calc_recref_mtz_lst, os.path.abspath(mtz_out))
                append_if_file_exist(Fextr_calc_recref_pdb_lst, os.path.abspath(pdb_rec))
                append_if_file_exist(Fextr_calc_recrealref_lst, os.path.abspath(pdb_rec_real))
                append_if_file_exist(Fextr_calc_realref_lst, os.path.abspath(pdb_real))
            elif maptype == 'kFextr_map':
                append_if_file_exist(kFextr_recref_mtz_lst, os.path.abspath(mtz_out))
                append_if_file_exist(kFextr_recref_pdb_lst, os.path.abspath(pdb_rec))
                append_if_file_exist(kFextr_recrealref_lst, os.path.abspath(pdb_rec_real))
                append_if_file_exist(kFextr_realref_lst, os.path.abspath(pdb_real))
            elif maptype == 'kFgenick_map':
                append_if_file_exist(kFgenick_recref_mtz_lst, os.path.abspath(mtz_out))
                append_if_file_exist(kFgenick_recref_pdb_lst, os.path.abspath(pdb_rec))
                append_if_file_exist(kFgenick_recrealref_lst, os.path.abspath(pdb_rec_real))
                append_if_file_exist(kFgenick_realref_lst, os.path.abspath(pdb_real))
            elif maptype == 'kFextr_calc_map':
                append_if_file_exist(kFextr_calc_recref_mtz_lst, os.path.abspath(mtz_out))
                append_if_file_exist(kFextr_calc_recref_pdb_lst, os.path.abspath(pdb_rec))
                append_if_file_exist(kFextr_calc_recrealref_lst, os.path.abspath(pdb_rec_real))
                append_if_file_exist(kFextr_calc_realref_lst, os.path.abspath(pdb_real))
            os.chdir(DH.outdir)

    print('-----------------------------------------')
    print("REFINE THE EXTRAPOLATED MODELS DONE")
    print('-----------------------------------------')

    if params.map_explorer.do_distance_analysis:
        print("ESTIMATE OPTIMAL OCCUPANCY USING THE DISTANCE METHOD")
        print('-----------------------------------------')

        print('-----------------------------------------', file=log)
        print("ESTIMATE OPTIMAL OCCUPANCY USING THE DISTANCE METHOD", file=log)
        print('-----------------------------------------', file=log)

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
        reversed_maptypes = final_maptypes[:]
        reversed_maptypes.reverse()
        occ_overview = {}
        for mp in reversed_maptypes:
            mp_type = mp.split("_map")[0]
            print("---%s---"%(mp_type))
            print("---%s---"%(mp_type), file= log)

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
            elif mp == 'qFgenick_map':
                recref_mtz_lst = qFgenick_recref_mtz_lst
                recref_pdb_lst = qFgenick_recref_pdb_lst
                recrealref_lst = qFgenick_recrealref_lst
                realref_lst    = qFgenick_realref_lst
            elif mp == 'qFextr_calc_map':
                recref_mtz_lst = qFextr_calc_recref_mtz_lst
                recref_pdb_lst = qFextr_calc_recref_pdb_lst
                recrealref_lst = qFextr_calc_recrealref_lst
                realref_lst    = qFextr_calc_realref_lst
            elif mp == 'Fextr_map':
                recref_mtz_lst = Fextr_recref_mtz_lst
                recref_pdb_lst = Fextr_recref_pdb_lst
                recrealref_lst = Fextr_recrealref_lst
                realref_lst    = Fextr_realref_lst
            elif mp == 'Fgenick_map':
                recref_mtz_lst = Fgenick_recref_mtz_lst
                recref_pdb_lst = Fgenick_recref_pdb_lst
                recrealref_lst = Fgenick_recrealref_lst
                realref_lst    = Fgenick_realref_lst
            elif mp == 'Fextr_calc_map':
                recref_mtz_lst = Fextr_calc_recref_mtz_lst
                recref_pdb_lst = Fextr_calc_recref_pdb_lst
                recrealref_lst = Fextr_calc_recrealref_lst
                realref_lst    = Fextr_calc_realref_lst
            elif mp == 'kFextr_map':
                recref_mtz_lst = kFextr_recref_mtz_lst
                recref_pdb_lst = kFextr_recref_pdb_lst
                recrealref_lst = kFextr_recrealref_lst
                realref_lst    = kFextr_realref_lst
            elif mp == 'kFgenick_map':
                recref_mtz_lst = kFgenick_recref_mtz_lst
                recref_pdb_lst = kFgenick_recref_pdb_lst
                recrealref_lst = kFgenick_recrealref_lst
                realref_lst    = kFgenick_realref_lst
            elif mp == 'kFextr_calc_map':
                recref_mtz_lst = kFextr_calc_recref_mtz_lst
                recref_pdb_lst = kFextr_calc_recref_pdb_lst
                recrealref_lst = kFextr_calc_recrealref_lst
                realref_lst    = kFextr_calc_realref_lst

            #Make a plot of the refinement R-factors, related to the specific maptype. The log-files should have the same prefix as the mtz-files.
            #This assumption is made in order to avoid storing the log-files in even another list
            #test = map(lambda fle: re.sub(r'mtz$','log', fle), recref_mtz_lst)
            #print(test)
            plot_Rfactors_per_alpha(map(lambda fle: re.sub(r'mtz$','log', fle), recref_mtz_lst), mp_type)
            print("", file=log)
            print("")

            #Check if the list with PDB files from real space refinement after reciprocal space refinement is complete
            #   If this is not the case (when reciprocal space or real space refinement had failed),
            #   the PDB files from real space refinment (direct real space refinement in the extrapolated maps) are used
            if ((len(recrealref_lst)==len(set(recrealref_lst)))):
                # print("Using recrealref_lst")
                pdb_list = recrealref_lst
            else:
                pdb_list = realref_lst
                # print("Using realref_lst")

            #Estimate alpha and occupancy based on the distances of the list with PDB files with the input model.
            #Don't use the water molecules to do the distance analysis because we don't know whether or not water was updated during the refinement
            # print("pdb_list", pdb_list)
            alpha, occ = Distance_analysis(pdb_list, Xtrapol8_params.occupancies.list_occ, resids_lst=
            DH.residue_list, use_waters = False, outsuffix = mp_type, log = log).extract_alpha()
            print("", file=log)
            print("")

            #add models and maps to pymol movie.
            #Not elegant with all the with hard coded regular expressions but avoids storing other lists
            pymol_pdb_list = pdb_list[:]
            pymol_pdb_list.remove(DH.X8_pdb_in)
            pymol_mtz_list = recref_mtz_lst[:]
            #Make Pymol movie with the reciprocal space refined maps if recrealref_lst is complete
            #Otherwise use the real space refined models + direct maps
            if pdb_list == recrealref_lst:
                ccp4_list = map(lambda fle: re.sub(r".mtz$", "_2mFo-DFc_filled.ccp4", fle), pymol_mtz_list)
                model_label = '%s_reciprocal_real_space'%(mp_type)
                ccp4_map_label = '%s_reciprocal_space'%(mp)
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
            if len(ccp4_list) == len(pymol_pdb_list) == len(Xtrapol8_params.occupancies.list_occ):
                #print("Pymol movie residue list: %s" %(DH.residue_list))
                PM = Pymol_movie(Xtrapol8_params.occupancies.list_occ, pdblst=pymol_pdb_list, ccp4_maps = ccp4_list,
                            resids_lst = DH.residue_list, model_label=model_label, ccp4_map_label=ccp4_map_label)
            else:
                #print("Pymol movie residue list: %s" %(DH.residue_list))
                PM = Pymol_movie(Xtrapol8_params.occupancies.list_occ, pdblst=pymol_pdb_list, resids_lst = DH.residue_list,
                            model_label=model_label)
            PM.write_pymol_script(outfile=pymol_script_out)
            
            #If estimated occupancy if not in list (will be case when using distance analysis or when plotalpha fails), take the closest occupancy from the list
            if occ not in Xtrapol8_params.occupancies.list_occ:
                occ = min(Xtrapol8_params.occupancies.list_occ, key = lambda x: abs(x-occ))
                alpha = 1/occ
            occ_dir = "%s/%s_%.3f" %(DH.outdir, dir_prefix, occ)

            #ddm calculation
            pdb_for_ddm = pdb_list[Xtrapol8_params.occupancies.list_occ.index(occ)+1]
            print("----Generate distance difference plot----")
            Difference_distance_analysis(DH.X8_pdb_in, pdb_for_ddm, outdir=occ_dir, scale=1.5).ddms()
            print('---------------------------')

            #Coot script
            mtzs_for_coot  = []
            if len(recref_mtz_lst) != len(Xtrapol8_params.occupancies.list_occ):
                occ_list_mtz = [(re.search(r'occupancy\_(.+?)\/',mtz).group(1)) for mtz in recref_mtz_lst]
                mtz_rec = recref_mtz_lst[occ_list_mtz.index(occ)]
            else:
                mtz_rec = recref_mtz_lst[Xtrapol8_params.occupancies.list_occ.index(occ)]
            append_if_file_exist(mtzs_for_coot, mtz_rec)
            if Xtrapol8_params.refinement.phenix_keywords.density_modification.density_modification:
                #mtz_dm = re.sub(".mtz$","_densitymod.mtz", mtz_rec)
                mtz_dm = re.sub(".mtz$","_dm.mtz", mtz_rec) #probably wrong because the name of the extrapolated structure factors is used and not the refined mtz
                append_if_file_exist(mtzs_for_coot, mtz_dm)

            #elif Xtrapol8_params.refinement.refmac_keywords.density_modification.density_modification:
                #mtz_dm = re.sub(".mtz$","_dm.mtz", mtz_rec) #probably wrong because the name of the extrapolated structure factors is used and not the refined mtz
                #append_if_file_exist(mtzs_for_coot, mtz_dm)
            mtz_extr = ["%s/%s"%(occ_dir,fle) for fle in os.listdir(occ_dir) if outname in fle and fle.endswith('m%s-DFc.mtz'%(mp_type))][0]
            append_if_file_exist(mtzs_for_coot,os.path.abspath(mtz_extr))

            pdbs_for_coot = [DH.X8_pdb_in,
                    recref_pdb_lst[Xtrapol8_params.occupancies.list_occ.index(occ)+1],
                    realref_lst[Xtrapol8_params.occupancies.list_occ.index(occ)+1],
                    recrealref_lst[Xtrapol8_params.occupancies.list_occ.index(occ)+1]]
            #if outname == 'triggered': #if dummy name applied, the files still contain the dummy name
                #mtzs_for_coot = map(lambda fle: re.sub(r"triggered",params.output.outname, fle), mtzs_for_coot)
                #pdbs_for_coot = map(lambda fle: re.sub(r"triggered",params.output.outname, fle), pdbs_for_coot)
            #find the mtz file of the FoFo:
            FoFo = [fle for fle in os.listdir(DH.outdir) if fle.endswith("FoFo.mtz")][0]
            script_coot = open_all_in_coot(DH.outdir+"/"+FoFo, pdbs_for_coot, mtzs_for_coot, DH.additional,
                                           occ_dir, mp_type)

            occ_overview[mp_type] = occ

            print("------------------------------------")
            print("------------------------------------", file=log)

        ##Add final lines to Pymol_script
        #if os.path.isfile('%s/pymol_movie.py' %(DH.outdir)):
            #Pymol_movie(Xtrapol8_params.occupancies.list_occ, resids_lst = DH.residue_list).write_pymol_appearance(
                #'%s/pymol_movie.py' %(
                #DH.outdir))

        print("Summary of occupancy determination:", file=log)
        print("Map type       Occupancy", file=log)
        for k in occ_overview:
            print("{:<15} {:>5.3f}".format(k, occ_overview[k]), file=log)
        print("-> Optimal occupancy of triggered state %.3f." %(occ), file=log)

        print("OCCUAPNCY DETERMINATION SUMMARY")
        print("Map type       Occupancy")
        for k in occ_overview:
            print("{:<15} {:>5.3f}".format(k, occ_overview[k]))
        print("-> Optimal occupancy of triggered state %.3f." %(occ))

        print('-----------------------------------------')
        print("ESTIMATE OPTIMAL OCCUPANCY USING THE DISTANCE METHOD DONE")
        print('-----------------------------------------')
    
    else:
 
        print("NO OCCUPANCY ESTIMATION USING DISTANCE METHOD")
        print('-----------------------------------------')

        ################################################################
        #No occupancy estimation but only generation of R-factor plot and pymol script
        #Therefore, here we loop over the maptypes (in a reversed order to keep with code above)
        #1) transfer the stored output files to a generic list
        #2) plot the refinement R-factors
        #3) append the refinement models and maps to the Pymol_movie script
        reversed_maptypes = final_maptypes[:]
        reversed_maptypes.reverse()
        occ_overview = {}
        for mp in reversed_maptypes:
            mp_type = mp.split("_map")[0]
            #print("---%s---"%(mp_type))

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
            elif mp == 'qFgenick_map':
                recref_mtz_lst = qFgenick_recref_mtz_lst
                recref_pdb_lst = qFgenick_recref_pdb_lst
                recrealref_lst = qFgenick_recrealref_lst
                realref_lst    = qFgenick_realref_lst
            elif mp == 'qFextr_calc_map':
                recref_mtz_lst = qFextr_calc_recref_mtz_lst
                recref_pdb_lst = qFextr_calc_recref_pdb_lst
                recrealref_lst = qFextr_calc_recrealref_lst
                realref_lst    = qFextr_calc_realref_lst
            elif mp == 'Fextr_map':
                recref_mtz_lst = Fextr_recref_mtz_lst
                recref_pdb_lst = Fextr_recref_pdb_lst
                recrealref_lst = Fextr_recrealref_lst
                realref_lst    = Fextr_realref_lst
            elif mp == 'Fgenick_map':
                recref_mtz_lst = Fgenick_recref_mtz_lst
                recref_pdb_lst = Fgenick_recref_pdb_lst
                recrealref_lst = Fgenick_recrealref_lst
                realref_lst    = Fgenick_realref_lst
            elif mp == 'Fextr_calc_map':
                recref_mtz_lst = Fextr_calc_recref_mtz_lst
                recref_pdb_lst = Fextr_calc_recref_pdb_lst
                recrealref_lst = Fextr_calc_recrealref_lst
                realref_lst    = Fextr_calc_realref_lst
            elif mp == 'kFextr_map':
                recref_mtz_lst = kFextr_recref_mtz_lst
                recref_pdb_lst = kFextr_recref_pdb_lst
                recrealref_lst = kFextr_recrealref_lst
                realref_lst    = kFextr_realref_lst
            elif mp == 'kFgenick_map':
                recref_mtz_lst = kFgenick_recref_mtz_lst
                recref_pdb_lst = kFgenick_recref_pdb_lst
                recrealref_lst = kFgenick_recrealref_lst
                realref_lst    = kFgenick_realref_lst
            elif mp == 'kFextr_calc_map':
                recref_mtz_lst = kFextr_calc_recref_mtz_lst
                recref_pdb_lst = kFextr_calc_recref_pdb_lst
                recrealref_lst = kFextr_calc_recrealref_lst
                realref_lst    = kFextr_calc_realref_lst
                
            #Make a plot of the refinement R-factors, related to the specific maptype. The log-files should have the same prefix as the mtz-files.
            #This assumption is made in order to avoid storing the log-files in even another list
            test = map(lambda fle: re.sub(r'mtz$','log', fle), recref_mtz_lst)
            #print(test)
            plot_Rfactors_per_alpha(map(lambda fle: re.sub(r'mtz$','log', fle), recref_mtz_lst), mp_type)
            print("", file=log)
            print("")

            #Check if the list with PDB files from real space refinement after reciprocal space refinement is complete
            #   If this is not the case (when reciprocal space or real space refinement had failed),
            #   the PDB files from real space refinment (direct real space refinement in the extrapolated maps) are used
            if ((len(recrealref_lst)==len(set(recrealref_lst)))):
                # print("Using recrealref_lst")
                pdb_list = recrealref_lst
            else:
                pdb_list = realref_lst
                # print("Using realref_lst")

            #add models and maps to pymol movie.
            #Not elegant with all the with hard coded regular expressions but avoids storing other lists
            pymol_pdb_list = pdb_list[:]
            pymol_mtz_list = recref_mtz_lst[:]
            #Make Pymol movie with the reciprocal space refined maps if recrealref_lst is complete
            #Otherwise use the real space refined models + direct maps
            if pdb_list == recrealref_lst:
                ccp4_list = map(lambda fle: re.sub(r".mtz$", "_2mFo-DFc_filled.ccp4", fle), pymol_mtz_list)
                model_label = '%s_reciprocal_real_space'%(mp_type)
                ccp4_map_label = '%s_reciprocal_space'%(mp)
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
            if len(ccp4_list) == len(pymol_pdb_list) == len(Xtrapol8_params.occupancies.list_occ):
                #print("Pymol movie residue list: %s" %(DH.residue_list))
                PM = Pymol_movie(Xtrapol8_params.occupancies.list_occ, pdblst=pymol_pdb_list, ccp4_maps = ccp4_list,
                             model_label=model_label, ccp4_map_label=ccp4_map_label, resids_lst = DH.residue_list)
            else:
                #print("Pymol movie residue list: %s" %(DH.residue_list))
                PM = Pymol_movie(Xtrapol8_params.occupancies.list_occ, pdblst=pymol_pdb_list, model_label=model_label, resids_lst = DH.residue_list)
            PM.write_pymol_script(outfile=pymol_script_out)
                
        #Add final lines to Pymol_script
        if os.path.isfile('%s/pymol_movie.py' %(DH.outdir)):
            PM.write_pymol_appearance(pymol_script_out)

    print("DONE! PLEASE INSPECT THE OUTPUT.")
    print('-----------------------------------------')

    print('-----------------------------------------', file=log)
    print("DONE! PLEASE INSPECT THE OUTPUT.", file=log)
    print('-----------------------------------------', file=log)


if __name__ == '__main__':
    run(sys.argv[1:])
    print('Hoera')
