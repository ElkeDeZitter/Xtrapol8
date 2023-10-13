# -*- coding: utf-8 -*-
"""
Run phenix for reciprocal and real space refinement.

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
import re, sys
import os
import glob
import iotbx.pdb
from mmtbx.scaling.matthews import p_vm_calculator
from libtbx import adopt_init_args
from cctbx import miller
from Fextr_utils import get_name
import subprocess
from iotbx.file_reader import any_file


class Phenix_refinements(object):
    def __init__(self,
                 mtz_in,
                 pdb_in,
                 additional                     = '',
                 F_column_labels                = 'QFEXTR',
                 strategy                       = ['individual_sites','individual_adp'],
                 rec_cycles                     = 5,
                 real_cycles                    = 5,
                 wxc_scale                      = 0.5,
                 wxu_scale                      = 1.0,
                 solvent                        = False,
                 sim_annealing                  = False,
                 sim_annealing_pars             = {},
                 map_sharpening                 = False,
                 weight_sel_crit                = {},
                 additional_reciprocal_keywords = [],
                 additional_real_keywords       = [],
                 log                            = sys.stdout):
        
        adopt_init_args(self, locals())
        
        self.strategy = "+".join(strategy)
        
        if map_sharpening:
            self.params = self.generate_map_params_bsharpening()
        else:
            self.params = ""
            
        self.mtz_name = get_name(self.mtz_in)

    def phenix_reciprocal_space_refinement(self):
        """
        use Bash line to run phenix.refine as usual (use of os.system is bad practice)
        """
        try:
            if self.F_column_labels.lower().startswith('q'):
                maptype = "q"+self.F_column_labels.lower()[1:].capitalize()
            elif self.F_column_labels.lower().startswith('k'):
                maptype = "k"+self.F_column_labels.lower()[1:].capitalize()
            else:
                maptype = self.F_column_labels.lower().capitalize()
            outprefix = re.sub(r"%s"%(maptype), "2m%s-DFc_reciprocal_space"%(maptype), self.mtz_name)
            if outprefix == self.mtz_name:
                raise AttributeError
        except AttributeError:
            outprefix = "%s_reciprocal_space"%(self.mtz_name)

        sim_annealing = ''
        if self.sim_annealing:
            sim_annealing += "simulated_annealing=True"
            for param in self.sim_annealing_pars.__dict__:
                if not param.startswith("_"):
                    sim_annealing += " simulated_annealing.%s=%s " %(param, str(self.sim_annealing_pars.__dict__[param]))
        else:
            sim_annealing += "simulated_annealing=False"
            
        additional_keywords_line = ''
        if len(self.additional_reciprocal_keywords) > 0:
            for keyword in self.additional_reciprocal_keywords:
                additional_keywords_line+= "%s " %(keyword)
            
        weight_selection_criteria = ""
        if self.weight_sel_crit.bonds_rmsd != None:
            weight_selection_criteria += "target_weights.weight_selection_criteria.bonds_rmsd=%.4f "%(self.weight_sel_crit.bonds_rmsd)
        if self.weight_sel_crit.angles_rmsd != None:
            weight_selection_criteria += "target_weights.weight_selection_criteria.angles_rmsd=%.4f "%(self.weight_sel_crit.angles_rmsd)
        if self.weight_sel_crit.r_free_minus_r_work != None:
            weight_selection_criteria += "target_weights.weight_selection_criteria.r_free_minus_r_work=%.4f "%(self.weight_sel_crit.r_free_minus_r_work)
            
        #TODO os.system -> subprocess.something (read the docs!)
        #cmd = "phenix.refine --overwrite %s %s  %s output.prefix=%s strategy=%s main.number_of_macro_cycles=%d refinement.output.write_model_cif_file=False refinement.input.xray_data.r_free_flags.disable_suitability_test=True refinement.input.xray_data.r_free_flags.ignore_pdb_hexdigest=True refinement.input.xray_data.r_free_flags.label='FreeR_flag' refinement.input.xray_data.r_free_flags.test_flag_value=1 nproc=4 wxc_scale=%f wxu_scale=%f ordered_solvent=%s write_maps=true %s %s %s" %(self.mtz_in, self.additional, self.pdb_in, outprefix, self.strategy, self.rec_cycles, self.wxc_scale, self.wxu_scale, self.solvent, self.params, weight_selection_criteria, sim_annealing)
        
        #reciprocal = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
        #p = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #cmd_output,_ = p.communicate() #get the output of the terminal
        #print(cmd_output)

        #TODO os.system -> subprocess.something (read the docs!)
        #cmd = "phenix.refine --overwrite %s %s  %s output.prefix=%s strategy=%s main.number_of_macro_cycles=%d refinement.output.write_model_cif_file=False refinement.input.xray_data.r_free_flags.disable_suitability_test=True refinement.input.xray_data.r_free_flags.ignore_pdb_hexdigest=True refinement.input.xray_data.r_free_flags.label='FreeR_flag' refinement.input.xray_data.r_free_flags.test_flag_value=1 nproc=4 wxc_scale=%f wxu_scale=%f ordered_solvent=%s write_maps=true %s %s %s" %(self.mtz_in, self.additional, self.pdb_in, outprefix, self.strategy, self.rec_cycles, self.wxc_scale, self.wxu_scale, self.solvent, self.params, weight_selection_criteria, sim_annealing)
    
        #reciprocal = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait() 
        #p = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #cmd_output,_ = p.communicate() #get the output of the terminal
        #print(cmd_output)

        reciprocal = os.system("phenix.refine --overwrite %s %s  %s output.prefix=%s strategy=%s "
                       "main.number_of_macro_cycles=%d refinement.output.write_model_cif_file=False "
                               "refinement.input.xray_data.r_free_flags.disable_suitability_test=True "
                               "refinement.input.xray_data.r_free_flags.ignore_pdb_hexdigest=True "
                               "refinement.input.xray_data.r_free_flags.label='FreeR_flag' "
                               "refinement.input.xray_data.r_free_flags.test_flag_value=1 nproc=4 wxc_scale=%f wxu_scale=%f ordered_solvent=%s write_maps=true %s %s %s %s" %(self.mtz_in, self.additional, self.pdb_in, outprefix, self.strategy, self.rec_cycles, self.wxc_scale, self.wxu_scale, self.solvent, self.params, weight_selection_criteria, sim_annealing, additional_keywords_line)) # wxc_scale=0.021 #target_weights.optimize_xyz_weight=True


        #Find output files, automatically
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
            mtz_out = "no_a_file"
            pdb_out = "refinement_did_not_finish_correcty"

        #Find output files: hardcoded appears easiest
        #mtz_out = [fle for fle in os.listdir(os.getcwd()) if outprefix in fle and fle.endswith('mtz')][0]
        #pdb_out = [fle for fle in os.listdir(os.getcwd()) if outprefix in fle and fle.endswith('pdb')][0]
        # mtz_out = "%s_001.mtz"%(outprefix)
        # pdb_out = "%s_001.pdb"%(outprefix)
        
        return mtz_out, pdb_out
    
    def get_phenix_version(self):
        """
        Weird construction to get the phenix version. This is required since some parameter names change between versions
        """
        try:
            phenix_version = int(re.search(r"phenix-1\.(.+?)\.", miller.__file__).group(1)) #This is not so robust. relies on the format being 'phenix.1.18.something' or 'phenix.1.18-something'
        except ValueError:
                phenix_version = int(re.search(r"phenix-1\.(.+?)\-", miller.__file__).group(1))
        except AttributeError:
            print('Update phenix! Verify that you are using at least Phenix.1.19.')
            phenix_version = 20 #let's assume then that the latest phenix is installed in case this fails for other reasons than a very old phenix version
        
        return phenix_version
    
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
        Some parameters have changed between version 1.17, 1.18 and 1.19 hence the weird construction to grap the version
        """       
        
        mtz_name = get_name(mtz_in)
            
        #Specify phenix version dependent parameters
        phenix_version = self.get_phenix_version()
        print("Phenix version 1.%d" %(phenix_version))
        
        if phenix_version >= 18:
            rotamer_restraints = 'rotamers.restraints.enabled=False' #rotamers.fit=all?
        else:
            rotamer_restraints = 'rotamer_restraints=False'
            
        if phenix_version >= 19:
            output_prefix = 'output.prefix=%s'%(mtz_name)
            model_format  = 'model_format=pdb'
            # outpdb        = "%s_real_space_refined_000.pdb"%(mtz_name)
            ramachandran_restraints = 'ramachandran_plot_restraints.enable=False'
        else:
            output_prefix = 'output.file_name_prefix=%s'%(mtz_name)
            model_format  = 'output.model_format=pdb'
            # outpdb        = "%s_real_space_refined.pdb"%(mtz_name)
            ramachandran_restraints = 'ramachandran_restraints=False'
            
        additional_keywords_line = ''
        if len(self.additional_real_keywords) > 0:
            for keyword in self.additional_real_keywords:
                additional_keywords_line+= "%s " %(keyword)
 
        
        real = os.system("phenix.real_space_refine %s %s %s "
                        "geometry_restraints.edits.excessive_bond_distance_limit=1000 refinement.run=minimization_global+adp scattering_table=n_gaussian c_beta_restraints=False %s refinement.macro_cycles=%d refinement.simulated_annealing=every_macro_cycle nproc=4 %s label='%s' %s %s ignore_symmetry_conflicts=True %s" %(mtz_in, self.additional, pdb_in, output_prefix, self.real_cycles, model_format, column_labels, rotamer_restraints, ramachandran_restraints, additional_keywords_line))

        #Find output file
        if real == 0 : #os.system has correctly finished. Then search for the last refined structure
            if phenix_version >= 19:
                try:
                    pdb_fles = glob.glob("%s_real_space_refined_???.pdb"%(mtz_name))
                    # [fle for fle in os.listdir(os.getcwd()) if "%s_independent_real_space_refined_0"%(mtz_name) in
                    #            fle and fle.endswith('pdb')]
                    pdb_fles.sort()
                    outpdb = pdb_fles[-1]
                except IndexError:
                    outpdb = "%s_real_space_refined_000.pdb" % (mtz_name)
            else:
                try:
                    pdb_fles = glob.glob("%s_real_space_refined.pdb"%(mtz_name))
                    # [fle for fle in os.listdir(os.getcwd()) if
                    #            "%s_independent_real_space_refined" % (mtz_name) in fle and fle.endswith('pdb')]
                    pdb_fles.sort()
                    outpdb = pdb_fles[-1]
                except IndexError:
                    outpdb = "%s_real_space_refined.pdb" % (mtz_name)
        else: # os.system has not correctly finished
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
        phenix_version = self.get_phenix_version()
        print("Phenix version 1.%d" %(phenix_version))
        
        if phenix_version >= 18:
            rotamer_restraints = 'rotamers.restraints.enabled=False' #rotamers.fit=all?
        else:
            rotamer_restraints = 'rotamer_restraints=False'
            
        if phenix_version >= 19:
            output_prefix = 'output.prefix=%s'%(ccp4_name)
            model_format  = 'model_format=pdb'
            # outpdb        = "%s_real_space_refined_000.pdb"%(ccp4_name)
            ramachandran_restraints = 'ramachandran_plot_restraints.enable=False'
        else:
            output_prefix = 'output.file_name_prefix=%s'%(ccp4_name)
            model_format  = 'output.model_format=pdb'
            # outpdb        = "%s_real_space_refined.pdb"%(ccp4_name)
            ramachandran_restraints = 'ramachandran_restraints=False'
            
        additional_keywords_line = ''
        if len(self.additional_real_keywords) > 0:
            for keyword in self.additional_real_keywords:
                additional_keywords_line+= "%s " %(keyword)
 
        
        real = os.system("phenix.real_space_refine %s %s %s "
                        "geometry_restraints.edits.excessive_bond_distance_limit=1000 refinement.run=minimization_global+adp scattering_table=n_gaussian c_beta_restraints=False %s refinement.macro_cycles=%d refinement.simulated_annealing=every_macro_cycle nproc=4 %s %s %s ignore_symmetry_conflicts=True resolution=%.2f %s" %(ccp4_in, self.additional, pdb_in, output_prefix, self.real_cycles, model_format, rotamer_restraints, ramachandran_restraints, resolution, additional_keywords_line))

        #Find output file
        if real == 0 : #os.system has correctly finished. Then search for the last refined structure
            if phenix_version >= 19:
                try:
                    pdb_fles = glob.glob("%s_real_space_refined_???.pdb"%(ccp4_name))
                    # [fle for fle in os.listdir(os.getcwd()) if "%s_independent_real_space_refined_0"%(ccp4_name) in
                    #            fle and fle.endswith('pdb')]
                    pdb_fles.sort()
                    outpdb = pdb_fles[-1]
                except IndexError:
                    outpdb = "%s_real_space_refined_000.pdb" % (ccp4_name)
            else:
                try:
                    pdb_fles = glob.glob("%s_real_space_refined.pdb"%(ccp4_name))
                    # [fle for fle in os.listdir(os.getcwd()) if
                    #            "%s_independent_real_space_refined" % (ccp4_name) in fle and fle.endswith('pdb')]
                    pdb_fles.sort()
                    outpdb = pdb_fles[-1]
                except IndexError:
                    outpdb = "%s_real_space_refined.pdb" % (ccp4_name)
        else: # os.system has not correctly finished
            outpdb = "not_a_file"

        return outpdb

        
    def get_solvent_content(self, pdb_in):
        """
        Extract solvent content
        """
        pdb_hier = iotbx.pdb.hierarchy.input(file_name=pdb_in)
        hier = pdb_hier.hierarchy
        overall_counts = hier.overall_counts()
        
        pdb_ini = iotbx.pdb.input(pdb_in)
        xray_structure = pdb_ini.xray_structure_simple()
        
        vm_calc = p_vm_calculator(xray_structure.crystal_symmetry(),
            n_residues=overall_counts.resname_classes.get("common_amino_acid", 0))
            #can be estended with n_bases=.overall_countsresname_classes.get("common_rna_dna", 0)
        return vm_calc.solc(vm=vm_calc.vm(copies=1))
    
    def generate_map_params_bsharpening(self):
        param_file = "map.params"
        o = open(param_file,"w")
        o.write('refinement {\n\
  electron_density_maps {\n\
    map_coefficients {\n\
      map_type = "2mFo-DFc"\n\
      mtz_label_amplitudes = "2FOFCWT"\n\
      mtz_label_phases = "PH2FOFCWT"\n\
      fill_missing_f_obs = True\n\
      sharpening = True\n\
      sharpening_b_factor = None\n\
    }\n\
    map_coefficients {\n\
      map_type = "2mFo-DFc"\n\
      mtz_label_amplitudes = "2FOFCWT_no_fill"\n\
      mtz_label_phases = "PH2FOFCWT_no_fill"\n\
      fill_missing_f_obs = False\n\
      sharpening = True\n\
      sharpening_b_factor = None\n\
    }\n\
    map_coefficients {\n\
      map_type = "mFo-DFc"\n\
      mtz_label_amplitudes = "FOFCWT"\n\
      mtz_label_phases = "PHFOFCWT"\n\
      fill_missing_f_obs = False\n\
      sharpening = False\n\
      sharpening_b_factor = None\n\
    }\n\
    map {\n\
      map_type = "2mFo-DFc"\n\
      fill_missing_f_obs = True\n\
      sharpening = True\n\
      format = ccp4\n\
      sharpening_b_factor = None\n\
    }\n\
    map {\n\
      map_type = "2mFo-DFc"\n\
      format = ccp4\n\
      fill_missing_f_obs = False\n\
      sharpening = True\n\
      sharpening_b_factor = None\n\
    }\n\
    map {\n\
      format = ccp4\n\
      map_type = "mFo-DFc"\n\
      fill_missing_f_obs = False\n\
      sharpening = False\n\
      sharpening_b_factor = None\n\
    }\n\
  }\n\
}')
        o.close()
        
        return param_file    

    
    def phenix_density_modification(self, mtz_in, pdb_in):
        """
        Run phenix.density_modification
        """
        
        solc = self.get_solvent_content(pdb_in)
        outname = re.sub(r".mtz$", "_densitymod.mtz", mtz_in)
        log_file = re.sub(r".mtz$", ".log", outname)
        
        print('Running density modification, output written to %s. Please wait...'%(log_file))
        os.system("phenix.density_modification %s %s input_files.map_coeffs_file=%s solvent_content=%.3f "
                  "denmod.mask_type=histograms output_files.output_mtz=%s > %s" %(mtz_in, pdb_in, mtz_in, solc, outname, log_file))
        
        return outname

    def write_refmac_for_dm(self, pdb_in):
        """
        Write and excecute a bash script to run Refmac in order to get an mtz-file that can be used by dm.
        """

        additional_lines = ''
        for cif in self.additional.split():
            if cif.endswith(".cif"):
                additional_lines += 'LIB_IN %s ' % (cif)

        mtz_out = "%s_for_dm.mtz" % (get_name(self.mtz_in))
        pdb_out = "%s_for_dm.pdb" % (get_name(self.mtz_in))
        log_file = "%s_refmac_for_dm.log" % (get_name(self.mtz_in))

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
eof" % (pdb_in, self.mtz_in, pdb_out, mtz_out, additional_lines, log_file, self.F_column_labels,
        self.F_column_labels))

        i.close()
        os.system("chmod +x %s" % (script_out))
        print('Running refmac with zero cycles to prepare files suitable for running dm afterwards, output written to '
              '%s. Please wait...' % (log_file))
        os.system("./%s" % (script_out))

        return mtz_out, pdb_out

    def write_density_modification_script(self, mtz_in, pdb_in, mtz_out, combine, cycles, log_file):
        """
        Write script to perform density modification with dm. In order to have the correct columns, the mtz-file
        should original from refmac.
        """
        solc = self.get_solvent_content(pdb_in)

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
eor\n' % (mtz_in, mtz_out, log_file, solc, combine, cycles, self.F_column_labels, self.F_column_labels))

        ccp4_map_name = re.sub(r".mtz$", ".ccp4", mtz_out)

        i.write('#generate map in ccp4 format\n\
fft hklin %s mapout %s <<eof > fft.log\n\
LABI F1=FDM PHI=PHIDM\n\
eof' % (mtz_out, ccp4_map_name))

        i.close()
        os.system("chmod +x %s" % (script_out))
        return script_out

    def ccp4_dm(self, pdb_in, combine,cycles):
        """
        Use dm to perform density modification. In order to get the correct columns, we first run refmac with 0
        cycles with the refined model from phenix.refine.
        """
        mtz_for_dm, pdb_for_dm = self.write_refmac_for_dm(pdb_in)
        mtz_out_dm = re.sub(r"for_dm.mtz$", "dm.mtz", mtz_for_dm)
        if (os.path.isfile(mtz_for_dm) and os.path.isfile(pdb_for_dm)):
            log_file = re.sub(r".mtz$", ".log", mtz_out_dm)
            script_dm = self.write_density_modification_script(mtz_for_dm, pdb_for_dm, mtz_out_dm, combine, cycles, log_file)
            print('Running density modification, output written to %s. Please wait...' % (log_file))
            os.system("./%s" % (script_dm))

        return mtz_out_dm
