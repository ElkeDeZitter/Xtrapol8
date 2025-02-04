# -*- coding: utf-8 -*-
"""
write bash script to run refmac and coot for reciprocal and real space refinement.

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
import iotbx.pdb
import iotbx.file_reader
from mmtbx.scaling.matthews import p_vm_calculator
from libtbx import adopt_init_args
from Fextr_utils import get_name
from iotbx.file_reader import any_file
from cctbx import miller

class Refmac_refinement(object):
    def __init__(self,
                 mtz_in,
                 pdb_in,
                 additional,
                 F_column_labels = 'QFEXTR',
                 rfree_col       = 'FreeR_flag',
                 fill_missing    = False,
                 add_rfree       = False,
                 refinement_weight         = 'AUTO',
                 refinement_weight_sigmas  = 'NOEX',
                 refinement_weighting_term = 0.2,
                 refinement_type       = 'REST',
                 TLS                   = False,
                 TLS_cycles            = 20,
                 bfac_set              = 30,
                 twinning              = False,
                 Brefinement           = 'ISOT',
                 cycles                = 5,
                 external_restraints   = None,
                 jelly_body_refinement = False,
                 jelly_body_sigma      = 0.03,
                 jelly_body_additional_restraints = None,
                 map_sharpening        = False,
                 additional_reciprocal_keywords = []):
        
        #Try this instead of huge repetition of the arguments
        adopt_init_args(self, locals())
                
        self.mtz_name = get_name(self.mtz_in)

    def add_missing_reflections(self):
        #add missing reflections with ccp4 and add R-free
        new_mtz = '%s_rfree.mtz' %(self.mtz_name)
        extra_line = 'uniqueify %s %s' %(self.mtz_in, new_mtz)
        self.rfree_col = 'FreeR_flag'
        self.mtz_in = new_mtz
        return extra_line

    def add_rfree_col(self):
        new_mtz = '%s_rfree.mtz' %(self.mtz_name)
        extra_line = 'freerflag hklin %s hklout %s << eof\nEND\neof' %(self.mtz_in, new_mtz)
        self.rfree_col = 'FreeR_flag'
        self.mtz_in = new_mtz
        return extra_line

    def write_refmac_input_reciprocal_space_refinement(self, mtz_out, pdb_out, log_file):
        """
        Prepare bash file to run refmac.
        BFAC SET <value> default 30
        """
        if self.fill_missing == True:
            extra_line = self.add_missing_reflections()
        elif self.add_rfree == True:
            extra_line = self.add_rfree_col()
        else:
            extra_line = ''
            
        additional_lines = ''
        for cif in self.additional.split():
            additional_lines+='LIB_IN %s '%(cif)
            
        if self.refinement_weight == 'AUTO':
            refinement_weight = '%s' %(self.refinement_weight)
        else:
            refinement_weight = '%s %s %f' %(self.refinement_weight_sigmas, self.refinement_weight, self.refinement_weighting_term)
            
        if self.TLS:
            TLS_line ='REFI TLSC %d\nBFAC SET %f' %(self.TLS_cycles, self.bfac_set)
        else:
            TLS_line = ''
            
        if self.twinning:
            twin_line = 'TWIN'
        else:
            twin_line = ''
            
        ext_restraints_line = ''
        for restraint in self.external_restraints:
            ext_restraints_line+= 'external %s\n' %(restraint)
        
        if self.jelly_body_refinement:
            if self.cycles < 25:
                self.cycles = self.cycles*5
            jelly_body_line = 'RIDG DIST %f\n' %(self.jelly_body_sigma)
            for restraint in self.jelly_body_additional_restraints:
                jelly_body_line+='RIDG %s\n' %(restraint)
        else:
            jelly_body_line = ''
            
        if self.map_sharpening:
            map_sharp_line = 'MAPC SHAR'
        else:
            map_sharp_line = ''
            
        additional_keywords_line = ''
        if len(self.additional_reciprocal_keywords) > 0:
            for keyword in self.additional_reciprocal_keywords:
                additional_keywords_line+= "%s " %(keyword)
        
        script_out = 'launch_refmac.sh'
        i = open(script_out,'w')
        i.write('#!/bin/sh \n\
%s\n\
#refmac:\n\
refmac5 HKLIN %s HKLOUT %s XYZIN %s XYZOUT %s %s << eor > %s\n\
MAKE HYDR No\n\
MAKE CHEC NONE\n\
MAKE SS Yes\n\
MAKE SYMM Yes\n\
MAKE SUGA Yes\n\
MAKE CONN No\n\
LABIN FP=%s SIGFP=SIG%s FREE=%s\n\
LABO FC=FC PHIC=PHIC FWT=2FOFCWT PHWT=PH2FOFCWT DELFWT=FOFCWT PHDELWT=PHFOFCWT\n\
WEIGHT %s\n\
%s\n\
REFI TYPE %s\n\
REFI RESI MLKF\n\
REFI BREF %s\n\
%s\
SCALe TYPE BULK\n\
NCYC %d\n\
MONI FEW\n\
NOHARVEST\n\
%s\
%s\
%s\n\
%s\n\
end\n\
eor\n' %(extra_line, self.mtz_in, mtz_out, self.pdb_in, pdb_out, additional_lines, log_file, self.F_column_labels, self.F_column_labels, self.rfree_col, refinement_weight, TLS_line, self.refinement_type, self.Brefinement, twin_line, self.cycles, ext_restraints_line, jelly_body_line, map_sharp_line, additional_keywords_line))
        
        ccp4_map_name = re.sub(r".mtz$", "_2mFo-DFc_filled.ccp4", mtz_out)
        
        i.write('#generate 2mFo-DFc in ccp4 format\n\
fft HKLIN %s MAPOUT %s << eof > fft.log\n\
LABIN F1=2FOFCWT PHI=PH2FOFCWT\n\
FILLIN\n\
end\n\
eof\n'%(mtz_out, ccp4_map_name))
        
        ccp4_diff_map_name = re.sub(r".mtz$", "_mFo-DFc.ccp4", mtz_out)
        
        i.write('#generate mFo-DFc in ccp4 format\n\
fft HKLIN %s MAPOUT %s << eof > fft.log\n\
LABIN F1=FOFCWT PHI=PHFOFCWT\n\
end\n\
eof\n'%(mtz_out, ccp4_diff_map_name))
        
        i.close()
        os.system("chmod +x %s" %(script_out))
        return script_out
    
    def get_solvent_content(self):
        """
        Extract solvent content
        """
        pdb_hier = iotbx.pdb.hierarchy.input(file_name=self.pdb_in)
        hier = pdb_hier.hierarchy
        overall_counts = hier.overall_counts()
        
        pdb_ini = iotbx.pdb.input(self.pdb_in)
        xray_structure = pdb_ini.xray_structure_simple()
        
        #scattering table should not be specified here because the structure is onluy used
        #to calculate the solvent content. No usage to calculate f_model
        
        vm_calc = p_vm_calculator(xray_structure.crystal_symmetry(),
            n_residues=overall_counts.resname_classes.get("common_amino_acid", 0))
            #can be estended with n_bases=.overall_countsresname_classes.get("common_rna_dna", 0)
        return vm_calc.solc(vm=vm_calc.vm(copies=1))
        
    def write_density_modification_script(self, mtz_in, mtz_out, combine, cycles, log_file):
        """
        Write script to perform density modification with dm
        """
        solc = self.get_solvent_content()
        
        script_out = 'launch_dm.sh'
        i = open(script_out,'w')
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
eor\n' %(mtz_in, mtz_out, log_file, solc, combine, cycles, self.F_column_labels, self.F_column_labels))
      
        ccp4_map_name = re.sub(r".mtz$", ".ccp4", mtz_out)
      
        i.write('#generate map in ccp4 format\n\
fft hklin %s mapout %s <<eof > fft.log\n\
LABI F1=FDM PHI=PHIDM\n\
eof' %(mtz_out, ccp4_map_name))

        i.close()
        os.system("chmod +x %s" %(script_out))
        return script_out        
        
    def reciprocal_space_refinement(self):
        try:
            if self.F_column_labels.lower().startswith('q'):
                maptype = "q"+self.F_column_labels.lower()[1:].capitalize()
            else:
                maptype = self.F_column_labels.lower().capitalize()
            mtz_out = re.sub(r"%s"%(maptype), "2m%s-DFc_refmac_reciprocal_space_001.mtz"%(maptype), self.mtz_name)
            pdb_out = re.sub(r"%s"%(maptype), "2m%s-DFc_refmac_reciprocal_space_001.pdb"%(maptype), self.mtz_name)
            log_file = re.sub(r"%s"%(maptype), "2m%s-DFc_refmac_reciprocal_space_001.log"%(maptype), self.mtz_name)
            if mtz_out == self.mtz_name:
                raise AttributeError
        except AttributeError:
            mtz_out = '%s_refmac_reciprocal_space_001.mtz' %(self.mtz_name)
            pdb_out = '%s_refmac_reciprocal_space_001.pdb' %(self.mtz_name)
            log_file = '%s_refmac_reciprocal_space_001.log' %(self.mtz_name)
            
        script_refmac = self.write_refmac_input_reciprocal_space_refinement(mtz_out, pdb_out, log_file)
        print('Running Refmac, output written to %s. Please wait...'%(log_file))
        reciprocal = os.system("./%s" %(script_refmac))
        
        if reciprocal != 0:
            mtz_out = "not_a_file"
            pdb_out = "refinement_did_not_finish_correcty"

        
        # if os.path.isfile(pdb_out) == False:
        #     print("pdb file from refinement not found. Refinement probably failed")
        #     pdb_out = self.pdb_in
        # if os.path.isfile(mtz_out) == False:
        #     print("mtz file from refinement not found. Refinement probably failed")
        #     mtz_out = self.mtz_in
            
        return mtz_out, pdb_out
    
    def ccp4_dm(self, pdb_in, combine, cycles):
        """
        Use dm to perform density modification
        """
        mtz_for_dm = re.sub(r".pdb$", ".mtz", pdb_in)
        mtz_out_dm = re.sub(r".pdb$", "_dm.mtz", pdb_in)
        if os.path.isfile(mtz_for_dm):
            log_file = re.sub(r".mtz$", ".log", mtz_out_dm)
            script_dm = self.write_density_modification_script(mtz_for_dm, mtz_out_dm, combine, cycles, log_file)
            print('Running density modification, output written to %s. Please wait...'%(log_file))
            os.system("./%s" %(script_dm))
        
        return mtz_out_dm


class Coot_refinement(object):
    def __init__(self,
                 ligands,
                 additional = ''):
        
        self.ligands    = ligands
        self.additional = additional
    
    def ligands_refinement(self, pdb_in):
        pdb_hier = iotbx.pdb.hierarchy.input(file_name=pdb_in)
        hier = pdb_hier.hierarchy
        
        if len(self.ligands) == 0:
            ligand_refinement = ""
        else:
            ligand_refinement="set_ligand_cluster_sigma_level(0.8)\nset_ligand_flexible_ligand_n_samples(10)\nset_matrix(120)\n" 
            for lig in self.ligands:
                for chain in hier.chains():
                    if chain.is_protein():
                        for res_group in chain.residue_groups():
                            for atom_group in res_group.atom_groups():
                                if atom_group.resname==lig:
                                    ligand_refinement+=('refine_zone(0,"%s",%d,%d,"%s")\naccept_regularizement()\n'%(chain.id, res_group.resseq_as_int(), res_group.resseq_as_int(), atom_group.altloc))
        return ligand_refinement
    
    def check_mtz_column(self, mtz_in, column_labels):
        """
        Check if the column is present in the mtz file, use the column labels. e.g. '2FOFCWT,PH2FOFCWT'
        """
        column_found = False
        hkl = any_file(mtz_in,force_type="hkl", raise_sorry_if_errors=False)
        for array in hkl.file_object.as_miller_arrays():
            if array.info().label_string() == column_labels:
                column_found = True
                
        return column_found
    
    def suggest_mtz_column(self, mtz_in, column_labels, column_label_strings_constraint=''):
        """
        Find the closest column labels in an mtz file.
        Additional constraints could be added to the string, e.g. should contain "2F" if searching for the 2FoFc type
        """
        from difflib import get_close_matches
        
        column_label_strings = []
        hkl = any_file(mtz_in,force_type="hkl", raise_sorry_if_errors=False)
        for array in hkl.file_object.as_miller_arrays():
            column_label_strings.append(array.info().label_string())
        
        column_suggestions = get_close_matches(column_labels, column_label_strings)
        column_suggestions.sort()
        
        close_column_found = False
        if len(column_suggestions) >=2:
            for suggestion in column_suggestions:
                if column_label_strings_constraint in suggestion:
                    close_column_found = True
                    column_labels_new = suggestion
                    break
                
        if close_column_found == False:
            column_labels_new = "2FOFCWT,PH2FOFCWT"
            
        return column_labels_new.split(",")

    def write_coot_input_real_space_refinement_mtz(self, pdb_in, mtz_in, column_labels):
        """
        Write input script for COOT based on the usage of an mtz file
        """
        mtz_name = get_name(mtz_in)
        #if "/" in mtz_in:
            #mtz_name = re.search(r"\/(.+?)\.mtz", mtz_in).group(1).split("/")[-1]
        #else:
            #mtz_name = re.sub("\.mtz","",mtz_in)

        column_labels_0_F = column_labels.split(",")[0].lstrip().rstrip()
        column_labels_0_P = column_labels.split(",")[1].lstrip().rstrip()
        column_labels_1_F = column_labels.split(",")[2].lstrip().rstrip()
        column_labels_1_P = column_labels.split(",")[3].lstrip().rstrip()
        
        if self.check_mtz_column(mtz_in, "{:s},{:s}".format(column_labels_0_F,column_labels_0_P)) == False:
            print("{:s}: column labels {:s},{:s} not found".format(mtz_in, column_labels_0_F,column_labels_0_P))
            column_labels_0_F,column_labels_0_P = self.suggest_mtz_column(mtz_in, "{:s},{:s}".format(column_labels_0_F,column_labels_0_P), column_label_strings_constraint="2")
            print("{:s}: columns {:s},{:s} will be used".format(mtz_in, column_labels_0_F, column_labels_0_P))

        if self.check_mtz_column(mtz_in, "{:s},{:s}".format(column_labels_1_F,column_labels_1_P)) == False:
            print("{:s}: column labels {:s},{:s} not found".format(mtz_in, column_labels_1_F,column_labels_1_P))
            column_labels_1_F,column_labels_1_P = self.suggest_mtz_column(mtz_in, "{:s},{:s}".format(column_labels_1_F,column_labels_1_P))
            print("{:s}: columns {:s},{:s} will be used".format(mtz_in, column_labels_1_F, column_labels_1_P))
        
        additional_lines = ""
        for cif in self.additional.split():
            if cif.endswith(".cif"):
                additional_lines+='read_cif_dictionary("%s")\n'%(cif)

        pdb_out = "%s_coot_real_space_refined_000.pdb"%(mtz_name)
        
        script_root = 'coot_real_space_refinement'
        script_out = '%s.py' %(script_root)
        j = 1
        while os.path.exists(script_out):
            script_root = "%s_%d" %(script_root, j)
            script_out = '%s.py' %(script_root)
            j += 1
            if j == 100: #to avoid endless loop
                break
        
        i = open(script_out,'w')
        i.write('handle_read_draw_molecule("%s")\n\
%s\
set_auto_read_column_labels("%s","%s",0)\n\
set_auto_read_column_labels("%s","%s",1)\n\
auto_read_make_and_draw_maps_from_mtz("%s")\n\
close_molecule(3)\n\
close_molecule(4)\n\
set_refine_ramachandran_angles(1)\n\
set_environment_distances_distance_limits(2.4,3.6)\n\
set_ligand_water_to_protein_distance_limits(2.4,3.6)\n\
set_refinement_immediate_replacement(1)\n\
%s\
set_matrix(20)\n\
stepped_refine_protein(0,10)\n\
accept_regularizement()\n\
fit_waters(0)\n\
accept_regularizement()\n\
write_pdb_file(0,"%s")\n\
coot_no_state_real_exit(1)' %(pdb_in, additional_lines, column_labels_0_F, column_labels_0_P, column_labels_1_F, column_labels_1_P, mtz_in, self.ligands_refinement(pdb_in), pdb_out))   
        
        i.close()
                
        #save_state_file_py("%s_coot_real_space_refined.py")\n\   
                
        return script_out, pdb_out
    
    def write_coot_input_real_space_refinement_ccp4(self, pdb_in, ccp4_in):
        """
        Write input script for COOT based on the usage of a ccp4 file
        """
        ccp4_name = get_name(ccp4_in)
        
        additional_lines = ""
        for cif in self.additional.split():
            if cif.endswith(".cif"):
                additional_lines+='read_cif_dictionary("%s")\n'%(cif)
                
                
        pdb_out = "%s_coot_real_space_refined_000.pdb"%(ccp4_name)
        
        script_root = 'coot_real_space_refinement'
        script_out = '%s.py' %(script_root)
        j = 1
        while os.path.exists(script_out):
            script_root = "%s_%d" %(script_root, j)
            script_out = '%s.py' %(script_root)
            j += 1
            if j == 100: #to avoid endless loop
                break
            
        i = open(script_out,'w')
        i.write('handle_read_draw_molecule("%s")\n\
%s\
handle_read_ccp4_map("%s", 0)\n\
set_refine_ramachandran_angles(1)\n\
set_environment_distances_distance_limits(2.4,3.6)\n\
set_ligand_water_to_protein_distance_limits(2.4,3.6)\n\
set_refinement_immediate_replacement(1)\n\
%s\
set_matrix(20)\n\
stepped_refine_protein(0,10)\n\
accept_regularizement()\n\
fit_waters(0)\n\
accept_regularizement()\n\
write_pdb_file(0,"%s")\n\
coot_no_state_real_exit(1)' %(pdb_in, additional_lines, ccp4_in, self.ligands_refinement(pdb_in), pdb_out))        
        
        i.close()
                
        #save_state_file_py("%s_coot_real_space_refined.py")\n\        
                
        return script_out, pdb_out
    
    def get_mtz_resolution(self, mtz_in):
        """
        Easy extraction of the resolution boundaries of an mtz file
        """
        reflections = any_file(mtz_in, force_type="hkl", raise_sorry_if_errors=True)
        low_res, high_res = reflections.file_content.file_content().max_min_resolution()
         
        return low_res, high_res
 
        
    def real_space_refinement_mtz(self,  mtz_in, pdb_in, column_labels):
        """
        Real space refinement within COOT based on an mtz file and spefic columns
        """
        coot_log = '%s_coot_real_space_refined_000.log' %(get_name(mtz_in))
        script_coot, pdb_out = self.write_coot_input_real_space_refinement_mtz(pdb_in, mtz_in, column_labels)
        print('Running Real space refinement in COOT, output written to %s. Please wait...' %(coot_log))
        os.system("coot --no-graphics --script %s > %s" %(script_coot, coot_log)) #Use a second coot install to avoid interference with other coot windows that might already be open
        
        return pdb_out
    
    def real_space_refinement_ccp4(self, ccp4_in, pdb_in, resolution):
        """
        Real space refinement within COOT based on a ccp4 file. No ambiguitiy stemming from multiple different columns
        
        Resolution is not required, but passed as an argument to remain consistent with real_space_refinement_ccp4 function within the Phenix_real_space_refinement class
        """
        coot_log = '%s_coot_real_space_refined_000.log' %(get_name(ccp4_in))
        script_coot, pdb_out = self.write_coot_input_real_space_refinement_ccp4(pdb_in, ccp4_in)
        print('Running Real space refinement in COOT, output written to %s. Please wait...' %(coot_log))
        os.system("coot --no-graphics --script %s > %s" %(script_coot, coot_log)) #Use a second coot install to avoid interference with other coot windows that might already be open
        
        return pdb_out
