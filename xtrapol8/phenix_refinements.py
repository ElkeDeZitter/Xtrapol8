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

import glob
import os
import re
import shlex
import subprocess
import sys
from difflib import get_close_matches

import iotbx.pdb
from iotbx.file_reader import any_file
from libtbx import adopt_init_args
from mmtbx.scaling.matthews import p_vm_calculator

from .Fextr_utils import get_name, get_phenix_version, get_phenix_version_number
from .programs.dm import dm
from .programs.fft import fft
from .programs.refmac import refmac_for_dm


class Phenix_reciprocal_space_refinement:
    def __init__(
        self,
        mtz_in,
        pdb_in,
        additional="",
        F_column_labels="QFEXTR",
        strategy=["individual_sites", "individual_adp"],
        rec_cycles=5,
        wxc_scale=0.5,
        wxu_scale=1.0,
        solvent=False,
        sim_annealing=False,
        sim_annealing_pars={},
        map_sharpening=False,
        scattering_table="n_gaussian",
        weight_sel_crit={},
        additional_reciprocal_keywords=[],
        log=sys.stdout,
    ):
        adopt_init_args(self, locals())

        self.strategy = "+".join(strategy)

        if map_sharpening:
            self.params = self.generate_map_params_bsharpening()
        else:
            self.params = ""

        self.mtz_name = get_name(self.mtz_in)
        
        self.phenix_version = get_phenix_version_number(get_phenix_version())

    def reciprocal_space_refinement(self):
        """
        run phenix.refine
        Support for version 1.18 and higher
        """
        try:
            if self.F_column_labels.lower().startswith("q"):
                maptype = "q" + self.F_column_labels.lower()[1:].capitalize()
            elif self.F_column_labels.lower().startswith("k"):
                maptype = "k" + self.F_column_labels.lower()[1:].capitalize()
            else:
                maptype = self.F_column_labels.lower().capitalize()
            outprefix = re.sub(
                r"%s" % (maptype),
                "2m%s-DFc_phenix_reciprocal_space" % (maptype),
                self.mtz_name,
            )
            if outprefix == self.mtz_name:
                raise AttributeError
        except AttributeError:
            outprefix = "%s_reciprocal_space" % (self.mtz_name)

        sim_annealing = ""
        if self.sim_annealing:
            sim_annealing += "simulated_annealing=True"
            for param in self.sim_annealing_pars.__dict__:
                if not param.startswith("_"):
                    sim_annealing += " simulated_annealing.%s=%s " % (
                        param,
                        str(self.sim_annealing_pars.__dict__[param]),
                    )
        else:
            sim_annealing += "simulated_annealing=False"

        additional_keywords_line = ""
        if len(self.additional_reciprocal_keywords) > 0:
            for keyword in self.additional_reciprocal_keywords:
                additional_keywords_line += "%s " % (keyword)

        weight_selection_criteria = ""
        if self.weight_sel_crit.bonds_rmsd != None:
            weight_selection_criteria += (
                "target_weights.weight_selection_criteria.bonds_rmsd=%.4f "
                % (self.weight_sel_crit.bonds_rmsd)
            )
        if self.weight_sel_crit.angles_rmsd != None:
            weight_selection_criteria += (
                "target_weights.weight_selection_criteria.angles_rmsd=%.4f "
                % (self.weight_sel_crit.angles_rmsd)
            )
        if self.weight_sel_crit.r_free_minus_r_work != None:
            weight_selection_criteria += (
                "target_weights.weight_selection_criteria.r_free_minus_r_work=%.4f "
                % (self.weight_sel_crit.r_free_minus_r_work)
            )

        if self.phenix_version <= 1.20:
            r_free_flag_parameters = "refinement.input.xray_data.r_free_flags.disable_suitability_test=True refinement.input.xray_data.r_free_flags.ignore_pdb_hexdigest=True refinement.input.xray_data.r_free_flags.label='FreeR_flag' refinement.input.xray_data.r_free_flags.test_flag_value=1"
        else:  # phenix 1.21 and higher
            r_free_flag_parameters = "data_manager.fmodel.xray_data.r_free_flags.ignore_pdb_hexdigest=True data_manager.fmodel.xray_data.r_free_flags.test_flag_value=1"
        # data_manager.fmodel.xray_data.r_free_flags.disable_suitability_test=True
        # Disable_suitability_test cannot be done. It keeps on giving an error message about the label and value. All combination have been tested, it seems that this does not work

        command = (
            "phenix.refine --overwrite %s %s  %s output.prefix=%s strategy=%s "
            "main.number_of_macro_cycles=%d refinement.output.write_model_cif_file=False "
            "refinement.main.scattering_table=%s "
            "%s refinement.main.nproc=4 wxc_scale=%f wxu_scale=%f ordered_solvent=%s write_maps=true %s %s %s %s"
            % (
                self.mtz_in,
                self.additional,
                self.pdb_in,
                outprefix,
                self.strategy,
                self.rec_cycles,
                self.scattering_table,
                r_free_flag_parameters,
                self.wxc_scale,
                self.wxu_scale,
                self.solvent,
                self.params,
                weight_selection_criteria,
                sim_annealing,
                additional_keywords_line,
            )
        )
        reciprocal = subprocess.call(shlex.split(command))

        # Find output files, automatically
        if reciprocal == 0:  # correctly finished, search for the last refined structure
            try:
                mtz_fles = glob.glob("%s_???.mtz" % (outprefix))
                # [fle for fle in os.listdir(os.getcwd()) if outprefix+"_0" in fle and fle.endswith('mtz') and not
                # fle.endswith("_densitymod.mtz")]
                mtz_fles.sort()
                mtz_out = mtz_fles[-1]
                pdb_fles = glob.glob("%s_???.pdb" % (outprefix))
                # [fle for fle in os.listdir(os.getcwd()) if outprefix+"_0" in fle and fle.endswith('pdb')]
                pdb_fles.sort()
                pdb_out = pdb_fles[-1]
            except IndexError:
                mtz_out = "%s_001.mtz" % (outprefix)
                pdb_out = "%s_001.pdb" % (outprefix)
        else:  # not correctly finished
            mtz_out = "not_a_file"
            pdb_out = "refinement_did_not_finish_correcty"

        return mtz_out, pdb_out

    def get_solvent_content(self):
        """
        Extract solvent content
        """
        pdb_hier = iotbx.pdb.hierarchy.input(file_name=self.pdb_in)
        hier = pdb_hier.hierarchy
        overall_counts = hier.overall_counts()

        pdb_ini = iotbx.pdb.input(self.pdb_in)
        xray_structure = pdb_ini.xray_structure_simple()

        # scattering table should not be specified here because the structure is onluy used
        # to calculate the solvent content. No usage to calculate f_model

        vm_calc = p_vm_calculator(
            xray_structure.crystal_symmetry(),
            n_residues=overall_counts.resname_classes.get("common_amino_acid", 0),
        )
        # can be estended with n_bases=.overall_countsresname_classes.get("common_rna_dna", 0)
        return vm_calc.solc(vm=vm_calc.vm(copies=1))

    def generate_map_params_bsharpening(self):
        param_file = "map.params"
        o = open(param_file, "w")
        o.write(
            'refinement {\n\
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
}'
        )
        o.close()

        return param_file

    def do_refmac_for_dm(self, pdb_in):
        "Run Refmac in order to get an mtz-file that can be used by dm."
        return refmac_for_dm(
            xyzin=pdb_in,
            hklin=self.mtz_in,
            libins=self.additional.split(),
            f_label=self.F_column_labels,
        )

    def do_density_modification(self, mtz_in, mtz_out, combine, cycles, log_file):
        """
        Perform density modification with dm.
        In order to have the correct columns, the mtz-file
        should original from refmac.
        """
        solc = self.get_solvent_content()
        dm(mtz_in, mtz_out, solc, combine, cycles, self.F_column_labels, log_file)
        ccp4_map_name = re.sub(r".mtz$", ".ccp4", mtz_out)
        fft(mtz_out, ccp4_map_name, "FDM", "PHIDM")

    def ccp4_dm(self, pdb_in, combine, cycles):
        """
        Use dm to perform density modification.
        In order to get the correct columns, we first run refmac with 0 cycles
        with the refined model from phenix.refine.
        """
        mtz_for_dm, _ = self.do_refmac_for_dm(pdb_in)
        mtz_out_dm = re.sub(r".pdb$", "_dm.mtz", pdb_in)
        if os.path.isfile(mtz_for_dm):
            log_file = re.sub(r".mtz$", ".log", mtz_out_dm)
            print(f"Running density modification, see {log_file}. Please wait...")
            self.do_density_modification(
                mtz_for_dm, mtz_out_dm, combine, cycles, log_file
            )
        return mtz_out_dm


class Phenix_real_space_refinement:
    def __init__(
        self,
        real_cycles=5,
        additional="",
        additional_real_keywords=[],
        scattering_table="n_gaussian",
        log=sys.stdout,
    ):
        adopt_init_args(self, locals())
        
        self.phenix_version = get_phenix_version_number(get_phenix_version())

    def check_mtz_column(self, mtz_in, column_labels):
        """
        Check if the column is present in the mtz file, use the column labels. e.g. '2FOFCWT,PH2FOFCWT'
        """
        column_found = False
        hkl = any_file(mtz_in, force_type="hkl", raise_sorry_if_errors=False)
        for array in hkl.file_object.as_miller_arrays():
            if array.info().label_string() == column_labels:
                column_found = True

        return column_found

    def suggest_mtz_column(
        self, mtz_in, column_labels, column_label_strings_constraint=""
    ):
        """
        Find the closest column labels in an mtz file.
        Additional constraints could be added to the string, e.g. should contain "2F" if searching for the 2FoFc type
        """
        column_label_strings = []
        hkl = any_file(mtz_in, force_type="hkl", raise_sorry_if_errors=False)
        for array in hkl.file_object.as_miller_arrays():
            column_label_strings.append(array.info().label_string())

        column_suggestions = get_close_matches(column_labels, column_label_strings)
        column_suggestions.sort()

        close_column_found = False
        if len(column_suggestions) >= 2:
            for suggestion in column_suggestions:
                if column_label_strings_constraint in suggestion:
                    close_column_found = True
                    column_labels_new = suggestion
                    break

        if close_column_found == False:
            column_labels_new = "2FOFCWT,PH2FOFCWT"

        return column_labels_new

    def get_mtz_resolution(self, mtz_in):
        """
        Easy extraction of the resolution boundaries of an mtz file
        """
        reflections = any_file(mtz_in, force_type="hkl", raise_sorry_if_errors=True)
        low_res, high_res = reflections.file_content.file_content().max_min_resolution()

        return low_res, high_res

    def real_space_refinement_mtz(self, mtz_in, pdb_in, column_labels):
        """
        Real space refinement based on mtz file and specified column labels
        run phenix.real_space_refine
        Support for phenix 1.18 and higher
        """

        mtz_name = get_name(mtz_in)

        # check if the expected columns in the mtz file can be found
        if self.check_mtz_column(mtz_in, column_labels) == False:
            print("{:s}: column labels {:s} not found".format(mtz_in, column_labels))
            column_labels = self.suggest_mtz_column(
                mtz_in, column_labels, column_label_strings_constraint="2"
            )
            print("{:s}: columns {:s} will be used".format(mtz_in, column_labels))

        rotamer_restraints = "rotamers.restraints.enabled=False"  # rotamers.fit=all?

        prefix = "%s_phenix" % (mtz_name)
        if self.phenix_version <= 1.18:
            output_prefix = 'output.file_name_prefix=%s'% (prefix)
            model_format  = 'output.model_format=pdb'
            # outpdb        = "%s_real_space_refined.pdb"%(mtz_name)
            ramachandran_restraints = 'ramachandran_restraints=False'
        else: # phenix version 1.19 and higher
            output_prefix = "output.prefix=%s" % (prefix)
            model_format = "model_format=pdb"
            ramachandran_restraints = "ramachandran_plot_restraints.enable=False"

        additional_keywords_line = ""
        if len(self.additional_real_keywords) > 0:
            for keyword in self.additional_real_keywords:
                additional_keywords_line += "%s " % (keyword)

        command = (
            "phenix.real_space_refine %s %s %s "
            "geometry_restraints.edits.excessive_bond_distance_limit=1000 refinement.run=minimization_global+adp scattering_table=%s c_beta_restraints=False %s refinement.macro_cycles=%d refinement.simulated_annealing=every_macro_cycle nproc=4 %s label='%s' %s %s ignore_symmetry_conflicts=True %s"
            % (
                mtz_in,
                self.additional,
                pdb_in,
                self.scattering_table,
                output_prefix,
                self.real_cycles,
                model_format,
                column_labels,
                rotamer_restraints,
                ramachandran_restraints,
                additional_keywords_line,
            )
        )
        real = subprocess.call(shlex.split(command))

        # Find output file
        if real == 0:  # correctly finished. search for the last refined structure
            pdb_pattern = (
                "%s_real_space_refined.pdb" % (prefix)
                if self.phenix_version <= 1.18
                else "%s_real_space_refined_???.pdb" % (prefix)  # phenix version 1.19 and higher
            )
            try:
                pdb_fles = glob.glob(pdb_pattern)
                pdb_fles.sort()
                outpdb = pdb_fles[-1]
            except IndexError:
                outpdb = "%s_real_space_refined_000.pdb" % (prefix)
        else:  # not correctly finished
            outpdb = "not_a_file"

        return outpdb

    def real_space_refinement_ccp4(self, ccp4_in, pdb_in, resolution):
        """
        Real space refinement based on ccp4 file
        run phenix.real_space_refine
        Some parameters have changed between version 1.7, 1.8 and 1.9 hence the weird construction to grap the version
        """

        ccp4_name = get_name(ccp4_in)

        # Specify phenix version dependent parameters

        rotamer_restraints = "rotamers.restraints.enabled=False"  # rotamers.fit=all?

        prefix = "%s_phenix" % (ccp4_name)
        if self.phenix_version <= 1.18:
            output_prefix = 'output.file_name_prefix=%s'% (prefix)
            model_format  = 'output.model_format=pdb'
            # outpdb        = "%s_real_space_refined.pdb"%(mtz_name)
            ramachandran_restraints = 'ramachandran_restraints=False'
        else: # phenix version 1.19 and higher
            output_prefix = "output.prefix=%s" % (prefix)
            model_format = "model_format=pdb"
            ramachandran_restraints = "ramachandran_plot_restraints.enable=False"

        additional_keywords_line = ""
        if len(self.additional_real_keywords) > 0:
            for keyword in self.additional_real_keywords:
                additional_keywords_line += "%s " % (keyword)

        command = (
            "phenix.real_space_refine %s %s %s "
            "geometry_restraints.edits.excessive_bond_distance_limit=1000 refinement.run=minimization_global+adp scattering_table=%s c_beta_restraints=False %s refinement.macro_cycles=%d refinement.simulated_annealing=every_macro_cycle nproc=4 %s %s %s ignore_symmetry_conflicts=True resolution=%.2f %s"
            % (
                ccp4_in,
                self.additional,
                pdb_in,
                self.scattering_table,
                output_prefix,
                self.real_cycles,
                model_format,
                rotamer_restraints,
                ramachandran_restraints,
                resolution,
                additional_keywords_line,
            )
        )
        real = subprocess.call(shlex.split(command))

        # Find output file
        if real == 0:  # correctly finished. search for the last refined structure
            pdb_pattern = (
                "%s_real_space_refined.pdb" % (prefix)
                if self.phenix_version <= 1.18
                else "%s_real_space_refined_???.pdb" % (prefix) # phenix version 1.19 and higher
            )
            try:
                pdb_fles = glob.glob(pdb_pattern)
                pdb_fles.sort()
                outpdb = pdb_fles[-1]
            except IndexError:
                outpdb = "%s_real_space_refined_000.pdb" % (prefix)
        else:  # not correctly finished
            outpdb = "not_a_file"

        return outpdb
