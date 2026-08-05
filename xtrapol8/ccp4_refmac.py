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

import os
import re
from difflib import get_close_matches

import coot_headless_api
import iotbx.pdb
from iotbx.file_reader import any_file
from libtbx import adopt_init_args
from mmtbx.scaling.matthews import p_vm_calculator

from .Fextr_utils import get_name, redirect_stdout_and_stderr
from .programs.dm import dm
from .programs.fft import fft
from .programs.refmac import refmac
from .programs.uniqueify import uniqueify

CHAPI = coot_headless_api.molecules_container_t(False)


class Refmac_refinement:
    def __init__(
        self,
        mtz_in,
        pdb_in,
        additional,
        F_column_labels="QFEXTR",
        fill_missing=False,
        refinement_weight="AUTO",
        refinement_weight_sigmas="NOEX",
        refinement_weighting_term=0.2,
        refinement_type="REST",
        TLS=False,
        TLS_cycles=20,
        bfac_set=30,
        twinning=False,
        Brefinement="ISOT",
        cycles=5,
        external_restraints=None,
        jelly_body_refinement=False,
        jelly_body_sigma=0.03,
        jelly_body_additional_restraints=None,
        map_sharpening=False,
        additional_reciprocal_keywords=[],
    ):
        # Try this instead of huge repetition of the arguments
        adopt_init_args(self, locals())

        self.mtz_name = get_name(self.mtz_in)

    def do_refmac_reciprocal_space_refinement(self, mtz_out, pdb_out, log_file):
        """
        Run refmac.
        BFAC SET <value> default 30
        """
        if self.fill_missing:
            self.mtz_in = uniqueify(self.mtz_in)

        if self.refinement_weight == "AUTO":
            refinement_weight = "AUTO"
        else:
            refinement_weight = "%s %s %f" % (
                self.refinement_weight_sigmas,
                self.refinement_weight,
                self.refinement_weighting_term,
            )

        if self.jelly_body_refinement and self.cycles < 25:
            self.cycles *= 5

        return_code = refmac(
            hklin=self.mtz_in,
            xyzin=self.pdb_in,
            hklout=mtz_out,
            xyzout=pdb_out,
            libins=self.additional.split(),
            log_file=log_file,
            f_label=self.F_column_labels,
            cycles=self.cycles,
            weight=refinement_weight,
            refi_type=self.refinement_type,
            refi_bref=self.Brefinement,
            twinning=self.twinning,
            use_tls=self.TLS,
            tls_cycles=self.TLS_cycles,
            bfac_set=self.bfac_set,
            use_jelly_body=self.jelly_body_refinement,
            jelly_body_sigma=self.jelly_body_sigma,
            additional_jelly_body_restraints=self.jelly_body_additional_restraints,
            external_restraints=self.external_restraints,
            map_sharpening=self.map_sharpening,
            additional_keywords=self.additional_reciprocal_keywords,
        )

        ccp4_map_name = re.sub(r".mtz$", "_2mFo-DFc_filled.ccp4", mtz_out)
        fft(mtz_out, ccp4_map_name, "2FOFCWT", "PH2FOFCWT")

        ccp4_diff_map_name = re.sub(r".mtz$", "_mFo-DFc.ccp4", mtz_out)
        fft(mtz_out, ccp4_diff_map_name, "FOFCWT", "PHFOFCWT")

        return return_code

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

    def do_density_modification(self, mtz_in, mtz_out, combine, cycles, log_file):
        "Perform density modification with dm"
        solc = self.get_solvent_content()
        dm(mtz_in, mtz_out, solc, combine, cycles, self.F_column_labels, log_file)
        ccp4_map_name = re.sub(r".mtz$", ".ccp4", mtz_out)
        fft(mtz_out, ccp4_map_name, "FDM", "PHIDM")

    def reciprocal_space_refinement(self):
        try:
            if self.F_column_labels.lower().startswith("q"):
                maptype = "q" + self.F_column_labels.lower()[1:].capitalize()
            else:
                maptype = self.F_column_labels.lower().capitalize()
            mtz_out = re.sub(
                r"%s" % (maptype),
                "2m%s-DFc_refmac_reciprocal_space_001.mtz" % (maptype),
                self.mtz_name,
            )
            pdb_out = re.sub(
                r"%s" % (maptype),
                "2m%s-DFc_refmac_reciprocal_space_001.pdb" % (maptype),
                self.mtz_name,
            )
            log_file = re.sub(
                r"%s" % (maptype),
                "2m%s-DFc_refmac_reciprocal_space_001.log" % (maptype),
                self.mtz_name,
            )
            if mtz_out == self.mtz_name:
                raise AttributeError
        except AttributeError:
            mtz_out = "%s_refmac_reciprocal_space_001.mtz" % (self.mtz_name)
            pdb_out = "%s_refmac_reciprocal_space_001.pdb" % (self.mtz_name)
            log_file = "%s_refmac_reciprocal_space_001.log" % (self.mtz_name)

        print("Running Refmac. Please wait...")
        reciprocal = self.do_refmac_reciprocal_space_refinement(
            mtz_out, pdb_out, log_file
        )

        if reciprocal != 0:
            mtz_out = "not_a_file"
            pdb_out = "refinement_did_not_finish_correcty"

        return mtz_out, pdb_out

    def ccp4_dm(self, pdb_in, combine, cycles):
        """
        Use dm to perform density modification
        """
        mtz_for_dm = re.sub(r".pdb$", ".mtz", pdb_in)
        mtz_out_dm = re.sub(r".pdb$", "_dm.mtz", pdb_in)
        if os.path.isfile(mtz_for_dm):
            log_file = re.sub(r".mtz$", ".log", mtz_out_dm)
            print(f"Running density modification, see {log_file}. Please wait...")
            self.do_density_modification(
                mtz_for_dm, mtz_out_dm, combine, cycles, log_file
            )
        return mtz_out_dm


class Coot_refinement:
    def __init__(self, ligands, additional=""):
        self.ligands = ligands
        self.additional = additional

    def ligands_refinement(self, pdb_in, imol):
        pdb_hier = iotbx.pdb.hierarchy.input(file_name=pdb_in)
        hier = pdb_hier.hierarchy

        if len(self.ligands) > 0:
            CHAPI.set_map_weight(120)
            for lig in self.ligands:
                for chain in hier.chains():
                    if chain.is_protein():
                        for res_group in chain.residue_groups():
                            for atom_group in res_group.atom_groups():
                                if atom_group.resname == lig:
                                    CHAPI.refine_residues(
                                        imol,
                                        chain.id,
                                        res_group.resseq_as_int(),
                                        "",  # insertion code
                                        atom_group.altloc,
                                        "SINGLE",  # refinement mode
                                        1000,  # number of cycles
                                    )

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

        return column_labels_new.split(",")

    def do_coot_real_space_refinement_mtz(self, pdb_in, mtz_in, column_labels):
        """
        Run real-space refinement with COOT (CHAPI) using an mtz file
        """
        f_label = column_labels.split(",")[0].strip()
        p_label = column_labels.split(",")[1].strip()
        labels = f"{f_label},{p_label}"
        if not self.check_mtz_column(mtz_in, labels):
            print(f"{mtz_in}: column labels {labels} not found")
            f_label, p_label = self.suggest_mtz_column(
                mtz_in, labels, column_label_strings_constraint="2"
            )
            print(f"{mtz_in}: columns {f_label},{p_label} will be used")
        imap = CHAPI.read_mtz(mtz_in, f_label, p_label, "", False, False)
        pdb_out = self.do_coot_real_space_refinement(pdb_in, imap, mtz_in)
        return pdb_out

    def do_coot_real_space_refinement_ccp4(self, pdb_in, ccp4_in):
        """
        Real-space refinement with COOT (CHAPI) using a ccp4 file
        """
        imap = CHAPI.read_ccp4_map(ccp4_in, False)
        pdb_out = self.do_coot_real_space_refinement(pdb_in, imap, ccp4_in)
        return pdb_out

    def do_coot_real_space_refinement(self, pdb_in, imap, map_source_filename):
        "Real-space refinement with COOT (CHAPI) with an already-loaded map"
        CHAPI.set_use_gemmi(False)
        imol = CHAPI.read_coordinates(pdb_in)
        for cif in self.additional.split():
            if cif.endswith(".cif"):
                CHAPI.import_cif_dictionary(cif, imol)

        CHAPI.set_imol_refinement_map(imap)

        self.ligands_refinement(pdb_in, imol)

        CHAPI.set_map_weight(20)
        CHAPI.set_use_rama_plot_restraints(True)
        real = CHAPI.refine_residues_using_atom_cid(imol, "//", "ALL", 1000)

        CHAPI.set_add_waters_water_to_protein_distance_lim_min(2.4)
        CHAPI.set_add_waters_water_to_protein_distance_lim_max(3.6)
        # CHAPI.add_waters(imol, imap)  # TODO: Fix memory error

        name = get_name(map_source_filename)
        pdb_out = f"{name}_coot_real_space_refined_000.pdb"
        CHAPI.write_coordinates(imol, pdb_out)

        CHAPI.close_molecule(imol)
        CHAPI.close_molecule(imap)
        
        if real != 1:  #real space refinement incorrectly finished
            print("real space refinement with coot failed")
        else:
            print("real space refinement with coot successful")

        return pdb_out

    def get_mtz_resolution(self, mtz_in):
        """
        Easy extraction of the resolution boundaries of an mtz file
        """
        reflections = any_file(mtz_in, force_type="hkl", raise_sorry_if_errors=True)
        low_res, high_res = reflections.file_content.file_content().max_min_resolution()

        return low_res, high_res

    def real_space_refinement_mtz(self, mtz_in, pdb_in, column_labels):
        """
        Real space refinement within COOT based on an mtz file and spefic columns
        """
        print("Running Real space refinement in COOT")
        coot_log = f"{get_name(mtz_in)}_coot_real_space_refined_000.log"
        with redirect_stdout_and_stderr(coot_log):
            pdb_out = self.do_coot_real_space_refinement_mtz(
                pdb_in, mtz_in, column_labels
            )

        return pdb_out

    def real_space_refinement_ccp4(self, ccp4_in, pdb_in, resolution):
        """
        Real space refinement within COOT based on a ccp4 file. No ambiguitiy stemming from multiple different columns

        Resolution is not required, but passed as an argument to remain consistent with real_space_refinement_ccp4 function within the Phenix_real_space_refinement class
        """
        print("Running Real space refinement in COOT. Please wait...")
        coot_log = f"{get_name(ccp4_in)}_coot_real_space_refined_000.log"
        with redirect_stdout_and_stderr(coot_log):
            pdb_out = self.do_coot_real_space_refinement_ccp4(pdb_in, ccp4_in)

        return pdb_out
