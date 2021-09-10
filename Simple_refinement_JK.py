from __future__ import division, print_function
import Fextr
import numpy as np
from Fextr_utils import *
from Fextr import SymManager, DataHandler, FobsFobs, Filesandmaps
import os,sys
from cctbx import miller, xray
from iotbx.file_reader import any_file
from mmtbx.scaling.matthews import p_vm_calculator
from libtbx import adopt_init_args
from cctbx import miller
from Fextr_utils import get_name
import glob


def run_simple_refinement(params, P, outdir_and_mtz_file_off_on_outname):

    # init
    # get output directory for Xtrapol8
    outdir = outdir_and_mtz_file_off_on_outname[0]
    os.chdir(outdir)
    # get mtz files
    mtz_on = outdir_and_mtz_file_off_on_outname[1]
    total = outdir_and_mtz_file_off_on_outname[3]
    # get log file
    logname = outdir_and_mtz_file_off_on_outname[2]
    log = open(logname, "w")
    # get outname
    outname = outdir_and_mtz_file_off_on_outname[4]
    # check mtz files
    check_all_files([mtz_on])
    #fix an occupancy
    occ=1
    # get reflections
   # reflections_off, reflections_on = open_mtz(mtz_on)
    reflections_on = any_file(mtz_on, force_type="hkl", raise_sorry_if_errors=True)


    print(
        '-----------------------------------------\nDATA PREPARATION\n-----------------------------------------')
    DH = DataHandler(params.Xtrapol8.input.reference_pdb, None,
                     params.Xtrapol8.input.additional_files,
                     outdir, mtz_on)

    print('CHECKING FILES\n==============================================================')
    # Check format of inputs files
    check_all_files([DH.pdb_in])
    DH.open_pdb_or_cif()

#converting I to F and adding Rfree
    # run truncate to get Rfree
    script_out = 'launch_truncate.sh'
    i = open(script_out, 'w')
    i.write("#!/bin/bash\n\
truncate hklin %s hklout tmp.mtz << EOF\n\
TITLE Truncate run on native data\n\
TRUNCATE YES\n\
ANOMALOUS NO\n\
LABIN IMEAN=IMEAN SIGIMEAN=SIGIMEAN \n\
SYMMETRY %s\n\
EOF\n\
uniqueify tmp.mtz\n\
out=`basename %s .mtz` \n\
mypath=`dirname %s` \n\
mv tmp-unique.mtz ${mypath}/${out}_FPFREE.mtz\n\
rm tmp.mtz tmp-unique.log" % (mtz_on, 'P6122', mtz_on, mtz_on))
    i.close()
    os.system('bash %s' % (script_out))
    new_mtz_on = (mtz_on.split('.mtz')[0])+ '_FPFREE.mtz'
    # Extract space group and unit cell from model pdb file
  #  DH.get_UC_and_SG()

    # Check if all cif files are given and extract the three-letter codes
    print("----Check additional files----")
    DH.check_additional_files(SG=P.spacegroup_on, UC= [P.a_on, P.b_on, P.c_on, P.alpha_on, P.beta_on, P.gamma_on])
    print('---------------------------')

    # extract columns from mtz files that needs to be substracted
    print("----Column extraction from reflection files----")
#    _, DH.fobs_on = DH.extract_fobs(None, reflections_on, params.Xtrapol8.input.low_resolution, params.Xtrapol8.input.high_resolution, log)
    print('---------------------------')

#    DH.generate_Rfree(DH.fobs_on, 0.05)
#    DH.generate_f_model(DH.fobs_on)

    JK_files_table = np.array(['occupancy', 'map_type', 'Refinement', 'map', 'model', 'labels', 'total', 'residlist']) #init the table
    print('------------- Running Refinements ---------------')
    # FoFo = FobsFobs(log, DH.fobs_on_scaled, DH.fobs_off_scaled)
    # FoFo.calculate_fdiff(kweight_scale=params.Xtrapol8.f_and_maps.kweight_scale)
    # FoFo.write_maps(DH.fmodel, DH.rfree, outname, qweighting=P.qFoFo_weight, kweighting=P.kFoFo_weight)
    Fextr = Fextrapolate(log,
                         None,
                         None,
                         None,
                         None,
                         None,
                         None,
                         None,
                         None,
                         None,
                         None,
                         occ,
                         name_out=outname,
                         neg_refl_handle=params.Xtrapol8.f_and_maps.negative_and_missing,
                         simple_refinement=True)

    mtz_out, pdb_rec, pdb_real, pdb_rec_real, refinement_rec, refinement_real, refinement_rec_real = Fextr.phenix_phenix_refinements(mtz_F=new_mtz_on,
                           mtz_map=mtz_on,pdb_in=DH.pdb_in,
                                                                                                   additional=DH.additional,
                                                                                                   keywords=params.Xtrapol8.refinement.phenix_keywords)
    refinement_labels = '2FOFCWT,PH2FOFCWT'
    newdir=os.getcwd()
    #create table to analyse the results of JK
    twoFextrFclist_refined = np.array([occ, 'F', refinement_rec, str(newdir + '/' + mtz_out), str(newdir + '/' +pdb_rec), refinement_labels, total, None])
    JK_files_table = np.vstack((JK_files_table, twoFextrFclist_refined))
    twoFextrFclist_refined = np.array([occ, 'F', refinement_real, None, str(newdir + '/' +pdb_real), None, total, None])
    JK_files_table = np.vstack((JK_files_table, twoFextrFclist_refined))
    twoFextrFclist_refined = np.array([occ, 'F', refinement_rec_real, None, str(newdir + '/' +pdb_rec_real), None, total, None])
    JK_files_table = np.vstack((JK_files_table, twoFextrFclist_refined))

    return (JK_files_table)

class Fextrapolate(object):
    """
    Class for the calculation, analysis and usage of extrapolated structure factors.
    """

    def __init__(self,
                 log,
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
                 occ=1,
                 name_out='Fextrapolate',
                 neg_refl_handle='fill_missing',
                 simple_refinement=False):
        self.log = log
        self.fdif = fdif
        self.fdif_q = fdif_q
        self.fdif_k = fdif_k
        self.sigf_diff = sigf_diff
        self.q_values = q_values
        self.k_values = k_values
        self.fobs_on = fobs_on
        self.fobs_off = fobs_off
        self.fmodel_fobs_off = fmodel_fobs_off
        self.rfree = rfree
        self.occ = occ
        self.alf = 1 / (self.occ)
        self.neg_refl_handle = neg_refl_handle
        if not simple_refinement:
            self.indices = fobs_off.indices()
            self.get_UC_and_SG()
        # print("Calculating Fextr and maps for occupancy %.3f." %(occ))
        if neg_refl_handle in ['no_fill', 'reject_no_fill', 'zero_no_fill', 'fcalc_no_fill', 'fref_no_fill',
                               'truncate_no_fill']:  # 'addconstant_no_fill', 'massage_no_fill'
            self.fill_missing = False
        else:
            self.fill_missing = True

        self.name_out = "%s_occ%.3f" % (name_out, occ)

    def phenix_phenix_refinements(self,
                                  mtz_F=None,
                                  mtz_map=None,
                                  pdb_in=None,
                                  additional=None,
                                  F_column_labels='QFEXTR',
                                  column_labels='2FOFCWT,PH2FOFCWT',
                                  keywords={}):

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

        ref = Phenix_refinements(mtz_F,
                                                    pdb_in,
                                                    additional=additional,
                                                    F_column_labels=F_column_labels,
                                                    strategy=keywords.refine.strategy,
                                                    rec_cycles=keywords.main.cycles,
                                                    real_cycles=keywords.real_space_refine.cycles,
                                                    wxc_scale=keywords.target_weights.wxc_scale,
                                                    wxu_scale=keywords.target_weights.wxu_scale,
                                                    solvent=keywords.main.ordered_solvent,
                                                    sim_annealing=keywords.main.simulated_annealing,
                                                    sim_annealing_pars=keywords.simulated_annealing,
                                                    map_sharpening=keywords.map_sharpening.map_sharpening,
                                                    weight_sel_crit=keywords.target_weights.weight_selection_criteria,
                                                    log=self.log)

        # print("Refinements:", file=self.log)
        # print("Refinements:")

        print("RECIPROCAL SPACE REFINEMENT WITH %s AND %s" % (mtz_F, pdb_in))
        mtz_out_rec, pdb_out_rec = ref.phenix_reciprocal_space_refinement()
        refinement_rec = False
        print("Output reciprocal space refinement:", file = self.log)
        print("----------------")
        print("Output reciprocal space refinement:")
        if os.path.isfile(pdb_out_rec):
            print("    pdb-file: %s" % (pdb_out_rec), file = self.log)
            print("    pdb-file: %s" % (pdb_out_rec))
        else:
            print("    pdb-file not found, %s incorrectly returned" % (pdb_in), file = self.log)
            print("    pdb-file not found, %s incorrectly returned" % (pdb_in))
            pdb_out_rec = pdb_in
        if os.path.isfile(mtz_out_rec):
            print("    mtz-file: %s" % (mtz_out_rec), file = self.log)
            print("    mtz-file: %s" % (mtz_out_rec))
            refinement_rec = 'reciprocal'
        else:
            print("    mtz-file not found. Refinement failed.", file = self.log)
            print("    mtz-file not found. Refinement failed.")
        print("----------------")

        if keywords.density_modification.density_modification:
            print("DENSITY MODIFICATION WITH %s AND %s" % (mtz_F, pdb_out_rec))
            # mtz_dm = ref.phenix_density_modification(mtz_out_rec, pdb_out_rec)
            mtz_dm = ref.ccp4_dm(pdb_out_rec, keywords.density_modification.combine,
                                 keywords.density_modification.cycles)
            print("Output density modification:", file=self.log)
            print("Output density modification:")
            if os.path.isfile(mtz_dm):
                print("    mtz-file: %s" % (mtz_dm), file=self.log)
                print("    mtz-file: %s" % (mtz_dm))
            else:
                print("    mtz-file not found. Density modification failed.", file=self.log)
                print("    mtz-file not found. Density modification failed.")

        print("REAL SPACE REFINEMENT WITH %s AND %s" % (mtz_map, pdb_in))
        pdb_out_real = ref.phenix_real_space_refinement(mtz_map, pdb_in, column_labels)
        refinement_real = False
        print("Output real space refinement:", file=self.log)
        print("----------------")
        print("Output real space refinement:")
        if os.path.isfile(pdb_out_real):
            print("    pdb-file: %s" % (pdb_out_real), file=self.log)
            print("    pdb-file: %s" % (pdb_out_real))
            refinement_real = 'real'
        else:
            print("    pdb-file not found, %s incorrectly returned" % (pdb_in), file=self.log)
            print("    pdb-file not found, %s incorrectly returned" % (pdb_in))
            pdb_out_real = pdb_in
        print("----------------")

        if (keywords.density_modification.density_modification and os.path.isfile(mtz_dm)):
            print("REAL SPACE REFINEMENT WITH %s AND %s" % (mtz_dm, pdb_out_rec))
            pdb_out_rec_real = ref.phenix_real_space_refinement(mtz_dm, pdb_out_rec, 'FWT,PHWT')
            refinement_rec_real = False
            print("Output real space refinement after reciprocal space refinement:", file=self.log)
            print("----------------")
            print("Output real space refinement after reciprocal space refinement:")
            if os.path.isfile(pdb_out_rec_real):
                print("    pdb-file: %s" % (pdb_out_rec_real), file=self.log)
                print("    pdb-file: %s" % (pdb_out_rec_real))
                refinement_rec_real = 'reciprocal + real'
            else:
                print("    pdb-file not found, %s incorrectly returned" % (pdb_out_rec), file=self.log)
                print("    pdb-file not found, %s incorrectly returned" % (pdb_out_rec))
                pdb_out_rec_real = pdb_out_rec
            print("----------------")
        else:
            print("REAL SPACE REFINEMENT WITH %s AND %s" % (mtz_out_rec, pdb_out_rec))
            pdb_out_rec_real = ref.phenix_real_space_refinement(mtz_out_rec, pdb_out_rec, '2FOFCWT,PH2FOFCWT')
            refinement_rec_real = False
            print("Output real space refinement after reciprocal space refinement:", file=self.log)
            print("----------------")
            print("Output real space refinement after reciprocal space refinement:")
            if os.path.isfile(pdb_out_rec_real):
                print("    pdb-file: %s" % (pdb_out_rec_real), file=self.log)
                print("    pdb-file: %s" % (pdb_out_rec_real))
                refinement_rec_real = 'reciprocal + real'
            else:
                print("    pdb-file not found, %s incorrectly returned" % (pdb_out_rec), file=self.log)
                print("    pdb-file not found, %s incorrectly returned" % (pdb_out_rec))
                pdb_out_rec_real = pdb_out_rec
            print("----------------")

        return mtz_out_rec, pdb_out_rec, pdb_out_real, pdb_out_rec_real, refinement_rec, refinement_real, refinement_rec_real


class Phenix_refinements(object):
    def __init__(self,
                 mtz_in,
                 pdb_in,
                 additional='',
                 F_column_labels='QFEXTR',
                 strategy=['individual_sites', 'individual_adp'],
                 rec_cycles=5,
                 real_cycles=5,
                 wxc_scale=0.5,
                 wxu_scale=1.0,
                 solvent=False,
                 sim_annealing=False,
                 sim_annealing_pars={},
                 map_sharpening=False,
                 weight_sel_crit={},
                 log=sys.stdout):

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
                maptype = "q" + self.F_column_labels.lower()[1:].capitalize()
            else:
                maptype = self.F_column_labels.lower().capitalize()
            outprefix = re.sub(r"%s" % (maptype), "2m%s-DFc_reciprocal_space" % (maptype), self.mtz_name)
            if outprefix == self.mtz_name:
                raise AttributeError
        except AttributeError:
            outprefix = "%s_reciprocal_space" % (self.mtz_name)

        sim_annealing = ''
        if self.sim_annealing:
            sim_annealing += "simulated_annealing=True"
            for param in self.sim_annealing_pars.__dict__:
                if not param.startswith("_"):
                    sim_annealing += " simulated_annealing.%s=%s " % (
                    param, str(self.sim_annealing_pars.__dict__[param]))
        else:
            sim_annealing += "simulated_annealing=False"

        weight_selection_criteria = ""
        if self.weight_sel_crit.bonds_rmsd != None:
            weight_selection_criteria += "target_weights.weight_selection_criteria.bonds_rmsd=%.4f " % (
                self.weight_sel_crit.bonds_rmsd)
        if self.weight_sel_crit.angles_rmsd != None:
            weight_selection_criteria += "target_weights.weight_selection_criteria.angles_rmsd=%.4f " % (
                self.weight_sel_crit.angles_rmsd)
        if self.weight_sel_crit.r_free_minus_r_work != None:
            weight_selection_criteria += "target_weights.weight_selection_criteria.r_free_minus_r_work=%.4f " % (
                self.weight_sel_crit.r_free_minus_r_work)

        reciprocal = os.system("phenix.refine --overwrite %s %s  %s output.prefix=%s strategy=%s "
                               "main.number_of_macro_cycles=%d refinement.output.write_model_cif_file=False "
                               "refinement.input.xray_data.r_free_flags.disable_suitability_test=True "
                               "refinement.input.xray_data.r_free_flags.ignore_pdb_hexdigest=True "
                               "refinement.input.xray_data.r_free_flags.label='FreeR_flag' "
                                "refinement.input.xray_data.labels='F,SIGF' "
                               "refinement.input.xray_data.r_free_flags.test_flag_value=1 nproc=4 wxc_scale=%f wxu_scale=%f ordered_solvent=%s write_maps=true %s %s %s" % (
                               self.mtz_in, self.additional, self.pdb_in, outprefix, self.strategy, self.rec_cycles,
                               self.wxc_scale, self.wxu_scale, self.solvent, self.params, weight_selection_criteria,
                               sim_annealing))  # wxc_scale=0.021 #target_weights.optimize_xyz_weight=True

        # Find output files, automatically
        if reciprocal == 0:  # os.system has correctly finished, then search for the last refined structure
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
        else:  # os.system has not correctly finished
            mtz_out = "no_a_file"
            pdb_out = "refinement_did_not_finish_correcty"

        # Find output files: hardcoded appears easiest
        # mtz_out = [fle for fle in os.listdir(os.getcwd()) if outprefix in fle and fle.endswith('mtz')][0]
        # pdb_out = [fle for fle in os.listdir(os.getcwd()) if outprefix in fle and fle.endswith('pdb')][0]
        # mtz_out = "%s_001.mtz"%(outprefix)
        # pdb_out = "%s_001.pdb"%(outprefix)

        return mtz_out, pdb_out

    def phenix_real_space_refinement(self, mtz_in, pdb_in, column_labels):
        """
        use Bash line to run phenix.real_space_refine as usual (use of os.system is bad practice).
        Some parameters have changed between version 1.7, 1.8 and 1.9 hence the weird construction to grap the version
        """

        mtz_name = get_name(mtz_in)

        # Spacify phenix version dependent parameters
        try:
            phenix_version = int(re.search(r"phenix-1\.(.+?)\.", miller.__file__).group(
                1))  # This is not so robust. relies on the format being 'phenix.1.18.something' or 'phenix.1.18-something'
        except ValueError:
            phenix_version = int(re.search(r"phenix-1\.(.+?)\-", miller.__file__).group(1))
        except AttributeError:
            print('Update phenix! Verify that you are using at least Phenix.1.17.')
            phenix_version = 19  # let's assume then that the latest phenix is installed in case this fails for other reasons than a very old phenix version

        print("Phenix version 1.%d" % (phenix_version))

        if phenix_version >= 18:
            rotamer_restraints = 'rotamers.restraints.enabled=False'  # rotamers.fit=all?
        else:
            rotamer_restraints = 'rotamer_restraints=False'

        if phenix_version >= 19:
            output_prefix = 'output.prefix=%s' % (mtz_name)
            model_format = 'model_format=pdb'
            # outpdb        = "%s_real_space_refined_000.pdb"%(mtz_name)
        else:
            output_prefix = 'output.file_name_prefix=%s' % (mtz_name)
            model_format = 'output.model_format=pdb'
            # outpdb        = "%s_real_space_refined.pdb"%(mtz_name)

        real = os.system("phenix.real_space_refine %s %s %s "
                         "geometry_restraints.edits.excessive_bond_distance_limit=1000 refinement.run=minimization_global+adp scattering_table=n_gaussian ramachandran_restraints=False c_beta_restraints=False %s refinement.macro_cycles=%d refinement.simulated_annealing=every_macro_cycle nproc=4 %s label='%s' %s ignore_symmetry_conflicts=True" % (
                         mtz_in, self.additional, pdb_in, output_prefix, self.real_cycles, model_format, column_labels,
                         rotamer_restraints))

        # Find output file
        if real == 0:  # os.system has correctly finished. Then search for the last refined structure
            if phenix_version >= 19:
                try:
                    pdb_fles = glob.glob("%s_real_space_refined_???.pdb" % (mtz_name))
                    # [fle for fle in os.listdir(os.getcwd()) if "%s_independent_real_space_refined_0"%(mtz_name) in
                    #            fle and fle.endswith('pdb')]
                    pdb_fles.sort()
                    outpdb = pdb_fles[-1]
                except IndexError:
                    outpdb = "%s_real_space_refined_000.pdb" % (mtz_name)
            else:
                try:
                    pdb_fles = glob.glob("%s_real_space_refined.pdb" % (mtz_name))
                    # [fle for fle in os.listdir(os.getcwd()) if
                    #            "%s_independent_real_space_refined" % (mtz_name) in fle and fle.endswith('pdb')]
                    pdb_fles.sort()
                    outpdb = pdb_fles[-1]
                except IndexError:
                    outpdb = "%s_real_space_refined.pdb" % (mtz_name)
        else:  # os.system has not correctly finished
            outpdb = "not_a_file"

        return outpdb
