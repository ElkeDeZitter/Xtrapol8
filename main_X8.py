
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

import JK_utils
import Parameters
from Parameters import Parameters


def run_X8(outdir_and_mtz_file_off_on, params, DH, P, master_phil, startdir, only_Xtrapol8=False):

#init
    #get output directory for Xtrapol8
    outdir = outdir_and_mtz_file_off_on[0]
    os.chdir(outdir)
    #get mtz files
    mtz_off = outdir_and_mtz_file_off_on[1]
    mtz_on = outdir_and_mtz_file_off_on[2]
    #check mtz files
    check_all_files([mtz_off, mtz_on])
    #get reflections
    reflections_off, reflections_on = open_mtz(mtz_off, mtz_on)
    #get log file
    logname = outdir_and_mtz_file_off_on[3]
    log = open(logname, "w")

#function from Fextr.DH
    def extract_fobs(low_res, high_res):
        """
        Extract the actual reflections from the data files and cut at resolution limits (if set)
        For now Friedel pairs will have to be merged.
        """
        DH.fobs_off, DH.fobs_on = Column_extraction(reflections_off,
                                                        reflections_on,
                                                        low_res,
                                                        high_res,
                                                        log=log).extract_columns()

        if DH.fobs_off.anomalous_flag():
            print(
                "I promised to keep the anomalous flags, but that was a lie. Xtrapol8 is not yet ready to handle anomalous data. For now, your Friedel pairs will be merged.",
                file=log)
            print(
                "I promised to keep the anomalous flags, but that was a lie. Xtrapol8 is not yet ready to handle anomalous data. For now, your Friedel pairs will be merged.")
            DH.fobs_off = DH.fobs_off.average_bijvoet_mates()
            DH.fobs_on = DH.fobs_on.average_bijvoet_mates()

        DH.fobs_off = DH.fobs_off.map_to_asu()
        # self.fobs_off = self.resolution_cutoff(self.fobs_off, low_res, high_res)
        DH.fobs_on = DH.fobs_on.map_to_asu()
        # self.fobs_on  = self.resolution_cutoff(self.fobs_on, low_res, high_res)

        # self.fobs_off = self.extract_colums(self.reflections_off, low_res, high_res)
        # self.fobs_on  = self.extract_colums(self.reflections_on, low_res, high_res)
        # self.fobs_on = []
        # for on in self.reflections_on:
        # f_on = self.extract_colums(on, res)
        # self.fobs_on.append(f_on)
        return (DH.fobs_off, DH.fobs_on)

#Almost original Xtrapol8 script: changed log in FobsFobs and Fextrapolate
    ##################################################################
    # extract columns from mtz files that needs to be substracted
    print("----Column extraction from reflection files----")
    DH.fobs_off, DH.fobs_on = extract_fobs(params.Xtrapol8.input.low_resolution, params.Xtrapol8.input.high_resolution)
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
            # err_on == 1 while err_off == 0 can only happen if pointless failed the reindexing
            print(
                "Changing space group and/or unit cell parameters of the triggered data in order to match with the model (if not already done by Pointless).",
                file=log)
            print(
                "Changing space group and/or unit cell parameters of the triggered data in order to match with the model (if not already done by Pointless).")
            DH.fobs_on = make_miller_array(DH.fobs_on.data(),
                                           DH.fobs_on.sigmas(),
                                           DH.SG,
                                           DH.UC,
                                           DH.fobs_on.indices())
        if err_off == 1:  # very weird if this would happen
            print(
                "Changing space group and/or unit cell parameters of the reference data in order to match with the model. This is weird! Was the model generated after refinement in the reference mtz??",
                file=log)
            print(
                "Changing space group and/or unit cell parameters of the reference data in order to match with the model. This is weird! Was the model generated after refinement in the reference mtz??")
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



    # Print reflection statistics to screen and log-file
    print("----Summary reference data:----")
    print(DH.fobs_on.show_comprehensive_summary())
    print("----Summary reference data:----", file=log)
    print(DH.fobs_on.show_comprehensive_summary(f=log), file=log)
    print("----Summary triggered data:----")
    print(DH.fobs_off.show_comprehensive_summary())
    print("------------------------------------")
    print("----Summary triggered data:----", file=log)
    print(DH.fobs_off.show_comprehensive_summary(f=log), file=log)
    print("------------------------------------", file=log)

    # Scale the data sets. Different steps taken:
    # 1) generate a set of Rfree reflections
    # 2) generate fmodel with reference reflections, model_in and Rfree flag
    # 3) update fmodel with bulk-solvent scaling, i.e. scale fobs_off with fcalc
    # 4) print R-factors and correlation coefficients between fobs_off and fcalc (this are the traditional Rwork/Rfree values)
    # 5) check the reflections agreement between the reference and triggered reflections (and all other generated arrays)
    # 6) scale triggered reflections with reference reflections
    # 7) check the reflections agreement between the reference and triggered reflections (and all other generated arrays)
    # 8) update fmodel with scaled reference and remove reflections that were rejected during the triggered scaling
    # 9) print R-factors and correlation coefficients between fobs_off and fcalc which an indicator for isomorfism (this is the famous Riso)
    print("----Scaling Fmodel with Freference----", file=log)
    print("----Scaling Fmodel with Freference----")
    DH.generate_Rfree(DH.fobs_off, 0.05)
    DH.generate_f_model(DH.fobs_off)
    print("fmodel generation, can take a few minutes")
    try:
        DH.fmodel.update_all_scales(show=True)  # , log=log)
    except RuntimeError:
        print("The model-based structure factors could not be scaled using the fast method. Try again with slow method")
        DH.fmodel.update_all_scales(show=True, fast=False)
    DH.fmodel.show()
    DH.fmodel.info().show_rfactors_targets_scales_overall(out=sys.stdout)
    print("Fobs,reference and Fcalc,reference scaled using mmtbx f_model")
    print('R_work:', DH.fmodel.r_work())
    print('R_free:', DH.fmodel.r_free())
    print("Fobs,reference and Fcalc,reference scaled using mmtbx f_model", file=log)
    print('R_work:', DH.fmodel.r_work(), file=log)
    print('R_free:', DH.fmodel.r_free(), file=log)
    DH.fobs_on = DH.get_common_indices_and_Fobs_off(
        DH.fobs_on, log)  # compare reflections and reassemble off_state data set after scaling of f_obs_off
    print("----Scaling Ftriggered with Freference----", file=log)
    print("----Scaling Ftriggered with Freference----")
    DH.scale_fobss(params.Xtrapol8.scaling.b_scaling, log)
    DH.fobs_on_scaled = DH.get_common_indices_and_Fobs_off(
        DH.fobs_on_scaled, log)  # compare reflections and reassemble off_state data set after scaling of f_obs_on
    DH.update_fmodel(DH.fobs_off_scaled)  # alter fmodel to remove reflections that were removed during fobs_on scaling
    riso, cciso = compute_r_factors(DH.fobs_off_scaled, DH.fobs_on_scaled, DH.rfree, log=log)
    print("Overall R_iso = %.4f " %(riso))
    print("Overall cc_iso = %.4f " %(cciso))
    print("Overall R_iso = %.4f " %(riso), file=log)
    print("Overall cc_iso = %.4f " %(cciso), file=log)

    P.get_parameters_DH_scaleit_X8(DH)

    # Print scaled reflections statistics to screen and log-file. Stupif that I'm printing everything to the screen and to the log-file. There must be something more elegant


    print("---------------------------", file=log)
    print("Resolution range= %2f - %2f A" % (P.dmax, P.dmin), file=log)
    print("---------------------------", file=log)
    print('---------------------------')
    print("Resolution range= %2f - %2f A" % (P.dmax, P.dmin))
    print('---------------------------')

    print("----Summary triggered scaled and sorted common reflections with reference data:----")
    print(DH.fobs_on_scaled.show_comprehensive_summary())
    print("----Summary triggered scaled and sorted common reflections with reference data:----", file=log)
    print(DH.fobs_on_scaled.show_comprehensive_summary(f=log), file=log)
    print("----Summary reference scaled and sorted common reflections with triggered data:----")
    print(DH.fobs_off_scaled.show_comprehensive_summary())
    print("----Summary reference scaled and sorted common reflections with triggered data:----", file=log)
    print(DH.fobs_off_scaled.show_comprehensive_summary(f=log), file=log)
    # time.sleep(3)

    print('-----------------------------------------')
    print('DATA PREPARATION DONE')
    print('-----------------------------------------')

    ################################################################
    print("CALCULATE (WEIGHTED) FO-FO FOURIER DIFFERENCE MAP")
    print('-----------------------------------------')

    print('-----------------------------------------', file=log)
    print("CALCULATE (WEIGHTED) FO-FO FOURIER DIFFERENCE MAP", file=log)
    print('-----------------------------------------', file=log)

    # calculate Fo-Fo and write map coefficients to mtz, ccp4 and xplor file
    FoFo = FobsFobs(log, DH.fobs_on_scaled, DH.fobs_off_scaled)
    FoFo.calculate_fdiff(kweight_scale=params.Xtrapol8.f_and_maps.kweight_scale)
    FoFo.write_maps(DH.fmodel, DH.rfree, P.outname, qweighting=P.qFoFo_weight, kweighting=P.kFoFo_weight)

    print("%s maps generated in mtz,ccp4 and xplor format"%(params.f_and_maps.fofo_type), file=log)
    #print("------------------------------------", file=log)
    print("%s maps generated in mtz,ccp4 and xplor format"%(params.f_and_maps.fofo_type))
    #print("------------------------------------")

    ################################################################

    # Use xplor map to find and integrate the peaks, annotate the peaks to residues and keep a list with the most important residues only based on a user-defined Z-score level
    print("\n************Map explorer************", file=log)
    print("\n************Map explorer************")
    if params.Xtrapol8.map_explorer.radius == None:
        params.Xtrapol8.map_explorer.radius = P.dmin
    map_expl_out_FoFo = map_explorer(FoFo.xplor_name, DH.pdb_in, params.Xtrapol8.map_explorer.radius,
                                     params.Xtrapol8.map_explorer.threshold, params.Xtrapol8.map_explorer.peak)
    residlist_zscore = Map_explorer_analysis(peakintegration_file=map_expl_out_FoFo, log=log).residlist_top(
        Z=params.Xtrapol8.map_explorer.z_score)

    print("FoFo map explored. Results in %s, residue list in residlist.txt and residues associated to highestpeaks in %s\n"
          %(map_expl_out_FoFo, residlist_zscore), file=log)
    print("FoFo map explored. Results in %s, residue list in residlist.txt and residues associated to highestpeaks in %s\n"
          %(map_expl_out_FoFo, residlist_zscore))
    #print('---------------------------')
    map_expl_out_FoFo = os.path.abspath(check_file_existance(map_expl_out_FoFo))
    residlist_zscore = os.path.abspath(check_file_existance(residlist_zscore))

    # Generate plot which indicates the secondary structure and integrated peak volume
    print("----Generate difference map plot----")
    Map_explorer_analysis(peakintegration_file=map_expl_out_FoFo, ligands=DH.extract_ligand_codes(), log=log).get_ss(
        DH.pdb_in)
    # print('---------------------------')

    # Set parameter for whether or not q_weighting is applied on the FoFo calculation. This will be used later in extrapolated map calculation.
    # k-weighting is not yet set, this might pose problems
    if P.qFoFo_weight == True:
        P.FoFo_type = 'qFo-Fo'
    elif P.kFoFo_weight == True:
        P.FoFo_type = 'kFo-Fo'
    else:
        P.FoFo_type = 'Fo-Fo'

    print(
        '-----------------------------------------\nCALCULATE (WEIGHTED) FO-FO FOURIER DIFFERENCE MAP DONE\n-----------------------------------------')

    if params.Xtrapol8.output.generate_fofo_only:
        print(
            '-----------------------------------------\nDONE! FOURIER DIFFERENCE MAP CALCULATED AND ANALYSIS PERFORMED.\n-----------------------------------------')
        print(
            '-----------------------------------------\nDONE! FOURIER DIFFERENCE MAP CALCULATED AND ANALYSIS PERFORMED.\n-----------------------------------------', file=log)
        # change names to real output name in case the dummy name was used
        if P.outname == 'triggered':
            # print("Replacing the dummpy outname ('%s') with true outname('%s')" %(outname, params.Xtrapol8.output.outname), file=log)
            print("Replacing the dummpy outname ('%s') with true outname('%s')" % (P.outname, params.output.outname))
            tr = [os.path.join(root, fle) for root, dirs, files in os.walk(outdir) for fle in files if P.outname in fle]
            _ = [os.rename(fle, fle.replace(P.outname, params.output.outname)) for fle in tr]
            FoFo.mtz_name = FoFo.mtz_name.replace(P.outname, params.output.outname)

        script_coot = open_all_in_coot("%s/%s" % (outdir, FoFo.mtz_name), [DH.pdb_in], [], DH.additional, outdir,
                                       P.FoFo_type)

        modified_phil = master_phil.format(python_object=params)
        modified_phil.show(out=open("Xtrapol8_out.phil", "w"))

        log.close()
        if (params.Xtrapol8.output.open_coot and check_program_path('coot')):
            os.system("coot --script %s" % (script_coot))

        sys.exit()

    ################################################################
    print(
        '-----------------------------------------\nCALCULATE (Q/K-WEIGHTED) FEXTRAPOLATED MAPS\n-----------------------------------------')
    print(
    '-----------------------------------------\nCALCULATE (Q/K-WEIGHTED) FEXTRAPOLATED MAPS\n-----------------------------------------', file=log)

    # For each map type, keep track of the map_explorer files. Therefore generate lists with output file names. Not elegant, but works.
    if P.qFextr_map:
        qFextr_map_expl_fles = [map_expl_out_FoFo]
        qFextr_recref_mtz_lst = []
        qFextr_recref_pdb_lst = [DH.pdb_in]
        qFextr_realref_lst = [DH.pdb_in]
        qFextr_recrealref_lst = [DH.pdb_in]
    if P.kFextr_map:
        kFextr_map_expl_fles = [map_expl_out_FoFo]
        kFextr_recref_mtz_lst = []
        kFextr_recref_pdb_lst = [DH.pdb_in]
        kFextr_realref_lst = [DH.pdb_in]
        kFextr_recrealref_lst = [DH.pdb_in]
    if P.Fextr_map:
        Fextr_map_expl_fles = [map_expl_out_FoFo]
        Fextr_recref_mtz_lst = []
        Fextr_recref_pdb_lst = [DH.pdb_in]
        Fextr_recrealref_lst = [DH.pdb_in]
        Fextr_realref_lst = [DH.pdb_in]
    if P.qFgenick_map:
        qFgenick_map_expl_fles = [map_expl_out_FoFo]
        qFgenick_recref_mtz_lst = []
        qFgenick_recref_pdb_lst = [DH.pdb_in]
        qFgenick_realref_lst = [DH.pdb_in]
        qFgenick_recrealref_lst = [DH.pdb_in]
    if P.kFgenick_map:
        kFgenick_map_expl_fles = [map_expl_out_FoFo]
        kFgenick_recref_mtz_lst = []
        kFgenick_recref_pdb_lst = [DH.pdb_in]
        kFgenick_realref_lst = [DH.pdb_in]
        kFgenick_recrealref_lst = [DH.pdb_in]
    if P.Fgenick_map:
        Fgenick_map_expl_fles = [map_expl_out_FoFo]
        Fgenick_recref_mtz_lst = []
        Fgenick_recref_pdb_lst = [DH.pdb_in]
        Fgenick_realref_lst = [DH.pdb_in]
        Fgenick_recrealref_lst = [DH.pdb_in]
    if P.qFextr_calc_map:
        qFextr_calc_map_expl_fles = [map_expl_out_FoFo]
        qFextr_calc_recref_mtz_lst = []
        qFextr_calc_recref_pdb_lst = [DH.pdb_in]
        qFextr_calc_realref_lst = [DH.pdb_in]
        qFextr_calc_recrealref_lst = [DH.pdb_in]
    if P.kFextr_calc_map:
        kFextr_calc_map_expl_fles = [map_expl_out_FoFo]
        kFextr_calc_recref_mtz_lst = []
        kFextr_calc_recref_pdb_lst = [DH.pdb_in]
        kFextr_calc_realref_lst = [DH.pdb_in]
        kFextr_calc_recrealref_lst = [DH.pdb_in]
    if P.Fextr_calc_map:
        Fextr_calc_map_expl_fles = [map_expl_out_FoFo]
        Fextr_calc_recref_mtz_lst = []
        Fextr_calc_recref_pdb_lst = [DH.pdb_in]
        Fextr_calc_realref_lst = [DH.pdb_in]
        Fextr_calc_recrealref_lst = [DH.pdb_in]
    # fast_and_furious mode: no refinement, but need to keep track of structure factor files
    if (params.Xtrapol8.f_and_maps.fast_and_furious): Fextr_mtz_lst = []

    # Remove pickle file with Fextr stats because otherwise plot will consist results from previous runs
    if os.path.isfile('%s/Fextr_binstats.pickle' % (outdir)):
        os.remove('%s/Fextr_binstats.pickle' % (outdir))
    # Remove pickle file with Number of negatove reflections because otherwise plot will consist results from previous runs
    if os.path.isfile('%s/Fextr_negative.pickle' % (outdir)):
        os.remove('%s/Fextr_negative.pickle' % (outdir))
    # Remove python movie movie because otherwise old files will be shown in the pymol session
    if os.path.isfile('%s/pymol_movie.py' % (outdir)):
        os.remove('%s/pymol_movie.py' % (outdir))

    ################################################################
    # Loop over occupancies and calculate extrapolated structure factors
    for occ in params.Xtrapol8.occupancies.list_occ:
        Fextr = Fextrapolate(log,
                             FoFo.fdif,
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
                             name_out=P.outname,
                             neg_refl_handle=params.Xtrapol8.f_and_maps.negative_and_missing)

        # Results stored in folder depending on occupancy and whether q-weighting is applied.
        new_dirpath_q, new_dirpath_k, new_dirpath = Fextr.create_output_dirs(outdir)

        # for each F_and_map type calculate structure factors and maps
        for mp in P.final_maptypes:
            if mp in ('qFextr_map', 'qFgenick_map', 'qFextr_calc_map'):
                os.chdir(new_dirpath_q)
            elif mp in ('kFextr_map', 'kFgenick_map', 'kFextr_calc_map'):
                os.chdir(new_dirpath_k)
            else:
                os.chdir(new_dirpath)

            # Depending on maptype, do following steps:
            # 1) calculate the structure factors, write out to mtz file, handle negative reflections and generate associated plots or write to pickle file for later usuage
            #   calculate map coefficients and write to mtz, xplor and ccp4 files (latter only for mFo-DFc type)
            # 2) get stats and write to pickle file in order to generate plot afterwards with all map types and occupancies
            # 3) compute signal to noise and plot
            # As same steps are repeated, but with a sifferent Fextr-function, I should think of a more clever way to reduce the redundancy here (and in following if clauses)
            if mp == 'qFextr_map':
                Fextr.fextr(qweight=True, kweight=False, outdir_for_negstats=outdir)
                get_Fextr_stats(occ, Fextr.fextr_ms, Fextr.maptype, FoFo.fdif_q_ms, P.FoFo_type, outdir)
                compute_f_sigf(Fextr.fextr_ms, '%s_occupancy%.3f.png' % (Fextr.maptype, Fextr.occ), log=log)
                # cc_list.append(plot_F1_F2(DH.fobs_off_scaled,Fextr.fextr_ms, F1_name = "Freference",F2_name = "Fextr"))
            elif mp == 'qFgenick_map':
                Fextr.fgenick(qweight=True, kweight=False, outdir_for_negstats=outdir)
                get_Fextr_stats(occ, Fextr.fgenick_ms, Fextr.maptype, FoFo.fdif_q_ms, P.FoFo_type, outdir)
                compute_f_sigf(Fextr.fgenick_ms, '%s_occupancy%.3f.png' % (Fextr.maptype, Fextr.occ), log=log)
            elif mp == 'qFextr_calc_map':
                Fextr.fextr_calc(qweight=True, kweight=False, outdir_for_negstats=outdir)
                get_Fextr_stats(occ, Fextr.fextr_calc_ms, Fextr.maptype, FoFo.fdif_q_ms, P.FoFo_type, outdir)
                compute_f_sigf(Fextr.fextr_calc_ms, '%s_occupancy%.3f.png' % (Fextr.maptype, Fextr.occ), log=log)
            if mp == 'kFextr_map':
                Fextr.fextr(qweight=False, kweight=True, outdir_for_negstats=outdir)
                get_Fextr_stats(occ, Fextr.fextr_ms, Fextr.maptype, FoFo.fdif_k_ms, P.FoFo_type, outdir)
                compute_f_sigf(Fextr.fextr_ms, '%s_occupancy%.3f.png' % (Fextr.maptype, Fextr.occ), log=log)
            elif mp == 'kFgenick_map':
                Fextr.fgenick(qweight=False, kweight=True, outdir_for_negstats=outdir)
                get_Fextr_stats(occ, Fextr.fgenick_ms, Fextr.maptype, FoFo.fdif_k_ms, P.FoFo_type, outdir)
                compute_f_sigf(Fextr.fgenick_ms, '%s_occupancy%.3f.png' % (Fextr.maptype, Fextr.occ), log=log)
            elif mp == 'kFextr_calc_map':
                Fextr.fextr_calc(qweight=False, kweight=True, outdir_for_negstats=outdir)
                get_Fextr_stats(occ, Fextr.fextr_calc_ms, Fextr.maptype, FoFo.fdif_k_ms, P.FoFo_type, outdir)
                compute_f_sigf(Fextr.fextr_calc_ms, '%s_occupancy%.3f.png' % (Fextr.maptype, Fextr.occ), log=log)
            elif mp == 'Fextr_map':
                Fextr.fextr(qweight=False, kweight=False, outdir_for_negstats=outdir)
                get_Fextr_stats(occ, Fextr.fextr_ms, Fextr.maptype, FoFo.fdif_c_ms, P.FoFo_type, outdir)
                compute_f_sigf(Fextr.fextr_ms, '%s_occupancy%.3f.png' % (Fextr.maptype, Fextr.occ), log=log)
            elif mp == 'Fgenick_map':
                Fextr.fgenick(qweight=False, kweight=False, outdir_for_negstats=outdir)
                get_Fextr_stats(occ, Fextr.fgenick_ms, Fextr.maptype, FoFo.fdif_c_ms, P.FoFo_type, outdir)
                compute_f_sigf(Fextr.fgenick_ms, '%s_occupancy%.3f.png' % (Fextr.maptype, Fextr.occ), log=log)
            elif mp == 'Fextr_calc_map':
                Fextr.fextr_calc(qweight=False, kweight=False, outdir_for_negstats=outdir)
                get_Fextr_stats(occ, Fextr.fextr_calc_ms, Fextr.maptype, FoFo.fdif_c_ms, P.FoFo_type, outdir)
                compute_f_sigf(Fextr.fextr_calc_ms, '%s_occupancy%.3f.png' % (Fextr.maptype, Fextr.occ), log=log)
            else:
                print("%s not recognised as extrapolated map type" % (mp))

            # Use xplor map of type mFo-DFc to find and integrate the peaks, annotate the peaks to residues

            print("\n************Map explorer************")
            print("\n************Map explorer************", file=log)
            map_expl_out = map_explorer(Fextr.xplor_name_FoFc, DH.pdb_in, params.Xtrapol8.map_explorer.radius,
                                        params.Xtrapol8.map_explorer.threshold, params.Xtrapol8.map_explorer.peak, maptype=Fextr.maptype)

            # depending on the map-type, append the output-file of mapexplorer to the correct list
            if mp == 'qFextr_map':
                append_if_file_exist(qFextr_map_expl_fles, os.path.abspath(map_expl_out))
            elif mp == 'qFgenick_map':
                append_if_file_exist(qFgenick_map_expl_fles, os.path.abspath(map_expl_out))
            elif mp == 'qFextr_calc_map':
                append_if_file_exist(qFextr_calc_map_expl_fles, os.path.abspath(map_expl_out))
            elif mp == 'kFextr_map':
                append_if_file_exist(kFextr_map_expl_fles, os.path.abspath(map_expl_out))
            elif mp == 'kFgenick_map':
                append_if_file_exist(kFgenick_map_expl_fles, os.path.abspath(map_expl_out))
            elif mp == 'kFextr_calc_map':
                append_if_file_exist(kFextr_calc_map_expl_fles, os.path.abspath(map_expl_out))
            elif mp == 'Fextr_map':
                append_if_file_exist(Fextr_map_expl_fles, os.path.abspath(map_expl_out))
            elif mp == 'Fgenick_map':
                append_if_file_exist(Fgenick_map_expl_fles, os.path.abspath(map_expl_out))
            elif mp == 'Fextr_calc_map':
                append_if_file_exist(Fextr_calc_map_expl_fles, os.path.abspath(map_expl_out))
            print("m%s-DFcalc map explored. Results in %s and residue list in associated residlist"
                               % (Fextr.maptype, map_expl_out))
            print("m%s-DFcalc map explored. Results in %s and residue list in associated residlist"
                               % (Fextr.maptype, map_expl_out), file=log)

            # In case of running in slow_and_rigorous / slow_and_curious (whatever you like to call the full way mode):
            # run refinement with phenix or refmac/coot
            if (params.Xtrapol8.f_and_maps.fast_and_furious == False and params.Xtrapol8.refinement.run_refinement):
                print("\n************Refinements************")
                print("\n************Refinements************", file=log)
                if params.Xtrapol8.refinement.use_refmac_instead_of_phenix:
                    mtz_out, pdb_rec, pdb_real, pdb_rec_real = Fextr.refmac_coot_refinements(pdb_in=DH.pdb_in,
                                                                                             additional=DH.additional,
                                                                                             ligands_list=DH.extract_ligand_codes(),
                                                                                             F_column_labels=
                                                                                             Fextr.FM.labels['data'],
                                                                                             map_column_labels='%s, PHI%s, %s, PHI%s'
                                                                                                               % (
                                                                                                                   Fextr.FM.labels[
                                                                                                                       'map_coefs_map'],
                                                                                                                   Fextr.FM.labels[
                                                                                                                       'map_coefs_map'],
                                                                                                                   Fextr.FM.labels[
                                                                                                                       'map_coefs_diff'],
                                                                                                                   Fextr.FM.labels[
                                                                                                                       'map_coefs_diff']),
                                                                                             keywords=params.Xtrapol8.refinement.refmac_keywords)
                else:
                    mtz_out, pdb_rec, pdb_real, pdb_rec_real = Fextr.phenix_phenix_refinements(pdb_in=DH.pdb_in,
                                                                                               additional=DH.additional,
                                                                                               F_column_labels=
                                                                                               Fextr.FM.labels['data'],
                                                                                               column_labels='%s,PHI%s' % (
                                                                                                   Fextr.FM.labels[
                                                                                                       'map_coefs_map'],
                                                                                                   Fextr.FM.labels[
                                                                                                       'map_coefs_map']),
                                                                                               keywords=params.Xtrapol8.refinement.phenix_keywords)
                print("--------------", file=log)
                # depending on the map-type, append the refinement output to the correct list
                # this is ugly, TODO: make an object to store the results in a clean and transparant way
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
                if mp == 'kFextr_map':
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

            # fast_and_furious mode: no refinement, but need to keep track of structure factor files
            elif params.Xtrapol8.f_and_maps.fast_and_furious:
                append_if_file_exist(Fextr_mtz_lst, os.path.abspath(Fextr.F_name))

            print("---> Results in %s" % (os.getcwd()))
            print("------------------------------------")
            print("---> Results in %s" % (os.getcwd()), file=log)
            print("------------------------------------", file=log)

        ################################################################
        # remove empty directories to avoid any confusion. secure because os.rmdir can only remove empty directories.
        if len(os.listdir(new_dirpath_q)) == 0:
            os.rmdir(new_dirpath_q)
        if len(os.listdir(new_dirpath)) == 0:
            os.rmdir(new_dirpath)

        # Go back to output directory and generate the plots from the pickle files
        os.chdir(outdir)
        plot_Fextr_sigmas()
        # plot_sigmas(maptype_lst=list(map(lambda x: re.sub(r"\_map$", "", x), final_maptypes)))
        plot_negative_reflections()

    ################################################################
    # Repeat the last step in order to make sure we have all plots correctly
    # Go back to output directory and generate the plots from the pickle files
    os.chdir(outdir)
    plot_Fextr_sigmas()
    # plot_sigmas(maptype_lst=list(map(lambda x: re.sub(r"\_map$", "", x), final_maptypes)))
    plot_negative_reflections()

    # plot_correlations(params.occupancies.list_occ, cc_list)

    print('-----------------------------------------')
    print("CALCULATE (Q/K-WEIGHTED) FEXTRAPOLATED MAPS DONE")
    print('-----------------------------------------')

    ################################################################

    print('-----------------------------------------')
    print("ESTIMATE OPTIMAL OCCUPANCY")
    print('-----------------------------------------')
    print('-----------------------------------------', file=log)
    print("ESTIMATE OPTIMAL OCCUPANCY", file=log)
    print('-----------------------------------------', file=log)

    # To estimate the optimal occupancy the integrated difference map peaks can be used and/or the structural movements in real space refinement (in slow_and_rigorous mode only)
    # In both cases, a file with residues to base the calculation on should be provided
    #   In case of faf, the Z-score based list will automatically be choses
    #   In sor mode, the user has the possibility to provide a list during a pauze of 5 minutes, but if nothong is provided or file is not found, the Z-score based list will be used.
    if params.Xtrapol8.f_and_maps.fast_and_furious:
        print(
            "BASED ON THE DIFFERENCE MAPS (Fo-Fo AND mFextr-DFc) AND THE THRESHOLD, RADIUS AND PEAK PARAMETERS YOU PROVIDED, WE SELECTED RESIDUES THAT ARE NEAR THE PEAKS. DEPENDING ON THE Z-SCORE WE SELECTED ONLY THE RESIDUES AROUND THE HIGHEST PEAKS. YOU'RE IN FAST AND FURIOUS MODE SO I WILL USE THESE. THIS ONLY WORKS WELL IF YOUR MAPEXPLORER PARAMETERS ARE WELL CHOSEN.")
        residlst = residlist_zscore
    else:
        print(
            "BASED ON THE DIFFERENCE MAPS AND THE THRESHOLD, RADIUS AND PEAK PARAMETERS, WE SELECTED RESIDUES THAT ARE NEAR THE PEAKS ('residlist.txt') DEPENDING ON THE Z-SCORE WE ALSO SELECTED ONLY THE RESIDUES AROUND THE HIGHEST PEAKS ('residlist_zscore.txt'). PLEASE INSPECT THE Fo-Fo DIFFERENCE MAP, MODIFY ONE OF THE RESIDUE ACCORDING TO THOSE RESIDUES THAT YOU BELIEVE HAVE REAL RELEVANT DIFFERENCE PEAKS AND STORE IT UNDER A NEW NAME.")
        if params.Xtrapol8.map_explorer.use_occupancy_from_distance_analysis:
            print(
                "WE WILL USE THE DIFFERENCE OF THE REAL-SPACE REFINED MODEL RESIDUES IN THE LIST TO SUGGEST A PROPER ALPHA-VALUE/OCCUPANCY OF THE TRIGGERED STATE. IF YOU DO NOT PROVIDE ANYTHING, THE Z-SCORE BASED LIST WILL BE USED AND THE RESULTING ALPHA-VALUE/OCCUPANCY MIGHT INFLUENCED BY ARTEFACTS.")
        else:
            print(
                "WE WILL USE THE DIFFERENCE MAPS PEAKS AROUND THOSE RESIDUES TO SUGGEST A PROPER ALPHA-VALUE/OCCUPANCY OF THE TRIGGERED STATE. IF YOU DO NOT PROVIDE ANYTHING, THE Z-SCORE BASED LIST WILL BE USED AND THE RESULTING ALPHA-VALUE/OCCUPANCY MIGHT INFLUENCED BY ARTEFACTS.")

        timeout = 300  # Time-out of 5 minutes to provide a different resid list than the z-score based list
        rlist, _, _ = select([sys.stdin], [], [], timeout)
        if rlist:
            s = sys.stdin.readline()
            residlst = s.strip("\n")
        else:
            residlst = residlist_zscore

        if os.path.isfile(residlst.lstrip().rstrip()):
            residlst = residlst.lstrip().rstrip()
        elif os.path.isfile(startdir + "/" + residlst.lstrip().rstrip()):
            residlst = startdir + "/" + residlst.lstrip().rstrip()
        else:
            residlst = residlist_zscore
        try:
            residlst = os.path.abspath(residlst)
            if os.path.isfile(residlst) == False:
                raise IOError
        except IOError:
            print("File %s not found. All peaks and residues will be used" % (str(residlst)))
            residlst = None
    print("Residue list used for estimation of occupancy of triggered state: %s\n" % (residlst), file=log)
    print("Residue list used for estimation of occupancy of triggered state: %s\n" % (residlst))

    ################################################################
    # occupancy estimation will be performed for all map types, but the final automatic descision will be based on a priority list.
    # This means that if maptypes Fextr and qFextr are selected, the output of qFextr will have priority.
    # Therefore, here we loop over the maptypes in a reversed order to
    # 1) transfer the stored output files to a generic list
    # 2) estimate occupacy based on the peakintegration volumes
    # If slow_and_rigorous:
    #   3) plot the refinement R-factors
    #   4) estimate the occupancy based on the distances of the real space refined models with the input model
    #   5) append the refinement models and maps to the Pymol_movie script
    #   6) make ddm plot
    #   7) write coot script
    reversed_maptypes = P.final_maptypes[:]
    reversed_maptypes.reverse()
    occ_overview = {}
    for mp in reversed_maptypes:
        mp_type = mp.split("_map")[0]
        print("---%s---" % (mp_type))
        print("---%s---" % (mp_type), file=log)

        if mp in ('qFextr_map', 'qFgenick_map', 'qFextr_calc_map'):
            dir_prefix = 'qweight_occupancy'
        elif mp in ('kFextr_map', 'kFgenick_map', 'kFextr_calc_map'):
            dir_prefix = 'kweight_occupancy'
        else:
            dir_prefix = 'occupancy'

        if mp == 'qFextr_map':
            recref_mtz_lst = qFextr_recref_mtz_lst
            recref_pdb_lst = qFextr_recref_pdb_lst
            recrealref_lst = qFextr_recrealref_lst
            realref_lst = qFextr_realref_lst
            map_expl_lst = qFextr_map_expl_fles
        elif mp == 'qFgenick_map':
            recref_mtz_lst = qFgenick_recref_mtz_lst
            recref_pdb_lst = qFgenick_recref_pdb_lst
            recrealref_lst = qFgenick_recrealref_lst
            realref_lst = qFgenick_realref_lst
            map_expl_lst = qFgenick_map_expl_fles
        elif mp == 'qFextr_calc_map':
            recref_mtz_lst = qFextr_calc_recref_mtz_lst
            recref_pdb_lst = qFextr_calc_recref_pdb_lst
            recrealref_lst = qFextr_calc_recrealref_lst
            realref_lst = qFextr_calc_realref_lst
            map_expl_lst = qFextr_calc_map_expl_fles
        if mp == 'kFextr_map':
            recref_mtz_lst = kFextr_recref_mtz_lst
            recref_pdb_lst = kFextr_recref_pdb_lst
            recrealref_lst = kFextr_recrealref_lst
            realref_lst = kFextr_realref_lst
            map_expl_lst = kFextr_map_expl_fles
        elif mp == 'kFgenick_map':
            recref_mtz_lst = kFgenick_recref_mtz_lst
            recref_pdb_lst = kFgenick_recref_pdb_lst
            recrealref_lst = kFgenick_recrealref_lst
            realref_lst = kFgenick_realref_lst
            map_expl_lst = kFgenick_map_expl_fles
        elif mp == 'kFextr_calc_map':
            recref_mtz_lst = kFextr_calc_recref_mtz_lst
            recref_pdb_lst = kFextr_calc_recref_pdb_lst
            recrealref_lst = kFextr_calc_recrealref_lst
            realref_lst = kFextr_calc_realref_lst
            map_expl_lst = kFextr_calc_map_expl_fles
        elif mp == 'Fextr_map':
            recref_mtz_lst = Fextr_recref_mtz_lst
            recref_pdb_lst = Fextr_recref_pdb_lst
            recrealref_lst = Fextr_recrealref_lst
            realref_lst = Fextr_realref_lst
            map_expl_lst = Fextr_map_expl_fles
        elif mp == 'Fgenick_map':
            recref_mtz_lst = Fgenick_recref_mtz_lst
            recref_pdb_lst = Fgenick_recref_pdb_lst
            recrealref_lst = Fgenick_recrealref_lst
            realref_lst = Fgenick_realref_lst
            map_expl_lst = Fgenick_map_expl_fles
        elif mp == 'Fextr_calc_map':
            recref_mtz_lst = Fextr_calc_recref_mtz_lst
            recref_pdb_lst = Fextr_calc_recref_pdb_lst
            recrealref_lst = Fextr_calc_recrealref_lst
            realref_lst = Fextr_calc_realref_lst
            map_expl_lst = Fextr_calc_map_expl_fles

            # Estimate alpha and occupancy based on the peakintegration area as stored in the peakintegration files
        alpha, occ = plotalpha(map_expl_lst, residlst, mp_type, log=log).estimate_alpha()

        if (params.Xtrapol8.f_and_maps.fast_and_furious == False and params.Xtrapol8.refinement.run_refinement):
            # If water molecules are updated during refinement, waters will be added and removed and their numbers are not in relation to the original waters in the input model
            #   therefore keeping them during the distance analysis would make no sense
            if params.Xtrapol8.refinement.phenix_keywords.main.ordered_solvent:
                distance_use_waters = False
            else:
                distance_use_waters = True

            # Make a plot of the refinement R-factors, related to the specific maptype. The log-files should have the same prefix as the mtz-files.
            # This assumption is made in order to avoid storing the log-files in even another list
            plot_Rfactors_per_alpha(map(lambda fle: re.sub(r'mtz$', 'log', fle), recref_mtz_lst), mp_type)
            print("", file=log)
            print("")

            # Check if the list with PDB files from real space refinement after reciprocal space refinement is complete
            #   If this is not the case (when reciprocal space or real space refinement had failed),
            #   the PDB files from real space refinment (direct real space refinement in the extrapolated maps) are used
            if ((len(recrealref_lst) == len(set(recrealref_lst)))):
                pdb_list = recrealref_lst
            else:
                pdb_list = realref_lst

            # Estimate alpha and occupancy based on the distances of the list with PDB files with the input model.
            # Overwrite the earlier found alpha and occupancy if use_occupancy_from_distance_analysis is True
            if params.Xtrapol8.map_explorer.use_occupancy_from_distance_analysis:
                alpha, occ = Distance_analysis(pdb_list, params.Xtrapol8.occupancies.list_occ, resids_lst=residlst,
                                               use_waters=distance_use_waters, outsuffix=mp_type,
                                               log=log).extract_alpha()
            else:
                _, _ = Distance_analysis(pdb_list, params.Xtrapol8.occupancies.list_occ, resids_lst=residlst,
                                         use_waters=distance_use_waters, outsuffix=mp_type, log=log).extract_alpha()
            # print("---------", file=log)
            print("", file=log)
            print("")

            # add models and maps to pymol movie.
            # Not elegant with all the with hard coded regular expressions but avoids storing other lists
            pymol_pdb_list = pdb_list[:]
            pymol_pdb_list.remove(DH.pdb_in)
            pymol_mtz_list = recref_mtz_lst[:]
            if P.outname == 'triggered':  # if dummy name applied, the files still contain the dummy name
                pymol_mtz_list = map(lambda fle: re.sub(r"triggered", params.output.outname, fle), pymol_mtz_list)
                pymol_pdb_list = map(lambda fle: re.sub(r"triggered", params.output.outname, fle), pymol_pdb_list)
            # Make Pymol movie with the reciprocal space refined maps if recrealref_lst is complete
            # Otherwise use the real space refined models + direct maps
            if pdb_list == recrealref_lst:
                ccp4_list = map(lambda fle: re.sub(r".mtz$", "_2mFo-DFc_filled.ccp4", fle), pymol_mtz_list)
                model_label = '%s_reciprocal_real_space' % (mp_type)
                ccp4_map_label = '%s_reciprocal_space' % (mp)
                # Pymol_movie(params.occupancies.list_occ, pdblst=pymol_pdb_list, ccp4_maps = ccp4_list, resids_lst = residlst, model_label='%s_reciprocal_real_space'%(mp_type), ccp4_map_label='%s_reciprocal_space'%(mp)).write_pymol_script()
            else:
                if mp == 'qFgenick_map':
                    ccp4_list = map(
                        lambda fle: re.search("(.+?)2mqFgenick-DFc_reciprocal", fle).group(1) + "mqFgenick-DFc.ccp4",
                        pymol_mtz_list)
                elif mp == 'kFgenick_map':
                    ccp4_list = map(
                        lambda fle: re.search("(.+?)2mkFgenick-DFc_reciprocal", fle).group(1) + "mkFgenick-DFc.ccp4",
                        pymol_mtz_list)
                elif mp == 'Fgenick_map':
                    ccp4_list = map(
                        lambda fle: re.search("(.+?)2mFgenick-DFc_reciprocal", fle).group(1) + "mFgenick-DFc.ccp4",
                        pymol_mtz_list)
                else:
                    ccp4_list = map(lambda fle: re.search("(.+?)\_reciprocal", fle).group(1) + ".ccp4", pymol_mtz_list)
                model_label = '%s_real_space' % (mp_type)
                ccp4_map_label = '%s' % (mp)
            if len(ccp4_list) == len(pymol_pdb_list) == len(params.Xtrapol8.occupancies.list_occ):
                Pymol_movie(params.Xtrapol8.occupancies.list_occ, pdblst=pymol_pdb_list, ccp4_maps=ccp4_list,
                            resids_lst=residlst, model_label=model_label,
                            ccp4_map_label=ccp4_map_label).write_pymol_script()
            else:
                Pymol_movie(params.Xtrapol8.occupancies.list_occ, pdblst=pymol_pdb_list, resids_lst=residlst,
                            model_label=model_label).write_pymol_script()

            # If estimated occupancy if not in list (will be case when using distance analysis or when plotalpha fails), take the closest occupancy from the list
            if occ not in params.Xtrapol8.occupancies.list_occ:
                occ = min(params.Xtrapol8.occupancies.list_occ, key=lambda x: abs(x - occ))
                alpha = 1 / occ
            occ_dir = "%s/%s_%.3f" % (outdir, dir_prefix, occ)

            # ddm calculation
            pdb_for_ddm = pdb_list[params.Xtrapol8.occupancies.list_occ.index(occ) + 1]
            print("----Generate distance difference plot----")
            Difference_distance_analysis(DH.pdb_in, pdb_for_ddm, ligands=DH.extract_ligand_codes(), outdir=occ_dir,
                                         scale=params.Xtrapol8.output.ddm_scale).ddms()
            print('---------------------------')

            # Coot script
            mtzs_for_coot = []
            if len(recref_mtz_lst) != len(params.Xtrapol8.occupancies.list_occ):
                occ_list_mtz = [(re.search(r'occupancy\_(.+?)\/', mtz).group(1)) for mtz in recref_mtz_lst]
                mtz_rec = recref_mtz_lst[occ_list_mtz.index(occ)]
            else:
                mtz_rec = recref_mtz_lst[params.Xtrapol8.occupancies.list_occ.index(occ)]
            append_if_file_exist(mtzs_for_coot, mtz_rec)
            if params.Xtrapol8.refinement.phenix_keywords.density_modification.density_modification:
                mtz_dm = re.sub(".mtz$", "_densitymod.mtz", mtz_rec)
                append_if_file_exist(mtzs_for_coot, mtz_dm)

            if params.Xtrapol8.refinement.refmac_keywords.density_modification.density_modification:
                mtz_dm = re.sub(".mtz$", "_dm.mtz", mtz_rec)
                append_if_file_exist(mtzs_for_coot, mtz_dm)
            mtz_extr = ["%s/%s" % (occ_dir, fle) for fle in os.listdir(occ_dir) if
                        P.outname in fle and fle.endswith('m%s-DFc.mtz' % (mp_type))][0]
            append_if_file_exist(mtzs_for_coot, os.path.abspath(mtz_extr))

            pdbs_for_coot = [DH.pdb_in,
                             recref_pdb_lst[params.Xtrapol8.occupancies.list_occ.index(occ) + 1],
                             realref_lst[params.Xtrapol8.occupancies.list_occ.index(occ) + 1],
                             recrealref_lst[params.Xtrapol8.occupancies.list_occ.index(occ) + 1]]
            # if outname == 'triggered': #if dummy name applied, the files still contain the dummy name
            # mtzs_for_coot = map(lambda fle: re.sub(r"triggered",params.output.outname, fle), mtzs_for_coot)
            # pdbs_for_coot = map(lambda fle: re.sub(r"triggered",params.output.outname, fle), pdbs_for_coot)
            script_coot = open_all_in_coot(outdir + "/" + FoFo.mtz_name, pdbs_for_coot, mtzs_for_coot, DH.additional,
                                           occ_dir, mp_type)

        elif (params.Xtrapol8.f_and_maps.fast_and_furious == False and params.refinement.run_refinement == False):
            # If estimated occupancy if not in list (will be case when using distance analysis or when plotalpha fails), take the closest occupancy from the list
            if occ not in params.Xtrapol8.occupancies.list_occ:
                occ = min(params.Xtrapol8.occupancies.list_occ, key=lambda x: abs(x - occ))
                alpha = 1 / occ
            occ_dir = "%s/%s_%.3f" % (outdir, dir_prefix, occ)
            mtzs_for_coot = []
            mtz_extr = ["%s/%s" % (occ_dir, fle) for fle in os.listdir(occ_dir) if
                        P.outname in fle and fle.endswith('m%s-DFc.mtz' % (mp_type))][0]
            append_if_file_exist(mtzs_for_coot, os.path.abspath(mtz_extr))
            pdbs_for_coot = [DH.pdb_in]
            script_coot = open_all_in_coot(outdir + "/" + FoFo.mtz_name, pdbs_for_coot, mtzs_for_coot, DH.additional,
                                           occ_dir, mp_type)

        else:
            # If estimated occupancy not in list (will be case when using distance analysis or when plotalpha fails), take the closest occupancy from the list
            if occ not in params.Xtrapol8.occupancies.list_occ:
                occ = min(params.Xtrapol8.occupancies.list_occ, key=lambda x: abs(x - occ))
                alpha = 1 / occ
            occ_dir = "%s/%s_%.3f" % (outdir, dir_prefix, occ)

        occ_overview[mp_type] = occ

        print("------------------------------------")
        print("------------------------------------", file=log)

    # Add final lines to Pymol_script
    if os.path.isfile('%s/pymol_movie.py' % (outdir)):
        Pymol_movie(params.Xtrapol8.occupancies.list_occ, resids_lst=residlst).write_pymol_appearance(
            '%s/pymol_movie.py' % (outdir))

    print("Summary of occupancy determination:", file=log)
    print("Map type       Occupancy", file=log)
    for k in occ_overview:
        print("{:<15} {:>5.3f}".format(k, occ_overview[k]), file=log)
    print("-> Optimal occupancy of triggered state %.3f." % (occ), file=log)

    print("OCCUAPNCY DETERMINATION SUMMARY")
    print("Map type       Occupancy")
    for k in occ_overview:
        print("{:<15} {:>5.3f}".format(k, occ_overview[k]))
    print("-> Optimal occupancy of triggered state %.3f." % (occ))

    print('-----------------------------------------')
    print("ESTIMATE OPTIMAL OCCUPANCY DONE")
    print('-----------------------------------------')

    ################################################################
    # Run refinements for chosen occupancy when in faf mode (as no refinement has yet run):
    # 1) move to directory where refinements have to be run
    # 2) find correct input files
    # 3) Run refinements
    # 4) make ddm plot
    # 5) write coot script
    if (params.Xtrapol8.f_and_maps.fast_and_furious and params.Xtrapol8.refinement.run_refinement):
        print('-----------------------------------------')
        print("FAST AND FURIOUS REFINEMENT")
        print('-----------------------------------------')
        print('-----------------------------------------', file=log)
        print("FAST AND FURIOUS REFINEMENT",file=log)
        print('-----------------------------------------',file=log)

        print("----Refinements with automatically found occupancy values----")
        print("----Refinements with automatically found occupancy values----", file=log)
        os.chdir(occ_dir)

        qFext_mtz_F = \
            [fle for fle in Fextr_mtz_lst if 'occupancy_%.3f' % (occ) in fle and fle.lower().endswith('_qfextr.mtz')][0]
        assert os.path.isfile(qFext_mtz_F)
        print("Extrapolated structure factors for reciprocal space refinement :%s" % (qFext_mtz_F))
        print("Extrapolated structure factors for reciprocal space refinement :%s" % (qFext_mtz_F), file=log)

        qFext_mtz_map = re.sub(r'qFextr.mtz$', '2mqFextr-DFc_mqFextr-DFc.mtz', qFext_mtz_F)
        assert os.path.isfile(qFext_mtz_map)
        print("Extrapolated map coefficients for real space refinement: %s\n" % (qFext_mtz_map))
        print("Extrapolated map coefficients for real space refinement: %s\n" % (qFext_mtz_map),file=log)

        Fextr.name_out = "%s_occ%.3f" % (P.outname, occ)

        if params.Xtrapol8.refinement.use_refmac_instead_of_phenix:
            mtz_out, pdb_rec, pdb_real, pdb_rec_real = Fextr.refmac_coot_refinements(mtz_F=qFext_mtz_F,
                                                                                     mtz_map=qFext_mtz_map,
                                                                                     pdb_in=DH.pdb_in,
                                                                                     additional=DH.additional,
                                                                                     ligands_list=DH.extract_ligand_codes(),
                                                                                     F_column_labels=Fextr.FM.labels[
                                                                                         'data'],
                                                                                     map_column_labels='%s, PHI%s, %s, PHI%s'
                                                                                                       % (
                                                                                                           Fextr.FM.labels[
                                                                                                               'map_coefs_map'],
                                                                                                           Fextr.FM.labels[
                                                                                                               'map_coefs_map'],
                                                                                                           Fextr.FM.labels[
                                                                                                               'map_coefs_diff'],
                                                                                                           Fextr.FM.labels[
                                                                                                               'map_coefs_diff']),
                                                                                     keywords=params.Xtrapol8.refinement.refmac_keywords)
        else:
            mtz_out, pdb_rec, pdb_real, pdb_rec_real = Fextr.phenix_phenix_refinements(mtz_F=qFext_mtz_F,
                                                                                       mtz_map=qFext_mtz_map,
                                                                                       pdb_in=DH.pdb_in,
                                                                                       additional=DH.additional,
                                                                                       F_column_labels=Fextr.FM.labels[
                                                                                           'data'],
                                                                                       column_labels='%s,PHI%s' % (
                                                                                           Fextr.FM.labels['map_coefs_map'],
                                                                                           Fextr.FM.labels[
                                                                                               'map_coefs_map']),
                                                                                       keywords=params.Xtrapol8.refinement.phenix_keywords)

        print("---> Results in %s_%.3f" % (dir_prefix, occ))
        print("---> Results in %s_%.3f" % (dir_prefix, occ), file=log)

        # Coot script
        mtzs_for_coot = [os.path.abspath(qFext_mtz_map), os.path.abspath(mtz_out)]
        if params.Xtrapol8.refinement.phenix_keywords.density_modification.density_modification:
            mtz_dm = re.sub(".mtz$", "_densitymod.mtz", mtz_out)
            if os.path.isfile(mtz_dm):
                mtzs_for_coot.append(os.path.abspath(mtz_dm))
        if params.Xtrapol8.refinement.phenix_keywords.density_modification.density_modification:
            mtz_dm = re.sub(".mtz$", "_dm.mtz", mtz_out)
            if os.path.isfile(mtz_dm):
                mtzs_for_coot.append(os.path.abspath(mtz_dm))
        append_if_file_exist(mtzs_for_coot, os.path.abspath(qFext_mtz_F))
        pdbs_for_coot = [DH.pdb_in, os.path.abspath(check_file_existance(pdb_rec)),
                         os.path.abspath(check_file_existance(pdb_real)),
                         os.path.abspath(check_file_existance(pdb_rec_real))]
        # if outname == 'triggered': #if dummy name applied, the files still contain the dummy name
        # mtzs_for_coot = map(lambda fle: re.sub(r"triggered",params.output.outname, fle), mtzs_for_coot)
        # pdbs_for_coot = map(lambda fle: re.sub(r"triggered",params.output.outname, fle), pdbs_for_coot)
        script_coot = open_all_in_coot(outdir + "/" + FoFo.mtz_name, pdbs_for_coot, mtzs_for_coot, DH.additional,
                                       occ_dir, "qFextr")

        print('-----------------------------------------')
        print("FAST AND FURIOUS REFINEMENT DONE")
        print('-----------------------------------------')
        print("----Generate distance difference plot----")
        Difference_distance_analysis(DH.pdb_in, pdb_rec_real, ligands=DH.extract_ligand_codes(), outdir=occ_dir,
                                     scale=params.Xtrapol8.output.ddm_scale).ddms()
        print('---------------------------')

    ################################################################
    # Make sure we are in the output directory
    if os.getcwd() != outdir:
        os.chdir(outdir)

    # change names to real output name in case the dummy name was used
    if P.outname == 'triggered':
        print('---------------------------', file=log)
        # print("Replacing the dummpy outname ('%s') with true outname('%s')" %(outname, params.output.outname), file=log)
        print("Replacing the dummpy outname ('%s') with true outname('%s')" % (P.outname, params.output.outname))
        tr = [os.path.join(root, fle) for root, dirs, files in os.walk(outdir) for fle in files if P.outname in fle]
        _ = [os.rename(fle, fle.replace(P.outname, params.output.outname)) for fle in tr]
        FoFo.mtz_name = FoFo.mtz_name.replace(P.outname, params.output.outname)

        coot_scripts = [os.path.join(root, fle) for root, dirs, files in os.walk(outdir) for fle in files if
                        fle.startswith("coot_all_")]
        for coot_script in coot_scripts:
            with open(coot_script, "r") as f:
                fle = f.read().split("\n")
            o = open(coot_script, "w")
            for lne in fle:
                if P.outname in lne:
                    lne = re.sub(P.outname, params.output.outname, lne)
                    o.write("%s\n" % lne)
                else:
                    o.write("%s\n" % lne)
            o.close()

    # Make sure we are in the output directory
    if os.getcwd() != outdir:
        os.chdir(outdir)
    # Make a list of all temporal files. Upon their removal disk space can be saved
    list_redundant_files(outdir)

    #################################################################
    # print('-----------------------------------------', file=log)
    # print("WRITE PYMOL AND COOT SCRIPTS", file=log)
    # print('-----------------------------------------', file=log)

    ##The pymol movie script is already made
    # print("pymol movie with models and maps in %s/pymol_movie.py. Run from terminal or pymol command line"%(outdir), file=log)

    ##The pymol script with just all models and maps is still missing
    ##This script is not working properly and models and maps are in wrong session subfolers or even missing
    # pymol_script = Pymol_visualization(DH.pdb_in, outdir).open_all_in_pymol()
    # print("pymol script with all models written in %s/%s. Can be run from a terminal with 'Pymol %s/%s' if Pymol is in your PATH or 'run %s/%s'from a pymol command line"%(outdir, pymol_script, outdir, pymol_script, outdir, pymol_script), file=log)

    # print('-----------------------------------------', file=log)
    # print("WRITE PYMOL AND COOT SCRIPTS DONE", file=log)
    # print('-----------------------------------------', file=log)

    ################################################################
    # Write a phil file. Parameters changed during the excecution of Xtrapol8 are possible when problems appeared
    modified_phil = master_phil.format(python_object=params)
    modified_phil.show(out=open("Xtrapol8_out.phil", "w"))

    print('-----------------------------------------')
    print("DONE! PLEASE INSPECT THE OUTPUT AND RERUN WITH DIFFERENT PARAMETERS IF NECESSARY.")
    print('-----------------------------------------')
    print('-----------------------------------------', file=log)
    print("DONE! PLEASE INSPECT THE OUTPUT AND RERUN WITH DIFFERENT PARAMETERS IF NECESSARY.", file=log)
    print('-----------------------------------------', file=log)
    if (params.Xtrapol8.f_and_maps.fast_and_furious == False and params.Xtrapol8.refinement.run_refinement):
        print("Pymol script with models and maps in %s/pymol_movie.py." % (outdir))
    print(
        "Coot script with models and maps in each output directory associated with estimated occupancy (coot_all_<maptype>.py).")
    print("Pymol script with models and maps in %s/pymol_movie.py." % (outdir), file=log)
    print(
        "Coot script with models and maps in each output directory associated with estimated occupancy (coot_all_<maptype>.py).", file=log)
    log.close()
    if (params.Xtrapol8.output.open_coot and check_program_path('coot')):
        os.system("coot --script %s" % (script_coot))

    ################################################################