# -*- coding: utf-8 -*-
"""
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

Script to extract an fmodel from a pdb file and mtz file. Used to extract phase info after refinement

"""
from __future__ import division, print_function
import sys, os
import argparse
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import math
from scipy.stats import pearsonr

from iotbx.file_reader import any_file
from cctbx import maptbx, miller, crystal, xray
from cctbx.array_family import flex
import mmtbx.f_model
import iotbx.pdb
import iotbx.phil

from Fextr_utils import get_name


master_phil = iotbx.phil.parse("""
input{
    triggered_mtz = None
        .type = path
        .help = Triggered data in mtz or mmcif format (merged).
        .expert_level = 0
    reference_pdb = None
        .type = path
        .help = Reference coordinates in pdb or mmcif format.
        .expert_level = 0
    }
occupancies{
    list_occ = None
        .type = floats(size_min=1, value_min=0, value_max=1)
        .help = List of occupancies to test (fractional). Will overwrite low_occ, high_occ and steps if defined.
        .expert_level = 0
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
    negative_and_missing = *truncate_and_fill truncate_no_fill fref_and_fill fref_no_fill fcalc_and_fill fcalc_no_fill keep_and_fill keep_no_fill reject_and_fill reject_no_fill zero_and_fill zero_no_fill fill_missing no_fill
        .type = choice(multi=False)
        .help = Handling of negative and missing extrapolated structure factor amplitudes (ESFAs) Note that this will not be applied on the Fourier difference map. If selected, filling of missing reflections is only carried on maps of the 2mFextr-DFcalc type. This parameters is NOT applicable for (q/k)Fgenick because negative reflections are rejected anyway. For refinement, default phenix.refine or refmac handling of negative/missing reflections is applied. keep_no_fill maps will be calculated in addition in all cases. keep_and_fill and keep_no_fill replace the old fill_missing and no_fill arguments which will become invalid keywords in future Xtrapol8 versions. Please check the manual for more information.
        .expert_level = 2
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
    }
    """)

class model(object):
    def __init__(self, pdb):
        self.pdb    = pdb
        self.prefix = get_name(self.pdb)
        self.fmodel = None
        self.fmodel_update = None

    def get_F(self, hkl):
        """
        Extract structure factor amplitudes
        """
        amplitude_types = [xray.observation_types.amplitude,
                        xray.observation_types.reconstructed_amplitude]
        intensity_types = [xray.observation_types.intensity,]

        f_obs =  None
        intensities_found = False

        for array in hkl.file_object.as_miller_arrays():
            if type(array.observation_type()) in amplitude_types:
                f_obs = array
                labels = f_obs.info().labels
                print("Found F's: %s" %(labels))
                break
            if type(array.observation_type()) in intensity_types:
                intensities_found = True

        if f_obs == None: 
            if intensities_found:
                print("No structure factors found. Please convert intensities to structure factors.")
            else:
                print("Could not find F's nor I's, please check the input file.")
                
        return f_obs

    def extract_f_model(self, mtz):
        """
        make an Fmodel based on an mtz file containing structure factor ampltitudes and pdb file
        """
        
        #read files
        reflections = any_file(mtz, force_type="hkl", raise_sorry_if_errors=True)
        model_in = any_file(self.pdb, force_type="pdb", raise_sorry_if_errors=True)
        
        #extrcat the data from the mtz file
        f_obs = self.get_F(reflections)
        if f_obs == None:
            print("Failed to extract structure factor amplitudes from {:s}".format(mtz))
            return
        print("\n------------------------------------")
        f_obs.show_summary()
        print("\n------------------------------------")
        
        #make new rfree array. Required to make f_model and easier than extracting from the input file
        rfree = f_obs.generate_r_free_flags(fraction=0.05)
        
        #get x-ray structure from input model. Required to make f_model
        pdb_ini = iotbx.pdb.input(self.pdb)
        xray_structure = pdb_ini.xray_structure_simple()
        
        #make f_model
        self.fmodel = mmtbx.f_model.manager(
                f_obs          = f_obs,
                r_free_flags   = rfree,
                xray_structure = xray_structure)
        
        #update the fmodel
        print("Updating scales")
        self.fmodel_update = self.fmodel.deep_copy()
        self.fmodel_update.update_all_scales(remove_outliers=False, log=sys.stdout)
        
        #give information about the updated fmodel
        print("fmodel without update of scales:")
        self.fmodel.show()
        print("\n------------------------------------")
        print("fmodel with update of scales:")
        self.fmodel_update.show()
        print("\n------------------------------------")
        
def plot_phases_and_fom(fmodel,
                        fmodel_prefix = "fmodel",
                        prefix = "Phase_info"):
    """
    plot the phases and fom for the updated and 
    """
    print("Extract and plot phase information: {:s}".format(fmodel_prefix))
    
    fmodel_phases = fmodel.f_model().phases()
    fmodel_phase_errors = flex.double(np.radians(fmodel.phase_errors()))
    fmodel_fom = fmodel.fom()
    fmodel_indices = fmodel_phases.indices()
    fmodel_alpha = fmodel.alpha_beta()[0]
    
    fmodel_bin_res_cent_lst     = []
    fmodel_phases_bin_lst       = []
    fmodel_phase_errors_bin_lst = []
    fmodel_fom_bin_lst          = []
    fmodel_alpha_bin_lst        = []
    
    fmodel_prefix_alt = "model1"
    print("{:s}: {:s}".format(fmodel_prefix_alt, fmodel_prefix))
    
    fmodel_phases.setup_binner(n_bins=20)
    print("bin  resolution range  #reflections        <fom>   <alpha> <phase> <phase-error>")
    for i_bin in fmodel_phases.binner().range_all():
        #get info from fmodel
        sel_fmodel = fmodel_phases.binner().selection(i_bin)
        fmodel_phases_bin = fmodel_phases.select(sel_fmodel).data()
        if fmodel_phases_bin.size() == 0 : continue
        fmodel_phases_bin_av = np.average(fmodel_phases_bin)
        fmodel_phases_bin_lst.append(fmodel_phases_bin_av)
        fmodel_phase_errors_bin = fmodel_phase_errors.select(sel_fmodel)
        fmodel_phase_errors_bin_av = np.average(fmodel_phase_errors_bin)
        fmodel_phase_errors_bin_lst.append(fmodel_phase_errors_bin_av)
        fmodel_fom_bin = fmodel_fom.select(sel_fmodel).data()
        fmodel_fom_bin_av = np.average(fmodel_fom_bin)
        fmodel_fom_bin_lst.append(fmodel_fom_bin_av)
        fmodel_alpha_bin = fmodel_alpha.select(sel_fmodel).data()
        fmodel_alpha_bin_av = np.average(fmodel_alpha_bin)
        fmodel_alpha_bin_lst.append(fmodel_alpha_bin_av)
        fmodel_bin_res_cent = np.average(fmodel_phases.binner().bin_d_range(i_bin))
        fmodel_bin_res_cent_lst.append(fmodel_bin_res_cent)
        #print info
        legend = fmodel_phases.binner().bin_legend(i_bin, show_counts=False)
        print("{:s} {:^10d} {:> 12.4f} {:> 6.4f} {:> 6.4f} {:> 6.4f}".format(legend, sel_fmodel.size(), fmodel_fom_bin_av, fmodel_alpha_bin_av, fmodel_phases_bin_av, fmodel_phase_errors_bin_av))

    fig, axs = plt.subplots(3, 2, figsize=(10, 10), squeeze=False)
    
    #plot the updated and f_model_fom
    #in function of reflection
    axs[0,0].scatter(range(fmodel_fom.data().size()), fmodel_fom.data(), color = "tab:blue", label="{:s} fom.\n Average: {:.2f}".format(fmodel_prefix_alt, np.mean(fmodel_fom.data())), marker = ".", zorder=0, alpha=0.1)
    axs[0,0].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    axs[0,0].plot(range(fmodel_fom.data().size()), flex.double(fmodel_fom.data().size(), np.mean(fmodel_fom.data())), color="blue", zorder=3)
    axs[0,0].set_xlabel("Reflection")
    axs[0,0].set_ylabel("fom")
    axs[0,0].set_ylim(-0.05, 1.05)
    axs[0,0].ticklabel_format(style="sci")
    
    #in function of resolution
    axs[0,1].scatter(fmodel_bin_res_cent_lst, fmodel_fom_bin_lst, color = "tab:blue", label="{:s} <fom>.".format(fmodel_prefix_alt), marker = ".", zorder=0)
    axs[0,1].set_xlim(np.max(fmodel_bin_res_cent_lst[1:]), np.min(fmodel_bin_res_cent_lst[1:]))
    axs[0,1].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    axs[0,1].set_xlabel("Resolution (A)")
    axs[0,1].set_ylabel("fom")
    axs[0,1].set_ylim(-0.05, 1.05)

    #plot the updated alpha (D)
    #in function of reflection
    axs[1,0].scatter(range(fmodel_alpha.data().size()), fmodel_alpha.data(), color = "tab:blue", label="{:s} alpha (D).\n Average: {:.2f}".format(fmodel_prefix_alt, np.mean(fmodel_alpha.data())), marker = ".", zorder=0, alpha=0.1)
    axs[1,0].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    axs[1,0].plot(range(fmodel_alpha.data().size()), flex.double(fmodel_alpha.data().size(), np.mean(fmodel_alpha.data())), color="blue", zorder=2)
    axs[1,0].set_xlabel("Reflection")
    axs[1,0].set_ylabel("alpha (D)")
    #axs[1,1].set_ylim(0, 1)
    axs[1,0].ticklabel_format(style="sci")
    
    #in function of resolution
    axs[1,1].scatter(fmodel_bin_res_cent_lst, fmodel_alpha_bin_lst, color = "tab:blue", label="{:s} <alpha (D)>.".format(fmodel_prefix_alt), marker = ".", zorder=0)
    axs[1,1].set_xlim(np.max(fmodel_bin_res_cent_lst[1:]), np.min(fmodel_bin_res_cent_lst[1:]))
    axs[1,1].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    axs[1,1].set_xlabel("Resolution (A)")
    axs[1,1].set_ylabel("alpha (D)")         
    
    
    #plot the phase and phase error of the f_model
    #in function of reflection
    axs[2,0].scatter(range(fmodel_phases.data().size()), fmodel_phases.data(), color = "tab:blue", label = "{:s} fmodel_phases.\n Average = {:.2f}".format(fmodel_prefix_alt, np.mean(fmodel_phases.data())), marker = ".")
    axs[2,0].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    axs[2,0].plot(range(fmodel_phases.data().size()), flex.double(fmodel_phases.data().size(), np.mean(fmodel_phases.data())), color="blue", zorder=3)
    axs[2,0].set_xlabel("Reflection")
    axs[2,0].set_ylabel("Phase")    
    axs[2,0].set_ylim(-2*math.pi, 2*math.pi)
    axs[2,0]
    ax2 = axs[2,0].twinx()
    ax2.fill_between(range(fmodel_phases.data().size()), (fmodel_phases.data()-fmodel_phase_errors), (fmodel_phases.data()+fmodel_phase_errors), color="tab:blue", alpha=0.2, label='phase error')
    #ax2.tick_params(axis='y')
    ax2.tick_params(top=False, labeltop=False, left=False, labelleft=False, right=False, labelright=False, bottom=False, labelbottom=False)
    ax2.set_ylabel('Phase error')
    ax2.set_ylim(-2*math.pi, 2*math.pi)
    
    #in function of resolution
    axs[2,1].scatter(fmodel_bin_res_cent_lst, fmodel_phases_bin_lst, color = "tab:blue", label = "{:s} <phases>.".format(fmodel_prefix_alt), marker = ".")
    axs[2,1].set_xlim(np.max(fmodel_bin_res_cent_lst[1:]), np.min(fmodel_bin_res_cent_lst[1:]))
    axs[2,1].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    axs[2,1].set_xlabel("Resolution (A)")
    axs[2,1].set_ylabel("Phase")
    #ax3 = axs[2,1].twinx()
    #ax3.fill_between(fmodel_bin_res_cent_lst,  np.subtract(fmodel_phases_bin_lst,fmodel_phase_errors_bin_lst),  np.subtract(fmodel_phases_bin_lst,fmodel_phase_errors_bin_lst), color="tab:red", alpha=0.2, label='phase error')
    #ax3.set_xlim(np.max(fmodel_bin_res_cent_lst[1:]), np.min(fmodel_bin_res_cent_lst[1:]))
    #ax3.tick_params(axis='y')
    #ax3.set_ylabel('Phase error {:s}'.format(fmodel_prefix_alt))
    #ax3.set_ylim(-2*math.pi, 2*math.pi)
    
    fig.suptitle("{:s}: {:s}".format(fmodel_prefix_alt, fmodel_prefix), ha='left', x=0.05)
    
    plt.subplots_adjust(hspace=0.35, wspace=0.5, left=0.09, right=0.88, top = 0.95)
    plt.savefig("{:s}.png".format(prefix))
    plt.close()
    
    print("{:s} fom. Average: {:.2f}".format(fmodel_prefix_alt, np.mean(fmodel_fom.data())))
    print("{:s} alpha (D). Average: {:.2f}".format(fmodel_prefix_alt, np.mean(fmodel_alpha.data())))
    print("{:s} phase. Average: {:.2f}".format(fmodel_prefix_alt, np.mean(fmodel_phases.data())))
    print("\n------------------------------------")
    
def compare_phases_and_fom(fmodel_1,
                           fmodel_2,
                           fmodel_1_prefix = "fmodel",
                           fmodel_2_prefix = "fmodel_2",
                           prefix = "Phase_info_correlation"):
    """
    Correlate the phases and fom of two fmodels
    """
    print("Get correlation between phase information: {:s} and {:s}".format(fmodel_1_prefix, fmodel_2_prefix))
    
    fmodel1_phases = fmodel_1.f_model().phases()
    fmodel1_phase_errors = flex.double(np.radians(fmodel_1.phase_errors()))
    fmodel1_fom = fmodel_1.fom()
    fmodel1_indices = fmodel1_phases.indices()
    fmodel1_alpha = fmodel_1.alpha_beta()[0]

    fmodel2_phases = fmodel_2.f_model().phases()
    fmodel2_phase_errors = flex.double(np.radians(fmodel_2.phase_errors()))
    fmodel2_fom = fmodel_2.fom()
    fmodel2_indices = fmodel2_phases.indices()
    fmodel2_alpha = fmodel_2.alpha_beta()[0]

    
    if fmodel1_fom.data().size() != fmodel2_fom.data().size():
        print("fmodels don't have the same number of Miller indices. Phase information cannot be compared.")
        CC_fom = CC_alpha = CC_phase = None
    elif np.all([fmodel1_indices != fmodel2_indices]):
        print("fmodels don't have the same Miller indices. Phase information cannot be compared.")
        CC_fom = CC_alpha = CC_phase = None
    else:
        fmodel1_bin_res_cent_lst     = []
        fmodel1_phases_bin_lst       = []
        fmodel1_phase_errors_bin_lst = []
        fmodel1_fom_bin_lst          = []
        fmodel1_alpha_bin_lst        = []
        
        fmodel2_bin_res_cent_lst     = []
        fmodel2_phases_bin_lst       = []
        fmodel2_phase_errors_bin_lst = []
        fmodel2_fom_bin_lst          = []
        fmodel2_alpha_bin_lst        = []

        fmodel1_phases.setup_binner(n_bins=20)
        fmodel2_phases.use_binning_of(fmodel1_phases)
        
        fmodel_1_prefix_alt = "model1"
        fmodel_2_prefix_alt = "model2"
        print("{:s}: {:s}\n{:s}: {:s}".format(fmodel_1_prefix_alt, fmodel_1_prefix, fmodel_2_prefix_alt, fmodel_2_prefix))
        print("                                              <fom>              <alpha>             <phase>         <phase error>")
        print("bin  resolution range  #reflections        {:s}  {:s}     {:s}  {:s}     {:s}  {:s}     {:s}  {:s}".format(fmodel_1_prefix_alt, fmodel_2_prefix_alt, fmodel_1_prefix_alt, fmodel_2_prefix_alt, fmodel_1_prefix_alt, fmodel_2_prefix_alt, fmodel_1_prefix_alt, fmodel_2_prefix_alt))

        for i_bin in fmodel1_phases.binner().range_all():
            #get info from fmodel1
            sel_fmodel1 = fmodel1_phases.binner().selection(i_bin)
            fmodel1_phases_bin = fmodel1_phases.select(sel_fmodel1).data()
            if fmodel1_phases_bin.size() == 0 : continue
            fmodel1_phases_bin_av = np.average(fmodel1_phases_bin)
            fmodel1_phases_bin_lst.append(fmodel1_phases_bin_av)
            fmodel1_phase_errors_bin = fmodel1_phase_errors.select(sel_fmodel1)
            fmodel1_phase_errors_bin_av = np.average(fmodel1_phase_errors_bin)
            fmodel1_phase_errors_bin_lst.append(fmodel1_phase_errors_bin_av)
            fmodel1_fom_bin = fmodel1_fom.select(sel_fmodel1).data()
            fmodel1_fom_bin_av = np.average(fmodel1_fom_bin)
            fmodel1_fom_bin_lst.append(fmodel1_fom_bin_av)
            fmodel1_alpha_bin = fmodel1_alpha.select(sel_fmodel1).data()
            fmodel1_alpha_bin_av = np.average(fmodel1_alpha_bin)
            fmodel1_alpha_bin_lst.append(fmodel1_alpha_bin_av)
            fmodel1_bin_res_cent = np.average(fmodel1_phases.binner().bin_d_range(i_bin))
            fmodel1_bin_res_cent_lst.append(fmodel1_bin_res_cent)
            #get info from fmodel1
            sel_fmodel2 = fmodel2_phases.binner().selection(i_bin)
            fmodel2_phases_bin = fmodel2_phases.select(sel_fmodel2).data()
            fmodel2_phases_bin_av = np.average(fmodel2_phases_bin)
            fmodel2_phases_bin_lst.append(fmodel2_phases_bin_av)
            fmodel2_phase_errors_bin = fmodel2_phase_errors.select(sel_fmodel2)
            fmodel2_phase_errors_bin_av = np.average(fmodel2_phase_errors_bin)
            fmodel2_phase_errors_bin_lst.append(fmodel2_phase_errors_bin_av)
            fmodel2_fom_bin = fmodel2_fom.select(sel_fmodel2).data()
            fmodel2_fom_bin_av = np.average(fmodel2_fom_bin)
            fmodel2_fom_bin_lst.append(fmodel2_fom_bin_av)
            fmodel2_alpha_bin = fmodel2_alpha.select(sel_fmodel2).data()
            fmodel2_alpha_bin_av = np.average(fmodel2_alpha_bin)
            fmodel2_alpha_bin_lst.append(fmodel2_alpha_bin_av)
            fmodel2_bin_res_cent = np.average(fmodel2_phases.binner().bin_d_range(i_bin))
            fmodel2_bin_res_cent_lst.append(fmodel2_bin_res_cent)
            
            #print info
            legend = fmodel1_phases.binner().bin_legend(i_bin, show_counts=False)
            print("{:s} {:^10d} {:> 12.4f} {:> 6.4f} {:> 10.4f} {:> 6.4f} {:> 10.4f} {:> 6.4f} {:> 10.4f} {:> 6.4f}".format(legend, sel_fmodel1.size(), fmodel1_fom_bin_av, fmodel2_fom_bin_av, fmodel1_alpha_bin_av, fmodel2_alpha_bin_av, fmodel1_phases_bin_av, fmodel2_phases_bin_av, fmodel1_phase_errors_bin_av, fmodel2_phase_errors_bin_av))
            
        #initiate plot
        fig, axs = plt.subplots(3, 3, figsize=(10, 10))
        
        #plot the updated and f_model_fom
        #in function of reflection
        axs[0,0].scatter(range(fmodel1_fom.data().size()), fmodel1_fom.data(), color = "tab:blue", label="{:s} fom.\n Average: {:.2f}".format(fmodel_1_prefix_alt,np.mean(fmodel1_fom.data())), marker = ".", alpha=0.1, zorder=1)
        axs[0,0].scatter(range(fmodel2_fom.data().size()), fmodel2_fom.data(), color = "tab:red", label="{:s} fom.\n Average: {:.2f}".format(fmodel_2_prefix_alt, np.mean(fmodel2_fom.data())), marker = ".", zorder=0, alpha=0.1)
        axs[0,0].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        axs[0,0].plot(range(fmodel1_fom.data().size()), flex.double(fmodel1_fom.data().size(), np.mean(fmodel1_fom.data())), color="blue", zorder=3)
        axs[0,0].plot(range(fmodel2_fom.data().size()), flex.double(fmodel2_fom.data().size(), np.mean(fmodel2_fom.data())), color="red", zorder=2)
        axs[0,0].set_xlabel("Reflection")
        axs[0,0].set_ylabel("fom")
        axs[0,0].set_ylim(-0.05, 1.05)
        axs[0,0].ticklabel_format(style="sci")
        
        #in function of resolution
        axs[0,1].scatter(fmodel1_bin_res_cent_lst, fmodel1_fom_bin_lst, color = "tab:blue", label="{:s} <fom>" .format(fmodel_1_prefix_alt), marker = ".", zorder=1)
        axs[0,1].scatter(fmodel2_bin_res_cent_lst, fmodel2_fom_bin_lst, color = "tab:red", label="{:s} <fom>".format(fmodel_2_prefix_alt), marker = ".", zorder=0)
        axs[0,1].set_xlim(np.max(fmodel1_bin_res_cent_lst[1:]), np.min(fmodel1_bin_res_cent_lst[1:]))
        axs[0,1].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        axs[0,1].set_xlabel("Resolution (A)")
        axs[0,1].set_ylabel("fom")
        axs[0,1].set_ylim(-0.05, 1.05)
        
        #correlation
        CC_fom = pearsonr(fmodel1_fom.data(), fmodel2_fom.data())[0]
        axs[0,2].scatter(fmodel1_fom.data(), fmodel2_fom.data(), color = "tab:blue", marker = ".", label="fom PearsonR = {:.2f}".format(CC_fom))
        axs[0,2].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        axs[0,2].set_xlabel("{:s} fom".format(fmodel_1_prefix_alt))
        axs[0,2].set_ylabel("{:s} fom".format(fmodel_2_prefix_alt))
        axs[0,2].set_ylim(-0.05, 1.05)
        axs[0,2].set_xlim(-0.05, 1.05)
                
        #plot the updated alpha (D)
        #in function of reflection
        axs[1,0].scatter(range(fmodel1_alpha.data().size()), fmodel1_alpha.data(), color = "tab:blue", label="{:s} alpha (D).\n Average: {:.2f}".format(fmodel_1_prefix_alt, np.mean(fmodel1_alpha.data())), marker = ".", alpha=0.1, zorder=1)
        axs[1,0].scatter(range(fmodel2_alpha.data().size()), fmodel2_alpha.data(), color = "tab:red", label="{:s} alpha (D).\n Average: {:.2f}".format(fmodel_2_prefix_alt, np.mean(fmodel2_alpha.data())), marker = ".", zorder=0, alpha=0.1)
        axs[1,0].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        axs[1,0].plot(range(fmodel1_alpha.data().size()), flex.double(fmodel1_alpha.data().size(), np.mean(fmodel1_alpha.data())), color="blue", zorder=3)
        axs[1,0].plot(range(fmodel2_alpha.data().size()), flex.double(fmodel2_alpha.data().size(), np.mean(fmodel2_alpha.data())), color="red", zorder=2)
        axs[1,0].set_xlabel("Reflection")
        axs[1,0].set_ylabel("alpha (D)")
        #axs[1,1].set_ylim(0, 1)
        axs[1,0].ticklabel_format(style="sci")
        
        #in function of resolution
        axs[1,1].scatter(fmodel1_bin_res_cent_lst, fmodel1_alpha_bin_lst, color = "tab:blue", label="{:s} <alpha (D)> ".format(fmodel_1_prefix_alt), marker = ".", zorder=1)
        axs[1,1].scatter(fmodel2_bin_res_cent_lst, fmodel2_alpha_bin_lst, color = "tab:red", label="{:s} <alpha (D)>".format(fmodel_2_prefix_alt), marker = ".", zorder=0)
        axs[1,1].set_xlim(np.max(fmodel1_bin_res_cent_lst[1:]), np.min(fmodel1_bin_res_cent_lst[1:]))
        axs[1,1].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        axs[1,1].set_xlabel("Resolution (A)")
        axs[1,1].set_ylabel("alpha (D)")

        #correlation
        CC_alpha = pearsonr(fmodel1_alpha.data(), fmodel2_alpha.data())[0]
        axs[1,2].scatter(fmodel1_alpha.data(), fmodel2_alpha.data(), color = "tab:blue", marker = ".", label="alpha (D) PearsonR = {:.2f}".format(CC_alpha))
        axs[1,2].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        axs[1,2].set_xlabel("{:s} alpha (D)".format(fmodel_1_prefix_alt))
        axs[1,2].set_ylabel("{:s} alpha (D)".format(fmodel_2_prefix_alt))
        #axs[1,2].set_ylim(0, 1)
        #axs[1,2].set_xlim(0, 1)
                
        #plot the phase and phase error
        #in function of reflection
        axs[2,0].scatter(range(fmodel1_phases.data().size()), fmodel1_phases.data(), color = "tab:blue", label = "{:s} phases.\n Average = {:.2f}".format(fmodel_1_prefix_alt, np.mean(fmodel1_phases.data())), marker = ".")
        axs[2,0].scatter(range(fmodel2_phases.data().size()), fmodel2_phases.data(), color = "tab:red", label = "{:s} phases.\n Average = {:.2f}".format(fmodel_2_prefix_alt, np.mean(fmodel2_phases.data())), marker = ".")
        axs[2,0].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        axs[2,0].plot(range(fmodel1_phases.data().size()), flex.double(fmodel1_phases.data().size(), np.mean(fmodel1_phases.data())), color="blue", zorder=3)
        axs[2,0].plot(range(fmodel2_phases.data().size()), flex.double(fmodel2_phases.data().size(), np.mean(fmodel2_phases.data())), color="blue", zorder=4)
        axs[2,0].set_xlabel("Reflection")
        axs[2,0].set_ylabel("Phase")
        axs[2,0].set_ylim(-2*math.pi, 2*math.pi)
        axs[2,0].ticklabel_format(axis="x", style="sci")
        ax2 = axs[2,0].twinx()
        ax2.fill_between(range(fmodel1_phases.data().size()), (fmodel1_phases.data()-fmodel1_phase_errors), (fmodel1_phases.data()+fmodel1_phase_errors), color="tab:blue", alpha=0.2, label='phase error')
        ax2.fill_between(range(fmodel2_phases.data().size()), (fmodel2_phases.data()-fmodel2_phase_errors), (fmodel2_phases.data()+fmodel2_phase_errors), color="tab:red", alpha=0.2, label='phase error')
        #ax2.tick_params(axis='y')
        ax2.tick_params(top=False, labeltop=False, left=False, labelleft=False, right=False, labelright=False, bottom=False, labelbottom=False)
        ax2.set_ylabel('Phase error')
        ax2.set_ylim(-2*math.pi, 2*math.pi)
        
        #in function of resolution
        axs[2,1].scatter(fmodel1_bin_res_cent_lst, fmodel1_phases_bin_lst, color = "tab:blue", label = "{:s} <phases>".format(fmodel_1_prefix_alt), marker = ".")
        axs[2,1].scatter(fmodel2_bin_res_cent_lst, fmodel2_phases_bin_lst, color = "tab:red", label = "{:s} <phases>".format(fmodel_2_prefix_alt), marker = ".")
        axs[2,1].set_xlim(np.max(fmodel1_bin_res_cent_lst[1:]), np.min(fmodel1_bin_res_cent_lst[1:]))
        axs[2,1].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        axs[2,1].set_xlabel("Resolution (A)")
        axs[2,1].set_ylabel("Phase")
        #ax3 = axs[2,1].twinx()
        #ax3.fill_between(fmodel_bin_res_cent_lst,  np.subtract(fmodel_phases_bin_lst,fmodel_phase_errors_bin_lst),  np.subtract(fmodel_phases_bin_lst,fmodel_phase_errors_bin_lst), color="tab:red", alpha=0.2, label='phase error')
        #ax3.set_xlim(np.max(fmodel1_bin_res_cent_lst[1:]), np.min(fmodel1_bin_res_cent_lst[1:]))
        #ax3.tick_params(axis='y')
        #ax3.set_ylabel('Phase error {:s}'.format(f_model_label))
        #ax3.set_ylim(-2*math.pi, 2*math.pi)
    
        #correlation
        CC_phase = pearsonr(fmodel1_phases.data(), fmodel2_phases.data())[0]
        axs[2,2].scatter(fmodel1_phases.data(), fmodel2_phases.data(), color = "tab:blue", marker = ".", label="phases PearsonR = {:.2f}".format(CC_phase))
        axs[2,2].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        axs[2,2].set_xlabel("{:s} phase".format(fmodel_1_prefix_alt))
        axs[2,2].set_ylabel("{:s} phase".format(fmodel_2_prefix_alt))
        axs[2,2].set_ylim(0, 1)
        axs[2,2].set_xlim(0, 1)
        
        fig.suptitle("{:s}: {:s}\n{:s}: {:s}".format(fmodel_1_prefix_alt, fmodel_1_prefix, fmodel_2_prefix_alt, fmodel_2_prefix), ha='left', x=0.05)
    
        plt.subplots_adjust(hspace=0.35, wspace=0.65, left=0.09, right=0.88, top = 0.88)
        plt.savefig("{:s}.png".format(prefix))
        plt.close()

        print("{:s} {:s} fom PearsonR = {:.2f}".format(fmodel_1_prefix_alt, fmodel_2_prefix_alt, CC_fom))
        print("{:s} {:s} alpha (d) PearsonR = {:.2f}".format(fmodel_1_prefix_alt, fmodel_2_prefix_alt, CC_alpha))
        print("{:s} {:s} phase PearsonR = {:.2f}".format(fmodel_1_prefix_alt, fmodel_2_prefix_alt, CC_phase))
        
    
    print("\n------------------------------------")    
    return CC_fom, CC_alpha, CC_phase
    
    
def plot_CCs(CC_array, model_names, prefix="Correlations"):
    """
    Plot the evolution of the CCs for the different models
    """
    CC_names = ["CC_fom", "CC_alpha", "CC_phase"]
    
    if len(CC_array.shape) > 1:
        n = CC_array.shape[0]
    else:
        n = len(CC_array.shape)
        
    m = len(model_names)
    
    if n != m:
        print("Provide a single model name per set of CCs.")
        return
        
    if n >= 1:
        fig, axs = plt.subplots(len(CC_names), 1, figsize=(10, 10))
        for i in range(len(CC_names)):
            if n == 1:
                axs[i].scatter(range(1,m+1), CC_array[0][i], color = "tab:blue", marker = "o", label = CC_names[i])
            else:
                axs[i].scatter(range(1,m+1), CC_array[:,i], color = "tab:blue", marker = "o", label = CC_names[i])
            #axs[i].legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
            axs[i].set_xlabel("Model")
            axs[i].set_ylabel(CC_names[i])
            axs[i].xaxis.set_major_locator(MaxNLocator(integer=True))
            
        title = ""
        for i in range(m):
            title += "model{:d}: {:s}\n".format(i+1, model_names[i])
        fig.suptitle(title, ha='left', x=0.05)
        plt.subplots_adjust(hspace=0.35, left=0.09, right=0.88, top = 0.88)
        plt.savefig("{:s}.png".format(prefix))
        plt.close()
        
    else:
        print("No CCs to plot")
        
def filter_other_models(other_models = [], maptype="qfextr"):
    """
    return a list with only those those models which correspond to the correct maptype
    possible maptypes: qfextr fextr kfextr qfgenick fgenick kfgenick qfextr_calc fextr_calc kfextr_calc
    """
    if maptype.lower() in ('qfextr', 'fextr', 'kfextr'):
        to_remove = [fle for fle in other_models if 'genick' in fle.lower() or 'fextr_calc' in fle.lower() or 'for_dm' in fle.lower()]
    elif maptype.lower() in ('qfgenick', 'fgenick', 'kfgenick'):
        to_remove = [fle for fle in other_models if 'fextr' in fle.lower() or 'fextr_calc' in fle.lower() or 'for_dm' in fle.lower()]
    elif maptype.lower() in ('qfextr_calc', 'fextr_calc', 'kfextr_calc'):
        to_remove = [fle for fle in other_models if 'genick' in fle.lower() or 'for_dm' in fle.lower()]
        to_remove_fextr = [fle for fle in other_models if 'fextr' in fle.lower() and not 'fextr_calc' in fle.lower()]
        to_remove += to_remove_fextr
    else:
        print("Invalid maptype, other models cannot be filtered.")
        to_remove = []
    other_models = list(set(other_models).symmetric_difference(set(to_remove)))
    
    return other_models
        

def run(pdb, mtz, other_models = [], other_data = [], phase_info = True, compare_phase_info = True):
    for fle in [pdb, mtz]:
        if os.path.isfile(fle) == False:
            print("File not found: {:s}".format(fle))
            sys.exit(1)
    
    models = []
    #get fmodel from the input model
    model_ref = model(pdb)
    model_ref.extract_f_model(mtz)
    if model_ref.fmodel != None and model_ref.fmodel_update != None:
        print("fmodel with and without update of scales extracted from {:s} and {:s}".format(pdb, mtz))
        models.append(model_ref)
        
    #check other models and get fmodels
    if len(other_models) == 0:
        compare_phase_info = False
        print("no other models provided")
    else:
        for i, pdb in enumerate(other_models):
            print(pdb)
            if os.path.isfile(pdb) == False:
                print("File not found: {:s}".format(pdb))
                continue
            else:
                model_other = model(pdb)
                if len(other_data) == 0:
                    model_other.extract_f_model(mtz)
                elif len(other_models) == len(other_data):
                    data = other_data[i]
                    print(data)
                    if os.path.isfile(data) == False:
                        print("File not found: {:s}".format(data))
                        continue
                    else:
                        model_other.extract_f_model(data)
                else:
                    print("Provide an other_data file for each PDB file, or provide no other_data files at all to use data")
                    compare_phase_info = False
                    break
                if model_other.fmodel != None and model_other.fmodel_update != None:
                    print("fmodel with and without update of scales extracted from {:s} and {:s}".format(pdb, mtz))
                    models.append(model_other)
                    
    if phase_info == False and compare_phase_info == False:
        print("No analysis option given. Use --phase_info or --compare_phase_info")

    if phase_info:
        for m in models:
            print("Plotting phase information from the fmodels without updating scales")
            plot_phases_and_fom(fmodel = m.fmodel,
                                fmodel_prefix = "fmodel_no_scale_update",
                                prefix = "{:s}_phase_info_no_scale_update".format(m.prefix))
            
            #plot updated fmodel phase info
            print("Plotting phase information from the fmodels after updating scales")
            plot_phases_and_fom(fmodel = m.fmodel_update,
                                fmodel_prefix = "fmodel_scale_update",
                                prefix = "{:s}_phase_info_scale_update".format(m.prefix))
    
    if compare_phase_info:
        CCs_fmodel        = []
        CCs_fmodel_update = []
        m_ref = models[0]
        for m in models[1:]:
            print("Comparing phase information from the fmodels without updating scales")
            CCs = compare_phases_and_fom(fmodel_1 = m_ref.fmodel,
                                    fmodel_2 = m.fmodel,
                                    fmodel_1_prefix = "{:s}".format(m_ref.prefix),
                                    fmodel_2_prefix = "{:s}".format(m.prefix),
                                    prefix = "{:s}_{:s}_no_scale_update".format(m_ref.prefix, m.prefix))
            CCs_fmodel.append(CCs)
            
            #plot updated fmodel phase info
            print("Comparing phase information from the fmodels after updating scales")
            CCs = compare_phases_and_fom(fmodel_1 = m_ref.fmodel_update,
                                    fmodel_2 = m.fmodel_update,
                                    fmodel_1_prefix = "{:s}".format(m_ref.prefix),
                                    fmodel_2_prefix = "{:s}".format(m.prefix),
                                    prefix = "{:s}_{:s}_scale_update".format(m_ref.prefix, m.prefix))
            CCs_fmodel_update.append(CCs)
        
        #plot the correlation per model, outnames need to be checked since they will be overwritten in case of an Xtrapol8 with multiple fextr and occupancies
        
        base = "correlations_{:s}_no_scale_update".format(m_ref.prefix)
        prefix_no_scale_update = base
        i = 1
        while os.path.isfile("{:s}.png".format(prefix_no_scale_update)):
            prefix_no_scale_update = "{:s}_{:d}".format(base,i)
            i+=1
            if i == 1000: 
                break
        base = "correlations_{:s}_scale_update".format(m_ref.prefix)
        prefix_scale_update = base
        i = 1
        while os.path.isfile("{:s}.png".format(prefix_scale_update)):
            prefix_scale_update = "{:s}_{:d}".format(prefix_scale_update, i)
            i+=1
            if i == 1000: 
                break
        
        model_names = [m.prefix for m in models[1:]]
        CCs_fmodel = np.array(CCs_fmodel)
        plot_CCs(CC_array = CCs_fmodel,
                 model_names = model_names,
                 prefix = prefix_no_scale_update)
        
        CCs_fmodel_update = np.array(CCs_fmodel_update)
        plot_CCs(CC_array = CCs_fmodel_update,
                 model_names = model_names,
                 prefix = prefix_scale_update)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "extract and fmodel object from a given pdb and mtzfile")
    parser.add_argument('-m', '--model', default='input.pdb', help='Coordinates in pdb format. Just a single model should be provided. Other models can be givven with -o argument.')
    parser.add_argument('-d', '--data', default='input.mtz', help='Data in mtz format.')
    parser.add_argument('-o', '--other_models', default=[], action = 'append', help='Data in mtz format. Optional, if not provided, then the data set given with the -d argument will be used to generate the fmodel. If provided, then exactly the same number of other data sets as other_models should be provided and in the same order.')
    parser.add_argument('-a', '--other_data', default=[], action = 'append', help='Coordinates in pdb format. At least one model should be provided to use the --compare_phase_info argument.')
    parser.add_argument('--phase_info', action = 'store_true', help='Plot phase information from the fmodel')
    parser.add_argument('--compare_phase_info', action = 'store_true', help='Compare the phase info in two models.')
    parser.add_argument('-x8', '--xtrapol8_phil', default= None, help="Xtrapol8 output phil file (Xtrapol8_out.phil) from which the input pdb and mtz files will be extracted. -m, -d, and -o will be ignored.")
    
    #print help if no arguments provided
    if len(sys.argv) < 2:
           parser.print_help()
           sys.exit(1)

    args = parser.parse_args()
    print(args)
    
    pdb = args.model
    mtz = args.data
    other_models = args.other_models
    other_data   = args.other_data
    phase_info   = args.phase_info
    compare_phase_info = args.compare_phase_info
    
    if args.xtrapol8_phil == None:
        run(pdb,
            mtz,
            other_models,
            other_data,
            phase_info = phase_info,
            compare_phase_info = compare_phase_info)
    
    else:
        #reset the input files to avoid that we'll start working with them
        pdb = None
        mtz = None
        other_models = []
        #Extract input from inputfile and command line
        #argument_interpreter = master_phil.command_line_argument_interpreter(home_scope="input")
        input_objects = iotbx.phil.process_command_line_with_files(
            args=[args.xtrapol8_phil],
            master_phil=master_phil
            )
        params = input_objects.work.extract()
        #modified_phil = master_phil.format(python_object=params)
        
        pdb = params.input.reference_pdb
        
        if params.output.generate_phil_only or params.output.generate_fofo_only:
            print("No input files from Xtrapol8 run with output.generate_phil_only or generate_fofo_only")
            sys.exit()
        if params.refinement.run_refinement == False and compare_phase_info == True:
            print("No refinement has been run in Xtrapol8. No models to compare")
            compare_phase_info = False
            
        mtz_and_pdb_files = {}
        maptypes = params.f_and_maps.f_extrapolated_and_maps
        for maptype in maptypes:
            for root, dirs, files in os.walk(params.output.outdir):
                for fle in files:
                    if fle.lower().endswith("{:s}.mtz".format(maptype)):
                        if root.endswith("maps-keep_no_fill"):
                            continue
                        mtz = os.path.join(root, fle)
                        other_models = [os.path.join(root, fle) for fle in os.listdir(root) if fle.lower().endswith(".pdb")]
                        other_models = filter_other_models(other_models, maptype)
                        other_models.sort()
                        mtz_and_pdb_files[mtz] = other_models
                        
        for mtz in mtz_and_pdb_files.keys():
            run(pdb, mtz, other_models = mtz_and_pdb_files[mtz], phase_info = phase_info, compare_phase_info = compare_phase_info)
        

