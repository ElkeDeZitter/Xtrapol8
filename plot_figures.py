# -*- coding: utf-8 -*-
"""
This file concerns plotting the Xtrapol8 figures from their respective pickle files.
Function names were given so that they are distinct from the function names in other Xtrapol8 scripts.
A return is called when the file does not exist. This might be replaced by an empty figure or so.
TODO: add the functionalities for the summed_difference_map and distance_analysis

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
from __future__ import division, print_function
import re, sys
import os
import numpy as np
import pickle
from cctbx.array_family import flex
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from Fextr_utils import get_name
import math


def plot_Riso_CCiso(pickle_file='Riso_CCiso.pickle'):
    """
    For Riso_CCiso plot
    Generated only once
    """
    if os.path.isfile(pickle_file) == False:
        return

    with open(pickle_file,'rb') as stats_file:
        bin_res_cent_lst, r_work_lst, cc_work_lst, r_work, cc_work = pickle.load(stats_file)
    
    fig,ax1 = plt.subplots(figsize=(10, 5))
    ax1.plot(bin_res_cent_lst[:], r_work_lst[:], marker = '.', color = 'red', label = 'Riso; overall %.4f' %(r_work))
    ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax1.set_xlabel('Resolution (A)')
    ax1.set_ylabel('Riso')
    ax1.yaxis.label.set_color('red')
    ax2 = ax1.twinx()
    ax2.plot(bin_res_cent_lst[:], cc_work_lst[:], marker = '^', markersize = 5, color = 'green', label = 'CCiso; overall %.4f' %(cc_work))
    ax2.set_ylabel('CCiso')
    ax2.yaxis.label.set_color('green')
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lne, []) for lne in zip(*lines_labels)]
    ax2.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    plt.title('Riso and CCiso for high resolution reflections', fontsize = 'medium',fontweight="bold")
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.80, top = 0.95)
    plt.savefig('Riso_CCiso.pdf', dpi=300, transparent=True)
    plt.savefig('Riso_CCiso.png', dpi=300)
    plt.close()
    
    
def plot_negativereflections(prefix):
    """
    For negative_reflections plot.
    This plot should be placed in the (q/k-weighted)_occupancy directory. If this script is not launched from the (q/k-weighted)_occupancy directory, an additional argument concerning the directory should be added, or pickle file should be directy provided.
    This plot should be prepared for each type of extrapolated structure factor
    prefix can be qFextr, kFextrm, Fextr, qFgenick, kFgenick, Fgenick, qFextr_calc, kFextr_calc, Fextr_calc
    """

    pickle_file = '%s_negative_reflections.pickle'%(prefix) 
    if os.path.isfile(pickle_file) == False:
        return

    
    with open(pickle_file,'rb') as stats_file:
        bin_res_cent_lst,neg_lst,neg_percent_lst,comp_lst, comp_true_lst,s = pickle.load(stats_file)
        
    #In Fextr_utils, "maptype" is used. In order to keep the plotting code the same, let's use the same name
    outname = '%s_negative_reflections' %(prefix)
    
    fig, (ax1, ax3) = plt.subplots(1,2, figsize=(10, 5))
    ax1.plot(bin_res_cent_lst[:], neg_lst[:], marker = ".", label='# Neg. ESFAs',color = 'red')
    ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax1.set_ylim(0, np.max(neg_lst))
    ax1.set_xlabel('Resolution (A)')
    ax1.set_ylabel('Absolute number of negative ESFAs')
    ax1.yaxis.label.set_color('red')
    ax1.set_title("Negative ESFAs for high resolution bins",fontsize = 'medium',fontweight="bold")
    ax2 = ax1.twinx()
    ax2.plot(bin_res_cent_lst[:], neg_percent_lst[:],  marker = 's', markersize = 3, label='% Neg. ESFAs', color = 'blue')
    ax2.set_ylim(0,100)
    ax2.set_ylabel('Negative ESFAs in resolution bin (%)')
    ax2.yaxis.label.set_color('blue')
    lines_labels_1 = [ax.get_legend_handles_labels() for ax in [fig.axes[0],fig.axes[2]]]
    lines_1, labels_1 = [sum(lne, []) for lne in zip(*lines_labels_1)]
    ax1.legend(lines_1, labels_1, loc='lower right', bbox_to_anchor=(0.82, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    
    s = np.full((bin_res_cent_lst.shape[0],1), 90)
    ax3.plot(bin_res_cent_lst[:], comp_lst[:], marker = '.', label='Completeness', color = 'red')
    ax3.plot(bin_res_cent_lst[:], comp_true_lst[:], marker = 's', markersize=3, label='True completeness', color = 'blue')
    ax3.plot(bin_res_cent_lst[:], s[:], linestyle = ':', label= '90 (%) Completeness', color = 'green')
    ax3.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax3.set_ylim(0,100)
    ax3.set_xlabel('Resolution (A)')
    ax3.set_ylabel('Completeness in resolution bin (%)')
    ax3.set_title("Completeness for high resolution bins",fontsize = 'medium',fontweight="bold")
    
    lines_labels_2 = fig.axes[1].get_legend_handles_labels()
    lines_2, labels_2 = [sum(lne, []) for lne in zip(lines_labels_2)]
    #ax3.legend(lines_2, labels_2, fontsize = 'x-small', framealpha=0.5, loc=3, bbox_to_anchor=(0.05, 0.05, 0.5, 0.5))
    ax3.legend(lines_2, labels_2, loc='lower right', bbox_to_anchor=(0.82, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    plt.subplots_adjust(hspace=0.25,wspace=0.5, left=0.09, right=0.88, top = 0.95)
    #plt.title('Completeness and Negative reflections')
    #fig.tight_layout()
    plt.savefig("%s.pdf"%(outname), dpi=300, transparent=True)
    plt.savefig("%s.png"%(outname), dpi=300)
    plt.close()

    
def plot_FsigF(prefix):
    """
    For FsigF plot.
    This plot should be placed in the (q/k-weighted)_occupancy directory
    This plot should be prepared for each type of extrapolated structure factor
    prefix can be qFextr, kFextrm, Fextr, qFgenick, kFgenick, Fgenick, qFextr_calc, kFextr_calc, Fextr_calc
    """
    
    pickle_file = '%s_FsigF.pickle'%(prefix)
    if os.path.isfile(pickle_file) == False:
        return

    
    with open(pickle_file,'rb') as stats_file:
        bin_res_cent_lst, f_sigf_lst, s, l, ids, idl = pickle.load(stats_file)
    
    fig,ax1 = plt.subplots(figsize=(10, 5))
    f_sigf_lst = np.asarray(f_sigf_lst)
    bin_res_cent_lst = np.asarray(bin_res_cent_lst)
    s = np.full((len(f_sigf_lst),1), 0.8)
    l = np.full((len(f_sigf_lst),1), 1.2)
    ax1.plot(bin_res_cent_lst[:], f_sigf_lst[:], marker = '.', label='<F/sig(F)>', color = 'red')
    ax1.plot(bin_res_cent_lst[:], s[:], linestyle = ':', label = '<F/sig(F)> = 0.8', color = 'blue') #(<I/sig(I)> = 2)
    ax1.plot(bin_res_cent_lst[:], l[:], linestyle = ':', label = '<F/sig(F)> = 1.2', color = 'green') #(<I/sig(I)> = 1.5)
    
    ax1.plot(np.array([ids]), np.array([0.8]), marker = 's', markersize=3, color='blue', label = 'estimation: %.2f A' %(ids))
    ax1.plot(np.array([idl]), np.array([1.2]), marker = '^', markersize=5, color = 'green', label = 'estimation: %.2f A' %(idl))
    ax1.set_xlabel('Resolution of bin center (A)')
    ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax1.set_ylabel('<F/sig(F)>')
    ax1.legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    plt.title('%s: <F/sig(F)> for high resolution bins' %(prefix), fontsize = 'medium', fontweight="bold")
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.82, top = 0.95)
    plt.savefig("%s_FsigF.pdf" %(prefix), dpi=300, transparent=True)
    plt.savefig("%s_FsigF.png" %(prefix), dpi=300)
    plt.close()


def plot_refinement_Rfactors_per_alpha(prefix):
    """
    For refinement_R-factors_per_alpha plot.
    This plot should be prepared for each type of extrapolated structure factor
    prefix can be qFextr, kFextrm, Fextr, qFgenick, kFgenick, Fgenick, qFextr_calc, kFextr_calc, Fextr_calc
    """
    #In Fextr_utils, "maptype" is used. In order to keep the plotting code the same, let's use the same name
    maptype = prefix
    
    pickle_file = '%s_refinement_R-factors_per_alpha.pickle' %(prefix)
    if os.path.isfile(pickle_file) == False:
        return

    
    with open(pickle_file,'rb') as stats_file:
        occ_lst, r_work_lst, r_free_lst, r_diff_lst = pickle.load(stats_file)
    
    fig,ax0 = plt.subplots(figsize=(10, 5))
    ax1 = ax0.twinx()
    
    ax0.plot(occ_lst, r_work_lst, color = 'red', marker = 'o', label = 'Rwork')
    ax0.plot(occ_lst, r_free_lst, color = 'blue', marker = 's', markersize = 5, label = 'Rfree')
    ax1.plot(occ_lst, r_diff_lst, color = 'green', marker = '^', label = 'Rfree-Rwork')
    ax0.set_xlabel('Occupancy of triggered state')
    ax0.set_ylabel('R-factor')
    ax1.set_ylabel('R-factor difference')
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lne, []) for lne in zip(*lines_labels)]
    ax1.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.75, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)

    plt.title('R-factors after reciprocal space refinement with %s' %(maptype),fontsize = 'medium',fontweight="bold")
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.82, top = 0.95)
    pltname = '%s_refinement_R-factors_per_alpha' %(maptype)
    plt.savefig("%s.pdf"%(pltname), dpi=300, transparent=True)
    plt.savefig("%s.png"%(pltname), dpi=300)
    plt.close()


def plot_k_estimation(pickle_file="k_estimation.pickle"):
    """
    For k_estimation plot
    Generated only once
    """
    if os.path.isfile(pickle_file) == False:
        return

    
    with open(pickle_file,'rb') as stats_file:
        bin_res_cent_lst, k_av_lst, k_max_lst, k_min_lst = pickle.load(stats_file)
        
    fig,ax1 = plt.subplots(figsize=(10, 5))
    ax1.set_xlabel('Resolution (A)')
    ax1.set_ylabel('Average k-weight in resolution bin')
    ax1.plot(bin_res_cent_lst[:], k_av_lst[:], marker = '.', label='Average k', color = 'red')
    ax1.tick_params(axis='y')
    ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax1.set_ylim(0,np.max(k_max_lst))
    ax2 = ax1.twinx()
    ax2.fill_between(bin_res_cent_lst[:], k_max_lst[:], k_min_lst[:], color='red', alpha=0.2, label='k range')
    ax2.set_ylim(0,np.max(k_max_lst))
    ax2.tick_params(axis='y')
    ax2.set_ylabel('k range within resolution bin')
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lne, []) for lne in zip(*lines_labels)]

    ax2.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.75, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    #fig.tight_layout()
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.82, top = 0.95)
    plt.title('Average k for high resolution reflections',fontsize = 'medium',fontweight="bold")
    plt.savefig('k_estimation.pdf', dpi=300, transparent=True)
    plt.savefig('k_estimation.png', dpi=300)
    plt.close()
    
def plot_q_estimation(pickle_file="q_estimation.pickle"):
    """
    For k_estimation plot
    Generated only once
    """
    if os.path.isfile(pickle_file) == False:
        return
    
    with open(pickle_file,'rb') as stats_file:
        bin_res_cent_lst, q_max_lst, q_min_lst = pickle.load(stats_file)
        
    fig,ax1 = plt.subplots(figsize=(10, 5))
    ax1.set_xlabel('Resolution (A)')
    ax1.set_ylabel('Average q in resolution bin')
    ax1.plot(bin_res_cent_lst[:], q_av_lst[:], marker = '.', label='Average q', color = 'red')
    ax1.tick_params(axis='y')
    ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax1.set_ylim(0,1)
    ax2 = ax1.twinx()
    ax2.fill_between(bin_res_cent_lst[:], q_max_lst[:], q_min_lst[:], color='red', alpha=0.2, label='q range')
    ax2.set_ylim(0,1)
    ax2.tick_params(axis='y')
    ax2.set_ylabel('q range within resolution bin')
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lne, []) for lne in zip(*lines_labels)]

    ax2.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.75, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    plt.title('Average q for high resolution reflections', fontsize = 'medium',fontweight="bold")
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.82, top = 0.95)
    plt.savefig('q_estimation.pdf', dpi=300, transparent=True)
    plt.savefig('q_estimation.png', dpi=300)
    plt.close()

def plot_Neg_Pos_reflections(suffix, pickle_file='Fextr_negative.pickle'):
    """
    For Neg_Pos_reflections plot.
    This plot should be prepared for each type of extrapolated structure factor
    This plot should be updated during run
    The input pickle file contains information on all types of extrapolation (the name is thus misleading)!!!!
    In this script, we reed the pickle and search specifically for maptype (qFextr, kFextrm, Fextr, qFgenick, kFgenick, Fgenick, qFextr_calc, kFextr_calc, Fextr_calc) and make the plot. In that sense it differs in the way the plots are generated in the command-line version
    """
    #In Fextr_utils, "maptype" is used. In order to keep the plotting code the same, let's use the same name
    maptype = suffix
    if os.path.isfile(pickle_file) == False:
        return
    
    with open(pickle_file,'rb') as stats_file:
         while True:
             try:
                 stats = np.array(pickle.load(stats_file))
                 #stats = np.array(tuple(stats), dtype='f8, S32, i4, i4, i4')
                 try:
                     alldata = np.vstack([alldata,stats[np.newaxis,...]])
                 except NameError:
                     alldata =stats[np.newaxis,...]
             except EOFError:
                 break

    indices = np.where(alldata[:,1] == maptype)[0]
    plt.close()
    fig, (ax0, ax2) = plt.subplots(1,2, figsize=(10, 5))
    ax1 = ax0.twinx()
    ax3 = ax2.twinx()
    
    for a in indices:
        ax0.plot(np.float(alldata[a][0]),  np.float(alldata[a][-1]), marker='o', color='red', label = maptype)
        ax2.plot(np.float(alldata[a][0]),  np.float(alldata[a][3]), marker='o', color='red', label = maptype)
            
    ax0.set_ylim(0,np.float(alldata[a][2]))
    ax1.set_ylim(0,100)
    ax2.set_ylim(0,np.float(alldata[a][2]))
    ax3.set_ylim(0,100)
    
    ax0.set_xlabel('Occupancy')
    ax0.set_ylabel('Absolute number')
    ax1.set_ylabel('Percentage of total number of reflections')
    ax0.set_title('%s: Negative reflections' %(maptype), fontsize = 'medium',fontweight="bold")
    
    ax2.set_xlabel('Occupancy')
    ax2.set_ylabel('Absolute number')
    ax3.set_ylabel('Percentage of total number of reflections')
    ax2.set_title('%s: Positive reflections'%(maptype), fontsize = 'medium',fontweight="bold")
    
    plt.subplots_adjust(hspace=0.25,wspace=0.5, left=0.09, right=0.88, top = 0.95)
    plt.savefig("Neg_Pos_reflections_%s.pdf" %(maptype), dpi=300, transparent=True)
    plt.savefig("Neg_Pos_reflections_%s.png" %(maptype), dpi=300)
    plt.close()

def plot_alpha_occupancy_determination(suffix):
    """
    For alpha_occupancy_determination plot.
    This plot should be prepared for each type of extrapolated structure factor
    prefix can be qFextr, kFextrm, Fextr, qFgenick, kFgenick, Fgenick, qFextr_calc, kFextr_calc, Fextr_calc
    """
    pickle_file = 'alpha_occupancy_determination_%s' %(suffix)
    if os.path.isfile(pickle_file) == False:
        return
    
    with open(pickle_file,'rb') as stats_file:
        alphas, occupancies, int1_norm, int2_norm, resids_lst_used, alphafound = pickle.load(stats_file)
        
    if resids_lst_used == False:
        ax1.plot(alphas, int1, 's', markersize = 5, color = 'blue', label='All peaks') #int1, int2 and conv will be the same, so we can plot only int1 and avoid the try catch
    else:
        try:
            ax1.plot(alphas, int1_norm, 'o', color = 'red', label='Selected residues')      
            ax1.plot(alphas, int2_norm, 's', markersize = 5, color = 'blue', label='All peaks')
            ax1.plot(alphas, results, '^', color="green", label='Selected residues with enhanced SNR')
        except: 
            ax1.plot(alphas, int1_norm, 'o', color = 'red')

    ax1.set_ylim([0.,1.1])
    ax1.set_xlim([np.min(alphas)*0.95,np.max(alphas)*1.05])
    ax1.set_xlabel('Alpha value = 1/occupancy')
    ax1.set_ylabel('Normalized difference map ratio')

    if resids_lst_used == False:
        ax2.plot(occupancies, int1_norm, 's', markersize = 5, color = 'blue', label='All peaks') #int1, int2 and conv will be the same, so we can plot only int1 and avoid the try catch
    else:
        try:
            ax2.plot(occupancies, int1_norm, 'o', color = 'red', label='Selected residues') 
            ax2.plot(occupancies, int2_norm, 's', markersize = 5, color = 'blue', label='All peaks')     
            ax2.plot(occupancies, results, '^', color="green", label='Selected residues with enhanced SNR')
        except:
            ax2.plot(alphas,int1, 'o', color = 'red')
            
    ax2.set_ylim([0.,1.1])
    ax2.set_xlim([np.min(occupancies)*0.95,np.max(occupancies)*1.05])
    ax2.set_xlabel('Triggered state occupancy')
    ax2.set_ylabel('Normalized difference map ratio')
    ax2.legend(loc='lower right', bbox_to_anchor=(0.92, -0.05, 0.45, 0.5), fontsize = 'x-small', framealpha=0.5)

    if alphafound:
        ax1.set_title('Alpha determination', fontsize = 'medium',fontweight="bold")
        ax2.set_title('Occupancy determination', fontsize = 'medium',fontweight="bold")
    else:
        ax1.set_title('Alpha determination IMPOSSIBLE', fontsize = 'medium',fontweight="bold")
        ax1.text(np.min(self.alphas), 0.5, 'no peaks found in at least one of the maps\n for the selected residues')
        ax2.set_title('Occupancy determination IMPOSSIBLE', fontsize = 'medium',fontweight="bold")
        ax2.text(np.min(self.occupancies), 0.5, 'no peaks found in at least one of the maps\n for the selected residues')
    plt.subplots_adjust(hspace=0.25,wspace=0.4, left=0.09, right=0.88, top = 0.95)
    plt.savefig("%s.pdf"%(outname), dpi=300, transparent=True,bbox_inches='tight', pad_inches = 0)
    plt.savefig("%s.png"%(outname), dpi=300,bbox_inches='tight', pad_inches = 0)
    plt.close()

def plot_FoFosigmas(pickle_file='Fextr_binstats.pickle'):
    """
    For FoFo-sigmas plot.
    Generated only once
    In this script, we reed the pickle and search only the info concerning the FoFo. In that sense it differs in the way the plots are generated in the command-line version
    """
    if os.path.isfile(pickle_file) == False:
        return
        
    with open(pickle_file,'rb') as stats_file:
        _, FoFo_type, _, bin_res_cent_lst, _, _, fdif_data_lst, fdif_sigmas_lst = np.array(pickle.load(stats_file))

    mn = 0
    mx = 0
    if min(fdif_data_lst[1:]) < mn:
        mn = min(fdif_data_lst[1:])
    if max(fdif_data_lst[1:]) > mx:
        mx = max(fdif_data_lst[1:])
    if min(fdif_sigmas_lst[1:]) < mn:
        mn = min(fdif_sigmas_lst[1:])
    if max(fdif_sigmas_lst[1:]) > mx:
        mx = max(fdif_sigmas_lst[1:])

    fig, ax0 = plt.subplots(1,1, figsize=(10, 5))
    ax1 = ax0.twinx()
    
    ax0.plot(bin_res_cent_lst[:], fdif_data_lst[:], marker = '.', color = 'tab:red', label = '%s' %(FoFo_type))
    ax1.plot(bin_res_cent_lst[:], fdif_sigmas_lst[:], marker = 'x', linestyle='--', color = 'tab:blue', label = 'sig(%s)' %(FoFo_type))
    ax0.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax0.set_xlabel('Resolution (A)')
    ax0.set_ylabel('%s' %(FoFo_type))
    ax0.yaxis.label.set_color('tab:red')
    ax1.set_ylabel('sig(%s)' %(FoFo_type))
    ax1.yaxis.label.set_color('tab:blue')
    lines_labels_1 = [ax.get_legend_handles_labels() for ax in [ax0,ax1]]
    lines_1, labels_1 = [sum(lne, []) for lne in zip(*lines_labels_1)]
    ax0.legend(lines_1, labels_1, loc='lower right', bbox_to_anchor=(0.73, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    ax0.set_title('%s for high resolution reflections' %(FoFo_type),fontsize = 'medium',fontweight="bold")
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.85, top = 0.95)
    plt.savefig('%s_sigmas.pdf' %(FoFo_type),dpi=300, transparent=True)
    plt.savefig('%s_sigmas.png' %(FoFo_type),dpi=300)
    plt.close()

def plot_Fextrsigmas(prefix, pickle_file='Fextr_binstats.pickle'):
    """
    For Fextr_sigmas plot.
    This plot should be prepared for each type of extrapolated structure factor
    This plot should be updated during run
    The input pickle file contains information on all types of extrapolation (the name is thus misleading)!!!!
    In this script, we reed the pickle and search specifically for maptype (qFextr, kFextrm, Fextr, qFgenick, kFgenick, Fgenick, qFextr_calc, kFextr_calc, Fextr_calc) and make the plot. In that sense it differs in the way the plots are generated in the command-line version
    """
    #In Fextr_utils, "maptype" is used. In order to keep the plotting code the same, let's use the same name
    maptype = prefix
    #Read the pickle file
    with open(pickle_file,'rb') as stats_file:
        while True:
            try:
                stats = np.array(pickle.load(stats_file))
                #stats = np.array(tuple(stats), dtype='f8, S32, i4, i4, i4')
                try:
                    alldata = np.vstack([alldata,stats[np.newaxis,...]])
                except NameError:
                    alldata =stats[np.newaxis,...]
            except EOFError:
                break
        
    #extract the occupancies
    occ_lst = list(set(alldata[:,0]))
    #get the indices concerning the specific maptype we are looking at
    indices = np.where(alldata[:,2] == maptype)[0]
    
    mn = 0
    mx = 0
    fig, ax0 = plt.subplots(1,1, figsize=(10, 5))
    ax1 = ax0.twinx()
    for a in indices:
        occ, _, _, bin_res_cent_lst, fextr_data_lst, fextr_sigmas_lst, _, _ = alldata[a]
        color_data  = plt.cm.Reds(int((np.log(1/occ))*100))
        color_sigma = plt.cm.Blues(int((np.log(1/occ))*100))
        ax0.plot(bin_res_cent_lst[:], fextr_data_lst[:], marker = '.', color = color_data, label = '%s, occ = %.3f' %(maptype, occ))
        ax1.plot(bin_res_cent_lst[:], fextr_sigmas_lst[:], marker = 'x', linestyle='--', color = color_sigma, label = 'sig(%s), occ = %.3f' %(maptype, occ))
        #Specify the minimum and maximum value
        if min(fextr_data_lst[1:]) < mn:
            mn = min(fextr_data_lst[1:])
        if max(fextr_data_lst[1:]) > mx:
            mx = max(fextr_data_lst[1:])
        if min(fextr_sigmas_lst[1:]) < mn:
            mn = min(fextr_sigmas_lst[1:])
        if max(fextr_sigmas_lst[1:]) > mx:
            mx = max(fextr_sigmas_lst[1:])
    
    ax0.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax0.set_xlabel('Resolution (A)')#, fontsize = 'small')
    ax0.set_ylabel('ESFAs')#, fontsize = 'small')
    ax0.yaxis.label.set_color('tab:red')
        
    ax0.set_ylim(mn, mx)
    ax1.set_ylim(mn, mx)

    ax0.legend(loc='lower right', bbox_to_anchor=(0.89, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    ax1.legend(loc='lower right', bbox_to_anchor=(1.17, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    ax1.set_ylabel('sig(ESFAs)')#, fontsize = 'small')
    ax1.yaxis.label.set_color('tab:blue')

    ax0.set_title('%s for high resolution reflections' %(maptype), fontsize = 'medium',fontweight="bold")
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.65, top = 0.95)
    plt.savefig('%s_sigmas.pdf' %(maptype),dpi=300, transparent=True)
    plt.savefig('%s_sigmas.png' %(maptype),dpi=300)
    plt.close()

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)


def plot_ddm(pickle_file = 'ddm.pickle', scale = 1.5):
    """
    For plotting the ddm.
    This plot should be prepared for each of the Fextr types but only with the model from the best estimated occupancy and should calculated only at the end
    The input pickle contains information on the chain id, ddm_residue and seq_info_unique
    The output filename should be the same as of the pickle file
    """
    with open(pickle_file,'rb') as stats_file:
        while True:
            try:
                stats = np.array(pickle.load(stats_file))
                try:
                    alldata = np.vstack([alldata,stats[np.newaxis,...]])
                except NameError:
                    alldata =stats[np.newaxis,...]
            except EOFError:
                break
    
    chains = np.unique(alldata[:,0])
    num_chains = chains.shape[0]
        
    if num_chains == 1:
        n_cols = 1
    else:
        n_cols = 2
        
    n_rows_prot = int(math.ceil(num_chains/n_cols))

    n_rows = n_rows_prot

    outname = get_name(pickle_file)
    outname_pdf = '%s.pdf'%(outname)
    outname_png = '%s.png'%(outname)
    
    if scale == None:
        scale = 0
        for ddm_residue in alldata[:,1]:
            mx_temp = np.max(ddm_residue)
            mn_temp = np.min(ddm_residue)
            if mx_temp > scale:
                scale = mx_temp
            if np.abs(mn_temp) > scale:
                scale = np.abs(mn_temp)
    else:
        scale = scale
    
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(10*n_rows, 10*n_cols), squeeze=False)#, constrained_layout=True)
    #fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)
    c = mcolors.ColorConverter().to_rgb
    rvb = make_colormap([c('blue'), c('white'), 0.40, c('white'),  0.60, c('white'), c('red')])

    col = -1
    row = -1
    old_chain_id = ''
    for chain in alldata:
        ID, ddm_residue,seq_info_unique = chain
        if ID != old_chain_id:
            if col == 0:
                col = 1
            else:
                col = 0
                row +=1
                        
        mask = np.zeros_like(ddm_residue, dtype=np.bool)
        mask[np.triu_indices_from(mask)] = True
        FINAL2 = np.ma.array(ddm_residue, mask=mask)
        
        tick_jump = int(np.round(len(seq_info_unique)/15,0)) #we want to add 15 seq_ticks
        seq_ticks = seq_info_unique[0::tick_jump] #residues for which we will show tick positions
        tick_pos = list(map(lambda x: x*tick_jump, range(len(seq_ticks))))
        
            
        img = axs[row, col].imshow(FINAL2, cmap=rvb, vmin=-scale, vmax=scale)
        
        axs[row, col].set_xticks(tick_pos)
        axs[row, col].set_xticklabels(seq_ticks)
        plt.setp(axs[row, col].xaxis.get_majorticklabels(), rotation=90, fontsize='small')
        
        axs[row, col].set_yticks(tick_pos)
        axs[row, col].set_yticklabels(seq_ticks)

        axs[row, col].set_title('Chain %s' %ID, fontsize = 'medium',fontweight="bold")
        axs[row, col].set_xlabel('Residues')
        axs[row, col].set_ylabel('Residues')
        
        axs[row, col].spines['top'].set_visible(False)
        axs[row, col].spines['right'].set_visible(False)
        
        fig.colorbar(img, ax=axs[row, col], fraction=0.046, pad=0.04)
            
        old_chain_id = ID

        
    fig.tight_layout()
    plt.savefig(outname_pdf, dpi=300, transparent=True)
    plt.savefig(outname_png, dpi=300)
    plt.close()
    
