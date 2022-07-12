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
"""

from __future__ import division, print_function
import re, sys
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cctbx import miller
from Fextr_utils import check_common_indices, make_miller_array
from cctbx.array_family import flex
import pickle

def calculate_k(f_obs_ref, f_obs_2, kweight_scale = 0.05, log=sys.stdout):
    """
    k= 1 / ( 1+ [(sigFlight^2+sigFdark^2)  / mean(sigFlight^2+sigFdark^2)-2) + (kweight_scale * (Flight-Fdark)^2 / mean(Flight-Fdark)^2))
    The mean values should be the mean from the total data set instead of the mean of the bin!
    """
    
    if check_common_indices([f_obs_ref, f_obs_2]) == False:
        print("Input files for k-weighthing don't have identical miller indices. Let's fix it")
        f_obs_ref, f_obs_2 = f_obs_ref.common_sets(f_obs_2)
        if check_common_indices([f_obs_ref, f_obs_2]) == False:
            sys.exit("Input files for k-weighthing don't have identical miller indices and I can't fix it",file=log)
            sys.exit("Input files for k-weighthing don't have identical miller indices and I can't fix it")
            
    sigmadsq_all = f_obs_ref.sigmas()**2+f_obs_2.sigmas()**2
    sigmasmeansq = flex.mean(sigmadsq_all)
    
    datad_all    = f_obs_2.data()-f_obs_ref.data()
    datameandsq  = flex.mean(datad_all**2)

    k        = flex.double()
    indices  = flex.miller_index()
    k_av_lst = []
    k_min_lst        = []
    k_max_lst        = []
    bin_res_cent_lst = []
    f_obs_ref.setup_binner(n_bins=20)
    f_obs_2.use_binning_of(f_obs_ref)
    print("\n************k-weighting statistics************", file=log)
    print("\n************k-weighting statistics************")
    #print("bin  resolution range  #reflections         <k>    Max.k   Min.k   #outliers", file=log)
    print("bin  resolution range  #reflections         <k>    Max.k   Min.k", file=log)
    print("bin  resolution range  #reflections         <k>    Max.k   Min.k")
    for i_bin in f_obs_ref.binner().range_all():
        sel_obs_ref   = f_obs_ref.binner().selection(i_bin)
        f_obs_ref_bin = f_obs_ref.select(sel_obs_ref)
        f_obs_2_bin   = f_obs_2.select(sel_obs_ref)
        if f_obs_ref_bin.size() == 0 : continue
    
        indices_bin = f_obs_ref_bin.indices()
        
        sigmadsq     = f_obs_ref_bin.sigmas()**2+f_obs_2_bin.sigmas()**2
        #sigmasmeansq = flex.mean(sigmadsq**0.5)**2
        sigmaterm    = sigmadsq/sigmasmeansq
        
        datad       = f_obs_2_bin.data()-f_obs_ref_bin.data()
        datadsq     = datad**2
        #datameandsq = flex.mean(datad)**2
        dataterm    = kweight_scale * (datadsq/datameandsq)
        
        k_bin = 1/(1 + sigmaterm + dataterm)
        
        #outlier rejection:
        #diffabs    = flex.abs(datad)
        #sel = k_bin * diffabs > 3* np.sqrt(flex.mean(datadsq))
        #n_outliers = sel.count(True)
        #k_bin.set_selected(sel, 0) 
        
        #Calculate statistics of this bin
        k_bin_av  = np.average(k_bin)
        k_av_lst.append(k_bin_av)
        k_bin_max = np.max(k_bin)
        k_max_lst.append(k_bin_max)
        k_bin_min = np.min(k_bin)
        k_min_lst.append(k_bin_min)
        bin_res_cent = np.average(f_obs_ref.binner().bin_d_range(i_bin))
        bin_res_cent_lst.append(bin_res_cent)
        legend = f_obs_ref.binner().bin_legend(i_bin, show_counts=False)
        #print("{:s} {:^10d} {:> 12.4f} {:> 6.4f} {:> 6.4f} {:^6d}".format(legend, f_obs_ref_bin.size(), k_bin_av, k_bin_max, k_bin_min, n_outliers), file=log)
        print("{:s} {:^10d} {:> 12.4f} {:> 6.4f} {:> 6.4f}".format(legend, f_obs_ref_bin.size(), k_bin_av, k_bin_max, k_bin_min), file=log)
        print("{:s} {:^10d} {:> 12.4f} {:> 6.4f} {:> 6.4f}".format(legend, f_obs_ref_bin.size(), k_bin_av, k_bin_max, k_bin_min))

        k = k.concatenate(k_bin)
        indices = indices.concatenate(indices_bin)
        
    #print('outliers have k set to 0')
    
    assert k.size() == f_obs_ref.sigmas().size(), "k has different size than data. Something went wrong during calculation..."
    
    #caluculate average k
    k_av = flex.mean(k)
    
    #Make Miller array in order to get same reflection order in K-array as in the data
    SG = re.search(r"(.+?)\(No",f_obs_ref.space_group_info().symbol_and_number()).group(1)
    UC = f_obs_ref.unit_cell()
    k_ms = make_miller_array(k, k, SG, UC, indices)
    
    k_ms,_ = k_ms.common_sets(f_obs_ref)
    if check_common_indices([f_obs_ref, k_ms]) == False:
        _,k_ms = f_obs_ref.common_sets(k_ms)
        if check_common_indices([f_obs_ref, k_ms]) == False:
            print('k-weight array has different indices than data sets. k-weighted maps will be nonsense', file=log)
            print('k-weight array has different indices than data sets. k-weighted maps will be nonsense')
    
    # Make graph
    k_av_lst         = np.asarray(k_av_lst)
    k_min_lst        = np.asarray(k_min_lst)
    k_max_lst        = np.asarray(k_max_lst)
    bin_res_cent_lst = np.asarray(bin_res_cent_lst)
    
    out=open('k_estimation.pickle' ,'wb') #write to pickle for GUI
    stats = [bin_res_cent_lst, k_av_lst, k_max_lst, k_min_lst]
    pickle.dump(stats,out)
    out.close()
    
    fig,ax1 = plt.subplots(figsize=(10, 5))
    ax1.set_xlabel('Resolution (A)')
    ax1.set_ylabel('Average k-weight in resolution bin')
    ax1.plot(bin_res_cent_lst[1:], k_av_lst[1:], marker = '.', label='Average k', color = 'red')
    ax1.tick_params(axis='y')
    ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax1.set_ylim(0,1)
    ax2 = ax1.twinx()
    ax2.fill_between(bin_res_cent_lst[1:], k_max_lst[1:], k_min_lst[1:], color='red', alpha=0.2, label='k range')
    ax2.set_ylim(0,1)
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
        
    return k_ms, k_av
