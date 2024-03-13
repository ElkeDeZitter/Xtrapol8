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

def outlier_rejection_only(f_obs_ref, f_obs_2, log=sys.stdout):
    """
    Outlier rejection only. Removes differences larger than 3 standard deviations.
    TODO: add correction for centric and acentric reflections
    """
    
    if check_common_indices([f_obs_ref, f_obs_2]) == False:
        print("Input files for outlier rejection don't have identical miller indices. Let's try fixing it")
        f_obs_ref, f_obs_2 = f_obs_ref.common_sets(f_obs_2)
        if check_common_indices([f_obs_ref, f_obs_2]) == False:
            sys.exit("Input files for outlier rejection don't have identical miller indices and I can't fix it",file=log)
            sys.exit("Input files for outlier rejection don't have identical miller indices and I can't fix it")
            
    c = flex.double() #to avoid confusion called 'c' instead of 'q'
    indices  = flex.miller_index()
    
    print("\n************outlier rejection only statistics************", file=log)
    print("\n************outlier rejection only statistics************")
    f_obs_ref.setup_binner(n_bins=20)
    f_obs_2.use_binning_of(f_obs_ref)
    for i_bin in f_obs_ref.binner().range_all():
        sel_obs_ref   = f_obs_ref.binner().selection(i_bin)
        f_obs_ref_bin = f_obs_ref.select(sel_obs_ref)
        f_obs_2_bin   = f_obs_2.select(sel_obs_ref)
        if f_obs_ref_bin.size() == 0 : continue

        c_bin = flex.double(f_obs_ref_bin.data().size(), 1)
        indices_bin = f_obs_ref_bin.indices()
        
        #outlier rejection:
        diffabs    = flex.abs(f_obs_2_bin.data()-f_obs_ref_bin.data())
        diffsquare = (f_obs_2_bin.data()-f_obs_ref_bin.data())**2
        sel = diffabs > 3* np.sqrt(flex.mean(diffsquare))
        legend = f_obs_ref.binner().bin_legend(i_bin, show_counts=False)
        print ('{:s} {:^10d} {:^6d} outliers for which correction factor c set to 0'.format(legend, f_obs_ref_bin.size(), sel.count(True)), file=log)
        print ('{:s} {:^10d} {:^6d} outliers for which correction factor c set to 0'.format(legend, f_obs_ref_bin.size(), sel.count(True)))
        c_bin.set_selected(sel, 0.0) 
    
        c = c.concatenate(c_bin)
        indices = indices.concatenate(indices_bin)
        
    #Make Miller array in order to get same reflection order in Q-array as in the data
    SG = re.search(r"(.+?)\(No",f_obs_ref.space_group_info().symbol_and_number()).group(1)
    UC = f_obs_ref.unit_cell()
    c_ms = make_miller_array(c, c, SG, UC, indices)
    
    c_ms,_ = c_ms.common_sets(f_obs_ref)
    if check_common_indices([f_obs_ref, c_ms]) == False:
        _,c_ms = f_obs_ref.common_sets(c_ms)
        if check_common_indices([f_obs_ref, c_ms]) == False:
            print('Outlier rejection array has different indices than data sets. Wrong reflections will be weighted down to zero', file=log)
            print('Outlier rejection array has different indices than data sets. Wrong reflections will be weighted down to zero')

    return c_ms

def calculate_q(f_obs_ref, f_obs_2, log=sys.stdout):
    
    if check_common_indices([f_obs_ref, f_obs_2]) == False:
        #print("Input files for q-weighthing don't have identical miller indices. Let's try fixing it", file=log)
        print("Input files for q-weighthing don't have identical miller indices. Let's try fixing it")
        f_obs_ref, f_obs_2 = f_obs_ref.common_sets(f_obs_2)
        if check_common_indices([f_obs_ref, f_obs_2]) == False:
            sys.exit("Input files for q-weighthing don't have identical miller indices and I can't fix it", file=log)
            sys.exit("Input files for q-weighthing don't have identical miller indices and I can't fix it")

    cent = f_obs_ref.centric_flags()
    mult = f_obs_ref.multiplicities()
    mult = flex.double(map(lambda x: float(x), mult.data()))
    mult = miller.array(miller_set= cent, data = mult)
    stol = f_obs_ref.unit_cell().stol(f_obs_ref.indices())
    res  = 1/(2*stol)
    res  = miller.array(miller_set= cent, data = res)
    num  = flex.double(map(lambda x: float(x+1),cent.data()))
    num  = miller.array(miller_set= cent, data = num)

    DECONV   = flex.double()
    NODECONV = flex.double()
    sigmadsq = flex.double()
    q        = flex.double()
    indices  = flex.miller_index()
    q_av_lst = []
    q_min_lst        = []
    q_max_lst        = []
    sigmadsq_lst     = []
    bin_res_cent_lst = []
    f_obs_ref.setup_binner(n_bins=20)
    f_obs_2.use_binning_of(f_obs_ref)
    cent.use_binning_of(f_obs_ref)
    mult.use_binning_of(f_obs_ref)
    res.use_binning_of(f_obs_ref)
    num.use_binning_of(f_obs_ref)
    
    print("\n************q-weighting statistics************", file=log)
    print("\n************q-weighting statistics************")
    print("bin  resolution range  #reflections <deltaFdeconv.est.> <deltaFnon-deconv.est.> <q>  Max.q  Min.q  #outliers  sDeconvolution", file=log)
    print("bin  resolution range  #reflections <deltaFdeconv.est.> <deltaFnon-deconv.est.> <q>  Max.q  Min.q  #outliers  sDeconvolution")
    for i_bin in f_obs_ref.binner().range_all():
        sel_obs_ref   = f_obs_ref.binner().selection(i_bin)
        f_obs_ref_bin = f_obs_ref.select(sel_obs_ref)
        f_obs_2_bin   = f_obs_2.select(sel_obs_ref)
        cent_bin      = cent.select(sel_obs_ref)
        mult_bin      = mult.select(sel_obs_ref)
        res_bin       = res.select(sel_obs_ref)
        num_bin       = num.select(sel_obs_ref)
        if f_obs_ref_bin.size() == 0 : continue
    
        indices_bin = f_obs_ref_bin.indices()
        
        diffsquare = (f_obs_2_bin.data()-f_obs_ref_bin.data())**2
        prod       = num_bin.data()*diffsquare
        qout       = prod/mult_bin.data()
        nodeconv   = flex.mean(qout)
        NODECONV.append(nodeconv)
        
        sigmasum   = flex.mean(f_obs_ref_bin.sigmas()**2+f_obs_2_bin.sigmas()**2)
        diffsquare = diffsquare-sigmasum
        prod       = num_bin.data()*diffsquare
        qout       = prod/mult_bin.data()
        deconv     = flex.mean(qout)
        DECONV.append(deconv)

        if deconv >= (nodeconv/2.0):
            sigmadsq_bin = deconv
            message ='deconvoluted estimate of difference is big enough...'
        else:
            sigmadsq_bin = nodeconv/2.0
            message ='deconvoluted estimate of differences too small...I have to correct...'
        sigmadsq_lst.append(sigmadsq_bin)
        enumerator_bin = mult_bin.data()*sigmadsq_bin/num_bin.data()
        divisor_bin    = (f_obs_ref_bin.sigmas()**2) + (f_obs_2_bin.sigmas()**2) + enumerator_bin
        q_bin          = enumerator_bin/divisor_bin
        
        #outlier rejection:
        diffabs    = flex.abs(f_obs_2_bin.data()-f_obs_ref_bin.data())
        diffsquare = (f_obs_2_bin.data()-f_obs_ref_bin.data())**2
        sel = q_bin * diffabs > 3* np.sqrt(flex.mean(diffsquare))
        n_outliers = sel.count(True)
        q_bin.set_selected(sel, 0)
        
        #Calculate statistics of this bin
        q_bin_av  = np.average(q_bin)
        q_av_lst.append(q_bin_av)
        q_bin_max = np.max(q_bin)
        q_max_lst.append(q_bin_max)
        q_bin_min = np.min(q_bin)
        q_min_lst.append(q_bin_min)
        bin_res_cent = np.average(f_obs_ref.binner().bin_d_range(i_bin))
        bin_res_cent_lst.append(bin_res_cent)
        legend = f_obs_ref.binner().bin_legend(i_bin, show_counts=False)
        print("{:s} {:^10d} {:> 10.4f} {:> 20.4f} {:> 12.4f} {:> 6.4f} {:> 6.4f} {:^6d} {:s}".format(legend, f_obs_ref_bin.size(), deconv, nodeconv, q_bin_av, q_bin_max, q_bin_min, n_outliers, message),file=log)
        print("{:s} {:^10d} {:> 10.4f} {:> 20.4f} {:> 12.4f} {:> 6.4f} {:> 6.4f} {:^6d} {:s}".format(legend, f_obs_ref_bin.size(), deconv, nodeconv, q_bin_av, q_bin_max, q_bin_min, n_outliers, message))

        q = q.concatenate(q_bin)
        indices = indices.concatenate(indices_bin)
        
    print('outliers have q set to 0', file=log)
    print('outliers have q set to 0')
    
    assert q.size() == f_obs_ref.sigmas().size(), "q has different size than data. Something went wrong during calculation..."
    #caluculate average q
    q_av = flex.mean(q)
    
    #Make Miller array in order to get same reflection order in Q-array as in the data
    SG = re.search(r"(.+?)\(No",f_obs_ref.space_group_info().symbol_and_number()).group(1)
    UC = f_obs_ref.unit_cell()
    q_ms = make_miller_array(q, q, SG, UC, indices)
    
    q_ms,_ = q_ms.common_sets(f_obs_ref)
    if check_common_indices([f_obs_ref, q_ms]) == False:
        _,q_ms = f_obs_ref.common_sets(q_ms)
        if check_common_indices([f_obs_ref, q_ms]) == False:
            print('Q-weight array has different indices than data sets. Reflections will get wrong weight...', file=log)
            print('Q-weight array has different indices than data sets. Reflections will get wrong weight...')
    
    # Make graph
    q_av_lst         = np.asarray(q_av_lst)
    q_min_lst        = np.asarray(q_min_lst)
    q_max_lst        = np.asarray(q_max_lst)
    bin_res_cent_lst = np.asarray(bin_res_cent_lst)
    sigmadsq_lst     = np.asarray(sigmadsq_lst)
    
    out=open('q_estimation.pickle' ,'wb') #write to pickle for GUI
    stats = [bin_res_cent_lst, q_av_lst, q_max_lst, q_min_lst]
    pickle.dump(stats,out)
    out.close()
    
    fig,ax1 = plt.subplots(figsize=(10, 5))
    ax1.set_xlabel('Resolution (A)')
    ax1.set_ylabel('q-weight')
    ax1.plot(bin_res_cent_lst[1:], q_av_lst[1:], marker = '.', label='Average q', color = 'red')
    #ax1.semilogy(bin_res_cent_lst[1:], q_av_lst[1:], marker = '.', label='Average q', color = 'red')
    #ax1.fill_between(bin_res_cent_lst[1:], q_max_lst[1:], q_min_lst[1:], color='red', alpha=0.2)
    ax1.tick_params(axis='y')
    ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax1.set_ylim(0,1)
    ax2 = ax1.twinx()
    #ax2.set_ylabel('Average difference F estimate in resolution bin')
    #ax2.semilogy(bin_res_cent_lst[1:], sigmadsq_lst[1:], marker = '.', label='Average deltaF estimate in resolution bin', color = 'blue')
    #ax2.tick_params(axis='y')
    ax2.fill_between(bin_res_cent_lst[1:], q_max_lst[1:], q_min_lst[1:], color='red', alpha=0.2, label='q range')
    ax2.set_ylim(0,1)
    ax2.tick_params(axis='y')
    ax2.set_ylabel('q-weight')
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lne, []) for lne in zip(*lines_labels)]

    #fig.legend(lines, labels, fontsize = 'x-small', framealpha=0.5, loc=6, bbox_to_anchor=(0.1, 0.1, 0.5, 0.5))
    #fig.tight_layout()
    ax2.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.75, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    plt.title('Average q for high resolution reflections', fontsize = 'medium',fontweight="bold")
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.82, top = 0.95)
    plt.savefig('q_estimation.pdf', dpi=300, transparent=True)
    plt.savefig('q_estimation.png', dpi=300)
    plt.close()
        
    return q_ms, q_av
