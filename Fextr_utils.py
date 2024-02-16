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
import os
import random
import numpy as np
import iotbx
import math
import pickle
from cctbx import miller, crystal, xray
from iotbx import mtz, symmetry
from iotbx.pdb import hierarchy
from iotbx.file_reader import any_file
from cctbx.array_family import flex
#import mmtbx.model
#import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy.stats import linregress
import uuid

def get_unique_id(id_length=20):
    """
    Function to get a unique id based on UUID with length id_length
    """
    if id_length > 36:
        id_length == 36
    return str(uuid.uuid4())[:id_length]

def generate_log_name(time_stamp):
    """
    Generate a unique name for the Xtrapol8 logfile.
    A short uuid of 20 characters is added to the logfile name.
    """
    uuid = get_unique_id(36)
    logname = "%s_Xtrapol8_%s.log" %(time_stamp, uuid)
    
    return logname

def remove_unique_id_from_log(log_name):
    """
    Remove the unqiue sequence from the log file
    """
    index = log_name.find("Xtrapol8")+len("Xtrapol8")
    new_name = log_name[:index]+".log"
    os.rename(log_name, new_name)
    
    return new_name

def update_file_paths_coot(outdir_old, outdir_new, coot_prefix = "coot_all_"):
    """
    Change the path in output files after changing the output directory.
    Needs to be run after renaming outdir_old to outdir_new
    """
    coot_scripts = [os.path.join(root, fle) for root, dirs, files in os.walk(outdir_new) for fle in files if fle.startswith(coot_prefix)]
    for coot_script in coot_scripts:
        with open(coot_script, "r") as f:
            fle = f.read().split("\n")
        o = open(coot_script,"w")
        for lne in fle:
            if outdir_old in lne:
                lne = re.sub(outdir_old, outdir_new, lne)
                o.write("%s\n"%lne)
            else:
                o.write("%s\n"%lne)
        o.close()
        
def update_file_paths_pymol(outdir_old, outdir_new, pymol_file = "pymol_movie.py"):
    """
    Change the path in output files after changing the output directory.
    Needs to be run after renaming outdir_old to outdir_new
    """
    pymol_script = '%s/%s' %(outdir_new, pymol_file)
    if os.path.isfile(pymol_script):
        with open(pymol_script, "r") as f:
            fle = f.read().split("\n")
        o = open(pymol_script,"w")
        for lne in fle:
            if outdir_old in lne:
                lne = re.sub(outdir_old, outdir_new, lne)
                o.write("%s\n"%lne)
            else:
                o.write("%s\n"%lne)
        o.close()
        
def update_file_paths_occupancy_recap(outdir_old, outdir_new, pickle_file = "occupancy_recap.pickle"):
    """
    Change the coot and ddm path in the occupancy_recap.pickle file
    """
    with open(pickle_file, "rb") as pf:
        results = pickle.load(pf)
         
    #change the paths of the coot and ddm files and safe in a new dictionary
    occ_overview = {}
    for fextr in results.keys():
        occ, coot, ddm = results[fextr]
        coot_corr = re.sub(outdir_old, outdir_new, coot)
        if ddm != None:
            ddm_corr = re.sub(outdir_old, outdir_new, ddm)
        else:
            ddm_corr = ddm
        occ_overview[fextr] = [occ,coot_corr,ddm_corr]

    #overwrite the pickle file with the new dictionary
    occ_pickle = open(pickle_file, "wb")
    pickle.dump(occ_overview, occ_pickle)
    occ_pickle.close()
    
def update_file_paths_log(log_name, outdir_old, outdir_new):
    """
    Change the paths in the log file after changing the output directory.
    Needs to be run after renaming outdir_old to outdir_new and after closing the log file.
    This might need to be done on a different moment than the update of other files.
    """
    full_log = "%s/%s" %(outdir_new, log_name)
    with open(full_log, "r") as f:
        fle = f.read().split("\n")
    o = open(full_log,"w")
    for lne in fle:
        if outdir_old in lne:
            lne = re.sub(outdir_old, outdir_new, lne)
            o.write("%s\n"%lne)
        else:
            o.write("%s\n"%lne)
    o.close()

def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def check_program_path(program):
    program_exist = False
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            program_exist = True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                program_exist = True
                break
                
    return exe_file, program_exist

def get_phenix_version():
    """
    Get the phenix version based on the full path of the phenix executable
    """
    phenix_path, program_exists = check_program_path("phenix")
    if program_exists:
        phenix_version = re.search(r"phenix-(.+?)\/", phenix_path).group(1)
    else:
        print("phenix could not be found")
        phenix_version = 0
        
    return phenix_version
    
def get_ccp4_version():
    """
    Get the ccp4 version based on the full path of the scaleit executable
    """
    scaleit_path, program_exists = check_program_path("scaleit")
    if program_exists:
        ccp4_version = re.search(r"ccp4-(.+?)\/", scaleit_path).group(1)
    else:
        print("scaleit (ccp4) could not be found")
        ccp4_version = 0
        
    return ccp4_version
    
def list_redundant_files(outdir):
    o = open('redundant_files.txt','w')
    for dirpath, dirnames, filenames in os.walk(os.getcwd()):
        for fle in filenames:
            if (fle.endswith("map") or fle.endswith("sh")):
                o.write("%s : %s\n" %(os.path.join(dirpath, fle), file_size(os.path.join(dirpath, fle))))
    o.close()

def convert_bytes(num):
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, unit)
        num /= 1024.0

def file_size(file_path):
    if os.path.isfile(file_path):
        file_info = os.stat(file_path)
        return convert_bytes(file_info.st_size)


def phenix_version_from_logfile(logfile):
    with open(logfile) as l:
        fle = l.read().split('\n')
    V_lne = [lne for lne in fle if 'Version' in lne][0]
    version = int(re.search(r"Version: 1\.(.+)\.", V_lne).group(1))
    return version

def get_name (fle):
    '''
    Get the name of the file without the path and format

    Parameters:
    -----------
    file (str)
        path of the file from which the name is searched
    file_format (format, example: .cif)
        format of the file input

    Returns:
    -----------
    name (str)
        name of the file without path or format
    '''
    try: #iotbx.file_reader.splitext from phenix 1.19
        _, file_format, _=iotbx.file_reader.splitext(fle)
    except ValueError: #iotbx.file_reader.splitext from phenix 1.18
        _, file_format =os.path.splitext(fle)
    #get format of the input file
    if "/" in fle:
        #if the path of the file contains '/'
        # name=re.split("/", file)[-1]
        # #the name of the file takes the name of the file in the last folder
        # if file[-4]+file[-3]+file[-2]+file[-1]== file_format:
        # #if the 4 last strings of the name are the file format
        #     name=name[:-4]
        #     #the file format is deducted from the name
        #find name independent of the lengt of the file_format
        name = re.search(r"\/(.+?)\%s$"%(file_format), fle).group(1).split("/")[-1]
    else:
        name = re.sub("\%s$"%(file_format), "", fle)

    #if len(name)>80:
        
        #new_name = raw_input("%s is a long name, please a shorter prefix to use. If you provide nothing (just press enter) I will use '%s' as prefix for the output files." %(name, name[:30]+name[-30:]))
        #if new_name == '' or new_name == ' ':
            #name = name[:30]+name[-30:]
        #else:
            #name = new_name
        #print("%s is a long name, let's call it '%s' in the output files." %(name, name[:30]+name[-30:]))
        #name = name[:30]+name[-30:]
    return name

def get_pdb_name(pdb_file):
    if "/" in pdb_file:
        name = re.search(r"\/(.+?)\.pdb$", pdb_file).group(1).split("/")[-1]
    else:
        name = re.sub("\.pdb$","",pdb_file)
    #if len(name)>80:
        #print("%s is a long name, let's call it %s" %(name, name[:30]+name[-30:]))
        #name = name[:30]+name[-30:]
    return name


def make_miller_array(data, sigma, SG, UC, indices):
    ms = miller.set(
        crystal_symmetry=crystal.symmetry(
            space_group_symbol=SG,
            unit_cell=UC),
        anomalous_flag=False,
        indices=indices)
    return ms.array(data=data, sigmas=sigma)

def phase_transfer(miller_array, phase_source): #taken from fobs_minus_fobs_map. Phases cannot be added to miller_array but have to be inserted with the data
    tmp = miller.array(miller_set = miller_array, data = flex.double(miller_array.indices().size(), 1)).phase_transfer(phase_source = phase_source)
    return miller.array(miller_set = miller_array, data = miller_array.data() * tmp.data() )
    #phased_data = miller_array.data() * tmp.data()
    #phased_ms = make_miller_array(phased_data, sigmas, SG, UC, indices)
    #return phased_ms
    
def check_file_existance(fle):
    if os.path.isfile(fle) == True:
        fle = fle
    else:
        print("%s cannot be found. Check if it is existing." %(fle))
        fle = ''
    return fle

def append_if_file_exist(lst, fle):
    if os.path.isfile(fle):
        lst.append(fle)
    return

def check_and_delete_hydrogen(pdb_file):
    pdb_hier = hierarchy.input(file_name=pdb_file)
    if pdb_hier.hierarchy.remove_hd() == 0:
        pdb_ed = pdb_file
    else:
        pdb_hier.hierarchy.remove_hd()
        p = open('model_edit.pdb', 'w')
        p.write(pdb_hier.hierarchy.as_pdb_string(crystal_symmetry=pdb_hier.input.crystal_symmetry()),output_break_records=False)
        p.close()
        pdb_ed = os.path.abspath(check_file_existance('model_edit.pdb'))
    return pdb_ed

def check_and_delete_altlocs(pdb_file):
    pdb_hier = hierarchy.input(file_name=pdb_file)
    if pdb_hier.hierarchy.altloc_indices().size() == 1:
        pdb_ed = pdb_file
    else:
        pdb_hier.hierarchy.remove_alt_confs(always_keep_one_conformer=True)
        p = open('model_NoAltlocs.pdb', 'w')
        p.write(pdb_hier.hierarchy.as_pdb_string(crystal_symmetry=pdb_hier.input.crystal_symmetry()),output_break_records=False)
        p.close()
        pdb_ed = os.path.abspath(check_file_existance('model_NoAltlocs.pdb'))
    return pdb_ed

def open_all_in_coot(FoFo, pdb_list, mtz_list, additional, outdir=os.getcwd(), suffix=''):
    """
    Not elegant at all, to be rewritten
    maps are not in logical order
    """
    
    script_coot = '%s/coot_all_%s.py'%(outdir, suffix)
    i = open(script_coot,'w')
    i.write('set_nomenclature_errors_on_read("auto-correct")\n') #choice between "auto-correct", "ignore", "prompt"
    
    additional_lines = ""
    if len(additional)>0:
        for cif in additional.split():
            if cif.endswith(".cif"):
                additional_lines+='read_cif_dictionary("%s")\n'%(cif)

    pdb_lines = ""
    for pdb in pdb_list:
        pdb_lines+='handle_read_draw_molecule("%s")\n' %(pdb)
    
    default_mtz_lines = ""
    special_mtz_lines = ""
    if len(mtz_list)>0:
        for mtz in mtz_list:
            if os.path.isfile(mtz):
                mtz_content = any_file(mtz,force_type="hkl").file_object.file_content()
                if '2FOFCWT' in mtz_content.column_labels(): #maps from refinement program
                    default_mtz_lines+='auto_read_make_and_draw_maps_from_mtz("%s")\n'%(mtz)
                elif 'FDM' in mtz_content.column_labels(): #maps from dm
                    default_mtz_lines+='auto_read_make_and_draw_maps_from_mtz("%s")\n'%(mtz)
                elif mtz_content.column_labels() == ['H', 'K', 'L', 'FP', 'SIGFP', 'PHIB', 'FOM', 'HLA', 'HLB', 'HLC', 'HLD', 'FWT', 'PHWT']: #from phenix.density_modification # not used anymore
                    default_mtz_lines+='auto_read_make_and_draw_maps_from_mtz("%s")\n'%(mtz)
                else:
                    if mtz_content.column_types() == ['H', 'H', 'H', 'F', 'P', 'F', 'P']:
                        special_mtz_lines += 'set_auto_read_column_labels("%s","%s",0)\nset_auto_read_column_labels("%s","%s",1)\nauto_read_make_and_draw_maps_from_mtz("%s")\n'%(mtz_content.column_labels()[3],mtz_content.column_labels()[4],mtz_content.column_labels()[5],mtz_content.column_labels()[6], mtz)

    FoFo_labels = any_file(FoFo,force_type="hkl").file_object.file_content().column_labels()
    FoFo_lines = 'make_and_draw_map("%s","%s","%s","",0,1)' %(FoFo,FoFo_labels[3],FoFo_labels[4])
    
    i.write('%s%s%s%s%s'%(pdb_lines, additional_lines, default_mtz_lines, special_mtz_lines, FoFo_lines))
    i.close()
    
    return script_coot
    #os.system("coot --script %s" %(script_coot))

def get_common_indices_dict(refl_dict): #Use at own risk, some columns might be lost
    #argument to be given should be like {'f_obs_off':f_obs_off, 'f_obs_on':f_obs_on, 'phase_info': phase_info}
    refl_dict_comm = {}
    for a in refl_dict: 
        for b in refl_dict:
            if b != a:
                a_comm = a+"_comm"
                b_comm = b+"_comm"
                a_comm, b_comm = refl_dict[a].common_sets(refl_dict[b])
                refl_dict_comm[a+"_comm"] = a_comm
                refl_dict_comm[b+"_comm"] = b_comm
    outlst = []
    for name in refl_dict_comm:
        outlst.append(refl_dict_comm[name])
    
    return outlst

def check_common_indices(refl_lst):
    for a in refl_lst:
        for b in refl_lst:
            if b != a:
                if a.indices().all_eq(b.indices()) == False or a.data().size() != b.data().size():
                    compatible = False
                    break
                else:
                    compatible = True
    return compatible

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def neg_neflecions_binning(miller_array, prefix, log=sys.stdout):
    """
    print and plot the number of negative reflections per bin.
    """
    outname = '%s_negative_reflections' %(prefix)
    bin_res_cent_lst = []
    neg_lst          = []
    neg_percent_lst  = []
    comp_lst         = []
    comp_true_lst    = []
    miller_array_noneg = miller_array.select(miller_array.data()>=0)
    miller_array.setup_binner(n_bins=20)
    miller_array_noneg.use_binning_of(miller_array)
    comp_bin_data = miller_array.completeness(use_binning=True).data
    #avoid problems when value is None, this is probably a consequence of a bad resolution cutoff
    comp_bin_data = [0.0 if v is None else v for v in comp_bin_data] 
    comp_true_bin_data = miller_array_noneg.completeness(use_binning=True).data
    #avoid problems when value is None, this is probably a consequence of a bad resolution cutoff
    comp_true_bin_data = [0.0 if v is None else v for v in comp_true_bin_data]
    print("\n************Completeness of extrapolated structure factors************", file=log)
    print("\n************Completeness of extrapolated structure factors************")
    print("bin  resolution range     #reflections  #neg.reflections  completeness(%)  compl.pos.reflections(%)",file=log)
    print("bin  resolution range     #reflections  #neg.reflections  completeness(%)  compl.pos.reflections(%)")
    for i_bin in miller_array.binner().range_all():
        sel = miller_array.binner().selection(i_bin)
        f_bin = miller_array.select(sel)
        if f_bin.size() == 0 : continue
        bin_res_cent = np.median(miller_array.binner().bin_d_range(i_bin))
        bin_res_cent_lst.append(bin_res_cent)
        neg = f_bin.select(~(f_bin.data()>=0)).data().size()
        neg_lst.append(neg)
        neg_percent = neg/f_bin.data().size() * 100
        neg_percent_lst.append(neg_percent)
        comp_lst.append(comp_bin_data[i_bin]*100)
        comp_true_lst.append(comp_true_bin_data[i_bin]*100)
        legend = miller_array.binner().bin_legend(i_bin, show_counts=False)
        print('{:s} {:^15d} {:>5d} ({:>5.2f}%) {:>15.2f} {:>15.2f}'.format(legend, f_bin.size(), neg, neg_percent, comp_bin_data[i_bin]*100, comp_true_bin_data[i_bin]*100), file=log)
        print('{:s} {:^15d} {:>5d} ({:>5.2f}%) {:>15.2f} {:>15.2f}'.format(legend, f_bin.size(), neg, neg_percent, comp_bin_data[i_bin]*100, comp_true_bin_data[i_bin]*100))
        #print("%s  %8d  %8d (%.2f %%)    %.2f%%      %.2f%%" % (legend, f_bin.size(), neg, neg_percent, comp_bin_data[i_bin]*100, comp_true_bin_data[i_bin]*100), file=log)

    neg_lst          = np.asarray(neg_lst)
    neg_percent_lst  = np.asarray(neg_percent_lst)
    bin_res_cent_lst = np.asarray(bin_res_cent_lst)
    comp_lst         = np.asarray(comp_lst)
    comp_true_lst    = np.asarray(comp_true_lst)
    assert neg_lst.shape==neg_percent_lst.shape==bin_res_cent_lst.shape==comp_lst.shape==comp_true_lst.shape
                
    plt.close()
    fig, (ax1, ax3) = plt.subplots(1,2, figsize=(10, 5))
    #ax1.plot(bin_res_cent_lst[1:], neg_lst[1:], linestyle = '-', label='# Neg. reflections',color = 'red')
    ax1.plot(bin_res_cent_lst[1:], neg_lst[1:], marker = ".", label='# Neg. ESFAs',color = 'red')
    ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax1.set_ylim(0, np.max(neg_lst))
    ax1.set_xlabel('Resolution (A)')
    ax1.set_ylabel('Absolute number of negative ESFAs')
    ax1.yaxis.label.set_color('red')
    ax1.set_title("Negative ESFAs for high resolution reflections",fontsize = 'medium',fontweight="bold")
    ax2 = ax1.twinx()
    #ax2.plot(bin_res_cent_lst[1:], neg_percent_lst[1:], linestyle = '-', label='% Neg. reflections', color = 'blue')
    ax2.plot(bin_res_cent_lst[1:], neg_percent_lst[1:], marker = 's', markersize = 3, label='% Neg. ESFAs', color = 'blue')
    ax2.set_ylim(0,100)
    ax2.set_ylabel('Percentage of negative ESFAs (%)')
    ax2.yaxis.label.set_color('blue')
    lines_labels_1 = [ax.get_legend_handles_labels() for ax in [fig.axes[0],fig.axes[2]]]
    lines_1, labels_1 = [sum(lne, []) for lne in zip(*lines_labels_1)]
    ax1.legend(lines_1, labels_1, loc='lower right', bbox_to_anchor=(0.82, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    
    s = np.full((bin_res_cent_lst.shape[0],1), 90)
    ax3.plot(bin_res_cent_lst[1:], comp_lst[1:], marker = '.', label='Completeness', color = 'red')
    ax3.plot(bin_res_cent_lst[1:], comp_true_lst[1:], marker = 's', markersize=3, label='True completeness', color = 'blue')
    ax3.plot(bin_res_cent_lst[1:], s[1:], linestyle = ':', label= '90 (%) Completeness', color = 'green')
    ax3.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax3.set_ylim(0,100)
    ax3.set_xlabel('Resolution (A)')
    ax3.set_ylabel('Completeness (%)')
    ax3.set_title("Completeness for high resolution reflections",fontsize = 'medium',fontweight="bold")
    
    #x_point1 = np.where(np.abs(comp_lst - 90) == np.partition(np.abs(comp_lst - 90), 0)[0])[0][0]
    #if x_point1 ==0: #lowest resolution bin, try second smallest value
        #x_point1 = np.where(np.abs(comp_lst - 90) == np.partition(np.abs(comp_lst - 90), 1)[1])[0][0]
    #if 0 < x_point1 < comp_lst.shape[0]-1:
        #if comp_lst[x_point1-1] > comp_lst[x_point1] and comp_lst[x_point1+1] > comp_lst[x_point1]:
            ##local minimum, try second smallest value
            #x_point1 = np.where(np.abs(comp_lst - 90) == np.partition(np.abs(comp_lst - 90), 1)[1])[0][0]
    #x1 = bin_res_cent_lst[x_point1]
    #y1 = comp_lst[x_point1]
    #try:
        #x_point2 = x_point1+1
        #x2 = bin_res_cent_lst[x_point2]
        #y2 = comp_lst[x_point2]
        #a = (y2-y1)/(x2-x1)
        #b = y1 - a*x1
        #idc = (4-b)/a
    #except IndexError:
        #idc = x1
    
    #x_point1 = np.where(np.abs(comp_true_lst - 90) == np.partition(np.abs(comp_true_lst - 90), 0)[0])[0][0]
    #if x_point1 ==0: #lowest resolution bin, try second smallest value
        #x_point1 = np.where(np.abs(comp_true_lst - 90) == np.partition(np.abs(comp_true_lst - 90), 1)[1])[0][0]
    #if 0 < x_point1 < comp_true_lst.shape[0]-1:
        #if comp_true_lst[x_point1-1] > comp_true_lst[x_point1] and comp_true_lst[x_point1+1] > comp_true_lst[x_point1]:
            ##local minimum, try second smallest value
            #x_point1 = np.where(np.abs(comp_true_lst - 90) == np.partition(np.abs(comp_true_lst - 90), 1)[1])[0][0]
    #x1 = bin_res_cent_lst[x_point1]
    #y1 = comp_true_lst[x_point1]
    #try:
        #x_point2 = x_point1+1
        #x2 = bin_res_cent_lst[x_point2]
        #y2 = comp_true_lst[x_point2]
        #a = (y2-y1)/(x2-x1)
        #b = y1 - a*x1
        #idt = (4-b)/a
    #except IndexError:
        #idt = x1
    
    #ax3.plot(np.array([idc]), np.array([90]), 'ro', label = '90 (%) Completeteness at {:.2f} A'.format(idc))
    #ax3.plot(np.array([idt]), np.array([90]), 'go', label = '90 (%) True completeteness at {:.2f} A'.format(idt))
    
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
    
    out=open('%s.pickle'%(outname) ,'wb') #write to pickle for GUI
    stats = [bin_res_cent_lst,neg_lst,neg_percent_lst,comp_lst, comp_true_lst,s]
    pickle.dump(stats,out)
    out.close()

    
    #fig, host = plt.subplots()
    #fig.subplots_adjust(right=0.55)
    #par1 = host.twinx()
    #par2 = host.twinx()
    #par3 = host.twinx()
    #par2.spines["right"].set_position(("axes", 1.25))
    #par3.spines["right"].set_position(("axes",1.5))
    #make_patch_spines_invisible(par2)
    #make_patch_spines_invisible(par3)
    #par2.spines["right"].set_visible(True)
    #par3.spines["right"].set_visible(True)
    #p1, = host.plot(bin_res_cent_lst[:], neg_lst[:], '-', label='# Neg. reflections', color = 'red')
    #p2, = par1.plot(bin_res_cent_lst[:], comp_lst[:], '-', label='Completeness', color = 'green')
    #p3, = par2.plot(bin_res_cent_lst[:],comp_true_lst[:], '-', label='True completeness', color = 'yellow')
    #p4, = par3.plot(bin_res_cent_lst[:], neg_percent_lst[:], '-', label='Neg. reflections', color = 'blue')
    #host.set_xlim(np.max(bin_res_cent_lst[:]), np.min(bin_res_cent_lst[:]))
    #host.set_ylim(0, np.max(neg_lst))
    #par1.set_ylim(0,100)
    #par2.set_ylim(0,100)
    #par3.set_ylim(0,100)
    #host.set_xlabel('Resolution (A)')
    #host.set_ylabel('Number of negative reflections in resolution bin')
    #par1.set_ylabel('Completeness in resolution bin (%)')
    #par2.set_ylabel('True completeness (%) = compl. of pos. refl. only')
    #par3.set_ylabel('Negative reflections in resolution bin (%)')
    #tkw = dict(size=4, width=1.5)
    #host.tick_params(axis='y', colors=p1.get_color(), **tkw)
    #par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    #par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    #par3.tick_params(axis='y', colors=p4.get_color(), **tkw)
    #host.tick_params(axis='x', **tkw)
    #lines = [p1, p2, p3, p4]
    #host.legend(lines, [l.get_label() for l in lines], framealpha=0.5, loc=3, bbox_to_anchor=(0.05, 0.05, 0.5, 0.5))
    #plt.title('Completeness and Negative reflections')
    #fig.tight_layout()
    #plt.savefig(outname, dpi=300, transparent=True)
    #plt.close()

def make_rfree_col(miller_array, fraction):
    tot = miller_array.data().size()
    num_free = int(tot*fraction)
    free_ind = random.sample(range(tot), num_free)
    free_col = np.zeros(tot, dtype=np.int32)
    for i in range(tot):
        if i in free_ind:
            free_col[i]=1
    return miller.array(miller_set = miller_array, data = flex.int(free_col))

    
def make_fwork_ffree(miller_array, r_free_flags):
    r_free_flags = r_free_flags.data().as_bool() #if r_free_flags are generated with make_rfree_col
    f_free = miller_array.select(r_free_flags)
    f_work = miller_array.select(~(r_free_flags))
    return f_work, f_free

def compute_f_sigf(miller_array, prefix, log=sys.stdout):
    """
    Calculate <F/sigma(F)> as equivalent of <I/sigma(I)> per resolution bin to estimate resolution of the (extrapolated) data. It also gives and estimation of the useful resolution, based on <F/sig(F)>
    """
    f_sigf_lst       = []
    bin_res_cent_lst = []
    miller_array.setup_binner(n_bins=20)
    print("\n************Extrapolated structure factors signal strength************", file=log)
    print("bin  resolution range  #reflections <F/sigma(F)>", file=log)
    print("\n************Extrapolated structure factors signal strength************")
    print("bin  resolution range  #reflections <F/sigma(F)>")
    for i_bin in miller_array.binner().range_all():
        sel = miller_array.binner().selection(i_bin)
        f_bin = miller_array.select(sel)
        if f_bin.size() == 0 : continue
        bin_res_cent = np.median(miller_array.binner().bin_d_range(i_bin))
        bin_res_cent_lst.append(bin_res_cent)
        f_sigf = np.average(f_bin.data()/f_bin.sigmas())
        f_sigf_lst.append(f_sigf)
        legend = miller_array.binner().bin_legend(i_bin, show_counts=False)
        print("%s  %8d  %.4f" % (legend, f_bin.size(), f_sigf), file=log)
        print("%s  %8d  %.4f" % (legend, f_bin.size(), f_sigf))
    
    plt.close()
    fig,ax1 = plt.subplots(figsize=(10, 5))
    f_sigf_lst = np.asarray(f_sigf_lst)
    bin_res_cent_lst = np.asarray(bin_res_cent_lst)
    s = np.full((len(f_sigf_lst),1), 0.8)
    l = np.full((len(f_sigf_lst),1), 1.2)
    ax1.plot(bin_res_cent_lst[1:], f_sigf_lst[1:], marker = '.', label='<F/sig(F)>', color = 'red')
    ax1.plot(bin_res_cent_lst[1:], s[1:], linestyle = ':', label = '<F/sig(F)> = 0.8', color = 'blue') #(<I/sig(I)> = 2)
    ax1.plot(bin_res_cent_lst[1:], l[1:], linestyle = ':', label = '<F/sig(F)> = 1.2', color = 'green') #(<I/sig(I)> = 1.5)
    
    x_point1 = np.where(np.abs(f_sigf_lst - 0.8) == np.partition(np.abs(f_sigf_lst - 0.8), 0)[0])[0][0]
    if x_point1 ==0: #lowest resolution bin, try second smallest value
        x_point1 = np.where(np.abs(f_sigf_lst - 0.8) == np.partition(np.abs(f_sigf_lst - 0.8), 1)[1])[0][0]
    if 0 < x_point1 < f_sigf_lst.shape[0]-1:
        if f_sigf_lst[x_point1-1] > f_sigf_lst[x_point1] and f_sigf_lst[x_point1+1] > f_sigf_lst[x_point1]:
            #local minimum, try second smallest value
            x_point1 = np.where(np.abs(f_sigf_lst - 0.8) == np.partition(np.abs(f_sigf_lst - 0.8), 1)[1])[0][0]
    x1 = bin_res_cent_lst[x_point1]
    y1 = f_sigf_lst[x_point1]
    try:
        x_point2 = x_point1+1
        x2 = bin_res_cent_lst[x_point2]
        y2 = f_sigf_lst[x_point2]
        a = (y2-y1)/(x2-x1)
        b = y1 - a*x1
        ids = (0.8-b)/a
    except IndexError:
        ids = x1
    
    x_point1 = np.where(np.abs(f_sigf_lst - 1.2) == np.partition(np.abs(f_sigf_lst - 1.2), 0)[0])[0][0]
    if x_point1 ==0: #lowest resolution bin, try second smallest value
        x_point1 = np.where(np.abs(f_sigf_lst - 1.2) == np.partition(np.abs(f_sigf_lst - 1.2), 1)[1])[0][0]
    if 0 < x_point1 < f_sigf_lst.shape[0]-1:
        if f_sigf_lst[x_point1-1] > f_sigf_lst[x_point1] and f_sigf_lst[x_point1+1] > f_sigf_lst[x_point1]:
            #local minimum, try second smallest value
            x_point1 = np.where(np.abs(f_sigf_lst - 1.2) == np.partition(np.abs(f_sigf_lst - 1.2), 1)[1])[0][0]
    x1 = bin_res_cent_lst[x_point1]
    y1 = f_sigf_lst[x_point1]
    try:
        x_point2 = x_point1+1
        x2 = bin_res_cent_lst[x_point2]
        y2 = f_sigf_lst[x_point2]
        a = (y2-y1)/(x2-x1)
        b = y1 - a*x1
        idl = (1.2-b)/a
    except IndexError:
        idl = x1
    
    ax1.plot(np.array([ids]), np.array([0.8]), marker = 's', markersize=3, color='blue', label = 'estimation: %.2f A' %(ids))
    ax1.plot(np.array([idl]), np.array([1.2]), marker = '^', markersize=5, color = 'green', label = 'estimation: %.2f A' %(idl))
    ax1.set_xlabel('Resolution (A)')
    ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax1.set_ylabel('<F/sig(F)>')
    ax1.legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    #fig.tight_layout()
    plt.title('%s: <F/sig(F)> for high resolution reflections' %(prefix), fontsize = 'medium', fontweight="bold")
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.82, top = 0.95)
    plt.savefig("%s_FsigF.pdf" %(prefix), dpi=300, transparent=True)
    plt.savefig("%s_FsigF.png" %(prefix), dpi=300)
    plt.close()
    
    out=open("%s_FsigF.pickle" %(prefix) ,'wb') #write to pickle for GUI
    stats = [bin_res_cent_lst, f_sigf_lst, s, l, ids, idl]
    pickle.dump(stats,out)
    out.close()

    
def plot_Rfactors_per_alpha(refine_log_lst, maptype):
    """
    Extract R-values from log-file phenix or refmac (Refmac with minimal output) and plot in function of the occupancy.
    """
    
    plt.close()
    fig,ax0 = plt.subplots(figsize=(10, 5))
    ax1 = ax0.twinx()
    
    r_work_lst = flex.double()
    r_free_lst = flex.double()
    occ_lst    = flex.double()
    for log_file in refine_log_lst:
        if os.path.isfile(log_file):
            occ = float(re.search('occupancy\_(.+?)\/', log_file).group(1))
            occ_lst.append(occ)
            with open(log_file) as fle:
                log = fle.readlines()
            if 'refmac' in log_file: #refinement performed in refmac
                r_work = float([line for line in log if  'R factor' in line][-1].split()[-1])
                r_free = float([line for line in log if  'R free' in line][-1].split()[-1])
            else: #refinement performed in phenix
                #r_line = log[-1]
                r_line = [lne for lne in log if lne.lstrip().startswith("Final R-work")][0]
                r_work,r_free = re.findall('Final R-work = (.+?), R-free = (.+?)\\n',r_line)[0]
            r_work_lst.append(float(r_work))
            r_free_lst.append(float(r_free))
    srt = flex.sort_permutation(occ_lst)
    occ_lst = occ_lst.select(srt)
    r_work_lst = r_work_lst.select(srt)
    r_free_lst = r_free_lst.select(srt)
    r_diff_lst = r_free_lst - r_work_lst
    ax0.plot(occ_lst, r_work_lst, color = 'red', marker = 'o', label = 'Rwork')
    ax0.plot(occ_lst, r_free_lst, color = 'blue', marker = 's', markersize = 5, label = 'Rfree')
    ax1.plot(occ_lst, r_diff_lst, color = 'green', marker = '^', label = 'Rfree-Rwork')
    ax0.set_xlabel('Triggered state occupancy')
    ax0.set_ylabel('R-factor')
    #ax0.legend(fontsize = 'xx-small', framealpha=0.5, loc='lower left', bbox_to_anchor=(0.0, 0., 0.5, 0.5))
    #ax1.set_xlabel('Occupancy of triggered state')
    ax1.set_ylabel('R-factor difference')
    ax1.yaxis.label.set_color('green')
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lne, []) for lne in zip(*lines_labels)]
    ax1.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.75, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)

    plt.title('R-factors after reciprocal space refinement with %s' %(maptype),fontsize = 'medium',fontweight="bold")
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.82, top = 0.95)
    pltname = '%s_refinement_R-factors_per_alpha' %(maptype)
    #fig.tight_layout()
    plt.savefig("%s.pdf"%(pltname), dpi=300, transparent=True)
    plt.savefig("%s.png"%(pltname), dpi=300)
    plt.close()
    
    out=open('%s_refinement_R-factors_per_alpha.pickle' %(maptype) ,'wb') #write to pickle for GUI
    stats = [occ_lst, r_work_lst, r_free_lst, r_diff_lst]
    pickle.dump(stats,out)
    out.close()

def get_Fextr_stats(occ, fextr_ms, Fextr_type, fdif_ms, FoFo_type, outdir, log=sys.stdout):
    
    bin_res_cent_lst       = flex.double()
    fextr_data_lst         = flex.double()
    fextr_sigmas_lst       = flex.double()
    fdif_data_lst          = flex.double()
    fdif_sigmas_lst        = flex.double()
    print("\n************Difference and extrapolated structure factors*************")
    if fextr_ms.data().size() == fdif_ms.data().size():
        fextr_ms.setup_binner(n_bins=20)
        fdif_ms.use_binning_of(fextr_ms)
        #print("Occupancy: %.3f, alpha: %2.3f" %(occ, 1/occ))
        print("bin  resolution range  #reflections <{0}> <sig({0})> <{1}> <sig({1})>".format(FoFo_type, Fextr_type))
        for i_bin in fextr_ms.binner().range_all():
            sel_fextr_ms     = fextr_ms.binner().selection(i_bin)
            fextr_ms_bin     = fextr_ms.select(sel_fextr_ms)
            fdif_ms_bin      = fdif_ms.select(sel_fextr_ms)
            if fextr_ms_bin.size() == 0 : continue
            legend = fextr_ms.binner().bin_legend(i_bin, show_counts=False)
            bin_res_cent_avg = np.median(fextr_ms.binner().bin_d_range(i_bin))
            bin_res_cent_lst.append(bin_res_cent_avg)
            fextr_data_avg = flex.mean(fextr_ms_bin.data())
            fextr_data_lst.append(fextr_data_avg)
            fextr_sigmas_avg = flex.mean(fextr_ms_bin.sigmas())
            fextr_sigmas_lst.append(fextr_sigmas_avg)
            fdif_data_avg = flex.mean(fdif_ms_bin.data())
            fdif_data_lst.append(fdif_data_avg)
            fdif_sigmas_avg = flex.mean(fdif_ms_bin.sigmas())
            fdif_sigmas_lst.append(fdif_sigmas_avg)
            print('{:s} {:^10d} {:> 10.4f} {:> 10.4f} {:> 10.4f} {:> 10.4f}'.format(legend, fextr_ms_bin.size(), fdif_data_avg, fdif_sigmas_avg, fextr_data_avg, fextr_sigmas_avg))

        assert bin_res_cent_lst.size()==fextr_data_lst.size()==fextr_sigmas_lst.size()==fdif_data_lst.size()==fdif_sigmas_lst.size(), 'list sizes to plot F and sigmas not equal'
    
    else: #If reflections are rejected, then fextr_ms.size() != fdif_ms.size()
        fdif_ms.setup_binner(n_bins=20)
        fextr_ms.setup_binner(n_bins=20)
        print("bin  resolution range_{0}  #reflections_{0}  <{0}> <sig({0})> | bin  resolution range_{1}  #reflections_{1} <{1}> <sig({1})>".format(FoFo_type, Fextr_type))
        for i_bin in fdif_ms.binner().range_all():
            sel_fdif_ms     = fdif_ms.binner().selection(i_bin)
            fdif_ms_bin     = fdif_ms.select(sel_fdif_ms)
            sel_fextr_ms    = fextr_ms.binner().selection(i_bin)
            fextr_ms_bin    = fextr_ms.select(sel_fextr_ms)            
            if fdif_ms_bin.size() == 0 : continue
            legend_fdif  = fdif_ms.binner().bin_legend(i_bin, show_counts=False)
            bin_res_cent_avg = np.median(fdif_ms.binner().bin_d_range(i_bin))
            bin_res_cent_lst.append(bin_res_cent_avg)
            fdif_data_avg = flex.mean(fdif_ms_bin.data())
            fdif_data_lst.append(fdif_data_avg)
            fdif_sigmas_avg = flex.mean(fdif_ms_bin.sigmas())
            fdif_sigmas_lst.append(fdif_sigmas_avg)
            legend_fextr = fextr_ms.binner().bin_legend(i_bin, show_counts=False)
            #bin_res_cent_avg = np.median(fextr_ms.binner().bin_d_range(i_bin))
            #bin_res_cent_lst.append(bin_res_cent_avg)
            fextr_data_avg = flex.mean(fextr_ms_bin.data())
            fextr_data_lst.append(fextr_data_avg)
            fextr_sigmas_avg = flex.mean(fextr_ms_bin.sigmas())
            fextr_sigmas_lst.append(fextr_sigmas_avg)
            print('{:s} {:^20d} {:> 10.4f} {:> 10.4f} | {:s} {:^20d} {:> 10.4f} {:> 10.4f}'.format(legend_fdif, fdif_ms_bin.size(), fdif_data_avg, fdif_sigmas_avg, legend_fextr, fextr_ms_bin.size(), fextr_data_avg, fextr_sigmas_avg))
        
        
    #fextr_ms.setup_binner(n_bins=20)
    #print("bin  resolution range  #reflections <%s> <sig(%s)>" %(Fextr_type, Fextr_type))
    #for i_bin in fextr_ms.binner().range_all():
        #sel_fextr_ms     = fextr_ms.binner().selection(i_bin)
        #fextr_ms_bin     = fextr_ms.select(sel_fextr_ms)
        #if fextr_ms_bin.size() == 0 : continue
        #legend = fextr_ms.binner().bin_legend(i_bin, show_counts=False)
        #bin_res_cent_avg = np.median(fextr_ms.binner().bin_d_range(i_bin))
        #bin_res_cent_lst.append(bin_res_cent_avg)
        #fextr_data_avg = flex.mean(fextr_ms_bin.data())
        #fextr_data_lst.append(fextr_data_avg)
        #fextr_sigmas_avg = flex.mean(fextr_ms_bin.sigmas())
        #fextr_sigmas_lst.append(fextr_sigmas_avg)
        #print('%s  %8d  %6.4f   %6.4f'%(legend, fextr_ms_bin.size(), fextr_data_avg, fextr_sigmas_avg))
        
    
    stats = [occ, FoFo_type, Fextr_type, bin_res_cent_lst, fextr_data_lst, fextr_sigmas_lst, fdif_data_lst, fdif_sigmas_lst]
    
    out=open('%s/Fextr_binstats.pickle' %(outdir),'ab') #write to pickle to avoid keeping in memory and for GUI
    pickle.dump(stats,out)
    out.close()
    print("Stats saved to %s/Fextr_binstats.pickle" %(outdir))
    
def plot_Fextr_sigmas(pickle_file='Fextr_binstats.pickle'):
    """
    To plot the Values of the difference map and extrapolated structure factors in a single figure for each extrapolated structure factor
    """
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
        
    #extract the maptypes
    maptype_lst = list(set(alldata[:,2]))
    #extract the occupancies
    occ_lst = list(set(alldata[:,0]))
    
    #For each maptype make the Fextr / sig(Fextr) plot for each maptype
    for maptype in maptype_lst:
        #print(maptype)
        mn = 0
        mx = 0
        plt.close()
        fig, ax0 = plt.subplots(1,1, figsize=(10, 5))
        ax1 = ax0.twinx()
        #get the indices concerning the specific maptype we are looking at
        indices = np.where(alldata[:,2] == maptype)[0]
        for a in indices:
            occ, _, _, bin_res_cent_lst, fextr_data_lst, fextr_sigmas_lst, _, _ = alldata[a]
            color_data  = plt.cm.Reds(int((np.log(1/occ))*100))
            color_sigma = plt.cm.Blues(int((np.log(1/occ))*100))
            ax0.plot(bin_res_cent_lst[1:], fextr_data_lst[1:], marker = '.', color = color_data, label = '%s, occ = %.3f' %(maptype, occ))
            ax1.plot(bin_res_cent_lst[1:], fextr_sigmas_lst[1:], marker = 'x', linestyle='--', color = color_sigma, label = 'sig(%s), occ = %.3f' %(maptype, occ))
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
        #ax0.tick_params(labelsize='x-small')
          
        ax0.set_ylim(mn, mx)
        ax1.set_ylim(mn, mx)

        #ax0.legend(loc='lower right', bbox_to_anchor=(-0.75, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        ax0.legend(loc='lower right', bbox_to_anchor=(0.89, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        ax1.legend(loc='lower right', bbox_to_anchor=(1.17, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        ax1.set_ylabel('sig(ESFAs)')#, fontsize = 'small')
        ax1.yaxis.label.set_color('tab:blue')
        #ax1.tick_params(labelsize='x-small')
        #ax1.legend(loc='lower right', bbox_to_anchor=(1.2, 0.05, 0.45, 2.5), fontsize = 'xx-small', framealpha=0.5)

        ax0.set_title('%s for high resolution reflections' %(maptype), fontsize = 'medium',fontweight="bold")
        #plt.show()
        plt.subplots_adjust(hspace=0.35, left=0.09, right=0.65, top = 0.95)
        #fig.tight_layout()
        plt.savefig('%s_sigmas.pdf' %(maptype),dpi=300, transparent=True)
        plt.savefig('%s_sigmas.png' %(maptype),dpi=300)
        plt.close()
        
    #Make the FoFo / sig(FoFo) plot
    mn = 0
    mx = 0
    plt.close()
    fig, ax0 = plt.subplots(1,1, figsize=(10, 5))
    ax1 = ax0.twinx()
    _, FoFo_type, _, bin_res_cent_lst, _, _, fdif_data_lst, fdif_sigmas_lst = alldata[0]
    ax0.plot(bin_res_cent_lst[1:], fdif_data_lst[1:], marker = '.', color = 'tab:red', label = '%s' %(FoFo_type))
    ax1.plot(bin_res_cent_lst[1:], fdif_sigmas_lst[1:], marker = 'x', linestyle='--', color = 'tab:blue', label = 'sig(%s)' %(FoFo_type))
    if min(fdif_data_lst[1:]) < mn:
        mn = min(fdif_data_lst[1:])
    if max(fdif_data_lst[1:]) > mx:
        mx = max(fdif_data_lst[1:])
    if min(fdif_sigmas_lst[1:]) < mn:
        mn = min(fdif_sigmas_lst[1:])
    if max(fdif_sigmas_lst[1:]) > mx:
        mx = max(fdif_sigmas_lst[1:])
                
    ax0.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax0.set_xlabel('Resolution (A)')#, fontsize = 'small')
    ax0.set_ylabel('%s' %(FoFo_type))#, fontsize = 'small')
    ax0.yaxis.label.set_color('tab:red')
    #ax0.tick_params(labelsize='x-small')
    #ax0.legend(loc='upper left', bbox_to_anchor=(-0.75, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    ax1.set_ylabel('sig(%s)' %(FoFo_type))#, fontsize = 'small')
    ax1.yaxis.label.set_color('tab:blue')
    #ax1.tick_params(labelsize='x-small')
    #ax1.legend(loc='center left', bbox_to_anchor=(-0.75, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    
    lines_labels_1 = [ax.get_legend_handles_labels() for ax in [ax0,ax1]]
    lines_1, labels_1 = [sum(lne, []) for lne in zip(*lines_labels_1)]
    ax0.legend(lines_1, labels_1, loc='lower right', bbox_to_anchor=(0.73, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    #ax0.legend(loc='lower right', bbox_to_anchor=(0.75, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)

    ax0.set_title('%s for high resolution reflections' %(FoFo_type),fontsize = 'medium',fontweight="bold")
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.85, top = 0.95)
    plt.savefig('%s_sigmas.pdf' %(FoFo_type),dpi=300, transparent=True)
    plt.savefig('%s_sigmas.png' %(FoFo_type),dpi=300)
    
    
def plot_Fextr_sigmas_allinone(pickle_file='Fextr_binstats.pickle', maptype_lst=['qFextr', 'Fextr', 'qFgenick', 'Fgenick', 'qFextr_calc', 'Fextr_calc'], occ_lst = None):
    """
    To plot the Values of the difference map and extrapolated structure factors in a single figure.
    Assuming that all data have the same resolution bin center. Should be clear whether this is the case from the tables generated by get_Fextr_stats.
    A list with maptypes and occupancies can be provided in order to gerenate custom plots if it would be run in an interactive mode.
    """
    
    plt.close()
    #fig, (ax0, ax2) = plt.subplots(1,2, figsize=(10, 5))
    fig = plt.figure(figsize=(10,5))
    ax0 = plt.subplot(211)
    ax1 = ax0.twinx()
    ax2 = plt.subplot(212)
    #ax3 = plt.subplot(224)
    ax3 = ax2.twinx()

    mn = 1000
    mx = 0

    with open(pickle_file,'rb') as stats_file:
        while True:
            try:
                stats = pickle.load(stats_file)
                try:
                    fdif_sigmas_lst
                    occ, FoFo_type, Fextr_type, _, fextr_data_lst, fextr_sigmas_lst, _,_ = stats
                except NameError:
                    occ, FoFo_type, Fextr_type, bin_res_cent_lst, fextr_data_lst, fextr_sigmas_lst, fdif_data_lst, fdif_sigmas_lst = stats
                if Fextr_type in maptype_lst:
                    if (occ_lst == None or occ in occ_lst):
                        if Fextr_type == 'qFextr':
                            color = plt.cm.Reds(int((np.log(1/occ))*100))
                        elif Fextr_type == 'Fextr':
                            color = plt.cm.Purples(int((np.log(1/occ))*100))
                        elif Fextr_type == 'qFgenick':
                            color = plt.cm.Blues(int((np.log(1/occ))*100))
                        elif Fextr_type == 'Fgenick':
                            color = plt.cm.Greens(int((np.log(1/occ))*100))
                        elif Fextr_type == 'qFextr_calc':
                            color = plt.cm.Oranges(int((np.log(1/occ))*100))
                        elif Fextr_type == 'Fextr_calc':
                            color = plt.cm.Greys(int((np.log(1/occ))*100))
                        else:
                            color = plt.cm.Greys(int((np.log(1/occ))*100))
                            
                        ax2.plot(bin_res_cent_lst[1:], fextr_data_lst[1:], marker = '.', color = color, label = '%s, occ = %.3f' %(Fextr_type, occ))
                        ax3.plot(bin_res_cent_lst[1:], fextr_sigmas_lst[1:], marker = 'x', linestyle='--', color = color, label = 'sig(%s), occ = %.3f' %(Fextr_type, occ))
                        
                        #Find min and max for plot, not very elegant
                        if min(fextr_data_lst[1:]) < mn:
                            mn = min(fextr_data_lst[1:])
                        if max(fextr_data_lst[1:]) > mx:
                            mx = max(fextr_data_lst[1:])
                        if min(fextr_sigmas_lst[1:]) < mn:
                            mn = min(fextr_sigmas_lst[1:])
                        if max(fextr_sigmas_lst[1:]) > mx:
                            mx = max(fextr_sigmas_lst[1:])
                        
            except EOFError:
                break
            
    ax0.plot(bin_res_cent_lst[1:], fdif_data_lst[1:], marker = '.', color = (0,1,0), label = '%s' %(FoFo_type))
    ax1.plot(bin_res_cent_lst[1:], fdif_sigmas_lst[1:], marker = 'x', linestyle='--', color = (0,0,1), label = 'sig(%s)' %(FoFo_type))
    
    ax0.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax2.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    #ax3.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    
    ax0.set_xlabel('Resolution (A)', fontsize = 'small')
    ax0.set_ylabel('%s' %(FoFo_type), fontsize = 'small')
    ax0.tick_params(labelsize='x-small')
    #ax0.legend(loc='upper left', bbox_to_anchor=(-0.75, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    ax1.set_ylabel('sig(%s)' %(FoFo_type), fontsize = 'small')
    ax1.tick_params(labelsize='x-small')
    #ax1.legend(loc='center left', bbox_to_anchor=(-0.75, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    ax0.set_title('%s for high resolution reflections' %(FoFo_type),fontsize = 'medium',fontweight="bold")
    
    lines_labels_1 = [ax.get_legend_handles_labels() for ax in [ax0,ax1]]
    lines_1, labels_1 = [sum(lne, []) for lne in zip(*lines_labels_1)]
    ax0.legend(lines_1, labels_1, loc='lower right', bbox_to_anchor=(0.82, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)

    ax2.set_xlabel('Resolution (A)', fontsize = 'small')
    ax2.set_ylabel('Extrapolated structure factors', fontsize = 'small')
    ax2.tick_params(labelsize='x-small')
    
    ax2.set_ylim(mn, mx)
    ax3.set_ylim(mn, mx)

    #ax2.legend(loc='lower right', bbox_to_anchor=(-0.75, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    ax2.legend(loc='lower right', bbox_to_anchor=(0.82, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    ax3.set_ylabel('sig(Extrapolated structure factors)', fontsize = 'small')
    ax3.tick_params(labelsize='x-small')
    #ax3.legend(loc='lower right', bbox_to_anchor=(1.2, 0.05, 0.45, 2.5), fontsize = 'xx-small', framealpha=0.5)

    ax2.set_title('%s for high resolution reflections' %(Fextr_type), fontsize = 'medium',fontweight="bold")
    #plt.show()
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.75, top = 0.95)
    #fig.tight_layout()
    plt.savefig('FoFo_Fextr_sigmas.pdf',dpi=300, transparent=True)
    plt.savefig('FoFo_Fextr_sigmas.png',dpi=300)
    plt.close()
    
def plot_sigmas_allinone(pickle_file='Fextr_binstats.pickle', maptype_lst=['qFextr', 'Fextr', 'qFgenick', 'Fgenick', 'qFextr_calc', 'Fextr_calc'], occ_lst = None):
    """
    To plot the Values of the sigmas. To understand when the error on the extrapolated structure factors becomes larger/smaller than the error in the difference structure factors.
    Assuming that all data have the same resolution bin center. Should be clear whether this is the case from the tables generated by get_Fextr_stats.
    A list with maptypes and occupancies can be provided in order to gerenate custom plots if it would be run in an interactive mode.
    """
    plt.close()
    #fig,ax0 = plt.subplots(figsize=(10, 5))

    mn = 1000
    mx = 0
    
    if len(maptype_lst)>1:
        n_cols = 2
    else:
        n_cols = 1 
    n_rows = int(math.ceil(len(maptype_lst)/2))
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(10*n_rows, 3*n_cols))

    with open(pickle_file,'rb') as stats_file:
        while True:
            try:
                stats = pickle.load(stats_file)
                try:
                    fdif_sigmas_lst
                    occ, FoFo_type, Fextr_type, _, _, fextr_sigmas_lst, _,_ = stats
                except NameError:
                    occ, FoFo_type, Fextr_type, bin_res_cent_lst, _, fextr_sigmas_lst, _, fdif_sigmas_lst = stats
                
                if Fextr_type in maptype_lst:
                    if (occ_lst == None or occ in occ_lst):
                        if Fextr_type == 'qFextr':
                            color = plt.cm.Reds(int((np.log(1/occ))*100))
                        elif Fextr_type == 'Fextr':
                            color = plt.cm.Purples(int((np.log(1/occ))*100))
                        elif Fextr_type == 'qFgenick':
                            color = plt.cm.Blues(int((np.log(1/occ))*100))
                        elif Fextr_type == 'Fgenick':
                            color = plt.cm.Greens(int((np.log(1/occ))*100))
                        elif Fextr_type == 'qFextr_calc':
                            color = plt.cm.Oranges(int((np.log(1/occ))*100))
                        elif Fextr_type == 'Fextr_calc':
                            color = plt.cm.Greys(int((np.log(1/occ))*100))
                        else:
                            color = plt.cm.Greys(int((np.log(1/occ))*100))
                            
                        if n_cols == 1:
                            axs.plot(bin_res_cent_lst[1:], fextr_sigmas_lst[1:], marker = 'x', linestyle='--', color = color, label = 'occ = %.3f' %(occ))
                        elif n_rows == 1:
                            axs[maptype_lst.index(Fextr_type)].plot(bin_res_cent_lst[1:], fextr_sigmas_lst[1:], marker = 'x', linestyle='--', color = color, label = 'occ = %.3f' %(occ))
                        else:
                            axs[int(math.floor(maptype_lst.index(Fextr_type)/2)), maptype_lst.index(Fextr_type)%2].plot(bin_res_cent_lst[1:], fextr_sigmas_lst[1:], marker = 'x', linestyle='--', color = color, label = 'occ = %.3f' %(occ))
                        
                        if min(fextr_sigmas_lst[1:]) < mn:
                            mn = min(fextr_sigmas_lst[1:])
                        if max(fextr_sigmas_lst[1:]) > mx:
                            mx = max(fextr_sigmas_lst[1:])
                        
            except EOFError:
                break
            
    if min(fdif_sigmas_lst[1:]) < mn:
        mn = min(fdif_sigmas_lst[1:])
    if max(fdif_sigmas_lst[1:]) > mx:
        mx = max(fdif_sigmas_lst[1:])
    
    if n_cols == 1:
        axs.plot(bin_res_cent_lst[1:], fdif_sigmas_lst[1:], marker = 'x', linestyle = '-',  linewidth = 2.5, color = (0,0,1), label = 'sig(%s)' %(FoFo_type))
        axs.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
        axs.set_xlabel('Resolution (A)', fontsize = 'small')
        axs.set_ylabel('sigma value', fontsize = 'small')
        axs.set_ylim(mn, mx)
        axs.tick_params(labelsize='x-small')
        axs.set_title('sigma(%s) and sigma(%s) for high resolution reflections' %(FoFo_type, maptype_lst[0]),fontsize = 'medium',fontweight="bold")
        axs.legend(loc='lower right', bbox_to_anchor=(0.70, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        plt.subplots_adjust(hspace=0.35, left=0.09, right=0.85, top = 0.85)
    elif n_rows == 1:
        for Fextr_type in maptype_lst:
            axs[maptype_lst.index(Fextr_type)].plot(bin_res_cent_lst[1:], fdif_sigmas_lst[1:], marker = 'x', linestyle = '-',  linewidth = 2.5, color = (0,0,1), label = 'sig(%s)' %(FoFo_type))
            axs[maptype_lst.index(Fextr_type)].set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
            axs[maptype_lst.index(Fextr_type)].set_xlabel('Resolution (A)', fontsize = 'small')
            axs[maptype_lst.index(Fextr_type)].set_ylabel('sigma value', fontsize = 'small')
            axs[maptype_lst.index(Fextr_type)].set_ylim(mn, mx)
            axs[maptype_lst.index(Fextr_type)].tick_params(labelsize='x-small')
            axs[maptype_lst.index(Fextr_type)].set_title('sigma(%s) and sigma(%s)\n for high resolution reflections' %(FoFo_type, Fextr_type),fontsize = 'medium',fontweight="bold")
            axs[maptype_lst.index(Fextr_type)].legend(loc='lower right', bbox_to_anchor=(0.82, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
            plt.subplots_adjust(hspace=0.35, left=0.09, right=0.85, top = 0.85, wspace=0.35)
    else:
        for Fextr_type in maptype_lst:
            axs[int(math.floor(maptype_lst.index(Fextr_type)/2)), maptype_lst.index(Fextr_type)%2].plot(bin_res_cent_lst[1:], fdif_sigmas_lst[1:], marker = 'x', linestyle = '-',  linewidth = 2.5, color = (0,0,1), label = 'sig(%s)' %(FoFo_type))
            axs[int(math.floor(maptype_lst.index(Fextr_type)/2)), maptype_lst.index(Fextr_type)%2].set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
            axs[int(math.floor(maptype_lst.index(Fextr_type)/2)), maptype_lst.index(Fextr_type)%2].set_xlabel('Resolution (A)', fontsize = 'small')
            #axs[int(math.floor(maptype_lst.index(Fextr_type)/2)), maptype_lst.index(Fextr_type)%2].tick_params(labelbottom=False)
            axs[int(math.floor(maptype_lst.index(Fextr_type)/2)), maptype_lst.index(Fextr_type)%2].set_ylabel('sigma value', fontsize = 'small')
            axs[int(math.floor(maptype_lst.index(Fextr_type)/2)), maptype_lst.index(Fextr_type)%2].set_ylim(mn, mx)
            axs[int(math.floor(maptype_lst.index(Fextr_type)/2)), maptype_lst.index(Fextr_type)%2].tick_params(labelsize='x-small')
            axs[int(math.floor(maptype_lst.index(Fextr_type)/2)), maptype_lst.index(Fextr_type)%2].set_title('sigma(%s) and sigma(%s) for high resolution reflections' %(FoFo_type, Fextr_type),fontsize = 'medium',fontweight="bold")
            axs[int(math.floor(maptype_lst.index(Fextr_type)/2)), maptype_lst.index(Fextr_type)%2].legend(loc='lower right', bbox_to_anchor=(0.65, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
            plt.subplots_adjust(hspace=0.55, left=0.05, right=0.95, top = 0.90, wspace=0.15, bottom = 0.1)
    ##plt.show()
    
    #fig.tight_layout()
    plt.savefig('FoFo_sigmas.pdf',dpi=300, transparent=True)
    plt.savefig('FoFo_sigmas.png',dpi=300)
    #plt.close()
    

def dump_negative_stats(occ, maptype, all_reflections, neg_reflections, outdir):
    """
    Store number of total, positive and negative reflections in pickle file to plot them afterwards. Useage of pickle file avoids storing stuff in memory and be able to acces the data afterwards.
    """
    pos_reflections = all_reflections - neg_reflections
    stats = [occ, maptype, all_reflections, pos_reflections, neg_reflections]
    
    out=open('%s/Fextr_negative.pickle' %(outdir),'ab')
    pickle.dump(stats,out)
    out.close()
    print("Number of negative reflections saved to %s/Fextr_negative.pickle" %(outdir))
    
def plot_negative_reflections(pickle_file='Fextr_negative.pickle'):
    """
    Read pickle file with number of negative reflections and plot them in a seperate figure for each type of extrapolated structure factors
    """
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

    maptype_lst=list(set(alldata[:,1]))
    for maptype in maptype_lst:
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
        
        ax0.set_xlabel('Triggered state occupancy')
        ax0.set_ylabel('Absolute number')
        ax1.set_ylabel('Percentage of total number of ESFAs')
        ax0.set_title('%s: Negative ESFAs' %(maptype), fontsize = 'medium',fontweight="bold")
        
        ax2.set_xlabel('Triggered state occupancy')
        ax2.set_ylabel('Absolute number')
        ax3.set_ylabel('Percentage of total number of ESFAs')
        ax2.set_title('%s: Positive ESFAs'%(maptype), fontsize = 'medium',fontweight="bold")
        
        #plt.title('%s' %(maptype),loc='center')
        plt.subplots_adjust(hspace=0.25,wspace=0.5, left=0.09, right=0.88, top = 0.95)
        #fig.tight_layout()
        plt.savefig("Neg_Pos_reflections_%s.pdf" %(maptype), dpi=300, transparent=True)
        plt.savefig("Neg_Pos_reflections_%s.png" %(maptype), dpi=300)
        plt.close()

        
def plot_negative_reflections_allinone(pickle_file='Fextr_negative.pickle'):
    """
    Read pickle file with number of negative reflections and plot them on the same figure.
    """
    
    plt.close()
    fig, (ax0, ax2) = plt.subplots(1,2, figsize=(10, 5))
    ax1 = ax0.twinx()
    ax3 = ax2.twinx()
    
    mn = 100000
    mx = 0
    
    with open(pickle_file,'rb') as stats_file:
        while True:
            try:
                stats = pickle.load(stats_file)
                occ, maptype, all_reflections, pos_reflections, neg_reflections = stats
                #print(occ, maptype, neg_reflections)
                if maptype == 'qFextr':
                    color='red'
                elif maptype == 'Fextr':
                    color='purple'
                elif maptype == 'qFgenick':
                    color='blue'
                elif maptype == 'Fgenick':
                    color='green'
                elif maptype == 'qFextr_calc':
                    color='orange'
                elif maptype == 'Fextr_calc':
                    color='gray'
                else:
                    color='black'
                ax0.plot(occ, neg_reflections, marker='o', color=color, label = maptype)
                ax1.plot(occ, (neg_reflections/all_reflections)*100, marker = 'x', color = color)
                ax2.plot(occ, pos_reflections, marker='o', color=color, label = maptype)
                ax3.plot(occ, (pos_reflections/all_reflections)*100, marker = 'x', color = color)
                if neg_reflections > mx:
                    mx = neg_reflections
                if pos_reflections > mx:
                    mx = pos_reflections
                if neg_reflections < mn:
                    mn = neg_reflections
                if pos_reflections < mn:
                    mn = pos_reflections
            except EOFError:
                break
    
    ax0.set_ylim(mn,mx)
    ax1.set_ylim(0,100)
    ax2.set_ylim(mn,mx)
    ax3.set_ylim(0,100)
    
    #ax0.legend(fontsize = 'x-small', framealpha=0.5, loc=1, bbox_to_anchor=(0.05, 0.5, 0.5, 0.45))
    #get unique labels
    lines, labels = ax0.get_legend_handles_labels()
    i = np.arange(len(labels))
    f = np.array([])
    unique_labels = list(set(labels))
    for ul in unique_labels:
        f = np.append(f,[i[np.array(labels)==ul][0]])
    lines = [lines[int(j)] for j in f]
    labels = [labels[int(j)] for j in f]
    ax2.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.92, -0.05, 0.45, 0.5), fontsize = 'x-small', framealpha=0.5)

    ax0.set_xlabel('Triggered state occupancy')
    ax0.set_ylabel('Number of negative reflections')
    ax1.set_ylabel('Percentage of total number of reflections (x)')
    ax0.set_title('Negative reflections', fontsize = 'medium',fontweight="bold")
    
    ax2.set_xlabel('Triggered state occupancy')
    ax2.set_ylabel('Number of positive reflections')
    ax3.set_ylabel('Percentage of total number of reflections (x)')
    ax2.set_title('Positive reflections', fontsize = 'medium',fontweight="bold")
    
    #plt.title('Number of negative and positive reflections',loc='center')
    plt.subplots_adjust(hspace=0.25,wspace=0.5, left=0.09, right=0.88, top = 0.95)
    #fig.tight_layout()
    plt.savefig("Neg_Pos_reflections.pdf", dpi=300, transparent=True)
    plt.savefig("Neg_Pos_reflections.png", dpi=300)
    plt.close()

def calculate_R(f_obs, f_calc): #to avoid a scale factor being taken into account
    #return abs(np.sum(abs(f_obs.data())-abs(f_calc.data())))/np.sum(abs(f_obs.data()))
    """
    Calculates normal R-factor frac(sum(||F| - k|F'||))(\sum(|F|)). For Riso, should be better to use frac(sum(||F| - k|F'||))(\sum((|F|+|F'|)/2))
    """
    return f_obs.r1_factor(f_calc)

def calculate_Riso(f_obs1, f_obs2):
    """
    Calculation of Riso, as done in phenix.fobs-fobs
    """
    num = flex.sum(flex.abs(f_obs1.data()-f_obs2.data()))
    den = flex.sum(flex.abs(f_obs1.data()+f_obs2.data())/2)
    r = None
    if(den!=0):
        r = num/den
    return r
    

def compute_r_factors(f_obs, f_calc, r_free_flags, log=sys.stdout):
    """
    R-factorcan can be calculated using miller-build in function .r1_factor (with or without emulate_sftools=True
    Using .r1_factor, R is calculated as R1 = frac(sum(||F| - k|F'||))(\sum(|F|)) where F is self.data() and F' is other.data() and k is the factor to put F' on the same scale as F.
    For calculation of Riso values:
        - Rfree and ccfree probably don't tell you anything
        - Riso the denominator is different as compared to standard R-factors. Nevertheless, the standard R-factor should be a proper approximation
    TODO: suggest resolution cutof based on Riso and CCiso values
    """
    if str(type(f_calc.data())) != "<class 'scitbx_array_family_flex_ext.double'>":
        f_calc = miller.array(miller_set = f_calc, data = f_calc.amplitudes().data())
            
    #f_obs_work, f_obs_test   = make_fwork_ffree(f_obs, r_free_flags)
    #f_calc_work, f_calc_test = make_fwork_ffree(f_calc, r_free_flags)
    
    #r_work = calculate_R(f_obs_work, f_calc_work)
    #r_free = calculate_R(f_obs_test, f_calc_test)
    #cc_work = f_obs_work.correlation(f_calc_work).coefficient()
    #cc_free = f_obs_test.correlation(f_calc_test).coefficient()
    
    f_obs_work  = f_obs
    f_calc_work = f_calc
    
    #R adn cc with all data. Just call it work to avoid to much changes below
    r_work = calculate_Riso(f_obs_work, f_calc_work)
    cc_work = f_obs_work.correlation(f_calc_work).coefficient()

    #print("r_work = %.4f  r_free = %.4f" %(r_work, r_free))
    #print("cc_work = %.4f cc_free = %.4f" %(cc_work, cc_free))
    
    bin_res_cent_lst = flex.double()
    r_work_lst  = flex.double()
    #r_free_lst  = flex.double()
    cc_work_lst = flex.double()
    #cc_free_lst = flex.double()
    f_obs_work.setup_binner(n_bins=20)
    #f_obs_test.use_binning_of(f_obs_work)
    f_calc_work.use_binner_of(f_obs_work)
    #f_calc_test.use_binning_of(f_obs_work)
    #print("bin  resolution range  #refl.work  #refl.test  riso-work riso-free  cciso-work  cciso-free", file=log)
    #print("bin  resolution range  #refl.work  #refl.test  riso-work riso-free  cciso-work  cciso-free")
    
    print("\n************isomorphism statistics************",file =log)
    print("\n************isomorphism statistics************")
    print("bin  resolution range         #refl  R_iso      cc_iso", file=log)
    print("bin  resolution range         #refl  R_iso      cc_iso")
    for i_bin in f_obs_work.binner().range_all():
        sel_work = f_obs_work.binner().selection(i_bin)
        #sel_test = f_obs_test.binner().selection(i_bin)
        fo_work_bin = f_obs_work.select(sel_work)
        fc_work_bin = f_calc_work.select(sel_work)
        #fo_test_bin = f_obs_test.select(sel_test)
        #fc_test_bin = f_calc_test.select(sel_test)
        if fc_work_bin.size() == 0 : continue
        bin_res_cent = np.median(f_obs_work.binner().bin_d_range(i_bin))
        bin_res_cent_lst.append(bin_res_cent)
        #r_work_bin = calculate_R(fo_work_bin, fc_work_bin)
        r_work_bin = calculate_Riso(fo_work_bin, fc_work_bin)
        r_work_lst.append(r_work_bin)
        #r_free_bin = calculate_R(fo_test_bin, fc_test_bin)
        #r_free_lst.append(r_free_bin)
        cc_work_bin = fo_work_bin.correlation(fc_work_bin).coefficient()
        cc_work_lst.append(cc_work_bin)
        #cc_free_bin = fo_test_bin.correlation(fc_test_bin).coefficient()
        #cc_free_lst.append(cc_free_bin)
        legend = f_obs_work.binner().bin_legend(i_bin, show_counts=False)
        #print("%s %8d  %8d   %.4f      %.4f       %.3f      %.3f" % (legend, fo_work_bin.size(),
        #fo_test_bin.size(), r_work_bin, r_free_bin, cc_work_bin, cc_free_bin))
        #print("%s %8d  %8d   %.4f      %.4f       %.3f      %.3f" % (legend, fo_work_bin.size(),
        #fo_test_bin.size(), r_work_bin, r_free_bin, cc_work_bin, cc_free_bin), file=log)
        print("%s %8d  %.4f       %.3f" % (legend, fo_work_bin.size(), r_work_bin, cc_work_bin), file=log)
        print("%s %8d  %.4f       %.3f" % (legend, fo_work_bin.size(), r_work_bin, cc_work_bin))
        
    #assert bin_res_cent_lst.size()==r_work_lst.size()==r_free_lst.size()==cc_work_lst.size()==cc_free_lst.size(), 'list sizes to plot Riso and CCiso not equal'
    assert bin_res_cent_lst.size()==r_work_lst.size()==cc_work_lst.size(), 'list sizes to plot Riso and CCiso not equal'
    
    out=open('Riso_CCiso.pickle' ,'wb') #write to pickle for GUI
    stats = [bin_res_cent_lst, r_work_lst, cc_work_lst, r_work, cc_work]
    pickle.dump(stats,out)
    out.close()

    
    fig,ax1 = plt.subplots(figsize=(10, 5))
    ax1.plot(bin_res_cent_lst[1:], r_work_lst[1:], marker = '.', color = 'red', label = 'Riso; overall %.4f' %(r_work))
    ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
    ax1.set_xlabel('Resolution (A)')
    ax1.set_ylabel('Riso')
    ax1.yaxis.label.set_color('red')
    #ax1.legend(fontsize = 'xx-small', framealpha=0.5, loc=0)
    ax2 = ax1.twinx()
    #ax2.plot(bin_res_cent_lst[1:], cc_work_lst[1:], marker = '.', color = 'green', label = 'CCiso,work; overall %.4f' %(cc_work))
    #ax2.plot(bin_res_cent_lst[1:], cc_free_lst[1:], marker = '.', color = 'yellow', label = 'CCiso,free; overall %.4f' %(cc_free))
    ax2.plot(bin_res_cent_lst[1:], cc_work_lst[1:], marker = '^', markersize = 5, color = 'green', label = 'CCiso; overall %.4f' %(cc_work))
    ax2.set_ylabel('CCiso')
    ax2.yaxis.label.set_color('green')
    #ax2.legend(fontsize = 'xx-small', framealpha=0.5, loc=0)
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lne, []) for lne in zip(*lines_labels)]
    #fig.legend(lines, labels, loc='lower right', fontsize = 'xx-small', framealpha=0.5) #bbox_to_anchor=(0.82, -0.05, 0.45, 0.5)
    ax2.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    plt.title('Riso and CCiso for high resolution reflections', fontsize = 'medium',fontweight="bold")
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.80, top = 0.95)
    #fig.tight_layout()
    plt.savefig('Riso_CCiso.pdf', dpi=300, transparent=True)
    plt.savefig('Riso_CCiso.png', dpi=300)
    plt.close()

    #return r_work, r_free
    return r_work, cc_work


def plot_F1_F2(F1, F2, F1_name = "Data_set_1", F2_name = "Data_set_2", log = sys.stdout):
    """
    Calculate the Pearson correlation coefficient between two data sets and plot them
    """
    
    #Calculate the correlation coefficient 
    #using the correlation function of cctbx.miller (avoids the need for equally sized data sets)
    cc = F1.correlation(F2).coefficient()
    
    #for plotting, only the common reflections can be used and the order should be the same
    if list(F1.indices()[-10:])==list(F2.indices()[-10:]) == False:
        #use cctbx.miller to maintain the common indices, this should also sort them
        F1, F2 = F1.common_sets(F2)
    if list(F1.indices()[-10:])==list(F2.indices()[-10:]) == False:
        #something went wrong in extracting the common indices
        print("Indices cannot be set equivalent, the scatter plot will not be made")
    else:
        ##get standard deviation in order to correct slope for linear fit
        #s1 = flex.mean_and_variance(F1.data()).unweighted_sample_standard_deviation()
        #s2 = flex.mean_and_variance(F2.data()).unweighted_sample_standard_deviation()
        ##calculate the slope of the linear fit
        #b  = cc * s2/s1
        ##calculate intercept of the linear fit
        ##a = F2.data()[0] - F1.data()[0] * b
        ##a = 0
        
        #use linear regression because a cannot be calculated
        b, a, _, _, _ = linregress(F1.data(), F2.data())
        
        plt.close()
        mn = 0
        mx = np.max([np.max(F1.data()), np.max(F2.data())])
        #plot the data
        fig,ax1 = plt.subplots(figsize=(10, 5))
        ax1.scatter(F1.data(), F2.data(), marker = '.', color = 'red', label = "Pearson correlation\n coefficient: %.4f" %(cc))
        #plot the linear fit, we only need two values for plotting a line
        lst = np.array([np.min(F1.data()), np.max(F1.data())])
        ax1.plot(lst, lst * b + a,  color = 'blue', label = "linear fit: y = %.2f x + %.2f" %(b, a))
        ax1.set_xlim(mn,mx)
        ax1.set_ylim(mn,mx)
        ax1.set_xlabel(F1_name)
        ax1.set_ylabel(F2_name)
        ax1.legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        plt.title('Correlation between %s and %s' %(F1_name, F2_name), fontsize = 'medium',fontweight="bold")
        plt.subplots_adjust(hspace=0.35, left=0.09, right=0.80, top = 0.95)
        plt.savefig("correlation_%s_%s.pdf"%(F1_name, F2_name), dpi=300, transparent=True)
        #plt.savefig("correlation_%s_%s.png"%(F1_name, F2_name), dpi=300, transparent=True)
        plt.close()
        
    return cc
    
def plot_correlations(occ_lst, correlation_list):
    """
    Plot the correlations calculated in plot_F1_F2 in function of the occupancy
    """
    assert len(occ_lst) == len(correlation_list), "list with occupancies and correlations not of equal length"
    
    alphas = list(map(lambda x: round(1/x, 3), occ_lst))
        
    plt.close()
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10, 5))
    ax1.plot(occ_lst, correlation_list, marker = 'o', color = 'red', label = "Correlation")
    #plot the linear fit, we only need two values for plotting a line
    ax1.set_xlabel("Triggered state occupancy")
    ax1.set_ylabel("Pearson correlalation coefficient")
    
    ax2.plot(alphas, correlation_list, marker = 'o', color = 'red', label = "Correlation")
    #plot the linear fit, we only need two values for plotting a line
    ax2.set_xlabel("Alpha")
    ax2.set_ylabel("Pearson correlalation coefficient")    
    
    #ax1.legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
    plt.title('Correlation between reference and extrapolated structure factors', fontsize = 'medium',fontweight="bold")
    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.80, top = 0.95)
    plt.savefig("correlations_per_alpha.pdf", dpi=300, transparent=True)
    #plt.savefig("correlations_per_alpha.png", dpi=300, transparent=True)
    plt.close()
    
    
    
        
