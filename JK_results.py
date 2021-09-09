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
import multiprocessing
from functools import partial
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.axislines import Subplot
import matplotlib.ticker as mtick
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

from Fextr_utils import *
import Fextr
from Log_file import print_terminal_and_log
import os
import numpy as np

def get_JK_results(tab_total, tab_list, outdir, ordered_solvent):
    '''
    Global function calculating the RMSD and CC to create the tables and the plot to compare the JK and total
    Args:
        tab_total: table for total
        tab_list: table for all JKs
        outdir: directory where to put the file with all JK results
        ordered_solvent: bool, variable from the phil

    Returns:
    file with CC and RMSD tables and plots
    '''
    os.chdir(outdir) #change directory to the output directory of the JK_results

#0. Use waters for comparison?
    use_waters = not ordered_solvent #if ordered solvent, the water in the models will not be used in comparison #TODO: check if this line is correct or if it has to be the opposit
    if use_waters:
        print_terminal_and_log('The water molecules will not be used in the comparison between models')
    else:
        print_terminal_and_log('The water molecules will be used in the comparison between models')

    tab_stack=np.stack(tab_list, axis=2) #stack all the table of files and data

#1. Create new file to print results
    #print JK results in a new log file
    JK_results_file = outdir + '/JK_results.log'
    JK_results_file_open = open(JK_results_file, 'w')
    print_terminal_and_log('The JK results can be found in %s'%(JK_results_file))

#2. Extract files from the X8 tables and get CC and RMSD
    #Get parameters and list of parameters
    for i in range (1,len(tab_stack[:,0,0])): #get the number of lines
        print_terminal_and_log('--- %s : Getting file names and variables ---' %(i))
        occ = float(tab_stack[i,0,0])#occupancy
        maptype = tab_stack[i,1,0]#map type
        refinement = tab_total[i][2]#variable for refinement:(None or False or real or rec or real+rec) of the total
        labels = tab_stack[i,5,0]#labels of the mtz file to get the average of

        #get total files
        if tab_total[i][2] != False: #if the refinement is not incomplete (None or real or rec or real+rec)
            total_map = tab_total[i][3] #mtz file created with all images
            total_model = tab_total[i][4] #pdb file created with all images
            totalresidlist_zscore = tab_total[1][7] #residue list created with all images
        else:
            total_map = 'None'
            total_model = 'None'
            totalresidlist_zscore = 'None'

        #init lists
        list_of_maps = [] #list of existing maps (without False refinement) filled if total_map exists
        list_of_models = [] #list of existing models (without False refinement) filled if total_model exists
        rmsd_total_fitted_list = [] #list of global RMSD (without False refinement) filled if total_model exists
        rmsd_residues_list = [] #list of residue RMSD (without False refinement) filled if total_model exists
        # assign variable with impossible value
        cc_total=-1

        if tab_total[i][2] != False:  # if the refinement is not incomplete (None or real or rec or real+rec)
#        if os.path.isfile(total_map) and os.path.isfile(total_model):  # if the total map and model to compare exist
            print_terminal_and_log('--- %s : Getting CC and RMSD ---' % (i))

            print_terminal_and_log('\nTable %s: occ=%s maptype=%s refinement=%s total_map=%s total_model=%s \n------------------------------------------------------------------------------------ ' % (
                    i, occ, maptype, refinement, total_map , total_model), log=JK_results_file_open)
            print_terminal_and_log('{0:<9s} {1:<27s} {2:<30s} {3:<31s} {4:<15s}'.format('JK_number',
                                                                    'total_RMSD_totalmodel-model',
                                                                    'residues_RMSD_totalmodel-model',
                                                                    'number_of_atoms_in_residue_list',
                                                                    'cc_totalmap-map'), log=JK_results_file_open)
            for j in range (0, len(tab_stack[0,0,:])): #get the number of JK
                #get parameters
                other_residlist_zscore = totalresidlist_zscore #the residue list used will be the one of the total, alternative=tab_stack[1,7,j] #residue list of the JK
                print_terminal_and_log('The list of residues of the total file will be used to calculate the RMSD of the atoms of the residue list')
                other_map = tab_stack[i,3,j] #map file
                other_model = tab_stack[i,4,j] #model file
                refinement_j = tab_stack[i,2,j] #variable for refinement:(None or False or real or rec or real+rec)
                JK_nb = tab_stack[i,6,j] #number of JK

                if refinement_j != False: #if the refinement is not incomplete (None or real or rec or real+rec)
                    #assign variables with impossible value
                    rmsd_total_fitted = -1
                    rmsd_residues = -1
                    nb_residlst = -1
                    cc = -1

                    # To avoid errors of type, changing None into str(None)
                    if other_map==None:
                        other_map='None'
                    if other_model == None:
                        other_model = 'None'
                    if total_map == None:
                        total_map = 'None'
                    if total_model == None:
                        total_model = 'None'

                #Get CC
                    if os.path.isfile(total_map):
                        if  os.path.isfile(other_map):
                            list_of_maps.append(other_map)
                            # get correlation coefficient between the map with all images and the other map
                            cc = float(get_cc_mtz_mtz(other_map, total_map, outdir)) #getting CC between the two maps
                        else:
                            print_terminal_and_log('The map for occ=%s, maptype=%s, JK=%s does not exist, the total RMSD can not be calculated' % (occ, maptype, JK_nb))
                    else:
                        print_terminal_and_log('The total map for occ=%s, maptype=%s does not exist, the total RMSD can not be calculated' % (occ, maptype))

                #Get RMSD
                    if os.path.isfile(total_model):
                        if os.path.isfile(other_model):
                            list_of_models.append(other_model)
                            # get the total rmsd between the total model and the JK model
                            rmsd_total_fitted, other_pdb_fitted = superpose_pdbs(total_model, other_model)  # get rmsd between other_pdb and total_pdb after other_pdb is fitted to total_pdb and get the fitted other pdb

                            # print('rmsd with superpose pdb', rmsd_total_fitted) #for information: the superpose pdb uses the waters
                            #
                            # #TODO: check if the rmsd calculations are correct
                            # residlist_all_total = get_residlist_from_all_atoms(total_model,
                            #                                                    use_waters=use_waters)
                            # coord_all_total, _ = get_coord_from_parser_selected(total_model,
                            #                                                     residlist_all_total)
                            # coord_all_total = np.array(coord_all_total)
                            # residlist_all = get_residlist_from_all_atoms(other_pdb_fitted,
                            #                                              use_waters=use_waters)
                            # coord_all, _ = get_coord_from_parser_selected(other_pdb_fitted,
                            #                                               residlist_all)
                            # coord_all = np.array(coord_all_total)
                            # rmsd_total_fitted = calculate_rmsd(coord_all_total, coord_all)
                            # print('rmsd with calcul', rmsd_total_fitted)


                            rmsd_total_fitted_list.append(rmsd_total_fitted)

                            if os.path.isfile(other_residlist_zscore) and os.path.isfile(totalresidlist_zscore): #if the residue list files exist
                                # get rmsd of the residue list between the total model and the JK model
                                total_residlst = get_residlist(totalresidlist_zscore, use_waters = use_waters)  # get unique residue list from the total_residue_list file
                                coord_total_pdb_residlst, info_total_pdb_residlst = get_coord_from_parser_selected(total_model, total_residlst)  # get coordinates of atoms of residlst from total_pdb
                                other_residlst = get_residlist(other_residlist_zscore, use_waters = use_waters)  # get unique residue list from the other_residue_list file
                                coord_other_pdb_residlst, info_other_pdb_residlst = get_coord_from_parser_selected(other_pdb_fitted, other_residlst)  # get coordinates of atoms of residlst from other_pdb
                                coord_total_common, _, coord_other_common, _, nb_residlst = select_common_coords(
                                    coord_total_pdb_residlst, info_total_pdb_residlst, coord_other_pdb_residlst,
                                    info_other_pdb_residlst)  # get list of coordinates of common atoms from other_pdb and total_pdb

                                rmsd_residues = calculate_rmsd(coord_total_common, coord_other_common)  # calculate rmsd between total_pdb and other_pdb for residue list
                                rmsd_residues_list.append(rmsd_residues)

                            else:
                                print_terminal_and_log('The residue list file %s does not exist, the RMSD can not be calculated for the residues' % (other_residlist_zscore))
                        else:
                            print_terminal_and_log('The model for occ=%s, maptype=%s, JK=%s does not exist, the total RMSD can not be calculated' % (occ, maptype, JK_nb))
                    else:
                        print_terminal_and_log('The total model for occ=%s, maptype=%s does not exist, the total RMSD can not be calculated' % (occ, maptype))

                #Print RESULTS
                    print_terminal_and_log('{0:<9s} {1:> 27.4f} {2:> 30.4f} {3:> 31d} {4:> 8.4f}'.format(str(JK_nb),
                                                                                         rmsd_total_fitted,
                                                                                         rmsd_residues,
                                                                                         nb_residlst,
                                                                                         float(cc)), log=JK_results_file_open)
                else:
                    print_terminal_and_log('The refinements did not work for occ=%s, maptype=%s, JK=%s. The file can not be used to calculate the CC, RMSD and the average map'%(occ, maptype, JK_nb))

        else:
            print_terminal_and_log('The refinement did not succeed for occ=%s, maptype=%s'%(occ, maptype))
            #print('The total map and model do not exist for occ=%s, maptype=%s'%(occ, maptype))

    #Get total-average CC
        if len(list_of_maps)<2 : #Check if there are enough maps to do the average
            print_terminal_and_log('There are less than 2 maps coming from Jack Knife for occ=%s, %s; the average of the maps can not be done'%(occ, maptype))
        else:
            # get the average map from list of maps
            avg_map = average_maps(list_of_maps, 'occ%s_%s_%srefine' % (occ, maptype, refinement), labels) #TODO: Average map that can be used for real space refinement and simulated annealing
            print_terminal_and_log('the average map of %s is %s' % (list_of_maps, avg_map))  # , file=JK_results_open)
            # get correlation coefficient between the map with all images and the average map
            cc_total = float(get_cc_mtz_mtz(avg_map, total_map, outdir))

    #Get mean and standard deviation on RMSD
        mean_rmsd_total_fitted_list = np.mean(list(rmsd_total_fitted_list))  # mean deviation of the global rmsd between total model and other models
        std_rmsd_total_fitted_list = np.std(list(rmsd_total_fitted_list))  # standard deviation of the global rmsd between total model and other models
        mean_rmsd_residues_list = np.mean(list(rmsd_residues_list)) # mean deviation of the residue rmsd between total model and other models
        std_rmsd_residues_list = np.std(list(rmsd_residues_list)) # standard deviation of the residue rmsd between total model and other models

    # Print RESULTS
        print_terminal_and_log('{0:^9s} {1:^27s} {2:^30s} {3:^31s} {4:^15s}'.format(' ',
                                                                     'mean',
                                                                     'mean',
                                                                     ' ',
                                                                     'cc_total-averagemap'), log=JK_results_file_open)
        print_terminal_and_log('{0:<9s} {1:^ 27.4f} {2:^ 30.4f} {3:^31s} {4:^8.4f}'.format(str(JK_nb),
                                                                           mean_rmsd_total_fitted_list,
                                                                           mean_rmsd_residues_list,
                                                                           ' ',
                                                                           float(cc_total)), log=JK_results_file_open)
        print_terminal_and_log('{0:^9s} {1:^27s} {2:^30s} {3:^31s} {4:^15s}'.format(' ',
                                                                     'StandDev:',
                                                                     'StandDev:',
                                                                     ' ',
                                                                     ' '), log=JK_results_file_open)
        print_terminal_and_log('{0:<9s} {1:^ 27.4f} {2:^ 30.4f} {3:^31s} {4:^15s}'.format(str(JK_nb),
                                                                           std_rmsd_total_fitted_list,
                                                                           std_rmsd_residues_list,
                                                                           ' ',
                                                                           ' '), log=JK_results_file_open)
    #Create the difference plot
        if os.path.isfile(total_model) and len(list_of_models)>0: #if the total model exists and there is at least one other JK model:
            calculate_difference_plot(total_model, list_of_models, use_waters, i, occ, maptype, refinement) #the values for the plot are calculated and the plot is created
            print_terminal_and_log(
                'see the plot : %s: occ=%s, maptype=%s, refinement=%s - Differences of mean coordinates of residues between the JK models and the total model' % (
                i, occ, maptype, refinement), log=JK_results_file_open) #under the table with CC and RMSD, print the corresponding plot
        else:
            print_terminal_and_log('No models no plot')

        #TODO:remove the temp_dir?         os.remove(outdir + '/temp_dir')

#functions to get CC of maps:
def average_maps(list_of_maps, mapoutname, labels):
    '''
    Get the average map of a list of maps

    Args:
        list_of_maps:
            list of maps (mtz files) to do the average of
        mapoutname:
            str name of the average map output (mtz file)
        labels:
            str,str name of the labels of the columns of the mtz files to get the average of

    Returns:
        avg_map:
            str average map name
    '''
    #get the list of maps in one str with spaces
    str_of_maps=''
    for mp in list_of_maps:
        str_of_maps = str_of_maps+ ' ' + mp

    avg_map = mapoutname + '_average-map.mtz' #name of the average map which will be created

    run_in_terminal('phenix.average_map_coeffs labels=%s %s map_coeffs_out=%s' % (labels, str_of_maps, avg_map), existing_files=[avg_map]) #run phenix.average_map_coeffs in the terminal

    return(avg_map)

def get_cc_mtz_mtz(mtz_avg_map, mtz, outdir):
    '''
    Get the correlation coefficient between the mtz file with all the images and the average mtz file from JK

    Args:
        mtz_avg_map:
            average mtz file
        mtz:
            mtz file to be compared
        outdir:
            output directory of the results (log file) and offset mtz file

    Returns:
        cc:
            float, Final Correlation Coefficient of maps given by get_cc_mtz_mtz

    '''

    mtz_outname = Fextr.get_name(mtz_avg_map) + '_offset.mtz' #offset (reset origin) version of mtz_avg_map to match mtz
    out_log_name=outdir+'/'+mtz_outname+'.log'

    # run phenix.get_cc_mtz_mtz in the terminal
    run_in_terminal('phenix.get_cc_mtz_mtz %s %s offset_mtz=%s output_dir=%s' %(mtz, mtz_avg_map, mtz_outname, outdir), existing_files=[mtz_outname])
    #change the name of the output log file
    run_in_terminal('mv %s %s'%(outdir+'/offset.log', out_log_name), wait=False, existing_files=[out_log_name])

    # get CC from the log file created by get_cc_mtz_mtz containing the final CC
    cc = -1 #In case of Error and CC does not exist, cc=-1
    cc_log_file_open = open(out_log_name, 'r')
    cc_log_file = cc_log_file_open.read().splitlines()  # get list of lines of the log file
    for line in cc_log_file:
        if ('Final CC of maps' in line): #if the line contains the final CC:
            cc = float(line[18:]) #get the CC as number

    return (cc)

#functions to get RMSD of models (some are used for the plot too):
def superpose_pdbs(total_pdb, other_pdb):
    """
    Get the total rmsd between the total model and an other model (fitted to the total model) and get the fitted other model
    Args:
        total_pdb:
            str, total model file
        other_pdb:
            str, other model file
    Returns:
        rmsd_total_fitted:
            float, rmsd between total model and fitted other model for all atomes
        other_pdb_fitted:
            str, other pdb file fitted to total model

    """
    #run phenix.superpose_pdbs in terminal
    _,cmd_output=run_in_terminal('phenix.superpose_pdbs %s %s'%(total_pdb, other_pdb))

    #get the total rmsd from the command output by taking the 5 values after 'RMSD (all matching atoms) (final):'
    for i in range(0, len(cmd_output)):
        if cmd_output[i:i+34] == 'RMSD (all matching atoms) (final):':
            rmsd_total_fitted = float(cmd_output[i+35:i+40]) #get the total rmsd between the models as a float

    other_pdb_fitted = os.getcwd() + '/' + get_name(other_pdb) + '.pdb_fitted.pdb' #name of the other_model file fitted to the total model

    return (rmsd_total_fitted, other_pdb_fitted)

def get_residlist_from_all_atoms(pdb_file, use_waters = False):
    """
    Author: EDZ
    Args:
        pdb_file:
    Returns:
    """
    pdb_hier = hierarchy.input(file_name=pdb_file)
    hier = pdb_hier.hierarchy

    residlst_all = []
    for chain in hier.chains():
        ch = chain.id
        for res_group in chain.residue_groups():
            resv = res_group.resseq
            for atom_group in res_group.atom_groups():
                resn = atom_group.resname
                if (use_waters == False and resn == 'HOH'): continue
                alt = atom_group.altloc
                line = "%4s %4s %4s %4s \n" % (tuple([resn, resv, ch, alt]))
                if line not in residlst_all:
                    residlst_all.append(line)
    residlst = np.array(residlst_all)

    return (residlst)

def get_residlist(resids_lst, use_waters = False):
    """
    Author:EDZ, adapted
    Args:
        resids_lst:
    Returns:
    """
    if resids_lst != None:
        with open(resids_lst) as rs_lst:
            residlst = rs_lst.readlines()
            residlst = [lne for lne in residlst if len(lne) > 0]
        if len(residlst) == 0:
            resids_lst = None
            return
        residlst = np.array(residlst)
        residlst_unique = list(set(residlst))
        residlst_unique.remove([lne for lne in residlst_unique if 'Resn Resv Chain Alt' in lne][0])
        if use_waters == False:
            residlst_unique = [lne for lne in residlst_unique if 'HOH' not in lne]
        residlst = np.array(residlst_unique)

        return(residlst)

def get_coord_from_parser_selected(pdb_file, residlst):
    """
    Author:EDZ
    Get the coordinates of the atoms in the pdb file
    """

    pdb_hier = hierarchy.input(file_name=pdb_file)
    hier = pdb_hier.hierarchy

    chains_to_check = list(set([resid.split()[2] for resid in residlst]))

    coord = []
    info = []
    for chain_to_check in chains_to_check:
        resids = [(resid.split()[0], resid.split()[1]) for resid in residlst if resid.split()[2] == chain_to_check]
        for chain in hier.chains():
            if chain.id == chain_to_check:
                for res_group in chain.residue_groups():
                    for atom_group in res_group.atom_groups():
                        if (atom_group.resname, res_group.resseq.lstrip().rstrip()) in resids:
                            # print (atom_group.resname, res_group.resseq.lstrip().rstrip())
                            for a in res_group.atoms():
                                # print a.name
                                coord.append(list(a.xyz))
                                i = a.fetch_labels()
                                # info.append((i.resname, i.resseq, i.chain_id, i.altloc, i.name, i.i_seq))
                                info.append([i.resname, i.resseq, i.chain_id, i.altloc, i.name]) #resseqq= num AA

    coord = np.asarray(coord)
    info = np.asarray(info)

    return coord, info

def select_common_coords(coord_total, info_total, coord_other, info_other):
    """
    Author:EDZ, adapted
    Compare the info with the common atom info and only retain the coordinates of the common atoms
    """
    # join the info to get a 1d array
    info_ar = np.array([",".join(i) for i in info_other])
    infototal_ar = np.array([",".join(i) for i in info_total])
    # get the the indices of the common atoms
    _, _, indices_retain = np.intersect1d(infototal_ar, info_ar, return_indices=True)
    _, _, indices_retain_total = np.intersect1d(info_ar, infototal_ar, return_indices=True)

    if indices_retain.shape[0] < info_other.shape[0] or indices_retain.shape[0] < info_total.shape[0]:
        print("Trimming atom selection to those common between all pdb files")

    # select the coordinates and info from the common atoms
    coord_other_common = coord_other[indices_retain]
    info_other_common = info_other[indices_retain]
    coord_total_common = coord_total[indices_retain_total]
    info_total_common = info_total[indices_retain_total]

    assert info_other.shape[0] == coord_other.shape[0]
    assert info_total.shape[0] == coord_total.shape[0]

    nb_residlst = coord_total.shape[0] #get the number of atoms in the residue list

    return coord_total_common, info_total_common, coord_other_common, info_other_common, nb_residlst

def calculate_rmsd(coord_total_common, coord_other_common):
    """
    Calculate the rmsd of two sets of atoms from the coordinates of the atoms
    Args:
        coord_total_common:
            list of coord of atoms, dim (3,X)
        coord_other_common:
            list of coord of atoms, dim (3,X)
        axis:

    Returns:
        rmsd:
            float, Root-mean-square deviation between the two sets of atoms

    """
    a= np.array(coord_total_common)
    b= np.array(coord_other_common)

    rmsd = np.sqrt(np.average(np.sum((a - b) ** 2, axis=1)))

    return (rmsd)


#functions to create the plot:
def calculate_difference_between_atoms(coord_listA, coord_listB, axis = 1):
    """
    Calculate distances between two list of atoms, gives the list of distances between each atoms
    Args:
        coord_listA: list of coordinates dim(X, 3)
        coord_listB: list of coordinates dim(X, 3)
        axis: nb of axis of variation of atoms from the lists

    Returns:
        distance_list: list of distances between each atoms
    """
    #get arrays
    A = np.array(coord_listA)
    B = np.array(coord_listB)

    #calculat the list of differences between the atoms of A and B
    distances_list = np.sqrt(np.sum((A - B) ** 2, axis=axis))

    return (distances_list)

def calculate_difference_per_coord(coord_listA, coord_listB):
    """
    Calculate distances between two coordinates of list of coordinates of atoms, gives the list of distances between the coordinate of each atoms
    Args:
        coord_listA: list of coordinates dim(X, 3)
        coord_listB: list of coordinates dim(X, 3)

    Returns:

    """
    A = np.array(coord_listA)
    B = np.array(coord_listB)

    # calculat the list of differences between the coord of A and B for x, y and z
    X_difference_list = (A[:,0] - B[:,0])
    Y_difference_list = (A[:,1] - B[:,1])
    Z_difference_list = (A[:,2] - B[:,2])

    return (X_difference_list, Y_difference_list, Z_difference_list)

def calculate_mean_point_of_residue(coord_list_residue):
    """
    Calculate the coordinates of the mean point of the atoms of a residue
    Args:
        coord_list_residue: list of coordinates of the atoms of the residue dim(X, 3)

    Returns:
        coord_mean_point: list of coordinates of the mean point of the residue dim(1, 3)
    """
    #not possible because the multiple conformations => more than 32 atoms to get the average of, which is not possible with np.ndarray
    # coord_list_residue= np.ndarray(coord_list_residue)
    # coord_mean_point=np.ndarray.mean(coord_list_residue,axis=0)

    coord_list_residue=np.array(coord_list_residue)
    coord_mean_point = [np.mean(tuple(coord_list_residue[:,0])),
                                  np.mean(tuple(coord_list_residue[:,1])),
                                  np.mean(tuple(coord_list_residue[:,2]))]

    return(coord_mean_point)

def cartesian_to_spherical(x, y, z):
    """
    Calculate the spherical coordinates from the cartesian coordinates
    Args:
        x, y, z: 3 lists of the cartesian coordinates

    Returns:
        r, theta, phi : 3 lists of the spherical coordinates
    """
    xy = x**2 + y**2
    r = np.sqrt(xy + z**2)
    theta = np.arctan2(np.sqrt(xy), z)
    phi = np.arctan2(y, x)

    return r, theta, phi


#main function that creates the plot:
def calculate_difference_plot(pdb_file_total, pdb_file_list, use_waters, i, occ, maptype, refinement):
    '''
    Calculate the differences between the JK models and the total model to create the final plot of the 'Differences of mean coordinates of residues between the JK models and the total model'
    Args:
        pdb_file_total: pdb file, created with all the images
        pdb_file_list: list of pdb files from JK for the specific occ, maptype, refinement
        use_waters: bool
        i: int, line of the tab_list (correspond to a specific occ, maptype, refinement)
        occ: float, occupancy
        maptype: str
        refinement: str

    Returns:
        Plot of 'Differences of mean coordinates of residues between the JK models and the total model'
    '''

#1. get the information from the total pdb file
    #get coordinates and info of all atoms of the total pdb file
    residlist_all_total = get_residlist_from_all_atoms(pdb_file_total, use_waters=use_waters)
    coord_all_total, info_all_total = get_coord_from_parser_selected(pdb_file_total, residlist_all_total)
    coord_all_total=np.array(coord_all_total)
    info_all_total=np.array(info_all_total)

    #get chains and number of chains from total pdb
    chain_id_list = list(set(info_all_total[:, 2])) #list of chain IDs of the total pdp file (=for each pdb file)
    nb_chains = len(chain_id_list) #number of chains in the pdb files

    def get_list_chainid_residue_begining_end_indice(info_all_total):
        #init list [chain_id, i_beginning_resid, i_end_resid]
        list_chainid_begeningresid_endresid = []

        i_resid_nb = 0
        i_begining_resid = 0
        for i in range(1, len(info_all_total[:, 1])): #for each line of the list of coord
            if info_all_total[i, 1] == info_all_total[i-1, 1]: #if the resid number is the same as the previous one
            #same residue
            #TODO: I do not know if we need to check if the atom is not water or part of a ligand without being in the chain S
                #if the atom is water ???
                #if the atom is part of a ligand???
                True

            else: #next residue
                #get the number of the next residue
                i_resid_nb += 1
                i_end_resid = i-1

                #add chain id, indices begining and end resid to list
                list_chainid_begeningresid_endresid.append([info_all_total[i-1, 2], i_begining_resid, i_end_resid])

                i_begining_resid=i #begining indice of the new resid

        return(list_chainid_begeningresid_endresid)

    #get the list of [chain_id, i_beginning_resid, i_end_resid] from the total pdb file
    list_chainid_begeningresid_endresid = get_list_chainid_residue_begining_end_indice(info_all_total)
    list_chainid_begeningresid_endresid = np.array(list_chainid_begeningresid_endresid)

#2. get the list [chain_id, resid_nb, coord_mean_point_of_resid] of the total
    def get_list_chainid_residnb_meanresidcoord(coord_all_total, info_all_total, list_chainid_begeningresid_endresid):

        #init list [chain_id, resid_nb, coord_mean_point_of_resid]
        list_chainid_residnb_meanresidcoord = []
        for i in range (0, len(list_chainid_begeningresid_endresid[:, 1])):
            chain_id = list_chainid_begeningresid_endresid[i, 0]
            i_begining_resid=list_chainid_begeningresid_endresid[i, 1]
            i_end_resid=list_chainid_begeningresid_endresid[i, 2]

            # get list of coordinates of atoms of resid
            coord_atoms_resid=[]

            for ind in range(int(i_begining_resid), int(i_end_resid)+1):
                coord_atoms_resid.append(coord_all_total[ind, :])

            # get the mean point of the resid
            coord_mean_point_of_resid = calculate_mean_point_of_residue(coord_atoms_resid)

            #resid number
            resid_nb = info_all_total[int(i_begining_resid)][1]

            list_chainid_residnb_meanresidcoord.append([chain_id, resid_nb, coord_mean_point_of_resid])

        return(list_chainid_residnb_meanresidcoord)

    list_chainid_residnb_meanresidcoord_total = get_list_chainid_residnb_meanresidcoord(coord_all_total, info_all_total, list_chainid_begeningresid_endresid) #get the list of [chain_id, resid_nb, coord_mean_point_of_resid] for the total pdb file
    list_chainid_residnb_meanresidcoord_total = np.array(list_chainid_residnb_meanresidcoord_total)


#Create initial plot
    plt.close()
    fig, axs = plt.subplots(4, nb_chains, figsize=(nb_chains*12, 3*5), squeeze=False)
#    figatm, axsatm = plt.subplots(3, nb_chains, figsize=(3*7, nb_chains*7), squeeze=False)

#3. for each JK pdb file, get the list [chain_id, resid_nb, coord_mean_point_of_resid]
    for pdb_file in pdb_file_list: #for each pdb file
        # get coordinates of all atoms of the pdb file
        residlist_all = get_residlist_from_all_atoms(pdb_file, use_waters=use_waters)
        coord_all, info_all = get_coord_from_parser_selected(pdb_file, residlist_all)
        coord_all=np.array(coord_all)
        list_chainid_residnb_meanresidcoord = get_list_chainid_residnb_meanresidcoord(coord_all, info_all_total, list_chainid_begeningresid_endresid) #get the list of [chain_id, resid_nb, coord_mean_point_of_resid] for the pdb file
        list_chainid_residnb_meanresidcoord = np.array(list_chainid_residnb_meanresidcoord)

#4. Separating the chains to plot and getting the coordinates to plot:
#X axis = AA_axis - list of residue numbers
#Y axis = atom_distances, X_distances, Y_distances, Z_distances - list of distances of atoms, Xcoord, Ycoord, Zcoord between the JK model and the total model
        col=0 #first column of the plot

    #Separating the chains
        for chain_id in chain_id_list: #for each chain_id in the list of chain_id for the pdb files
            #init for each chain
            AA_axis=[] #residue list = X axis
            mean_coord=[] #list of the mean coordinates of each residue for the pdb file
            mean_coord_total=[] #list of the mean coordinates of each residue for the total pdb file
            for indice in range(0, len(list_chainid_residnb_meanresidcoord[:,0])): #for each indice of the list of [chain_id, resid_nb, coord_mean_point_of_resid]
                if chain_id == list_chainid_residnb_meanresidcoord[indice,0]: #if it is the same chain
                    AA_axis.append(list_chainid_residnb_meanresidcoord[indice,1]) #add residue nb to list of residue (X axis)
                    mean_coord_total.append(list_chainid_residnb_meanresidcoord_total[indice,2]) #add coordinates of the total mean residue to the list of total mean coordinates for each residue
                    mean_coord.append(list_chainid_residnb_meanresidcoord[indice,2]) #add coordinates of the mean residue to the list of mean coordinates for each residue

            #calculating the differences between JK model and total model
#            meanresidue_distances = calculate_difference_between_atoms(mean_coord_total, mean_coord) #getting the distances between the mean of the residues of the two models (=r_distances)
            X_distances, Y_distances, Z_distances = calculate_difference_per_coord(mean_coord_total, mean_coord)
            r_distances, theta_distances, phi_distances = cartesian_to_spherical(X_distances, Y_distances, Z_distances)
            rxthetaxphi_distances= r_distances*theta_distances*phi_distances

            #Alternative: calculating the differences between JK model and total model for each atom
#            atom_distances = calculate_difference_between_atoms(coord_all, coord_all_total)
#            X_distances_atom, Y_distances_atom, Z_distances_atom = calculate_difference_per_coord(coord_all, coord_all_total)
#            r_distances_atom, theta_distances_atom, phi_distances_atom = cartesian_to_spherical(X_distances_atom, Y_distances_atom, Z_distances_atom)


#5. Ploting #TODO: Also plotting the B factors somehow?
            #plotting distances between mean point of residues of JK model and total model
            ax1= axs[0, col].scatter(AA_axis, r_distances, marker='.', label='JK')
            axs[1, col].scatter(AA_axis, theta_distances, marker='.', label='JK')
            axs[2, col].scatter(AA_axis, phi_distances, marker='.', label='JK')
            axs[3, col].scatter(AA_axis, rxthetaxphi_distances, marker='.', label='JK')

        # Alternative: Ploting differences between every atoms for each residue
            # plotting distances between mean point of residues of JK model and total model
            #            axs[0, col].scatter(info_all[:,1], r_distances_atom, marker='.', label='JK')
            #            axs[1, col].scatter(info_all[:,1], theta_distances_atom, marker='.', label='JK')
            #            axs[2, col].scatter(info_all[:,1], phi_distances_atom,  marker='.', label='JK')

        # trying to plot the atoms gradually per color: did not work
#            no_points = len(r_distances)
#            axs[0, col].set_color_cycle([ax1(1. * i / (no_points - 1)) for i in range(0,no_points - 1)])

            #additions to plot for every plot
            for row in range(0,4):
                #putting x ticks every factor 10 and adding minor ticks at every factor 5
                axs[row, col].xaxis.set_major_locator(MultipleLocator(10))
                axs[row, col].xaxis.set_minor_locator(MultipleLocator(5))
                #adding line at y = 0, corresponding to the model
                axs[row, col].axhline(0, color='grey', linewidth=0.8)
                #adding legend
                axs[row, col].legend()

            #adding y labels
            axs[0, col].set_ylabel('r (A)')
            axs[1, col].set_ylabel('theta (radian)')
            axs[2, col].set_ylabel('phi (radian)')
            axs[3, col].set_ylabel('r * theta * phi (A*radian2)')

            #adding x label and chain
            axs[3, col].set_xlabel('Residues')
            axs[0, col].set_title('Chain %s' % chain_id, fontsize='medium', fontweight="bold")
            col+=1

    #global title of the plot
    fig.suptitle('%s: occ=%s, maptype=%s, refinement=%s - Differences of mean coordinates of residues between the JK models and the total model'%(i, occ, maptype, refinement))

    plt.subplots_adjust(hspace=0.35, left=0.09, right=0.65, top=0.95)
    plt.savefig('%s_JK_Difference_plot.pdf' %(i), dpi=300, transparent=False)
    plt.savefig('%s_JK_Difference_plot.png' % (i), dpi=300)

    plt.show()
    plt.close()


#for test
# total_tab=[['occupancy', 'map_type', 'Refinement', 'map', 'model', 'labels', 'total',
#   'residlist'],
#  [1, 'qFoFo', None,
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/Cyt1A_mqFoFo.mtz',
#   None, 'QFOFOWT,PHIQFOFOWT', 'total',
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/residlist_Zscore2.00.txt'],
#  [0.1, 'qFextr_map', None,
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc.mtz',
#   None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 'total', None],
#  [0.1, 'qFextr_map', 'reciprocal',
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.mtz',
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#   '2FOFCWT,PH2FOFCWT', 'total', None],
#  [0.1, 'qFextr_map', False, None,
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.100//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#   None, 'total', None],
#  [0.1, 'qFextr_map', 'reciprocal + real', None,
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#   None, 'total', None],
#  [0.30000000000000004, 'qFextr_map', None,
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc.mtz',
#   None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 'total', None],
#  [0.30000000000000004, 'qFextr_map', 'reciprocal',
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.mtz',
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#   '2FOFCWT,PH2FOFCWT', 'total', None],
#  [0.30000000000000004, 'qFextr_map', 'real' ,None,
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#   None, 'total', None],
#  [0.30000000000000004, 'qFextr_map', False, None,
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#   None, 'total', None],
#  [0.5, 'qFextr_map', None,
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc.mtz',
#   None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 'total', None],
#  [0.5, 'qFextr_map', 'reciprocal',
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.mtz',
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#   '2FOFCWT,PH2FOFCWT', 'total', None],
#  [0.5, 'qFextr_map', 'real', None,
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#   None, 'total', None],
#  [0.5, 'qFextr_map', False, None,
#   '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_total/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#   None, 'total', None]]
# tab_list=[np.array([['occupancy', 'map_type', 'Refinement', 'map', 'model', 'labels',
#         'total', 'residlist'],
#        [1, 'qFoFo', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/Cyt1A_mqFoFo.mtz',
#         None, 'QFOFOWT,PHIQFOFOWT', 0,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/residlist_Zscore2.00.txt'],
#        [0.1, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 0,
#         None],
#        [0.1, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 0, None],
#        [0.1, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 0, None],
#        [0.1, 'qFextr_map', 'reciprocal + real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#         None, 0, None],
#        [0.30000000000000004, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 0,
#         None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 0, None],
#        [0.30000000000000004, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.300//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 0, None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal + real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#         None, 0, None],
#        [0.5, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 0,
#         None],
#        [0.5, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 0, None],
#        [0.5, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 0, None],
#        [0.5, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_0/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 0, None]], dtype=object), np.array([['occupancy', 'map_type', 'Refinement', 'map', 'model', 'labels',
#         'total', 'residlist'],
#        [1, 'qFoFo', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/Cyt1A_mqFoFo.mtz',
#         None, 'QFOFOWT,PHIQFOFOWT', 1,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/residlist_Zscore2.00.txt'],
#        [0.1, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 1,
#         None],
#        [0.1, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 1, None],
#        [0.1, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 1, None],
#        [0.1, 'qFextr_map', 'reciprocal + real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#         None, 1, None],
#        [0.30000000000000004, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 1,
#         None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 1, None],
#        [0.30000000000000004, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 1, None],
#        [0.30000000000000004, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 1, None],
#        [0.5, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 1,
#         None],
#        [0.5, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 1, None],
#        [0.5, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 1, None],
#        [0.5, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_1/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 1, None]], dtype=object), np.array([['occupancy', 'map_type', 'Refinement', 'map', 'model', 'labels',
#         'total', 'residlist'],
#        [1, 'qFoFo', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/Cyt1A_mqFoFo.mtz',
#         None, 'QFOFOWT,PHIQFOFOWT', 2,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/residlist_Zscore2.00.txt'],
#        [0.1, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 2,
#         None],
#        [0.1, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 2, None],
#        [0.1, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.100//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 2, None],
#        [0.1, 'qFextr_map', 'reciprocal + real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#         None, 2, None],
#        [0.30000000000000004, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 2,
#         None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 2, None],
#        [0.30000000000000004, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.300//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 2, None],
#        [0.30000000000000004, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 2, None],
#        [0.5, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 2,
#         None],
#        [0.5, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 2, None],
#        [0.5, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.500//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 2, None],
#        [0.5, 'qFextr_map', 'reciprocal + real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_2/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#         None, 2, None]], dtype=object), np.array([['occupancy', 'map_type', 'Refinement', 'map', 'model', 'labels',
#         'total', 'residlist'],
#        [1, 'qFoFo', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/Cyt1A_mqFoFo.mtz',
#         None, 'QFOFOWT,PHIQFOFOWT', 3,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/residlist_Zscore2.00.txt'],
#        [0.1, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 3,
#         None],
#        [0.1, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 3, None],
#        [0.1, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.100//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 3, None],
#        [0.1, 'qFextr_map', 'reciprocal + real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#         None, 3, None],
#        [0.30000000000000004, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 3,
#         None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 3, None],
#        [0.30000000000000004, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.300//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 3, None],
#        [0.30000000000000004, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 3, None],
#        [0.5, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 3,
#         None],
#        [0.5, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 3, None],
#        [0.5, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 3, None],
#        [0.5, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_3/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 3, None]], dtype=object), np.array([['occupancy', 'map_type', 'Refinement', 'map', 'model', 'labels',
#         'total', 'residlist'],
#        [1, 'qFoFo', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/Cyt1A_mqFoFo.mtz',
#         None, 'QFOFOWT,PHIQFOFOWT', 4,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/residlist_Zscore2.00.txt'],
#        [0.1, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 4,
#         None],
#        [0.1, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 4, None],
#        [0.1, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 4, None],
#        [0.1, 'qFextr_map', 'reciprocal + real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#         None, 4, None],
#        [0.30000000000000004, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 4,
#         None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 4, None],
#        [0.30000000000000004, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.300//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 4, None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal + real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#         None, 4, None],
#        [0.5, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 4,
#         None],
#        [0.5, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 4, None],
#        [0.5, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.500//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 4, None],
#        [0.5, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_4/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 4, None]], dtype=object), np.array([['occupancy', 'map_type', 'Refinement', 'map', 'model', 'labels',
#         'total', 'residlist'],
#        [1, 'qFoFo', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/Cyt1A_mqFoFo.mtz',
#         None, 'QFOFOWT,PHIQFOFOWT', 5,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/residlist_Zscore2.00.txt'],
#        [0.1, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 5,
#         None],
#        [0.1, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 5, None],
#        [0.1, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 5, None],
#        [0.1, 'qFextr_map', 'reciprocal + real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#         None, 5, None],
#        [0.30000000000000004, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 5,
#         None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 5, None],
#        [0.30000000000000004, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.300//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 5, None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal + real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#         None, 5, None],
#        [0.5, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 5,
#         None],
#        [0.5, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 5, None],
#        [0.5, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 5, None],
#        [0.5, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_5/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 5, None]], dtype=object), np.array([['occupancy', 'map_type', 'Refinement', 'map', 'model', 'labels',
#         'total', 'residlist'],
#        [1, 'qFoFo', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/Cyt1A_mqFoFo.mtz',
#         None, 'QFOFOWT,PHIQFOFOWT', 6,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/residlist_Zscore2.00.txt'],
#        [0.1, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 6,
#         None],
#        [0.1, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 6, None],
#        [0.1, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 6, None],
#        [0.1, 'qFextr_map', 'reciprocal + real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#         None, 6, None],
#        [0.30000000000000004, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 6,
#         None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 6, None],
#        [0.30000000000000004, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.300//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 6, None],
#        [0.30000000000000004, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 6, None],
#        [0.5, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 6,
#         None],
#        [0.5, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 6, None],
#        [0.5, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.500//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 6, None],
#        [0.5, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_6/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 6, None]], dtype=object), np.array([['occupancy', 'map_type', 'Refinement', 'map', 'model', 'labels',
#         'total', 'residlist'],
#        [1, 'qFoFo', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/Cyt1A_mqFoFo.mtz',
#         None, 'QFOFOWT,PHIQFOFOWT', 7,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/residlist_Zscore2.00.txt'],
#        [0.1, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 7,
#         None],
#        [0.1, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 7, None],
#        [0.1, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.100//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 7, None],
#        [0.1, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 7, None],
#        [0.30000000000000004, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 7,
#         None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 7, None],
#        [0.30000000000000004, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.300//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 7, None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal + real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#         None, 7, None],
#        [0.5, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 7,
#         None],
#        [0.5, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 7, None],
#        [0.5, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.500//mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/data_X8/data/6t14_modifiedChainID.pdb',
#         None, 7, None],
#        [0.5, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_7/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 7, None]], dtype=object), np.array([['occupancy', 'map_type', 'Refinement', 'map', 'model', 'labels',
#         'total', 'residlist'],
#        [1, 'qFoFo', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/Cyt1A_mqFoFo.mtz',
#         None, 'QFOFOWT,PHIQFOFOWT', 8,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/residlist_Zscore2.00.txt'],
#        [0.1, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 8,
#         None],
#        [0.1, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 8, None],
#        [0.1, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 8, None],
#        [0.1, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 8, None],
#        [0.30000000000000004, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 8,
#         None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 8, None],
#        [0.30000000000000004, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 8, None],
#        [0.30000000000000004, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 8, None],
#        [0.5, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 8,
#         None],
#        [0.5, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 8, None],
#        [0.5, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 8, None],
#        [0.5, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_8/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 8, None]], dtype=object), np.array([['occupancy', 'map_type', 'Refinement', 'map', 'model', 'labels',
#         'total', 'residlist'],
#        [1, 'qFoFo', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/Cyt1A_mqFoFo.mtz',
#         None, 'QFOFOWT,PHIQFOFOWT', 9,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/residlist_Zscore2.00.txt'],
#        [0.1, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 9,
#         None],
#        [0.1, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 9, None],
#        [0.1, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 9, None],
#        [0.1, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.100/Cyt1A_occ0.100_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 9, None],
#        [0.30000000000000004, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 9,
#         None],
#        [0.30000000000000004, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 9, None],
#        [0.30000000000000004, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 9, None],
#        [0.30000000000000004, 'qFextr_map', False, None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.300/Cyt1A_occ0.300_2mqFextr-DFc_reciprocal_space_001.pdb',
#         None, 9, None],
#        [0.5, 'qFextr_map', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc.mtz',
#         None, '2QFEXTRFCWT,PHI2QFEXTRFCWT QFEXTRFCWT,PHIQFEXTRFCWT', 9,
#         None],
#        [0.5, 'qFextr_map', 'reciprocal',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.mtz',
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001.pdb',
#         '2FOFCWT,PH2FOFCWT', 9, None],
#        [0.5, 'qFextr_map', 'real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_mqFextr-DFc_real_space_refined_000.pdb',
#         None, 9, None],
#        [0.5, 'qFextr_map', 'reciprocal + real', None,
#         '/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_CYT1A_WT_PH10_ALL_triggered/CYT1A_WT_PH10_ALL_JK_90percent_9/Xtrapol8/qweight_occupancy_0.500/Cyt1A_occ0.500_2mqFextr-DFc_reciprocal_space_001_real_space_refined_000.pdb',
#         None, 9, None]], dtype=object)]
# outdir='/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/Xtrapol8_test_1_Cyt1A/Cyt1A_JackKnife_Xtrapol8_20/JK_average_and_comparison_results'
#
# ordered_solvent=True
# #
# get_JK_results(total_tab, tab_list, outdir, ordered_solvent)