# -*- coding: utf-8 -*-
"""
Script to calculate the occupancy based on the real space refined models of a reference model in Extrapolated map
coefficients calculated with different occupancy (alpha).

Automatically run in Xtrapol8 routine but can be run on a standalone basis.

What do you need?
- a complete Xtrapol8 run in Calm_and_curious or Fast_and_furious mode
- an ESFA type (this should have been run in the Xtrapol8 run that precedes the differencemap_analysis
- Optional:
    - Residue list for which the difference map peaks will be used to estimate the occupancy (in the format as an Xtrapol8 residuelist)
        If no residue list provided, a residue list will be calculated based on the Z-score
    - suffix to be added to the output files (will be used as a prefix for some files)
    - Log-file
    - Output directory
    
usage
-----
To get the help message:
$ phenix.python distance_analysis.py
or
$ phenix.python distance_analysis.py --help
or
$ phenix.python distance_analysis.py -h

To run with an input file:
$ phenix.python distance_analysis.py distance_analysis.phil

To run with command line arguments:
$ phenix.python distance_analysis.py input.Xtrapol8=Xtrapol8/Xtrapol8_out.phil input.f_extrapolated_and_maps=qfextr
 map_explorer.residue_list=residlist_adapted.txt output.outdir=distance_test

To run with an input file and command line arguments:
$ phenix.python distance_analysis.py distance_analysis.phil input.f_extrapolated_and_maps=fextr_calc
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
TODO: find elegant alternative for global variables

"""

from __future__ import division, print_function
import sys, re, os
import string
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from iotbx.pdb import hierarchy
from scipy.optimize import curve_fit
import scipy.stats
import glob
from Fextr_utils import check_file_existance

colorlib=['purple','indigo','rebeccapurple', 'midnightblue', 'darkblue', 'mediumblue', 'blue', 'royalblue', 'dodgerblue', 'cornflowerblue', 'deepskyblue', 'lightskyblue', 'cadetblue','darkcyan', 'darkturquoise', 'mediumturquoise', 'turquoise', 'aqua', 'mediumaquamarine', 'aquamarine', 'mediumspringgreen', 'springgreen', 'green', 'lime', 'lawngreen', 'chartreuse', 'greenyellow', 'yellow', 'gold', 'goldenrod', 'orange', 'darkorange', 'chocolate', 'darksalmon', 'orangered', 'red', 'firebrick', 'maroon', 'darkred', 'black']

def sigmoid_fit(x, L, k, x0):
    """
    logistic function = sigmoidal
    L = maximum
    k = steepness
    x0 = sigmoids midpoint
    add additional term of +1 so that the function starts at alpha of 1
    """
    #print("L, k, x0", L, k, x0)
    return L / (1 + np.exp(-k*(x-x0)-1))

def sigmoid_fit_2(fact,x):
    """
    variant of previous function for easier use with np.apply_along_axis
    """
    L, k, x0 = fact
    return L / (1 + np.exp(-k*(x-x0)-1))
    
    
def logfit(x, a, b, c):
    """
    logistic function = logfit
    """
    #c = 0
    return a * (1- np.exp(-b * (x-1))) + c * x

def logfit_2(fact,x):
    """
    variant of previous function for easier use with np.apply_along_axis
    """
    a, b, c = fact
    #c = 0
    return a * (1- np.exp(-b * (x-1))) + c * x

def sumofsquares(arr):
    return np.sum(arr**2)

def totalsumsquares(arr):
    return np.sum((arr-np.mean(arr))**2)

def fitting(sel, x):
    """
    fit with logfit model. Not used anymore.
    """
    if np.all(sel ==0):
        return np.array([0.,0.,0.])
    else:
        try:
            #Extend the bmax-bkmin range for 0.05 on each side
            bmax = val +0.05 #see below for val definition and why bmax=val
            bmin = val/(maxalpha-1) -0.05#see below for val definition and why bmin=val/100
            if bmax < bmin:
                return np.array([0.,0.,0.])
            
            popt,_ = curve_fit(logfit, x,sel,sigma=np.sqrt(x),absolute_sigma=False,
                                bounds=((0.95*np.max(sel), bmin, -0.001), (1.05*np.max(sel), bmax, 0.001)))
        except ValueError: 
            popt = np.array([0.,0.,0.])
        #print(np.array([popt[0], popt[1], popt[2]]))
        
        #Excplude based on alpha. 
        alpha = (val/popt[1])+1
        if alpha < 1.0:
            #print("Minimum alpha reached. alpha and occupancy:", alpha, 1/alpha)
            return np.array([0.,0.,0.])
        if alpha > maxalpha:
            #print("Maximum alpha reached. alpha and occupancy:", alpha, 1/alpha)
            return np.array([0.,0.,0.])

        return np.array(popt)
    
def sigmoid_fitting(sel, x):
    """
    fit via the sigmoidal model instead of the logfit
    p0 not defined
    Bounds: L: 0.80*np.max(sel), 1.15*np.max(sel)
            k x (alpha - x0) = val
            alpha has to lie between 1 and 100, and thus also the x0 value has to lie between 1 and 100
            k: val - val/100 or val - val/(maximum of probed alphas)
            x0: 1 - 100 or 1 - maximum of probed alphas
    """
    if np.all(sel ==0):
        return np.array([0.,0.,0.])
    else:
        try:
            kmax = val/1 +0.05 #Extend the kmax-kmin range for 0.05 on each side
            kmin = val/maxalpha - 0.05
            if kmax < kmin:
                return np.array([0.,0.,0.])
                #kmin = val
                #kmax = val/maxalpha
            x0max =maxalpha
            Lmin = 0.95*np.max(sel)
            Lmax = 1.05*np.max(sel)
            popt,_ = curve_fit(sigmoid_fit, x, sel, sigma=np.sqrt(x),absolute_sigma=False, bounds = ((Lmin, kmin, 1), (Lmax, kmax, x0max)))
        except ValueError: 
            popt = np.array([0.,0.,0.])
        except RuntimeError:
            try:
                #print("loosen fitting boundaries")
                kmax = val/1 +0.15 #Extend the kmax-kmin range for 0.15 on each side
                kmin = val/maxalpha - 0.15
                if kmax < kmin:
                    return np.array([0.,0.,0.])
                    #kmin = val
                    #kmax = val/maxalpha
                x0max =maxalpha
                Lmin = 0.80*np.max(sel)
                Lmax = 1.20*np.max(sel)
                popt,_ = curve_fit(sigmoid_fit, x, sel, sigma=np.sqrt(x),absolute_sigma=False, bounds = ((Lmin, kmin, 1), (Lmax, kmax, x0max)))
            except ValueError: 
                popt = np.array([0.,0.,0.])
            except RuntimeError:
                popt = np.array([0.,0.,0.])
            
        #Exclude based on alpha. 
        alpha = ((val)/popt[1])+ popt[2]
        if alpha < 1.0:
            #print("Minimum alpha reached. alpha and occupancy:", alpha, 1/alpha)
            return np.array([0.,0.,0.])
        if alpha > maxalpha:
            #print("Maximum alpha reached. alpha and occupancy:", alpha, 1/alpha)
            return np.array([0.,0.,0.])
        #print(np.array([popt[0], popt[1], popt[2]]))
        return np.array(popt)
    
class Distance_analysis(object):
    """
    Class to carry out the distance analysis. Fitting is done with a sigmoidal function. In the past we used an exponential function but the sigmoidal appears to work better.
    Needs a list of pdb files, with the first one being the reference file. 
    The occupancy list should be in the same order as the pdb files and no value should be given for the first pdb file (we'll add an alpha value of 1)
    """
    def __init__(self, pdblst, occupancies, resids_lst = None, plateau_fraction=0.99, use_waters = True, outsuffix = '', log = sys.stdout):
        self.resids_lst       = resids_lst
        if plateau_fraction < 0.0:
            print("Distance analysis plateau fraction should be set between 0.0 and 1.0")
            self.plateau_fraction = 0.99
        elif plateau_fraction > 1.0:
            print("Distance analysis plateau fraction should be set between 0.0 and 1.0")
            self.plateau_fraction = 0.99
        else:
            self.plateau_fraction = plateau_fraction
        self.occupancies      = occupancies
        self.use_waters       = use_waters
        self.outsuffix        = outsuffix
        self.log              = log
        
        print("DISTANCE ANALYSIS")
        print("DISTANCE ANALYSIS", file=self.log)
        alphas = list(map(lambda x: round(1/x, 3), occupancies))
            
        if (1.0 not in alphas and len(alphas) != len(pdblst)):
            alphas = [1]+alphas

        if len(alphas) != len(pdblst):
            print("Different number of models and occupancies. The distance analysis cannot be carried out.")
            sys.exit()
            
        ##exponential fitting:
        ##Search alpha for 0.95 = 1*(1-np.exp(-b*alpha))
        ##   0.05 = 1/(np.exp(b*alpha)
        ##   ln(1/0.05) = b*alpha
        ##   ln(20)/b = alpha
        ##For 98% of plateau level, need ln(50)
        ##If a != 1, we need to use the gerenal case:
        ##find x for 0.95 = a*(1-np.exp(-b*x))
        ##   ln(-a/(0.95-a))/b = x
        ##for any plateau-fraction (np.log is natural logarithm):
        ##Just as done for the sigmoidal fit, a shift to the right has to be done. However, it cannot be brought into the constant val and needs to be taken into account upon calculation
        #val = np.log(1/(1-self.plateau_fraction))
        #global val
        ##val/b = alpha
        ##To find max and min b-values:
        ##   b = val/(alpha-1)
        ##alpha between 1 (but deviding by 1-1=0 is impossible) and maxalpha
        ##bmax = val
        ##bmin = val/(maxalpha-1)
        
        
        #Sigmoidal fitting: plateau_value = 1/ (1+ exp(-k x (alpha - alphainflection))
        #-k x (alpha - alphainflection) = ln[(1/plateau_fraction) -1]
        #k x (alpha - alphainflection) = - ln[(1/plateau_fraction) -1] = val
        #To estimate alpha however, we should take into account that alpha between 0 and 1 is impossible, hence the curve should actually start at 1.
        #   This comes down to shifting the complete sigmoidal to the right. Hence to calcalute alpha, one should use the following formula
        #   k x (alpha - alphainflection) -1 = val
        global val
        val = -1 * np.log((1/self.plateau_fraction)-1) + 1
        
     
        #Sort the pdb files and alphas in order to have alpha from small to large, this is important for fitting
        #This might not work in pyhton3
        zipped = zip(alphas, pdblst)
        zipped_sorted = sorted(zipped, key = lambda x:x[0])
        alphas, pdblst =zip(*zipped_sorted)
        self.alphas = list(alphas)
        self.pdblst = list(pdblst)
        
        global maxalpha
        maxalpha = np.max(self.alphas)
     
        print("Reference pdb file:", self.pdblst[0], file=log)
        print("Reference pdb file:", self.pdblst[0])
        self.get_residlist()
        if self.resids_lst == None:
            print("Using all atoms for distance analysis is computationally heavy and can take a while.")
            self.get_residlist_from_all_atoms(self.pdblst[0])
        print("Distance analysis based on %d residues" %(self.residlst.shape[0]), file=self.log)
        print("Distance analysis based on %d residues" %(self.residlst.shape[0]))
        #print("self.residlst",self.residlst)
        self.initialize_dicts(self.pdblst[0]) #Not used anymore?
        self.check_common_atoms()
                
        
    def get_residlist(self):
        """
        Make a residue list based on an input map-explorer like file
        """
        if self.resids_lst != None:
            with open(self.resids_lst) as rs_lst:
                residlst = rs_lst.read().split("\n")
                residlst = [lne for lne in residlst if len(lne)>0]
            if len(residlst) == 0: 
                self.resids_lst = None
                return
            residlst = np.array(residlst)
            residlst_unique = list(set(residlst))
            residlst_unique.remove([lne for lne in residlst_unique if 'Resn Resv Chain Alt' in lne][0])
            if self.use_waters == False:
                residlst_unique = [lne for lne in residlst_unique if 'HOH' not in lne]
            self.residlst = np.array(residlst_unique)
            
    def get_residlist_from_all_atoms(self, pdb_file):
        """
        Make a residue list based on all atoms in a PDB file
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
                    if (self.use_waters == False and resn == 'HOH'):continue
                    alt = atom_group.altloc
                    line = "%4s %4s %4s %4s \n" %(tuple([resn, resv, ch, alt]))
                    if line not in residlst_all:
                        residlst_all.append(line)
        self.residlst = np.array(residlst_all)
        
    def check_common_atoms(self):
        """
        Get list with the array atom info in order to make sure that we only compare those . Genereally, this should be performed only between the reference (first pdb from the list) and second pdb file as the others should contain the same atoms as the second one, unless use_waters was wrongly set to True. Hence, check all of them
        """
        
        if len(self.pdblst) > 1:
            #Check the common atoms between the first the 2 PDB-files of the list:
            #get info of first pdb
            _, info_0 = self.get_coord_from_parser_selected(self.pdblst[0])
            #join the info to get a 1d array
            info_0_ar = np.array([",".join(i) for i in info_0])
            #get info of second pdb
            _, info_1 = self.get_coord_from_parser_selected(self.pdblst[1])
            #join the info to get a 1d array
            info_1_ar = np.array([",".join(i) for i in info_1])
            #get the common info
            comm = np.intersect1d(info_0_ar, info_1_ar,return_indices=False)
            
            #Check common atoms between the common info and the info of the other pdb files
            for pdb in self.pdblst[2:]:
                #get the info
                _, info = self.get_coord_from_parser_selected(pdb)
                #join the info to get a 1d array
                info_ar = np.array([",".join(i) for i in info])
                #get the common info
                comm = np.intersect1d(comm, info_ar, return_indices=False)
                
        self.common_info = comm

    def initialize_dicts(self, pdb_file):
        """
        From pdb-file, should be the reference model so launch with pdblst[0], extract all residue types and their atoms
        """
        pdb_hier = hierarchy.input(file_name=pdb_file)
        hier = pdb_hier.hierarchy

        residatoms = {}
        for chain in hier.chains():
            #if chain.is_protein():
                for res_group in chain.residue_groups():
                    for atom_group in res_group.atom_groups():
                        if atom_group.resname not in residatoms:
                            atoms = list(set([a.name for a in atom_group.atoms()]))
                            residatoms[atom_group.resname]=atoms
        self.residatoms = residatoms
    
    def get_coord_from_parser_selected(self, pdb_file):
        """
        Get the coordinates of the atoms in the pdb file
        """
        
        pdb_hier = hierarchy.input(file_name=pdb_file)
        hier = pdb_hier.hierarchy
        
        chains_to_check = list(set([resid.split()[2] for resid in self.residlst]))            
        
        coord = []
        info  = []
        for chain_to_check in chains_to_check:
            resids = [(resid.split()[0], resid.split()[1]) for resid in self.residlst if resid.split()[2]==chain_to_check]
            for chain in hier.chains():
                if chain.id == chain_to_check:
                    for res_group in chain.residue_groups():
                        for atom_group in res_group.atom_groups():
                            if (atom_group.resname, res_group.resseq.lstrip().rstrip()) in resids:
                                #print (atom_group.resname, res_group.resseq.lstrip().rstrip())
                                for a in res_group.atoms():
                                    #print a.name
                                    coord.append(list(a.xyz))
                                    i = a.fetch_labels()
                                    #info.append((i.resname, i.resseq, i.chain_id, i.altloc, i.name, i.i_seq))
                                    info.append((i.resname, i.resseq, i.chain_id, i.altloc, i.name))
                         
        coord = np.asarray(coord)
        info  = np.asarray(info)
 
        return coord, info
    
    def select_common_coords(self, coord, info):
        """
        Compare the info with the common atom info and only retain the coordinates of the common atoms
        """
        #join the info to get a 1d array
        info_ar = np.array([",".join(i) for i in info])
        #get the the indices of the common atoms
        _,_,indices_retain = np.intersect1d(self.common_info, info_ar, return_indices=True)
        
        if indices_retain.shape[0] < info.shape[0]:
            print("Trimming atom selection to those common between all pdb files")
        
        #select the coordinates and info from the common atoms
        coord = coord[indices_retain]
        info  = info[indices_retain]
        #print("info", info)
        assert info.shape[0] == coord.shape[0]
        
        return coord, info
    
    def get_d(self, p1, p2, axis=0):
        """
        :param p1: numpy array of dim (X,3) or (3)
        :param p2: numpy array of dim (3) or (3)
        :return:
        """
        p1 = np.array(p1)
        p2 = np.array(p2)
        
        return np.sqrt(np.sum((p1-p2)**2, axis=axis))
    
    def get_dist_matrix(self, coords):
        N = coords.shape[0]
        #print arr.shape
        diff = np.zeros((N, N))
        for i in range(N):
            CA = coords[i]
        #print CA
            diff[:,i] = self.get_d(coords, CA, axis=1)
        #print(diff)
        return diff
    
    def find_peaks(self, a,b):
        x = np.array(a)
        y = np.array(b)
        max = np.max(x)
        lenght = len(a)
        ret = []
        for i in range(lenght):
            ispeak = True
            if i-1 > 0:
                ispeak &= (x[i] > 1.5 * x[i-1]) # peak value is at least 150 percent of neighbour peak
            if i+1 < lenght:
                ispeak &= (x[i] > 1.5 * x[i+1])

            ispeak &= (x[i] > 0.25 * max)   # peak value is at least 25 percent of max value
            if ispeak:
                ret.append(y[i])
        if ret == []:
            ret = list(y[np.where(x==max)])
        return ret
    
    def remove_brackets_from_string(self, string):
        """
        Just remove all brackets from string
        """
        while '(' in string:
            string = re.sub('\(', '', string)
        while ')' in string:
            string = re.sub('\)','', string)
        return string

    def get_all_distances(self):#, mindiff = 0.05):
        """
        Calculate distances for atoms in of the residlist (or all atoms in case no residlist provided)
        """

        for pdb in self.pdblst:
            #Extract coordinates from residues in the residue list
            print("Extracting coordinates from %s" %(pdb))
            print("Extracting coordinates from %s" %(pdb), file=self.log)
            coord, info = self.get_coord_from_parser_selected(pdb)
            #only use coord and info from the common atoms
            coord, info = self.select_common_coords(coord, info)
            #Calculate the distance matric
            print("calculating distances")
            distances = np.array(self.get_dist_matrix(coord))
            #We actually only needs the upper triangle because all distances are double
            distances = np.triu(distances)
            #print("max distance:", np.max(distances))
            #Assemble the difference matrices from the different pdb-files
            #alldistances.append(distances)
            try:
                alldistances = np.vstack([alldistances, distances[np.newaxis,...]])
            except NameError:
                alldistances = distances[np.newaxis,...]
            #print("alldistances.shape", alldistances.shape)
            #Substract the first pdb-file (which is the reference)
            print("calculate difference")
            difference = np.subtract(distances, alldistances[0])
            #difference=distances-alldistances[0]
            #difference=np.array(difference)
            #Only keep the differences that are larger dan mindiff: this is not done here anymore: maxids and max_proj not used
            #print("Filter distances > %.2f" %(mindiff))
            #maxids = np.where(np.abs(difference)>mindiff)
            ##minids = np.where(difference<-mindiff)
            #max_proj=np.zeros(np.shape(difference))
            #max_proj[maxids] = difference[maxids] 
            #max_proj[minids] = difference[minids] 
            try:
                alldifferences= np.vstack([alldifferences, difference[np.newaxis,...]])
            except NameError:
                alldifferences = difference[np.newaxis,...]
            #print("alldifferences.shape", alldifferences.shape)
        self.info = info 
        self.alldistances   = alldistances
        self.alldifferences = alldifferences
        #self.alldistances = np.asarray(alldistances)
        #self.alldifferences = np.asarray(alldifferences)
        
        assert self.info.shape[0] == self.alldistances.shape[1] == self.alldistances.shape[2] == self.alldifferences.shape[1] == self.alldifferences.shape[2]

        try:
            print("max sampled distance: %s" %(np.max(alldistances)), file=self.log)
            print("properties of differences (max, mean, median): %s /// %s /// %s " %(np.max(alldifferences), np.mean(alldifferences), np.median(alldifferences)), file=self.log)
            print("max sampled distance: %s" %(np.max(alldistances)))
            print("properties of differences (max, mean, median): %s /// %s /// %s " %(np.max(alldifferences), np.mean(alldifferences), np.median(alldifferences)))
        except ValueError:
            #This happens in case the alldistances matrix is empty
            print("No distance differences observed", file=self.log)
            print("No distance differences observed")
                
    def plot_distance_differences_and_get_possible_b(self, maxdist = 10.0, mindist = 2.0, maxdiff = 10, mindiff=0.25):
        
        info_for_print = [[i[0],i[2],i[1].lstrip().rstrip(), '(%s%s)'%(i[4].lstrip().rstrip(), i[3].lstrip().rstrip())] for i in self.info]
        
        #Set all differences to a positive value
        self.alldifferences = np.abs(self.alldifferences)

        # get matrix with all interesting combinations distances for alpha calculation
        print("Filter distance differences between %.2f and %.2f" %(mindiff, maxdiff))
        for a in range(1,len(self.alphas)):
            s = np.where(self.alldifferences[a] <= maxdiff, 1, 0)
            t = np.where(self.alldifferences[a] >= mindiff, 1, 0)
            u = np.where(self.alldistances[a] < maxdist, 1, 0)
            v = np.where(self.alldistances[a] > mindist, 1, 0)
            w = np.add(np.add(s, t), np.add(u, v))
            toplot_alpha = np.triu(np.where(w==4, 1, 0))
            try:
                toplot_matrix = np.add(toplot_matrix, toplot_alpha)
            except NameError:
                toplot_matrix = toplot_alpha
              
        #we don't need self.alldistances anymore, so let's clear it already
        del self.alldistances
        
        toplot_matrix = np.where(toplot_matrix >=1, 1, 0)
        assert self.info.shape[0] == toplot_matrix.shape[0] == toplot_matrix.shape[1]
        
        plt.close()
        fig, (ax0,ax1) = plt.subplots(1,2, figsize=(10, 5))

        x = np.array(self.alphas)
        
        #In case fitting fails, we still have a matrix to continue with (but script will fail further downstream)
        #fitting_matrix = np.zeros((3,self.alldifferences.shape[1], self.alldifferences.shape[2]))
        
        #print("--------------------")
        print("Fitting distances.")
        #generate the fitting_matrix which contains on every n,m the popt from fitting all
        # n,m elements for each alpha. Each element is thus a list of 3 values.
        # Assert that with every operation the matrices have the correct shape
        
        #First set distances in alldistances which are too small to fit to 0 by multiplying with toplot_matrix:
        self.alldifferences = np.multiply(self.alldifferences, toplot_matrix)
        
        #set waring and error to ignore because we will to divide by zero
        np.seterr(divide='ignore', invalid='ignore')
        
        ##exponential fitting:
        #fitting_matrix = np.apply_along_axis(fitting, axis=0, arr=self.alldifferences, x =x)
        #sigmoidal fitting:
        fitting_matrix = np.apply_along_axis(sigmoid_fitting, axis=0, arr=self.alldifferences, x =x)
        
        assert fitting_matrix.shape[0] == 3
        assert fitting_matrix.shape[1] == fitting_matrix.shape[2] == self.alldifferences.shape[1]
        
        #calculate the residuals
        ##exponential fitting:
        #res_matrix = np.subtract(self.alldifferences, np.apply_along_axis(logfit_2, axis=0, arr=fitting_matrix, x=x))
        #sigmoidal fitting:
        res_matrix = np.subtract(self.alldifferences, np.apply_along_axis(sigmoid_fit_2, axis=0, arr=fitting_matrix, x=x))
        assert res_matrix.shape == self.alldifferences.shape
        
        #calculate sum of squares of the residuals
        sumofsquares_res_matrix = np.apply_along_axis(sumofsquares, axis=0, arr=res_matrix)
        assert sumofsquares_res_matrix.shape == self.alldifferences.shape[1:]
        
        #calculate total sum of squares of the differences
        totalsumsquares_matrix = np.apply_along_axis(totalsumsquares, axis=0, arr=self.alldifferences)
        assert totalsumsquares_matrix.shape == self.alldifferences.shape[1:]
        
        #calculate chi square
        chisq_matrix = np.nansum(np.divide(res_matrix, np.std(self.alldifferences, axis=0))**2, axis=0)
        assert chisq_matrix.shape == self.alldifferences.shape[1:]
        
        #calculate r square
        r_squared_matrix = np.subtract(1,np.divide(sumofsquares_res_matrix,totalsumsquares_matrix))
        assert r_squared_matrix.shape == self.alldifferences.shape[1:]
        r_squared_matrix = np.where(np.isnan(r_squared_matrix), 0, r_squared_matrix)
        assert r_squared_matrix.shape  == self.alldifferences.shape[1:]
        
        #reset normal warning messages upon zero-divisions
        np.seterr(divide='warn', invalid='warn')
        
        #we don't need intermediate matrices anymore, so let's clear them already
        del res_matrix, sumofsquares_res_matrix, totalsumsquares_matrix
        
        #update toplot_matrix with only the row-column pairs where r_square > 0.95 and 0 < chisq < 1.5:
        a = np.where((1.5 - chisq_matrix) > 0, 1, 0)
        b = np.where(r_squared_matrix > 0.95, 1, 0)
        toplot_matrix = np.multiply(a, b)
        #print("np.max(toplot_matrix)", np.max(toplot_matrix))
    
        #Extraction of possible b and amplitude. This assumes that both of them are non-zero
        #for simoidal plot, keep similar names as for the exponential fit avoid too much changes:
            #k = possible_b
            #L = possible_amp
        possible_b = np.multiply(fitting_matrix[1], toplot_matrix)
        possible_b = possible_b.ravel(order='C')[np.flatnonzero(possible_b)]
        
        possible_amp = np.multiply(fitting_matrix[0], toplot_matrix)
        possible_amp = possible_amp.ravel(order='C')[np.flatnonzero(possible_amp)]
        assert possible_b.shape[0] == possible_amp.shape[0]
        
        #sigmoidal fit need to extract the inflection point too:
        possible_x0 = np.multiply(fitting_matrix[2], toplot_matrix)
        possible_x0 = possible_x0.ravel(order='C')[np.flatnonzero(possible_x0)]
        assert possible_b.shape[0] == possible_amp.shape[0] == possible_x0.shape[0]
        
        #Calculate the average relevant distance
        fullist = np.multiply(self.alldifferences, toplot_matrix)
        fulllist_av = np.array([0])
        for a in range(1, len(self.alphas)):
            fulllist_av = np.append(fulllist_av,
                      np.mean(fullist[a].ravel(order='C')[np.flatnonzero(fullist[a])]))
        fulllist_av = np.where(np.isnan(fulllist_av), 0, fulllist_av)
        assert fulllist_av.shape[0] == len(self.alphas)
        #print("fulllist_av", fulllist_av)
        
        #only loop over all distances to plot if less than 400
        if np.where(np.any(toplot_matrix==1, axis =1))[0].shape[0] < 400:
            print("Less than 400 atoms involved in distance calculation. Let's plot all distances:", file=self.log)
            print("     Distance                                  R2   chi2   occupancy", file=self.log)
            print("Less than 400 atoms involved in distance calculation. Let's plot all distances:")
            print("     Distance                                  R2   chi2   occupancy")

            for n in np.where(np.any(toplot_matrix==1, axis =1))[0].astype(int):
            #n = int(n)
                toplot2 = np.where(toplot_matrix[n]==1)[0].astype(int)
                for j in toplot2:
                    a=0
                    title = "--".join(["-".join(info_for_print[n]),"-".join(info_for_print[j])])
                    ##exponential fitting:
                    #print('{:^40s} {:>10.2f} {:>5.2f} {:>5.2f}'.format(title, r_squared_matrix[n,j], chisq_matrix[n,j], 1/(val/fitting_matrix[1, n,j] + 1)), file=self.log)
                    #print('{:^40s} {:>10.2f} {:>5.2f} {:>5.2f}'.format(title, r_squared_matrix[n,j], chisq_matrix[n,j], 1/(val/fitting_matrix[1, n,j] + 1)))
                    #ax0.plot(x,self.alldifferences[:, n, j],color="%s"%(colorlib[a]), linestyle=':', linewidth=0.30, label='Distance')# ,marker='s', label='Method 2, {:.0%} occ.'.format(occ))          #+'; '+str(int(wavenumber[0]))+r' cm$^{-1}$')
                    #ax0.plot(x, np.abs(logfit_2(fitting_matrix[:,n, j], x)), 'r--', linewidth=0.20, label = 'Exp. fit')#label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
                    
                    #sigmoidal fitting:
                    print('{:^40s} {:>10.2f} {:>5.2f} {:>5.2f}'.format(title, r_squared_matrix[n,j], chisq_matrix[n,j], 1/((val/fitting_matrix[1, n,j])+ fitting_matrix[2, n,j])), file=self.log)
                    print('{:^40s} {:>10.2f} {:>5.2f} {:>5.2f}'.format(title, r_squared_matrix[n,j], chisq_matrix[n,j], 1/((val/fitting_matrix[1, n,j])+ fitting_matrix[2, n,j])))
                    ax0.plot(x,self.alldifferences[:, n, j],color="%s"%(colorlib[a]), linestyle=':', linewidth=0.30, label='Individual distance')
                    ax0.plot(x, np.abs(sigmoid_fit_2(fitting_matrix[:,n, j], x)), 'r--', linewidth=0.20, label = 'Individual fit')
                a+=1
        else:
            print("Too much distances, only plot averages", file=self.log)
            print("Too much distances, only plot averages")
                
        #plot the average distance
        ax0.plot(x,fulllist_av, label = "Average distance")
        xx = np.arange(1,max(self.alphas),0.1)
        #xx = np.arange(0,100,0.1)
        
        #Plot the fit of the average
        #If not all items of fulllist_av == 0
        if np.min(fulllist_av) < np.max(fulllist_av):
            ##exponential fitting
            #poptave,pcovave = curve_fit(logfit,
                                        #x,
                                        #fulllist_av,
                                        #sigma=np.sqrt(x)+x,
                                        #absolute_sigma=False,
                                        #bounds=((0.95*np.min(fulllist_av), val/maxalpha, -0.001), (1.05*np.max(fulllist_av),val, 0.001)))
            #Sigmoidal fitting
            #poptave,_ = curve_fit(sigmoid_fit,
                                        #x,
                                        #fulllist_av,
                                        #sigma=np.sqrt(x)+x,
                                        #absolute_sigma=False,
                                        #bounds=((0.95*np.min(fulllist_av), val/maxalpha -0.05, 1), (1.05*np.max(fulllist_av),val/1 +0.05, maxalpha)))
            poptave = sigmoid_fitting(fulllist_av, x)

            ##exponential fitting:
            #ax0.plot(xx, logfit(xx,poptave[0],poptave[1],0), 'b-', label='Fit of average fit')
            #sigmoidal fitting:
            ax0.plot(xx, sigmoid_fit(xx,poptave[0],poptave[1],poptave[2]), 'b-', label='Fit of average distance')

        else:
            poptave = np.array([0.,0.,0.])
            
        #Plot the average of all fits
        ##exponential fitting:
        #ax0.plot(xx, logfit(xx,np.median(possible_amp),(np.mean(possible_b)),0), 'b--', label='Average of all fits')
        #sigmoidal fitting:
        ax0.plot(xx, sigmoid_fit(xx,np.mean(possible_amp),np.mean(possible_b),np.mean(possible_x0)), 'b--', label='Average of all fits')


        #Set names of the plot axes 
        ax0.set_xlabel('alpha [1/occupancy]')
        ax0.set_ylabel('Distance [A]')
        
        #Set legend of the plot
        lines, labels = ax0.get_legend_handles_labels()
        i = np.arange(len(labels))
        f = np.array([])
        unique_labels = list(set(labels))
        for ul in unique_labels:
            f = np.append(f,[i[np.array(labels)==ul][0]])
        lines = [lines[int(j)] for j in f]
        labels = [labels[int(j)] for j in f]
        #ax.legend(lines, labels, fontsize = 'x-small', framealpha=0.5, loc=1, bbox_to_anchor=(0.05, 0.5, 0.5, 0.45))
        ax0.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.85, -0.05, 0.45, 0.5), fontsize = 'xx-small', framealpha=0.5)
        
        #Set title of the plot
        ax0.set_title("Distance differences", fontsize = 'medium')#,fontweight="bold")
        
        
        #Plot the histogram with possible occupancies
        ##exponential fitting:
        #possible_occ = np.round(1/(val/possible_b +1 ), decimals=3)
        #sigmoidal fitting:
        possible_occ = np.round(1/(((val)/possible_b)+ possible_x0), decimals=3)
        #print(possible_occ)
        counts, bins = np.histogram(possible_occ)
        ax1.hist(bins[:-1], bins, weights=counts)
        ax1.set_xlabel("Occupancy")
        if min(bins[:-1])-0.05 > 0:
            xmin = min(bins[:-1])-0.05
        else:
            xmin = 0
            
        if max(bins[:-1])+0.05 < 1:
            xmax = max(bins[:-1])+0.05
        else:
            xmax = 1
        ax1.set_xlim(xmin, xmax)
        #ax1.set_xlim([0, 1])
        ax1.set_ylabel("Occurence")
        ax1.set_title("Occupancy histogram", fontsize = 'medium')
        
        #Set figure title
        fig.suptitle("Distance analysis between reference model and real-space-refined models", fontsize = 'medium',fontweight="bold")
        #Set figure layour
        plt.subplots_adjust(hspace=0.35, wspace=0.40, left=0.09, right=0.95, top = 0.85)
        
        #plt.show()
        plt.savefig("Distance_difference_plot_%s.pdf" %(self.outsuffix), dpi=300, transparent=True)
        plt.savefig("Distance_difference_plot_%s.png" %(self.outsuffix), dpi=300)
        plt.close()
        
        ##exponential fitting:
        #return possible_b, poptave
        #sigmoidal fitting:
        return possible_b, poptave, possible_x0
        
    def extract_alpha(self):
        
        self.get_all_distances()
        if self.alldistances.shape[1] == 0:
            print("Occupancy cannot be determined",file=self.log)
            print("Occupancy cannot be determined")
            return 1.0, 1.0
        
        if self.resids_lst == None:
            ##exponential fitting:
            #possible_b, poptave = self.plot_distance_differences_and_get_possible_b()
            #sigmoidal fitting:
            possible_b, poptave, possible_x0 = self.plot_distance_differences_and_get_possible_b()
        else:
            ##exponential fitting:
            #possible_b, poptave = self.plot_distance_differences_and_get_possible_b()
            #sigmoidal fitting:
            possible_b, poptave, possible_x0 = self.plot_distance_differences_and_get_possible_b()
            
        if np.all(poptave == np.array([0.,0.,0.])):
            #print("poptave", poptave)
            print("Distance difference fitting failed", file=self.log)
            print("Distance difference fitting failed")
            return 1.0, 1.0
        
        if possible_b.shape[0] > 0:
            #Calculate b statistics
            b_mean = np.mean(possible_b)
            b_stdev = np.std(possible_b)
            counts, bins = np.histogram(possible_b)
            peaks = self.find_peaks(counts,bins)
            if len(peaks)>1:
                peaks = np.average(peaks)
            mode = scipy.stats.mode(np.round(possible_b,2))[0][0]
            
            #if sigmoidal fitting:
            #calculate x0 statistics
            x0_mean = np.mean(possible_x0)
            x0_stdev = np.std(possible_x0)
            counts_x0, bins_x0 = np.histogram(possible_x0)
            peaks_x0 = self.find_peaks(counts_x0,bins_x0)
            if len(peaks_x0)>1:
                peaks_x0 = np.average(peaks_x0)
            mode_x0 = scipy.stats.mode(np.round(possible_x0,2))[0][0]
            
            #set waring and error to ignore because possibility of division by zero
            np.seterr(divide='ignore', invalid='ignore')
            
            #Print results to log file and to stdout
            try:
                ##exponential fitting:
                #alp_mean = val/b_mean + 1
                #occ_mean = 1/alp_mean
                #occ_stdev = 1/(val/b_stdev)
                #occ_interval_max = occ_mean + occ_stdev
                #if occ_interval_max > 1:
                    #occ_interval_max = 1.0
                #occ_interval_min = occ_mean - occ_stdev
                #if occ_interval_min < 0:
                    #occ_interval_min= 0.0
                #occ_histogram = 1/(val/np.asarray(peaks, dtype='float64') + 1)
                #occ_mode = 1/(val/mode + 1)
                #occ_av_distance = 1/(val/poptave[1]+ 1)
                
                #sigmoidal fitting:
                alp_mean = ((val)/b_mean)+ x0_mean
                occ_mean = 1/alp_mean
                occ_stdev = 1/(val/b_stdev)
                occ_interval_max = occ_mean + occ_stdev
                if occ_interval_max > 1:
                    occ_interval_max = 1.0
                occ_interval_min = occ_mean - occ_stdev
                if occ_interval_min < 0:
                    occ_interval_min= 0.0
                
                occ_histogram = 1/(((val)/np.asarray(peaks, dtype='float64'))+ np.asarray(peaks_x0, dtype='float64'))
                occ_mode = 1/((val/mode)+ mode_x0)
                occ_av_distance = 1/((val/poptave[1])+ poptave[2])
                
                print("Occupancy based on average of fit results for %d distances: %.3f +/- %.3f" %(np.sum(counts), occ_mean, occ_stdev), file=self.log)
                print("  Occupancy intervals : %.3f  > %.3f >  %.3f " %(occ_interval_max, occ_mean, occ_interval_min), file=self.log)
                print("Other:", file=self.log)
                print("  Occupancy based on histogram of fit results for %d distances: %.3f" %(np.sum(counts),occ_histogram), file=self.log)
                print("  Occupancy based on most occuring value in histogram: %.3f" %(occ_mode), file=self.log)
                print("  Occupancy based on fit of average distance: %.3f" %(occ_av_distance), file=self.log)
                
                print("Occupancy based on average of fit results for %d distances: %.3f +/- %.3f" %(np.sum(counts), occ_mean, occ_stdev))
                print("  Occupancy intervals : %.3f  > %.3f >  %.3f " %(occ_interval_max, occ_mean, occ_interval_min))
                print("Other:")
                print("  Occupancy based on histogram of fit results for %d distances: %.3f" %(np.sum(counts),occ_histogram))
                print("  Occupancy based on most occuring value in histogram: %.3f" %(occ_mode))
                print("  Occupancy based on fit of average distance: %.3f" %(occ_av_distance))


            except TypeError:
                print("Problem with occupancy extraction", file=self.log)
                print("Problem with occupancy extraction")
                
            #reset warnind and error behavior
            np.seterr(divide='warn', invalid='warn')
                
        else:
            print("No distance differences within allowed range, alpha and occupancy cannot be estimated.", file=self.log)
            print("No distance differences within allowed range, alpha and occupancy cannot be estimated.")
            
            alp_mean = 1.0
            occ_mean = 1/alp_mean
        
        return alp_mean, occ_mean


class Filefinder(object):
    """
    Find input files for the standalone version of the distance_analysis, based on the output of an Xtrapol8 run.
    """
    def __init__(self,
                 X8_outdir ="Xtrapol8",
                 X8_outname = "Xtrapol8",
                 X8_list_occ = [0.1],
                 X8_f_extrapolated_and_maps = "qFextr",
                 distance_f_extrapolated_and_maps = "qFextr"):
        self.X8_outdir                  = X8_outdir
        self.X8_outname                 = X8_outname
        self.X8_list_occ                = X8_list_occ
        self.X8_f_extrapolated_and_maps = X8_f_extrapolated_and_maps
        self.distance_f_extrapolated_and_maps = distance_f_extrapolated_and_maps
    
    def find_pdbs(self):
        """
        Find the PDB files  given the X8_outdir, X8_outname, the X8_list_occ list and X8_f_extrapolated_and_maps 
        """
        if self.distance_f_extrapolated_and_maps not in self.X8_f_extrapolated_and_maps:
            print("ESFA file type not found in Xtrapol8 output (this might be a bug)")
        
        if self.distance_f_extrapolated_and_maps.startswith("q"):
            first_part = "qweight_"
        elif self.distance_f_extrapolated_and_maps.startswith("k"):
            first_part = "kweight_"
        else:
            first_part = ''
                        
        maptype = re.sub("f","F", self.distance_f_extrapolated_and_maps)
        last_part = "2m{:s}-DFc".format(maptype)
                
        pdb_list = []
        for occ in self.X8_list_occ:
            d = "{:s}/{:s}occupancy_{:.3f}/{:s}_occ{:.3f}_{:s}".format(self.X8_outdir, first_part, occ, self.X8_outname, occ, last_part)
            
            try:
                pdb_fles = glob.glob("{:s}*real_space*.pdb".format(d))
                pdb_fles.sort(key=lambda x: os.path.getmtime(x))
                pdb_out = pdb_fles[-1]
            except IndexError:
                pdb_out = glob.glob("{:s}*.pdb".format(d))
                pdb_fles.sort(key=lambda x: os.path.getmtime(x))
                pdb_out = pdb_fles[-1]
                
            pdb_list.append(os.path.abspath(check_file_existance(pdb_out)))
            print("occ: {:.3f} pdb found: {:s}".format(occ, pdb_out))

        return pdb_list
    
if __name__ == '__main__':

    import iotbx.phil
    from libtbx.utils import Usage
    
    from master import master_phil
    Xtrapol8_master_phil = master_phil

    master_phil = iotbx.phil.parse("""
    input{
        Xtrapol8_out = None
            .type = path
            .help = Xtrapol8_out.phil which can be found in the Xtrapol8 output directory
            .expert_level = 0
        f_extrapolated_and_maps = *qfextr fextr kfextr qfgenick fgenick kfgenick qfextr_calc fextr_calc kfextr_calc
            .type = choice(multi=False)
            .help = The type of ESFAs for which the difference map analysis will be carried out. The Xtrapol8 run prior to these analysis should include the ESFA type of choice. You can only specify one, launch mutliple runs if you want to repeat on with different ESFA types.
            .expert_level = 0
        reference_pdb = None
            .type = path
            .help = Reference model for difference map analysis
            .expert_level = 0
        }
    map_explorer{
        residue_list = None
            .type = path
            .help = list with residues to take into account for the occupancy estimation in same style as the output from map-explorer (e.g. residlist_Zscore2.00.txt). All residues in the reference model will be used if no residue_list is provided.
            .expert_level = 0
        use_waters = False
            .type = bool
            .help = Also use the water molecules in the analysis. This requires that waters have not been removed or added in comparison to the reference model
            .expert_level = 0
        }
    output{
        outdir = Distance_analysis
            .type = str
            .help = Output directory. 'Distance_analysis' will be used if not specified.
            .expert_level = 0
        suffix = None
            .type = str
            .help = suffix/prefix to be added to the output files (e.g. the Fextrapoled map type).
            .expert_level = 0
        log_file = None
            .type = str
            .help = write results to a file.
            .expert_level = 0
    }
    """, process_includes=True)

    #print help if no arguments provided or "--help" or "-h"
    if len(sys.argv) < 2:
           master_phil.show(attributes_level=1)
           raise Usage("phenix.python distance_analysis.py + [.phil] + [arguments]\n arguments only overwrite .phil if provided last")
           sys.exit(1)
    if "--help" in sys.argv or "-h" in sys.argv:
           master_phil.show(attributes_level=1)
           raise Usage("phenix.python distance_analysis.py + [.phil] + [arguments]\n arguments only overwrite .phil if provided last")
           sys.exit(1)

    #Extract input from inputfile and command line
    input_objects = iotbx.phil.process_command_line_with_files(
        args=sys.argv[1:],
        master_phil=master_phil
        )
    params = input_objects.work.extract()
        
    #Extract info form Xtrapol8 run
    Xtrapol8_input_objects = iotbx.phil.process_command_line_with_files(
        args = [params.input.Xtrapol8_out],
        master_phil = Xtrapol8_master_phil
        )
    Xtrapol8_params = Xtrapol8_input_objects.work.extract()

    model_pdb = os.path.abspath(check_file_existance(params.input.reference_pdb))
    
    #extract and search for distance-analysis parameters and input files
    pdbs = Filefinder(X8_outdir = Xtrapol8_params.output.outdir,
                      X8_outname = Xtrapol8_params.output.outname,
                      X8_list_occ = Xtrapol8_params.occupancies.list_occ,
                      X8_f_extrapolated_and_maps = Xtrapol8_params.f_and_maps.f_extrapolated_and_maps,
                      distance_f_extrapolated_and_maps = params.input.f_extrapolated_and_maps).find_pdbs()
                      
    pdblst = [model_pdb]+pdbs
    occupancies = Xtrapol8_params.occupancies.list_occ
    
    assert len(pdblst) == len(occupancies)+1, "Please provide a model PDB and an occupancy value for each pdb in the PDB list"

    if params.output.suffix != None:
        suffix = params.output.suffix
    else:
        suffix = params.input.f_extrapolated_and_maps
        
    use_waters = params.map_explorer.use_waters

    if params.map_explorer.residue_list != None:
        residue_list = os.path.abspath(params.map_explorer.residue_list)
    else:
        residue_list = None
        
    outdir = params.output.outdir
    i = 1
    while os.path.exists(outdir):
        if os.path.isdir(outdir):
            if len(os.listdir(outdir)) ==0:
                break
        outdir = "%s_%d" %(params.output.outdir, i)
        i += 1
        if i == 1000: #to avoid endless loop, but this leads to a max of 1000 runs
            break
    try:
        os.mkdir(outdir)
        print('Output directory being created: %s'%(outdir))
    except OSError:
        try:
            os.makedirs(outdir)
            print('Output directory being created: %s'%(outdir))
        except OSError:
            print("Output directory: %s" %(outdir))
            
    os.chdir(outdir)
    
    if params.output.log_file == None:
        log = open('distance_analysis_%s.log' %(suffix), 'w')
    else:
        log = open(params.output.log_file, 'w')
        
    #Run the distance analysis
    Distance_analysis(pdblst, occupancies, residue_list, use_waters = use_waters, outsuffix = suffix, log = log).extract_alpha()

    #Close the log file
    log.close()
