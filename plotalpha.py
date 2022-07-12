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

from __future__ import print_function
import sys, os, re
import glob
import numpy as np
from matplotlib import pyplot as plt
import pickle

class plotalpha(object):
    """
    Calculation of the alpha and occupancy based on the integrated difference maps peaks stored in map_explorer output files.
    """
    def __init__(self, map_expl_files, resids_lst=None, outsuffix = 'qFextr', log=sys.stdout):
        self.map_expl_files = map_expl_files
        self.resids_lst     = resids_lst
        self.outsuffix      = outsuffix
        self.log            = log
        self.alphafound     = False
        
        occupanies = [(re.search(r'occupancy\_(.+?)\/',map_expl_fle).group(1)) for map_expl_fle in self.map_expl_files[1:]]
        occupancies = map(lambda x: float(x), occupanies)
        self.occupancies = np.array(occupancies)
        self.alphas = np.reciprocal(occupancies)
        
        #open residlist to calculated alpha using selected residues
        self.get_residlist()
        
        print("MAP EXPLORER ANALYSIS")
        print("MAP EXPLORER ANALYSIS", file=self.log)

    
    def get_residlist(self):
        """
        Extract the residues from the residue list
        """
        if self.resids_lst != None:
            with open(self.resids_lst) as rs_lst:
                residlst = rs_lst.readlines()
                residlst = [lne for lne in residlst if len(lne)>0]
            if len(residlst) == 0: 
                self.resids_lst = None
            residlst = np.array(residlst)
            residlst_unique = list(set(residlst))
            residlst_unique.remove([lne for lne in residlst_unique if 'Resn Resv Chain Alt' in lne][0])
            residlst = np.array(residlst_unique)
            self.residlst = residlst
            
    def get_final_list(self, map_expl_fle):
        """
        Extract the lines from the map explorer file and store those those only have the residues that are in the in the residlist in an additional list
        """
        with open(map_expl_fle) as fle:
            lines = fle.readlines()
        lines = [line for line in lines if len(line.split())>1 and "Resn" not in line]
        lines = np.array(lines) #lines contains the usefull lines from the peakintegration file
        finallst = [] #finallist will contains all the lines from "lines" for the residues in residuelist
        if self.resids_lst == None: #If no residlist provided, then all residues will be used
            finallst = lines 
        else:
            for resid in self.residlst:
                test =  [line for line in lines if ','.join(resid.split()) in ','.join(line.split())] #make sure all residues in residue list are in peakintegration file. This avoids addition of residues for which no peaks are found because this would mess up the analysis
                for lne in test:
                    if len(lne)>0:
                        finallst.append(lne) #make sure no empty lines are present
        return lines, finallst
    
    def get_integrated_values(self, map_expl_fle):
        """
        Get the integrated values for the selected and all residues.
        Use the integrated values of all residues as an estimate for the noise in the map
        """
        print("Extracting peak information from %s"%(map_expl_fle), file=self.log)
        print("Extracting peak information from %s"%(map_expl_fle))
        lines, finallst = self.get_final_list(map_expl_fle)
        tmp   = 0
        tmp2  = 0
        intnoise = 0
        for line in finallst:
            tmp  += np.abs(float(line.split()[-3])) #sum all the peak-values to tmp for residues in finallist
        for line in lines:
            tmp2 += np.abs(float(line.split()[-3])) #sum all the peak-values to tmp2
            intnoise+=np.abs(float(line.split()[-3])) #all peaks will be used to calculate the noise
        return tmp, tmp2, intnoise
            
    def normalize_array(self, array):
        return array/np.max(array)
    
    def normalize_array_double(self, array):
        return (array-np.min(array))/(np.max(array)-np.min(array))
    
    def get_integration_lists(self):
        """
        Store the integrated values as signal-to-noise values in same order as alphas and occupancies in appropriate arrays
        Take care that the noise is calculated as the sqrt of the sum of the intnoise instead of average of the intnoise
        Usage of the arrays depends on the method for the alphadetermination
        """
        self.int1 = np.zeros(self.alphas.shape[0])
        self.int2 = np.zeros(self.alphas.shape[0])
        self.exp1 = np.zeros(1)
        self.exp2 = np.zeros(1)
        for i,expl_fle in enumerate(self.map_expl_files):
            tmp, tmp2, intnoise = self.get_integrated_values(expl_fle)
            if i == 0: #i==0 means that 'expl_fle' is the peaklist from the Fo-Fo map
                try:
                    noise = np.sqrt(intnoise)
                    #Devision by noise to do some kind of normalization, usefull when map_explorer parameters were not superwill chosen
                    self.exp1[0]  = tmp / noise #'experimental' peak-height = sum(all selected)/sqrt(all peaks)
                    self.exp2[0]  = tmp2 / noise #'experimental' peak-height = sum(all)/sqrt(all peaks)
                except: #exception in which case?
                    self.exp1[0] = tmp
                    self.exp2[0] = tmp2
                    
            else: #i > 0 means that expl_fle is from the extrapolated difference maps
                try: 
                    noise = np.sqrt(intnoise)
                    self.int1[i-1] = tmp / noise #'integrated' peak-height = sum(all selected)/sqrt(all peaks)
                    self.int2[i-1] = tmp2 / noise #''integrated' peak-height = sum(all)/sqrt(all peaks)
                except: #exception in which case? ZeroDivisionError?
                    self.int1[i-1] = tmp 
                    self.int2[i-1] = tmp2
                print("occ: %.3f alpha:%.3f integrated peak area selected: %.3f integrated peak area all: %.3f" %(self.occupancies [i-1], self.alphas[i-1], self.int1[i-1], self.int2[i-1]))
                
        self.int1_norm = self.normalize_array(self.int1)
        self.int2_norm = self.normalize_array(self.int2)
        self.exp1_norm = self.normalize_array(np.array(self.exp1))
        self.exp2_norm = self.normalize_array(np.array(self.exp2))
        
    def alphadetermination(self):
        """
        Alphadetermination with signal enhancement. Uses selected peaks only
        """
        try:
            results = self.int1_norm * self.normalize_array(1 - np.sqrt(((self.int1_norm - self.exp1_norm)/(self.int1_norm+self.exp1_norm))**2))
            self.results = self.normalize_array(results) #normalization
        except ZeroDivisionError:
            print("no peaks found in at least one of the maps")
            self.results = self.int1_norm 
        gooda = np.where(self.results==np.max(self.results))
            
        try:
            a = np.max(self.alphas[gooda])
            o = np.max(self.occupancies[gooda])
            self.alphafound = True
        except:
            print("No good alpha value could be found. Best alpha will be set to 1, occupancy will be set to 1 and subsequently corrected to the highest value in the input list. Consider changing map_explorer parameters and run again.", file=self.log)
            print("No good alpha value could be found. Best alpha will be set to 1, occupancy will be set to 1 and subsequently corrected to the highest value in the input list. Consider changing map_explorer parameters and run again.")
            a = 1
            o = 1
            self.alphafound = False
        return a, o
    
    def alphadetermination_special(self):
        """
        Alternative alpha determnation from original plotalpha special script.
        Alphadetermination with deconvolution of all peaks with selected peaks
        """
        try:
            results = self.int1_norm * self.int2_norm
            self.results = self.normalize_array(results)
        except ZeroDivisionError:
            print("no peaks found in at least one of the maps")
            self.results = self.int1_norm
        gooda = np.where(self.results==np.nanmax(self.results)) #ignore nan if present
        
        try:
            a = np.max(self.alphas[gooda])
            o = np.max(self.occupancies[gooda])
            self.alphafound = True
        except:
            print("No good alpha value could be found. Best alpha will be set to 1, occupancy will be set to 1 and subsequently corrected to the highest value in the input list. Consider changing map_explorer parameters and run again.", file=self.log)
            print("No good alpha value could be found. Best alpha will be set to 1, occupancy will be set to 1 and subsequently corrected to the highest value in the input list. Consider changing map_explorer parameters and run again.")
            a = 1
            o = 1
        return a, o
    
    def alphadetermination_simple(self):
        """
        Alternative alpha determination from initial plotalpha script
        Alphadetermination with selected peaks only
        """
        try:
            results = (self.int1 / self.exp1) + np.average(self.int1)
            self.results = self.normalize_array(results)
        except ZeroDivisionError:
            print("no peaks found in at least one of the maps")
            self.results = self.int1_norm
        gooda = np.where(self.results==np.nanmax(self.results))
        
        try:
            a = np.max(self.alphas[gooda])
            o = np.max(self.occupancies[gooda])
            self.alphafound = True
        except:
            print("No good alpha value could be found. Best alpha will be set to 1, occupancy will be set to 1 and subsequently corrected to the highest value in the input list. Consider changing map_explorer parameters and run again.", file=self.log)
            print("No good alpha value could be found. Best alpha will be set to 1, occupancy will be set to 1 and subsequently corrected to the highest value in the input list. Consider changing map_explorer parameters and run again.")
            a = 1
            o = 1
        return a, o
    
    def alphadetermination_kyp(self):
        """
        Alternative alpha determination with double normalization according to Kyprianos suggestion.
        Alphadermination with selected peaks only
        """
        results = self.int1 - self.exp1
        try:
            self.results = self.normalize_array_double(results)
        except ZeroDivisionError:
            print("no peaks found in at least one of the maps")
            self.results = self.int1_norm
        gooda = np.where(self.results==np.nanmax(self.results))
        
        try:
            a = np.max(self.alphas[gooda])
            o = np.max(self.occupancies[gooda])
            self.alphafound = True
        except:
            print("No good alpha value could be found. Best alpha will be set to 1, occupancy will be set to 1 and subsequently corrected to the highest value in the input list. Consider changing map_explorer parameters and run again.", file=self.log)
            print("No good alpha value could be found. Best alpha will be set to 1, occupancy will be set to 1 and subsequently corrected to the highest value in the input list. Consider changing map_explorer parameters and run again.")
            a = 1
            o = 1
        return a, o
    
    def plot_alpha_determination(self, outname):
        
        #write pickle for GUI
        out=open('alpha_occupancy_determination_%s.pickle' %(self.outsuffix) ,'wb') #write to pickle for GUI
        stats = [self.alphas,self.occupancies, self.int1_norm, self.int2_norm, self.results]
        if self.resids_lst == None:
            stats.append(False)
        else:
            stats.append(True)
        
        if self.alphafound:
            stats.append(True)
        else:
            stats.append(False)
        pickle.dump(stats,out)
        out.close()

        fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10, 5))
        
        if self.resids_lst == None:
            ax1.plot(self.alphas, self.int1, 's', markersize = 5, color = 'blue', label='All peaks') #int1, int2 and conv will be the same, so we can plot only int1 and avoid the try catch
        else:
            try:
                ax1.plot(self.alphas, self.int1_norm, 'o', color = 'red', label='Selected residues')      
                ax1.plot(self.alphas, self.int2_norm, 's', markersize = 5, color = 'blue', label='All peaks')
                ax1.plot(self.alphas, self.results, '^', color="green", label='Selected residues with enhanced SNR')
            except: #in case int1 and int2 don't exist, when would this be the case?
                ax1.plot(self.alphas, self.int1_norm, 'o', color = 'red')

        ax1.set_ylim([0.,1.1])
        ax1.set_xlim([np.min(self.alphas)*0.95,np.max(self.alphas)*1.05])
        ax1.set_xlabel('Alpha value = 1/occupancy')
        ax1.set_ylabel('Normalized difference map ratio')

        if self.resids_lst == None:
            ax2.plot(self.occupancies, self.int1_norm, 's', markersize = 5, color = 'blue', label='All peaks') #int1, int2 and conv will be the same, so we can plot only int1 and avoid the try catch
        else:
            try:
                ax2.plot(self.occupancies, self.int1_norm, 'o', color = 'red', label='Selected residues') 
                ax2.plot(self.occupancies, self.int2_norm, 's', markersize = 5, color = 'blue', label='All peaks')     
                ax2.plot(self.occupancies, self.results, '^', color="green", label='Selected residues with enhanced SNR')
            except:
                ax2.plot(self.alphas,self.int1, 'o', color = 'red')
                
        ax2.set_ylim([0.,1.1])
        ax2.set_xlim([np.min(self.occupancies)*0.95,np.max(self.occupancies)*1.05])
        ax2.set_xlabel('Triggered state occupancy')
        ax2.set_ylabel('Normalized difference map ratio')
        ax2.legend(loc='lower right', bbox_to_anchor=(0.92, -0.05, 0.45, 0.5), fontsize = 'x-small', framealpha=0.5)

        if self.alphafound:
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
    
    def estimate_alpha(self):
        self.get_integration_lists()
        print("mFext-DFc maps          | Sum of all integrated peaks | Sum of selected integrated peaks", file=self.log)
        print("Occupancy       alpha         absolute   normalized         absolute   normalized", file=self.log)
        print("mFext-DFc maps          | Sum of all integrated peaks | Sum of selected integrated peaks")
        print("Occupancy       alpha         absolute   normalized         absolute   normalized")

        for i in range(self.alphas.shape[0]):
            print("{:>5.3f} {:>15.3f} {:>15.3f} {:>10.3f} {:>18.3f} {:>10.3f}".format(self.occupancies[i],self.alphas[i], self.int2[i], self.int2_norm[i], self.int1[i], self.int1_norm[i]), file=self.log)
            print("{:>5.3f} {:>15.3f} {:>15.3f} {:>10.3f} {:>18.3f} {:>10.3f}".format(self.occupancies[i],self.alphas[i], self.int2[i], self.int2_norm[i], self.int1[i], self.int1_norm[i]))

        print("Fobs-Fobs map {:>23.3f} {:>10.3f} {:>18.3f} {:>10.3f}".format(self.exp2[0], self.exp2_norm[0], self.exp1[0], self.exp1_norm[0]), file=self.log)
        print("Fobs-Fobs map {:>23.3f} {:>10.3f} {:>18.3f} {:>10.3f}".format(self.exp2[0], self.exp2_norm[0], self.exp1[0], self.exp1_norm[0]))

        #print("Standard alphadetermination", file=self.log)
        self.alpha, self.occ = self.alphadetermination()
        #alternatives:
        #print("Special alphadetermination", file=self.log)
        #self.alpha, self.occ = self.alphadetermination_special()
        #print("Simple alphadetermination", file=self.log)
        #self.alpha, self.occ = self.alphadetermination_simple()
        #print("Kyprianos alphadetermination", file=self.log)
        #self.alpha, self.occ = self.alphadetermination_kyp()
        print("Best alpha value is %.3f, meaning an occupancy of %.3f"%(self.alpha, self.occ),file=self.log)
        print("Best alpha value is %.3f, meaning an occupancy of %.3f"%(self.alpha, self.occ))
        
        self.plot_alpha_determination('alpha_occupancy_determination_%s' %(self.outsuffix))
        
        return self.alpha, self.occ

