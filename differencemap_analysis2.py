# -*- coding: utf-8 -*-
"""
Script to calculate the occupancy based on the Fourier difference map (FoFo) and extrapolated difference maps (mFextr-DFcalc).
This is automatically run in Xtrapol8 routine but can be run on a standalone basis using this script.

What do you need?
- An FoFo map
- A reference model
- list with mFextr-DFcalc map files.
- list with occupancies in the same order as the associated mFextr-DFcalc map files
- Optional:
    - Residue list for which the difference map peaks will be used to estimate the occupancy (in the format as an Xtrapol8 residuelist)
        If no residue list provided, a residue list will be calculated based on the Z-score
    - Outsuffix to be added to the output files (will be used as a prefix for some files)
    - Log-file. If not provided, then output will written to a predifined file
    
    
usage
-----
To get the help message:
$ phenix.python difference_analysis.py

To run with a specific set of arguments:
$ phenix.python difference_analysis.py -f my_fofo.map -m my_model.pdb -x occ1/mfextr-Dfcalc_with_occ0.1.map,occ2/mfextr-Dfcalc_with_occ0.2.map -o 0.1,0.2

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

from __future__ import print_function
import sys
import os
import re
import argparse
import numpy as np
import scipy.stats
from libtbx import adopt_init_args
import mmtbx
from mmtbx import utils
from cctbx import sgtbx, crystal
from iotbx import pdb, ccp4_map
import iotbx.xplor.map
from iotbx.file_reader import any_file
from map_explorer2 import map_explorer
#from map_explorer_analysis import Map_explorer_analysis
from plotalpha2 import plotalpha
from Fextr_utils import check_file_existance


class Difference_analysis(object):
    def __init__(self,
                 pdb_in,
                 fofo,
                 additional_files = [],
                 fextrfcalc_list = [],
                 occupancies = [],
                 residue_list = None,
                 radius = 10,
                 threshold = 3.5,
                 peak = 4.0,
                 z_score= 2.0,
                 prefix = '',
                 log = sys.stdout):
    
        #Try this instead of  repetition of the arguments
        adopt_init_args(self, locals())
        assert len(fextrfcalc_list) == len(occupancies), 'the number of difference maps and occupancies is not equal, this will be nonsense and end with an error somewhere'
        
        #Sort the map files and occupancies in order to have occupancies from small to large, probably just for cosmethics
        #This might not work in pyhton3
        zipped = zip(occupancies, fextrfcalc_list)
        zipped_sorted = sorted(zipped, key = lambda x:x[0])
        occupancies, fextrfcalc_list =zip(*zipped_sorted)
        self.occupancies     = list(occupancies)
        self.fextrfcalc_list = list(fextrfcalc_list)
        
        #initiate an empty list with the peakintegration files:
        self.map_exp_files = []
        
        #check if all additional files are provided, this is needed to extract the ligand codes:
        self.check_additional_files()
        
    
    def check_additional_files(self):
        """
        Check for all additional files if they exist, if the minimum of additional files is present for phenix
        """
        #Loop over the additional_files and check if it is a cif-file
        self.cif_objects = []
        if len(self.additional_files)>0:
            for fle in self.additional_files:
                if 'cif' in fle:
                    try:
                        cif_object=mmtbx.monomer_library.server.read_cif(file_name=fle)
                    except Exception:
                        #raise  AssertionError,"Unable to read the cif file %s" %(fle)
                        print("Unable to read the cif file: %s, this file will be ignored. This might affect some other parts in the analysis." %(fle), file=self.log)
                    else:
                        if (len(cif_object) > 0):
                            self.cif_objects.append((fle,cif_object))
                            
        #Check if all cif-files are provided to interpret the pdb-file
        model_in = any_file(self.pdb_in, force_type="pdb", raise_sorry_if_errors=True)
        SpaceGroup=sgtbx.space_group_info(symbol=str(model_in.crystal_symmetry().space_group_info()))
        crystal_symmetry=crystal.symmetry(unit_cell=model_in.crystal_symmetry().unit_cell(),
                                          space_group_info=SpaceGroup)
        processed_pdb_files_srv = utils.process_pdb_file_srv(
            crystal_symmetry          = crystal_symmetry,
            pdb_parameters            = pdb.input(self.pdb_in),
            cif_objects               = self.cif_objects)
        processed_pdb_files_srv.process_pdb_files(pdb_file_names = [self.pdb_in])

        
    def extract_ligand_codes(self):
        """
        Get list with three-letter codes for ligands. Needed in further analysis.
        """
        #Read the additional files and check if these are sufficient to completely interprete the pdb file (aka, are all ligand codes provided)
        #self.check_additional_files()
        
        #Extact the ligand codes from the cif-files
        cif_list = []
        for cif in self.cif_objects:
            for comp in cif[1]:
                try:
                    ligand = re.search(r'comp_(.+?)$', comp).group(1)
                    if len(ligand) == 3 and ligand not in cif_list:
                        cif_list.append(ligand)
                except AttributeError:
                    continue
            for ligand in cif[1]['comp_list']['_chem_comp.id']:
                if ligand not in cif_list:
                    cif_list.append(ligand)
        return cif_list
    
    def check_and_make_dir(self, outdir):
        """
        Create directory if it doesn't exist yet
        """
        
        if os.path.exists(outdir) == False:
            os.mkdir(outdir)
            print('Output directory: %s'%(outdir))
            print('Output directory: %s'%(outdir), file=self.log)

    def run_map_explorer(self, difference_map, map_type = ""):
        """
        Run map_explorer and obtain the file containing the peakintegration information
        """
        
        map_expl_file, residue_list, mask, ref = map_explorer(difference_map, self.pdb_in, self.radius, self.threshold, self.peak, self.z_score)
        
        return map_expl_file, residue_list, mask, ref
    
    def get_map_explorer_results_FoFo(self):
        """
        Run map_explorer on the FoFo map and append the peakintegration file to a list
        If no residue list is provided: search for the highest peaks in the map using the Z-score
        Make the secondary structure plot (ss-plot) with the difference map peaks
        """
        print("************Map explorer: %s ************" %(self.fofo), file=self.log)
        print("************Map explorer: %s ************" %(self.fofo))
        
        #Run map_explorer
        map_expl_out_FoFo, self.residue_list, self.mask, self.ref = self.run_map_explorer(self.fofo)
        map_expl_out_FoFo = os.path.abspath(check_file_existance(map_expl_out_FoFo))
        print("FoFo map explored. Results in %s" %(map_expl_out_FoFo), file=self.log)
        print("FoFo map explored. Results in %s" %(map_expl_out_FoFo))
        
        #append the absolute path to the list of map_explorer files
        self.map_exp_files.append(map_expl_out_FoFo)
        
        #Get the Z-score list
        #if (self.residue_list == None or os.path.abspath(self.residue_list) == False):
        #    residlist_zscore  = Map_explorer_analysis(peakintegration_file = map_expl_out_FoFo,log=self.log).residlist_top(Z=self.z_score)
        #    residlist_zscore  = os.path.abspath(check_file_existance(residlist_zscore))
        #    print("Residue list in residlist.txt and residues associated to highestpeaks in %s\n" %(residlist_zscore), file=self.log)
        #    print("Residue list in residlist.txt and residues associated to highestpeaks in %s\n" %(residlist_zscore))
        self.residue_list = os.path.abspath(check_file_existance(self.residue_list))
        #else:
        #    self.residue_list = os.path.abspath(self.residue_list)
        #    print("Residue list provided: %s" %(self.residue_list), file=self.log)
        #    print("Residue list provided: %s" %(self.residue_list))

        # TODO: Generate plot which indicates the secondary structure and integrated peak volume
        #Map_explorer_analysis(peakintegration_file = map_expl_out_FoFo, ligands = self.extract_ligand_codes(),log=self.log).get_ss(self.pdb_in)
        #print("Difference map plot generated")


    def get_map_explorer_results_FextrFcalc(self):
        """
        Run map_eplorer on each mFextr-DFcalc map and append results to the list with peakintegration files
        """
        
        #start_dir = os.getcwd()
        self.integrated_values = []
        self.pearsonCC = []
        try:
            self.fofo_data = ccp4_map.map_reader(file_name=self.fofo).data.as_numpy_array()
        except:
            self.fofo_data = iotbx.xplor.map.reader(file_name=self.fofo).data.as_numpy_array()
        for i, fextrfcalc in enumerate(self.fextrfcalc_list):
            print("\n************Map explorer: %s ************" %(fextrfcalc), file=self.log)
            print("\n************Map explorer: %s ************" %(fextrfcalc))
            
            fextrfcalc = os.path.abspath(fextrfcalc)
            
            #Need to work in different directories because
            #1) map_explorer writes always the same filename
            #2) plotalpha used the occupancy determined from the directory
            #occ = self.occupancies[i]
            #out_dir = "%s_occupancy_%.3f" %(self.prefix, occ)
            #self.check_and_make_dir(out_dir)
            #os.chdir(out_dir)
            if fextrfcalc.endswith('ccp4'):
                density_map = ccp4_map.map_reader(file_name=fextrfcalc)
            else:
                density_map = iotbx.xplor.map.reader(file_name=fextrfcalc)
            data = density_map.data.as_numpy_array()

            pos = 0
            neg = 0
            for i in range(self.mask.shape[1]):
                tmp = data[self.mask[0, i]].sum()
                if tmp > 0: pos+= tmp
                else: neg -= tmp
            self.integrated_values.append([pos, neg, pos+neg])
            self.pearsonCC.append(scipy.stats.pearsonr(self.fofo_data.flatten(), data.flatten())[0])
            #Run map_explorer
            #map_expl_out = self.run_map_explorer(fextrfcalc, map_type = self.prefix)
            #map_expl_out = os.path.abspath(check_file_existance(map_expl_out))
            #print("FoFo map explored. Results in %s" %(map_expl_out), file=self.log)
            #print("FoFo map explored. Results in %s" %(map_expl_out))
            
            #self.map_exp_files.append(map_expl_out)
            
            #os.chdir(start_dir)
        self.integrated_values = np.array(self.integrated_values)
        self.pearsonCC = np.array(self.pearsonCC)
        #print(self.integrated_values.size)
        #print(self.pearsonCC.size)

            
    def run_difference_map_analysis(self):
        """
        Perform the different steps of the difference map analys:
        1) map explorer
        2) plotalpha
        """
        
        print("---------------------------------------------", file=self.log)
        #Map_explorer and map_explorer_analysis on the FoFo map
        self.get_map_explorer_results_FoFo()
        
        #Map_explorer om the FextrFcalc maps
        self.get_map_explorer_results_FextrFcalc()
        
        #Estimate alha and occupancy based on the peakintegration area as stored in the peakintegration files
        print("---------------------------------------------", file=self.log)
        alpha, occ = plotalpha(self.occupancies, self.integrated_values, self.pearsonCC, self.ref, self.prefix, log=self.log).estimate_alpha()
    
        print("---------------------------------------------", file=self.log)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = 'Standalone version of the difference map analysis for occupancy estimation.')
    
    #input files
    parser.add_argument('-f', '--fofo_map', default="my_input.map", help="Fourier difference map (FoFo) in xplor format (.map)")
    parser.add_argument('-m', '--model_pdb', default='input.pdb', help='Reference coordinates in pdb format.')
    parser.add_argument('-a', '--additional_files', default = None, help='Additional files required for correct interpration of the pdb file, e.g ligand cif file, restraints file. Comma-seperated, no spaces.')
    parser.add_argument('-x', '--fextrfcalc_list', default='mfextr-Dfcalc_with_occ0.1.map,mfextr-Dfcalc_with_occ0.2.map', help='list of mfextr-Dfcalc map files in xplor format (.map) to be analysed. Comma-seperated, no spaces.')
    parser.add_argument('-o', '--occupancies', default ='0.1, 0.2', help='list of occupancies, in the same order as the mfextr-Dfcalc map files. Comma-seperated, no spaces.')
    
    #map_explorer parameters
    parser.add_argument('-rl', '--residue_list', default = None, help='list with residues to take into account for the occupancy estimation in same style as the output from map-explorer (e.g. residlist_Zscore2.00.txt). If no file is provided, a residue list will be generated based on the map_explorer analysis and Z-score value.')
    parser.add_argument('-r', '--radius', default = 2.0, type=float, help="maximum radius in Angstrom to allocate a density blob to a protein atom (e.g. the resolution of the extrapolated maps)")
    parser.add_argument('-t', '--peak_integration_floor', default = 3.5, type=float,  help="integration threshold in sigma")
    parser.add_argument('-p', '--peak_detection_threshold', default = 4.0, type=float, help="Peak detection threshold in sigma")
    parser.add_argument('-z', '--z_score', default = 2.0, type=float, help='Z-score to determine residue list with only highest peaks')

    #auxiliary paramters
    parser.add_argument('-s', '--suffix', default = '', help='suffix/prefix to be added to the output files (e.g. the Fextrapoled map type).')
    parser.add_argument('-l', '--log_file', default=None, help='write results to a file.')
    
    
    #print help if no arguments provided
    if len(sys.argv) < 2:
           parser.print_help()
           sys.exit(1)

    #interprete arguments
    args = parser.parse_args()
    
    fofo_map         = os.path.abspath(args.fofo_map)
    model_pdb        = os.path.abspath(args.model_pdb)
    additional_files = args.additional_files.split(",")
    fextrfcalcs      = args.fextrfcalc_list.split(",")
    occupancies      = list(map(lambda x : float(x), args.occupancies.split(",")))
    if len(fextrfcalcs) != len(occupancies):
        print("Number of occupancies and mFextr-DFcalc maps is not equal. Please from a single occupancy for each map.")
        sys.exit()
        
    residue_list = args.residue_list
    radius       = args.radius
    threshold    = args.peak_integration_floor
    peak         = args.peak_detection_threshold
    z_score      = args.z_score
    suffix       = args.suffix
    
    
    if args.log_file == None:
        log = open('differencemap_analysis_%s.log' %(suffix), 'w')
    else:
        log = open(args.log_file, 'w')
        
        
    #Run the analysis
    Difference_analysis(model_pdb,
                        fofo_map,
                        additional_files,
                        fextrfcalcs,
                        occupancies,
                        residue_list,
                        radius,
                        threshold,
                        peak,
                        z_score,
                        suffix,
                        log).run_difference_map_analysis()
            
    #Close the log file
    log.close()


    