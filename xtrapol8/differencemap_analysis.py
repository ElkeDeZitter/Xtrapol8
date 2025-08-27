# -*- coding: utf-8 -*-
"""
Script to calculate the occupancy based on the Fourier difference map (FoFo) and extrapolated difference maps (mFextr-DFcalc).
This is automatically run in Xtrapol8 routine but can be run on a standalone basis using this script.

What do you need?
- a complete Xtrapol8 run in Calm_and_curious or Fast_and_furious mode
- an ESFA type (this should have been run in the Xtrapol8 run that precedes the differencemap_analysis
- map-explorer parameters
- Optional:
    - Residue list for which the difference map peaks will be used to estimate the occupancy (in the format as an Xtrapol8 residuelist)
        If no residue list provided, a residue list will be calculated based on the Z-score
    - suffix to be added to the output files (will be used as a prefix for some files)
    - Log-file
    - Output directory

usage
-----
To get the help message:
$ phenix.python difference_analysis.py
or
$ phenix.python difference_analysis.py --help
or
$ phenix.python difference_analysis.py -h

To run with an input file:
$ phenix.python difference_analysis.py difference_analysis.phil

To run with command line arguments:
$ phenix.python difference_analysis.py input.Xtrapol8=Xtrapol8/Xtrapol8_out.phil input.f_extrapolated_and_maps=qfextr
 map_explorer.residue_list=residlist_adapted.txt map_explorer.peak_integration_floor=4 map_explorer.peak_detection_threshold=4 output.outdir=diffmap_test

To run with an input file and command line arguments:
$ phenix.python difference_analysis.py difference_analysis.phil input.f_extrapolated_and_maps=fextr_calc
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
from libtbx import adopt_init_args
import mmtbx
from mmtbx import utils
from cctbx import sgtbx, crystal
from iotbx import pdb
from iotbx.file_reader import any_file
from iotbx import ccp4_map
import iotbx.xplor.map
from scipy.stats import pearsonr
from map_explorer import map_explorer
from map_explorer_analysis import Map_explorer_analysis
from plotalpha import plotalpha
from Fextr_utils import check_file_existance
from libtbx.utils import Usage

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
        
        #initiate the mask as None
        self.mask = None
        
    
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

    
    def do_map_explorer_and_get_data_FoFo(self):
        """
        Run map_explorer on the FoFo map and append the peakintegration results to a list
        If no residue list is provided: search for the highest peaks in the map using the Z-score
        Make the secondary structure plot (ss-plot) with the difference map peaks
        Write the actual map data to an object as it will be required to calculate the Pearson correlation coef.
        """
        print("************Map explorer: %s ************" %(self.fofo), file=self.log)
        print("************Map explorer: %s ************" %(self.fofo))
        
        #Run map_explorer
        map_expl_out_FoFo, residlist_zscore, mask, FoFo_ref = map_explorer(self.fofo, self.pdb_in, self.radius, self.threshold, self.peak, self.z_score, self.residue_list)
                
        map_expl_out_FoFo = os.path.abspath(check_file_existance(map_expl_out_FoFo))

        print("FoFo map explored. Results in %s" %(map_expl_out_FoFo), file=self.log)
        print("FoFo map explored. Results in %s" %(map_expl_out_FoFo))
        
        #Generate plot which indicates the secondary structure and integrated peak volume
        Map_explorer_analysis(peakintegration_file = map_expl_out_FoFo, ligands = self.extract_ligand_codes(),log=self.log).get_ss(self.pdb_in)
        print("Difference map plot generated")
        
        #Extract the FoFo data
        if self.fofo.endswith('ccp4'):
            self.fofo_data = ccp4_map.map_reader(file_name=self.fofo).data.as_numpy_array()
        else:
            self.fofo_data = iotbx.xplor.map.reader(file_name=self.fofo).data.as_numpy_array()
        
        self.map_exp_files.append(FoFo_ref)
        self.mask = mask
        
    def get_map_explorer_results_FextrFcalc(self, mask = None):
        """
        Run map_eplorer on each mFextr-DFcalc map and append results to the list with peakintegration files
        """
        
        start_dir = os.getcwd()
        
        for i,fextrfcalc in enumerate(self.fextrfcalc_list):
            print("\n************Map explorer: %s ************" %(fextrfcalc), file=self.log)
            print("occupancy: %.3f" %(self.occupancies[i]), file=self.log)
            print("\n************Map explorer: %s ************" %(fextrfcalc))
            print("occupancy: %.3f" %(self.occupancies[i]))
            
            #fextrfcalc = os.path.abspath(fextrfcalc)
            
            #Need to work in different directories because
            #1) map_explorer writes always the same filename
            #2) plotalpha used the occupancy determined from the directory
            #occ = self.occupancies[i]
            #out_dir = "%s_occupancy_%.3f" %(self.prefix, occ)
            #self.check_and_make_dir(out_dir)
            #os.chdir(out_dir)
            
            ##Run map_explorer
            #map_expl_out = self.run_map_explorer(fextrfcalc, map_type = self.prefix)
            #map_expl_out = os.path.abspath(check_file_existance(map_expl_out))
            #print("FoFo map explored. Results in %s" %(map_expl_out), file=self.log)
            #print("FoFo map explored. Results in %s" %(map_expl_out))
            
            if fextrfcalc.endswith('ccp4'):
                data = ccp4_map.map_reader(file_name=fextrfcalc).data.as_numpy_array()
            else:
                data = iotbx.xplor.map.reader(file_name=fextrfcalc).data.as_numpy_array()
            pos = 0
            neg = 0
            #print(mask[0,0])
            
            try:
                for i in range(self.mask.shape[1]):
                    tmp = data[self.mask[0, i]].sum()
                    if tmp > 0: pos+= tmp
                    else: neg -= tmp
            except AttributeError: #when self.mask = None  as set in the __init__
                print("Mask not found. Run get_map_explorer_results_FoFo first and make sure it ran correctly.")
                return
                
            try:
                CC = pearsonr(self.fofo_data.flatten(), data.flatten())[0]
            except ValueError:
                #the fft_map function might change the cystal gridding. In that case the maps do not have the same shape
                #and thus the CC cannot be calculated. Not nice, but at least Xtrapol8 can continue.
                #in this case, also the output from plotalpha is wrong since the mask will be incorrectly projected!!!
                print("Pearson correlation factor could not be calculated. The CC will be set to zero.")
                CC = 0

            self.map_exp_files.append([CC, pos, neg, pos+neg])
            
            #os.chdir(start_dir)
            
    def run_difference_map_analysis(self):
        """
        Perform the different steps of the difference map analys:
        1) map explorer
        2) plotalpha
        """
        print("map_explorer parameters:", file=self.log)
        print("  peak_integration_floor: %f" %(self.threshold), file=self.log)
        print("  peak_detection_threshold: %f" %(self.peak), file=self.log)
        print("  radius: %f" %(self.radius), file=self.log)
        if self.residue_list == None:
            print("  z_score: %f" %(self.z_score), file=self.log)
        else:
            print("  residue list: %s" %(self.residue_list), file=self.log)
        
        print("---------------------------------------------", file=self.log)
        #Map_explorer and map_explorer_analysis on the FoFo map
        self.do_map_explorer_and_get_data_FoFo()
        
        #Map_explorer om the FextrFcalc maps
        self.get_map_explorer_results_FextrFcalc()
        
        #Estimate alha and occupancy based on the peakintegration area as stored in the peakintegration files
        print("---------------------------------------------",file=self.log)
        _, _ , _, _= plotalpha(self.occupancies, self.map_exp_files[1:], self.map_exp_files[0], self.prefix, log=self.log).estimate_alpha()
    
        print("---------------------------------------------", file=self.log)
        
class Filefinder(object):
    def __init__(self,
                 X8_outdir ="Xtrapol8",
                 X8_outname = "Xtrapol8",
                 X8_list_occ = [0.1],
                 X8_fofo_type = "qFoFo",
                 X8_f_extrapolated_and_maps = "qFextr",
                 diffmap_f_extrapolated_and_maps = "qFextr"):
        self.X8_outdir                  = X8_outdir
        self.X8_outname                 = X8_outname
        self.X8_list_occ                = X8_list_occ
        self.X8_fofo_type               = X8_fofo_type
        self.X8_f_extrapolated_and_maps = X8_f_extrapolated_and_maps
        self.diffmap_f_extrapolated_and_maps = diffmap_f_extrapolated_and_maps
        
    def find_fofo(self):
        """
        Find the FoFo map based on the X8_outdir and X8_fofo_type and X8_outname
        """
        maptype = re.sub("f","F", self.X8_fofo_type)
        last_part = "m{:s}.ccp4".format(maptype)

        f = "{:s}/{:s}_{:s}".format(self.X8_outdir, self.X8_outname,last_part)
        fofo = os.path.abspath(check_file_existance(f))
        
        return fofo
    
    def find_fextfc(self):
        """
        Find the Fextr-Fc type of maps given the X8_outdir, X8_outname, the X8_list_occ list and X8_f_extrapolated_and_maps 
        """
        if self.diffmap_f_extrapolated_and_maps not in self.X8_f_extrapolated_and_maps:
            print("ESFA file type not found in Xtrapol8 output (this might be a bug)")
        
        if self.diffmap_f_extrapolated_and_maps.startswith("q"):
            first_part = "qweight_"
        elif self.diffmap_f_extrapolated_and_maps.startswith("k"):
            first_part = "kweight_"
        else:
            first_part = ''
            
        maptype = re.sub("f","F", self.diffmap_f_extrapolated_and_maps)
        last_part = "m{:s}-DFc.ccp4".format(maptype)
                
        fextrfcalc_list = []
        for occ in self.X8_list_occ:
            f = "{:s}/{:s}occupancy_{:.3f}/{:s}_occ{:.3f}_{:s}".format(self.X8_outdir, first_part, occ, self.X8_outname, occ, last_part)
            fextrfcalc_list.append(os.path.abspath(check_file_existance(f)))
            
        return fextrfcalc_list
            
if __name__ == "__main__":
    
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
        }
    map_explorer{
        residue_list = None
            .type = str
            .help = list with residues to take into account for the occupancy estimation in same style as the output from map-explorer in an Xtrapol8 run (e.g. residlist_Zscore2.00.txt). If no file is provided, a residue list will be generated based on the map_explorer analysis and Z-score value.
            .expert_level = 0
        peak_integration_floor = 3.5
            .type = float
            .help = Floor value for peak integration (sigma). Peaks will be integrated from their maximum value towards this lower bound to avoid integration of noise.
            .expert_level = 0
        peak_detection_threshold = 4.0
            .type = float
            .help = Peak detection threshold (sigma). Only peaks with an absolute value equal or above this value will be integrated.
        radius = None
            .type = float
            .help = Maximum radius (A) to allocate a density blob to a protein atom in map explorer. Resolution will be used if not specified.
            .expert_level = 0
        z_score = 2.0
            .type = float
            .help = Z-score to determine residue list with only highest peaks.
            .expert_level = 0
        }
    output{
        outdir = Differencemap_analysis
            .type = str
            .help = Output directory. 'Differencemap_analysis' be used if not specified.
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
           raise Usage("phenix.python differencemap_analysis.py + [.phil] + [arguments]\n arguments only overwrite .phil if provided last")
           sys.exit(1)
    if "--help" in sys.argv or "-h" in sys.argv:
           master_phil.show(attributes_level=1)
           raise Usage("phenix.python differencemap_analysis.py + [.phil] + [arguments]\n arguments only overwrite .phil if provided last")
           sys.exit(1)

    #Extract input from inputfile and command line
    input_objects = iotbx.phil.process_command_line_with_files(
        args=sys.argv[1:],
        master_phil=master_phil
        )
    params = input_objects.work.extract()
    
    #Extract info form Xtrapol8 run
    if params.input.Xtrapol8_out == None:
        print("input.Xtrapol8_out not defined")
        sys.exit(1)
    if os.path.isfile(params.input.Xtrapol8_out) == False:
        print("File not found: {:s}". format(params.input.Xtrapol8_out))
        sys.exit(1)
              
    Xtrapol8_input_objects = iotbx.phil.process_command_line_with_files(
        args = [params.input.Xtrapol8_out],
        master_phil = Xtrapol8_master_phil
        )
    Xtrapol8_params = Xtrapol8_input_objects.work.extract()
    
    #extract and search for difference-map-analysis parameters and input files
    fofo_map = Filefinder(X8_outdir = Xtrapol8_params.output.outdir,
                          X8_outname = Xtrapol8_params.output.outname,
                          X8_fofo_type = Xtrapol8_params.f_and_maps.fofo_type).find_fofo()
    
    model_pdb = Xtrapol8_params.input.reference_pdb
    
    additional_files = Xtrapol8_params.input.additional_files

    fextrfcalcs = Filefinder(X8_outdir = Xtrapol8_params.output.outdir,
                             X8_outname = Xtrapol8_params.output.outname,
                             X8_f_extrapolated_and_maps = Xtrapol8_params.f_and_maps.f_extrapolated_and_maps,
                             X8_list_occ = Xtrapol8_params.occupancies.list_occ,
                             diffmap_f_extrapolated_and_maps = params.input.f_extrapolated_and_maps).find_fextfc()
    
    occupancies = Xtrapol8_params.occupancies.list_occ
    
    if len(fextrfcalcs) != len(occupancies):
        print("Number of occupancies and mFextr-DFcalc maps is not equal. Please provide a single occupancy for each map.")
        sys.exit(1)

        
    if params.map_explorer.residue_list != None:
        if os.path.isfile(params.map_explorer.residue_list) == False:
            print("File not found: {:s}". format(params.map_explorer.residue_list))
            sys.exit(1)
        residue_list = os.path.abspath(params.map_explorer.residue_list)
    else:
        residue_list = None
        
    radius    = params.map_explorer.radius
    threshold = params.map_explorer.peak_integration_floor
    peak      = params.map_explorer.peak_detection_threshold
    z_score   = params.map_explorer.z_score
    
    if params.output.suffix != None:
        suffix = params.output.suffix
    else:
        suffix = params.input.f_extrapolated_and_maps
            
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
        log = open('differencemap_analysis_%s.log' %(suffix), 'w')
    else:
        log = open(params.output.log_file, 'w')
        
        
    #Run the difference-map analysis
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


    
