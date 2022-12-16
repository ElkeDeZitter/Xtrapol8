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

import os, re
import numpy as np

class Pymol_visualization(object):
    """
    class to open all output from Xtrapol8 calculations. Will search in outdir for the FoFo map, extrapolated map coefficients and output of reciprocal and real space refinement.
    It will only use ccp4 format map files as mtz cannot be loaded in open source pymol. Maps will be loaded but no isomesh or isosurface will be drawn.
    """
    def __init__(self, pdb_in, outdir):
        self.pdb_in = pdb_in
        self.outdir = outdir
    
    def get_basename(self, fle, extention):
        if "/" in fle:
            name = re.search(r"\/(.+?)\.%s" %(extention), fle).group(1).split("/")[-1]
        else:
            name = re.sub("\.%s" %(extention),"",fle)
        return name

        
    def find_fofo(self):
        """
        In outdir search for the FoFo maps in ccp4 format
        """
        maps = [self.outdir+"/"+fle for fle in os.listdir(self.outdir) if fle.lower().endswith("fofo.ccp4")]
        if len(maps) > 1:
            print("multiple FoFo difference maps found. All will be loaded")
        return maps

    def find_subdirectories(self):
        """
        in outdir find the names of the subd_subdirectories, these will contain the output of the extrapolated structure factor calculations
        """
        dirs =[dr for dr in os.listdir(self.outdir) if ('occupancy' in dr and os.path.isdir(self.outdir+"/"+dr))]
        return dirs
        
    def find_extrapolated_maps(self, dr):
        """
        In a subdirectory dr from outdir, search for the extrapolated maps in ccp4 format. If we also want to show the output of refinement, then the output should be written into ccp4 format
        """
        maps = [self.outdir+"/"+dr+"/"+fle for fle in os.listdir(self.outdir+"/"+dr) if 'ccp4' in fle]
        return maps
    
    def extract_maptype_from_ccp4_file(self, file_name):
        """
        From the ccp4 file name, extract the maptype. !!! Warning: This is very sensitive to changes to file names!!!!
        For (q)Fextr and 9q)Fextr_calc names:
        ccp4_name_2FoFc  = '%s_2m%s-DFc.ccp4' %(prefix, maptype)
        ccp4_name_FoFc   = '%s_m%s-DFc.ccp4' %(prefix, maptype)
        for (q)Fgenick names:
        ccp4_name_2FoFc  = '%s_m%s.ccp4' %(prefix, maptype)
        ccp4_name_FoFc   = '%s_m%s-DFc.ccp4' %(prefix, maptype)
        """
        if 'extr' in file_name.lower():
            try:
                maptype = re.search(r"2m(.+?)-dfc", file_name.lower()).group(1)
            except AttributeError:
                try:
                    maptype = re.search(r"m(.+?)-dfc", file_name.lower()).group(1)
                except AttributeError:
                    return None
        elif 'genick' in file_name.lower():
            try:
                maptype = re.search(r"m(.+?)-dfc", file_name.lower()).group(1)
            except AttributeError:
                try:
                    maptype = re.search(r"m(.+?)\.ccp4", file_name.lower()).group(1)
                except AttributeError:
                    return None
        else: #if not extr or genick in filename
            maptype = None
        return maptype
                    
    def extract_maptype_from_pdb_file(self, file_name):
        """
        From the pdb file name, extract the maptype. !!! Warning: This is very sensitive to changes to file names!!!!
        """
        
        if 'extr' in file_name.lower():
            try:
                maptype = re.search(r"_m(.+?)-dfc", file_name.lower()).group(1)
            except AttributeError:
                return None
        elif 'genick' in file_name.lower():
            try:
                maptype = re.search(r"fgenick_m(.+?)-dfc", file_name.lower()).group(1)
            except AttributeError:
                return None
        else: #if not extr or genick in filename
            maptype = None
        return maptype
    
    def find_real_space_refined_models(self, dr):
        """
        In a subdirectory dr from outdir, search for the pdb files which come from the real space refinement into the extrapolated map (no reciprocal space or real space after reciprocal space refinement)
        """
        models = [self.outdir+"/"+dr+"/"+fle for fle in os.listdir(self.outdir+"/"+dr) if 'pdb' in fle and ('reciprocal' and 'refmac') not in fle]
        return models
    
    def open_all_in_pymol(self):
        """
        Write script that opens all files defined above. Script can be run using 'Pymol <where/it/is>/pymol_all.py
        """
        
        script_pymol = self.outdir+'/pymol_all.py'
        i = open(script_pymol, 'w')
        i.write("from pymol import cmd\ncmd.set('group_auto_mode', 1)\n")
        
        #load dark model
        name_dark = self.get_basename(self.pdb_in, 'pdb')
        i.write("cmd.load('%s','%s')\n" %(self.pdb_in, name_dark))
        
        #load Fo-Fo difference maps
        i.write("cmd.group('FoFo_map')\n")
        for mp in self.find_fofo():
            name = self.get_basename(mp, 'ccp4')
            i.write("cmd.load('%s', 'FoFo_map.%s')\n"%(mp, name))
        
        #load extrapolated maps and models
        for dr in self.find_subdirectories():
            i.write("cmd.group('%s')\n" %(dr))
            maptypes = []
            for mp in self.find_extrapolated_maps(dr):
                name = self.get_basename(mp, 'ccp4')
                maptype = self.extract_maptype_from_ccp4_file(name)
                if maptype != None:
                    if maptype not in maptypes: #if maptype not in the list yet
                        i.write("cmd.group('%s.%s')\n" %(dr, maptype)) #create the subgroup
                        maptypes.append(maptype) #append to map so that the subgroup will not be re-created
                    i.write("cmd.load('%s', '%s.%s.%s_map')\n"%(mp, dr, maptype, name))
                else: #maptype could not be extracted
                    i.write("cmd.load('%s', '%s.%s_map')\n"%(mp, dr, name))
            for model in self.find_real_space_refined_models(dr):
                name = self.get_basename(model, 'pdb')
                maptype = self.extract_maptype_from_pdb_file(name)
                if maptype != None:
                    if maptype not in maptypes: #if maptype not in the list yet
                        i.write("cmd.group('%s.%s')\n" %(dr, maptype)) #create the subgroup
                        maptypes.append(maptype) #append to map so that the subgroup will not be re-created
                    i.write("cmd.load('%s', '%s.%s.%s')\n"%(model, dr, maptype, name))
                else: #maptype could not be extracted
                    i.write("cmd.load('%s', '%s.%s')\n"%(model, dr, name))
                    
        i.write("cmd.center('%s')"%(name_dark))
        
        i.close()
        return script_pymol


class Pymol_movie(object):
    """
    Aim is to merge Pymol visualization 
    """
    def __init__(self,
                 occupancies,
                 pdblst=None,
                 ccp4_maps=None,
                 resids_lst=None,
                 model_label="",
                 ccp4_map_label=""):
        
        self.resids_lst     = resids_lst
        self.model_label    = model_label
        self.ccp4_map_label = ccp4_map_label
        
        if pdblst != None:
            #Sort the pdb files and alphas in order to have alpha from small to large, this is important for fitting
            #This might not work in pyhton3
            alphas = list(map(lambda x: round(1/x, 3), occupancies))
            if ccp4_maps!=None:
                zipped = zip(alphas, pdblst, ccp4_maps)
                zipped_sorted = sorted(zipped, key = lambda x:x[0])
                alphas, pdblst, ccp4_maps =zip(*zipped_sorted)
                alphas = list(alphas)
                self.pdblst = list(pdblst)
                self.ccp4_maps = list(ccp4_maps)
            else:
                zipped = zip(alphas, pdblst)
                zipped_sorted = sorted(zipped, key = lambda x:x[0])
                alphas, pdblst =zip(*zipped_sorted)
                alphas = list(alphas)
                self.pdblst = list(pdblst)    
                self.ccp4_maps = None
        else:
            self.pdblst = None
            self.ccp4_maps = None
            print("No models selected for pymol visualization.")
        
    def get_residlist(self):
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
            self.residlst = np.array(residlst_unique)
            
    def get_syntax(self):
        self.resids = False
        self.resid_lines = []
        self.waters = False
        self.water_lines = []
        for resid in self.residlst:
            try:
                resn, resv, chain, alt = resid.split("\n")[0].split()[:4]
                if resn == 'HOH':
                    self.waters = True
                    self.water_lines.append('(chain %s and resid %s and altloc %s)'%(chain, resv, alt))
                else:
                    self.resids = True
                    self.resid_lines.append('(chain %s and resid %s and altloc %s)'%(chain, resv, alt))
            except ValueError:
                resn, resv, chain= resid.split("\n")[0].split()[:3]
                if resn == 'HOH':
                    self.waters = True
                    self.water_lines.append('(chain %s and resid %s)'%(chain, resv))
                else:
                    self.resids = True
                    self.resid_lines.append('(chain %s and resid %s)'%(chain, resv))

            
    def load_models(self,outscript):
        for i,pdb in enumerate(self.pdblst):
            outscript.write('cmd.load("%s","models_%s",state=%d)\n'%(pdb, self.model_label, i+1))
    
    def show_default(self, outscript):
        outscript.write('cmd.hide()\n')
        outscript.write('cmd.show("lines")\n')
        outscript.write('cmd.show("nonbonded")\n')
                        
    def show_sticks_and_spheres(self, outscript):
        if self.resids_lst == None:
            print("Showing all residues in pymol movie.")
            outscript.write('cmd.show("sticks")\n')
            outscript.write('cmd.show("nb_spheres")\n')
        else:
            if self.resids:
                outscript.write('cmd.show("sticks", "%s")\n'%(" or ".join(self.resid_lines)))
            if self.waters:
                outscript.write('cmd.show("nb_spheres", "%s")\n'%(" or ".join(self.water_lines)))
                
    def load_maps(self, outscript):
        for i, ccp4 in enumerate(self.ccp4_maps):
            outscript.write('cmd.load("%s","maps_%s",%d,"ccp4")\n'%(ccp4,self.ccp4_map_label, i+1))

    def show_mesh(self, outscript):
        if self.resids_lst == None:
            print("No residues selected for isomesh.")
            # if len(self.residlist)<50:
            #     print("No residues selected for isomesh. The molecule is small thus showing isomesh for all residues")
            #     outscript.write('cmd.isomesh("%s", "maps_%s", 1.0, "all", carve=1.6)\n'%(self.ccp4_map_label, self.ccp4_map_label))
            # else:
            #     print("No residues selected for isomesh.")
        else:
            if (self.waters and self.resids):
                selection = " or ".join([" or ".join(self.resid_lines)," or ".join(self.water_lines)])
            elif self.resids:
                selection = " or ".join(self.resid_lines)
            else:
                selection = " or ".join(self.water_lines)
            outscript.write('cmd.isomesh("%s", "maps_%s", 1.0, "%s", carve=1.6)\n'%(self.ccp4_map_label, self.ccp4_map_label, selection))
            
    def additional_settings(self, outscript):
        outscript.write('cmd.set("movie_fps", 5)\n')
        outscript.write('cmd.show("mesh", "F*")\n')
        outscript.write('cmd.show("mesh", "qF*")\n')
        outscript.write('cmd.set("mesh_width", 0.3)\n')
        
    def remove_old_outscript(self, outscript):
        try:
            os.remove(outscript)
        except OSError:
            print("%s doen not exist, hence cannot be removed.")
    
    def write_pymol_script(self, outfile="pymol_movie.py"):
        self.get_residlist()
        if self.resids_lst != None:
            self.get_syntax()
        o = open(outfile, "a")
        self.load_models(o)
        if self.ccp4_maps != None:
            self.load_maps(o)
            self.show_mesh(o)
        o.close()
        #return outfile
        
    def write_pymol_appearance(self, outfile="pymol_movie.py"):
        self.get_residlist()
        if self.resids_lst != None:
            self.get_syntax()
        o = open(outfile, "a")
        self.show_default(o)
        self.show_sticks_and_spheres(o)
        self.additional_settings(o)
        o.close()
        #return outfile
    
    
    
