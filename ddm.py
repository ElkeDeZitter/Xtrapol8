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
import numpy as np
import argparse
import os, sys
import math
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from Fextr_utils import get_name
from iotbx.pdb import hierarchy
from cctbx.array_family import flex

class Difference_distance_analysis(object):
    """
    Class to calculate the ddm between two structures.
    Ligands are not (yet) taken into account
    """
    def __init__(self, pdb1, pdb2, ligands = [], outdir=None, scale=1.5, log=sys.stdout):
        self.ligands = ligands
        self.scale   = scale
        self.log     = log
        
        if outdir==None:
            self.outdir = os.getcwd()
        else:
            self.outdir = outdir
        
        self.pdb1_name = get_name(pdb1)
        print("reference pdb file: %s" %(pdb1), file=self.log)
        pdb1_hier = hierarchy.input(file_name=pdb1)
        self.hier1 = pdb1_hier.hierarchy
        print("number of atoms at start: %s" %(self.hier1.atoms_size()), file=self.log)
        self.ligand_coords1, self.ligand_info1 = self.remove_ligands_and_get_coord(self.hier1)
        print("number of atoms after ligand and water removal: %d" %(self.hier1.atoms_size()), file=self.log)
        #self.remove_altlocs(self.hier1)
        self.hier1.remove_alt_confs(True)
        print("number of atoms after altloc removal: %d" %(self.hier1.atoms_size()), file=self.log)
        print("number of protein atoms for ddm calculation: %d" %(self.hier1.atoms_size()), file=self.log)
        
        self.pdb2_name = get_name(pdb2)
        print("refined pdb file: %s" %(pdb2), file=self.log)
        pdb2_hier = hierarchy.input(file_name=pdb2)
        self.hier2 = pdb2_hier.hierarchy
        print("number of atoms at start: %d" %(self.hier2.atoms_size()), file=self.log)
        self.ligand_coords2, ligand_info2 = self.remove_ligands_and_get_coord(self.hier2)
        print("number of atoms after ligand and water removal: %d" %(self.hier2.atoms_size()), file=self.log)
        #self.remove_altlocs(self.hier2)
        self.hier2.remove_alt_confs(True)
        print("number of atoms after altloc removal: %d" %(self.hier2.atoms_size()), file=self.log)
        print("number of protein atoms for ddm calculation: %d" %(self.hier2.atoms_size()), file=self.log)
        
    def remove_ligands_and_get_coord(self, pdb_hierarchy):
        """
        Function to remove ligands from pdb_hierarchy
        remove ligands because they might mess up with numbers
        remove water molecules
        """
        ligand_coords = []
        ligand_info   = []
        
        hetatm_residues = []
        for chain in pdb_hierarchy.chains():
            for res_group in chain.residue_groups():
                for atom_group in res_group.atom_groups():
                    # for each atom of the pdb file:
                    for atom in atom_group.atoms():
                        # 1.add the HETATM residues to the list of ligands
                        if atom.hetero:
                            # in the pdb file, if the line of the atom starts with HETATM
                            hetatm_residues.append(atom_group.resname)
                    # add the residue name to the list hetatm_residues
                    self.ligands = hetatm_residues + self.ligands
                    # get rid of the repetition of residue name of the list of ligands
                    self.ligands = list(set(self.ligands))
                    # 2.remove ligands
                    if (atom_group.resname in self.ligands and atom_group.resname!= 'HOH'):
                        #print(chain.id, atom_group.resname)
                        for a in res_group.atoms():
                            ligand_coords.append(list(a.xyz))
                            i = a.fetch_labels()
                            ligand_info.append((i.resname, i.resseq, i.chain_id, i.altloc, i.name, i.i_seq))
                        res_group.remove_atom_group(atom_group)
                    if atom_group.resname == 'HOH':
                        res_group.remove_atom_group(atom_group)
                    if res_group.atoms_size() == 0:
                        chain.remove_residue_group(res_group)
        return ligand_coords, ligand_info
    
    #def remove_altlocs(self, pdb_hierarchy):
        #"""
        #Function to remove laternative conformations
        #"""
        #for chain in pdb_hierarchy.chains():
            #if chain.is_protein():
                #for res_group in chain.residue_groups():
                    #if len(res_group.conformers())>1:
                        #first_altloc = res_group.conformers()[0].altloc
                        #for atom_group in res_group.atom_groups():
                            #if (atom_group.altloc != "" and atom_group.altloc != first_altloc):
                                ##print("remove:",atom_group.id_str())
                                #res_group.remove_atom_group(atom_group)
    
    def get_offset(self, chain):
        """
        Function to get the offset of a certain chain
        """

        try:
            offset = chain.residues()[0].resseq_as_int()
        except AssertionError:
            offset = chain.conformers()[0].residues()[0].resseq_as_int()
                    
        return offset
    
    def get_last_residue_number(self, pdb_hierarchy, chain):
        
        for c in pdb_hierarchy.chains():
            if c.id == chain.id:
                try:
                    last = c.residues()[-1].resseq_as_int()
                except AssertionError:
                    last = c.conformers()[0].residues()[-1].resseq_as_int()
                    
        return last

                    
    def get_coord(self, pdb_hierarchy, chain):
        """
        Get the coordinates of the atoms in the pdb file.
        Because of badly placed TER cards, a chain may be read as multiple chains, hence need to loop over all chains again instead of working with the input chain. Need to find a proper way to merge chains with same ID or get rid of TER cards.
        """
        
        coord = []
        info  = []
        for c in pdb_hierarchy.chains():
            if c.id == chain.id:
                for res_group in c.residue_groups():
                    for a in res_group.atoms():
                        #print a.name
                        coord.append(list(a.xyz))
                        i = a.fetch_labels()
                        info.append(i.resseq)
                         
        coord = np.asarray(coord)
 
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

    
    def get_diff_matrix(self, arr):
        N = arr.shape[0]
        #print(arr.shape)
        diff = np.zeros((N, N))
        for i in range(N):
            CA = arr[i]
            diff[:,i] = self.get_d(arr, CA, axis=1)
        return diff

                    
    def calculate_ddm(self, chain1, chain2):

        #print "Processing first pdb" # Reference
        coords1, seq_info = self.get_coord(self.hier1, chain1)
        diff_m1 = self.get_diff_matrix(coords1)

        #"Processing second pdb" -
        coords2,_ = self.get_coord(self.hier2, chain2)
        diff_m2 = self.get_diff_matrix(coords2)
        
        try:
            diff = diff_m2 - diff_m1
        except ValueError:
            print('Different number of Calpha between pdb1 (%s, %i) and pdb2 (%s, %i)\n'
                'Please check your pdb files.'%(self.pdb1_name, diff_m1.shape[0], self.pdb2_name, diff_m2.shape[0]), file=self.log)
            diff = np.array([0,0])
        return diff, seq_info
    
    def ddm_residue(self, ddm, seq_info):
        """
        From the all atom ddm, calculate a new ddm that returns only the average distance per residue pair
        """
        seq_info_unique, indices = np.unique(np.array(map(lambda x: int(x), seq_info)), return_inverse=True)
        
        n_residues = seq_info_unique.shape[0]
        ddm_residue = np.zeros((n_residues, n_residues))
        for i in range(n_residues):
            for j in range(n_residues):
                ddm_residue[i,j] = np.average(ddm[np.where(indices==i)[0],:][:,np.where(indices==j)[0]])
                
        return ddm_residue, seq_info_unique
    
    def make_colormap(self,seq):
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
    
    def ddms(self):
        """
        Calculate ddm
        """
        chains = [chain.id for chain in self.hier1.chains() if chain.is_protein()]
        num_chains = len(set(chains))
            
        if num_chains == 1:
            n_cols = 1
        else:
            n_cols = 2
            
        n_rows_prot = int(math.ceil(num_chains/n_cols))

        #if num_lig == 0:
            #n_rows_lig = 0
        #else:
            #n_rows_lig = 1
        #n_rows = n_rows_prot + n_rows_lig
        
        n_rows = n_rows_prot

        outname = '%s/ddm_%s_vs_%s'%(self.outdir, self.pdb1_name,self.pdb2_name)
        outname_pdf = '%s.pdf'%(outname)
        outname_png = '%s.png'%(outname)
        
        #Also write the ddm_residue and seq_info_unique to a pickle file
        outname_pickle = open("%s.pickle"%(outname), 'ab')

        fig, axs = plt.subplots(n_rows, n_cols, figsize=(10*n_rows, 10*n_cols), squeeze=False)#, constrained_layout=True)
        #fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)
        c = mcolors.ColorConverter().to_rgb
        rvb = self.make_colormap([c('blue'), c('white'), 0.40, c('white'),  0.60, c('white'), c('red')])


        col = -1
        row = -1
        old_chain_id = ''
        for chain in self.hier1.chains():
            if chain.is_protein():
                ID = chain.id
                if ID != old_chain_id:
                    if col == 0:
                        col = 1
                    else:
                        col = 0
                        row +=1
                    chain1 = chain.detached_copy()
                    offset = self.get_offset(chain1)
                    
                    for c in self.hier2.chains():
                        if c.is_protein():
                            if (c.id == ID and self.get_offset(c) == offset):
                                chain2 = c.detached_copy()
                            
                    ddm, seq_info = self.calculate_ddm(chain1, chain2)
                    
                    try:
                        n,m = ddm.shape
                        if ddm.shape[0] == 0:
                            raise ValueError
                    except ValueError:
                        print("ddm cannot be generated.")
                        #Need to return the png outname so that it can be easily found by the GUI. This can be much more elegant though
                        return outname_png
                    
                    assert n==m, "Difference matrix for chain %s is incorrect"
                    
                    #For all atom ddm
                    #mask = np.zeros_like(ddm, dtype=np.bool)
                    #mask[np.triu_indices_from(mask)] = True
                    #FINAL2 = np.ma.array(ddm, mask=mask)
                    
                    #last = self.get_last_residue_number(self.hier1, chain1)
                    #total = last+offset+1
                    ##steps = int(total/10)
                    #seq_info = np.array(list(map(lambda x: int(x), seq_info)))
                    #seq_ticks = seq_info[0::int(len(seq_info)/15)]
                    #tick_pos  = [np.where(seq_info==i)[0][0] for i in seq_ticks]
                    #steps = len(seq_ticks)
                    #if self.scale is None:
                        #scale = max(np.abs(np.min(ddm)), np.abs(np.max(ddm)))
                    #else:
                        #scale = self.scale
                    
                    #For per residue ddm
                    ddm_residue, seq_info_unique = self.ddm_residue(ddm, seq_info)
                    
                    #write info to pickle file
                    stats = [chain.id, ddm_residue, seq_info_unique]
                    pickle.dump(stats, outname_pickle)
                    
                    
                    mask = np.zeros_like(ddm_residue, dtype=np.bool)
                    mask[np.triu_indices_from(mask)] = True
                    FINAL2 = np.ma.array(ddm_residue, mask=mask)
                    
                    tick_jump = int(np.round(len(seq_info_unique)/15,0)) #we want to add 15 seq_ticks
                    seq_ticks = seq_info_unique[0::tick_jump] #residues for which we will show tick positions
                    tick_pos = list(map(lambda x: x*tick_jump, range(len(seq_ticks))))

                    
                    if self.scale == None:
                        scale = max(np.abs(np.min(ddm_residue)), np.abs(np.max(ddm_residue)))
                    else:
                        scale = self.scale
                        
                    #if (n_rows > 1 and n_cols > 1):
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
                        
                    #elif (n_rows > 1 and n_cols == 1):
                        #img = axs[row].imshow(FINAL2, cmap=rvb, vmin=-scale, vmax=scale)
                        
                        #axs[row].set_xticks(tick_pos)
                        #axs[row].set_xticklabels(seq_ticks)
                        #plt.setp(axs[row].xaxis.get_majorticklabels(), rotation=90, fontsize='small')
                        
                        #axs[row].set_yticks(tick_pos)
                        #axs[row].set_yticklabels(seq_ticks)

                        #axs[row].set_title('Chain %s' %ID, fontsize = 'medium',fontweight="bold")
                        #axs[row].set_xlabel('Residues')
                        #axs[row].set_ylabel('Residues')
                        
                        #axs[row].spines['top'].set_visible(False)
                        #axs[row].spines['right'].set_visible(False)
                        
                        #fig.colorbar(img, ax=axs[row], fraction=0.046, pad=0.04)
                        
                    #elif (n_rows == 1 and n_cols > 1):
                        #img = axs[col].imshow(FINAL2, cmap=rvb, vmin=-scale, vmax=scale)
                        
                        #axs[col].set_xticks(tick_pos)
                        #axs[col].set_xticklabels(seq_ticks)
                        #plt.setp(axs[col].xaxis.get_majorticklabels(), rotation=90, fontsize='small')
                        
                        #axs[col].set_yticks(tick_pos)
                        #axs[col].set_yticklabels(seq_ticks)

                        #axs[col].set_title('Chain %s' %ID, fontsize = 'medium',fontweight="bold")
                        #axs[col].set_xlabel('Residues')
                        #axs[col].set_ylabel('Residues')
                        
                        #axs[col].spines['top'].set_visible(False)
                        #axs[col].spines['right'].set_visible(False)
                        
                        #fig.colorbar(img, ax=axs[col], fraction=0.046, pad=0.04)
                    
                    #else: #situation of n_rows == 1 and n_cols ==1
                        #img = axs.imshow(FINAL2, cmap=rvb, vmin=-scale, vmax=scale)
                        #axs.set_xticks(tick_pos)
                        #axs.set_xticklabels(seq_ticks)
                        #plt.setp(axs.xaxis.get_majorticklabels(), rotation=90, fontsize='small')
                        
                        #axs.set_yticks(tick_pos)
                        #axs.set_yticklabels(seq_ticks)

                        #axs.set_title('Chain %s' %ID, fontsize = 'medium',fontweight="bold")
                        #axs.set_xlabel('Residues')
                        #axs.set_ylabel('Residues')
                        
                        #axs.spines['top'].set_visible(False)
                        #axs.spines['right'].set_visible(False)
                        
                        #fig.colorbar(img, ax=axs, fraction=0.046, pad=0.04)
                        
                    old_chain_id = ID
        
        #if len(ligands)>0:            
            #ligands.sort()
            #row += 1
            #col = 0
            
        fig.tight_layout()
        plt.savefig(outname_pdf, dpi=300, transparent=True)
        plt.savefig(outname_png, dpi=300)
        plt.close()
        
        outname_pickle.close()
        
        #Need to return the png outname so that it can be easily found by the GUI. This can be much more elegant though
        return outname_png
