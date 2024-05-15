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
import sys
import os
import numpy as np
import iotbx.xplor.map
from iotbx import ccp4_map
from scipy import ndimage
from iotbx.pdb import hierarchy
import scipy.stats

# TODO: convert to usage of iotbx instead of reading pdb as text file
def get_atom_info_all(pdb):
    """
    This function return a list of the coordinates and info of each atom
    coordinates: a.xyz
    residue name: i.resname
    residue number: i.resseq
    residue chain ID: i.chain_id
    residue altloc : i.altloc
    residue atomname : i.name
    residue serialnumber: i.i_seq+1
    """
    
    pdb_hier = hierarchy.input(file_name=pdb)
    hier = pdb_hier.hierarchy
    
    coord = []
    info  = []
    
    for chain in hier.chains():
        for res_group in chain.residue_groups():
            for atom_group in res_group.atom_groups():
                for a in res_group.atoms():
                    coord.append(a.xyz)
                    i = a.fetch_labels()
                    info.append((i.resname, i.resseq, i.chain_id, i.altloc, i.name, i.i_seq+1))
                    
    return np.array(coord), info

def get_atom_info_residlist(pdb, resid_list_in):
    """
    This function return a list of the coordinates and info of each atom
    coordinates: a.xyz
    residue name: i.resname
    residue number: i.resseq
    residue chain ID: i.chain_id
    residue altloc : i.altloc
    residue atomname : i.name
    residue serialnumber: i.i_seq+1
    but only for the residues specified in the resid_list_in
    """
    
    with open(resid_list_in) as r:
        residlst = r.read().split("\n")
    residlst = [lne for lne in residlst if "Resn" not in lne] #remove header lines
    residlst = [lne.split() for lne in residlst] #make list of list
    
    residlst_corr = []
    for lne in residlst:
        #append an empty item if no altloc provided
        #remove from residlist if not 3 or 4 items in lne
        if len(lne) == 4:
            residlst_corr.append(lne)
        elif len(lne) == 3:
            residlst_corr.append(lne+[''])

    pdb_hier = hierarchy.input(file_name=pdb)
    hier = pdb_hier.hierarchy
    
    #chains = [chain.id for chain in hier.chains()]
    
    coord = []
    info  = []
    
    for chain in hier.chains():
        for res_group in chain.residue_groups():
            for atom_group in res_group.atom_groups():
                for a in res_group.atoms():
                    i = a.fetch_labels()
                    if [i.resname.lstrip().rstrip(), i.resseq.lstrip().rstrip(), i.chain_id.lstrip().rstrip(), i.altloc.lstrip().rstrip()] in residlst_corr:
                        coord.append(a.xyz)
                        info.append((i.resname, i.resseq, i.chain_id, i.altloc, i.name, i.i_seq+1))
                        
                    
    return np.array(coord), info


def get_d(P1, P2, axis=0):
    """
    :param P1: tuple of coordinates x1,y1,z1 (float)
    :param P2: tuple of coordinates x2,y2,z2 (float)
    :param axis: axis along which computing the sum (0 by default, should not be changed)
    :return: the distance between these two points
    """
    return np.sqrt(np.sum((P1-P2)**2, axis=axis))


def blob_detection(map_object, threshold, peak, coord, radius, info, out=sys.stdout, residlst_out = None):
    #neg = []
    map_object.open_map()
    arr_pos = do_blob_search(map_object, threshold, peak, coord, radius, info)
    print_results([arr_pos], out, residlst_out)
    arr_neg = do_blob_search(map_object, -threshold, peak, coord, radius, info)
    print_results([arr_neg], out, residlst_out)
    # print(arr_pos.shape, arr_neg.shape)
    # print(arr_pos, arr_neg)
    if arr_pos.shape == (0,):
        arr_final = arr_neg
    elif arr_neg.shape == (0,):
        arr_final = arr_pos
    else:
        arr_final = np.concatenate((arr_pos, arr_neg), axis=0)
    return arr_final

def do_blob_search(map_object, threshold, peak, coord_pdb, radius, info):
    """
    1) Find peaks >= the peak (peak_detection_threshold) and asing blobs to all vaxels around the peak that have a value >= threshold (peak_integration_floor)
    2) Asign a unique number to each of voxels of each blob
    3) Only keep those blobs which have the peak voxel wihtin the radius from a protein atom
    4) Asign the blobs to the atoms and store the atom associated peak/blob information
    """

    s = ndimage.morphology.generate_binary_structure(3,3)
    if threshold < 0:
        copy = np.where(map_object.data <= threshold, 1, 0)
        copy_peak = np.where(map_object.data <= -peak, 1, 0)
        bla = ('negative', 'below')
    else:
        copy = np.where(map_object.data >= threshold, 1, 0)
        copy_peak = np.where(map_object.data >= peak, 1, 0)
        bla = ('positive', 'above')

    labeled_array, num_features = ndimage.label(copy, structure=s)
    labeled_array_peak, num_features_peak = ndimage.label(copy_peak, structure=s)

    all_atoms = []
    mask_indices = []
    blobs = []
    bla = (num_features_peak,) + bla
    print('Found  %i %s blobs %s peak threshold... Now integrating and allocating them.'%bla)

    for i in range(1, num_features_peak+1):
        atom = []

        feature_to_integrate = labeled_array[np.nonzero(labeled_array_peak == i)][0]
        indices = np.nonzero(labeled_array == feature_to_integrate)
        blob = map_object.data[indices]
        #print(any([np.array_equal(blob, y) for y in blobs]))
        if np.abs(blob).max() < peak or blob.size < 2 or any([np.array_equal(blob, y) for y in blobs]):
            continue
        blobs.append(blob)
        sm = blob.sum()
        size = blob.size
        if threshold < 0.:
            mx = blob.min()
            argmax = blob.argmin()
        else:
            mx = blob.max()
            argmax = blob.argmax()
        index_max = np.transpose(np.nonzero(labeled_array == feature_to_integrate))[argmax]
        coord_map = map_object.get_coord_from_map_indices(index_max)

        D = get_d(coord_pdb, coord_map, axis=1)

        if coord_pdb is not None: #if no pdb is provided (should not happen within Xtrapol8)
            if (radius is not None and D.min() <= radius) or radius is None: #only select blobs which are close enough to atoms
                atom.extend(info[D.argmin()])
                atom.extend([D.min(), mx, sm, size, indices])
        else:
            atom.extend((mx, sm, size))

        if atom: #should normally always be the case

            #if coord_pdb is not None and len(atom) == 10 and atom not in all_atoms: all_atoms.append(tuple(atom))
            #else: all_atoms.append(tuple(atom))
            all_atoms.append(tuple(atom))

    del copy
    del copy_peak
    #print("done")
    # Nasty but remove duplicate / Use of a decorator ?
    #arr_all_atoms = np.array(list(set(tuple(all_atoms))))
    all_atoms = np.array(all_atoms)

    if all_atoms.size > 0:
        sorted_indices = np.argsort(all_atoms[:, 8].astype(np.float32))
        if threshold < 0:
            return all_atoms[sorted_indices]
        else:
            all_atoms = all_atoms[sorted_indices]
            return all_atoms[::-1]

    else:
        return all_atoms
    
def check_inputs(map_name, pdb, radius, threshold, peak, resid_list_in=None, log=sys.stdout):
    """
    peak: peak_detection_threshold (sigma)
    hreshold: peak_integration_floor (sigma)
    """

    if pdb is not None:
        #pdb = pdb[0]
        if not os.path.isfile(pdb):
            print("Sorry, no such pdb file:%s"%pdb)
            sys.exit(1)
        elif resid_list_in!=None:
            if not os.path.isfile(resid_list_in):
                print("Sorry, no such residue list file:%s"%pdb)
                sys.exit(1)
            else:
                coord, info = get_atom_info_residlist(pdb, resid_list_in)
        else:
            coord, info = get_atom_info_all(pdb)
    else: coord = None

    if not os.path.isfile(map_name):
        print("Sorry, no such map file:%s" % map_name, file=log)
        print("Sorry, no such map file:%s" % map_name)
        sys.exit(1)
    else:
        if map_name.endswith('ccp4'):
            map_object = CCP4_Maps(map_name)
            #xplor = iotbx.xplor.map.reader(file_name=map)
        else:
            try:
              map_object = XPLOR_Maps(map_name)
              #    = ccp4_map.map_reader(file_name=map)
            except:
              print("Sorry, %s map is not a valid XPLOR or CCP4 map. Aborting map explorer." % map, file=log)
              print("Sorry, %s map is not a valid XPLOR or CCP4 map. Aborting map explorer." % map)
              sys.exit(1)

    if radius is not None and coord is None:
        print('You cannot provide a radius without a valid pdb file', file=log)
        print('You cannot provide a radius without a valid pdb file')
        radius = None

    if radius is not None:
        radius = float(radius)

    return map_object, threshold, peak, coord, radius, info

def print_results(blobs, out=sys.stdout, residlst_out=None):
    residIDlst = []
    for blob in blobs:
        out.write("%4s %4s %4s %3s %4s %4s %7s %7s %8s %6s" %("Resn", "Resv", "Chain" ,"Alt", "Atom", "ID", "d" ,"Peak" ,"Int" ,"#Voxels\n"))
        if residlst_out != None:
            residlst_out.write("%4s %4s %4s %3s\n" %("Resn", "Resv", "Chain" ,"Alt"))
        print("%4s %4s %4s %3s %4s %4s %7s %7s %7s %6s" %("Resn", "Resv", "Chain" ,"Alt", "Atom", "ID", "d" ,"Peak" ,"Int" ,"#Voxels"))
        for b in blob:
            b= list(b)
            b[6] = float(b[6])
            b[7] = float(b[7])
            b[8] = float(b[8])

            out.write('%4s %4s %4s %3s %4s  %4s %7.3f %7.3f %8.2f %6s\n'%tuple(b[:-1]))
            print('%4s %4s %4s %3s %4s  %4s %7.3f %7.3f %8.2f %6s'%tuple(b[:-1]))
            residID = '%4s %4s %4s %3s\n' %(b[0], b[1], b[2], b[3])
            if residID not in residIDlst and residlst_out != None:
                residIDlst.append(residID)
                residlst_out.write(residID)
        print("")
        #out.write("\n")

def map_explorer(map_name, pdb, radius, peak_integration_floor, peak_detection_threshold, zscore= 2, resid_list_in=None, log=sys.stdout):
    
    args  = check_inputs(map_name, pdb, radius, peak_integration_floor, peak_detection_threshold, resid_list_in=resid_list_in, log=log)
    map_object, threshold, peak, coord, radius, info = args
    outname = "peakintegration.txt"
    out = open(outname, 'w')
        
    if resid_list_in == None:
        print("Zscore: %4.2f" % zscore)
        residlst_out = open("residlist.txt" ,'w')
        blobs = blob_detection(map_object, threshold, peak, coord, radius, info, out, residlst_out=residlst_out)
        residlst_out.close()
    else:
        blobs = blob_detection(map_object, threshold, peak, coord, radius, info, out)
    out.close()
    
    if len(blobs) >=1 :
        peaks = blobs[:, 8].astype(np.float32)
    else:
        mask = np.ndarray([0,0])
        return outname, outname, mask, [0, 0, 0]

    if resid_list_in == None:
        peaks_corr = np.where(peaks < 0, peaks - (np.max(np.where(peaks > 0, -np.infty, peaks))), peaks - (np.min(np.where(peaks < 0, np.infty,peaks))))
        idx = np.where(scipy.stats.zscore(np.abs(peaks_corr)) > zscore)
        if idx[0].shape[0] >= 1:
            mask = blobs[idx, -1] #Keep the coordinates of the blob rather than storing the whole mask 
        else:
            #zscoring did not work, no blobs retained, use all blobs
            mask = blobs[:, -1]
            shape1 = mask.shape[0]
            mask = mask.reshape((1,shape1))

        outname_zscore = "peakintegration_Zscore%4.2f.txt" % (zscore)
        out_zscore = open(outname_zscore, 'w')
        residlst_zscore = open("residlist_Zscore%4.2f.txt" % (zscore), 'w')

        print("Peaks kept after Z-scoring")
        print_results(blobs[idx, :], out= out_zscore, residlst_out=residlst_zscore)
        out_zscore.close()
        residlst_zscore.close()
        
        pos = 0
        neg = 0
        for peak in blobs[idx, 8][0]:

            if peak > 0: pos += peak
            else:
                neg -= peak

    else:
        #No Z-scoring, make mask using the blobs around the residues in the resid list
        print("Residue list is provided, Z-scoring will not be applied. All peaks around the specified residue(s) will be included in the analysis")
        mask = blobs[:, -1]
        shape1 = mask.shape[0]
        mask = mask.reshape((1,shape1)) #to make compatible with the mask generated based on Z-scoring
        outname_zscore = None
        
        pos = 0
        neg = 0
        
        for i in range(shape1): #can use shape 1 because all blobs are used in the mask
            peak = blobs[i, 8] #use integrated value, not peak height
            if peak > 0: pos += peak
            else:
                neg -= peak
            
    return outname, outname_zscore, mask, [pos, neg, pos+neg]
    
class Maps(object):

    def __init__(self, map_name):
        self.map_name = map_name
        # map origin
        self.origin = [0, 0, 0]
        # grid is the total number of points along each unit cell axis
        self.grid = [0, 0, 0]
        # unit cell parameters
        self.alpha = 90.
        self.beta = 90.
        self.gamma = 90.
        self.a = 0.
        self.b = 0.
        self.c = 0.

    def open_map(self):
        """
        Overwritten in sub-classes
        """
        return

    def coord_transform_setup(self):
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = self.unit_cell
        self.alpha = np.deg2rad(self.alpha)
        self.beta = np.deg2rad(self.beta)
        self.gamma = np.deg2rad(self.gamma)
        self.cosa_star = (np.cos(self.beta) * np.cos(self.gamma) - np.cos(self.alpha)) / (np.sin(self.beta) * np.sin(self.gamma))

        self.sina_star = np.sqrt(1 - self.cosa_star ** 2)
        self.transfo_matrix = np.array([[ self.a, self.b * np.cos(self.gamma), self.c * np.cos(self.beta)],
                                        [0,       self.b * np.sin(self.gamma), - self.c * self.cosa_star * np.sin(self.beta)],
                                        [0,       0,                           self.c * np.sin(self.beta)*self.sina_star]
                                        ]).T

    def get_coord_from_map_indices(self, indices, verbose=False):
        """
        Will return cartesian coordinates from grid map indices
        first is the map origin (int array (3,))
        grid is the number of map points covering a unit cell (int array (3,))
        cell_axes is the unit cell parameters a, b and c (float array (3,))
        WARNING: Does not work properl
        """
        index = indices + self.origin
        cartesian = np.matmul(index / self.grid, self.transfo_matrix)
        if verbose:
            print('Indices %3i %3i %3i converted to x,y,z: %8.2f %8.2f %8.2f' % (tuple(indices) + tuple(cartesian)))
        return cartesian

    # not used in this script
    def get_indices_from_coord(self, coord, origin_xyz, steps_xyz, verbose=False):
        indices = np.abs(coord - origin_xyz) / steps_xyz
        return np.rint(indices).astype(np.int32)


class CCP4_Maps(Maps):

    def __init__(self, map_name):
        super(CCP4_Maps, self).__init__(map_name)

    def open_map(self):
        self.map_object = ccp4_map.map_reader(file_name=self.map_name)
        self.grid = np.array(self.map_object.unit_cell_grid, dtype=np.float32)
        self.origin = np.array(self.map_object.data.as_double().origin(), dtype=np.float32)
        self.unit_cell = np.array(self.map_object.unit_cell_parameters, dtype=np.float32)
        self.data = self.map_object.data.as_numpy_array()
        self.coord_transform_setup()


class XPLOR_Maps(Maps):

    def __init__(self, map_name):
        super(XPLOR_Maps, self).__init__(map_name)

    def open_map(self):
        print(self.map_name)
        self.map_object = iotbx.xplor.map.reader(file_name=self.map_name)
        self.grid = np.array(self.map_object.gridding.n, dtype=np.float32)
        self.origin = np.array(self.map_object.gridding.first, dtype=np.float32)
        self.unit_cell = np.array(self.map_object.unit_cell.parameters(), dtype=np.float32)
        self.data = self.map_object.data.as_numpy_array()
        self.coord_transform_setup()


