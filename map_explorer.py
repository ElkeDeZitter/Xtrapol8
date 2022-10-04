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


# TODO: convert to usage of iotbx instead of reading pdb as text file
def get_atom_info(pdb):
    """
    This function return a list of the coordinates and info of each atom
    coordinates: a.xyz
    residue name: i.resname
    residue number: i.resseq
    residue chain ID: i.chain_id
    residue altloc : i.altloc
    residue atomname : i.name
    rsidue serialnumber: i.i_seq+1
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


def get_d(P1, P2, axis=0):
    """
    :param P1: tuple of coordinates x1,y1,z1 (float)
    :param P2: tuple of coordinates x2,y2,z2 (float)
    :param axis: axis along which computing the sum (0 by default, should not be changed)
    :return: the distance between these two points
    """
    return np.sqrt(np.sum((P1-P2)**2, axis=axis))


def blob_detection(out, residlst, map_object, threshold, peak, coord, radius,info):
    pos = []
    neg = []
    map_object.open_map()
    #try:
    #  grid  = np.array(xplor.gridding.n, dtype=np.float32)
    #  first = np.array(xplor.gridding.first, dtype=np.float32)
    #  unit_cell = np.array(xplor.unit_cell.parameters(), dtype=np.float32)
    #  data = xplor.data.as_numpy_array()
    #except:
    #  grid = np.array(xplor.unit_cell_grid, dtype=np.float32)
    #  first = np.array(xplor.data.as_double().origin(), dtype=np.float32)
    #  unit_cell = np.array(xplor.unit_cell_parameters, dtype=np.float32)
    #  data = xplor.data.as_numpy_array()
    #map_param = (grid, first, unit_cell)
    pos.append(do_blob_search(map_object, threshold, peak, coord, radius, info))
    print_results(pos, out, residlst)
    neg.append(do_blob_search(map_object, -threshold, peak, coord, radius, info))
    print_results(neg, out, residlst)
    return pos+neg

def do_blob_search(map_object, threshold, peak, coord_pdb, radius, info):

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
        bla = (num_features_peak,) + bla
        print('Found  %i %s blobs %s peak threshold... Now integrating and allocating them.'%bla)

        for i in range(1, num_features_peak+1):
            atom = []

            feature_to_integrate = labeled_array[np.nonzero(labeled_array_peak == i)][0]
            indices = np.nonzero(labeled_array == feature_to_integrate)
            blob = map_object.data[indices]

            if np.abs(blob).max() < peak or blob.size < 2:
                continue

            sum = blob.sum()
            size = blob.size
            if threshold < 0.:
                max = blob.min()
                argmax = blob.argmin()
            else:
                max = blob.max()
                argmax = blob.argmax()
            indices = np.transpose(np.nonzero(labeled_array == feature_to_integrate))[argmax]
            coord_map = map_object.get_coord_from_map_indices(indices)

            D = get_d(coord_pdb, coord_map, axis=1)

            if coord_pdb is not None:
                if (radius is not None and D.min() <= radius) or radius is None:
                    atom.extend(info[D.argmin()])
                    print(info[D.argmin()], indices)
                    atom.extend([D.min(), max, sum, size])
            else:
                atom.extend((max, sum, size))

            if atom:
                if coord_pdb is not None and len(atom) == 10: all_atoms.append(tuple(atom))
                else: all_atoms.append(tuple(atom))

        del copy
        del copy_peak
        # Nasty but remove duplicate / Use of a decorator ?
        all_atoms = np.array(list(set(tuple(all_atoms))))
        if all_atoms.size > 0:
            sorted_indices = np.argsort(all_atoms[:, 8].astype(np.float32))
            if threshold < 0:
                return all_atoms[sorted_indices]
            else:
                all_atoms = all_atoms[sorted_indices]
                return all_atoms[::-1]

        else:
            return all_atoms
    
def check_inputs(map_name, pdb, radius, threshold, peak, log):
    """
    peak: peak_detection_threshold (sigma)
    hreshold: peak_integration_floor (sigma)
    """

    if pdb is not None:
        #pdb = pdb[0]
        if not os.path.isfile(pdb):
            print("Sorry, no such pdb file:%s"%pdb)
            sys.exit(1)
        else:
            coord, info = get_atom_info(pdb)
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
              map_object = XPLOR_Maps(map_name) \
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

def print_results(blobs, out, residlst):
    residIDlst = []
    for blob in blobs:
        out.write("%4s %4s %4s %3s %4s %4s %7s %7s %7s %6s" %("Resn", "Resv", "Chain" ,"Alt", "Atom", "ID", "d" ,"Peak" ,"Int" ,"#Voxels\n"))
        residlst.write("%4s %4s %4s %3s\n" %("Resn", "Resv", "Chain" ,"Alt"))
        print("%4s %4s %4s %3s %4s %4s %7s %7s %7s %6s" %("Resn", "Resv", "Chain" ,"Alt", "Atom", "ID", "d" ,"Peak" ,"Int" ,"#Voxels"))
        for b in blob:
            b= list(b)
            b[6] = float(b[6])
            b[7] = float(b[7])
            b[8] = float(b[8])

            out.write('%4s %4s %4s %3s %4s  %4s %7.3f %7.3f %7.3f %6s\n'%tuple(b))
            print('%4s %4s %4s %3s %4s  %4s %7.3f %7.3f %7.3f %6s'%tuple(b))
            residID = '%4s %4s %4s %3s\n' %(b[0],b[1],b[2],b[3])
            if residID not in residIDlst:
                residIDlst.append(residID)
                residlst.write(residID)
        print("")
        out.write("\n")

def map_explorer(map, pdb, radius, peak_integration_floor, peak_detection_threshold, maptype='', log=sys.stdout):
    if len(maptype)>0:
        maptype = maptype+'_'
    outname = "%speakintegration.txt" %(maptype)
    out = open(outname, 'w')
    residlst = open("%sresidlist.txt" %(maptype),'w')
    args  = check_inputs(map, pdb, radius, peak_integration_floor, peak_detection_threshold, log=log)
    blobs = blob_detection(out, residlst, *args)
    out.close()
    residlst.close()
    return outname
    
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
        print(self.a)
        self.cosa_star = (np.cos(self.beta) * np.cos(self.gamma) - np.cos(self.alpha)) / (np.sin(self.beta) * np.sin(self.gamma))
        print(self.cosa_star)
        print(self.c)
        print(self.beta)

        self.sina_star = np.sqrt(1 - self.cosa_star ** 2)
        self.transfo_matrix = np.array([[ self.a, self.b * np.cos(self.gamma), self.c * np.cos(self.beta)],
                                        [0,       self.b * np.sin(self.gamma), - self.c * self.cosa_star * np.sin(self.beta)],
                                        [0,       0,                           self.c * np.sin(self.beta)*self.sina_star]
                                        ])

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
        self.map_object = ccp4_map.map_reader(file_name=self.map_name)
        self.grid = np.array(self.map_object.unit_cell_grid, dtype=np.float32)
        self.origin = np.array(self.map_object.data.as_double().origin(), dtype=np.float32)
        self.unit_cell = np.array(self.map_object.unit_cell_parameters, dtype=np.float32)
        self.data = self.map_object.data.as_numpy_array()
        self.coord_transform_setup()


