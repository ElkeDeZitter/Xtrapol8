"""
authors and contact information
-------
Elke De Zitter - elke.de-zitter@ibs.fr
Nicolac Coquelle - nicolas.coquelle@esrf.fr
Thomas Barends - Thomas.Barends@mpimf-heidelberg.mpg.de
Jacques Philippe Colletier - jacques-Philippe.colletier@ibs.fr

"""

from __future__ import print_function
import sys
import os
import numpy as np
import iotbx.xplor.map
import argparse
from scipy import ndimage
from iotbx.pdb import hierarchy

#TODO: convert to usage of iotbx instead of reading pdb as text file

from iotbx import ccp4_map
from iotbx import file_reader

def map_data(map_obj):
    m_data = map_obj.data.as_double()
    n_real = map_obj.unit_cell_grid
    if(n_real == m_data.all()):
        return map_obj.data.as_numpy_array() / map_obj.header_rms
    else:
      # XXX hideously SLOW! MOVE TO C++
      # map_new = flex.double(flex.grid(n_real), 0)
      map_new = np.empty(n_real)
      o = m_data.origin()
      f = m_data.focus()
      for i in range(o[0],f[0]):
        for j in range(o[1],f[1]):
          for k in range(o[2],f[2]):
            map_new[i%n_real[0], j%n_real[1], k%n_real[2]] = m_data[i, j, k] / map_obj.header_rms
      return map_new

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


#Following definitions have to be changed towards iotbx.pdb.hierarchy as to read non-standard pdb files (and mmcif)
def getAtomSerialNumber(line):
    """
    integer
    """
    return int(line[6:11])

def getAtomName(line):
    return line[12:16]

def getAltLoc(line):
    return line[16]

def getResidueName(line):
    return line[17:20]

def getChainID(line):
    if line[20] == ' ':
        chainID = line[21]
    else:
        chainID = line[20]
    return chainID

def getResidueNumber(line):
    return line[22:26]


#Not used in this script
def getICode(line):
    return line[26]

def getCoordinates(line):
    """
    Each float is of format 8.3
    """
    return float(line[30:38]), float(line[38:46]), float(line[46:54])

#not used in this script
def getOccupancy(line):
    """
    f6.2
    """
    return float(line[54:60])

#Function in which it is used is not used in this script
def getBFactor(line):
    """
    f6.2
    """
    return float(line[60:66])

#Not used in this script
def getElement(line):
    return line[76:78]

#Not used in this script
def getCharge(line):
    return line[78:80]

def get_d(P1,P2, axis=0):
    """
    :param P1: tuple of coordinates x1,y1,z1 (float)
    :param P2: tuple of coordinates x2,y2,z2 (float)
    :return: the distance between these two points
    """

    return np.sqrt(np.sum((P1-P2)**2,axis=axis))

#Not used in this script
def get_B(pdb):
    lines = open(pdb).readlines()
    Bs = []
    for line in lines:
        if line.split()[0] in ['ATOM', 'HETATM']:
            Bs.append(getBFactor(line))
    return np.array(Bs)

#Not used in this script
def get_coord_via_PyMOL(pdb):
    cmd.load("%s"%pdb,"mymol")
    #cmd.select("A", "mymol and alt A")
    stored.coord = []
    stored.info = []
    res = []
    cmd.iterate_state(1, selector.process("mymol"), "stored.coord.append([x,y,z])")
    cmd.iterate_state(1, selector.process("mymol"), "stored.info.append((resn, resv, chain, alt, name, ID))")

    cmd.iterate("mymol and ( (not alt B and name CA) or (resn HOH))", "res.append((resv, resn))", space={'res': res} )
    return np.array(stored.coord), stored.info

def get_coord_from_parser(pdb):
    lines = open(pdb).readlines()
    coord = []
    info = []
    for line in lines:
        if line.split()[0] in ['ATOM','HETATM']:
            coord.append(getCoordinates(line))
            info.append((getResidueName(line),getResidueNumber(line),getChainID(line),getAltLoc(line),getAtomName(line),str(getAtomSerialNumber(line))))
    return np.array(coord), info

def get_coord_from_map_indices(indices, first, grid, cell_axes,verbose=False ):
    index = indices + first
    cartesian = index / grid * cell_axes
    if verbose:
        print('Indices %3i %3i %3i converted to x,y,z: %8.2f %8.2f %8.2f'%(tuple(indices)+tuple(cartesian)))
    return cartesian

#not used in this script
def get_indices_from_coord(coord,origin_xyz,steps_xyz, verbose=False):
    indices =  np.abs( coord - origin_xyz ) / steps_xyz
    return np.rint(indices).astype(np.int32)

def blob_detection(out, residlst, xplor, threshold, peak, coord, radius,info):
    pos = []
    neg = []
    try:
      grid  = np.array(xplor.gridding.n, dtype=np.float32)
      first = np.array(xplor.gridding.first, dtype=np.float32)
      unit_cell = np.array(xplor.unit_cell.parameters(), dtype=np.float32)
      data = xplor.data.as_numpy_array()  
    except:
      grid = np.array(xplor.unit_cell_grid, dtype=np.float32)
      first = np.array(xplor.data.as_double().origin(), dtype=np.float32)
      unit_cell = np.array(xplor.unit_cell_parameters, dtype=np.float32) 
      data = map_data(xplor)
    map_param = (grid, first, unit_cell)
    pos.append(do_blob_search(data, threshold, peak, map_param, coord, radius,info))
    print_results(pos, out, residlst)
    neg.append(do_blob_search(data, -threshold, peak, map_param, coord, radius,info))
    print_results(neg, out, residlst)
    return pos+neg

def do_blob_search(data, threshold, peak, map_param, coord_pdb, radius, info):

        grid, first, unit_cell = map_param
        s=ndimage.morphology.generate_binary_structure(3,3)
        if threshold < 0:
            copy = np.where(data <= threshold, 1, 0)
            copy_peak = np.where(data <= -peak, 1, 0)
            bla = ('negative','below')
        else:
            copy = np.where(data >= threshold, 1, 0)
            copy_peak = np.where(data >= peak, 1, 0)
            bla = ('positive','above')

        labeled_array, num_features = ndimage.label(copy, structure=s)
        labeled_array_peak, num_features_peak = ndimage.label(copy_peak, structure=s)

        all_atoms = []
        bla = (num_features_peak,) + bla
        print('Found  %i %s blobs %s peak threshold... Now integrating and allocating them.'%bla)

        for i in range(1, num_features_peak+1):
            atom = []

            feature_to_integrate = labeled_array[np.nonzero(labeled_array_peak == i)][0]
            blob = data[np.nonzero(labeled_array == feature_to_integrate)]

            if np.abs(blob).max() < peak or blob.size < 2 : continue

            sum = blob.sum()
            size = blob.size
            if threshold < 0.:
                max = blob.min()
                argmax = blob.argmin()
            else:
                max = blob.max()
                argmax = blob.argmax()
            indices = np.transpose(np.nonzero(labeled_array == feature_to_integrate))[argmax]
            coord_map = get_coord_from_map_indices(indices, first, grid,unit_cell[0:3])

            D = get_d(coord_pdb, coord_map, axis=1)

            if coord_pdb is not None:
                if (radius is not None and D.min() <= radius) or radius is None:
                    atom.extend(info[D.argmin()])
                    print(info[D.argmin()])
                    #cmd.iterate('mymol and ID %i'%D.argmin(), 'atom.extend((resn, resv, chain, alt, name, ID))', space = {'D': D, 'atom': atom})
                    atom.extend([D.min(), max, sum, size])
                    #print max, coord_map


            else: atom.extend((max, sum, size))
            if atom:
                if coord_pdb is not None and len(atom) == 10: all_atoms.append(tuple(atom))
                else: all_atoms.append(tuple(atom))

        del copy
        del copy_peak
        # Nasty but remove duplicate / Use of a decorator ?
        all_atoms = np.array(list(set(tuple(all_atoms))))
        #all_atoms = np.array(all_atoms)
        if all_atoms.size > 0:
            sorted_indices = np.argsort(all_atoms[:,8].astype(np.float32))
            if threshold < 0:
                return all_atoms[sorted_indices]
            else:
                all_atoms = all_atoms[sorted_indices]
                return all_atoms[::-1]

        else: return all_atoms
    
def check_inputs(map, pdb, radius, threshold, peak, log):

    if pdb is not None:
        #pdb = pdb[0]
        if not os.path.isfile(pdb):
            print("Sorry, no such pdb file:%s"%pdb)
            sys.exit(1)
        else:
            coord, info = get_atom_info(pdb)
    else: coord = None

    if not os.path.isfile(map):
        print("Sorry, no such map file:%s"%(map),file=log)
        print("Sorry, no such map file:%s"%(map))
        sys.exit(1)
    else:
        try:
            xplor = iotbx.xplor.map.reader(file_name=map)
        except:
            try:
              xplor = ccp4_map.map_reader(file_name=map)
            except:
              print("Sorry, %s map is not a valid XPLOR or CCP4 map. Aborting map explorer."%(map), file=log)
              print("Sorry, %s map is not a valid XPLOR or CCP4 map. Aborting map explorer."%(map))
              sys.exit(1)

    if radius is not None and coord is None:
        print('You cannot provide a radius without a valid pdb file', file=log)
        print('You cannot provide a radius without a valid pdb file')
        radius = None

    if radius is not None:
        radius = float(radius)
        #try: radius = float(radius[0])
        #except TypeError: radius = float(radius)
    return xplor, threshold, peak, coord, radius, info

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

#Modified by Elke to be called from Fextr.py and get rid of argparser:
#if __name__ == '__main__':
def map_explorer(map, pdb, radius, threshold, peak, maptype='', log=sys.stdout):
    if len(maptype)>0:
        maptype = maptype+'_'
    outname = "%speakintegration.txt" %(maptype)
    out = open(outname, 'w')
    residlst = open("%sresidlist.txt" %(maptype),'w')
    args  = check_inputs(map, pdb, radius, threshold, peak, log=log)
    blobs = blob_detection(out, residlst, *args)
    out.close()
    residlst.close()
    return outname
    
    #parser = ArgumentParser()
    #args= check_inputs(parser)
    #results = np.zeros((10,10))    
    #i = 20
    #j = 10
    #radius = 1.5
    #peak = 4.0
    #integration = 3.0
    #pdb = 'INT_int%i_on%i.pdb'%(i,j)
    #coord, res, info = get_coord_via_PyMOL(pdb)
    #map = 'LCLS_399_all_int%i_on%i_Fcalc.map'%(i,j)
    #xplor = iotbx.xplor.map.reader(file_name=map)
    #args = (xplor, integration, peak, coord, radius, info)
    #blobs = blob_detection(*args)
    #pos, neg = blobs
    #pos1 = pos[pos[:,0] == 'PIA']
    #neg1 = neg[neg[:,0] == 'PIA']
    #for alt in ['A', 'B', 'C', 'D']:
    #    print 'Alt %s'%alt
    #    print pos1[pos1[:,3] == alt][:,8].astype(np.float32).sum()    
    #    print neg1[neg1[:,3] == alt][:,8].astype(np.float32).sum()
    #if pos1.size == 0: print 'Sum = 0'
    #else:
    #s = 0.
    #s += pos[pos[:,0] == 'PIA'][:,8].astype(np.float32).sum()
    #s += np.abs(neg[neg[:,0] == 'PIA'][:,8].astype(np.float32).sum())
    #print '%5i %5i %8.2f'%(i,j,s)
    #results[i/5,j/5] = s
    #np.save('3ps_FoFc_LI56.npy',results)        
