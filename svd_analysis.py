# coding=utf-8
from __future__ import print_function
import re
import struct
import iotbx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.ndimage
import scipy
import scipy.stats
import scipy.signal
from libtbx.utils import Usage
from Fextr_utils import check_file_existance
from libtbx import adopt_init_args
import pickle
from iotbx import ccp4_map
import iotbx.map_tools

import glob
import sys
import os
#######################READCCP4MAP################
#
def readccp4map(mapinname):
    """
    Read a ccp4 map file
    Input:  mapinname - name of the ccp4 map file
    Output: size - number of columns, rows, slices
            start - starting point in columns, rows, slices
            intervals - number of intervals in columns, rows, slices
            uc - unit cell parameters a, b, c, alpha, beta, gamma
            order - order of columns, rows, slices
            skew - the skew matrix
            skew_trn - the skew translation vector
            symops - the symmetry operators
            mrc_start - start of the map in the original mrc file
            maparray - the 3D array of the map values
    This function could probably be replaced by a cctbx function:
        from iotbx import ccp4_map
        ccp4_map_obj = ccp4_map.map_reader(file_name=mapinname)
    struct.unpack could probably be replaced by np.frombuffer (should be faster and more elegant)
    """
    # with open(mapinname, "rb") as mapinfile:
    #     metadata=mapinfile.read(208)
    mapinfile=open(mapinname, "rb")  # Open the map file in binary read mode
    data_in=mapinfile.read(208) # Read the first 208 bytes (the header section)
    print("data_in attribute",data_in.__getattribute__)
    # mapinfile.close() # Closing as to correctly re-open later. This seem to be required to read in correctly data further.
    mapinfo=struct.unpack("10i 6f 3i 3f 3i 12f 15f", data_in) # Unpack the header
    n_col=mapinfo[0] #Number of columns (fastest changing in map array)
    n_row=mapinfo[1] #Number of rows (second fastest changing in map array)
    n_sli=mapinfo[2] #Number of slices (slowest changing in map array)
    size=[n_col,n_row,n_sli] #ccp4_map_obj.data.all()
    # MODE=mapinfo[3] #Data type (2 = float, 0 = int8, 1 = int16) #Not used
    c_start=mapinfo[4] #starting point of the map in columns
    r_start=mapinfo[5] #starting point of the map in rows
    s_start=mapinfo[6] #starting point of the map in slices
    start=[c_start,r_start,s_start] #ccp4_map_obj.get_origin() but different order
    n_x=mapinfo[7] #Number of intervals along X (columns)
    n_y=mapinfo[8] #Number of intervals along Y (rows)
    n_z=mapinfo[9] #Number of intervals along Z (slices)
    intervals=[n_x,n_y,n_z]
    # uc_a=mapinfo[10] #Unit cell parameter a
    # uc_b=mapinfo[11] #Unit cell parameter b
    # uc_c=mapinfo[12] #Unit cell parameter c
    # uc_alpha=mapinfo[13] #Unit cell parameter alpha
    # uc_beta=mapinfo[14] #Unit cell parameter beta
    # uc_gamma=mapinfo[15] #
    # CELL=np.asarray([uc_a,uc_b,uc_c,uc_alpha,uc_beta,uc_gamma])
    uc = np.asarray(mapinfo[10:16]) #ccp4_map_obj.unit_cell_parameters
    map_c=mapinfo[16] #Which axis corresponds to columns
    map_r=mapinfo[17] #Which axis corresponds to rows
    map_s=mapinfo[18] #Which axis corresponds to slices
    order=[map_c,map_r,map_s]
    #print 'Order of columns, rows, slices:',order
    # a_min=mapinfo[19] #Minimum density value #Not used #ccp4_map_obj.header_min
    # a_max=mapinfo[20] #Maximum density value #Not used #ccp4_map_obj.header_max
    # a_mean=mapinfo[21] #Mean density value #Not used #ccp4_map_obj.header_mean
    # i_sg=mapinfo[22] #Space group number (0 or 1 if unknown) #Not used
    n_symbt=mapinfo[23] #Number of bytes used for storing symmetry operators
    # skew_flag=mapinfo[24] #Flag for skew matrix (0 if not present, 1 if present) #Not used
    skew=[[mapinfo[25],mapinfo[26],mapinfo[27]],
          [mapinfo[28],mapinfo[29],mapinfo[30]],
          [mapinfo[31],mapinfo[32],mapinfo[33]]] # The skew matrix
    skew=np.asarray(skew) 
    skew_trn=mapinfo[34:37] # The skew translation vector
    skew_trn=np.asarray(skew_trn)
    # FUTURE=mapinfo[37:52] #Future expansion space #Not used
    mrc_start=[mapinfo[49],mapinfo[50],mapinfo[51]] #start of the map in the original mrc file

    data_in=mapinfile.read(8) #Read the word 'MAP ' and the 4-byte machine stamp
    mapinfo=struct.unpack("fi",data_in)
    a_rms=mapinfo[0] #RMS deviation of map #Not used #ccp4_map_obj.header_rms
    NLABL=mapinfo[1] #Meaning? #Not used
    
    LABEL=mapinfile.read(800) #Meaning? #Not used
    # with open(mapinname, "rb") as mapinfile:
    symops=mapinfile.read(n_symbt)

    n_voxels=n_col*n_row*n_sli
    formatstring=str(n_voxels)+'f'
    # with open(mapinname, "rb") as mapinfile:
    data_in=mapinfile.read(n_voxels*4)
    print("data_in attribute",data_in.__getattribute__)
    maparray=struct.unpack(formatstring,data_in) #ccp4_map_obj.data or ccp4_map_obj.map_data()
    mapinfile.close()
    maparray=np.asarray(maparray)
    maparray=np.reshape(maparray, (n_sli,n_row,n_col))

    return size, start, intervals, uc, order, skew, skew_trn, symops, mrc_start, maparray

class CCP4_Maps(object):

    def __init__(self, map_name):
        self.map_name = map_name

    def open_map(self):
        self.map_object = ccp4_map.map_reader(file_name=self.map_name)
        self.grid = self.map_object.unit_cell_grid
        self.origin = np.array(self.map_object.data.as_double().origin(), dtype=np.float32)
        self.unit_cell = self.map_object.unit_cell_parameters
        self.data = self.map_object.data.as_numpy_array()

# def skew(a,b,c,alpha,beta,gamma):
#     """
#     Calculate the skew or scale matrix from unit cell parameters
#     Input:  a,b,c - unit cell lengths
#             alpha,beta,gamma - unit cell angles in degrees
#     Output: s - the 3x3 skew matrix
#     Seems like it is not userd an_ywhere
#     """
#     pi=3.14159265359
#     alpha=pi*alpha/180.0
#     beta=pi*beta/180.0
#     gamma=pi*gamma/180.0

#     vsq=1.0 - (np.cos(alpha)*np.cos(alpha)) - (np.cos(beta)*np.cos(beta)) - (np.cos(gamma)*np.cos(gamma)) - (2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
#     v=np.sqrt(vsq) # volume of unit parallelepiped

    
#     s=np.zeros((3,3))
#     s[0,0]=1.0/a
#     s[0,1]=-1.0*np.cos(gamma)/(a*np.sin(gamma))
#     s[0,2]=(np.cos(alpha)*np.cos(gamma)-np.cos(beta)) / (a*v*np.sin(gamma))
#     s[1,1]=1.0/(b*np.sin(gamma))
#     s[1,2]=(np.cos(beta)*np.cos(gamma)-np.cos(alpha)) / (b*v*np.sin(gamma))
#     s[2,2]=np.sin(gamma)/(c*v)
#     return s

#
def writeccp4map_full(mapoutname, size, start, intervals, uc, order, skew, skew_trn, symops, maparray):
    """
    Write a ccp4 map file
    Input:  mapoutname - name of the ccp4 map file to be written
            size - number of columns, rows, slices
            start - starting point in columns, rows, slices
            intervals - number of intervals in columns, rows, slices
            uc - unit cell parameters a, b, c, alpha, beta, gamma
            order - order of columns, rows, slices
            skew - the skew matrix
            skew_trn - the skew translation vector
            symops - the symmetry operators # Not used
            maparray - the 3D array of the map values
    Output: a ccp4 map file
    This function could probably be replaced by a cctbx function
    """
    header=bytes()
    for value in size:
        header += struct.pack('i',value) #Number of columns, rows, slices
    header += struct.pack('i',int(2)) #Set the mode to 2 (i.e. float)
    for value in start:
        header += struct.pack('i',value) #Starting point in columns, rows, slices
    for value in intervals:
        header += struct.pack('i',value) #Number of intervals in columns, rows, slices
    for value in uc:
        header += struct.pack('f',value) #Unit cell parameters a, b, c, alpha, beta, gamma
    for value in order:
        header += struct.pack('i',value)
    a_min=np.min(maparray) #Minimum density value
    a_max=np.max(maparray) #Maximum density value
    a_mean=np.average(maparray) #Mean density value
    header += struct.pack('3f',a_min,a_max,a_mean)
    i_sg=1 #Space group number (0 or 1 if unknown)
    n_symbt=0
    #skew_flag=1 # IF skew IS TO BE WRITTEN OUT 
    skew_flag=0 # IF NO skew E.G. FOR PYMOL
    header += struct.pack('3i',i_sg,n_symbt,skew_flag)
    #skew=np.reshape(skew,(9)) # IF skew IS TO BE WRITTEN OUT
    skew=np.zeros((9)) # IF NO skew E.G. FOR PYMOL
    skew_trn=np.zeros((3))
    for value in skew:
        header += struct.pack('f',value) #The skew matrix
    for value in skew_trn:
        header += struct.pack('f',value) #The skew translation vector
    for n in range(15):
        header += struct.pack('i',int(0)) #Future expansion space: add zeros
    header += struct.pack('4c',b'M',b'A',b'P',b' ') #Adding the word 'MAP '
    #header += struct.pack('4c',chr(0x44), chr(0x41), chr(0x00), chr(0x00))
    header += struct.pack('4c',b'D',b'A',b'0',b'0') #Adding the little-endian machine stamp 'D','A',0,0
    a_rms=np.std(maparray)
    header += struct.pack('f',a_rms) #RMS deviation of map
    header += struct.pack('i',3) 
    for n in range(800):
        header += struct.pack('c',b' ') #Adding 800 bytes of blank space instead of LABEL

    mapoutfile=open(mapoutname, "wb")
    mapoutfile.write(header)
    mapsize=np.shape(maparray)
    maplength=(mapsize[0]*mapsize[1]*mapsize[2])
    maparray=np.reshape(maparray,(maplength))
    for n in range(maplength):
        voxelout=struct.pack('f',maparray[n])
        mapoutfile.write(voxelout)
    mapoutfile.close()


class SVD_analysis(object):
    def __init__(self,
                    map_2mFextr_DFc_list = [],
                    occupancies=[],
                    xray_structure = None,
                    numvec = 5, #Need to implement function that automatically detects the number of useful number of vectors
                    prefix = '',
                    log = sys.stdout):
        
        #Try this instead of  repetition of the arguments
        adopt_init_args(self, locals())
        assert len(map_2mFextr_DFc_list) == len(occupancies), 'the number of difference maps and occupancies is not equal, this will be nonsense and end with an error somewhere'
        
        #Sort the map files and occupancies in order to have occupancies from small to large, probably just for cosmethics
        #This might not work in pyhton3
        zipped = zip(occupancies, map_2mFextr_DFc_list)
        zipped_sorted = sorted(zipped, key = lambda x:x[0])
        occupancies, map_2mFextr_DFc_list =zip(*zipped_sorted)
        self.occupancies     = list(occupancies)
        self.map_2mFextr_DFc_list = list(map_2mFextr_DFc_list)
        
    # def check_orthonormality(self, matrix):
    #     """
    #     Check if matrix is orthonormal
    #     """
    #     return np.allclose(np.dot(matrix,matrix.T), np.eye(matrix.shape[0]))
    
    def check_orthonormality(self, matrix, tol=1e-3):
        """
        Check if matrix is orthonormal (should be faster than the implementation above).
        The set tolerance should be evaluated
        """
        prod = np.dot(matrix, matrix.T)
        return np.linalg.norm(prod - np.eye(prod.shape[0]), ord='fro') < tol

    def run_svd(self, dataset):
        """
        Run SVD to obtain the follinf vectors:
        u: Unitary matrix having left singular vectors as columns
        s: The singular values, sorted in non-increasing order
        vh: Unitary matrix having right singular vectors as rows
        """
        print('Now performing Singular Value Decomposition...')
        #u,s,v=np.linalg.svd(dataset,full_matrices=False)
        err = 0
        try:
            u, s, vh =scipy.linalg.svd(dataset, full_matrices=False)
        except LinAlgError:
            print('SVD computation did not converge')
            err += 1
        
        #Check SVD
        if u.shape != dataset.shape: #left singular values #shape: (n_voxels, n_maps)
            print("u shape not ok")
            err += 1
        elif s.shape != (dataset.shape[1],): #diagonal matrix with singular values #shape: (min(n_voxels, n_maps),), so here (n_maps,)
            print("s shape not ok")
            err += 1
        elif vh.shape != (dataset.shape[1],dataset.shape[1]): #right singular values #shape: (n_maps, n_maps)
            print("vh shape not ok")
            err += 1
        # elif self.check_orthonormality(u) == False: #This takes a lot of memory
        #     print("u not orthonormal")
        #     err += 1
        elif self.check_orthonormality(vh.T) == False:
            print("vh not orthonormal")
            err += 1
        
        if err >= 1:
            print("SVD failed! Results will be nonsense")
            success = False
        else:
            success = True
            
        return u, s, vh, success
        
    def plot_singular_values(self,s):
        """
        Plot the singular values from the diagonal of s
        """    
        plt.bar(np.arange(len(self.map_2mFextr_DFc_list)),s)
        plt.title('Singular values for w1.00')
        axes = plt.gca()
        axes.set_xlim([-2,12])
        axes.set_ylim([0,5000])
        #plt.show()
        #outfig1=outdir+'Singular_values.png'
        plt.savefig('singular_values_w1.00.png', dpi=300)
        plt.close()
    
    def plot_right_singular_values(self, vh):
        
        #for v in range(numvec):
        #    leg="vector #"+str(v)
        #    plt.semilogx((times[:]),vh[v,:],label=leg)
        #    print('****** Right Singular Vector #'+str(v)+' ******')
        #    print(vh[v,:])

        #plt.title('Time dependence (right singular vector) #'+str(v))
        #plt.xlabel('Time')
        #plt.ylabel('Amplitude')
        #plt.legend()
        #splt.show()

        colorlist=['xkcd:purple','xkcd:royal blue','xkcd:blue','xkcd:blue','xkcd:aqua','xkcd:lime green','xkcd:neon green','xkcd:green','xkcd:gold','xkcd:golden rod','xkcd:light orange','xkcd:orange','xkcd:red orange','xkcd:red','xkcd:dark red',]

        # leg="vector #"+str(v)
        #plt.rc('xtick', labelsize=8)
        #plt.xticks(rotation=90, fontsize='8')

        matplotlib.rc('xtick', labelsize=8) 
        matplotlib.rc('ytick', labelsize=8) 
            
        print(np.shape(self.occupancies))
        print(self.occupancies[:])
    
        # initiate plot
        n_cols = 2
        n_rows = int(math.ceil(self.numvec/n_cols))
        fig, axs = plt.subplots(n_rows, n_cols, figsize=(5 * n_rows, 5 * n_cols),
                                squeeze=False)  # , constrained_layout=True)
        col = -1
        row = -1
        for vec in range(self.numvec):
            # go the next position in the plot
            if col == 0:
                col = 1
            else:
                col = 0
                row += 1
            # fig.subplots_adjust(left=0.1, bottom=0.25, right=0.9, top=0.95,wspace=0.6, hspace=0.5)
            axs[(row,col)].bar((self.occupancies[:]),np.abs(vh[vec,:]), color=colorlist[vec]) #,label=leg,marker="o",linewidth=1, markersize=6)
            axs[(row,col)].set(xlabel='dataset', ylabel='Amplitude')
            #ax.set_xlim([0.001,1000000])
            print((self.occupancies[:]))
            print(np.abs(vh[vec,:]))
            print((vh[vec,:]))
            print("")
            axs[(row,col)].set_ylim([0,1])
            axs[(row,col)].axhline(y=0.4, color='gray' , linestyle='--')
            #ax.axhline(y=0.3, color='gray' , linestyle='--')
            for tick in axs[(row,col)].xaxis.get_ticklabels():
                tick.set_fontsize('medium')
                tick.set_rotation(90)

            for tick in axs[((row,col))].yaxis.get_ticklabels():
                tick.set_fontsize('medium')
            
        fig.tight_layout()
        outname = "{s}_right_singular_values.png".format(self.prefix)
        plt.savefig(outname, dpi=300)
    
    def read_ccp4_map(self, mapinname):
        map_object = CCP4_Maps(mapinname)
        map_object.open_map()
        print('Read map file {} with size {}'.format(mapinname, map_object.data.shape))
        return map_object
    
    def write_ccp4_map(self, maparrray, uc, grid):
        
        if self.xray_structure != None:
            sites_cart = self.xray_structure.sites_cart()
            
        outname = "{:s}".format(self.prefix)
            
        iotbx.map_tools.write_ccp4_map(
        sites_cart=sites_cart,                     # Cartesian coordinates of atoms (usually from your model)
        unit_cell=uc,                              # Unit cell parameters (a, b, c, alpha, beta, gamma)
        map_data=maparrray,                        # 3D map as a flex array or numpy array
        n_real=grid,                               # Grid dimensions (tuple of 3 ints)
        buffer=5.0,                                # Buffer size (default 5.0)
        file_name=outname)                         # Output file name
        
        return outname
        
    def run_svd_analysis(self):
    
        for n in range(len(self.map_2mFextr_DFc_list)):
            # size, start, intervals, uc, order, skew, skew_trn, symops, _, totalmap = readccp4map(self.map_2mFextr_DFc_list[n])
            maparray = self.read_ccp4_map(self.map_2mFextr_DFc_list[n])
            totalmap=maparray.data/np.std(maparray.data) #Scale the maps to have the same standard deviation #Is this necessary ?
            try:
                dataset[:,n] = np.reshape(totalmap, (length))
            except NameError:
                length = totalmap.size
                nmaps=len(self.map_2mFextr_DFc_list)
                dataset = np.zeros((length, nmaps), dtype='float32')
                dataset[:,n] = np.reshape(totalmap, (length))
                
        print('Done. Read', np.shape(dataset),'voxels into memory.')
        
        print('Now performing Singular Value Decomposition...')
        u, s, vh, success = self.run_svd(dataset)
        
        #reflect = float(str(float(scipy.linalg.det(u) * scipy.linalg.det(v))))
        #if reflect == -1.0:
        #    s[-1] = -s[-1]
        #    u[:,-1] = -u[:,-1]

        # implement ? Eo should be the inital residuals
        #RMSD = E0 - (2.0 * sum(s))
        #RMSD = numpy.sqrt(abs(RMSD / L))

        #print('checking for reflections')
        #print(u)
        #print(s)
        #print(v)

        self.plot_singular_values(s)

        print('')
        print('Now cleaning up the maps...')
        print('There are',len(self.map_2mFextr_DFc_list),'vectors of which we are using the first',self.numvec,'.')
        s_prime=np.zeros((len(self.map_2mFextr_DFc_list),len(self.map_2mFextr_DFc_list)))
        for vec in range(self.numvec):
            s_prime[vec,vec]=s[vec]
            
        #####s_prime[1,1]=0.0 ##################REMOVE WHEN DONE#########################
        dataset_cleaned=np.dot(u,np.dot(s_prime,vh))
        print('')
        print('Writing out the cleaned maps...')
        for d in range(len(self.map_2mFextr_DFc_list)):
            cleanedmap=dataset_cleaned[:,d]
            originalmap=dataset[:,d]
            C=np.corrcoef(originalmap,cleanedmap)
            print('Correlations between original and cleaned for map',d,'at occupancy',self.occupancies[d],':')
            print(C[0,1])
            cleanedmap=np.reshape(cleanedmap,totalmap.shape)
            std=np.std(cleanedmap)
            cleanedmap=cleanedmap/std
            mapoutname = self.write_ccp4_map(cleanedmap, maparray.unit_cell, maparray.grid)

        self.plot_right_singular_values()

        print("Code still in development")
        sys.exit()

        pymolloadername=os.path.join(outdir,'vectormaps_w1.00.pml')
        pymolloader=open(pymolloadername,"a")
        print('Writing out the left singular vectors... ')
        print('pymolcommands: ')
        print('################')
        pymolloader.write('reinitialize \n')
        pymolloader.write('load %s/darkmodel.pdb \n' %outdir)
        pymolloader.write('select PIA, resn PIA \n')

        #vectorzero= -1 * np.reshape(u[:,0],(x,y,z))
        #vectorzero=vectorzero/np.std(vectorzero)


        #if (np.mean(vectorzero)) < 0 :
        #    invertsvd0 = False
        #else:
        #    invertsvd0 = True


        #if (np.max(vectorzero)) - (np.min(vectorzero)) < 0 :
        #    invertsvd = True
        #else:
        #    invertsvd = False

        for d in range(numvec):
            vector=u[:,d]
            m, bins, patches = plt.hist(vector, 100, density=1, facecolor='green', alpha=0.5)
            print("%s ---- %.5f  %.5f  %.5f ----" %(d, np.max(vector), np.median(vector), np.min(vector)))
            
            vector=np.reshape(vector,totalmap.shape) # * -1
            print("%s ---- %.5f  %.5f  %.5f ----" %(d, np.max(vector), np.median(vector), np.min(vector)))
            #vector=vector/np.std(vector)

            #what actually work for negative values
            #testit= -1 * vector  + np.reshape(u[:,0],(x,y,z)) 

            
            #if invertsvd0: 
                #testit = -1 * testit 
                #vectorzero = -1 * vectorzero
            #if (vectorzero+vector) > (vectorzero-vector) 
                
            #for positive values the subtraction addition should work
            
                
            
            outname=os.path.join(outdir,'SVD_w1.00_original'+str(d)+'.ccp4')
            #outtestit=outdir+'testit_'+str(d)+'.ccp4'
            vector=vector/np.std(vector)
            #testit=testit/np.std(testit)

            #print("%s ---- %.5f  %.5f  %.5f ---- %.5f  %.5f  %.5f" %(d, np.max(vector), np.mean(vector), np.min(vector), np.max(testit), np.mean(testit), np.min(testit)))

            writeccp4map_full(outname,size,start,intervals,uc,order,skew,skew_trn,symops,vector)
            #writeccp4map_full(outtestit,size,start,intervals,uc,order,skew,skew_trn,symops,testit)

            pymolloader.write('load %s, vec_%s \n'%(outname,d))
            pymolloader.write('isomesh pos_%s, vec_%s ,3.0 \n' %(d,d))
            pymolloader.write('isomesh neg_%s, vec_%s ,3.0 \n' %(d,d))
        pymolloader.write('color green, pos_* \n')
        pymolloader.write('color red, pos_* \n')

        pymolloader.close()
    


class Filefinder(object):
    def __init__(self,
                 X8_outdir ="Xtrapol8",
                 X8_outname = "Xtrapol8",
                 X8_list_occ = [0.1],
                 X8_f_extrapolated_and_maps = "qFextr",
                 f_extrapolated_and_maps = "qFextr"):
        self.X8_outdir                  = X8_outdir
        self.X8_outname                 = X8_outname
        self.X8_list_occ                = X8_list_occ
        self.X8_f_extrapolated_and_maps = X8_f_extrapolated_and_maps
        self.f_extrapolated_and_maps    = f_extrapolated_and_maps
    
    def find_2fextfc(self):
        """
        Find the 2Fextr-Fc type of maps given the X8_outdir, X8_outname, the X8_list_occ list and X8_f_extrapolated_and_maps 
        """
        if self.f_extrapolated_and_maps not in self.X8_f_extrapolated_and_maps:
            print("ESFA file type not found in Xtrapol8 output (this might be a bug)")
        
        if self.f_extrapolated_and_maps.startswith("q"):
            first_part = "qweight_"
        elif self.f_extrapolated_and_maps.startswith("k"):
            first_part = "kweight_"
        else:
            first_part = ''
            
        maptype = re.sub("f","F", self.f_extrapolated_and_maps)
        last_part = "2m{:s}-DFc.ccp4".format(maptype)
                
        map_2fextrfcalc_list = []
        for occ in self.X8_list_occ:
            f = "{:s}/{:s}occupancy_{:.3f}/{:s}_occ{:.3f}_{:s}".format(self.X8_outdir, first_part, occ, self.X8_outname, occ, last_part)
            map_2fextrfcalc_list.append(os.path.abspath(check_file_existance(f)))
            
        return map_2fextrfcalc_list

if __name__ == "__main__":
    # run(sys.argv[1:])
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
            .help = The type of ESFAs for which the SVD map analysis will be carried out. The Xtrapol8 run prior to these analysis should include the ESFA type of choice. You can only specify one, launch mutliple runs if you want to repeat on with different ESFA types.
            .expert_level = 0
        }
    output{
        outdir = SVD_analysis
            .type = str
            .help = Output directory. 'SVD_analysis' be used if not specified.
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
           raise Usage("phenix.python svd_analysis.py + [.phil] + [arguments]\n arguments only overwrite .phil if provided last")
           sys.exit(1)
    if "--help" in sys.argv or "-h" in sys.argv:
           master_phil.show(attributes_level=1)
           raise Usage("phenix.python svd_analysis.py + [.phil] + [arguments]\n arguments only overwrite .phil if provided last")
           sys.exit(1)

    #Extract input from inputfile and command line
    input_objects = iotbx.phil.process_command_line_with_files(
        args=sys.argv[1:],
        master_phil=master_phil
        )
    params = input_objects.work.extract()
    
    #Extract info from Xtrapol8 run
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
    
    model_pdb = Xtrapol8_params.input.reference_pdb
    pdb_ini = iotbx.pdb.input(model_pdb)
    xray_structure = pdb_ini.xray_structure_simple()
    
    additional_files = Xtrapol8_params.input.additional_files

    map_2fextrfcalc_list = Filefinder(X8_outdir = Xtrapol8_params.output.outdir,
                             X8_outname = Xtrapol8_params.output.outname,
                             X8_f_extrapolated_and_maps = Xtrapol8_params.f_and_maps.f_extrapolated_and_maps,
                             X8_list_occ = Xtrapol8_params.occupancies.list_occ,
                             f_extrapolated_and_maps = params.input.f_extrapolated_and_maps).find_2fextfc()
    
    occupancies = Xtrapol8_params.occupancies.list_occ
    
    if len(map_2fextrfcalc_list) != len(occupancies):
        print("Number of occupancies and mFextr-DFcalc maps is not equal. Please provide a single occupancy for each map.")
        sys.exit(1)
        
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
        log = open('svd_analysis_%s.log' %(suffix), 'w')
    else:
        log = open(params.output.log_file, 'w')
        
    SVD_analysis(map_2mFextr_DFc_list = map_2fextrfcalc_list,
                 occupancies = occupancies,
                 xray_structure=xray_structure,
                 prefix = suffix,
                 log = log).run_svd_analysis()
    
