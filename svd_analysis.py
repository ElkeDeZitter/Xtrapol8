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

import glob
import sys
import os
#######################READCCP4MAP################
#
def readccp4map(mapinname):
    """
    Read a ccp4 map file
    Input:  mapinname - name of the ccp4 map file
    Output: SIZE - number of columns, rows, slices
            START - starting point in columns, rows, slices
            INTERVALS - number of intervals in columns, rows, slices
            CELL - unit cell parameters a,b,c,alpha,beta,gamma
            ORDER - order of columns, rows, slices
            SKEW - the skew matrix
            SKEWTRN - the skew translation vector
            SYMOPS - the symmetry operators
            MRC_START - start of the map in the original mrc file
            maparray - the 3D array of the map values
    """
    mapinfile=open(mapinname, "rb")
    datain=mapinfile.read(208)

    mapinfo=struct.unpack("10i 6f 3i 3f 3i 12f 15f", datain)
    NC=mapinfo[0]
    NR=mapinfo[1]
    NS=mapinfo[2]
    SIZE=[NC,NR,NS]
    MODE=mapinfo[3]
    NCSTART=mapinfo[4]
    NRSTART=mapinfo[5]
    NSSTART=mapinfo[6]
    START=[NCSTART,NRSTART,NSSTART]
    NX=mapinfo[7]
    NY=mapinfo[8]
    NZ=mapinfo[9]
    INTERVALS=[NX,NY,NZ]
    A=mapinfo[10]
    B=mapinfo[11]
    C=mapinfo[12]
    ALPHA=mapinfo[13]
    BETA=mapinfo[14]
    GAMMA=mapinfo[15]
    CELL=np.asarray([A,B,C,ALPHA,BETA,GAMMA])
    MAPC=mapinfo[16]
    MAPR=mapinfo[17]
    MAPS=mapinfo[18]
    ORDER=[MAPC,MAPR,MAPS]
    #print 'Order of columns, rows, slices:',ORDER
    AMIN=mapinfo[19]
    AMAX=mapinfo[20]
    AMEAN=mapinfo[21]
    ISPG=mapinfo[22]
    NSYMBT=mapinfo[23]
    LSKFLG=mapinfo[24]
    SKEW=[[mapinfo[25],mapinfo[26],mapinfo[27]],
          [mapinfo[28],mapinfo[29],mapinfo[30]],
          [mapinfo[31],mapinfo[32],mapinfo[33]]]
    SKEW=np.asarray(SKEW)
    SKEWTRN=mapinfo[34:37]
    SKEWTRN=np.asarray(SKEWTRN)
    FUTURE=mapinfo[37:52]

    MRC_START=[mapinfo[49],mapinfo[50],mapinfo[51]]

    datain=mapinfile.read(8) # read the word 'MAP ' and the 4-byte machine stamp

    datain=mapinfile.read(8)

    mapinfo=struct.unpack("fi",datain)
    ARMS=mapinfo[0]
    NLABL=mapinfo[1]


    LABEL=mapinfile.read(800)
    SYMOPS=mapinfile.read(NSYMBT)

    NVOXELS=NC*NR*NS
    formatstring=str(NVOXELS)+'f'
    datain=mapinfile.read(NVOXELS*4)
    maparray=struct.unpack(formatstring,datain)
    mapinfile.close()
    maparray=np.asarray(maparray)
    maparray=np.reshape(maparray, (NS,NR,NC))
    return SIZE,START,INTERVALS,CELL,ORDER,SKEW,SKEWTRN,SYMOPS,MRC_START,maparray




#######################skew########################
# Calculate the skew or scale matrix from unit cell
# parameters
#
def skew(a,b,c,alpha,beta,gamma):
    """
    Calculate the skew or scale matrix from unit cell parameters
    Input:  a,b,c - unit cell lengths
            alpha,beta,gamma - unit cell angles in degrees
    Output: s - the 3x3 skew matrix
    """
    pi=3.14159265359
    alpha=pi*alpha/180.0
    beta=pi*beta/180.0
    gamma=pi*gamma/180.0

    vsq=1.0 - (np.cos(alpha)*np.cos(alpha)) - (np.cos(beta)*np.cos(beta)) - (np.cos(gamma)*np.cos(gamma)) - (2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
    v=np.sqrt(vsq) # volume of unit parallelepiped

    
    s=np.zeros((3,3))
    s[0,0]=1.0/a
    s[0,1]=-1.0*np.cos(gamma)/(a*np.sin(gamma))
    s[0,2]=(np.cos(alpha)*np.cos(gamma)-np.cos(beta)) / (a*v*np.sin(gamma))
    s[1,1]=1.0/(b*np.sin(gamma))
    s[1,2]=(np.cos(beta)*np.cos(gamma)-np.cos(alpha)) / (b*v*np.sin(gamma))
    s[2,2]=np.sin(gamma)/(c*v)
    return s

###################################################
#######################WRITECCP4MAP_FULL####################
#
def writeccp4map_full(mapoutname,SIZE,START,INTERVALS,CELL,ORDER,SKEW,SKEWTRN,SYMOPS,maparray):
    """
    Write a ccp4 map file
    Input:  mapoutname - name of the ccp4 map file to be written
            SIZE - number of columns, rows, slices
            START - starting point in columns, rows, slices
            INTERVALS - number of intervals in columns, rows, slices
            CELL - unit cell parameters a,b,c,alpha,beta,gamma
            ORDER - order of columns, rows, slices
            SKEW - the skew matrix
            SKEWTRN - the skew translation vector
            SYMOPS - the symmetry operators
            maparray - the 3D array of the map values
    Output: a ccp4 map file
    """
    HEADER=bytes()

    #print 'Now in writeccp4map function'
    #print 'SIZE looks like:',SIZE
    for value in SIZE:
        HEADER += struct.pack('i',value)
    #print 'Setting the mode to 2...'    
    HEADER += struct.pack('i',int(2)) #SET THE MODE TO 2
    #print 'I think START looks like:',START
    for value in START:
        HEADER += struct.pack('i',value)
    #print 'I think INTERVALS looks like:',INTERVALS
    for value in INTERVALS:
        HEADER += struct.pack('i',value)
    #print 'I think CELL looks like:',CELL

    for value in CELL:
        HEADER += struct.pack('f',value)
    #print 'I think ORDER looks like:',ORDER
    for value in ORDER:
        HEADER += struct.pack('i',value)
        
    AMIN=np.min(maparray)
    AMAX=np.max(maparray)
    AMEAN=np.average(maparray)
    #print 'Min,max,av:',AMIN,AMAX,AMEAN
    HEADER += struct.pack('3f',AMIN,AMAX,AMEAN)
    #print 'I think INTERVALS looks like:',INTERVALS
    #for value in INTERVALS:
    #    HEADER += struct.pack('i',value)
    
    ISPG=1
    NSYMBT=0
    #LSKFLG=1 # IF SKEW IS TO BE WRITTEN OUT
    LSKFLG=0 # IF NO SKEW E.G. FOR PYMOL
    
    HEADER += struct.pack('3i',ISPG,NSYMBT,LSKFLG)

    #SKEW=np.reshape(SKEW,(9)) # IF SKEW IS TO BE WRITTEN OUT
    SKEW=np.zeros((9)) # IF NO SKEW E.G. FOR PYMOL
    SKEWTRN=np.zeros((3))
    #print 'I think the SKEW looks like:',SKEW
    for value in SKEW:
        #print value
        HEADER += struct.pack('f',value)

    #print 'I think the SKEWTRN looks like:',SKEWTRN
    for value in SKEWTRN:
        HEADER += struct.pack('f',value)
    #print 'adding 15 zeros...'
    for n in range(15):
        #print n,
        HEADER += struct.pack('i',int(0))
    #print ''
    #print 'Adding the word MAP_'
    HEADER += struct.pack('4c',b'M',b'A',b'P',b' ')
    # little-endian machine stamp 'D','A',0,0
    #HEADER += struct.pack('4c',chr(0x44), chr(0x41), chr(0x00), chr(0x00))
    HEADER += struct.pack('4c',b'D',b'A',b'0',b'0')

    ARMS=np.std(maparray)
    HEADER += struct.pack('f',ARMS)
    HEADER += struct.pack('i',3)
    for n in range(800):
        HEADER += struct.pack('c',b' ')

    mapoutfile=open(mapoutname, "wb")
    mapoutfile.write(HEADER)
    mapsize=np.shape(maparray)
    #print 'The map size is',mapsize
    maplength=(mapsize[0]*mapsize[1]*mapsize[2])
    #print 'The map length is',maplength
    maparray=np.reshape(maparray,(maplength))
    for n in range(maplength):
        voxelout=struct.pack('f',maparray[n])
        mapoutfile.write(voxelout)
    mapoutfile.close()
############################################################################## 


class SVD_analysis(object):
    def __init__(self,
                    map_2mFextr_DFc_list = [],
                    occupancies=[],
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
        
    def run_svd_analysis(self):
    
        size, start, interval, cell, order, skew, skewtrn, SYMOPS, mrc_start, totalmap = readccp4map(self.map_2mFextr_DFc_list[0])

        s = np.shape(totalmap)
        x = s[0]
        y = s[1]
        z = s[2]
        length = x*y*z
        nmaps=len(self.map_2mFextr_DFc_list)
        assert length == totalmap.size
        dataset = np.zeros((length, nmaps), dtype='float32')


        for n in range(len(self.map_2mFextr_DFc_list)):
            size,start,interval,cell,order,skew,skewtrn,SYMOPS,mrc_start, totalmap = readccp4map(self.map_2mFextr_DFc_list[n])

            #If this scale necessary ?
            totalmap=totalmap/np.std(totalmap)

            print('Read map file {} with size {}'.format(self.map_2mFextr_DFc_list[n], np.shape(totalmap)))
            dataset[:,n] = np.reshape(totalmap, (length))

        print('Done. Read', np.shape(dataset),'voxels into memory.')
        print('Now performing Singular Value Decomposition...')

        #U,S,V=np.linalg.svd(dataset,full_matrices=False)
        U, S, V=scipy.linalg.svd(dataset, full_matrices=False)
        print('Done. The prefactors for the singular values are:')
        print(U.shape)
        print(S.shape)
        print(V.shape)


        # we already have our solution, in the results from SVD.
        # we just need to check for reflections. U and V are orthonormal,
        # so their det's are +/-1.
        #print(U.shape)
        #print(V.shape)



        print(scipy.linalg.det(V))

        #reflect = float(str(float(scipy.linalg.det(U) * scipy.linalg.det(V))))
        #if reflect == -1.0:
        #    S[-1] = -S[-1]
        #    U[:,-1] = -U[:,-1]

        # implement ? Eo should be the inital residuals
        #RMSD = E0 - (2.0 * sum(S))
        #RMSD = numpy.sqrt(abs(RMSD / L))

        #print('checking for reflections')
        #print(U)
        #print(S)
        #print(V)




        plt.bar(np.arange(len(self.map_2mFextr_DFc_list)),S)
        plt.title('Singular values for w1.00')
        axes = plt.gca()
        axes.set_xlim([-2,12])
        axes.set_ylim([0,5000])
        #plt.show()
        #outfig1=outdir+'Singular_values.png'
        plt.savefig(outdir+'Singular_values_w1.00.png', dpi=300)
        plt.close()

        #plt.bar(np.arange(len(map_2mFextr_DFc_list)),S)
        #plt.title('Singular values')
        #plt.show()

        #numvec=input('How many vectors do you want to use ? ')
        #numvec=int(numvec)

        numvec=5
        #numvec=sys.argv[1]


        print('')
        print('Now cleaning up the maps...')
        print('There are',len(self.map_2mFextr_DFc_list),'vectors of which we are using the first',numvec,'.')
        Sprime=np.zeros((len(self.map_2mFextr_DFc_list),len(self.map_2mFextr_DFc_list)))
        for v in range(numvec):
            Sprime[v,v]=S[v]
            
            
        #####Sprime[1,1]=0.0 ##################REMOVE WHEN DONE#########################
        dataset_cleaned=np.dot(U,np.dot(Sprime,V))
        print('')
        print('Writing out the cleaned maps...')
        print('pymolcommands:')
        print('################')
        print('load dark.pdb, dark')
        for d in range(len(self.map_2mFextr_DFc_list)):
            cleanedmap=dataset_cleaned[:,d]
            #m, bins, patches = plt.hist(vector, 100, normed=1, facecolor='green', alpha=0.5)
            #plt.show()

            cleanedmap=np.reshape(cleanedmap,(x,y,z))
            std=np.std(cleanedmap)
            cleanedmap=cleanedmap/std
            mapoutname=os.path.join(outdir,'cleaned_time'+str(d)+'.ccp4')
            writeccp4map_full(mapoutname,size,start,interval,cell,order,skew,skewtrn,SYMOPS,cleanedmap)

            print('load',mapoutname,', cleanedmap,'+str(d+1))
        print('isomesh cleaned'+str(d)+', cleanedmap ,1.5')
            #print('isomesh neg'+str(d)+', cleanedmap'+str(d)+' ,-3.0')
        print('color slate, cleaned')
            #print('color red, neg'+str(d))


        print(len(self.map_2mFextr_DFc_list))
        print(dataset.shape)
        print(len(self.occupancies))
        for d in range(len(self.map_2mFextr_DFc_list)):
            cleanedmap=dataset_cleaned[:,d]
            originalmap=dataset[:,d]
            C=np.corrcoef(originalmap,cleanedmap)
            print('Correlations between original and cleaned for map',d,'at time',self.occupancies[d],':')
            print(C[0,1])

        #for v in range(numvec):
        #    leg="vector #"+str(v)
        #    plt.semilogx((times[:]),V[v,:],label=leg)
        #    print('****** Right Singular Vector #'+str(v)+' ******')
        #    print(V[v,:])

        #plt.title('Time dependence (right singular vector) #'+str(v))
        #plt.xlabel('Time')
        #plt.ylabel('Amplitude')
        #plt.legend()
        #splt.show()

        colorlist=['xkcd:purple','xkcd:royal blue','xkcd:blue','xkcd:blue','xkcd:aqua','xkcd:lime green','xkcd:neon green','xkcd:green','xkcd:gold','xkcd:golden rod','xkcd:light orange','xkcd:orange','xkcd:red orange','xkcd:red','xkcd:dark red',]
        #colorlist=['xkcd:purple','xkcd:blue','xkcd:aqua','xkcd:lime green','xkcd:wheat','xkcd:tangerine','xkcd:red','xkcd:dark red',]

        leg="vector #"+str(v)
        #plt.rc('xtick', labelsize=8)
        #plt.xticks(rotation=90, fontsize='8')

        matplotlib.rc('xtick', labelsize=8) 
        matplotlib.rc('ytick', labelsize=8) 

        print(V[v,:])
        for d in range(len(self.map_2mFextr_DFc_list)):
            print(V[v,:d])
            
        print(np.shape(occupancies))
        print (occupancies[:])
        for v in range(numvec):
            #fig, ax = plt.subplots(1,1)
            #ax.semilogx((times[:]),np.abs(V[v,:]), c=colorlist[v],label=leg,marker="o",linewidth=1, markersize=6)
            #ax.set(xlabel='Time (ps)', ylabel='Amplitude')
            #ax.set_xlim([0.001,1000000])
            #ax.set_ylim([0,1])
            #ax.axhline(y=0.5, color='gray' , linestyle='--')
            #ax.axhline(y=0.3, color='gray' , linestyle='--')
        
            fig, ax = plt.subplots(1,1)
            fig.subplots_adjust(left=0.1, bottom=0.25, right=0.9, top=0.95,wspace=0.6, hspace=0.5)
            ax.bar((occupancies[:]),np.abs(V[v,:]), color=colorlist[v]) #,label=leg,marker="o",linewidth=1, markersize=6)
            ax.set(xlabel='dataset', ylabel='Amplitude')
            #ax.set_xlim([0.001,1000000])
            print((occupancies[:]))
            print(np.abs(V[v,:]))
            print((V[v,:]))
            print("")
            ax.set_ylim([0,1])
            ax.axhline(y=0.4, color='gray' , linestyle='--')
            #ax.axhline(y=0.3, color='gray' , linestyle='--')
            for tick in ax.xaxis.get_ticklabels():
                tick.set_fontsize('medium')
                tick.set_rotation(90)

            for tick in ax.yaxis.get_ticklabels():
                tick.set_fontsize('medium')
            #tick.set_rotation(90)
            
        
            plt.savefig(os.path.join(outdir,'Waves_vector%i.png' %v), dpi=300)


        pymolloadername=os.path.join(outdir,'vectormaps_w1.00.pml')
        pymolloader=open(pymolloadername,"a")
        print('Writing out the left singular vectors... ')
        print('pymolcommands: ')
        print('################')
        pymolloader.write('reinitialize \n')
        pymolloader.write('load %s/darkmodel.pdb \n' %outdir)
        pymolloader.write('select PIA, resn PIA \n')

        #vectorzero= -1 * np.reshape(U[:,0],(x,y,z))
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
            vector=U[:,d]
            m, bins, patches = plt.hist(vector, 100, density=1, facecolor='green', alpha=0.5)
            print("%s ---- %.5f  %.5f  %.5f ----" %(d, np.max(vector), np.median(vector), np.min(vector)))
            
            vector=np.reshape(vector,(x,y,z)) # * -1
            print("%s ---- %.5f  %.5f  %.5f ----" %(d, np.max(vector), np.median(vector), np.min(vector)))
            #vector=vector/np.std(vector)

            #what actually work for negative values
            #testit= -1 * vector  + np.reshape(U[:,0],(x,y,z)) 

            
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

            writeccp4map_full(outname,size,start,interval,cell,order,skew,skewtrn,SYMOPS,vector)
            #writeccp4map_full(outtestit,size,start,interval,cell,order,skew,skewtrn,SYMOPS,testit)

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
                 prefix = suffix,
                 log = log).run_svd_analysis()
    
