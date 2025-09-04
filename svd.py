from __future__ import print_function
import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.ndimage
import scipy
import scipy.stats
import scipy.signal

    import glob
    import sys, os
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

#####################################################################################################

def run(args):
    """
    Perform Singular Value Decomposition on a series of maps
    Input:  args[0] - path to the maps
            args[1] - output directory
    Output: cleaned maps, singular values, singular vectors as png files, pymol loader script
        """
    
    #best SVD ALL 1.8 A: 8,15,8,10,8,8,15
    #best SVD ALL 1.6 A: 10, 15, 10, 15, 8, 8, 10  
    #best SVD TUNNEL 1.6 A: 10, 15, 8, 8, 8, 8, 10   or 8 10 8 8 8 8 10 or 8 10 8 10 8 8 10
    #best SVD MINIMAL 1.6 A: 10, 15, 8, 8, 8, 8, 10   or 8 10 8 8 8 8 10 or 8 10 8 10 8 8 10   - nothing really good

    #bestbest SVD ALL 1.6 A BSHARP: 8 15 8 10 8 8 15

    path=args[0]
    output_dir = args[1]

    mapinnames = sorted(glob.glob(os.path.join(path,f'kweight_occupancy_0.???/*_occ0.???_2mkFextr-DFc.ccp4')))
    #mapinnames = sorted(glob.glob(os.path.join(path,f'occupancy_0.???/*_occ0.???_2mFextr-DFc.ccp4')))
    times = [float(f.split('occupancy_')[1][:5]) for f in mapinnames]

    labels = [str(x) for x in times] #['0.0', '0.05', '0.075', '0.09', '0.10', '0.11', '0.12', '0.13', '0.14', '0.15', '0.16', '0.17', '0.20', '0.30', '0.40', '0.50', '0.60']

    print(labels)
    print(times)

    import os
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    size, start, interval, cell, order, skew, skewtrn, SYMOPS, mrc_start, totalmap = readccp4map(mapinnames[0])

    s = np.shape(totalmap)
    x = s[0]
    y = s[1]
    z = s[2]
    length = x*y*z
    nmaps=len(mapinnames)
    assert length == totalmap.size
    dataset = np.zeros((length, nmaps), dtype='float32')


    for n in range(len(mapinnames)):
        size,start,interval,cell,order,skew,skewtrn,SYMOPS,mrc_start, totalmap = readccp4map(mapinnames[n])

        #If this scale necessary ?
        totalmap=totalmap/np.std(totalmap)

        print(f'Read map file {mapinnames[n]} with size {np.shape(totalmap)}')
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




    plt.bar(np.arange(len(mapinnames)),S)
    plt.title('Singular values for w1.00')
    axes = plt.gca()
    axes.set_xlim([-2,12])
    axes.set_ylim([0,5000])
    #plt.show()
    #outfig1=output_dir+'Singular_values.png'
    plt.savefig(output_dir+'Singular_values_w1.00.png', dpi=300)
    plt.close()

    #plt.bar(np.arange(len(mapinnames)),S)
    #plt.title('Singular values')
    #plt.show()

    #numvec=input('How many vectors do you want to use ? ')
    #numvec=int(numvec)

    numvec=5
    #numvec=sys.argv[1]


    print('')
    print('Now cleaning up the maps...')
    print('There are',len(mapinnames),'vectors of which we are using the first',numvec,'.')
    Sprime=np.zeros((len(mapinnames),len(mapinnames)))
    for v in range(numvec):
        Sprime[v,v]=S[v]
        
        
    #####Sprime[1,1]=0.0 ##################REMOVE WHEN DONE#########################
    dataset_cleaned=np.dot(U,np.dot(Sprime,V))
    print('')
    print('Writing out the cleaned maps...')
    print('pymolcommands:')
    print('################')
    print('load dark.pdb, dark')
    for d in range(len(mapinnames)):
        cleanedmap=dataset_cleaned[:,d]
        #m, bins, patches = plt.hist(vector, 100, normed=1, facecolor='green', alpha=0.5)
        #plt.show()

        cleanedmap=np.reshape(cleanedmap,(x,y,z))
        std=np.std(cleanedmap)
        cleanedmap=cleanedmap/std
        mapoutname=os.path.join(output_dir,'cleaned_time'+str(d)+'.ccp4')
        writeccp4map_full(mapoutname,size,start,interval,cell,order,skew,skewtrn,SYMOPS,cleanedmap)

        print('load',mapoutname,', cleanedmap,'+str(d+1))
    print('isomesh cleaned'+str(d)+', cleanedmap ,1.5')
        #print('isomesh neg'+str(d)+', cleanedmap'+str(d)+' ,-3.0')
    print('color slate, cleaned')
        #print('color red, neg'+str(d))


    print(len(mapinnames))
    print(dataset.shape)
    print(len(times))
    for d in range(len(mapinnames)):
        cleanedmap=dataset_cleaned[:,d]
        originalmap=dataset[:,d]
        C=np.corrcoef(originalmap,cleanedmap)
        print('Correlations between original and cleaned for map',d,'at time',times[d],':')
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
    """
    fig, axs = plt.subplots(3, 5, figsize=(20,12))
    #axs[0, 0].semilogx((times[:]),np.abs(V[0,:]), c='xkcd:purple',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[0, 0].set_title('Vector 0')
    #axs[0, 1].semilogx((times[:]),np.abs(V[1,:]), c='xkcd:blue',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[0, 1].set_title('Vector 1')
    #axs[0, 2].semilogx((times[:]),np.abs(V[2,:]), c='xkcd:aqua',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[0, 2].set_title('Vector 2')
    #axs[0, 3].semilogx((times[:]),np.abs(V[3,:]), c='xkcd:neon green',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[0, 3].set_title('Vector 3')
    #axs[0, 4].semilogx((times[:]),np.abs(V[4,:]), c='xkcd:lime green',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[0, 4].set_title('Vector 4')



    #axs[1, 0].semilogx((times[:]),np.abs(V[5,:]), c='xkcd:golder rod green',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[1, 0].set_title('Vector 5')
    #axs[1, 1].semilogx((times[:]),np.abs(V[6,:]), c='xkcd:light orange',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[1, 1].set_title('Vector 6')
    #axs[1, 2].semilogx((times[:]),np.abs(V[7,:]), c='xkcd:green',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[1, 2].set_title('Vector 7')
    #axs[1, 3].semilogx((times[:]),np.abs(V[8,:]), c='xkcd:gold',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[1, 3].set_title('Vector 8')
    #axs[1, 4].semilogx((times[:]),np.abs(V[9,:]), c='xkcd:golden rod',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[1, 4].set_title('Vector 10')


    #axs[2, 0].semilogx((times[:]),np.abs(V[10,:]), c='xkcd:light orange',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[2, 0].set_title('Vector 11')
    #axs[2, 1].semilogx((times[:]),np.abs(V[11,:]), c='xkcd:orange',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[2, 1].set_title('Vector 12')
    #axs[2, 2].semilogx((times[:]),np.abs(V[12,:]), c='xkcd:red orange',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[2, 2].set_title('Vector 13')
    #axs[2, 3].semilogx((times[:]),np.abs(V[13,:]), c='xkcd:red',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[2, 3].set_title('Vector 14')
    #axs[2, 4].semilogx((times[:]),np.abs(V[14,:]), c='xkcd:dark red',label=leg,marker="o",linewidth=1, markersize=6)
    #axs[2, 4].set_title('Vector 15')



    axs[0, 0].bar((labels[:]),np.abs(V[0,:]), color=colorlist)#, c='xkcd:purple',label=leg,marker="o",linewidth=1, markersize=6)
    axs[0, 0].set_title('Vector 0')
    axs[0, 1].bar((labels[:]),np.abs(V[1,:]), color=colorlist)#, c='xkcd:royal blue',label=leg,marker="o",linewidth=1, markersize=6)
    axs[0, 1].set_title('Vector 1')
    axs[0, 2].bar((labels[:]),np.abs(V[2,:]), color=colorlist)#, c='xkcd:blue',label=leg,marker="o",linewidth=1, markersize=6)
    axs[0, 2].set_title('Vector 2')
    axs[0, 3].bar((labels[:]),np.abs(V[3,:]), color=colorlist)#, c='xkcd:baby blue',label=leg,marker="o",linewidth=1, markersize=6)
    axs[0, 3].set_title('Vector 3')
    axs[0, 4].bar((labels[:]),np.abs(V[4,:]), color=colorlist)#, c='xkcd:aqua',label=leg,marker="o",linewidth=1, markersize=6)
    axs[0, 4].set_title('Vector 4')



    axs[1, 0].bar((labels[:]),np.abs(V[5,:]), color=colorlist)#, c='xkcd:lime green',label=leg,marker="o",linewidth=1, markersize=6)
    axs[1, 0].set_title('Vector 5')
    axs[1, 1].bar((labels[:]),np.abs(V[6,:]), color=colorlist)#, c='xkcd:neon green',label=leg,marker="o",linewidth=1, markersize=6)
    axs[1, 1].set_title('Vector 6')
    axs[1, 2].bar((labels[:]),np.abs(V[7,:]), color=colorlist)#, c='xkcd:green',label=leg,marker="o",linewidth=1, markersize=6)
    axs[1, 2].set_title('Vector 7')
    axs[1, 3].bar((labels[:]),np.abs(V[8,:]), color=colorlist)#, c='xkcd:gold',label=leg,marker="o",linewidth=1, markersize=6)
    axs[1, 3].set_title('Vector 8')
    axs[1, 4].bar((labels[:]),np.abs(V[9,:]), color=colorlist)#, c='xkcd:golden rod',label=leg,marker="o",linewidth=1, markersize=6)
    axs[1, 4].set_title('Vector 10')


    axs[2, 0].bar((labels[:]),np.abs(V[10,:]), color=colorlist)#, c='xkcd:light orange',label=leg,marker="o",linewidth=1, markersize=6)
    axs[2, 0].set_title('Vector 11')
    axs[2, 1].bar((labels[:]),np.abs(V[11,:]), color=colorlist)#, c='xkcd:orange',label=leg,marker="o",linewidth=1, markersize=6)
    axs[2, 1].set_title('Vector 12')
    axs[2, 2].bar((labels[:]),np.abs(V[12,:]), color=colorlist)#, c='xkcd:red orange',label=leg,marker="o",linewidth=1, markersize=6)
    axs[2, 2].set_title('Vector 13')
    axs[2, 3].bar((labels[:]),np.abs(V[13,:]), color=colorlist)#, c='xkcd:red',label=leg,marker="o",linewidth=1, markersize=6)
    axs[2, 3].set_title('Vector 14')
    axs[2, 4].bar((labels[:]),np.abs(V[14,:]), color=colorlist)#, c='xkcd:dark red',label=leg,marker="o",linewidth=1, markersize=6)
    axs[2, 4].set_title('Vector 15')

    fig.subplots_adjust(left=0.1, bottom=0.15, right=0.9, top=0.95,wspace=0.6, hspace=0.5)

    for h in range(0:14):
        print(labels[:]),np.abs(V[h,:])
        print("ok")
        
    #fig.canvas.draw()

    #axs.subplots_adjust(left=0.125, bottom=0.15, right=0.9, top=0.9,wspace=0.2, hspace=0.2)


    for ax in axs.flat:
        #ax.set(xlabel='Time (ps)', ylabel='Amplitude')
        #ax.set_xlim([0.001,1000000])
        #ax.set_ylim([0,1])
        #ax.axhline(y=0.5, color='gray' , linestyle='--')
        #ax.axhline(y=0.3, color='gray' , linestyle='--')
        
        ax.set(xlabel='dataset', ylabel='Amplitude')
        
        #ax.set_xlim([0.001,1000000])
        ax.set_ylim([0,1])
        ax.axhline(y=0.45, color='gray' , linestyle='--')
        #ax.axhline(y=0.3, color='gray' , linestyle='--')
        for tick in ax.xaxis.get_ticklabels():
            tick.set_fontsize(8)
        tick.set_rotation(90)

        #plt.show()
        #outfig1=output_dir+'Singular_values.png'
    plt.savefig(output_dir+'Waves_w1.00.png', dpi=300)
    """

    print(V[v,:])
    for d in range(len(mapinnames)):
        print(V[v,:d])
        
    print(np.shape(labels))
    print (labels[:])
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
        ax.bar((labels[:]),np.abs(V[v,:]), color=colorlist[v]) #,label=leg,marker="o",linewidth=1, markersize=6)
        ax.set(xlabel='dataset', ylabel='Amplitude')
        #ax.set_xlim([0.001,1000000])
        print((labels[:]))
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
        
    
        plt.savefig(os.path.join(output_dir,'Waves_vector%i.png' %v), dpi=300)




    pymolloadername=os.path.join(output_dir,'vectormaps_w1.00.pml')
    pymolloader=open(pymolloadername,"a")
    print('Writing out the left singular vectors... ')
    print('pymolcommands: ')
    print('################')
    pymolloader.write('reinitialize \n')
    pymolloader.write('load %s/darkmodel.pdb \n' %output_dir)
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
        
            
        
        outname=os.path.join(output_dir,'SVD_w1.00_original'+str(d)+'.ccp4')
        #outtestit=output_dir+'testit_'+str(d)+'.ccp4'
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

if __name__ == "__main__":
    run(sys.argv[1:])
