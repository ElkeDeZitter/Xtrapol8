import JK_utils
import multiprocessing
import JK_image_merging_and_create_mtz_file
import os
from functools import partial

import Fextr_utils
#QUESTIONS???
#Xtrapol8 with which dark mtz and light mtz?

# TODO:
### verif si argument d'entree pas overwritten  si ajoute dans other pour fonction crystfel
### dans JK_utils run_in_terminal ameliorer subprocess.Popen (voir si fichiers crees 2 fois)
### check if programms executable OK
### check symmetry unit cell (+gauss) OK
#get unit cell from pdb if also X8 and choose what to get in case there are different
### get pointgroup, system... from spacegroup
### Stream module to clean
### muliple process: get right order in log file
### adapt for light and dark state
### if spacegroup written with spaces, get rid of them
#ajout processor
#ajout high and low resolution

def image_merging_and_create_mtz_file(dir_streamfile, pointgroup, other_process_hkl, other_partialator, cell,
             other_stats_compare_hkl, a, b, c, alpha, beta, gamma, dir_cryst_prog, spacegroup, method_process_hkl, method_partialator):
    '''
    For each stream file, merging intensities, getting figure of merit, creating mtz file in the corresponding directory
    Args:
        dir_streamfile:
            directory of the stream file
        pointgroup:
        other_process_hkl:
        other_partialator:
        cell:
        other_stats_compare_hkl:
        a:
        b:
        c:
        alpha:
        beta:
        gamma:
        dir_cryst_prog:
        spacegroup:
        logdir:
        method_process_hkl:
        method_partialator:

    Returns:
        outdirfiles:
            directory where to fid the mtz file
        mtzoutdir:
            mtz file full directory
    '''
    outdirfiles = dir_streamfile[0] #getting directory of the output files and the input stream file
    stream_fle = dir_streamfile[1] #getting the stream file

    os.chdir(outdirfiles) #changing directory to folder with stream file and where to put output files
    output_file_name = Fextr_utils.get_name(stream_fle) #get name of output files = stream file name

    #1. Merging intensities (creating hkl files)
    if method_process_hkl: #if the method to use is montecarlo
        hkl_file, do_hkl12_statistics = JK_image_merging_and_create_mtz_file.Image_merging_and_create_mtz_file(stream_fle, output_file_name, pointgroup, dir_cryst_prog).merge_I_montecarlo(other_process_hkl)
        #merging intensities with montecarlo

    if method_partialator: #if the method to use is partialator
        hkl_file, do_hkl12_statistics = JK_image_merging_and_create_mtz_file.Image_merging_and_create_mtz_file(stream_fle, output_file_name, pointgroup, dir_cryst_prog).merge_I_partialator(other_partialator)
        #merging intensities with partialator

    #2. Get figures of merit
    statistics_file_name = JK_image_merging_and_create_mtz_file.Image_merging_and_create_mtz_file(stream_fle, output_file_name, pointgroup, dir_cryst_prog).statistics( cell, other_stats_compare_hkl, outdirfiles, do_hkl12_statistics)
    JK_utils.print_terminal_and_log('The figures of merit are regrouped in the file : %s' % (statistics_file_name))

    #3. Create mtz files
    mtzoutfile = JK_image_merging_and_create_mtz_file.Image_merging_and_create_mtz_file(stream_fle, output_file_name, pointgroup, dir_cryst_prog).create_mtz(spacegroup, a, b, c, alpha, beta, gamma, hkl_file)
    JK_utils.print_terminal_and_log('mtz file %s created in %s' % (mtzoutfile, outdirfiles))
    mtzoutdir = outdirfiles + '/' + mtzoutfile #get directory of mtz file

    return [outdirfiles, mtzoutdir] #returning directory of mtz file

def run_JK(P, outdir,  stream_file, stream_file_name, n_frames_to_keep, system, pointgroup, unique_axis, a, b, c, alpha, beta, gamma,
                                             spacegroup, state='off', total=False):
    '''

    Args:
        P:
            class of parameters
        outdir:
            output directory for JK
        stream_file:
        stream_file_name:
        n_frames_to_keep:
        system:
        pointgroup:
        unique_axis:
        a:
        b:
        c:
        alpha:
        beta:
        gamma:
        spacegroup:
        state:
            'off' or 'on'
            state of the crystal (dark or triggered)
        total:
            bool
            the function is run with a complete stream file

    Returns:
        mtzoutdirs:
            list of [directory where to find mtz file, mtz file complete directory]
    '''

    JK_utils.print_terminal_and_log('CREATING NEW DIRECTORIES AND MOVING LOG FILE TO NEW DIRECTORY\n==============================================================')
    #Create new directories
    newoutdir, table_dir_streamfile = JK_utils.create_files(stream_file_name, P.percentage, n_frames_to_keep,
                                                                P.repeats, stream_file, outdir, state=state, total=total)
    JK_utils.print_terminal_and_log('The output directory of the files for Jack Knife is : %s' %(newoutdir))

    JK_utils.print_terminal_and_log('CREATING CELL FILE\n==============================================================')
    #creating cell file containing the symmetry of the crystal
    cell= JK_utils.create_cell_file(system, unique_axis, a, b, c, alpha, beta, gamma, newoutdir,
                                         spacegroup)

    JK_utils.print_terminal_and_log('MERGING IMAGES AND CREATING MTZ FILES\n==============================================================')
    #for each section of images (new stream file): merging images, getting figure of merit and creating mtz file - simultaniously
    pool_process = multiprocessing.Pool(P.processors)  # Create a multiprocessing Pool with the number of processors defined
    mtzoutdirs=pool_process.map(partial(image_merging_and_create_mtz_file, pointgroup=pointgroup, other_process_hkl=P.other_process_hkl, other_partialator=P.other_partialator, cell=cell, other_stats_compare_hkl=P.other_stats_compare_hkl, a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma, dir_cryst_prog=P.dir_cryst_prog,  spacegroup=spacegroup, method_process_hkl=P.method_process_hkl, method_partialator=P.method_partialator), table_dir_streamfile) # process the function image_merging_and_create_mtz_file with dir_streamfile iterable with pool in table_dir_streamfile

    return (mtzoutdirs) #list of [directory where to find mtz file, mtz file complete directory]
