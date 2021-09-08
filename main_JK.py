from __future__ import division, print_function
import multiprocessing
import os
from functools import partial

from Stream_EDZ import Stream
import JK_image_merging_and_create_mtz_file
import Log_file
from Log_file import print_terminal_and_log
import Fextr_utils

def create_files (stream_file_name, percentage, n_frames_to_keep, repeats, stream_file, outdir, state=None, total=False):
    '''
    Create the folder for the state dark/light and in it for the number of iterations of the percentage wanted and in them create the stream files with random selection
    Args:
        stream_file_name:
        percentage:
        n_frames_to_keep:
        repeats:
        stream_file:
        outdir:
        state:
            'off' or 'on'
            state of the crystal (dark or triggered)
        total:
            bool
            the function is run with a complete stream file

    Returns:
        newoutdir:
            str
            new output path
        table_dir_streamfile:
            list of lists of 4 arguments [[outdirfiles, newstreamfiledir, log, total],[outdirfiles, newstreamfiledir, log, i]...]
            list of output directory for every portion of images
        table_dir_streamfile_total:
            list of 4 arguments [outdirfiles, newstreamfiledir, log, 'total'] or empty list
            list of output directory for all images
    '''

    print_terminal_and_log('>>> Creating directories <<<')

    #create new folder for state dark/light or nothing
    if state == 'off':
        newoutdir = outdir + '/JK_' + stream_file_name +'_reference'
    elif state == 'on':
        newoutdir = outdir + '/JK_' + stream_file_name +'_triggered'
    else:
        newoutdir = outdir + '/JK_' + stream_file_name

    if os.path.isdir(newoutdir) == False:
        os.mkdir(newoutdir)
    else:
        i=1
        while os.path.isdir(newoutdir):
            if state == 'off':
                newoutdir = outdir + '/JK_' + stream_file_name + '_reference' + str(i)
            elif state == 'on':
                newoutdir = outdir + '/JK_' + stream_file_name + '_triggered' + str(i)
            else:
                newoutdir = outdir + '/JK_' + stream_file_name + str(i)

            i += 1
        os.mkdir(newoutdir)
    print_terminal_and_log('new directory created: %s' % (newoutdir))

    table_dir_streamfile_total = []  # initiate list for total
    table_dir_streamfile = [] #initiate list for other JK

    #create all folders and stream files
    i = 0 #index
    while i < repeats: #index for the number of divisions wanted

    # 1.create folders
        outdirfiles = "%s/%s_JK_%ipercent_%i" % (newoutdir, stream_file_name, percentage, i)
        if os.path.isdir(outdirfiles) == False:
            os.mkdir(outdirfiles)
        print('%s directory created' % (outdirfiles))

    # 1,5.create log file
        logname, log = Log_file.create_log(outdirfiles)

    #2.create stream files
        newstreamfiledir = '%s/%s_%iimages.stream' % (outdirfiles, stream_file_name, n_frames_to_keep) #new stream file directory with new name

        Stream(stream_file).save_random_selection(n_frames_to_keep, newstreamfiledir) #creates a new stream file with only the number of random images to keep
        print_terminal_and_log('%s stream file created in %s' % (newstreamfiledir, outdirfiles))

        table_dir_streamfile.append([outdirfiles, newstreamfiledir, logname, i]) #add list of output directory and new stream file directory to the list

        # print the summary of the stream file created
        print_terminal_and_log("----> %s : %s_%iimages.stream <---- " % (i, stream_file_name,  n_frames_to_keep))
        Stream(newstreamfiledir).get_stream_summary()
        print_terminal_and_log("------------------")

        i += 1 #index

    #IF TOTAL : Create folder and list with total stream file (all images)
    if total == True:
        total_outdirfiles = newoutdir + '/' + stream_file_name + '_total'
        os.mkdir(total_outdirfiles)
        logname, log = Log_file.create_log(total_outdirfiles)
        table_dir_streamfile_total = [total_outdirfiles, stream_file, logname, 'total']  # add list of total output directory and stream file directory to the total list

    return (newoutdir, table_dir_streamfile, table_dir_streamfile_total)

def create_cell_file (system, unique_axis, a, b, c, alpha, beta, gamma, newoutdir, spacegroup):
    '''
    Adapted from: Thomas Barends
    Create the cell file containing the crystal symmetry

    Args:
        system:
        unique_axis:
        a:
        b:
        c:
        alpha:
        beta:
        gamma:
        newoutdir:

    Returns:
        cell
            file .cell
            cell file directory containing the symmetry of the crystal

    '''
    cell_name = newoutdir + '/' + 'crystal_symmetry.cell'
    cell = open(cell_name, 'w')
    cell.write('CrystFEL unit cell file version 1.0\n')
    cell.write('\n')
    cell.write('lattice_type = ' + str(system) + '\n')
    spg = str(spacegroup)
    centering = spg[0]
    if unique_axis!=None:
        cell.write('unique_axis = ' + str(unique_axis) + '\n')
    cell.write('centering = ' + centering + '\n')
    cell.write('a = ' + str(a) + ' A\n')
    cell.write('b = ' + str(b) + ' A\n')
    cell.write('c = ' + str(c) + ' A\n')
    cell.write('al = ' + str(alpha) + ' deg\n')
    cell.write('be = ' + str(beta) + ' deg\n')
    cell.write('ga = ' + str(gamma) + ' deg\n')
    cell.close()
    print_terminal_and_log('Cell file created: spacegroup=%s, UnitCell=%s %s %s %s %s %s' %(spacegroup, a, b, c, alpha, beta, gamma))

    return (cell_name)

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
            directory where to find the mtz file
        mtzoutdir:
            mtz file full directory
        total:
            Bool, file is with total data
    '''
    outdirfiles = dir_streamfile[0] #getting directory of the output files and the input stream file
    stream_fle = dir_streamfile[1] #getting the stream file
    logname = dir_streamfile[2] #getting the log file
    total = dir_streamfile[3] #'total'if the files are for the total data and number of JK if not
    log = open(logname, "w")
    os.chdir(outdirfiles) #changing directory to folder with stream file and where to put output files
    output_file_name = Fextr_utils.get_name(stream_fle) #get name of output files = stream file name

    #1. Merging intensities (creating hkl files)
    if method_process_hkl: #if the method to use is montecarlo
        hkl_file, do_hkl12_statistics = JK_image_merging_and_create_mtz_file.Image_merging_and_create_mtz_file(stream_fle, output_file_name, pointgroup, dir_cryst_prog, log).merge_I_montecarlo(other_process_hkl)
        #merging intensities with montecarlo

    if method_partialator: #if the method to use is partialator
        hkl_file, do_hkl12_statistics = JK_image_merging_and_create_mtz_file.Image_merging_and_create_mtz_file(stream_fle, output_file_name, pointgroup, dir_cryst_prog, log).merge_I_partialator(other_partialator)
        #merging intensities with partialator

    #2. Get figures of merit
    statistics_file_name = JK_image_merging_and_create_mtz_file.Image_merging_and_create_mtz_file(stream_fle, output_file_name, pointgroup, dir_cryst_prog, log).statistics( cell, other_stats_compare_hkl, outdirfiles, do_hkl12_statistics)
    print_terminal_and_log('The figures of merit are regrouped in the file : %s' % (statistics_file_name), log=log)

    #3. Create mtz files
    mtzoutdir = JK_image_merging_and_create_mtz_file.Image_merging_and_create_mtz_file(stream_fle, output_file_name, pointgroup, dir_cryst_prog, log).create_mtz(spacegroup, a, b, c, alpha, beta, gamma, hkl_file)#create mtz and get full directory of mtz file
    print_terminal_and_log('mtz file %s created in %s' % (mtzoutdir, outdirfiles), log=log)

    return [outdirfiles, mtzoutdir, total] #returning directory of mtz file and total or number of JK

def run_JK(P, outdir,  stream_file, stream_file_name, n_frames_to_keep, system, pointgroup, unique_axis, a, b, c, alpha, beta, gamma,
                                             spacegroup, log, state='off', total=False):
    '''
    Creating directories and cell file for JK.
    Then, merging images getting figure of merit of the reflections and creating mtz files.

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
            list of [directory where to find mtz file, mtz file complete directory, True/False:the file is for total data]
        mtzoutdirs_total:
            [directory where to find total mtz file, total mtz file complete directory, 'total']
    '''
#1. Create new directories
    print_terminal_and_log('CREATING NEW DIRECTORIES\n==============================================================')
    #Create new directories
    newoutdir, table_dir_streamfile, table_dir_streamfile_total = create_files(stream_file_name, P.percentage, n_frames_to_keep,
                                                                P.repeats, stream_file, outdir, state=state, total=total)
    print_terminal_and_log('The output directory of the files for Jack Knife is : %s' %(newoutdir), log=log)

#2. Create CELL file
    print_terminal_and_log('CREATING CELL FILE\n==============================================================')
    #creating cell file containing the symmetry of the crystal
    cell= create_cell_file(system, unique_axis, a, b, c, alpha, beta, gamma, newoutdir,
                                         spacegroup)

#3. Merging images and creating mtz files
    print_terminal_and_log('MERGING IMAGES AND CREATING MTZ FILES\n==============================================================')
    #for each section of images (new stream file): merging images, getting figure of merit and creating mtz file - simultaniously
    pool_process = multiprocessing.Pool(P.processors)  # Create a multiprocessing Pool with the number of processors defined
    mtzoutdirs = pool_process.map(partial(image_merging_and_create_mtz_file, pointgroup=pointgroup, other_process_hkl=P.other_process_hkl, other_partialator=P.other_partialator, cell=cell, other_stats_compare_hkl=P.other_stats_compare_hkl, a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma, dir_cryst_prog=P.dir_cryst_prog,  spacegroup=spacegroup, method_process_hkl=P.method_process_hkl, method_partialator=P.method_partialator), table_dir_streamfile) # process the function image_merging_and_create_mtz_file with dir_streamfile iterable with pool in table_dir_streamfile

    #Loop to test without multiprocessing
    # mtzoutdirs=[]
    # for lst in table_dir_streamfile:
    #     mtzoutdir=image_merging_and_create_mtz_file(lst, pointgroup = pointgroup, other_process_hkl = P.other_process_hkl, other_partialator = P.other_partialator, cell = cell, other_stats_compare_hkl = P.other_stats_compare_hkl, a = a, b = b, c = c, alpha = alpha, beta = beta, gamma = gamma, dir_cryst_prog = P.dir_cryst_prog, spacegroup = spacegroup, method_process_hkl = P.method_process_hkl, method_partialator = P.method_partialator)
    #     mtzoutdirs.append(mtzoutdir)

    if total:
        mtzoutdirs_total = image_merging_and_create_mtz_file(table_dir_streamfile_total, pointgroup, P.other_process_hkl, P.other_partialator, cell,
                                          P.other_stats_compare_hkl, a, b, c, alpha, beta, gamma, P.dir_cryst_prog, spacegroup,
                                          P.method_process_hkl, P.method_partialator) #run function image_merging_and_create_mtz_file for the total stream file (all images)
    else: mtzoutdirs_total = None

    return (mtzoutdirs, mtzoutdirs_total) #list of [directory where to find mtz file, mtz file complete directory, nb_JK], [directory where to find total mtz file, total mtz file complete directory, 'total']...

