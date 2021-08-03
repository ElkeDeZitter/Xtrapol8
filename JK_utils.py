from __future__ import division, print_function
import re
import os
import sys
import subprocess
import shutil
from datetime import datetime
from iotbx.file_reader import any_file

import time

from Stream_EDZ import Stream

#class Fextr DH
# def from_cif_create_pdb_file(file, outdir):
#     """
#     If the input model file (pdb_in) has a format cif, a new pdb file of the model is created and replaces the cif file (pdb_in)
#
#     Parameters:
#     -----------
#     self.pdb_in (file)
#         model file input by the user
#     self.outdir (string)
#         output directory
#     file_format (format, example: .cif)
#         format of the model file input by the user
#
#     Returns:
#     -----------
#     self.pdb_in (pdb file)
#         Created model pdb file from the model cif input file
#     """
#
#     try:  # iotbx.file_reader.splitext from phenix 1.19
#         _, file_format, _ = iotbx.file_reader.splitext(file)
#     except ValueError:  # iotbx.file_reader.splitext from phenix 1.18
#         _, file_format = os.path.splitext(file)
#
#     # get the file path and type
#     if file_format == '.cif':
#         # if the file type is cif
#         pdb_hier = hierarchy.input(file_name=file)
#         p = open(file + '/' + Fextr_utils.get_name(file) + '.pdb', 'w')
#         # create and open a new pdb file to write in, the file will be in the new outdir and have the same name as the input file with a pdb format
#         p.write(pdb_hier.hierarchy.as_pdb_string(crystal_symmetry=pdb_hier.input.crystal_symmetry()))
#         # write the columns of the cif file into the pdb file
#         p.close()
#         if check_single_file(outdir + '/' + Fextr_utils.get_name(file) + '.pdb'):
#             # check if the created pdb file exists
#             file = os.path.abspath(outdir + '/' + Fextr_utils.get_name(file) + '.pdb')
#             # the input model used is the created pdb file with the info from cif file
#     return file
#
# def check_and_delete_hydrogen(pdb_in, outdir):
#     """
#     If hydrogen exists, make new pdb-file with removed hydrogens.
#     Possibility to achieve extrapolated structure factors with such high precision as to see hydrogens is very low.
#     """
#     pdb_hier = hierarchy.input(file_name=pdb_in)
#     outname = 'model_edit.pdb'
#     if pdb_hier.hierarchy.remove_hd() != 0:
#         print("Model contains hydrogen atoms. Create a new model without these atoms in the output directory: %s" %(outname))
#         print("Model contains hydrogen atoms. Create a new model without these atoms in the output directory: %s" % (outname), file=log)
#         pdb_hier.hierarchy.remove_hd()
#         p = open('%s/%s' %(outdir, outname), 'w')
#         p.write(pdb_hier.hierarchy.as_pdb_string(crystal_symmetry=pdb_hier.input.crystal_symmetry()))
#         p.close()
#         if check_single_file('%s/%s' %(outdir, outname)):
#             pdb_in = os.path.abspath('%s/%s' %(outdir, outname))
#         else:
#             pdb_in = os.path.abspath(pdb_in)
#     else:
#         pdb_in = os.path.abspath(pdb_in) #We will need absolute path of pdb file for later refinements in subdirectories
#     return pdb_in
#
# def get_UC_and_SG(model_in):
#     """
#     Extract unit cell and space group from the model.
#     """
#     #self.SG = re.search(r"(.+?)\(No",self.fobs_off.space_group_info().symbol_and_number()).group(1)
#     #self.UC = self.fobs_off.unit_cell()
#     self.SG = str(self.model_in.crystal_symmetry().space_group_info())
#     self.UC = self.model_in.crystal_symmetry().unit_cell()

#class JK_0
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
            list of lists of 2 argument [[outdirfiles, newstreamfiledir],[outdirfiles, newstreamfiledir]...]
            list of output directory for every portion of images
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

    table_dir_streamfile = [] #initiate list

    #create all folders and stream files
    i = 0 #index
    while i < repeats: #index for the number of divisions wanted

    # 1.create folders
        outdirfiles = "%s/%s_JK_%ipercent_%i" % (newoutdir, stream_file_name, percentage, i)
        if os.path.isdir(outdirfiles) == False:
            os.mkdir(outdirfiles)
        print('%s directory created' % (outdirfiles))

    #2.create stream files
        newstreamfiledir = '%s/%s_%iimages.stream' % (outdirfiles, stream_file_name, n_frames_to_keep) #new stream file directory with new name

        Stream(stream_file).save_random_selection(n_frames_to_keep, newstreamfiledir) #creates a new stream file with only the number of random images to keep
        print_terminal_and_log('%s stream file created in %s' % (newstreamfiledir, outdirfiles))

        table_dir_streamfile.append([outdirfiles, newstreamfiledir]) #add list of output directory and new stream file directory to the list

        # print the summary of the stream file created
        print_terminal_and_log("----> %s : %s_%iimages.stream <---- " % (i, stream_file_name,  n_frames_to_keep))
        Stream(newstreamfiledir).get_stream_summary()
        print_terminal_and_log("------------------")

        i += 1 #index

    #add total stream file to the list
    if total == True:
        total_outdirfiles = newoutdir + '/' + stream_file_name + '_total'
        os.mkdir(total_outdirfiles)
        table_dir_streamfile.append([total_outdirfiles, stream_file])  # add list of total output directory and stream file directory to the list

    return (newoutdir, table_dir_streamfile)

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

#class global
def print_file_content(fle_in_dir, fle_out_dir):
    '''
    Copy the content of one file into a new file

    Args:
        fle_in_dir:
            file directory
            input file directory
        fle_out_dir:
            open file with 'w' or file directory
            output file
    Returns:

    '''
    if os.path.isfile(fle_in_dir) == False:
        filelines = None
        print('%s does not exist' %(fle_in_dir))
    else:
        infile = open(fle_in_dir, 'r') # open file
        filelines = infile.read().splitlines() # read lines
        infile.close() # close file

    if filelines != None:

        #if the fle_out_dir is a file directory
        if type(fle_out_dir) == str:

            #create file if doesn't exist
            if os.path.isfile(fle_out_dir) == False:
                os.mknod(fle_out_dir)

            # lines of input file printed in output file
            fle_out_dir=open(fle_out_dir, 'w')
            print('------------------------------------------------------------------------------------------', file=fle_out_dir)
            for line in filelines:
                print(line, file=fle_out_dir)
            print('------------------------------------------------------------------------------------------', file=fle_out_dir)
            fle_in_dir.close()
            fle_out_dir.close()

        # if the fle_out_dir is an open file
        elif type(fle_out_dir == file):
            if fle_out_dir.mode == 'w':
                #lines of input file printed in output file
                print('------------------------------------------------------------------------------------------', file=fle_out_dir)
                for line in filelines:
                    print(line, file=fle_out_dir)
                print('------------------------------------------------------------------------------------------', file=fle_out_dir)
            else:
                # lines of input file printed in output file
                fle_out_dir = open(fle_out_dir, 'w')
                print('------------------------------------------------------------------------------------------', file=fle_out_dir)
                for line in filelines:
                    print(line, file=fle_out_dir)
                print('------------------------------------------------------------------------------------------', file=fle_out_dir)
                fle_out_dir.close()

        print_terminal_and_log('%s content copied into %s' % (fle_in_dir, fle_out_dir))
    return (filelines)

def run_in_terminal (cmd, existing_files=False, wait=True):
    '''
    Run the command in the terminal and wait until it has completed

    Args:
        cmd:
            str
            command to execute in terminal
        existing_files:
            list of file(s) that has to exist before continuing the code, if False no file has to exist
        wait:
            bool
            the code will not be continued until the command is finished if True
    Returns:
        cmd_done:
            bool
            the command was executed in terminal
        cmd_output:
            str
            str printed in the terminal after execution of the command line
    '''

    cmd_done = False #the command hasn't been executed
    print_terminal_and_log(cmd) #print the command in the terminal and log file

    if wait == True:
        #??? How to improve code and not do command twice?
        P = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait() #launch the command and wait until finished
        p = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE)  #launch the command
        _, cmd_output = p.communicate() #get the output of the terminal
        print_terminal_and_log(cmd_output) #print the output of the terminal in terminal and log file
#        p.wait() #wait until the command finished ???

        i=0
        while P != 0: #the code will not continue while the command isn't finished
            time.sleep(10)  # time delay of 10s
            #to avoid unended run
            i+=1
            if i==100:
                print_terminal_and_log("this is not working or you are running heavy stuff on a slow machine")
                sys.exit()
                break

        if existing_files != False:
            i=0
            for file in existing_files:
                 while not os.path.isfile(file): #the code will not continue while the file doesn't exist
                    time.sleep(10) #time delay of 10s
                    #to avoid unended run
                    i += 1
                    if i == 100:
                        print_terminal_and_log("this is not working or you are running heavy stuff on a slow machine")
                        sys.exit()
                        break

                 print_terminal_and_log('%s created' % (file)) #print the file was created

        cmd_done = True #the command has been executed


    else:
        P = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE) # launch the command
        cmd_done = True  # the command has been executed
        _, cmd_output = P.communicate() #get the output of the terminal
        print_terminal_and_log(cmd_output) #print the output of the terminal in terminal and log file

    return (cmd_done, cmd_output)


#class log():

def create_log(outdir):
    '''
    Generate log-file in the output directory

    Args:
        outdir:
            dierctory
            directory of the created log file
    Returns:
    '''

    now = datetime.now().strftime('%Y-%m-%d_%Hh%M') #get date
    logname = "%s/%s_Xtrapol8.log" % (outdir, now) #directory and name of the log file
    #if a log file with that name already exists, change the log name to log name_i (with i a number)
    i = 1
    while os.path.isfile(logname):
        logname = "%s/%s_Xtrapol8_%d.log" % (outdir, now, i)
        i += 1
        if i == 50:
            break

    #write the header for log file
    log = open(logname, "w")
    global log
    print("Xtrapol8 -- version 0.9.5 -- run date: %s" % (now), file=log)
    print('-----------------------------------------')
    print("Xtrapol8 -- version 0.9.5 -- run date: %s" % (now))
    return (logname, log)

def move_log (logname, outdir, newoutdir):
    '''
    Move log file to new output directory
    Args:
        logname:
            str
            log file name
        outdir:
            path
            current directory of the log file
        newoutdir:
            path
            new directory where to move the log file
    Returns:

    '''
    if os.path.isfile(logname): #if the log file exists
        shutil.move(logname, logname.replace(outdir, newoutdir))

def print_terminal_and_log(x, log1=None):
    '''
    Print x in the log file and the terminal
    Args:
        x:
            str, float; int, bool ...
            element to print
    Returns:
    '''
    if log1!=None and os.path.isfile(log1):
        print(x, file=log1)
        print(x)
    else:
        print(x, file=log)
        print(x)



#From Fextr DH
def check_single_file(fle):
    if fle == None:
        return False
    else:
        return os.path.isfile(fle)



#Class get params
#From Fextr DH
def get_UC_and_SG(self):
    """
    Extract unit cell and space group from the reference data set
    """
    SG = re.search(r"(.+?)\(No",self.fobs_off.space_group_info().symbol_and_number()).group(1)
    #What is fobs_off???
    UC = self.fobs_off.unit_cell()
    return UC, SG

    model_in = any_file("DATA/darkmodel.pdb", force_type="pdb", raise_sorry_if_errors=True)
    test = model_in.crystal_symmetry().unit_cell()
    a, b, c, alf, bet, gam = test.parameters()

# class Stream():
# #From Stream_EDZ
#
#
#     def save_random_selection(self, n, outstreamfilename):
#         #??? indexed_images, frames and header? What are those?
#         '''
#         Get a random selection of images of a stream file
#
#         Args:
#             outstreamfilename: stream_file_directory_and_name
#             n: number of samples wanted
#         indexed_images
#         header
#
#         Returns:
#             creates new stream file with number of samples wanted named 'root_nimages.stream'
#         '''
#         total = self.indexed_images
#         sele = random.sample(range(total),n)
#
#         out = open(outstreamfilename, 'w')
#         print(self.header, file=out)
#         print ('Saving %d frames to %s' %(n, outstreamfilename))
#         for i in sele:
#             print >> out, ''.join(self.frames[i].all), #print(''.join(self.frames[i].all), file=out)
#         out.close()