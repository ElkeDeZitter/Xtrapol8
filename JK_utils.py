from __future__ import division, print_function
import re
import os
import sys
import subprocess
import shutil
from datetime import datetime
from iotbx.file_reader import any_file
import time

#log_main =open('/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/2021-08-03_08h29_Xtrapol8.log', "w")

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

def run_in_terminal (cmd, existing_files=False, wait=True, log=None, Crystfel_program=False):
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
        p = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE)  # launch the command
        P = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait() #launch the command and wait until finished

        if Crystfel_program==True: #special case of the crystfel programs to get the output in the terminal
            _, cmd_output = p.communicate() #get the output of the terminal

    #        p.wait() #wait until the command finished ??? TODO: Improve the code for running a command in the terminal with subprocess.Popen
             #???How to improve code and not do command twice?
        else:
            cmd_output, _=p.communicate() #get the output of the terminal

        i=0
        while P != 0: #the code will not continue while the command isn't finished
            time.sleep(10)  # time delay of 10s TODO: Choose a correct time delay or make it a parameter in the function, or change 'if i==100'
            #to avoid unended run
            i+=1
            if i==100:
                print_terminal_and_log("this is not working or you are running heavy stuff on a slow machine. Did you check your inputs? Unit Cell? Space Group? Unique_Axis? Other_Statistics?")
                sys.exit()
                break

        if existing_files == False:
            cmd_done=True # the command has been executed

    else:
        P = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE) # launch the command
        if existing_files == False:
            cmd_done=True # the command has been executed
        if Crystfel_program==True:
            _, cmd_output = P.communicate() #get the output of the terminal
        else:
            cmd_output, _ = P.communicate() #get the output of the terminal

    if existing_files != False:
        i=0
        for file in existing_files:
             while not os.path.isfile(file) == True: #the code will not continue while the file doesn't exist
                time.sleep(10) #time delay of 10s TODO: Choose a correct time delay or make it a parameter in the function
                #to avoid unended run
                i += 1
                if i == 100:
                    print_terminal_and_log("this is not working or you are running heavy stuff on a slow machine. Did you check your inputs? Unit Cell? Space Group? Unique_Axis? Other_Statistics?")
                    sys.exit()
                    break
             output_outdir=os.getcwd()
             print_terminal_and_log('%s created in %s' % (file, output_outdir), log=log) #print the file was created

             cmd_done = True  # the command has been executed

#    print_terminal_and_log(cmd_output, log=log)  # print the output of the terminal in terminal and log file

    return (cmd_done, cmd_output)

#class log():
def create_log(outdir, global_log=False):
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
    if global_log == True:
        log_main=log
        global log_main
    print("Xtrapol8 -- version 0.9.5 -- run date: %s" % (now), file=log)
    print('-----------------------------------------')
    print("Xtrapol8 -- version 0.9.5 -- run date: %s" % (now))
    return (logname, log)

def print_terminal_and_log(x, log=log_main):
    '''
    Print x in the log file and the terminal
    Args:
        x:
            str, float; int, bool ...
            element to print
    Returns:
    '''
    print(x, file = log)
    print(x)

