from __future__ import division, print_function
import os
from datetime import datetime

#for test of part of code
#log_main =open('/mnt/ibs-equipe-weik.04/edezitter/PaulaOeser/2021-08-03_08h29_Xtrapol8.log', "w")

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

def print_terminal_and_log(x, log=None): #needs to be done that way instead of log=log_main because the program will not run it, giving an error on the assignement of log_main
    '''
    Print x in the log file and the terminal
    Args:
        x:
            str, float; int, bool ...
            element to print
    Returns:
    '''
    if log!=None:
        print(x, file = log)
        print(x)
    else:
        print(x, file = log_main)
        print(x)

def update_log_main(new_log_main):
    '''
    Update the global log_main into the new_log_main
    Args:
        new_log_main: directory of the log file to update as log_main

    Returns:

    '''
    log_main = open(new_log_main, "w")
    global log_main

