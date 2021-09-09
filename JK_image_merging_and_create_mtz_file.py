from __future__ import division, print_function
import os
import sys

from Log_file import print_terminal_and_log
from Fextr_utils import get_name, run_in_terminal, print_file_content

class Image_merging_and_create_mtz_file(object):

    def __init__(self, stream_file, output_file_name, pointgroup, dir_cryst_prog, log):
        self.stream_file = stream_file
        self.output_file_name = output_file_name
        self.pointgroup = pointgroup
        self.dir_cryst_prog = dir_cryst_prog
        self.log = log


    def merge_I_montecarlo(self, other_process_hkl):
        '''
        Merge intensities of the stream file with the montecarlo method

        Args:
            stream_file:
            output_file_name:
            pointgroup:
            other_process_hkl:
                other options to add to the command process_hkl
            dir_cryst_prog:
                directory of the CrystFEL programs

        Returns:
            [hkl_file_name]:
                hkl file
            do_hkl12_statistics
                bool
                the statistics with hkl1 and hkl1 are done

        '''

        print_terminal_and_log('>>>MERGING INTENSITIES WITH MONTE CARLO METHOD<<<\n--------------------%s----------------------' % (self.stream_file))

        # merge intensities with montecarlo:
        process_hkl_done, _ = run_in_terminal("%s/process_hkl -i %s -o %s.hkl -y %s %s" % (
        self.dir_cryst_prog, self.stream_file, self.output_file_name, self.pointgroup, other_process_hkl),
                                                       existing_files=[self.output_file_name + '.hkl'], log=self.log)
        # merging of all the images
        process_hkl1_done, _ = run_in_terminal("%s/process_hkl -i %s -o  %s.hkl1 --even-only -y %s %s" % (
        self.dir_cryst_prog, self.stream_file, self.output_file_name, self.pointgroup, other_process_hkl),
                                                        existing_files=[self.output_file_name + '.hkl1'], log=self.log)
        # merging of the even half of the images
        process_hkl2_done, _ = run_in_terminal("%s/process_hkl -i %s -o %s.hkl2 --odd-only -y %s %s" % (
        self.dir_cryst_prog, self.stream_file, self.output_file_name, self.pointgroup, other_process_hkl),
                                                        existing_files=[self.output_file_name + '.hkl2'], log=self.log)
        # merging of the odd half of the images

        # file names
        hkl_file_name = self.output_file_name + '.hkl'
        even_hkl_file_name = self.output_file_name + '.hkl1'
        odd_hkl_file_name = self.output_file_name + '.hkl2'

        # check if process_hkl done
        if process_hkl1_done == False or process_hkl2_done == False:
            do_hkl12_statistics = False  # parameters for figure of merit false
            print_terminal_and_log('process_hkl did not succeed for the two halfs, the compared figure of merit will not be calculated', log=self.log)
        else:
            do_hkl12_statistics = True  # parameters for figure of merit true

        if process_hkl_done == False:
            print_terminal_and_log('Process_hkl did not succeed', log=self.log)
            sys.exit()  # stop program

        return (hkl_file_name, do_hkl12_statistics)


    def merge_I_partialator(self, other_partialator):
        '''
        Merge intensities of the stream file with the partialator

        Args:
            stream_file:
            output_file_name:
            pointgroup:
            other_partialator:
                other options to add to the command partialator
            dir_cryst_prog:
                directory of the CrystFEL programs

        Returns:
            [hkl_file_name]
        '''

        print_terminal_and_log(
            '>>>MERGING INTENSITIES WITH PARTIALATOR<<<\n--------------------%s----------------------' % (self.stream_file), log=self.log)

        # merging of intensities with partialator
        partialator, _ = run_in_terminal("%s/partialator -i %s -o %s.hkl -y %s %s" % (
        self.dir_cryst_prog, self.stream_file, self.output_file_name, self.pointgroup, other_partialator),
                                                  existing_files=[self.output_file_name + '.hkl', self.output_file_name + '.hkl1',
                                                                  self.output_file_name + '.hkl2'], log=self.log)

        # file names
        hkl_file_name = self.output_file_name + '.hkl'
        hkl1_file_name = self.output_file_name + '.hkl1'
        hkl2_file_name = self.output_file_name + '.hkl2'

        # check if process_hkl done
        if partialator == False:
            print_terminal_and_log('Partialator did not succeed', log=self.log)
            sys.exit()  # stop program
        if partialator == True: do_hkl12_statistics = True  # if the partialator worked, the code can continue and statistics can be calculated

        return (hkl_file_name, do_hkl12_statistics)


    def statistics(self, cell_symmetry, other_stats_compare_hkl, outdir, do_hkl12_statistics):
        '''
        Print figure of merit (check_hkl and compare_hkl Rsplit CC CCstar) in one file

        Args:
            output_file_name:
            pointgroup:
            cell_symmetry:
                cell file with cell symmetry
            other_stats_compare_hkl:
                other options to add to the command compare_hkl
            dir_cryst_prog:
                directory of the CrystFEL programs
            outdir:
                output directory
            do_hkl12_statistics:
                bool
                the statistics with hkl1 and hkl1 are done
        Var:
            statistics_file:
                open file with all the results of the statistics

        Returns:
            statistics_file_name:
                directory of the file containing all figures of merit
        '''

        print_terminal_and_log('>>>CALCULATING STATISTICS<<<\n------------------------------------------')

        # calculate figure of merit check_hkl
        print_terminal_and_log('---check_hkl---', log=self.log)
        _, check_hkl_output = run_in_terminal(
            "%s/check_hkl %s.hkl -y %s -p %s" % (self.dir_cryst_prog, self.output_file_name, self.pointgroup, cell_symmetry),
            wait=True, existing_files=[outdir + '/shells.dat'], log=self.log, Crystfel_program=True)  # run command in terminal

        statistics_file_name = '%s/statistics_%s.hkl.log' % (outdir, self.output_file_name)
        statistics_file = open(statistics_file_name, 'w')  # create file for all figures of merit
        print('---check_hkl---', file=statistics_file)
        print(check_hkl_output, file=statistics_file)  # print global output of check_hkl to figure of merit file
        print_file_content('%s/shells.dat' % (outdir),
                                    statistics_file)  # copy results of figure of merits (shells.dat) into figure of merit file
        os.remove('%s/shells.dat' % (outdir))  # delete shells.dat

        if do_hkl12_statistics == True:
            # calculate figure of merit compare_hkl for 'Rsplit', 'CC','CCstar'
            for stat_type in ['Rsplit', 'CC', 'CCstar']:
                print_terminal_and_log('---compare_hkl--- %s' % (stat_type), log=self.log)
                _, compare_hkl_output = run_in_terminal(
                    "%s/compare_hkl %s.hkl1 %s.hkl2 -y %s -p %s --fom=%s %s" % (
                    self.dir_cryst_prog, self.output_file_name, self.output_file_name, self.pointgroup, cell_symmetry, stat_type,
                    other_stats_compare_hkl), wait=True, existing_files=[outdir + '/shells.dat'], log=self.log, Crystfel_program=True)  # run command in terminal

                print('---compare_hkl--- %s' % (stat_type), file=statistics_file)
                print(compare_hkl_output,
                      file=statistics_file)  # print global output of compare_hkl to figure of merit file
                print_file_content('%s/shells.dat' % (outdir),
                                            statistics_file)  # copy results of figure of merits (shells.dat) into figure of merit file
                os.remove('%s/shells.dat' % (outdir))  # delete shells.dat
        statistics_file.close()

        return (statistics_file_name)

        # def plot_stats():
        #     #stats plot
        #     infile=open('shells.dat','r')
        #     inlines=infile.readlines()
        #     infile.close()
        #     print inlines
        #     D=[]
        #     SNR=[]
        #     for line in inlines[1:]:
        #         elements=line.split()
        #         D.append(float(elements[9]))
        #         SNR.append(float(elements[6]))
        #     D=np.asarray(D)
        #     SNR=np.asarray(SNR)
        #
        #     for n in range(len(SNR)):
        #         if SNR[n]<1.0:
        #             break
        #     resolution=D[n-1]
        #     print ('I think the resolution is',resolution,'Angstrom')
        #
        #     Lest=SNR[0]
        #     kest=0.1
        #     x0est=resolution
        #
        #     Dfit=np.append(D,0)
        #     SNRfit=np.append(SNR,0)
        #
        #
        #     p0=[Lest,kest,x0est]
        #     coeff, var_matrix=curve_fit(logistic,Dfit,SNRfit,p0=p0)
        #
        #     Dplot=np.linspace(0,D[0],100)
        #
        #     linex=[resolution, resolution]
        #     liney=[np.min(SNR), np.max(SNR)]
        #     plt.plot(D,SNR,label='SNR')
        #     plt.plot(Dplot,logistic(Dplot,*coeff),label='fit (not used)')
        #     plt.plot(linex,liney,'k--',label='limit')
        #     plt.title('SNR vs resolution')
        #     plt.xlabel(r'Resolution [$\AA$]')
        #     plt.ylabel('SNR')
        #     plt.legend()
        #     plt.show()


    def create_mtz(self, spacegroup, a, b, c, alpha, beta, gamma, hkl_file):
        '''
        Create the mtz file from the hkl file

        Args:
            pointgroup:
            a:
            b:
            c:
            alpha:
            beta:
            gamma:
            output_file_name:
            hkl_file:
                file .hkl
                hkl file containing the data to transfer (type output_file_name.hkl(1 or 2 or nothing)

        Returns:
            mtzoutfile:
                mtz file created
        '''

        hkl_file_name = get_name(hkl_file)
        print_terminal_and_log('>>>CREATE MTZ for %s<<<\n------------------------------------------' % (hkl_file_name), log=self.log)

        # get file names with format
        dir = os.getcwd()
        mtzoutfile = dir + '/' + self.output_file_name + '.mtz'
        tmphkl = dir + '/' + self.output_file_name + '.temp.hkl'

        # create tmphkl file = hkl file without the lines from 'End of reflections' until end
        hkl_file_open = open(hkl_file, 'r')
        tmphkl_open = open(tmphkl, 'w')
        temp = hkl_file_open.read().splitlines()  # get list of lines of the hkl file

        # print lines of hkl file in tmphkl file until the 'End of reflections' line
        for line in temp:
            if not ('End of reflections' in line):
                print(line, file=tmphkl_open)
            else:
                break
        hkl_file_open.close()
        tmphkl_open.close()

        print_terminal_and_log('Running f2mtz...', log=self.log)
        # run f2mtz in terminal to create mtz file
        run_in_terminal('unset SYMINFO\n'
                                 'f2mtz HKLIN %s HKLOUT %s > out.html << EOF\n'
                                 'TITLE Reflections from CrystFEL\n'
                                 'NAME PROJECT wibble CRYSTAL wibble DATASET wibble\n'
                                 'CELL %s %s %s %s %s %s\n'
                                 'SYMM %s\n'
                                 'SKIP 3\n'
                                 'LABOUT H K L IMEAN SIGIMEAN\n'
                                 'CTYPE  H H H J     Q\n'
                                 'FORMAT "(3(F4.0,1X),F10.2,10X,F10.2)"\n'
                                 'EOF' % (tmphkl, mtzoutfile, a, b, c, alpha, beta, gamma, spacegroup),
                                 existing_files=[mtzoutfile], log=self.log)

        # remove tmphkl file???    os.remove(tmphkl)
        return mtzoutfile
