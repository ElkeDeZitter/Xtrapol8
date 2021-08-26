from __future__ import division, print_function
import re
import os
import sys
import random
import subprocess
import shutil
from select import select
from datetime import datetime
import numpy as np
from iotbx.file_reader import any_file
from iotbx import symmetry
from iotbx.pdb import hierarchy
from cctbx.array_family import flex
from cctbx import maptbx, miller, crystal, xray
import iotbx.phil
import iotbx.map_tools
from libtbx.utils import Usage
import mmtbx.f_model
import mmtbx.map_tools
import mmtbx.maps.utils
from mmtbx import utils
from cctbx import sgtbx
from iotbx import pdb
from mmtbx.scaling.matthews import p_vm_calculator
import scipy.stats
from cctbx import crystal
import multiprocessing
import time

from Stream_EDZ import Stream
import Fextr_utils
import JK_utils
from JK_utils import print_terminal_and_log as print_in_T_and_log

def check_programs_JK():
    '''
    Check if necessary programs are executable

    Args:

    Returns:
    '''

    # check programs
    # check if ccp4 is installed and f2mtz is executable
    if Fextr_utils.check_program_path('f2mtz') == False:
        print_in_T_and_log("ccp4 not found. Make sure ccp4 is correctly instaled and added to PATH.")
        sys.exit()

    # check if CrystFEL is installed and programs are excutable
    if Fextr_utils.check_program_path('process_hkl') == False or Fextr_utils.check_program_path(
            'partialator') == False or Fextr_utils.check_program_path(
        'check_hkl') == False or Fextr_utils.check_program_path('compare_hkl') == False:
        print_in_T_and_log(
            'CrystFEL not found. Make sure CrystFEL is correctly instaled and added to PATH.')
        sys.exit()

def check_programs_X8(params):

    # Check if non-phenix programs can be found:
    if Fextr_utils.check_program_path('coot') == False:
        print_in_T_and_log("COOT not found.")
    if Fextr_utils.check_program_path('scaleit') == False:
        print_in_T_and_log("scaleit not found. Data will not be scaled.")
        params.Xtrapol8.scaling.b_scaling = 'no'
    if params.Xtrapol8.refinement.use_refmac_instead_of_phenix:
        if (Fextr_utils.check_program_path('refmac5') and Fextr_utils.check_program_path('coot')) == False:
            print_in_T_and_log("refmac and/or COOT not found. Phenix will be used for refinement.")
            params.Xtrapol8.refinement.use_refmac_instead_of_phenix = False

class Parameters():
    def __init__(self, params):
        self.params = params

    def get_parameters_JK_X8(self):

        print_in_T_and_log('>>>Getting input parameters from phil file<<<')

        # getting programs to run
        self.run_JackKnife = False
        self.run_Xtrapol8 = False
        if 'run_JackKnife' in self.params.options.programs:
            self.run_JackKnife = True
        if 'run_Xtrapol8' in self.params.options.programs:
            self.run_Xtrapol8 = True
        if self.params.options.programs == []:
            print("%s not defined. Nothing will be run." % (
                self.params.options.programs.__phil_path__()))
            sys.exit()

        self.JK_one_stream_file = False
        if self.params.JackKnife.Off_state.stream_file_off!=None and self.params.JackKnife.On_state.stream_file_on==None and self.run_JackKnife and (not self.run_Xtrapol8):
            self.JK_one_stream_file = True

        #getting number of processors to use
        self.processors = self.params.options.processors
        if self.processors > multiprocessing.cpu_count():
            self.processors = multiprocessing.cpu_count()



    def get_parameters_JK(self):
        '''
        Get the values of all the parameters from the phil file for JackKnife and deduced from them with simple calculs

        Args:
            params:

        Returns:
            repeats, stream_file, stream_file_name, fraction, percentage, n_frames_to_keep, pointgroup, other_process_hkl, other_partialator, other_stats_compare_hkl, a, b ,c, alpha, beta, gamma, system, unique_axis, dir_cryst_prog, method_process_hkl, method_partialator
        '''

    # getting input parameters
        self.repeats = self.params.JackKnife.input.repeats
        if self.repeats == None:
            print_in_T_and_log('If you want to run Jack Knife, please give a value to repeats (Number of times JackKnife will be repeated)')
            sys.exit()

        self.fraction = self.params.JackKnife.input.fraction
        if self.fraction == None:
            print_in_T_and_log('If you want to run Jack Knife, please give a value to fraction (Percentage of images used for JackKnife (fractional))')
            sys.exit()
        # check if the fraction is correct
        if self.fraction == 0:
            print_in_T_and_log("please select fraction between 0.0 and 1.0")
            sys.exit()
        # get percentage
        self.percentage = int(self.fraction * 100)

        self.dir_cryst_prog = self.params.JackKnife.input.directory_crystfel_programs
        if self.dir_cryst_prog == None :
            print_in_T_and_log("%s not defined. No path for crystfel programs will be applied and might cause problems if they are not correctly sourced." % (
                    self.params.JackKnife.input.directory_crystfel_programs.__phil_path__()))
            self.dir_cryst_prog = ''

    # getting Off state parameters
        self.stream_file_off = self.params.JackKnife.Off_state.stream_file_off
        if self.run_JackKnife and self.stream_file_off == None:
            print('JackKnife can not be run with no stream file.')
            sys.exit()
        #check file
        Fextr_utils.check_all_files([self.stream_file_off])

        # read the stream file and get the number of indexed images, and the number of frames to keep
        print_in_T_and_log("---Reading stream file of the reference state---")
        self.S_off = Stream(self.stream_file_off)
        self.stream_file_name_off = Fextr_utils.get_name(self.stream_file_off)
        self.a_off, a_array_off, astdev_off, self.b_off, b_array_off, bstdev_off, self.c_off, c_array_off, cstdev_off, self.alpha_off, alpha_array_off, alphastdev_off, self.beta_off, beta_array_off, betastdev_off, self.gamma_off, gamma_array_off, gammasdtev_off = Stream(self.stream_file_off).get_cell_stats()

        # test for normal distribution
        print_in_T_and_log(' - Checking normal distribution of unit cell parameters - ')

        def normal_distribution_test_unit_cell(array, variable_tested):
            try:
                statistic, p = scipy.stats.normaltest(array)
                print_in_T_and_log('%s : p= %s , statistic= %s' %(variable_tested, p, statistic))
                if p >= 0.05:
                    print_in_T_and_log("Peaks not normal distributed (p-value for normality test = %.3f). The %s value might not be correct. The process will continue but the unit cell might be wrong." % (p, variable_tested))
                    time.sleep(15)  # time delay of 15s
            except ValueError:
                print_in_T_and_log("Test for normality could not be performed. Probably not enough peaks. Jack Knife will not be run because there might not be enough images.")
                sys.exit()

        normal_distribution_test_unit_cell(a_array_off, 'a_off')
        normal_distribution_test_unit_cell(b_array_off, 'b_off')
        normal_distribution_test_unit_cell(c_array_off, 'c_off')
        normal_distribution_test_unit_cell(alpha_array_off, 'alpha_off')
        normal_distribution_test_unit_cell(beta_array_off, 'beta_off')
        normal_distribution_test_unit_cell(gamma_array_off, 'gamma_off')

        # getting unit cell parameters
        print_in_T_and_log('---Getting unit cell parameters---')

        if self.params.JackKnife.Off_state.unit_cell.use_UC_and_SG_from_pdb and self.run_Xtrapol8: #get pdb unit cell parameters ??? if not pdb but cif
            # check file
            Fextr_utils.check_all_files([self.params.Xtrapol8.input.reference_pdb])
            #get unit cell parameters from pdb file
            model_in = any_file(self.params.Xtrapol8.input.reference_pdb, force_type="pdb", raise_sorry_if_errors=True)
            unitcell = model_in.crystal_symmetry().unit_cell()
            self.a_off, self.b_off, self.c_off, self.alpha_off, self.beta_off, self.gamma_off = unitcell.parameters()
            spacegroup_off = str(model_in.crystal_symmetry().space_group_info())
            self.spacegroup_off = ''
            for i in spacegroup_off:
                if i != ' ':
                    self.spacegroup_off += i

        else: #get input unit cell parameters if given

            if self.params.JackKnife.Off_state.unit_cell.a != None:
                self.a_off = self.params.JackKnife.Off_state.unit_cell.a
            else:
                print('Since a was not given, the average a will be taken')
            if self.params.JackKnife.Off_state.unit_cell.b != None:
                self.b_off = self.params.JackKnife.Off_state.unit_cell.b
            else:
                print('Since b was not given, the average b will be taken')
            if self.params.JackKnife.Off_state.unit_cell.c != None:
                self.c_off = self.params.JackKnife.Off_state.unit_cell.c
            else:
                print('Since c was not given, the average c will be taken')
            if self.params.JackKnife.Off_state.unit_cell.alpha != None:
                self.alpha_off = self.params.JackKnife.Off_state.unit_cell.alpha
            else:
                print('Since alpha was not given, the average alpha will be taken')
            if self.params.JackKnife.Off_state.unit_cell.beta != None:
                self.beta_off = self.params.JackKnife.Off_state.unit_cell.beta
            else:
                print('Since beta was not given, the average beta will be taken')
            if self.params.JackKnife.Off_state.unit_cell.gamma != None:
                self.gamma_off = self.params.JackKnife.Off_state.unit_cell.gamma
            else:
                print('Since gamma was not given, the average gamma will be taken')

            # get symmetry parameters
            spacegroup_off = self.params.JackKnife.Off_state.spacegroup
            if spacegroup_off == None:
                print('Please give the space group for JackKnife')
                sys.exit()
            else:
                self.spacegroup_off = ''
                for i in spacegroup_off:
                    if i != ' ':
                        self.spacegroup_off += i

        self.unique_axis_off = self.params.JackKnife.Off_state.unique_axis #??? comment recuperer unique axis du pdb?
        if self.unique_axis_off == None or self.unique_axis_off == 'default':
            self.unique_axis_off = None


        # get pointgroup and system HORRILE BUT SHOULD WORK
        def get_system_and_pointgroup(a, b, c, alpha, beta, gamma, spacegroup, unique_axis):
            '''
            Checking if spacegroup exists
            From spacegroup, getting system and pointgroup
            Checking unit cell correct for spacegroup, if not change according to system
            Args:
                a:
                b:
                c:
                alpha:
                beta:
                gamma:
                spacegroup:
                unique_axis:

            Returns:
                spacegroup, system, pointgroup, a, b, c, alpha, beta, gamma
            '''

            merge_friedel_pairs = True #var to give a spacegroup with which the fliedel pairs will be merged when merging intensities

            # 1. Checking if spacegroup exists
            unit_cell = [a, b, c, alpha, beta, gamma]
            try:
                symmetry = crystal.symmetry(space_group=spacegroup, unit_cell=unit_cell)
            except RuntimeError:
                print('The space group given does not exist, please make sure it is written correctly')
                sys.exit()

            nb_spacegroup = crystal.symmetry.space_group_number(symmetry)
            print('spacegroup number:', nb_spacegroup)

            # 2.  From spacegroup, getting system and pointgroup
            if nb_spacegroup in [1, 2]:
                system = 'triclinic'
                if nb_spacegroup == 1:
                    if merge_friedel_pairs:
                        pointgroup = '-1'
                    else:
                        pointgroup = '1'
                if nb_spacegroup == 2:
                    pointgroup = '-1'

            elif 3 <= nb_spacegroup <= 15:
                system = 'monoclinic'
                if 3 <= nb_spacegroup <= 5:
                    if merge_friedel_pairs:
                        pointgroup = '2/m'
                    else:
                        pointgroup = '2'
                if 6 <= nb_spacegroup <= 9:
                    pointgroup = 'm'
                if 10 <= nb_spacegroup <= 15:
                    pointgroup = '2/m'

            elif 16 <= nb_spacegroup <= 74:
                system = 'orthogonal'
                if 16 <= nb_spacegroup <= 24:
                    if merge_friedel_pairs:
                        pointgroup = 'mmm'
                    else:
                        pointgroup = '222'
                if 25 <= nb_spacegroup <= 46:
                    pointgroup = 'mm2'
                if 47 <= nb_spacegroup <= 74:
                    pointgroup = 'mmm'

            elif 75 <= nb_spacegroup <= 142:
                system = 'tetratgonal'
                if 75 <= nb_spacegroup <= 80:
                    if merge_friedel_pairs:
                        pointgroup = '-4'
                    else:
                        pointgroup = '4'
                if 81 <= nb_spacegroup <= 82:
                    pointgroup = '-4'
                if 83 <= nb_spacegroup <= 88:
                    pointgroup = '4/m'
                if 89 <= nb_spacegroup <= 98:
                    if merge_friedel_pairs:
                        pointgroup = '4/mmm'
                    else:
                        pointgroup = '422'
                if 99 <= nb_spacegroup <= 110:
                    pointgroup = '4mm'
                if 111 <= nb_spacegroup <= 122:
                    pointgroup = '-42m'
                if 123 <= nb_spacegroup <= 142:
                    pointgroup = '4/m'

            elif 143 <= nb_spacegroup <= 167:
                system = 'trigonal'
                if 143 <= nb_spacegroup <= 146:
                    if merge_friedel_pairs:
                        pointgroup = '-3'
                    else:
                        pointgroup = '3'
                if 147 <= nb_spacegroup <= 148:
                    pointgroup = '-3'
                if 149 <= nb_spacegroup <= 155:
                    if merge_friedel_pairs:
                        pointgroup = '-31m'
                    else:
                        pointgroup = '312'
                if 156 <= nb_spacegroup <= 161:
                    pointgroup = '3m'
                if 162 <= nb_spacegroup <= 167:
                    pointgroup = '-3m'

            elif 168 <= nb_spacegroup <= 194:
                system = 'hexagonal'
                if 168 <= nb_spacegroup <= 173:
                    if merge_friedel_pairs:
                        pointgroup = '-6'
                    else:
                        pointgroup = '6'
                if nb_spacegroup <= 174:
                    pointgroup = '-6'
                if 175 <= nb_spacegroup <= 176:
                    pointgroup = '6/m'
                if 177 <= nb_spacegroup <= 182:
                    if merge_friedel_pairs:
                        pointgroup = '6/mmm'
                    else:
                        pointgroup = '622'
                if 183 <= nb_spacegroup <= 186:
                    pointgroup = '6mm'
                if 187 <= nb_spacegroup <= 190:
                    pointgroup = '-6m2'
                if 191 <= nb_spacegroup <= 194:
                    pointgroup = '6/m'

            elif 195 <= nb_spacegroup <= 230:
                system = 'cubic'
                if 195 <= nb_spacegroup <= 199:
                    if merge_friedel_pairs:
                        pointgroup = 'm-3'
                    else:
                        pointgroup = '23'
                if 200 <= nb_spacegroup <= 206:
                    pointgroup = 'm-3'
                if 207 <= nb_spacegroup <= 214:
                    if merge_friedel_pairs:
                        pointgroup = 'm-3m'
                    else:
                        pointgroup = '432'
                if 215 <= nb_spacegroup <= 220:
                    pointgroup = '-43m'
                if 221 <= nb_spacegroup <= 230:
                    pointgroup = 'm-3m'

            #get the default unique axis corresponding to the system ??? correct?
            if unique_axis==None:
                if system=='monoclinic':
                    unique_axis='b'
                if system=='hexagonal' or system=='tetragonal':
                    unique_axis='c'

            if unique_axis != None and (system == 'monoclinic' or system == 'hexagonal' or system == 'tetragonal'):
                pointgroup = pointgroup + '_ua' + unique_axis
            else:
                print('the unique axis will not be taken into account because the system is not appropriate')

            # 3. Checking unit cell correct for spacegroup, if not change according to system
            adapt_unit_cell_to_system=True
            if adapt_unit_cell_to_system:
#            if not crystal.symmetry.is_compatible_unit_cell(symmetry):  # the unit cell is not compatible with the spacegroup ???
                if system == 'monoclinic' and unique_axis == 'a':
                    if beta != 90:
                        beta = 90.0
                        print('the beta angle was changed from %s to 90 since the system is monoclinic' % (beta))
                    if gamma != 90:
                        print('the gamma angle was changed from %s to 90 since the system is monoclinic' % (gamma))
                        gamma = 90.0

                if system == 'monoclinic' and unique_axis == 'b':
                    if alpha != 90:
                        print('the beta angle was changed from %s to 90 since the system is monoclinic' % (alpha))
                        alpha = 90.0
                    if gamma != 90:
                        print('the gamma angle was changed from %s to 90 since the system is monoclinic' % (gamma))
                        gamma = 90.0

                if system == 'monoclinic' and (unique_axis == 'c' or unique_axis == None):
                    if alpha != 90:
                        print('the beta angle was changed from %s to 90 since the system is monoclinic' % (alpha))
                        alpha = 90.0
                    if beta != 90:
                        print('the beta angle was changed from %s to 90 since the system is monoclinic' % (beta))
                        beta = 90.0

                if system == 'orthorhombic':
                    if alpha != 90:
                        print('the beta angle was changed from %s to 90 since the system is monoclinic' % (alpha))
                        alpha = 90.0
                    if beta != 90:
                        print('the beta angle was changed from %s to 90 since the system is monoclinic' % (beta))
                        beta = 90.0
                    if gamma != 90:
                        print('the gamma angle was changed from %s to 90 since the system is monoclinic' % (gamma))
                        gamma = 90.0

                if system == 'cubic':
                    if alpha != 90:
                        print('the beta angle was changed from %s to 90 since the system is cubic' % (alpha))
                        alpha = 90.0
                    if beta != 90:
                        print('the beta angle was changed from %s to 90 since the system is cubic' % (beta))
                        beta = 90.0
                    if gamma != 90:
                        print('the gamma angle was changed from %s to 90 since the system is cubic' % (gamma))
                        gamma = 90.0
                    if a != b or a != c or b != c:
                        mean_abc = (a + b + c) / 3
                        print(
                            'a, b and c were changed from a=%s, b=%s, c=%s to a=%s, b=%s, c=%s since the system is cubic' % (
                                a, b, c, mean_abc, mean_abc, mean_abc))

                if system == 'hexagonal':
                    if alpha != 90:
                        print('the beta angle was changed from %s to 90 since the system is hexagonal' % (alpha))
                        alpha = 90.0
                    if beta != 90:
                        print('the beta angle was changed from %s to 90 since the system is hexagonal' % (beta))
                        beta = 90.0
                    if gamma != 120:
                        print('the gamma angle was changed from %s to 120 since the system is hexagonal' % (gamma))
                        gamma = 120.0
                    if unique_axis == 'c' or unique_axis == None:
                        if a != b:
                            mean_ab = (a + b) / 2
                            print(
                                'a anb b were changed from a=%s, b=%s to a=%s, b=%s since the system is trigonal or hexagonal' % (
                                    a, b, mean_ab, mean_ab))
                            a = mean_ab
                            b = mean_ab
                    if unique_axis == 'a':
                        if c != b:
                            mean_bc = (b + c) / 2
                            print(
                                'a anb b were changed from b=%s, c=%s to a=%s, b=%s since the system is trigonal or hexagonal' % (
                                    b, c, mean_bc, mean_bc))
                            b = mean_bc
                            c = mean_bc
                    if unique_axis == 'b':
                        if a != c:
                            mean_ac = (a + c) / 2
                            print(
                                'a anb b were changed from a=%s, c=%s to a=%s, b=%s since the system is trigonal or hexagonal' % (
                                    a, c, mean_ac, mean_ac))
                            a = mean_ac
                            c = mean_ac

                if system == 'trigonal':
                    if alpha != 90:
                        print('the beta angle was changed from %s to 90 since the system is trigonal' % (alpha))
                        alpha = 90.0
                    if beta != 90:
                        print('the beta angle was changed from %s to 90 since the system is trigonal' % (beta))
                        beta = 90.0
                    if a != b or a != c or b != c:
                        mean_abc = (a + b + c) / 3
                        print(
                            'a, b and c were changed from a=%s, b=%s, c=%s to a=%s, b=%s, c=%s since the system is trigonal' % (
                                a, b, c, mean_abc, mean_abc, mean_abc))

                if system == 'tetragonal':
                    if alpha != 90:
                        print('the beta angle was changed from %s to 90 since the system is tetragonal' % (alpha))
                        alpha = 90.0
                    if beta != 90:
                        print('the beta angle was changed from %s to 90 since the system is tetragonal' % (beta))
                        beta = 90.0
                    if gamma != 90:
                        print('the gamma angle was changed from %s to 90 since the system is tetragonal' % (gamma))
                        gamma = 90.0
                    if unique_axis == 'c' or unique_axis == None:
                        if a != b:
                            mean_ab = (a + b) / 2
                            print(
                                'a anb b were changed from a=%s, b=%s to a=%s, b=%s since the system is trigonal or hexagonal' % (
                                    a, b, mean_ab, mean_ab))
                            a = mean_ab
                            b = mean_ab
                    if unique_axis == 'a':
                        if c != b:
                            mean_bc = (b + c) / 2
                            print(
                                'a anb b were changed from b=%s, c=%s to a=%s, b=%s since the system is trigonal or hexagonal' % (
                                    b, c, mean_bc, mean_bc))
                            b = mean_bc
                            c = mean_bc
                    if unique_axis == 'b':
                        if a != c:
                            mean_ac = (a + c) / 2
                            print(
                                'a anb b were changed from a=%s, c=%s to a=%s, b=%s since the system is trigonal or hexagonal' % (
                                    a, c, mean_ac, mean_ac))
                            a = mean_ac
                            c = mean_ac

            unit_cell = [a, b, c, alpha, beta, gamma]

            print_in_T_and_log('space group = %s\n'
                                            'system = %s\n'
                                            'unique axis = %s\n'
                                            'point group = %s\n'
                                            'unit cell = %s' % (spacegroup, system, unique_axis, pointgroup, unit_cell))
            return (spacegroup, system, pointgroup, a, b, c, alpha, beta, gamma, unique_axis)

        self.spacegroup_off, self.system_off, self.pointgroup_off, self.a_off, self.b_off, self.c_off, self.alpha_off, self.beta_off, self.gamma_off, self.unique_axis_off = get_system_and_pointgroup(self.a_off, self.b_off, self.c_off, self.alpha_off, self.beta_off, self.gamma_off, self.spacegroup_off, self.unique_axis_off)

        # get nb of indexed images and number of images to use for JK
        self.n_frames_off = len(self.S_off.frames)  # search number of indexed images (frame is an indexed image)
        self.n_frames_to_keep_off = int(self.n_frames_off * self.fraction)

    # getting On state parameters
        if not self.JK_one_stream_file:
            self.stream_file_on = self.params.JackKnife.On_state.stream_file_on

            # check file
            Fextr_utils.check_all_files([self.stream_file_on])

            # read the stream file and get the number of indexed images, and the number of frames to keep
            print_in_T_and_log("---Reading stream file of the triggered state---")
            S_on = Stream(self.stream_file_on)
            self.stream_file_name_on = Fextr_utils.get_name(self.stream_file_on)
            self.a_on, a_array_on, astdev_on, self.b_on, b_array_on, bstdev_on, self.c_on, c_array_on, cstdev_on, self.alpha_on, alpha_array_on, alphastdev_on, self.beta_on, beta_array_on, betastdev_on, self.gamma_on, gamma_array_on, gammasdtev_on = Stream(
                self.stream_file_on).get_cell_stats()

            # Convert unit cell parameters (nm to A)
            self.a_on = self.a_on * 10
            self.b_on = self.b_on * 10
            self.c_on = self.c_on * 10

            # test for normal distribution
            print_in_T_and_log(' - Checking normal distribution of unit cell parameters - ')

            normal_distribution_test_unit_cell(a_array_on, 'a_on')
            normal_distribution_test_unit_cell(b_array_on, 'b_on')
            normal_distribution_test_unit_cell(c_array_on, 'c_on')
            normal_distribution_test_unit_cell(alpha_array_on, 'alpha_on')
            normal_distribution_test_unit_cell(beta_array_on, 'beta_on')
            normal_distribution_test_unit_cell(gamma_array_on, 'gamma_on')

            # getting unit cell parameters
            print_in_T_and_log('---Getting unit cell parameters---')
            if self.params.JackKnife.On_state.unit_cell.a != None:
                self.a_on = self.params.JackKnife.On_state.unit_cell.a
            else:
                print('Since a was not given, the average a will be taken')
            if self.params.JackKnife.On_state.unit_cell.b != None:
                self.b_on = self.params.JackKnife.On_state.unit_cell.b
            else:
                print('Since b was not given, the average b will be taken')
            if self.params.JackKnife.On_state.unit_cell.c != None:
                self.c_on = self.params.JackKnife.On_state.unit_cell.c
            else:
                print('Since c was not given, the average c will be taken')
            if self.params.JackKnife.On_state.unit_cell.alpha != None:
                self.alpha_on = self.params.JackKnife.On_state.unit_cell.alpha
            else:
                print('Since alpha was not given, the average alpha will be taken')
            if self.params.JackKnife.On_state.unit_cell.beta != None:
                self.beta_on = self.params.JackKnife.On_state.unit_cell.beta
            else:
                print('Since beta was not given, the average beta will be taken')
            if self.params.JackKnife.On_state.unit_cell.gamma != None:
                self.gamma_on = self.params.JackKnife.On_state.unit_cell.gamma
            else:
                print('Since gamma was not given, the average gamma will be taken')

            # get symmetry parameters
            spacegroup_on = self.params.JackKnife.On_state.spacegroup
            if spacegroup_on == None:
                print('Please give the space group for JackKnife')
                sys.exit()
            else:
                self.spacegroup_on = ''
                for i in spacegroup_on:
                    if i != ' ':
                        self.spacegroup_on += i

            self.unique_axis_on = self.params.JackKnife.On_state.unique_axis
            if self.unique_axis_on == None or self.unique_axis_on == 'default':
                self.unique_axis_on = None

            # get pointgroup and system
            self.spacegroup_on, self.system_on, self.pointgroup_on, self.a_on, self.b_on, self.c_on, self.alpha_on, self.beta_on, self.gamma_on, self.unique_axis_on = get_system_and_pointgroup(
                self.a_on, self.b_on, self.c_on, self.alpha_on, self.beta_on, self.gamma_on, self.spacegroup_on,
                self.unique_axis_on)

            # get nb of indexed images and number of images to use for JK
            self.n_frames_on = len(S_on.frames)  # search number of indexed images (frame is an indexed image)
            self.n_frames_to_keep_on = int(self.n_frames_on * self.fraction)


    # getting method for intensity merging, pointgroup and other parameters
        self.method_process_hkl = self.method_partialator = False

        if self.params.JackKnife.Scaling_and_merging.algorithm == ['process_hkl']:
            self.other_process_hkl = self.params.JackKnife.Scaling_and_merging.process_hkl.other_process_hkl
            self.other_partialator = None
            self.method_process_hkl = True

        elif self.params.JackKnife.Scaling_and_merging.algorithm == ['partialator']:
            self.other_partialator = self.params.JackKnife.Scaling_and_merging.partialator.other_partialator
            self.other_process_hkl = None
            self.method_partialator = True

        else:
            print("%s not defined. Monte Carlo method for intensity merging will be applied with no additional arguments." % (
                self.params.JackKnife.Scaling_and_merging.__phil_path__()))
            self.method_process_hkl = True
            self.other_process_hkl = None
            self.other_partialator = None

        if self.other_process_hkl==None: self.other_process_hkl=''
        if self.other_partialator==None: self.other_partialator=''

    #getting statistics
        self.other_stats_compare_hkl = self.params.JackKnife.Statistics.other_stats_compare_hkl
        if self.other_stats_compare_hkl == None:
            self.other_stats_compare_hkl = ''

    #getting low and high resolution if specified only if run_Xtrapol8 too #??? add to check_hkl too?
        self.lowres = self.highres = None
        if self.run_Xtrapol8:
            if self.params.Xtrapol8.input.low_resolution != None:
                self.lowres = self.params.Common_X8.low_resolution
                if not '--rmin' in self.other_stats_compare_hkl and not '--lowres' in self.other_stats_compare_hkl:
                    self.other_stats_compare_hkl += ' --lowres ' + str(self.lowres) #low resolution added to the other_stats_compare_hkl
            if self.params.Xtrapol8.input.high_resolution != None:
                self.highres = self.params.Common_X8.high_resolution
                if not '--rmax' in self.other_stats_compare_hkl and not '--highres' in self.other_stats_compare_hkl:
                    self.other_stats_compare_hkl += ' --highres ' + str(self.highres) #high resolution added to the other_stats_compare_hkl

#        return (
#            self.repeats, self.stream_file, self.stream_file_name, self.fraction, self.percentage, self.n_frames_to_keep, self.pointgroup,
#            self.other_process_hkl,
#            self.other_partialator, self.other_stats_compare_hkl, self.a, self.b, self.c, self.alpha, self.beta, self.gamma, self.system, self.unique_axis,
#            self.dir_cryst_prog,
#            self.method_process_hkl, self.method_partialator, self.spacegroup)

    def get_parameters_X8(self):
            # specify extrapolated structure factors and map types
            self.qFextr_map = self.qFgenick_map = self.qFextr_calc_map = self.Fextr_map = self.Fgenick_map = self.Fextr_calc_map = self.kFextr_map = self.kFgenick_map = self.kFextr_calc_map = False
            # calculate only the specified map types
            if 'qfextr' in self.params.Xtrapol8.f_and_maps.f_extrapolated_and_maps:
                self.qFextr_map = True
            if 'fextr' in self.params.Xtrapol8.f_and_maps.f_extrapolated_and_maps:
                self.Fextr_map = True
            if 'kfextr' in self.params.Xtrapol8.f_and_maps.f_extrapolated_and_maps:
                self.kFextr_map = True
            if 'qfgenick' in self.params.Xtrapol8.f_and_maps.f_extrapolated_and_maps:
                self.qFgenick_map = True
            if 'kfgenick' in self.params.Xtrapol8.f_and_maps.f_extrapolated_and_maps:
                self.kFgenick_map = True
            if 'fgenick' in self.params.Xtrapol8.f_and_maps.f_extrapolated_and_maps:
                self.Fgenick_map = True
            if ('qfextr_calc') in self.params.Xtrapol8.f_and_maps.f_extrapolated_and_maps:
                self.qFextr_calc_map = True
            if ('kfextr_calc') in self.params.Xtrapol8.f_and_maps.f_extrapolated_and_maps:
                self.kFextr_calc_map = True
            if ('fextr_calc') in self.params.Xtrapol8.f_and_maps.f_extrapolated_and_maps:
                self.Fextr_calc_map = True

            if self.params.Xtrapol8.f_and_maps.fofo_type == 'fofo':
                self.qFoFo_weight = False
                self.kFoFo_weight = False
            elif self.params.Xtrapol8.f_and_maps.fofo_type == 'qfofo':
                self.qFoFo_weight = True
                self.kFoFo_weight = False
            elif self.params.Xtrapol8.f_and_maps.fofo_type == 'kfofo':
                self.qFoFo_weight = False
                self.kFoFo_weight = True
            else:
                print_in_T_and_log("%s not defined. Q-weighting will be applied." % (self.params.Xtrapol8.f_and_maps.fofo_type.__phil_path__()))
                self.qFoFo_weight = True
                self.kFoFo_weight = False

            # Only run the non-weighted
            if self.params.Xtrapol8.f_and_maps.only_no_weight:
                self.qFextr_map = qFgenick_map = qFextr_calc_map = kFextr_map = kFgenick_map = kFextr_calc_map = False
                self.Fextr_map = Fgenick_map = Fextr_calc_map = True
            # Only run the k-weighted
            if self.params.Xtrapol8.f_and_maps.only_kweight:
                self.qFextr_map = qFgenick_map = qFextr_calc_map = Fextr_map = Fgenick_map = Fextr_calc_map = False
                self.kFextr_map = kFgenick_map = kFextr_calc_map = True
            # Only run the q-weighted
            if self.params.Xtrapol8.f_and_maps.only_qweight:
                self.Fextr_map = Fgenick_map = Fextr_calc_map = kFextr_map = kFgenick_map = kFextr_calc_map = False
                self.qFextr_map = qFgenick_map = qFextr_calc_map = True
            # Run all maps
            if self.params.Xtrapol8.f_and_maps.all_maps:  # calculate all Fextr map types
                self.qFextr_map = self.qFgenick_map = self.qFextr_calc_map = self.Fextr_map = self.Fgenick_map = self.Fextr_calc_map = self.kFextr_map = self.kFgenick_map = self.kFextr_calc_map = True

            # if all map types being false:
            if self.qFextr_map == self.qFgenick_map == self.qFextr_calc_map == self.Fextr_map == self.Fgenick_map == self.Fextr_calc_map == self.kFextr_map == self.kFgenick_map == self.kFextr_calc_map == False and self.params.Xtrapol8.f_and_maps.fast_and_furious == False:
                print_in_T_and_log(
                    'The combination of arguments used to define extrapolated structure factors and maps leads to no calculations at all. The default will be applied: qFextr')
                self.qFextr_map = True
            # if fast_and_furious mode: overwrite all F and map setting to default:
            if self.params.Xtrapol8.f_and_maps.fast_and_furious:
                # change parameters for Xtrapol8_out.phil
                self.params.Xtrapol8.f_and_maps.fofo_type = 'qfofo'
                self.params.Xtrapol8.f_and_maps.f_extrapolated_and_maps = ['qfextr']
                self.params.Xtrapol8.f_and_maps.only_no_weight = self.params.Xtrapol8.f_and_maps.all_maps = self.params.Xtrapol8.f_and_maps.only_kweight = self.params.Xtrapol8.f_and_maps.only_qweight = False
                # change working parameters
                self.qFoFo_weight = self.qFextr_map = True
                self.kFoFo_weight = self.qFgenick_map = self.qFextr_calc_map = self.Fextr_map = self.Fgenick_map = self.Fextr_calc_map = self.kFextr_map = self.kFgenick_map = self.kFextr_calc_map = False
                self.params.Xtrapol8.map_explorer.use_occupancy_from_distance_analysis = False
                self.params.Xtrapol8.f_and_maps.negative_and_missing = 'truncate_and_fill'
            #If JK is run, calm and curious method will be used
            if self.run_JackKnife:
                self.params.Xtrapol8.f_and_maps.fast_and_furious = False
                print_in_T_and_log('Jack Knife can not be run with Fast and Furious. Calm and Curious will be run with the few parameters of Fast and Furious.')

            # Bring all maptypes to be calculated together in list instead of using loose variables:
            self.all_maptypes = ['qFextr_map', 'Fextr_map', 'qFgenick_map', 'Fgenick_map', 'qFextr_calc_map',
                            'Fextr_calc_map',
                            'kFextr_map', 'kFgenick_map', 'kFextr_calc_map']
            self.all_maps = [self.qFextr_map, self.Fextr_map, self.qFgenick_map, self.Fgenick_map, self.qFextr_calc_map, self.Fextr_calc_map, self.kFextr_map,
                        self.kFgenick_map, self.kFextr_calc_map]
            self.maptypes_zip = zip(self.all_maptypes, self.all_maps)
            self.final_maptypes = [mp[0] for mp in self.maptypes_zip if mp[1] == True]

            print("final_maptypes", self.final_maptypes)

            if self.qFoFo_weight:
                self.params.Xtrapol8.f_and_maps.fofo_type = 'qfofo'
            elif self.kFoFo_weight:
                self.params.Xtrapol8.f_and_maps.fofo_type = 'kfofo'
            else:
                self.params.Xtrapol8.f_and_maps.fofo_type = 'fofo'

            # get list with occupancies from start-end-steps or list
            self.occ_step = (self.params.Xtrapol8.occupancies.high_occ - self.params.Xtrapol8.occupancies.low_occ) / self.params.Xtrapol8.occupancies.steps
            if self.params.Xtrapol8.occupancies.list_occ == None:
                self.occ = self.params.Xtrapol8.occupancies.low_occ
                self.occ_lst = []
                while self.occ <= self.params.Xtrapol8.occupancies.high_occ:
                    self.occ_lst.append(self.occ)
                    self.occ += self.occ_step
            else:
                self.occ_lst = self.params.Xtrapol8.occupancies.list_occ
            self.occ_lst.sort()
            if len(self.occ_lst) == 0:
                print("No input occupancies found. Xtrapol8 will stop after the FoFo calculation.")
                self.params.Xtrapol8.output.generate_fofo_only = True
                self.params.Xtrapol8.occupancies.list_occ = None
            else:
                self.params.Xtrapol8.occupancies.list_occ = self.occ_lst

    def get_parameters_JK_X8_output(self, params):
        '''
        Get the values of the output parameters from the phil file (for JK and/or X8)
        Args:
            params
        Returns:
            outname, output
        '''

        # get output name, or take prefix of triggered mtz or take dummy name if name is too long.
        if self.params.output.outname == None:
            if self.run_Xtrapol8 and not self.run_JackKnife:
                self.params.output.outname = Fextr_utils.get_name(params.Xtrapol8.input.triggered_mtz)
            else:
                self.params.output.outname = self.stream_file_name_on

        if len(self.params.output.outname) > 50:
            self.outname = 'triggered'
            print_in_T_and_log(
                "output.outname='%s' is too long. It will be substituted by '%s' during the excecution of Xtrapol8 and at the end the file names will be converted to the real output names. Please take this into account when inspecting log files." % (
                    self.params.output.outname, self.outname))
            print_in_T_and_log('---------------------------')
        else:
            self.outname = self.params.output.outname

        print('GETTING OUTPUT DIRECTORY\n==============================================================')
        #get output and change to output directory
        self.startdir = os.getcwd()
        if self.params.output.outdir != None and os.path.isdir(self.params.output.outdir):
            self.outdir = self.params.output.outdir
        elif self.params.output.outdir != None and os.path.exists(self.params.output.outdir)==False:
            try:
                os.mkdir(self.params.output.outdir)
                print('Output directory not present thus being created: %s' % (self.params.output.outdir))
            except OSError:
                os.makedirs(self.params.output.outdir)
            self.outdir = os.path.abspath(self.params.output.outdir)
        else:
            self.outdir = os.getcwd()

        if not self.JK_one_stream_file:
            outdir_0 = self.outdir + '/' + self.outname
            if self.run_JackKnife:
                outdir_0 = outdir_0 + '_JackKnife'
            if self.run_Xtrapol8:
                outdir_0 = outdir_0 + '_Xtrapol8'

            outdir_i=outdir_0
            i = 1
            while os.path.isdir(outdir_i):
                outdir_i = outdir_0 + '_' + str(i)
                i += 1
            self.outdir=outdir_i

        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)

    def get_parameters_DH_scaleit_X8(self, DH):

        # Extract the minimum and maximum resolution, as defined by the input parameters and/or the common reflections between the two scaled datasets
        self.dmax, self.dmin = DH.fobs_off_scaled.d_max_min()
        # set the params.resolution boundaries so that they appear in the output.phil
        self.params.Xtrapol8.input.high_resolution = self.dmin
        self.params.Xtrapol8.input.low_resolution = self.dmax

        return(self.dmin, self.dmax)














