"""
 panelIO.py

 Created 10/21/2020

--------
-authors and contact information
--------
-Elke De Zitter - elke.de-zitter@ibs.fr
-Nicolas Coquelle - nicolas.coquelle@esrf.fr
-Thomas Barends - Thomas.Barends@mpimf-heidelberg.mpg.de
-Jacques Philippe Colletier - jacques-Philippe.colletier@ibs.fr
--------
-license information
--------
-Copyright (c) 2021 Elke De Zitter, Nicolas Coquelle, Thomas Barends and Jacques-Philippe Colletier
-see https://github.com/ElkeDeZitter/Xtrapol8/blob/main/LICENSE
--------
"""

from matplotlib.backends.backend_wxagg import (
    FigureCanvasWxAgg as FigureCanvas,
    NavigationToolbar2WxAgg as NavigationToolbar,
)
from matplotlib.figure import Figure
import matplotlib.cm as cm
import matplotlib.colors as mcolors

import numpy as np
import wx
import os, re
from wx.lib.scrolledpanel import ScrolledPanel
from wxtbx import metallicbutton
from wx.lib.pubsub import pub
import pickle
import glob
import math
import sys
sys.path.append("..")
from Fextr_utils import get_name

script_dir = os.path.dirname(os.path.abspath(__file__))

class GradientButton (metallicbutton.MetallicButton):
    def __init__ (self, parent, label='', label2='', bmp=None,
                  size=wx.DefaultSize, style=metallicbutton.MB_STYLE_DEFAULT,
                  handler_function=None, user_data=None, start_color=(218,218,218),
                  gradient_percent=0, highlight_color=(160,190,210),
                  label_size=11, caption_size=8, button_margin=5,
                  disable_after_click=0, bmp2=None) :
        if (isinstance(bmp, str) or isinstance(bmp, unicode)):
            img = wx.Image(bmp, type=wx.BITMAP_TYPE_ANY, index=-1)
            #if size is not None :
            #        (w,h) = size
            #        img.Rescale(w,h)
            bmp = img.ConvertToBitmap()
        #    bmp = StandardBitmap(to_str(bmp))
        #if (isinstance(bmp2, str) or isinstance(bmp, unicode)):
        #    bmp2 = StandardBitmap(to_str(bmp2))
        if (user_data is None):
            user_data = 'None'
        metallicbutton.MetallicButton.__init__(self,
                                               parent=parent,
                                               label=label,
                                               label2=label2,
                                               bmp=bmp,
                                               size=size,
                                               style=style,
                                               name=user_data,     # AutoBuild can pass unicode path
                                               start_color=start_color,
                                               gradient_percent=gradient_percent,
                                               highlight_color=highlight_color,
                                               label_size=label_size,
                                               caption_size=caption_size,
                                               button_margin=button_margin,
                                               disable_after_click=disable_after_click,
                                               bmp2=bmp2)
        self.Enable(False)

    def OnEnter(self, evt):
        return

    def OnLeave(self,evt):
      return

    def OnHIGH(self, evt):
        self.SetState(2)

    def OnNormal(self, evt):
      self.SetState(0)

    def OnFocus(self, evt):
      return


class TabLog(wx.Panel):
    """
    This will be the first tab of the Results Notebook
    (which holds the gui and the different processing steps)
    """
    # ----------------------------------------------------------------------
    LOG_ENTRIES_CNC = ["Data Preparation",
                   "Fourier Difference Map Calculation",
                   "Fourier Difference Map Analysis",
                   "Fextrapolated Map computation",
                   "Optimal occupancy estimation"]

    LOG_ENTRIES_FNF = ["Data Preparation",
                   "Fourier Difference Map Calculation",
                   "Fourier Difference Map Analysis",
                   "Fextrapolated Map computation",
                   "Optimal occupancy estimation",
                   "Refinements",]

    LOG_STEPS_INIT = ["DATA PREPARATION",
                      "CALCULATE (WEIGHTED) FO-FO FOURIER DIFFERENCE MAP",
                      "************Map explorer************"]

    LOG_IDX_INIT = [0, 1, 2]
    LOG_STEPS_CORE = ["type of structure factors and maps for occupancy"]
    LOG_IDX_CORE_CNC = [3]
    LOG_IDX_CORE_FNF = [3]

    LOG_STEPS_OCC = ["ESTIMATE OPTIMAL OCCUPANCY"]
    LOG_IDX_OCC_CNC = [4]
    LOG_IDX_OCC_FNF = [4]

    LOG_STEPS_REF_FNF = ["FAST AND FURIOUS REFINEMENT"]
    LOG_IDX_REF_FNF = [5]

    def __init__(self, parent, options):

        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

        self.count = 0
        self.options = options
        self.fNf = self.options.f_and_maps.fast_and_furious
        self.FoOnly = self.options.output.generate_fofo_only
        self.indexLOG = 0
        self.LOG_STEPS = []
        self.LOG_IDX = []
        self.buildLogSteps()
        self.buttons = []

        for i in range(len(self.LOG_ENTRIES)):
            label2 = ''
            if i == 3:
                label2 = ' '
            self.buttons.append(GradientButton(self, label=self.LOG_ENTRIES[i], label2=label2))

        self.ButtonSizer = wx.BoxSizer(wx.VERTICAL)
        for button in self.buttons:
            self.ButtonSizer.Add(button, proportion=0, flag=wx.EXPAND | wx.ALL, border=10)
            self.ButtonSizer.AddSpacer(10)
        self.mainSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.LogTextCtrl = wx.TextCtrl(self, wx.ID_ANY, size=(800, 600),
                  style=wx.TE_MULTILINE | wx.TE_READONLY)

        self.mainSizer.Add(self.LogTextCtrl)
        self.mainSizer.Add(self.ButtonSizer, proportion=1, flag=wx.EXPAND, border=10)
        self.SetSizer(self.mainSizer)

    def buildLogSteps(self):

        try:
            N_occ = len(self.options.occupancies.list_occ)
            self.occ = self.options.occupancies.list_occ[0]
        except TypeError:
            occs = np.linspace(self.options.occupancies.low_occ,
                                      self.options.occupancies.high_occ,
                                      self.options.occupancies.steps+1, endpoint=True)
            N_occ = np.size(occs)
            self.occ = occs[0]

        self.LOG_STEPS.extend(self.LOG_STEPS_INIT)
        self.LOG_IDX.extend(self.LOG_IDX_INIT)

        if self.FoOnly:
            qftype = 0
            self.qftype = ''
            self.LOG_IDX_CORE = []
            self.LOG_ENTRIES = self.LOG_ENTRIES_CNC[0:3]

        elif self.fNf:
            self.qftype='qFextr'
            qftype = 1
            self.LOG_IDX_CORE = self.LOG_IDX_CORE_FNF
            self.LOG_ENTRIES = self.LOG_ENTRIES_FNF

        else:
            qftype = len(self.options.f_and_maps.f_extrapolated_and_maps)
            self.qftype = self.options.f_and_maps.f_extrapolated_and_maps
            self.LOG_IDX_CORE = self.LOG_IDX_CORE_CNC
            self.LOG_ENTRIES = self.LOG_ENTRIES_CNC

        N = qftype * N_occ
        self.LOG_STEPS.extend(self.LOG_STEPS_CORE * N)
        self.LOG_IDX.extend(self.LOG_IDX_CORE * N)

        if not self.FoOnly:
            self.LOG_STEPS.extend(self.LOG_STEPS_OCC)
            if self.fNf:
                self.LOG_IDX.extend(self.LOG_IDX_OCC_FNF)
                self.LOG_STEPS.extend(self.LOG_STEPS_REF_FNF)
                self.LOG_IDX.extend(self.LOG_IDX_REF_FNF)
            else:
                self.LOG_IDX.extend(self.LOG_IDX_OCC_CNC)

    def updateLog(self, line):
        self.LogTextCtrl.WriteText(line)
        try:

            if self.LOG_STEPS[self.indexLOG] in line:
                if self.indexLOG > 0:
                    self.buttons[self.LOG_IDX[self.indexLOG - 1]].OnNormal(None)
                self.buttons[self.LOG_IDX[self.indexLOG]].OnHIGH(None)
                if self.LOG_STEPS_CORE[0] in line:
                    qftype = line.split()[1]
                    occ = line.split()[10]
                    caption = '%s - occ %s' % (qftype, occ)
                    self.buttons[self.LOG_IDX[self.indexLOG]]._label2 = caption

                if self.indexLOG < len(self.LOG_STEPS) - 1: self.indexLOG += 1

        except IndexError:
            pass

    def CreateCoot(self):
        Button2 = wx.Button(self, wx.ID_ANY, label="Open in Coot")
        bmp = wx.Image(os.path.join(script_dir, 'pngs/coot_logo_small.png'), type=wx.BITMAP_TYPE_ANY,
                       index=-1).ConvertToBitmap()

        # if size is not None :
        #        (w,h) = size
        #        img.Rescale(w,h)
        Button2.SetBitmapMargins((2, 2))
        Button2.SetBitmap(bmp)
        Button2.Bind(wx.EVT_BUTTON, self.onCoot)

        self.ButtonSizer.Add(Button2, proportion=0, flag=wx.EXPAND | wx.ALL, border=10)
        self.mainSizer.Layout()

    def onCoot(self, evt):
        script_coot = os.path.join(self.options.output.outdir,'coot_all_qFo-Fo.py')
        os.system("coot --script %s &" % (script_coot))


class TabMainImg(ScrolledPanel):
    """
    """
    # ----------------------------------------------------------------------
    def __init__(self, parent):
        ScrolledPanel.__init__(self, parent=parent, style=wx.VSCROLL | wx.HSCROLL)
        self.SetupScrolling()
        self.parent = parent
        self.mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.plotSizer = wx.BoxSizer(wx.VERTICAL)
        self.mainSizer.Add(self.plotSizer, 1, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL)
        self.SetSizer(self.mainSizer)
        self.SetAutoLayout(1)
        self.photoMaxSize = 1000
        self.mapping = {
            "Riso_CCiso.pickle": self.plot_Riso_CCiso,
            "q_estimation.pickle": self.plot_q_estimation,
            "k_estimation.pickle": self.plot_k_estimation,
            "Fextr_negative.pickle": self.plot_Neg_Pos_reflections,
            "Fextr_binstats.pickle": self.plot_sigmas
        }
        for fextr in ['Fextr',  'qFextr', 'kFextr',
                      'Fgenick', 'qFgenick', 'kFgenick',
                      'Fextr_calc', 'qFextr_calc', 'kFextr_calc']:
            self.mapping["alpha_occupancy_determination_%s.pickle"  % fextr] = self.plot_alpha_occupancy_determination
            self.mapping["%s_refinement_R-factors_per_alpha.pickle" % fextr] = self.plot_refinement_Rfactors_per_alpha




    def addPlot(self, pickle_file):

        _, pickle_name = os.path.split(pickle_file)
        canvas = self.mapping[pickle_name](pickle_file)
        self.plotSizer.Add(canvas, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_CENTER_HORIZONTAL)

        toolbar = NavigationToolbar(canvas)
        toolbar.Realize()
        # By adding toolbar in sizer, we are able to put it at the bottom
        # of the frame - so appearance is closer to GTK version.
        self.plotSizer.Add(toolbar, 0, wx.LEFT | wx.EXPAND)

        # update the axes menu on the toolbar
        toolbar.update()
        self.plotSizer.AddSpacer(60)
        #self.SetSizer(self.sizer)
        self.Fit()

    def plot_k_estimation(self, pickle_file):
        with open(pickle_file, 'rb') as stats_file:
            bin_res_cent_lst, k_av_lst, k_max_lst, k_min_lst = pickle.load(stats_file)

        self.figure = Figure(figsize=(10, 5))
        ax1 = self.figure.add_subplot(111)
        ax1.set_xlabel('Resolution (A)')
        ax1.set_ylabel('Average k-weight in resolution bin')
        ax1.plot(bin_res_cent_lst[:], k_av_lst[:], marker='.', label='Average k', color='red')
        ax1.tick_params(axis='y')
        ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
        ax1.set_ylim(0, 1)
        ax2 = ax1.twinx()
        ax2.fill_between(bin_res_cent_lst[:], k_max_lst[:], k_min_lst[:], color='red', alpha=0.2, label='k range')
        ax2.set_ylim(0, 1)
        ax2.tick_params(axis='y')
        ax2.set_ylabel('k range within resolution bin')
        lines_labels = [ax.get_legend_handles_labels() for ax in self.figure.axes]
        lines, labels = [sum(lne, []) for lne in zip(*lines_labels)]

        ax2.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.75, -0.05, 0.45, 0.5), fontsize='small',
                   framealpha=0.5)
        #self.figure.tight_layout()
        ax1.set_title('Average k for high resolution reflections', fontsize='medium', fontweight="bold")
        self.figure.subplots_adjust(hspace=0.35, left=0.09, right=0.82, top=0.95)
        canvas = FigureCanvas(self, -1, self.figure)
        return canvas

    def plot_Riso_CCiso(self, pickle_file):
        with open(pickle_file, 'rb') as stats_file:
            bin_res_cent_lst, r_work_lst, cc_work_lst, r_work, cc_work = pickle.load(stats_file)

        self.figure = Figure(figsize=(10, 5))#, tight_layout=True)#, dpi=100)
        ax1 = self.figure.add_subplot(111)
        ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
        ax1.set_xlabel('Resolution (A)')
        ax1.set_ylabel('Riso')
        ax1.yaxis.label.set_color('red')
        ax1.plot(bin_res_cent_lst[:], r_work_lst[:], marker='.', color='red', linewidth=2,
                 label='Riso; overall %.4f' % (r_work))

        ax2 = ax1.twinx()
        ax2.plot(bin_res_cent_lst[:], cc_work_lst[:],  marker = '^', markersize = 5, color = 'green', linewidth=2,
                     label='CCiso; overall %.4f' % (cc_work))
        ax2.set_ylabel('CCiso')
        ax2.yaxis.label.set_color('green')
        lines_labels = [ax.get_legend_handles_labels() for ax in self.figure.axes]
        lines, labels = [sum(lne, []) for lne in zip(*lines_labels)]
        ax2.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.81, -0.05, 0.45, 0.5), fontsize='small',
                       framealpha=0.5)
        ax1.set_title('Riso and CCiso for high resolution reflections', fontsize='medium', fontweight="bold")
        self.figure.subplots_adjust(hspace=0.35, left=0.09, right=0.82, top=0.95)
        canvas = FigureCanvas(self, -1, self.figure)
        return canvas

    def plot_q_estimation(self, pickle_file="q_estimation.pickle"):
        """
        For k_estimation plot
        Generated only once
        """
        if os.path.isfile(pickle_file) == False:
            return

        try:
            with open(pickle_file, 'rb') as stats_file:
                bin_res_cent_lst, q_av_lst, q_max_lst, q_min_lst = pickle.load(stats_file)
        except:
            with open(pickle_file, 'rb') as stats_file:
                bin_res_cent_lst, q_max_lst, q_min_lst = pickle.load(stats_file)
                q_av_lst = None


        self.figure = Figure(figsize=(10, 5))#, tight_layout=True)
        ax1 = self.figure.add_subplot(111)
        ax1.set_xlabel('Resolution (A)')
        ax1.set_ylabel('Average q in resolution bin')
        if q_av_lst is not None: ax1.plot(bin_res_cent_lst[:], q_av_lst[:], marker='.', label='Average q', color='red')
        ax1.tick_params(axis='y')
        ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
        ax1.set_ylim(0, 1)
        ax2 = ax1.twinx()
        ax2.fill_between(bin_res_cent_lst[:], q_max_lst[:], q_min_lst[:], color='red', alpha=0.2, label='q range')
        ax2.set_ylim(0, 1)
        ax2.tick_params(axis='y')
        ax2.set_ylabel('q range within resolution bin')
        lines_labels = [ax.get_legend_handles_labels() for ax in self.figure.axes]
        lines, labels = [sum(lne, []) for lne in zip(*lines_labels)]
        ax2.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.75, -0.05, 0.45, 0.5), fontsize='small',
                   framealpha=0.5)
        ax1.set_title('Average q for high resolution reflections', fontsize='medium', fontweight="bold")
        self.figure.subplots_adjust(hspace=0.35, left=0.09, right=0.82, top=0.95)
        canvas = FigureCanvas(self, -1, self.figure)
        return canvas

    def plot_Neg_Pos_reflections(self, suffix, pickle_file):
        """
        For Neg_Pos_reflections plot.
        This plot should be prepared for each type of extrapolated structure factor
        This plot should be updated during run
        The input pickle file contains information on all types of extrapolation (the name is thus misleading)!!!!
        In this script, we reed the pickle and search specifically for maptype (qFextr, kFextrm, Fextr, qFgenick, kFgenick, Fgenick, qFextr_calc, kFextr_calc, Fextr_calc) and make the plot. In that sense it differs in the way the plots are generated in the command-line version
        """
        # In Fextr_utils, "maptype" is used. In order to keep the plotting code the same, let's use the same name
        maptype = suffix
        if os.path.isfile(pickle_file) == False:
            return

        with open(pickle_file, 'rb') as stats_file:
            while True:
                try:
                    stats = np.array(pickle.load(stats_file))
                    # stats = np.array(tuple(stats), dtype='f8, S32, i4, i4, i4')
                    try:
                        alldata = np.vstack([alldata, stats[np.newaxis, ...]])
                    except NameError:
                        alldata = stats[np.newaxis, ...]
                except EOFError:
                    break

        indices = np.where(alldata[:, 1] == maptype)[0]

        self.figure = Figure(figsize=(10, 5))
        ax0, ax2 = self.figure.subplots(1,2)
        #fig, (ax0, ax2) = plt.subplots(1, 2, figsize=(10, 5))
        ax1 = ax0.twinx()
        ax3 = ax2.twinx()

        for a in indices:
            ax0.plot(np.float(alldata[a][0]), np.float(alldata[a][-1]), marker='o', color='red', label=maptype)
            ax2.plot(np.float(alldata[a][0]), np.float(alldata[a][3]), marker='o', color='red', label=maptype)

        ax0.set_ylim(0, np.float(alldata[a][2]))
        ax1.set_ylim(0, 100)
        ax2.set_ylim(0, np.float(alldata[a][2]))
        ax3.set_ylim(0, 100)

        ax0.set_xlabel('Occupancy')
        ax0.set_ylabel('Absolute number')
        ax1.set_ylabel('Percentage of total number of reflections')
        ax0.set_title('%s: Negative reflections' % (maptype), fontsize='medium', fontweight="bold")

        ax2.set_xlabel('Occupancy')
        ax2.set_ylabel('Absolute number')
        ax3.set_ylabel('Percentage of total number of reflections')
        ax2.set_title('%s: Positive reflections' % (maptype), fontsize='medium', fontweight="bold")

        self.figure.subplots_adjust(hspace=0.25, wspace=0.5, left=0.09, right=0.88, top=0.95)
        canvas = FigureCanvas(self, -1, self.figure)
        return canvas

    def plot_refinement_Rfactors_per_alpha(self, prefix, pickle_file):
        """
        For refinement_R-factors_per_alpha plot.
        This plot should be prepared for each type of extrapolated structure factor
        prefix can be qFextr, kFextrm, Fextr, qFgenick, kFgenick, Fgenick, qFextr_calc, kFextr_calc, Fextr_calc
        """
        # In Fextr_utils, "maptype" is used. In order to keep the plotting code the same, let's use the same name
        maptype = prefix

        #pickle_file = '%s_refinement_R-factors_per_alpha.pickle' % (prefix)
        if os.path.isfile(pickle_file) == False:
            return

        with open(pickle_file, 'rb') as stats_file:
            occ_lst, r_work_lst, r_free_lst, r_diff_lst = pickle.load(stats_file)

        self.figure = Figure(figsize=(10, 5))
        ax0 = self.figure.add_subplot(111)
        ax1 = ax0.twinx()

        ax0.plot(occ_lst, r_work_lst, color='red', marker='o', label='Rwork')
        ax0.plot(occ_lst, r_free_lst, color='blue', marker = 's', markersize = 5, label='Rfree')
        ax1.plot(occ_lst, r_diff_lst, color='green', marker = '^', label='Rfree-Rwork')
        ax0.set_xlabel('Occupancy of triggered state')
        ax0.set_ylabel('R-factor')
        ax1.set_ylabel('R-factor difference')
        lines_labels = [ax.get_legend_handles_labels() for ax in self.figure.axes]
        lines, labels = [sum(lne, []) for lne in zip(*lines_labels)]
        ax1.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.75, -0.05, 0.45, 0.5), fontsize='x-small',
                   framealpha=0.5)

        ax0.set_title('R-factors after reciprocal space refinement with %s' % (maptype), fontsize='medium',
                  fontweight="bold")
        self.figure.subplots_adjust(hspace=0.35, left=0.09, right=0.82, top=0.95)
        canvas = FigureCanvas(self, -1, self.figure)
        return canvas


    def plot_alpha_occupancy_determination(self, suffix, pickle_file):
        """
        For alpha_occupancy_determination plot.
        This plot should be prepared for each type of extrapolated structure factor
        prefix can be qFextr, kFextrm, Fextr, qFgenick, kFgenick, Fgenick, qFextr_calc, kFextr_calc, Fextr_calc
        """
        #pickle_file = 'alpha_occupancy_determination_%s.pickle' % (suffix)
        if not os.path.isfile(pickle_file):
            return

        with open(pickle_file, 'rb') as stats_file:
            alphas, occupancies, int1_norm, int2_norm, results, resids_lst_used, alphafound = pickle.load(stats_file)

        self.figure = Figure(figsize=(10, 5))
        ax1, ax2 = self.figure.subplots(1, 2)

        if resids_lst_used == False:
            ax1.plot(alphas, int1_norm, 's', markersize = 5, color = 'blue',
                     label='All peaks')  # int1, int2 and conv will be the same, so we can plot only int1 and avoid the try catch
        else:
            try:
                ax1.plot(alphas, int1_norm, 'o', color = 'red', label='Selected residues')
                ax1.plot(alphas, int2_norm, 's', markersize = 5, color = 'blue', label='All peaks')
                ax1.plot(alphas, results, '^', color="green", label='Selected residues with enhanced SNR')
            except:
                ax1.plot(alphas, int1_norm, 'o', color = 'red')

        ax1.set_ylim([0., 1.1])
        ax1.set_xlim([np.min(alphas) * 0.95, np.max(alphas) * 1.05])
        ax1.set_xlabel('Alpha value = 1/occupancy')
        ax1.set_ylabel('Normalized difference map ratio')

        if resids_lst_used == False:
            ax2.plot(occupancies, int1_norm, 's', markersize = 5, color = 'blue',
                     label='All peaks')  # int1, int2 and conv will be the same, so we can plot only int1 and avoid the try catch
        else:
            try:
                ax2.plot(occupancies, int1_norm, 'o', color = 'red', label='Selected residues')
                ax2.plot(occupancies, int2_norm, 's', markersize = 5, color = 'blue', label='All peaks')
                ax2.plot(occupancies, results, '^', color="green", label='Selected residues with enhanced SNR')
            except:
                ax2.plot(occupancies, int1_norm, 'o', color = 'red')

        ax2.set_ylim([0., 1.1])
        ax2.set_xlim([np.min(occupancies) * 0.95, np.max(occupancies) * 1.05])
        ax2.set_xlabel('Triggered state occupancy')
        ax2.set_ylabel('Normalized difference map ratio')
        ax2.legend(loc='lower right', bbox_to_anchor=(0.92, -0.05, 0.45, 0.5), fontsize='x-small', framealpha=0.5)

        if alphafound:
            ax1.set_title('Alpha determination', fontsize='medium', fontweight="bold")
            ax2.set_title('Occupancy determination', fontsize='medium', fontweight="bold")
        else:
            ax1.set_title('Alpha determination IMPOSSIBLE', fontsize='medium', fontweight="bold")
            ax1.text(np.min(alphas), 0.5, 'no peaks found in at least one of the maps')
            ax2.set_title('Occupancy determination IMPOSSIBLE', fontsize='medium', fontweight="bold")
            ax2.text(np.min(occupancies), 0.5, 'no peaks found in at least one of the maps')
        self.figure.subplots_adjust(hspace=0.25, wspace=0.4, left=0.09, right=0.88, top=0.95)
        canvas = FigureCanvas(self, -1, self.figure)
        return canvas

    def addImg(self, filepath):
        img = wx.Image(filepath, wx.BITMAP_TYPE_ANY)
        W = img.GetWidth()
        H = img.GetHeight()
        if W > H:
            NewW = self.photoMaxSize
            NewH = self.photoMaxSize * H / W
        else:
            NewH = self.photoMaxSize
            NewW = self.photoMaxSize * W / H
        img = img.Scale(NewW, NewH)

        self.newimg = wx.StaticBitmap(self, wx.ID_ANY,
                                      wx.BitmapFromImage(img))
        self.mainSizer.Add(self.newimg, proportion=0,  flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.mainSizer.Add(wx.StaticLine(self, wx.ID_ANY))
        self.mainSizer.AddSpacer(60)
        self.FitInside()

    def addChoices(self, selection):
        self.FextrSelection = wx.Choice(self, wx.ID_ANY, choices=selection)
        self.mainSizer.Add(self.FextrSelection, 0, wx.ALIGN_CENTER)
        self.FextrSelection.SetSelection(0)
        self.ImgSizer = wx.BoxSizer(wx.VERTICAL)

        self.mainSizer.AddSpacer(30)
        self.mainSizer.Add(self.ImgSizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.FitInside()
        index = self.parent.GetSelection()
        pub.sendMessage("updateFextr", evt=None)

    def addFextrPlot(self, fextr, pickle_file):
        _, pickle_name = os.path.split(pickle_file)
        canvas = self.mapping[pickle_name](fextr, pickle_file)
        self.ImgSizer.Add(canvas, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_CENTER_HORIZONTAL)

        toolbar = NavigationToolbar(canvas)
        toolbar.Realize()
        # By adding toolbar in sizer, we are able to put it at the bottom
        # of the frame - so appearance is closer to GTK version.
        self.ImgSizer.Add(toolbar, 0, wx.LEFT | wx.EXPAND)

        # update the axes menu on the toolbar
        toolbar.update()
        self.ImgSizer.AddSpacer(60)
        #self.SetSizer(self.sizer)
        self.FitInside()

    def plot_sigmas(self, prefix=None, pickle_file='Fextr_binstats.pickle'):
        if prefix.endswith('pickle'):
            return self.plot_FoFosigmas(pickle_file=prefix)
        else:
            return self.plot_Fextrsigmas(prefix=prefix, pickle_file=pickle_file)

    def plot_Fextrsigmas(self, prefix, pickle_file='Fextr_binstats.pickle'):
        """
        For Fextr_sigmas plot.
        This plot should be prepared for each type of extrapolated structure factor
        This plot should be updated during run
        The input pickle file contains information on all types of extrapolation (the name is thus misleading)!!!!
        In this script, we reed the pickle and search specifically for maptype (qFextr, kFextrm, Fextr, qFgenick, kFgenick, Fgenick, qFextr_calc, kFextr_calc, Fextr_calc) and make the plot. In that sense it differs in the way the plots are generated in the command-line version
        """
        # In Fextr_utils, "maptype" is used. In order to keep the plotting code the same, let's use the same name
        maptype = prefix
        # Read the pickle file
        with open(pickle_file, 'rb') as stats_file:
            while True:
                try:
                    stats = np.array(pickle.load(stats_file))
                    # stats = np.array(tuple(stats), dtype='f8, S32, i4, i4, i4')
                    try:
                        alldata = np.vstack([alldata, stats[np.newaxis, ...]])
                    except NameError:
                        alldata = stats[np.newaxis, ...]
                except EOFError:
                    break

        # extract the occupancies
        occ_lst = list(set(alldata[:, 0]))
        # get the indices concerning the specific maptype we are looking at
        indices = np.where(alldata[:, 2] == maptype)[0]

        mn = 0
        mx = 0
        self.figure = Figure(figsize=(10, 5))
        ax0 = self.figure.add_subplot(111)
        #fig, ax0 = plt.subplots(1, 1, figsize=(10, 5))
        ax1 = ax0.twinx()
        for a in indices:
            occ, _, _, bin_res_cent_lst, fextr_data_lst, fextr_sigmas_lst, _, _ = alldata[a]
            color_data = cm.Reds(int((np.log(1 / occ)) * 100))
            color_sigma = cm.Blues(int((np.log(1 / occ)) * 100))
            ax0.plot(bin_res_cent_lst[:], fextr_data_lst[:], marker='.', color=color_data,
                     label='%s, occ = %.3f' % (maptype, occ))
            ax1.plot(bin_res_cent_lst[:], fextr_sigmas_lst[:], marker='x', linestyle='--', color=color_sigma,
                     label='sig(%s), occ = %.3f' % (maptype, occ))
            # Specify the minimum and maximum value
            if min(fextr_data_lst[1:]) < mn:
                mn = min(fextr_data_lst[1:])
            if max(fextr_data_lst[1:]) > mx:
                mx = max(fextr_data_lst[1:])
            if min(fextr_sigmas_lst[1:]) < mn:
                mn = min(fextr_sigmas_lst[1:])
            if max(fextr_sigmas_lst[1:]) > mx:
                mx = max(fextr_sigmas_lst[1:])

        ax0.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
        ax0.set_xlabel('Resolution (A)')  # , fontsize = 'small')
        ax0.set_ylabel('ESFAs')  # , fontsize = 'small')
        ax0.yaxis.label.set_color('tab:red')

        ax0.set_ylim(mn, mx)
        ax1.set_ylim(mn, mx)

        ax0.legend(loc='lower right', bbox_to_anchor=(0.89, -0.05, 0.45, 0.5), fontsize='x-small', framealpha=0.5)
        ax1.legend(loc='lower right', bbox_to_anchor=(1.17, -0.05, 0.45, 0.5), fontsize='x-small', framealpha=0.5)
        ax1.set_ylabel('sig(ESFAs)')  # , fontsize = 'small')
        ax1.yaxis.label.set_color('tab:blue')

        ax0.set_title('%s for high resolution reflections' % (maptype), fontsize='medium', fontweight="bold")
        self.figure.subplots_adjust(hspace=0.35, left=0.09, right=0.65, top=0.95)
        canvas = FigureCanvas(self, -1, self.figure)
        return canvas

        #plt.savefig('%s_sigmas.pdf' % (maptype), dpi=300, transparent=True)
        #plt.savefig('%s_sigmas.png' % (maptype), dpi=300)
        #plt.close()

    def plot_FoFosigmas(self, pickle_file='Fextr_binstats.pickle'):
        """
        For FoFo-sigmas plot.
        Generated only once
        In this script, we reed the pickle and search only the info concerning the FoFo. In that sense it differs in the way the plots are generated in the command-line version
        """
        if os.path.isfile(pickle_file) == False:
            return

        with open(pickle_file, 'rb') as stats_file:
            _, FoFo_type, _, bin_res_cent_lst, _, _, fdif_data_lst, fdif_sigmas_lst = np.array(pickle.load(stats_file))

        mn = 0
        mx = 0
        if min(fdif_data_lst[1:]) < mn:
            mn = min(fdif_data_lst[1:])
        if max(fdif_data_lst[1:]) > mx:
            mx = max(fdif_data_lst[1:])
        if min(fdif_sigmas_lst[1:]) < mn:
            mn = min(fdif_sigmas_lst[1:])
        if max(fdif_sigmas_lst[1:]) > mx:
            mx = max(fdif_sigmas_lst[1:])

        self.figure = Figure(figsize=(10, 5))
        ax0 = self.figure.add_subplot(111)

        ax1 = ax0.twinx()

        ax0.plot(bin_res_cent_lst[:], fdif_data_lst[:], marker='.', color='tab:red', label='%s' % (FoFo_type))
        ax1.plot(bin_res_cent_lst[:], fdif_sigmas_lst[:], marker='x', linestyle='--', color='tab:blue',
                 label='sig(%s)' % (FoFo_type))
        ax0.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
        ax0.set_xlabel('Resolution (A)')
        ax0.set_ylabel('%s' % (FoFo_type))
        ax0.yaxis.label.set_color('tab:red')
        ax1.set_ylabel('sig(%s)' % (FoFo_type))
        ax1.yaxis.label.set_color('tab:blue')
        lines_labels_1 = [ax.get_legend_handles_labels() for ax in [ax0, ax1]]
        lines_1, labels_1 = [sum(lne, []) for lne in zip(*lines_labels_1)]
        ax0.legend(lines_1, labels_1, loc='lower right', bbox_to_anchor=(0.73, -0.05, 0.45, 0.5), fontsize='xx-small',
                   framealpha=0.5)
        ax0.set_title('%s for high resolution reflections' % (FoFo_type), fontsize='medium', fontweight="bold")
        self.figure.subplots_adjust(hspace=0.35, left=0.09, right=0.85, top=0.95)
        canvas = FigureCanvas(self, -1, self.figure)
        return canvas

    def addFextrImg(self, filepath):
        img = wx.Image(filepath, wx.BITMAP_TYPE_ANY)
        W = img.GetWidth()
        H = img.GetHeight()
        if W > H:
            NewW = self.photoMaxSize
            NewH = self.photoMaxSize * H / W
        else:
            NewH = self.photoMaxSize
            NewW = self.photoMaxSize * W / H
        img = img.Scale(NewW, NewH)

        self.newimg = wx.StaticBitmap(self, wx.ID_ANY, wx.BitmapFromImage(img))
        self.ImgSizer.Add(self.newimg, proportion=0,  flag=wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, border=20)
        self.FitInside()

    def Clear(self, evt):
        while not self.ImgSizer.IsEmpty():
            item = self.ImgSizer.GetItem(0)
            if item.IsSpacer:
                self.ImgSizer.Detach(0)
            else:
                window = self.ImgSizer.GetItem(0).GetWindow()
                window.Destroy()
        pub.sendMessage("updateFextr", evt=None)#,tabindex=index)

class TabOccResults(ScrolledPanel):
    """
     This will be the first tab of the Configure Notebook (which holds all X8 inputs)
     """

    # ----------------------------------------------------------------------
    def __init__(self, parent, options):
        ScrolledPanel.__init__(self, parent=parent, style=wx.VSCROLL | wx.HSCROLL)
        self.SetupScrolling()
        self.mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.photoMaxSize = 1000
        self.options = options
        self.finished = False
        self.best_occ = {}
        self.ddm = {}
        self.coot_scripts = {}
        self.setupUI()

        self.setupBindings()


    def setupUI(self):
        try:
            self.occ_list = [str(choice) for choice in self.options.occupancies.list_occ]

        except TypeError:
            self.occ_list = ['%.3f' %choice for choice in np.linspace(self.options.occupancies.low_occ,
                                                                      self.options.occupancies.high_occ,
                                                                      self.options.occupancies.steps+1, endpoint=True)]
        self.OccChoice = wx.Choice(self, wx.ID_ANY, choices=self.occ_list)
        self.OccChoice.SetSelection(0)
        self.occ = self.occ_list[0]
        occStatic = wx.StaticText(self, wx.ID_ANY, "Occupancies : ")

        if self.options.f_and_maps.fast_and_furious:
            choices = ['qfextr']
        else:
            choices = self.options.f_and_maps.f_extrapolated_and_maps
        for choice in choices:
            self.best_occ[choice] = None
        self.FextrChoice = wx.Choice(self, wx.ID_ANY, choices=choices)
        self.FextrChoice.SetSelection(0)
        self.fextr = self.options.f_and_maps.f_extrapolated_and_maps[0]

        FextrStatic = wx.StaticText(self, wx.ID_ANY, "Fextr type : ")
        self.best_occ_Static = wx.StaticText(self, wx.ID_ANY, "best estimation @ ...")
        self.coot_button = wx.Button(self, wx.ID_ANY, label="Open in Coot")
        bmp = wx.Image(os.path.join(script_dir, 'pngs/coot_logo_small.png'), type=wx.BITMAP_TYPE_ANY,
                       index=-1).ConvertToBitmap()
        self.coot_button.SetBitmapMargins((2, 2))
        self.coot_button.SetBitmap(bmp)
        self.coot_button.Bind(wx.EVT_BUTTON, self.onCoot)
        
        self.occNfextrSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.occNfextrSizer.Add(occStatic, 0, wx.ALIGN_CENTER_VERTICAL, border=5)
        self.occNfextrSizer.AddSpacer(10)
        self.occNfextrSizer.Add(self.OccChoice, 0, wx.ALIGN_CENTER_VERTICAL, border=5)
        self.occNfextrSizer.AddSpacer(15)
        self.occNfextrSizer.Add(self.best_occ_Static, 0, wx.ALIGN_CENTER_VERTICAL, border=5)
        self.occNfextrSizer.Hide(self.best_occ_Static)
        self.occNfextrSizer.AddSpacer(30)
        self.occNfextrSizer.Add(FextrStatic, 0, wx.ALIGN_CENTER_VERTICAL, border=5)
        self.occNfextrSizer.AddSpacer(10)                
        self.occNfextrSizer.Add(self.FextrChoice, 0, wx.ALIGN_CENTER_VERTICAL, border=5)
        self.occNfextrSizer.AddSpacer(30)
        self.occNfextrSizer.Add(self.coot_button, 0, wx.ALIGN_CENTER_VERTICAL, border=5)
        self.occNfextrSizer.Layout()
        self.occNfextrSizer.Hide(self.coot_button)
        self.occNfextrSizer.Layout()
        self.ImgSizer = wx.BoxSizer(wx.VERTICAL)

        self.mainSizer.AddSpacer(30)
        self.mainSizer.Add(self.occNfextrSizer, 0,  wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)
        self.mainSizer.AddSpacer(30)
        self.mainSizer.Add(self.ImgSizer, 1, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL)
        self.SetSizer(self.mainSizer)
        self.SetAutoLayout(1)

    def setupBindings(self):
        self.Bind(wx.EVT_CHOICE, self.onSelection, self.OccChoice)
        self.Bind(wx.EVT_CHOICE, self.onSelection, self.FextrChoice)

    def onSelection(self, evt):
        self.occ = self.OccChoice.GetStringSelection()
        self.fextr = self.FextrChoice.GetStringSelection()

        while not self.ImgSizer.IsEmpty():
            item = self.ImgSizer.GetItem(0)
            if item.IsSpacer:
                self.ImgSizer.Detach(0)
            else:
                window = self.ImgSizer.GetItem(0).GetWindow()
                window.Destroy()
        self.ImgSizer.Layout()
        self.mainSizer.Layout()
        path_occ = ''
        if self.fextr.startswith('q'):
            path_occ += 'qweight_'
            Fextr = self.fextr[0] + self.fextr[1].upper() + self.fextr[2:]
        elif self.fextr.startswith('k'):
            path_occ += 'kweight_'
            Fextr = self.fextr[0] + self.fextr[1].upper() + self.fextr[2:]
        else:
            Fextr =self.fextr[0].upper() + self.fextr[1:]

        self.Fextr_png_name = Fextr
        path_occ += 'occupancy_%.3f' % float(self.occ)
        path = os.path.join(self.options.output.outdir, path_occ)


        fn_neg = os.path.join(path, '%s_negative_reflections.pickle' % Fextr)
        if os.path.isfile(fn_neg):
            self.addPlot(fn_neg)
        else:
            self.addEmptyPlot()

        FsigF = os.path.join(path, '%s_FsigF.pickle' % Fextr)
        if os.path.isfile(FsigF):
            self.addPlot(FsigF)
        else:
            self.addEmptyPlot()

        if self.finished:
            self.best_occ_Static.SetLabel("best estimation @ %s"%self.best_occ[self.fextr])#    self.OccChoice.FindString(s)
            if float(self.occ) == float(self.best_occ[self.fextr]):
                try:
                    fn_ddm = re.sub(".png$", ".pickle", self.ddm[self.fextr]) #get name of pickle file for ddm. Not elegant
                    #fn_ddm = get_name(self.ddm[self.fextr])+".pickle"
                    if os.path.isfile(fn_ddm):
                        self.addPlot(fn_ddm)
                    else:
                        self.addImg(self.ddm[self.fextr])
                except AttributeError:
                    pass
                except TypeError:
                    pass
                    
                if not self.coot_button.IsShown():
                    self.occNfextrSizer.Show(self.coot_button)
            else:
                if self.coot_button.IsShown():
                    self.occNfextrSizer.Hide(self.coot_button)
        self.mainSizer.Layout()
        evt.Skip()

    def onCoot(self, evt):
        os.system("coot --script %s &" % self.coot_scripts[self.fextr])

    def plot_FsigF(self, pickle_file):
        """
        For FsigF plot.
        This plot should be placed in the (q/k-weighted)_occupancy directory
        This plot should be prepared for each type of extrapolated structure factor
        prefix can be qFextr, kFextrm, Fextr, qFgenick, kFgenick, Fgenick, qFextr_calc, kFextr_calc, Fextr_calc
        """
        #prefix = pickle_file.strip('_FsigF.pickle')
        # pickle_file = '%s_FsigF.pickle' % (prefix)
        prefix = get_name(pickle_file)
        if os.path.isfile(pickle_file) == False:
            return

        with open(pickle_file, 'rb') as stats_file:
            bin_res_cent_lst, f_sigf_lst, s, l, ids, idl = pickle.load(stats_file)

        #fig, ax1 = plt.subplots(figsize=(10, 5))
        self.figure = Figure(figsize=(10, 5))
        ax1 = self.figure.add_subplot(111)

        f_sigf_lst = np.asarray(f_sigf_lst)
        bin_res_cent_lst = np.asarray(bin_res_cent_lst)
        s = np.full((len(f_sigf_lst), 1), 0.8)
        l = np.full((len(f_sigf_lst), 1), 1.2)
        ax1.plot(bin_res_cent_lst[:], f_sigf_lst[:], marker = '.', label='<F/sig(F)>', color='red')
        ax1.plot(bin_res_cent_lst[:], s[:], linestyle=':', label='<F/sig(F)> = 0.8', color='blue')  # (<I/sig(I)> = 2)
        ax1.plot(bin_res_cent_lst[:], l[:], linestyle=':', label='<F/sig(F)> = 1.2',
                 color='green')  # (<I/sig(I)> = 1.5)

        ax1.plot(np.array([ids]), np.array([0.8]), marker = 's', markersize=3, color='blue', label='estimation: %.2f A' % (ids))
        ax1.plot(np.array([idl]), np.array([1.2]), marker = '^', markersize=5, color = 'green', label='estimation: %.2f A' % (idl))
        ax1.set_xlabel('Resolution of bin center (A)')
        ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
        ax1.set_ylabel('<F/sig(F)>')
        ax1.legend(loc='lower right', bbox_to_anchor=(0.79, -0.05, 0.45, 0.5), fontsize='xx-small', framealpha=0.5)
        ax1.set_title('%s: <F/sig(F)> for high resolution bins' % (prefix), fontsize='medium', fontweight="bold")
        self.figure.subplots_adjust(hspace=0.35, left=0.09, right=0.82, top=0.95)
        canvas = FigureCanvas(self, -1, self.figure)
        return canvas

    def plot_negativereflections(self, pickle_file):
        """
        For negative_reflections plot.
        This plot should be placed in the (q/k-weighted)_occupancy directory. If this script is not launched from the (q/k-weighted)_occupancy directory, an additional argument concerning the directory should be added, or pickle file should be directy provided.
        This plot should be prepared for each type of extrapolated structure factor
        prefix can be qFextr, kFextrm, Fextr, qFgenick, kFgenick, Fgenick, qFextr_calc, kFextr_calc, Fextr_calc
        """

        if os.path.isfile(pickle_file) == False:
            return

        with open(pickle_file, 'rb') as stats_file:
            bin_res_cent_lst, neg_lst, neg_percent_lst, comp_lst, comp_true_lst, s = pickle.load(stats_file)

        # In Fextr_utils, "maptype" is used. In order to keep the plotting code the same, let's use the same name

        self.figure = Figure(figsize=(10, 5))
        ax1, ax3 = self.figure.subplots(1, 2)

        ax1.plot(bin_res_cent_lst[:], neg_lst[:], marker = ".", label='# Neg. ESFAs', color='red')
        ax1.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
        ax1.set_ylim(0, np.max(neg_lst))
        ax1.set_xlabel('Resolution (A)')
        ax1.set_ylabel('Absolute number of negative ESFAs')
        ax1.yaxis.label.set_color('red')
        ax1.set_title("Negative ESFAs for high resolution bins", fontsize='medium', fontweight="bold")
        ax2 = ax1.twinx()
        ax2.plot(bin_res_cent_lst[:], neg_percent_lst[:],  marker = 's', markersize = 3, label='% Neg. ESFAs', color='blue')
        ax2.set_ylim(0, 100)
        ax2.set_ylabel('Negative ESFAs in resolution bin (%)')
        ax2.yaxis.label.set_color('blue')
        lines_labels_1 = [ax.get_legend_handles_labels() for ax in [self.figure.axes[0], self.figure.axes[2]]]
        lines_1, labels_1 = [sum(lne, []) for lne in zip(*lines_labels_1)]
        ax1.legend(lines_1, labels_1, loc='lower right', bbox_to_anchor=(0.82, -0.05, 0.45, 0.5), fontsize='xx-small',
                   framealpha=0.5)

        s = np.full((bin_res_cent_lst.shape[0], 1), 90)
        ax3.plot(bin_res_cent_lst[:], comp_lst[:], marker = '.', label='Completeness', color='red')
        ax3.plot(bin_res_cent_lst[:], comp_true_lst[:], marker = 's', markersize=3, label='True completeness', color='blue')
        ax3.plot(bin_res_cent_lst[:], s[:], linestyle=':', label='90 (%) Completeness', color='green')
        ax3.set_xlim(np.max(bin_res_cent_lst[1:]), np.min(bin_res_cent_lst[1:]))
        ax3.set_ylim(0, 100)
        ax3.set_xlabel('Resolution (A)')
        ax3.set_ylabel('Completeness in resolution bin (%)')
        ax3.set_title("Completeness for high resolution bins", fontsize='medium', fontweight="bold")

        lines_labels_2 = self.figure.axes[1].get_legend_handles_labels()
        lines_2, labels_2 = [sum(lne, []) for lne in zip(lines_labels_2)]
        # ax3.legend(lines_2, labels_2, fontsize = 'x-small', framealpha=0.5, loc=3, bbox_to_anchor=(0.05, 0.05, 0.5, 0.5))
        ax3.legend(lines_2, labels_2, loc='lower right', bbox_to_anchor=(0.82, -0.05, 0.45, 0.5), fontsize='xx-small',
                   framealpha=0.5)
        self.figure.subplots_adjust(hspace=0.25, wspace=0.5, left=0.09, right=0.88, top=0.95)
        canvas = FigureCanvas(self, -1, self.figure)
        return canvas
    
    def make_colormap(self,seq):
        """Return a LinearSegmentedColormap
        seq: a sequence of floats and RGB-tuples. The floats should be increasing
        and in the interval (0,1).
        """
        seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
        cdict = {'red': [], 'green': [], 'blue': []}
        for i, item in enumerate(seq):
            if isinstance(item, float):
                r1, g1, b1 = seq[i - 1]
                r2, g2, b2 = seq[i + 1]
                cdict['red'].append([item, r1, r2])
                cdict['green'].append([item, g1, g2])
                cdict['blue'].append([item, b1, b2])
        return mcolors.LinearSegmentedColormap('CustomMap', cdict)


    def plot_ddm(self, pickle_file, scale = None):
        """
        For plotting the ddm.
        This plot should be prepared for each of the Fextr types but only with the model from the best estimated occupancy and should calculated only at the end
        The input pickle contains information on the chain id, ddm_residue and seq_info_unique
        The output filename should be the same as of the pickle file
        """
        
        if os.path.isfile(pickle_file) == False:
            return

        with open(pickle_file,'rb') as stats_file:
            while True:
                try:
                    stats = np.array(pickle.load(stats_file))
                    try:
                        alldata = np.vstack([alldata,stats[np.newaxis,...]])
                    except NameError:
                        alldata =stats[np.newaxis,...]
                except EOFError:
                    break
        
        chains = np.unique(alldata[:,0])
        num_chains = chains.shape[0]
            
        if num_chains == 1:
            n_cols = 1
        else:
            n_cols = 2
            
        n_rows_prot = int(math.ceil(num_chains/n_cols))

        n_rows = n_rows_prot

        outname = get_name(pickle_file)
        outname_pdf = '%s.pdf'%(outname)
        outname_png = '%s.png'%(outname)
        
        if scale == None:
            scale = 0
            for ddm_residue in alldata[:,1]:
                mx_temp = np.max(ddm_residue)
                mn_temp = np.min(ddm_residue)
                if mx_temp > scale:
                    scale = mx_temp
                if np.abs(mn_temp) > scale:
                    scale = np.abs(mn_temp)
        else:
            scale = scale
            
        if n_rows == n_cols == 1:
            self.figure = Figure(figsize=(10, 10), tight_layout=True)
        else:
            self.figure = Figure(figsize=(10*n_rows, 5*n_cols), tight_layout=True)
        #self.figure = Figure(figsize=(10,5), tight_layout=True)
        axs = self.figure.subplots(n_rows, n_cols, squeeze=False)
        
        #fig, axs = plt.subplots(n_rows, n_cols, figsize=(10*n_rows, 10*n_cols), squeeze=False)#, constrained_layout=True)
        c = mcolors.ColorConverter().to_rgb
        rvb = self.make_colormap([c('blue'), c('white'), 0.40, c('white'),  0.60, c('white'), c('red')])

        col = -1
        row = -1
        old_chain_id = ''
        for chain in alldata:
            ID, ddm_residue,seq_info_unique = chain
            if ID != old_chain_id:
                if col == 0:
                    col = 1
                else:
                    col = 0
                    row +=1
                            
            mask = np.zeros_like(ddm_residue, dtype=np.bool)
            mask[np.triu_indices_from(mask)] = True
            FINAL2 = np.ma.array(ddm_residue, mask=mask)
            
            tick_jump = int(np.round(len(seq_info_unique)/15,0)) #we want to add 15 seq_ticks
            seq_ticks = seq_info_unique[0::tick_jump] #residues for which we will show tick positions
            tick_pos = list(map(lambda x: x*tick_jump, range(len(seq_ticks))))
            
                
            img = axs[row, col].imshow(FINAL2, cmap=rvb, vmin=-scale, vmax=scale)
            
            axs[row, col].set_xticks(tick_pos)
            axs[row, col].set_xticklabels(seq_ticks)
            
            for tick in axs[row, col].get_xticklabels():
                tick.set_rotation(90)
            #self.figure.setp(axs[row, col].xaxis.get_majorticklabels(), rotation=90, fontsize='small')
            
            axs[row, col].set_yticks(tick_pos)
            axs[row, col].set_yticklabels(seq_ticks)

            axs[row, col].set_title('Chain %s' %ID, fontsize = 'medium',fontweight="bold")
            axs[row, col].set_xlabel('Residues')
            axs[row, col].set_ylabel('Residues')
            
            axs[row, col].spines['top'].set_visible(False)
            axs[row, col].spines['right'].set_visible(False)
            
            self.figure.colorbar(img, ax=axs[row, col], fraction=0.046, pad=0.04)
                
            old_chain_id = ID

            
        #fig.tight_layout()
        #plt.savefig(outname_pdf, dpi=300, transparent=True)
        #plt.savefig(outname_png, dpi=300)
        #plt.close()
        
        canvas = FigureCanvas(self, -1, self.figure)
        return canvas


    def addEmptyPlot(self):
        self.figure = Figure(figsize=(10, 5))
        ax1 = self.figure.subplots(1, 1)
        ax1.text(0.3, 0.5, "Data not available yet")
        #ax1.set_title("Negative reflections for high resolution bins", fontsize='medium', fontweight="bold")
        self.figure.subplots_adjust(hspace=0.25, wspace=0.5, left=0.09, right=0.88, top=0.95)
        canvas = FigureCanvas(self, -1, self.figure)
        self.ImgSizer.Add(canvas, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_CENTER_HORIZONTAL)

        toolbar = NavigationToolbar(canvas)
        toolbar.Realize()
        # By adding toolbar in sizer, we are able to put it at the bottom
        # of the frame - so appearance is closer to GTK version.
        self.ImgSizer.Add(toolbar, 0, wx.LEFT | wx.EXPAND)

        # update the axes menu on the toolbar
        toolbar.update()
        self.ImgSizer.AddSpacer(60)
        # self.SetSizer(self.sizer)]
        self.ImgSizer.Layout()
        self.mainSizer.Layout()
        self.FitInside()


    def addPlot(self, pickle_file):

        _, pickle_name = os.path.split(pickle_file)
        if 'negative_reflections'in pickle_name:
            canvas = self.plot_negativereflections(pickle_file)
        elif pickle_name.startswith("ddm_"):
            canvas = self.plot_ddm(pickle_file)
        else:
            canvas = self.plot_FsigF(pickle_file)
        self.ImgSizer.Add(canvas, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_CENTER_HORIZONTAL)

        toolbar = NavigationToolbar(canvas)
        toolbar.Realize()
        # By adding toolbar in sizer, we are able to put it at the bottom
        # of the frame - so appearance is closer to GTK version.
        self.ImgSizer.Add(toolbar, 0, wx.LEFT | wx.EXPAND)

        # update the axes menu on the toolbar
        toolbar.update()
        self.ImgSizer.AddSpacer(60)
        # self.SetSizer(self.sizer)
        self.ImgSizer.Layout()
        self.mainSizer.Layout()
        self.FitInside()

        return

    def addImg(self, filepath):
            img = wx.Image(filepath, wx.BITMAP_TYPE_ANY)
            W = img.GetWidth()
            H = img.GetHeight()
            if W > H:
                NewW = self.photoMaxSize
                NewH = self.photoMaxSize * H / W
            else:
                NewH = self.photoMaxSize
                NewW = self.photoMaxSize * W / H
            img = img.Scale(NewW, NewH)

            self.newimg = wx.StaticBitmap(self, wx.ID_ANY,
                                          wx.BitmapFromImage(img))
            self.ImgSizer.Add(self.newimg, proportion=0, flag=wx.ALIGN_CENTER_HORIZONTAL)
            self.ImgSizer.Layout()
            self.mainSizer.Layout()
            self.FitInside()

    def onFinished(self):
        pickle_fn = os.path.join(self.options.output.outdir,"occupancy_recap.pickle")
        results = pickle.load(open(pickle_fn, "rb"))
        fextrs = self.options.f_and_maps.f_extrapolated_and_maps
        for fextr in fextrs:
            Fextr=''
            for s in fextr:
                if s == 'f':
                    Fextr+= s.upper()
                else:
                    Fextr += s
            occ, coot, ddm = results[Fextr]
            self.best_occ[fextr] = occ
            self.coot_scripts[fextr] = coot
            self.ddm[fextr] = ddm
        fextr = self.FextrChoice.GetStringSelection()
        self.occNfextrSizer.Show(self.best_occ_Static)
        self.best_occ_Static.SetLabel("best estimation @ %s"%self.best_occ[fextr])
        self.finished = True











