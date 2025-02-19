# -*- coding: utf-8 -*-
"""
-------
authors and contact information
-------
Elke De Zitter - elke.de-zitter@ibs.fr
Nicolas Coquelle - nicolas.coquelle@esrf.fr
Thomas Barends - Thomas.Barends@mpimf-heidelberg.mpg.de
Jacques Philippe Colletier - jacques-Philippe.colletier@ibs.fr
-------
license information
-------
Copyright (c) 2021 Elke De Zitter, Nicolas Coquelle, Thomas Barends and Jacques-Philippe Colletier
see https://github.com/ElkeDeZitter/Xtrapol8/blob/main/LICENSE
-------
"""
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import (
    FigureCanvasWxAgg as FigureCanvas,
    NavigationToolbar2WxAgg as NavigationToolbar,
)
from matplotlib.figure import Figure

import os.path
import glob
import threading
from threading import Thread
import sys
import wx
import subprocess
from wx.lib.pubsub import pub
from gui import panelIO, panelExtrapolation, panelRefinement, panelLog
from gui.panelLog import TabLog, TabMainImg, TabOccResults

from libtbx.phil import parse

from Fextr import master_phil
from wx.aui import AuiNotebook
import version
#from wx.lib.agw.flatnotebook import FlatNotebook as AuiNotebook


script_dir = os.path.dirname(os.path.abspath(__file__))
# Notebook styles to have, or not, a widget to close the tab
bookStyleNO = wx.aui.AUI_NB_DEFAULT_STYLE & ~(wx.aui.AUI_NB_CLOSE_ON_ACTIVE_TAB)
bookStyleYES = wx.aui.AUI_NB_DEFAULT_STYLE & wx.aui.AUI_NB_CLOSE_ON_ACTIVE_TAB

# Get access to file types on Mac
wx.SystemOptions.SetOption(u"osx.openfiledialog.always-show-types","1")

########################################################################
class NotebookConfigure(AuiNotebook):
    """
    Notebook for all X8 parameters
    tabIO: Inputs and Outputs
    tabExt X8 settings
    Third tab: Refinement
    """

    # ----------------------------------------------------------------------
    def __init__(self, parent):
        AuiNotebook.__init__(self, parent, id=wx.ID_ANY, style=bookStyleNO)
        #wx.BK_DEFAULT
                             # wx.BK_TOP
                             # wx.BK_BOTTOM
                             # wx.BK_LEFT
                             # wx.BK_RIGHT

        # Create the first tab and add it to the notebook
        self.tabIO = panelIO.TabIO(self)
        self.AddPage(self.tabIO, "Input / Output")

        # Create and add the second tab
        self.tabExt = panelExtrapolation.TabExtrapolation(self)
        self.AddPage(self.tabExt, "FoFo / Extrapolation")

        # Create and add the third tab
        self.tabRefine = panelRefinement.TabRefinement(self)
        self.AddPage(self.tabRefine, "Refinement")


class NoteBookResults(AuiNotebook):
    # ----------------------------------------------------------------------
    def __init__(self, parent, options):
        AuiNotebook.__init__(self, parent, id=wx.ID_ANY, style=bookStyleNO)
        self.tabLog = TabLog(self, options)
        self.options = options
        self.AddPage(self.tabLog, "Log")
        font = wx.Font(10, wx.TELETYPE, wx.NORMAL, wx.NORMAL, False, 'Consolas')
        self.tabLog.LogTextCtrl.SetFont(font)

        self.tabImg = TabMainImg(self)
        self.AddPage(self.tabImg, "Main")

        self.tabOcc = TabOccResults(self, self.options)
        if not self.options.output.generate_fofo_only:
            self.AddPage(self.tabOcc,"Occupancies")
        self.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, self.OnPageChanged, self)

    def OnPageChanged(self, evt):
        index = self.GetSelection()
        if index == 2:
            self.tabOcc.onSelection(evt)

    def onFinished(self):
        # Update statusbar
        if self.options.output.generate_fofo_only:
            self.tabLog.CreateCoot()
        else:
            self.tabOcc.onFinished()


class MainNotebook(AuiNotebook):

    def __init__(self, parent):
        AuiNotebook.__init__(self, parent, id=wx.ID_ANY, style=bookStyleNO)
        self.Configure = NotebookConfigure(self)
        self.AddPage(self.Configure, "Configure")
        self.Runs = -1
        self.ResultsBooks = []
        self.threads = []
        self.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, self.OnPageChanged, self)
        pub.subscribe(self.OnUpdateLog, 'update_log')
        pub.subscribe(self.onFinished, 'END')

    def StopRun(self):
        pageIdx = self.GetSelection()
        if pageIdx > 0:
            thread = self.threads[pageIdx - 1]
            if thread is not None:
                if not thread.stopped():
                    thread.stop()
                    print("Run %i stopped" % pageIdx)
                else:
                    print("Run %i already stopped" % pageIdx)


    def OnUpdateLog(self, Nlog, line):
        self.ResultsBooks[Nlog].tabLog.updateLog(line)

    def onFinished(self, Nlog):
        self.ResultsBooks[Nlog].onFinished()

    def OnPageChanged(self, evt):
        pg = self.GetCurrentPage()
        self.idx = self.GetPageIndex(pg)

        if self.idx > 0 : self.SetWindowStyleFlag(bookStyleYES)
        else : self.SetWindowStyleFlag(bookStyleNO)

    def OnPageClose(self, idx):
        SelectedThread = self.threads[idx]
        if SelectedThread.is_alive():
            Stop = wx.MessageDialog(None, 'Job is not finished!\n Do you want to stop it ?', 'WorkStatus', wx.YES_NO | wx.NO_DEFAULT).ShowModal()
            if Stop == wx.ID_YES:
                SelectedThread.stop()
            else:
                evt.Veto()
                return
        self.threads.pop(idx)


class X8Thread(Thread):
    """This is the thread which will run the code"""
    # ----------------------------------------------------------------------
    def __init__(self, input_phil, Nlog):
        self.input = input_phil
        self.Nlog = Nlog
        Thread.__init__(self)
        self.daemon = True
        self._stop = threading.Event()

    #----------------------------------------------------------------------
    def run(self):
        count = 0
        # For debugging purpose
        #pub.sendMessage("END", Nlog=self.Nlog)
        #return
        p = subprocess.Popen(['phenix.python', os.path.join(script_dir, 'Fextr.py'), 'tmp.phil'],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        while True:

            line = p.stdout.readline()
            if not line:
                pub.sendMessage("END", Nlog=self.Nlog)
                return
            if self.stopped():
                print("Run stopped")
                p.terminate()
                return

            wx.CallAfter(pub.sendMessage, "update_log", line=line, Nlog=self.Nlog)


    def stop(self):
        self._stop.set()

    def stopped(self):
        return self._stop.isSet()

########################################################################
########################################################################
class MainFrame(wx.Frame):
    """
    Main Frame holding all widgets.
    This Frame also holds part of the model (design to be improved)
    """
    # ----------------------------------------------------------------------
    def __init__(self,args):
        """Constructor"""
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "XtrapolG8 -- version %s" %(version.VERSION),
                          size=(1100, 1000)
                          )


        default_font = wx.Font(11, wx.MODERN, wx.NORMAL, wx.NORMAL, False, 'MS Shell Dlg 2')

        # Adding a MenuBar
        menubar = wx.MenuBar()
        fileMenu = wx.Menu()
        fileItemOpenPhil = fileMenu.Append(wx.ID_ANY, 'Open phil', 'Open input phil')
        fileItemLoadResults = fileMenu.Append(wx.ID_ANY, 'Load Results', 'Load Results from previous Xtrapol8 runs')
        fileItemSavePhil = fileMenu.Append(wx.ID_ANY, 'Save phil', 'Save all inputs in a phil file')
        fileItemQuit = fileMenu.Append(wx.ID_EXIT, 'Quit', 'Quit application')
        menubar.Append(fileMenu, "&File")
        self.SetFont(default_font)
        self.SetMenuBar(menubar)
        self.Bind(wx.EVT_MENU, self.OnClose, fileItemQuit)
        self.Bind(wx.EVT_MENU, self.OnOpenPhil, fileItemOpenPhil)
        self.Bind(wx.EVT_MENU, self.OnLoadResults, fileItemLoadResults)
        self.Bind(wx.EVT_MENU, self.OnSavePhil, fileItemSavePhil)
        self.Bind(wx.EVT_CLOSE, self.OnClose)

        # Adding the ToolBar
        self.ToolBar = wx.ToolBar(self, -1)
        self.ToolBar.SetToolBitmapSize(size=(1, 1))
        # self.ToolBar.AddTool(101, wx.Bitmap(os.path.join(script_dir,"gui/pngs/settings_scaled.png")))
        self.ToolBar.AddTool(102, wx.Bitmap(os.path.join(script_dir,"gui/pngs/run_scaled.png")))
        self.ToolBar.AddTool(103, wx.Bitmap(os.path.join(script_dir,"gui/pngs/cancel_scaled.png")))
        self.ToolBar.Bind(wx.EVT_TOOL, self.OnToolBar)
        self.ToolBar.Realize()

        #The panel will hold the MainNotebook
        panel = wx.Panel(parent=self)
        self.notebook = MainNotebook(panel)
        self.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CLOSE, self.OnPageClose, self.notebook) #Closgin a run tab

        # This list will holds all inputs (type phil objects) of the different runs (only Run tabs in the Gui which are still visible)
        self.inputs = []

        # Timers will be used to update the gui with the different figures
        self.timer = wx.Timer(self)
        self.timerFextr = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.update, self.timer)
        self.Bind(wx.EVT_TIMER, self.updateFextr, self.timerFextr)
        pub.subscribe(self.updateFextr, "updateFextr")
        #pub.subscribe(self.X8ModeChanged, "X8Mode")

        # These variables will be used by the timer callbacks to update the gui accordingly
        self.pngs = ["Riso_CCiso.pickle",
                     "q_estimation.pickle",
                     "k_estimation.pickle",
                     'summed_difference_peaks.png',
                     'Fextr_binstats.pickle'] # First set of pngs displayed in the gui (used by the first timer)
        self.png = self.pngs[0]
        self.pngs_idx = []  # This list will hold a reference for each run tab in the gui

        self.Fextr_pngs = ['Fextr_binstats.pickle',
                           'Fextr_negative.pickle',
                           'alpha_occupancy_determination_tmp.pickle',
                           'tmp_refinement_R-factors_per_alpha.pickle',
                           'Distance_difference_plot_tmp.png'] #Second set of pngs (used by the second timer)
                           # - note the tmp par of the string which will be replaced by the appropriate Fextr type

        # Sizer holding the notebbok
        self.fsizer = wx.BoxSizer(wx.VERTICAL)
        self.fsizer.Add(self.notebook, 1, wx.ALL | wx.EXPAND, 10)
        panel.SetSizer(self.fsizer)
        self.Layout()
        self.Show()
        if len(args) > 1:
            phil = args[1]
            #try:
            self.OnOpenPhil(event=None,phil_file=phil)
            #except:
            #    pass


    #def X8ModeChanged(self ,mode):
        #old, new = mode
        #modes = ["FoFo", "FNF", "CNC"]
        #new_mode = modes[new]
        ## Saving phil for future restauration
        #input_phil = self.extract_phil()
        #modified_phil = master_phil.format(python_object=input_phil)
        #modified_phil.show(out=open(".%s.phil"%old, "w"))

        ## Restauration if possible
        #self.notebook.Configure.tabExt.currentX8Mode = new_mode
        #phil_file = ".%s.phil" % new_mode
        #if os.path.exists(phil_file):
            #user_params = self.extract_debug_phil(open(phil_file).read())
            #self.SetWidgetsTabExt(user_params,SetX8=False)
            ##self

    def update(self, event):
        """
        This function is called by the timer
        It looks for results in the output directory to display in the gui.
        (Only if the thread is active, and only for the run the user is currently looking at)
        :param event: wx.EVT_TIMER
        :return: None
        """
        run = self.notebook.GetSelection() - 1
        if run >= 0:
            if self.notebook.threads[run] is not None and self.notebook.threads[run].is_alive():
                if self.pngs_idx[run] < len(self.pngs):
                    png = self.pngs[self.pngs_idx[run]]
                    
                    filepath = os.path.join(self.inputs[run].output.outdir, png)

                    if os.path.isfile(filepath):
                        if filepath.endswith('pickle'):
                            self.notebook.ResultsBooks[run].tabImg.addPlot(filepath)
                        else:
                            self.notebook.ResultsBooks[run].tabImg.addImg(filepath)

                        if self.pngs_idx[run] <= len(self.pngs) - 1:
                            self.pngs_idx[run] += 1
                else:
                    if not self.timerFextr.IsRunning():
                        self.timerFextr.Start(2000)
                    if not hasattr(self.notebook.ResultsBooks[run].tabImg, 'FextrSelection'):
                        self.notebook.ResultsBooks[run].tabImg.addChoices(self.inputs[run].f_and_maps.f_extrapolated_and_maps)
                        self.notebook.ResultsBooks[run].tabImg.FextrSelection.Bind(wx.EVT_CHOICE, self.notebook.ResultsBooks[run].tabImg.Clear)

    def updateFextr(self, evt):
        run = self.notebook.GetSelection() - 1
        if run >= 0:
            if hasattr(self.notebook.ResultsBooks[run].tabImg, 'FextrSelection'):
                tab = self.notebook.ResultsBooks[run].tabImg
                Fextr = tab.FextrSelection.GetStringSelection()
                Total = tab.ImgSizer.GetItemCount()
                for j in range(Total, len(self.Fextr_pngs)):
                    if Fextr[0] in ['q','k']:
                        Fextr = Fextr[0] + Fextr[1].upper() + Fextr[2:]
                    else:
                        Fextr = Fextr[0].upper() + Fextr[1:]
                    png = self.Fextr_pngs[j].replace('tmp', Fextr)
                    filepath = os.path.join(self.inputs[run].output.outdir, png)

                    if os.path.isfile(filepath):
                        if filepath.endswith('pickle'):
                            tab.addFextrPlot(Fextr, filepath)
                        else:
                            tab.addFextrImg(filepath)
                    else:
                        #return
                        print("%s does not exists" %filepath)

    def OnPageClose(self, evt):
        # will check that the run is not running - will clean its thread list accordingly
        # Variables to update after deleting the run tab of the main notebook

        run = self.notebook.GetSelection() - 1
        SelectedThread = self.notebook.threads[run]
        if SelectedThread is not None:
            if SelectedThread.is_alive():
                Stop = wx.MessageDialog(None, 'Job is not finished!\n Do you want to stop it ?', 'WorkStatus',
                                        wx.YES_NO | wx.NO_DEFAULT).ShowModal()
                # print Stop
                if Stop == wx.ID_YES:
                    print("Clicked YES")
                    SelectedThread.stop()

                else:
                    evt.Veto()
                    return
        self.notebook.threads.pop(run)
        self.notebook.ResultsBooks.pop(run)
        self.inputs.pop(run)
        self.pngs_idx.pop(run)
        self.notebook.Runs -= 1


    def OnClose(self, event):
        # Should check if any job is running
        if self.timer.IsRunning(): self.timer.Stop()
        if self.timerFextr.IsRunning(): self.timer.Stop()
        phils = glob.glob(".*phil")
        for phil in phils:
            os.remove(phil)
        self.Destroy()

    def OnToolBar(self, event):
        id = event.GetId()
        # if id == 101:
        #     return
        if id == 102:
            self.OnrunX8()
        elif id == 103:
            self.OnStopRun()

    def AddResultsTab(self):
        self.notebook.Runs += 1
        self.notebook.ResultsBooks.append(NoteBookResults(self.notebook, self.input_phil))
        self.notebook.AddPage(self.notebook.ResultsBooks[self.notebook.Runs], "Run #%i" % (self.notebook.Runs + 1))
        

    def OnrunX8(self,):
        if self.check_user_input():

            self.input_phil = self.extract_phil()
            modified_phil = master_phil.format(python_object=self.input_phil)
            modified_phil.show(out=open("tmp.phil", "w"))
            self.AddResultsTab()
            if not self.timer.IsRunning():
                self.timer.Start(2000)
            thread = X8Thread(self.input_phil, self.notebook.Runs)
    
            ## All this should be in a class that deals with the outputs
            self.notebook.threads.append(thread)
            self.pngs_idx.append(0)
            self.inputs.append(self.input_phil)
            thread.start()

    def OnLoadResults(self, evt):
        PathResults = self.onBrowseDir(evt=None)
        Phil = os.path.join(PathResults, 'Xtrapol8_out.phil')
        if os.path.exists(Phil):
            print(Phil)
            self.OnOpenPhil(event=None, phil_file=Phil)
            self.input_phil = self.extract_phil()
            self.AddResultsTab()
            log = glob.glob(os.path.join(PathResults, "*Xtrapol8*.log"))[0]
            self.pngs_idx.append(0)
            self.inputs.append(self.input_phil)
            self.notebook.threads.append(None)
            run = self.notebook.Runs
            self.notebook.ResultsBooks[run].tabLog.LogTextCtrl.WriteText(open(log).read())
            for png in self.pngs:
                filepath = os.path.join(PathResults, png)
                if os.path.isfile(filepath):
                    #print(png)
                    if filepath.endswith("pickle"):
                        self.notebook.ResultsBooks[run].tabImg.addPlot(filepath)
                    else:
                        self.notebook.ResultsBooks[run].tabImg.addImg(filepath)
            self.notebook.ResultsBooks[run].tabImg.addChoices(self.inputs[run].f_and_maps.f_extrapolated_and_maps)
            self.notebook.ResultsBooks[run].tabImg.FextrSelection.Bind(wx.EVT_CHOICE,
                                                                       self.notebook.ResultsBooks[run].tabImg.Clear)

            #if hasattr(self.notebook.ResultsBooks[run].tabImg, 'FextrSelection'):
            tab = self.notebook.ResultsBooks[run].tabImg
            Fextr = tab.FextrSelection.GetStringSelection()
            Total = tab.ImgSizer.GetItemCount()
            for j in range(Total, len(self.Fextr_pngs)):
                if Fextr[0] in ['q','k']:
                    Fextr = Fextr[0] + Fextr[1].upper() + Fextr[2:]
                else:
                    Fextr = Fextr[0].upper() + Fextr[1:]
                png = self.Fextr_pngs[j].replace('tmp', Fextr)
                filepath = os.path.join(PathResults, png)
                if os.path.isfile(filepath):
                    if filepath.endswith('pickle'):
                        tab.addFextrPlot(Fextr, filepath)
                    else:
                        tab.addFextrImg(filepath)
                else:
                    print("%s does not exists" %filepath)
            if self.inputs[run].output.generate_fofo_only:
                self.notebook.ResultsBooks[run].tabLog.CreateCoot()
            else:
                self.notebook.ResultsBooks[run].tabOcc.onFinished()

        return
    
    def check_user_input(self):
        tabIO = self.notebook.Configure.tabIO
        message_err="Error with your input files.\nXtrapol8 needs:\n"
        err = 0
        if len(tabIO.files["Reference model"]) == 0:
            message_err += "\n- a reference model (pdb or cif)"
            err = 1
        if len(tabIO.files["Reference model"]) > 1:
            message_err += "\n- a single reference model (pdb or cif)"
            err = 1
        if len(tabIO.files["Reference mtz"]) == 0:
            message_err += "\n- a reference mtz (mtz or cif)"
            err = 1
        if len(tabIO.files["Reference mtz"]) > 1:
            message_err += "\n- a single reference mtz (mtz or cif)"
            err = 1

        if len(tabIO.files["Triggered mtz"]) == 0:
            message_err += "\n- at least one triggered mtz (mtz or cif)"

        if len(tabIO.outdir_sizer.TextCtrl.GetValue()) == 0:
            tabIO.outdir_sizer.TextCtrl.SetValue(os.getcwd()+'/Xtrapol8')
        else:
            path = tabIO.outdir_sizer.TextCtrl.GetValue()
            if os.path.exists(path):
                #Keep outdir given by user if it empty
                if len(os.listdir(path)) == 0:
                    path = path
                ##Keep outdir given by user if it only contains Xtrapol8 log-files:
                #elif len([fle for fle in os.listdir(path) if fle.endswith("Xtrapol8.log")]) == len(os.listdir(path)):
                    #path = path
                else:
                    path = self.get_new_path(path)

            tabIO.outdir_sizer.TextCtrl.SetValue(path)
        if err == 1:
            message_err += '.'
            _ = wx.MessageDialog(self, message=message_err, style=wx.OK).ShowModal()
            return False
        else:
            return True

    def get_new_path(self, path):
        root = path.split('_')
        if len(root) == 1:
            path += '_1'
        else:
            try:
                N = int(root[-1]) + 1
                path = '_'.join(root[:-1]) + '_%i' % N
            except ValueError:
                path += '_1'
        if os.path.exists(path):
            return self.get_new_path(path)
        else:
            return path

    def OnStopRun(self):
        self.notebook.StopRun()

    def extract_debug_phil(self, phil_str):
        user_phil = parse(phil_str)
        user_params = master_phil.fetch(source=user_phil).extract()
        #modified_phil = master_phil.format(python_object=user_params)
        #modified_phil.show(out=open("tmp.phil", "w"))
        return user_params

    def OnOpenPhil(self, event,phil_file=None):
        wildcard = "Phil files (*.phil)|*.phil|" \
                   "All files (*.*)|*.*"
        if phil_file is None:
            phil_file = self.onBrowse(wildcard=wildcard)
        user_params = self.extract_debug_phil(open(phil_file).read())
        self.SetWidgets(user_params)

    def OnSavePhil(self, event):
        filename = self.onBrowse(style=wx.FD_SAVE)
        input_phil = self.extract_phil()
        modified_phil = master_phil.format(python_object=input_phil)
        modified_phil.show(out=open(filename, "w"))

    def SetWidgets(self, user_params):
        self.SetWidgetsTabIO(user_params)
        self.SetWidgetsTabExt(user_params)
        self.SetWidgetsTabRef(user_params)
    
    def SetWidgetsTabIO(self, user_params):
        tabIO = self.notebook.Configure.tabIO

        # Clear List
        tabIO.list.DeleteAllItems()
        tabIO.files = {'Reference model': [],
                      'Reference mtz': [],
                      'Triggered mtz': [],
                      'Restraints': []}

        # Fill listCtrl with input files
        if user_params.input.reference_mtz is not None:
            index = tabIO.list.InsertStringItem(sys.maxint, user_params.input.reference_mtz)
            tabIO.list.SetStringItem(index, 1, "Reference mtz")
            tabIO.files["Reference mtz"].append(user_params.input.reference_mtz)
            #tabIO.extract_dmin_dmax(user_params.input.reference_mtz)
        if user_params.input.triggered_mtz is not None:
            index = tabIO.list.InsertStringItem(sys.maxint, user_params.input.triggered_mtz)
            tabIO.list.SetStringItem(index, 1, "Triggered mtz")
            tabIO.files["Triggered mtz"].append(user_params.input.triggered_mtz)
            #tabIO.extract_dmin_dmax(user_params.input.triggered_mtz)
        if user_params.input.reference_pdb is not None:
            index = tabIO.list.InsertStringItem(sys.maxint, user_params.input.reference_pdb)
            tabIO.list.SetStringItem(index, 1, "Reference model")
            tabIO.files["Reference model"].append(user_params.input.reference_pdb)

        if user_params.input.additional_files is not None:
            for cif in user_params.input.additional_files:
                index = tabIO.list.InsertStringItem(sys.maxint, cif)
                tabIO.list.SetStringItem(index, 1, "Restraints")
                tabIO.files["Restraints"].append(cif)
        # Resolution
        if user_params.input.high_resolution is not None:
            try:
                tabIO.highRes.SetValue("%4.2f" % float(user_params.input.high_resolution))
            except ValueError:
                pass
        if user_params.input.low_resolution is not None:
            try:
                tabIO.lowRes.SetValue("%4.2f" % float(user_params.input.low_resolution))
            except ValueError:
                pass

        #if user_params.input.high_resolution is not "None":
        if user_params.output.outdir is not None:
            tabIO.outdir_sizer.TextCtrl.SetValue(user_params.output.outdir)
        if user_params.output.outname is not None:
            tabIO.outname.SetValue(user_params.output.outname)

    def SetWidgetsTabExt(self, user_params,SetX8=True):
        #####################
        ### Occ - Ext_tab ###
        #####################
        if SetX8:
            self.setX8_mode(user_params)

        tabExt = self.notebook.Configure.tabExt
        tabExt.LowTextCtrl.SetValue(str(user_params.occupancies.low_occ))
        tabExt.HighTextCtrl.SetValue(str(user_params.occupancies.high_occ))
        tabExt.StepsTextCtrl.SetValue(str(user_params.occupancies.steps))

        if user_params.occupancies.list_occ is not None:
            tabExt.ListTextCtrl.SetValue(', '.join([str(x) for x in user_params.occupancies.list_occ]))
        #######################
        ### \Occ - Ext_tab\ ###
        #######################

        ################################
        ### Maps N Scaling - Ext_tab ###
        ###############################
        tabExt.FoChoice.SetStringSelection(user_params.f_and_maps.fofo_type)
        tabExt.ScalingChoice.SetStringSelection(user_params.scaling.b_scaling)
        try:
            tabExt.kscale.SetValue("%4.3f" % float(user_params.f_and_maps.kweight_scale))
        except ValueError:
            pass
        tabExt.updateKScale(None)

        self.fill_map_types(user_params)
        neg_N_missing = user_params.f_and_maps.negative_and_missing
        neg, fill = neg_N_missing.split('_')[0:2]
        if neg in ['fill', 'no']:
            tabExt.negChoice.SetStringSelection('keep')
        else:
            tabExt.negChoice.SetStringSelection(neg)
        if (fill == 'no' or neg == 'no'):
            tabExt.missChoice.SetSelection(1)
        else:
            tabExt.missChoice.SetSelection(0)

        tabExt.peak_detection_thresholdTextCtrl.SetValue(str(user_params.map_explorer.peak_detection_threshold))
        tabExt.peak_integration_floorTextCtrl.SetValue(str(user_params.map_explorer.peak_integration_floor))
        tabExt.RadiusTextCtrl.SetValue(str(user_params.map_explorer.radius))
        tabExt.ZscoreTextCtrl.SetValue(str(user_params.map_explorer.z_score))
        tabExt.OccEstimation.SetStringSelection(user_params.map_explorer.occupancy_estimation)
        if user_params.map_explorer.use_occupancy_from_distance_analysis:
            tabExt.OccEstimation.SetSelection(2)
        if user_params.f_and_maps.fast_and_furious is True:
            if (user_params.map_explorer.occupancy_estimation == 'distance_analysis' or 
            user_params.map_explorer.use_occupancy_from_distance_analysis == True):
                tabExt.OccEstimation.SetSelection(0)
        #tabExt.DistanceAnalysis.SetValue(user_params.map_explorer.use_occupancy_from_distance_analysis)
        
        #Scaling resolution boundaries
        if user_params.scaling.high_resolution is not None:
            try:
                tabExt.ScalingHighRes.SetValue("%4.2f" % float(user_params.scaling.high_resolution))
            except ValueError:
                pass
        if user_params.scaling.low_resolution is not None:
            try:
                tabExt.ScalingLowRes.SetValue("%4.2f" % float(user_params.scaling.low_resolution))
            except ValueError:
                pass
        
    def SetWidgetsTabRef(self, user_params):
        # Refinement        
        tabRef = self.notebook.Configure.tabRefine

        tabRef.RunRef.SetValue(user_params.refinement.run_refinement)
        tabRef.onRefChanged(None)
        
        tabRef.SoftChoiceReci.SetStringSelection(user_params.refinement.reciprocal_space)
        tabRef.SoftChoiceReal.SetStringSelection(user_params.refinement.real_space)
        
        if user_params.refinement.use_refmac_instead_of_phenix:
            tabRef.SoftChoiceReci.SetSelection(1)
            tabRef.SoftChoiceReal.SetSelection(1)

        tabRef.onSoftReciChanged(None)
        tabRef.onSoftRealChanged(None)

        tabRef.wxc_scale_TextCtrl.SetValue(str(user_params.refinement.phenix_keywords.target_weights.wxc_scale))
        tabRef.wxu_scale_TextCtrl.SetValue(str(user_params.refinement.phenix_keywords.target_weights.wxu_scale))
        tabRef.bonds_rmsd_TextCtrl.SetValue(str(user_params.refinement.phenix_keywords.target_weights.weight_selection_criteria.bonds_rmsd))
        tabRef.angle_rmsd_TextCtrl.SetValue(str(user_params.refinement.phenix_keywords.target_weights.weight_selection_criteria.angles_rmsd))
        tabRef.rf_minus_rw.SetValue(str(user_params.refinement.phenix_keywords.target_weights.weight_selection_criteria.r_free_minus_r_work))

        self.set_strategy(user_params.refinement.phenix_keywords.refine.strategy)

        tabRef.NCyclesReciprocal_TextCtrl.SetValue(str(user_params.refinement.phenix_keywords.main.cycles))
        tabRef.ordered_solvent.SetValue(user_params.refinement.phenix_keywords.main.ordered_solvent)
        tabRef.sim_ann.SetValue(user_params.refinement.phenix_keywords.main.simulated_annealing)

        tabRef.start_T.SetValue(str(user_params.refinement.phenix_keywords.simulated_annealing.start_temperature))
        tabRef.final_T.SetValue(str(user_params.refinement.phenix_keywords.simulated_annealing.final_temperature))
        tabRef.cooling_rate.SetValue(str(user_params.refinement.phenix_keywords.simulated_annealing.cool_rate))
        tabRef.mode.SetStringSelection(user_params.refinement.phenix_keywords.simulated_annealing.mode)

        tabRef.map_sharpening.SetValue(user_params.refinement.phenix_keywords.map_sharpening.map_sharpening)

        tabRef.NCyclesReal_TextCtrl.SetValue(str(user_params.refinement.phenix_keywords.real_space_refine.cycles))

        tabRef.density_modification.SetValue(user_params.refinement.phenix_keywords.density_modification.density_modification)
        tabRef.combine.SetStringSelection(user_params.refinement.phenix_keywords.density_modification.combine)
        tabRef.cycles.SetValue(str(user_params.refinement.phenix_keywords.density_modification.cycles))

        tabRef.AUTO.SetStringSelection(user_params.refinement.refmac_keywords.target_weights.weight)
        tabRef.NOEX.SetStringSelection(user_params.refinement.refmac_keywords.target_weights.experimental_sigmas)
        tabRef.weighting_term.SetValue(str(user_params.refinement.refmac_keywords.target_weights.weighting_term))

        tabRef.jelly_body_refinement.SetStringSelection(str(user_params.refinement.refmac_keywords.restraints.jelly_body_refinement))
        tabRef.jbs_TextCtrl.SetValue(str(user_params.refinement.refmac_keywords.restraints.jelly_body_sigma))

        #print(len(user_params.refinement.refmac_keywords.restraints.jelly_body_additional_restraints))
        jb_add_restraints_lst = user_params.refinement.refmac_keywords.restraints.jelly_body_additional_restraints
        if len(jb_add_restraints_lst) == 0:
            jb_add_restraints = "None"
        else:
            jb_add_restraints = ' '.join(str(jb_add_restraints_lst))
        tabRef.jbar_TextCtrl.SetValue(jb_add_restraints)
        external_res_lst = user_params.refinement.refmac_keywords.restraints.external_restraints
        if len(external_res_lst) == 0:
            external_res = "None"
        else:
            external_res = ' '.join(str(external_res_lst))
        tabRef.extr_TextCtrl.SetValue(external_res)
        tabRef.REFType.SetStringSelection(user_params.refinement.refmac_keywords.refine.type)
        tabRef.TLS.SetStringSelection(str(user_params.refinement.refmac_keywords.refine.TLS))
        tabRef.tls_cycles.SetValue(str(user_params.refinement.refmac_keywords.refine.TLS_cycles))
        tabRef.bfac_set.SetValue(str(user_params.refinement.refmac_keywords.refine.bfac_set))
        tabRef.Bref.SetStringSelection(user_params.refinement.refmac_keywords.refine.Brefinement)
        tabRef.refmac_cycles.SetValue(str(user_params.refinement.refmac_keywords.refine.cycles))
        tabRef.twinning.SetValue(user_params.refinement.refmac_keywords.refine.twinning)
        tabRef.REFMAC_map_sharpening.SetValue(user_params.refinement.refmac_keywords.map_sharpening.map_sharpening)
        tabRef.refmac_DM.SetValue(user_params.refinement.refmac_keywords.density_modification.density_modification)
        tabRef.refmac_combine.SetStringSelection(user_params.refinement.refmac_keywords.density_modification.combine)
        tabRef.refmac_DM_cycles.SetValue(str(user_params.refinement.refmac_keywords.density_modification.cycles))

    def set_strategy(self,strategy):
        checkboxes = ['individual_sites', 'individual_sites_real_space', 'rigid_body',
                      'individual_adp', 'group_adp', 'tls', 'occupancies', 'group_anomalous']
        # Resetting all checkboxes to False
        for checkbox in checkboxes:
            getattr(self.notebook.Configure.tabRefine, checkbox).SetValue(False)
        # Applying user strategy
        for step in strategy:
            getattr(self.notebook.Configure.tabRefine, step).SetValue(True)

    def setX8_mode(self, user_params):
        if user_params.output.generate_fofo_only is True:
            self.notebook.Configure.tabExt.X8Modes.SetSelection(0)
            self.notebook.Configure.tabExt.onFoFo()
            self.notebook.Configure.tabExt.currentX8Mode = "FoFo" 
            return
        if user_params.f_and_maps.fast_and_furious is True:
            self.notebook.Configure.tabExt.X8Modes.SetSelection(1)
            self.notebook.Configure.tabExt.onFastNFurious()
            self.notebook.Configure.tabExt.currentX8Mode = "FNF"
            return
        else:
            self.notebook.Configure.tabExt.X8Modes.SetSelection(2)
            self.notebook.Configure.tabExt.onCalmNCurious()
            self.notebook.Configure.tabExt.currentX8Mode = "CNC"

    def fill_map_types(self, user_params):
        checkBoxes = ['qfextr', 'fextr', 'kfextr', 'qfgenick', 'kfgenick', 'fgenick', 'qfextr_calc', 'kfextr_calc',
                      'fextr_calc']
        tabExt = self.notebook.Configure.tabExt

        # Reinitialize all checkboxes
        for checkBox in checkBoxes:
            getattr(tabExt, checkBox).SetValue(False)

        # All maps
        if user_params.f_and_maps.all_maps is True:
            for checkBox in checkBoxes:
                getattr(tabExt, checkBox).SetValue(True)
            return

        only = []

        if user_params.f_and_maps.only_qweight is True:
            for checkBox in checkBoxes:
                if checkBox.startswith('q'):
                    getattr(tabExt, checkBox).SetValue(True)
            only.append(True)

        if user_params.f_and_maps.only_kweight is True:
            for checkBox in checkBoxes:
                if checkBox.startswith('k'):
                    getattr(tabExt, checkBox).SetValue(True)
            only.append(True)

        if user_params.f_and_maps.only_no_weight is True:
            for checkBox in checkBoxes:
                if checkBox.startswith('f'):
                    getattr(tabExt, checkBox).SetValue(True)
            only.append(True)

        if True in only:
            return
        else:
            for fextr in checkBoxes:
                if fextr in user_params.f_and_maps.f_extrapolated_and_maps:
                    getattr(self.notebook.Configure.tabExt, fextr).SetValue(True)
                #else:
                #    getattr(self.notebook.Configure.tabExt, fextr).SetValue(False)

    def onBrowse(self,style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST,wildcard=''):
            """
            File selection.
            :param evt: wx.EVT_BUTTON (guessed)
            :param text_static: str, modified TextStatic widget
            :param multi: bool, is text_static multiline (This could be determined by text_static attribute)
            :param key: str, used to set the value of this key in self.inputs dict
            :param ext: Could be used to filter file extensions (Not a the moment
            :return: None
            """
            dlg = wx.FileDialog(self,
                                message='Open File',
                                defaultDir = '',
                                defaultFile='',
                                wildcard=wildcard,
                                style=style)

            if dlg.ShowModal() == wx.ID_CANCEL:
                return

            return dlg.GetPath()

    def onBrowseDir(self, evt):
        """
        Directory selection
        :param evt: wx.EVT_BUTTON (guessed)
        :param text_static: str, modified TextStatic widget
        :param multi: bool, is text_static multiline (This could be determined by text_static attribute)
        :param key: str, used to set the value of this key in self.inputs dict
        :param ext: Could be used to filter file extensions (Not a the moment
        :return: None
        """
        dlg = wx.DirDialog(self,
                           'Choose Dir',
                           '',
                           wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_CANCEL:
            return

        path = dlg.GetPath()
        return path
    
    def extract_phil(self):
        """
        All inputs will be collected from widgets to build the phil object master_phil
        :return: phil input file
        """
        #master_phil = iotbx.phil.parse(input_string,process_includes=True)

        tabIO = self.notebook.Configure.tabIO

        #TODO: Here check that you have a single ref model, ref mtz and at least a trigg mtz
        try:
            tobeparsed = "input.reference_mtz = %s \n" % tabIO.files["Reference mtz"][0]
        except IndexError:
            tobeparsed = "input.reference_mtz = None \n"
        
        try:
            tobeparsed += "input.triggered_mtz = %s \n" % tabIO.files["Triggered mtz"][0]
        except IndexError:
            tobeparsed += "input.triggered_mtz = None \n"
        try:
            tobeparsed += "input.reference_pdb = %s\n" % tabIO.files["Reference model"][0]
        except IndexError:
            tobeparsed += "input.reference_pdb = None\n" 
        
        for _add_file in tabIO.files["Restraints"]:
            if len(_add_file) > 0: tobeparsed += "input.additional_files = %s\n" % _add_file
        if len(tabIO.highRes.GetValue()) > 0:
            highRes = tabIO.highRes.GetValue()
        else:
            highRes = "None"

        if len(tabIO.lowRes.GetValue()) > 0:
            lowRes = tabIO.lowRes.GetValue()
        else:
            lowRes = "None"

        tobeparsed += "input.high_resolution = %s \n" % highRes + \
                      "input.low_resolution = %s \n" % lowRes

        #######################
        ### Output - IO_tab ###
        #######################
        tobeparsed += "output.outdir = %s\n" % self.get_txtctrl_values(tabIO.outdir_sizer.TextCtrl) + \
                      "output.outname = %s\n" % self.get_txtctrl_values(tabIO.outname) + \
                      "output.GUI = True\noutput.open_coot = False\n"

        #####################
        ### Occ - Ext_tab ###
        #####################
        tabExt = self.notebook.Configure.tabExt
        tobeparsed += "occupancies.low_occ = %s\n" % self.get_txtctrl_values(tabExt.LowTextCtrl) + \
                      "occupancies.high_occ = %s\n" % self.get_txtctrl_values(tabExt.HighTextCtrl) +\
                      "occupancies.steps = %s\n" % self.get_txtctrl_values(tabExt.StepsTextCtrl) +\
                      "occupancies.list_occ = %s\n" % self.get_txtctrl_values(tabExt.ListTextCtrl)
        #######################
        ### \Occ - Ext_tab\ ###
        #######################

        ################################
        ### Maps N Scaling - Ext_tab ###
        ################################
        fofo = ["qfofo", "fofo", "kfofo"]
        choice = tabExt.FoChoice.GetStringSelection()
        fofo_str = self.buildSTR(fofo, choice)

        scaling = ["no", "isotropic", "anisotropic"]
        choice = tabExt.ScalingChoice.GetStringSelection()
        scaling_str = self.buildSTR(scaling, choice)
        tobeparsed += "f_and_maps.fofo_type = %s\n" % fofo_str + \
                      "scaling.b_scaling = %s\n" % scaling_str +\
                      "f_and_maps.kweight_scale = %s\n" % self.get_txtctrl_values(tabExt.kscale)

        X8Mode = tabExt.X8Modes.GetSelection()
        if X8Mode == 1:
            tobeparsed += "f_and_maps.fast_and_furious = True\n"
        if X8Mode == 0:
            tobeparsed += "output.generate_fofo_only = True\n"


        ################################

        #####################
        ### Type of Fextr ###
        #####################
        checkBoxes = ['qfextr', 'fextr', 'kfextr',
                      'qfgenick',  'fgenick',  'kfgenick',
                      'qfextr_calc', 'fextr_calc', 'kfextr_calc']
        FextrSelection =  []
        for checkBox in checkBoxes:
            if getattr(self.notebook.Configure.tabExt, checkBox).IsChecked():
                FextrSelection.append('*%s'%checkBox)
            else:
                FextrSelection.append(checkBox)

        if len(FextrSelection) > 0:
            tobeparsed += "f_and_maps.f_extrapolated_and_maps = %s\n" % ' '.join(FextrSelection)
        else:
            tobeparsed += "f_and_maps.f_extrapolated_and_maps = None\n"
        #    print("You need to at least select one type of extrapolated structure factors to compute")
        #    return
        ################################

        ##############################
        ### Map explorer - Ext_tab ###
        ##############################
        tobeparsed += "map_explorer.peak_integration_floor = %s\n" % self.get_txtctrl_values(tabExt.peak_integration_floorTextCtrl) + \
                      "map_explorer.peak_detection_threshold = %s\n" % self.get_txtctrl_values(tabExt.peak_detection_thresholdTextCtrl) + \
                      "map_explorer.radius = %s\n" % self.get_txtctrl_values(tabExt.RadiusTextCtrl) + \
                      "map_explorer.z_score = %s\n" % self.get_txtctrl_values(tabExt.ZscoreTextCtrl) + \
                      "map_explorer.occupancy_estimation = %s\n" %(tabExt.OccEstimation.GetStringSelection())
                      #"map_explorer.use_occupancy_from_distance_analysis = %s\n" % tabExt.DistanceAnalysis.GetValue()
        #TODO: user_params.map_explorer.use_occupancy_from_distance_analysis = ?
        ##############################

        ###########################################
        ### Neg N Missing Refelctions - Ext tab ###
        ###########################################
        neg = self.notebook.Configure.tabExt.negChoice.GetStringSelection()
        missing = self.notebook.Configure.tabExt.missChoice.GetStringSelection()
        if missing == 'fill': missing = 'and_fill'
        #if neg == '--':
            #if missing == 'and_fill':
                #neg = 'fill'
                #missing = 'missing'
            #else:
                #neg = 'no'
                #missing = 'fill'


        tobeparsed += "f_and_maps.negative_and_missing = %s_%s\n" % (neg, missing)
        ###########################################
        
        ###########################################
        ### Scaling Resolution ###
        ###########################################
        if len(tabExt.ScalingLowRes.GetValue()) > 0:
            SR_low = tabExt.ScalingLowRes.GetValue()
        else:
            SR_low = "None"
        if len(tabExt.ScalingHighRes.GetValue()) > 0:
            SR_high = tabExt.ScalingHighRes.GetValue()
        else:
            SR_high = "None"
        
        tobeparsed += "scaling.high_resolution = %s\n" %(SR_low) +\
            "scaling.low_resolution = %s\n" %(SR_high)
        
        #tobeparsed += "scaling.high_resolution = %s\n" %(tabExt.ScalingHighRes.GetStringSelection()) +\
            #"scaling.low_resolution = %s\n" %(tabExt.ScalingLowRes.GetStringSelection())

        ########################
        ### Phenix - Ref_tab ###
        ########################
        tabRefine = self.notebook.Configure.tabRefine
        
        tobeparsed += "refinement.run_refinement = %s\n" % tabRefine.RunRef.IsChecked() + \
                      "refinement.use_refmac_instead_of_phenix = False\n" +\
                      "refinement.reciprocal_space = %s\n" %tabRefine.SoftChoiceReci.GetStringSelection() +\
                      "refinement.real_space = %s\n" %tabRefine.SoftChoiceReal.GetStringSelection() + \
                      "refinement.phenix_keywords.target_weights.wxc_scale = %s\n" % self.get_txtctrl_values(tabRefine.wxc_scale_TextCtrl) +\
                      "refinement.phenix_keywords.target_weights.wxu_scale = %s\n" % self.get_txtctrl_values(tabRefine.wxu_scale_TextCtrl) +\
                      "refinement.phenix_keywords.target_weights.weight_selection_criteria.bonds_rmsd = %s\n" % self.get_txtctrl_values(tabRefine.bonds_rmsd_TextCtrl) + \
                      "refinement.phenix_keywords.target_weights.weight_selection_criteria.angles_rmsd = %s\n" % self.get_txtctrl_values(tabRefine.angle_rmsd_TextCtrl) + \
                      "refinement.phenix_keywords.target_weights.weight_selection_criteria.r_free_minus_r_work = %s\n" % self.get_txtctrl_values(tabRefine.rf_minus_rw) + \
                      "refinement.phenix_keywords.main.cycles = %s\n" % self.get_txtctrl_values(tabRefine.NCyclesReciprocal_TextCtrl) +\
                      "refinement.phenix_keywords.main.ordered_solvent = %s\n" % tabRefine.ordered_solvent.IsChecked()
        Ref_checkBoxes = ['individual_sites', 'individual_sites_real_space', 'rigid_body', 'individual_adp', 'group_adp', 'tls',
                          'occupancies', 'group_anomalous']
        strategy = ""
        for checkBox in Ref_checkBoxes:
            if getattr(self.notebook.Configure.tabRefine, checkBox).IsChecked():
                strategy += " *%s" % checkBox
            else:
                strategy += " %s" % checkBox

        tobeparsed += "refinement.phenix_keywords.refine.strategy =%s\n" % strategy + \
                      "refinement.phenix_keywords.main.cycles = %s\n" % self.get_txtctrl_values(
            tabRefine.NCyclesReciprocal_TextCtrl) + \
                      "refinement.phenix_keywords.main.ordered_solvent = %s\n" % tabRefine.ordered_solvent.IsChecked() + \
                      "refinement.phenix_keywords.main.simulated_annealing = %s\n" % tabRefine.sim_ann.IsChecked() + \
                      "refinement.phenix_keywords.simulated_annealing.start_temperature = %s\n" % self.get_txtctrl_values(
            tabRefine.start_T) + \
                      "refinement.phenix_keywords.simulated_annealing.final_temperature = %s\n" % self.get_txtctrl_values(
            tabRefine.final_T) + \
                      "refinement.phenix_keywords.simulated_annealing.cool_rate = %s\n" % self.get_txtctrl_values(
            tabRefine.cooling_rate) + \
                      "refinement.phenix_keywords.map_sharpening.map_sharpening = %s\n" % tabRefine.map_sharpening.IsChecked() + \
                      "refinement.phenix_keywords.real_space_refine.cycles = %s\n" % self.get_txtctrl_values(
            tabRefine.NCyclesReal_TextCtrl) + \
                      "refinement.phenix_keywords.density_modification.density_modification = %s\n" % tabRefine.density_modification.IsChecked() + \
                      "refinement.phenix_keywords.density_modification.combine = %s\n" % tabRefine.combine.GetStringSelection() + \
                      "refinement.phenix_keywords.density_modification.cycles = %s\n" % self.get_txtctrl_values(
            tabRefine.cycles)

        SA_mode = ["every_macro_cycle", "second_and_before_last", "once", "first",  "first_half"]
        choice = tabRefine.mode.GetStringSelection()
        mode_str = self.buildSTR(SA_mode, choice)
        tobeparsed += "refinement.phenix_keywords.simulated_annealing.mode = %s\n" % mode_str

        DM_combine = ["PERT", "OMIT"]
        choice = tabRefine.combine.GetStringSelection()
        DM_combine_str = self.buildSTR(DM_combine, choice)
        tobeparsed += "refinement.phenix_keywords.density_modification.combine = %s\n" %DM_combine_str
        ########################
        ### Refmac - Ref_tab ###
        ########################
        weights = ["AUTO",  "MATRIx"]
        choice = tabRefine.AUTO.GetStringSelection()
        AUTO_str = self.buildSTR(weights, choice)

        EXP_SIG = ["NOEX", "EXPE"]
        choice = tabRefine.NOEX.GetStringSelection()
        EXP_SIG_str = self.buildSTR(EXP_SIG, choice)
        tobeparsed += "refinement.refmac_keywords.target_weights.experimental_sigmas = %s \n" %EXP_SIG_str + \
                      "refinement.refmac_keywords.target_weights.weight = %s\n" % AUTO_str + \
                      "refinement.refmac_keywords.target_weights.weighting_term = %s\n" % self.get_txtctrl_values(tabRefine.weighting_term) + \
                      "refinement.refmac_keywords.restraints.external_restraints = %s\n" % self.get_txtctrl_values(tabRefine.extr_TextCtrl) +\
                      "refinement.refmac_keywords.restraints.jelly_body_additional_restraints = %s\n" % self.get_txtctrl_values(tabRefine.jbar_TextCtrl) + \
                      "refinement.refmac_keywords.restraints.jelly_body_refinement = %s\n" % tabRefine.jelly_body_refinement.GetStringSelection() +\
                      "refinement.refmac_keywords.restraints.jelly_body_sigma = %s\n" % self.get_txtctrl_values(tabRefine.jbs_TextCtrl)

        ref_type = ["RESTrained",  "UNREstrained",  "RIGId"]
        choice = tabRefine.REFType.GetStringSelection()
        ref_str = self.buildSTR(ref_type, choice)

        Bref = ["OVERall",  "ISOTropic"]
        choice = tabRefine.Bref.GetStringSelection()
        Bref_str = self.buildSTR(Bref, choice)
        tobeparsed += "refinement.refmac_keywords.refine.type = %s\n" % ref_str + \
                      "refinement.refmac_keywords.refine.cycles = %s\n" % self.get_txtctrl_values(tabRefine.refmac_cycles) + \
                      "refinement.refmac_keywords.refine.TLS =  %s\n" % tabRefine.TLS.GetStringSelection() +\
                      "refinement.refmac_keywords.refine.TLS_cycles = %s\n" % self.get_txtctrl_values(tabRefine.tls_cycles) +\
                      "refinement.refmac_keywords.refine.Brefinement = %s\n" % Bref_str +\
                      "refinement.refmac_keywords.refine.bfac_set = %s\n" % self.get_txtctrl_values(tabRefine.bfac_set) +\
                      "refinement.refmac_keywords.refine.twinning = %s\n" %tabRefine.twinning.IsChecked() +\
                      "refmac_keywords.refine.cycles = %s\n" % self.get_txtctrl_values(tabRefine.refmac_cycles) +\
                      "refinement.refmac_keywords.map_sharpening.map_sharpening = %s\n" %tabRefine.REFMAC_map_sharpening.IsChecked() +\
                      "refinement.refmac_keywords.density_modification.density_modification = %s\n" % tabRefine.refmac_DM.IsChecked() + \
                      "refinement.refmac_keywords.density_modification.combine = %s\n" % tabRefine.refmac_combine.GetStringSelection() + \
                      "refinement.refmac_keywords.density_modification.cycles = %s\n" % self.get_txtctrl_values(tabRefine.refmac_DM_cycles)

            ########################
        user_phil = parse(tobeparsed)
        user_params = master_phil.fetch(source=user_phil).extract()
        return user_params

    def buildSTR(self, lst, choice):
        mode_str = ''
        for mode in lst:
            if mode == choice:
                mode_str += '*%s ' %mode
            else:
                mode_str += mode+' '
        return mode_str
        ########################

    def get_txtctrl_values(self,TxtCtrl):
        value = TxtCtrl.GetValue()
        if len(value) > 0:
            return value
        else: return "None"

# ----------------------------------------------------------------------
if __name__ == "__main__":
    app = wx.App(False)
    frame = MainFrame(sys.argv)
    frame.Show(True)
    app.MainLoop()
