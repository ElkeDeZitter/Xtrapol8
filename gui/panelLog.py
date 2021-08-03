#----------------------------------------------------------------------
# panelIO.py
#
# Created 10/21/2020
#
# Author: Nicolas Coquelle - nicolas.coquelle@esrf.fr
#
#----------------------------------------------------------------------
import numpy as np
import wx
import os
from wx.lib.scrolledpanel import ScrolledPanel
from wxtbx import metallicbutton
from wx.lib.pubsub import pub
import glob

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
        bmp = wx.Image('gui/pngs/coot_logo_small.png', type=wx.BITMAP_TYPE_ANY, index=-1).ConvertToBitmap()

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
        ScrolledPanel.__init__(self, parent=parent, style=wx.VSCROLL)
        self.SetupScrolling()
        self.parent = parent
        self.mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.mainSizer)
        self.SetAutoLayout(1)
        self.photoMaxSize = 1000

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

        self.newimg = wx.StaticBitmap(self, wx.ID_ANY,
                                      wx.BitmapFromImage(img))
        self.ImgSizer.Add(self.newimg, proportion=0,  flag=wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, border=20)
        self.FitInside()

    def Clear(self, evt):
        while not self.ImgSizer.IsEmpty():
            window = self.ImgSizer.GetItem(0).GetWindow()
            window.Destroy()
        pub.sendMessage("updateFextr", evt=None)#,tabindex=index)


class TabOccResults(ScrolledPanel):
    """
     This will be the first tab of the Configure Notebook (which holds all X8 inputs)
     """

    # ----------------------------------------------------------------------
    def __init__(self, parent, options):
        ScrolledPanel.__init__(self, parent=parent, style=wx.VSCROLL)
        self.SetupScrolling()
        self.mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.mainSizer)
        self.SetAutoLayout(1)
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
            self.occ_list = [str(choice) for choice in np.linspace(self.options.occupancies.low_occ,
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
        bmp = wx.Image('gui/pngs/coot_logo_small.png', type=wx.BITMAP_TYPE_ANY, index=-1).ConvertToBitmap()
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
        self.mainSizer.Add(self.ImgSizer, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL)
        self.mainSizer.Layout()
        
    def setupBindings(self):
        self.Bind(wx.EVT_CHOICE, self.onSelection, self.OccChoice)
        self.Bind(wx.EVT_CHOICE, self.onSelection, self.FextrChoice)

    def onSelection(self, evt):
        self.occ = self.OccChoice.GetStringSelection()
        self.fextr = self.FextrChoice.GetStringSelection()
        
        while not self.ImgSizer.IsEmpty():
            item = self.ImgSizer.GetItem(0)
            window = item.GetWindow()
            try:
                window.Destroy()
            except: item.Destroy()

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


        fn_neg = os.path.join(path, '%s_negative_reflections.png' % Fextr)
        if os.path.isfile(fn_neg):
            self.addImg(fn_neg)

        FsigF = os.path.join(path, '%s_occupancy%.3f.png_FsigF.png' % (Fextr, float(self.occ)))
        if os.path.isfile(FsigF):
            self.addImg(FsigF)
        
        if self.finished:
            self.best_occ_Static.SetLabel("best estimantion @ %s"%self.best_occ[self.fextr])#    self.OccChoice.FindString(s)
            if self.occ == self.best_occ[self.fextr]:
                self.addImg(self.ddm[self.fextr])
                if not self.coot_button.IsShown():
                    self.occNfextrSizer.Show(self.coot_button)
            else:
                if self.coot_button.IsShown():
                    self.occNfextrSizer.Hide(self.coot_button)
                    
        evt.Skip()

    def onCoot(self, evt):
        script_coot = os.path.join(self.coot_scripts[self.fextr], 'coot_all_%s.py' % self.Fextr_png_name)
        os.system("coot --script %s &" % (script_coot))

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
            self.FitInside()

    def onFinished(self):
        fextrs = self.options.f_and_maps.f_extrapolated_and_maps
        for fextr in fextrs:
            if fextr.startswith('q'):
                paths = glob.glob(os.path.join(self.options.output.outdir,'qweight_occupancy*'))
                self.best_occ[fextr] = self.get_best_occupancy(paths, fextr)
            elif fextr.startswith('k'):
                paths = glob.glob(os.path.join(self.options.output.outdir,'kweight_occupancy*'))
                self.best_occ[fextr] = self.get_best_occupancy(paths, fextr)
            else:
                paths = glob.glob(os.path.join(self.options.output.outdir,'occupancy*'))
                self.best_occ[fextr] = self.get_best_occupancy(paths, fextr)
        fextr = self.FextrChoice.GetStringSelection()
        self.occNfextrSizer.Show(self.best_occ_Static)
        self.best_occ_Static.SetLabel("best estimantion @ %s"%self.best_occ[fextr])

    def get_best_occupancy(self, paths, fextr):
        for path in paths:
            ddms = glob.glob(os.path.join(path, 'ddm*png'))
            if len(ddms) > 0:
                for ddm in ddms:
                    if fextr[2:] in ddm:
                        self.ddm[fextr] = ddm
                        self.coot_scripts[fextr] = os.path.dirname(ddm)
                        s = ddm.split('occupancy_')[1][0:5].strip()
                        return s                        
                        #self.OccChoice.SetString(index, "* %s" %s)
                        #Add ddm image











