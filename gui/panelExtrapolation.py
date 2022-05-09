"""
 panelExtrapolation.py


 Created 10/21/2020


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

import wx
from wx.lib.pubsub import pub

from utils import CharValidator
from wx.lib.scrolledpanel import ScrolledPanel

class TabExtrapolation(ScrolledPanel):
    """
    This will be the second notebook tab with settings for extrapolation calculation and analysis
    """
    #----------------------------------------------------------------------
    def __init__(self, parent):
        """"""
        ScrolledPanel.__init__(self, parent=parent, style=wx.VSCROLL | wx.HSCROLL)
        self.SetupScrolling()
        self.createAndLayout()
        self.parent = parent
        
    def createAndLayout(self):

        defont = wx.Font(11, wx.MODERN, wx.NORMAL, wx.NORMAL, False, 'MS Shell Dlg 2')
        bigfont = wx.Font(10, wx.MODERN, wx.NORMAL, wx.NORMAL, False, 'MS Shell Dlg 2')
        self.SetFont(defont)
        #####################
        ###   X8 Modes    ###
        #####################
        self.X8Modes = wx.RadioBox(self, 1, "X8 Modes", size=(800, -1),choices=["FoFo only", "Fast N Furious", "Calm N Curious" ])
        self.X8Modes.Bind(wx.EVT_RADIOBOX, self.onRadioBox)
        self.X8Modes.SetSelection(2)
        self.currentX8Mode = "CNC"

        #####################
        ###  Occupancies  ###
        #####################
        width_TextCtrl = 50
        height_TextCtrl = 30
        Occ = wx.StaticBox(self, 1, "Occupancies", size=(800, 200))
        Occ.SetFont(defont)
        occ_sizer = wx.BoxSizer(wx.HORIZONTAL)

        LowTxt = wx.StaticText(self, wx.ID_ANY, "Lowest :", size=(60, -1))
        self.LowTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "0.1", style= wx.TE_PROCESS_ENTER,
                                    size=(width_TextCtrl, height_TextCtrl), validator=CharValidator('no-alpha'))
        LowTxt.SetFont(defont)
        self.LowTextCtrl.SetFont(defont)
        HighTxt = wx.StaticText(self, wx.ID_ANY, "Highest :", size=(60, -1))
        self.HighTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "1", style= wx.TE_PROCESS_ENTER,
                                   size=(width_TextCtrl, height_TextCtrl), validator=CharValidator('no-alpha'))
        HighTxt.SetFont(defont)
        self.HighTextCtrl.SetFont(defont)
        StepsTxt = wx.StaticText(self, wx.ID_ANY, "Steps :", size=(60, -1))
        self.StepsTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "3", style= wx.TE_PROCESS_ENTER,
                                   size=(width_TextCtrl, height_TextCtrl), validator=CharValidator('no-alpha'))
        StepsTxt.SetFont(defont)
        self.StepsTextCtrl.SetFont(defont)

        occ_sizer.Add(LowTxt, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        occ_sizer.Add(self.LowTextCtrl, 0, wx.ALL  | wx.ALIGN_CENTER_VERTICAL, 10)
        occ_sizer.Add(HighTxt, 0, wx.LEFT| wx.ALIGN_CENTER_VERTICAL, 10)
        occ_sizer.Add(self.HighTextCtrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        occ_sizer.Add(StepsTxt, 0, wx.ALL| wx.ALIGN_CENTER_VERTICAL, 10)
        occ_sizer.Add(self.StepsTextCtrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)

        list_occ_sizer = wx.BoxSizer(wx.HORIZONTAL)
        ListTxt = wx.StaticText(self, wx.ID_ANY, "Occupancies List : ", size=(150, -1))
        self.ListTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "", style=wx.TE_PROCESS_ENTER,
                                   size=(120, height_TextCtrl))
        list_occ_sizer.Add(ListTxt, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        list_occ_sizer.Add(self.ListTextCtrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)

        occ_sizer_final = wx.StaticBoxSizer(Occ, wx.VERTICAL)
        occ_sizer_final.Add(occ_sizer, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 0)
        #occ_sizer_final.AddSpacer(20)
        occ_sizer_final.Add(list_occ_sizer, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 0)
        #######################
        ###  \Occupancies\  ###
        #######################


        ########################
        ###  Maps N Scaling  ###
        ########################
        Maps =  wx.StaticBox(self, 1, "Maps And Scaling", size=(800, 200))

        FoTxt = wx.StaticText(self, wx.ID_ANY, "Type of FoFo Map : ")
        FoTxt.SetFont(defont)
        self.FoChoice = wx.Choice(self, wx.ID_ANY, choices = ["qfofo", "fofo", "kfofo"])
        self.FoChoice.SetSelection(0)
        kscale = wx.StaticText(self, wx.ID_ANY, "K weight Scale : ")
        self.kscale = wx.TextCtrl(self, wx.ID_ANY, "0.05")

        #MNS_fgs = wx.FlexGridSizer(2, 4, 10, 10)

        MNS_sizer = wx.BoxSizer(wx.HORIZONTAL)
        MNS_sizer.Add(FoTxt, 0, wx.LEFT  | wx.ALIGN_CENTER_VERTICAL, 10)
        MNS_sizer.Add(self.FoChoice, 0, wx.LEFT  | wx.ALIGN_CENTER_VERTICAL, 10)
        MNS_sizer.AddSpacer(30)
        MNS_sizer.Add(kscale, 0, wx.LEFT | wx.ALIGN_CENTER_VERTICAL, 10)
        MNS_sizer.Add(self.kscale, 0, wx.LEFT | wx.ALIGN_CENTER_VERTICAL, 10)
        ScalingTxt = wx.StaticText(self, wx.ID_ANY, "B-factor Scaling : ", size=(140, -1))
        self.ScalingChoice = wx.Choice(self, wx.ID_ANY, choices=["no", "isotropic","anisotropic"])
        self.ScalingChoice.SetSelection(2)
        ScalingTxt.SetFont(defont)

        S_sizer = wx.BoxSizer(wx.HORIZONTAL)
        S_sizer.Add(ScalingTxt, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        S_sizer.Add(self.ScalingChoice, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)

        MNS_final = wx.StaticBoxSizer(Maps, wx.VERTICAL)
        MNS_final.AddSpacer(10)
        MNS_final.Add(MNS_sizer)
        MNS_final.AddSpacer(20)
        MNS_final.Add(S_sizer)
        ##########################
        ###  \Maps N Scaling\  ###
        ##########################


        self.ExtSF_box = wx.StaticBox(self, 1, "Extrapolated structure factors", size=(800, 200))
        self.ExtSF = wx.StaticBoxSizer(self.ExtSF_box, wx.VERTICAL)

        fgs = wx.FlexGridSizer(3, 5, 10, 10)
        self.qfextr = wx.CheckBox(self, id=wx.ID_ANY,label="qfextr")
        self.fextr =  wx.CheckBox(self, id=wx.ID_ANY, label="fextr")
        self.kfextr = wx.CheckBox(self, id=wx.ID_ANY, label="kfextr")
        self.qfgenick = wx.CheckBox(self, id=wx.ID_ANY, label="qfgenick")
        self.fgenick = wx.CheckBox(self, id=wx.ID_ANY, label="fgenick")
        self.kfgenick = wx.CheckBox(self, id=wx.ID_ANY, label="kfgenick")
        self.qfextr_calc = wx.CheckBox(self, id=wx.ID_ANY, label="qfextr_calc")
        self.fextr_calc = wx.CheckBox(self, id=wx.ID_ANY, label="fextr_calc")
        self.kfextr_calc = wx.CheckBox(self, id=wx.ID_ANY, label="kfextr_calc")
        self.qfextr.SetValue(True)

        self.ID1 = wx.NewId()
        self.ID2 = wx.NewId()
        self.ID3 = wx.NewId()
        self.ID4 = wx.NewId()
        self.ID5 = wx.NewId()

        self.qwFextr = [self.qfextr, self.qfgenick, self.qfextr_calc]
        self.non_qwFextr = [self.fextr, self.fgenick, self.fextr_calc]
        self.kFextr = [self.kfextr, self.kfgenick, self.kfextr_calc]
        self.allfextr_btns = self.qwFextr + self.non_qwFextr + self.kFextr

        self.allfextr = wx.Button(self, id=self.ID1,label="Select all")
        self.none = wx.Button(self, id=self.ID2, label="Clear all")
        self.qonly = wx.Button(self, id=self.ID3, label="All q-weighted")
        self.nonq_only =  wx.Button(self, id=self.ID4, label="All non-weighted")
        self.konly = wx.Button(self, id=self.ID5, label="All k-weighted")

        self.selectButtons = [self.allfextr, self.none, self.qonly, self.nonq_only, self.konly]

        fgs.AddMany([self.qfextr,      self.kfextr,      self.fextr,      self.nonq_only, self.allfextr,
                     self.qfgenick,    self.kfgenick,    self.fgenick,    self.qonly, self.none,
                     self.qfextr_calc, self.kfextr_calc, self.fextr_calc, self.konly])

        self.ExtSF.Add(fgs, 0, wx.EXPAND, border=5)


        self.ExpertMode = wx.Button(self, id=wx.ID_ANY,label="Expert Mode ?")
        
        
        MapExplorer = wx.StaticBox(self, 1, "Map Explorer", size=(800,250))

        self.MES = wx.StaticBoxSizer(MapExplorer,wx.VERTICAL)
        peak_detection_threshold = wx.StaticText(self, wx.ID_ANY, "Peak_detection_threshold (sigma)", size=(90, 20))
        peak_detection_threshold.SetFont(defont)
        peak_integration_floor = wx.StaticText(self, wx.ID_ANY, "Peak_integration_floor (sigma)", size=(120, 20))
        peak_integration_floor.SetFont(defont)
        radius = wx.StaticText(self, wx.ID_ANY, "Radius (A)", size=(90, -1))
        radius.SetFont(defont)
        zscore = wx.StaticText(self, wx.ID_ANY, "z-score", size=(60, -1))
        zscore.SetFont(defont)
        blank = wx.StaticText(self, wx.ID_ANY, "", size=(60, -1))
        blank2 = wx.StaticText(self, wx.ID_ANY, "", size=(60, -1))

        width_TextCtrl = 100
        self.peak_detection_thresholdTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "4.0", style=wx.TE_PROCESS_ENTER,
                                    size=(width_TextCtrl, 20))
        self.peak_integration_floorTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "3.5", style=wx.TE_PROCESS_ENTER,
                                   size=(width_TextCtrl, 20))

        ### Should be changed by highest resolution
        self.RadiusTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "1.5", style=wx.TE_PROCESS_ENTER,
                                   size=(width_TextCtrl, 20))
        self.ZscoreTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "2.0", style=wx.TE_PROCESS_ENTER,
                                     size=(width_TextCtrl, 20))

        MES_fgs = wx.FlexGridSizer(rows=2, cols=5, vgap=10, hgap=10)
        MES_fgs.AddMany([peak_detection_threshold, self.peak_detection_thresholdTextCtrl, blank, peak_integration_floor, self.peak_integration_floorTextCtrl, radius, self.RadiusTextCtrl, blank2, zscore, self.ZscoreTextCtrl])
        self.MES.AddSpacer(5)
        self.MES.Add(MES_fgs)
        self.MES.AddSpacer(10)
        self.DistanceAnalysis = wx.CheckBox(self, wx.ID_ANY, "Use occupancy from distance analysis")
        self.MES.Add(self.DistanceAnalysis)

        NM = wx.StaticBox(self, 1, "Negative and Missing reflections", size=(800, 200))
        self.NM = wx.StaticBoxSizer(NM, wx.VERTICAL)
        neg = wx.StaticText(self, wx.ID_ANY, "Negative", size=(60, -1))
        neg.SetFont(defont)
        missing = wx.StaticText(self, wx.ID_ANY, "Missing", size=(60, -1))
        missing.SetFont(defont)
        self.negChoice = wx.Choice(self, wx.ID_ANY, choices=["truncate", "keep", "reject", "zero", "foff", "fcalc", "--"])
        self.negChoice.SetFont(defont)
        self.negChoice.SetSelection(0)
        self.missChoice = wx.Choice(self, wx.ID_ANY, choices=["fill", "no_fill"])
        self.missChoice.SetFont(defont)
        self.missChoice.SetSelection(0)
        NM_fgs = wx.FlexGridSizer(1, 5, 10, 10)
        NM_fgs.AddMany([neg, self.negChoice, blank, missing, self.missChoice])
        self.NM.AddSpacer(5)
        self.NM.Add(NM_fgs)

        self.FinalSizer = wx.BoxSizer(wx.VERTICAL)
        self.FinalSizer.Add(self.X8Modes, 0, wx.ALL, 10)
        self.FinalSizer.Add(occ_sizer_final, 0, wx.ALL, 10)
        self.FinalSizer.Add(MNS_final, 0, wx.ALL, 10)
        self.FinalSizer.Add(self.ExtSF, 0, wx.ALL, 10)
        self.FinalSizer.Add(self.ExpertMode, 0, wx.ALL, 10)
        self.FinalSizer.Add(self.MES, 0, wx.ALL, 10)
        self.FinalSizer.Hide(self.MES)
        self.FinalSizer.Add(self.NM, 0, wx.ALL, 10)
        self.FinalSizer.Hide(self.NM)
        self.SetSizer(self.FinalSizer)
        self.SetAutoLayout(True)

        self.allfextr.Bind(wx.EVT_BUTTON, self.onSelectFextr)
        self.none.Bind(wx.EVT_BUTTON, self.onSelectFextr)
        self.qonly.Bind(wx.EVT_BUTTON, self.onSelectFextr)
        self.nonq_only.Bind(wx.EVT_BUTTON, self.onSelectFextr)
        self.konly.Bind(wx.EVT_BUTTON, self.onSelectFextr)
        self.ExpertMode.Bind(wx.EVT_BUTTON, self.OnExpertMode)
        self.FoChoice.Bind(wx.EVT_CHOICE, self.updateKScale)
        self.kscale.Disable()

    def onSelectFextr(self, evt):
        id = evt.GetId()
        if id == self.ID1:
            for checkbox in self.allfextr_btns:
                checkbox.SetValue(True)
        elif id == self.ID2:
            for checkbox in self.allfextr_btns:
                checkbox.SetValue(False)
        elif id == self.ID3:
            for checkbox in self.qwFextr:
                checkbox.SetValue(True)
            #for checkbox in self.non_qwFextr:
            #    checkbox.SetValue(False)
            #for checkbox in self.kFextr:
            #    checkbox.SetValue(False)
        elif id == self.ID4:
            #for checkbox in self.qwFextr:
            #    checkbox.SetValue(False)
            for checkbox in self.non_qwFextr:
                checkbox.SetValue(True)
            #for checkbox in self.kFextr:
            #    checkbox.SetValue(False)
        elif id == self.ID5:
            #for checkbox in self.qwFextr:
            #    checkbox.SetValue(False)
            #for checkbox in self.non_qwFextr:
            #    checkbox.SetValue(False)
            for checkbox in self.kFextr:
                checkbox.SetValue(True)

        #self.GetFextr()

    #def GetFextr(self):
    #    fextr = []
    #    for checkbox in self.allfextr_btns:
    #        fextr.append(checkbox.GetValue())

    def OnExpertMode(self,evt):
        if not self.FinalSizer.IsShown(self.MES):
            self.FinalSizer.Show(self.MES)
            self.FinalSizer.Show(self.NM)
            self.ExpertMode.SetLabel("Expert Mode")
            self.FinalSizer.Layout()
        else:
            self.FinalSizer.Hide(self.MES)
            self.FinalSizer.Hide(self.NM)
            self.ExpertMode.SetLabel("Expert Mode ?")
        self.FinalSizer.Layout()

    def onRadioBox(self, evt):
        idx = self.X8Modes.GetSelection()
        pub.sendMessage("X8Mode", mode=(self.currentX8Mode, idx))
        if idx == 2:
            self.onCalmNCurious()
        elif idx == 1:
            self.onFastNFurious()
        elif idx == 0:
            self.onFoFo()

    def onCalmNCurious(self):
        if not self.FinalSizer.IsShown(self.ExtSF):
            self.FinalSizer.Show(self.ExtSF)
            self.FinalSizer.Layout()
        self.FoChoice.Enable()
        self.updateKScale(None)
        for checkbox in self.allfextr_btns:
            checkbox.Enable()
        for btn in self.selectButtons:
            btn.Enable()

    def onFastNFurious(self):
        if not self.FinalSizer.IsShown(self.ExtSF):
            self.FinalSizer.Show(self.ExtSF)
            self.FinalSizer.Layout()
        self.FoChoice.SetStringSelection("qfofo")
        self.FoChoice.Disable()
        self.kscale.Disable()
        for checkbox in self.allfextr_btns:
            checkbox.SetValue(False)
            checkbox.Disable()
        for btn in self.selectButtons:
            btn.Disable()
        self.qfextr.SetValue(True)

    def onFoFo(self):
        if self.FinalSizer.IsShown(self.ExtSF):
            self.FinalSizer.Hide(self.ExtSF)
            self.FinalSizer.Layout()
        self.FoChoice.Enable()
        self.updateKScale(None)

    def updateKScale(self, evt):
        if self.FoChoice.GetStringSelection() == 'kfofo':
            self.kscale.Enable()
        else:
            self.kscale.Disable()

