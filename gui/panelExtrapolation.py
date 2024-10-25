# -*- coding: utf-8 -*-
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
        
        width_TextCtrl = 100
        height_TextCtrl = 20
        
        width_ChoiceCtrl = 100
        height_ChoiceCtrl = 30
        
        
        #blanks for filling spaces in fgs.AddMany
        blank = wx.StaticText(self, wx.ID_ANY, "", size=(1, -1))
        blank1 = wx.StaticText(self, wx.ID_ANY, "", size=(1, -1))
        blank2 = wx.StaticText(self, wx.ID_ANY, "", size=(60, -1))
        blank3 = wx.StaticText(self, wx.ID_ANY, "", size=(60, -1))
        blank4 = wx.StaticText(self, wx.ID_ANY, "", size=(60, -1))
        blank5 = wx.StaticText(self, wx.ID_ANY, "", size=(60, -1))

        
        #####################
        ###   X8 Modes    ###
        #####################
        
        self.X8Modes = wx.RadioBox(self, 1, "Xtrapol8 Modes", size=(800, -1),choices=["FoFo only ", "Fast and Furious ", "Calm and Curious " ])
        self.X8Modes.Bind(wx.EVT_RADIOBOX, self.onRadioBox)
        self.X8Modes.SetSelection(2)
        self.currentX8Mode = "CNC"

        #####################
        ###  Occupancies  ###
        #####################
        

        Occ = wx.StaticBox(self, 1, "Occupancies", size=(800, 200))
        Occ.SetFont(defont)
        occ_sizer = wx.BoxSizer(wx.HORIZONTAL)

        LowTxt = wx.StaticText(self, wx.ID_ANY, "Lowest :", size=(60, -1))
        self.LowTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "0.1", style= wx.TE_PROCESS_ENTER,
                                    size=(50, height_TextCtrl), validator=CharValidator('no-alpha'))
        LowTxt.SetFont(defont)
        self.LowTextCtrl.SetFont(defont)
        HighTxt = wx.StaticText(self, wx.ID_ANY, "Highest :", size=(60, -1))
        self.HighTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "0.5", style= wx.TE_PROCESS_ENTER,
                                   size=(50, height_TextCtrl), validator=CharValidator('no-alpha'))
        HighTxt.SetFont(defont)
        self.HighTextCtrl.SetFont(defont)
        StepsTxt = wx.StaticText(self, wx.ID_ANY, "Steps :", size=(60, -1))
        self.StepsTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "5", style= wx.TE_PROCESS_ENTER,
                                   size=(50, height_TextCtrl), validator=CharValidator('no-alpha'))
        StepsTxt.SetFont(defont)
        self.StepsTextCtrl.SetFont(defont)

        occ_sizer.Add(LowTxt, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        occ_sizer.Add(self.LowTextCtrl, 0, wx.ALL  | wx.ALIGN_CENTER_VERTICAL, 10)
        occ_sizer.Add(HighTxt, 0, wx.LEFT| wx.ALIGN_CENTER_VERTICAL, 10)
        occ_sizer.Add(self.HighTextCtrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        occ_sizer.Add(StepsTxt, 0, wx.ALL| wx.ALIGN_CENTER_VERTICAL, 10)
        occ_sizer.Add(self.StepsTextCtrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)

        list_occ_sizer = wx.BoxSizer(wx.HORIZONTAL)
        ListTxt = wx.StaticText(self, wx.ID_ANY, "Occupancies List :", size=(150, -1))
        self.ListTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "", style=wx.TE_PROCESS_ENTER,
                                   size=(120, height_TextCtrl))
        list_occ_sizer.Add(ListTxt, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        list_occ_sizer.Add(self.ListTextCtrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)

        self.occ_sizer_final = wx.StaticBoxSizer(Occ, wx.VERTICAL)
        self.occ_sizer_final.Add(occ_sizer, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 0)
        #self.occ_sizer_final.AddSpacer(20)
        self.occ_sizer_final.Add(list_occ_sizer, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 0)

        ########################
        ###  FoFo (earlier called Maps and Scaling)  ###
        ########################
        
        Maps =  wx.StaticBox(self, 1, "Fourier difference map (FoFo)", size=(800, 200))

        FoTxt = wx.StaticText(self, wx.ID_ANY, "Type of FoFo Map :")
        FoTxt.SetFont(defont)
        self.FoChoice = wx.Choice(self, wx.ID_ANY, choices = ["qfofo", "fofo", "kfofo"])
        self.FoChoice.SetSelection(0)
        kscale = wx.StaticText(self, wx.ID_ANY, "K weight Scale :")
        self.kscale = wx.TextCtrl(self, wx.ID_ANY, "0.05")

        #MNS_fgs = wx.FlexGridSizer(2, 4, 10, 10)

        MNS_sizer = wx.BoxSizer(wx.HORIZONTAL)
        MNS_sizer.Add(FoTxt, 0, wx.LEFT  | wx.ALIGN_CENTER_VERTICAL, 10)
        MNS_sizer.Add(self.FoChoice, 0, wx.LEFT  | wx.ALIGN_CENTER_VERTICAL, 10)
        MNS_sizer.AddSpacer(30)
        MNS_sizer.Add(kscale, 0, wx.LEFT | wx.ALIGN_CENTER_VERTICAL, 10)
        MNS_sizer.Add(self.kscale, 0, wx.LEFT | wx.ALIGN_CENTER_VERTICAL, 10)
        
        # ScalingTxt = wx.StaticText(self, wx.ID_ANY, "B-factor Scaling : ", size=(140, -1))
        # self.ScalingChoice = wx.Choice(self, wx.ID_ANY, choices=["no", "isotropic","anisotropic"])
        # self.ScalingChoice.SetSelection(2)
        # ScalingTxt.SetFont(defont)

        # S_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # S_sizer.Add(ScalingTxt, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        # S_sizer.Add(self.ScalingChoice, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        
        MNS_final = wx.StaticBoxSizer(Maps, wx.VERTICAL)
        # MNS_final.AddSpacer(10)
        MNS_final.Add(MNS_sizer)
        # MNS_final.AddSpacer(20)
        # MNS_final.Add(S_sizer)

        ##########################
        ###  ESFAs  ###
        ##########################
        
        self.ExtSF_box = wx.StaticBox(self, 1, "Extrapolated structure factor amplitudes", size=(800, 200))
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

        self.ExtSF.AddSpacer(5)
        self.ExtSF.Add(fgs, 0, wx.EXPAND, border=5)

        ##########################
        ###  Map-explorer  ###
        ##########################
        
        MapExplorer = wx.StaticBox(self, 1, "Map-explorer", size=(800,250))

        self.MES = wx.StaticBoxSizer(MapExplorer,wx.VERTICAL)
        peak_detection_threshold = wx.StaticText(self, wx.ID_ANY, "Peak_detection_threshold (sigma) :", size=(250, 20))
        peak_detection_threshold.SetFont(defont)
        peak_integration_floor = wx.StaticText(self, wx.ID_ANY, "Peak_integration_floor (sigma) :", size=(250, 20))
        peak_integration_floor.SetFont(defont)
        radius = wx.StaticText(self, wx.ID_ANY, "Radius (A) :", size=(250, 20))
        radius.SetFont(defont)
        zscore = wx.StaticText(self, wx.ID_ANY, "z-score :", size=(250, 20))
        zscore.SetFont(defont)

        self.peak_detection_thresholdTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "4.0", style=wx.TE_PROCESS_ENTER,
                                    size=(width_TextCtrl, height_TextCtrl))
        self.peak_integration_floorTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "3.5", style=wx.TE_PROCESS_ENTER,
                                   size=(width_TextCtrl, height_TextCtrl))

        ### Should be changed by highest resolution
        self.RadiusTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "1.5", style=wx.TE_PROCESS_ENTER,
                                   size=(width_TextCtrl, height_TextCtrl))
        self.ZscoreTextCtrl = wx.TextCtrl(self, wx.ID_ANY, "2.0", style=wx.TE_PROCESS_ENTER,
                                     size=(width_TextCtrl, height_TextCtrl))
        
        self.occ_est = wx.StaticText(self, wx.ID_ANY, "Occupancy estimation :", size=(180, 20))
        self.occ_est.SetFont(defont)
        self.occ_list_all = ["difference_map_maximization", "difference_map_PearsonCC", "distance_analysis"]
        self.occ_list_noref = ["difference_map_maximization", "difference_map_PearsonCC"]
        self.OccEstimation = wx.Choice(self, wx.ID_ANY, choices=self.occ_list_all, name = "Occupancy estimation")
        self.OccEstimation.SetFont(defont)
        self.OccEstimation.SetSelection(0)

        MES_fgs = wx.FlexGridSizer(rows=3, cols=5, vgap=10, hgap=10)
        MES_fgs.AddMany([peak_detection_threshold, self.peak_detection_thresholdTextCtrl, blank, peak_integration_floor, self.peak_integration_floorTextCtrl,
                         radius, self.RadiusTextCtrl, blank1, zscore, self.ZscoreTextCtrl, 
                         self.occ_est, self.OccEstimation])
        self.MES.AddSpacer(5)
        self.MES.Add(MES_fgs, 0, wx.EXPAND, border=5)
        # self.MES.AddSpacer(10)
        #self.MES.Add(self.OccEstimation)
        #self.DistanceAnalysis = wx.CheckBox(self, wx.ID_ANY, "Use occupancy from distance analysis")
        #self.MES.Add(self.DistanceAnalysis)
        
        ##########################
        ###  Expert mode  ###
        ##########################
        
        self.ExpertMode = wx.Button(self, id=wx.ID_ANY,label="Expert Mode ?")
        self.expert = False

        ##########################
        ###  Negative and missing  ###
        ##########################
        
        NM = wx.StaticBox(self, 1, "Negative and Missing reflections", size=(800, 200))
        self.NM = wx.StaticBoxSizer(NM, wx.VERTICAL)
        neg = wx.StaticText(self, wx.ID_ANY, "Negative :", size=(200, 20))
        neg.SetFont(defont)
        missing = wx.StaticText(self, wx.ID_ANY, "Missing :", size=(200, 20))
        missing.SetFont(defont)
        self.negChoice = wx.Choice(self, wx.ID_ANY, choices=["truncate", "reject", "zero", "fref", "fcalc", "keep", "absolute"], size=(width_ChoiceCtrl, height_ChoiceCtrl))
        self.negChoice.SetFont(defont)
        self.negChoice.SetSelection(0)
        self.missChoice = wx.Choice(self, wx.ID_ANY, choices=["fill", "no_fill"], size=(width_ChoiceCtrl, height_ChoiceCtrl) )
        self.missChoice.SetFont(defont)
        self.missChoice.SetSelection(0)
        NM_fgs = wx.FlexGridSizer(rows=1, cols=5, vgap=10, hgap=10)
        NM_fgs.AddMany([neg, self.negChoice, blank2, missing, self.missChoice])
        self.NM.AddSpacer(5)
        self.NM.Add(NM_fgs, 0, wx.EXPAND, border=5)
        
        ########################
        ### Scaling (earlier called Scaling resolution) ###
        ########################
        
        SR = wx.StaticBox(self, 1, "Scaling", size=(800, 200))
        self.SR = wx.StaticBoxSizer(SR, wx.VERTICAL)
        
        ScalingTxt = wx.StaticText(self, wx.ID_ANY, "B-factor Scaling :", size=(200, 20))
        self.ScalingChoice = wx.Choice(self, wx.ID_ANY, choices=["no", "isotropic","anisotropic"], size=(width_ChoiceCtrl, height_ChoiceCtrl))
        self.ScalingChoice.SetSelection(2)
        ScalingTxt.SetFont(defont)
        
        ScalingHighResTxt = wx.StaticText(self, wx.ID_ANY, label="Scaling high resolution :", size=(200, 20))
        ScalingHighResTxt.SetFont(defont)
        self.ScalingHighRes = wx.TextCtrl(self, wx.ID_ANY, "", style=wx.TE_PROCESS_ENTER, size=(width_TextCtrl, height_TextCtrl))
        ScalingLowResTxt = wx.StaticText(self, wx.ID_ANY, label="Scaling low resolution :", size=(200, 20))
        ScalingLowResTxt.SetFont(defont)
        self.ScalingLowRes = wx.TextCtrl(self, wx.ID_ANY, "", style=wx.TE_PROCESS_ENTER, size=(width_TextCtrl, height_TextCtrl))
        
        SR_sizer = wx.BoxSizer(wx.HORIZONTAL)
        SR_sizer.Add(ScalingHighResTxt, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        SR_sizer.Add(self.ScalingHighRes, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        # SR_sizer.AddSpacer(30)
        SR_sizer.Add(ScalingLowResTxt, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        SR_sizer.Add(self.ScalingLowRes, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        
        SR_fgs = wx.FlexGridSizer(rows=2, cols=5, vgap=10, hgap=10)
        SR_fgs.AddMany([ScalingTxt, self.ScalingChoice, blank2, blank3, blank4, 
                        ScalingLowResTxt, self.ScalingLowRes, blank5, ScalingHighResTxt, self.ScalingHighRes])
        self.SR.AddSpacer(5)
        self.SR.Add(SR_fgs, 0, wx.EXPAND, border=5)

                
        ########################
        ### Final layout ###
        ########################

        self.FinalSizer = wx.BoxSizer(wx.VERTICAL)
        self.FinalSizer.Add(self.X8Modes, 0, wx.ALL, 10)
        self.FinalSizer.Add(self.occ_sizer_final, 0, wx.ALL, 10)
        self.FinalSizer.Add(MNS_final, 0, wx.ALL, 10)
        self.FinalSizer.Add(self.ExtSF, 0, wx.ALL, 10)
        self.FinalSizer.Add(self.MES, 0, wx.ALL, 10)
        self.FinalSizer.Add(self.ExpertMode, 0, wx.ALL, 10)
        # self.FinalSizer.Hide(self.MES)
        self.FinalSizer.Add(self.NM, 0, wx.ALL, 10)
        self.FinalSizer.Hide(self.NM)
        self.FinalSizer.Add(self.SR, 0, wx.ALL, 10)
        self.FinalSizer.Hide(self.SR)
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
        if not self.FinalSizer.IsShown(self.SR):
            # self.FinalSizer.Show(self.MES)
            self.FinalSizer.Show(self.NM)
            self.FinalSizer.Show(self.SR)
            self.ExpertMode.SetLabel("Expert Mode")
            #self.FinalSizer.Layout()
            self.expert = True
        else:
            # self.FinalSizer.Hide(self.MES)
            self.FinalSizer.Hide(self.NM)
            self.FinalSizer.Hide(self.SR)
            self.ExpertMode.SetLabel("Expert Mode ?")
            self.expert = False
        self.FinalSizer.Layout()
        
        idx = self.X8Modes.GetSelection()
        if idx == 2:
            self.onCalmNCurious()
        elif idx == 1:
            self.onFastNFurious()
        elif idx == 0:
            self.onFoFo()

    def onRadioBox(self, evt):
        idx = self.X8Modes.GetSelection()
        #pub.sendMessage("X8Mode", mode=(self.currentX8Mode, idx))
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
        if self.expert == True:
            if not self.FinalSizer.IsShown(self.NM):
                self.FinalSizer.Show(self.NM)
                self.FinalSizer.Layout()
        else:
            self.FinalSizer.Hide(self.NM)
            self.FinalSizer.Layout()
        self.negChoice.Enable()
        self.missChoice.Enable()
        #self.DistanceAnalysis.Enable() 
        OccEst_ini = self.OccEstimation.GetStringSelection()
        self.OccEstimation.SetItems(self.occ_list_all)
        if OccEst_ini in self.occ_list_all:
            self.OccEstimation.SetStringSelection(OccEst_ini)
        if self.FinalSizer.IsShown(self.MES):
            self.OccEstimation.Show()
            self.occ_est.Show()
        if not self.FinalSizer.IsShown(self.occ_sizer_final):
            self.FinalSizer.Show(self.occ_sizer_final)
            self.FinalSizer.Layout()

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
        if self.expert == True:
            if not self.FinalSizer.IsShown(self.NM):
                self.FinalSizer.Show(self.NM)
                self.FinalSizer.Layout()
        else:
            self.FinalSizer.Hide(self.NM)
            self.FinalSizer.Layout()
        self.negChoice.SetStringSelection("truncate")
        self.negChoice.Disable()
        self.missChoice.SetStringSelection("fill")
        self.missChoice.Disable()
        #if not self.DistanceAnalysis.IsShown():
            #self.DistanceAnalysis.Show()
        #self.DistanceAnalysis.SetValue(False)
        #self.DistanceAnalysis.Disable()
        OccEst_ini = self.OccEstimation.GetStringSelection()
        self.OccEstimation.SetItems(self.occ_list_noref)
        if OccEst_ini in self.occ_list_noref:
            self.OccEstimation.SetStringSelection(OccEst_ini)
        if self.FinalSizer.IsShown(self.MES):
            self.OccEstimation.Show()
            self.occ_est.Show()
        if not self.FinalSizer.IsShown(self.occ_sizer_final):
            self.FinalSizer.Show(self.occ_sizer_final)
            self.FinalSizer.Layout()

    def onFoFo(self):
        if self.FinalSizer.IsShown(self.ExtSF):
            self.FinalSizer.Hide(self.ExtSF)
            self.FinalSizer.Layout()
        self.FoChoice.Enable()
        self.updateKScale(None)
        if self.FinalSizer.IsShown(self.NM):
            self.FinalSizer.Hide(self.NM)
            self.FinalSizer.Layout()
        #self.DistanceAnalysis.SetValue(False)
        self.OccEstimation.Hide()
        self.occ_est.Hide()
        
        if self.FinalSizer.IsShown(self.occ_sizer_final):
            self.FinalSizer.Hide(self.occ_sizer_final)
            self.FinalSizer.Layout()

    def updateKScale(self, evt):
        if self.FoChoice.GetStringSelection() == 'kfofo':
            self.kscale.Enable()
        else:
            self.kscale.Disable()

