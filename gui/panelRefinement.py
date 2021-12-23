"""
 panelRefinement.py

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
from wx.lib.scrolledpanel import ScrolledPanel

class TabRefinement(ScrolledPanel):
    """
    This will be the third notebook tab
    """

    # ----------------------------------------------------------------------
    def __init__(self, parent):
        """"""
        ScrolledPanel.__init__(self, parent=parent, style=wx.VSCROLL | wx.HSCROLL)
        self.SetupScrolling()
        self.createAndLayout()

    def createAndLayout(self):

        width_TextCtrl = 50
        defont = wx.Font(11, wx.MODERN, wx.NORMAL, wx.NORMAL, False, 'MS Shell Dlg 2')
        self.SetFont(defont)

        # Perform refinement and soft selection
        soft_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.RunRef = wx.CheckBox(self, wx.ID_ANY, "Perform Refinement with software")
        self.SoftChoice = wx.Choice(self, wx.ID_ANY, choices=['Phenix','Refmac'])
        self.RunRef.SetValue(True)  # Default
        self.SoftChoice.SetSelection(0)  # Default
        soft_sizer.Add(self.RunRef, 0, wx.ALIGN_CENTER_VERTICAL)
        soft_sizer.Add(self.SoftChoice, 0, wx.ALIGN_CENTER_VERTICAL)


        # Refinement Strategy
        strategy = wx.StaticBox(self, 1, "Refinement Strategy", size=(800, 300))
        self.strategy = wx.StaticBoxSizer(strategy, wx.VERTICAL)

        strategy_fgs = wx.FlexGridSizer(rows=4, cols=4, hgap=10, vgap=10)
        self.individual_sites = wx.CheckBox(self, wx.ID_ANY, "XYZ (Reciprocal space)")
        self.individual_sites_real_space = wx.CheckBox(self, wx.ID_ANY, "XYZ (Real space)")
        self.rigid_body = wx.CheckBox(self, wx.ID_ANY, "Rigid Body")
        self.individual_adp = wx.CheckBox(self, wx.ID_ANY, "Individual ADP")
        self.group_adp = wx.CheckBox(self, wx.ID_ANY, "Group ADP")
        self.tls = wx.CheckBox(self, wx.ID_ANY, "TLS")
        self.occupancies = wx.CheckBox(self, wx.ID_ANY, "Occupancies")
        self.group_anomalous = wx.CheckBox(self, wx.ID_ANY, "Group anomalous")
        
        NCyclesReciprocal = wx.StaticText(self, wx.ID_ANY, "Cycles (Reciprocal Space) :")
        NCyclesReal = wx.StaticText(self, wx.ID_ANY, "Cycles (Real Space) :")
        self.NCyclesReciprocal_TextCtrl = wx.TextCtrl(self, wx.ID_ANY, "3", style=wx.TE_PROCESS_ENTER,
                                        size=(width_TextCtrl, 30))
        self.NCyclesReal_TextCtrl = wx.TextCtrl(self, wx.ID_ANY, "3", style=wx.TE_PROCESS_ENTER,
                                                      size=(width_TextCtrl, 30))


        border=10
        strategy_fgs.AddMany([(self.individual_sites, 0, wx.ALIGN_CENTER_VERTICAL, border),
                              (self.individual_sites_real_space, 0, wx.ALIGN_CENTER_VERTICAL, border),
                              (self.rigid_body, 0, wx.ALIGN_CENTER_VERTICAL, border),
                              (self.individual_adp, 0, wx.ALIGN_CENTER_VERTICAL, border),
                              (self.group_adp, 0, wx.ALIGN_CENTER_VERTICAL, border),
                              (self.tls, 0, wx.ALIGN_CENTER_VERTICAL, border),
                              (self.occupancies, 0, wx.ALIGN_CENTER_VERTICAL, border),
                              (self.group_anomalous, 0, wx.ALIGN_CENTER_VERTICAL, border),
                              (NCyclesReciprocal, 0, wx.ALIGN_CENTER_VERTICAL, border),
                              (self.NCyclesReciprocal_TextCtrl, 0, wx.ALIGN_CENTER_VERTICAL, border),
                              (NCyclesReal, 0, wx.ALIGN_CENTER_VERTICAL, border),
                              (self.NCyclesReal_TextCtrl, 0, wx.ALIGN_CENTER_VERTICAL, border)])
        self.strategy.Add(strategy_fgs, 0, wx.ALL, border)
        self.individual_adp.SetValue(True) #default
        self.individual_sites.SetValue(True) #default

        # Weights and targets
        weights = wx.StaticBox(self, 1, "Weights and targets", size=(800, -1))
        self.weights = wx.StaticBoxSizer(weights, wx.VERTICAL)
        weights_GridBag = wx.GridBagSizer(5, 10)#, vgap=10, hgap=10)

        wxc_scale = wx.StaticText(self, wx.ID_ANY, "wxc_scale :", size=(100, -1))
        wxu_scale = wx.StaticText(self, wx.ID_ANY, "wxu_scale :", size=(100, -1))
        bonds_rmsd = wx.StaticText(self, wx.ID_ANY, "Bonds rmsd :", size=(100, -1))
        angle_rmsd = wx.StaticText(self, wx.ID_ANY, "Angles rmsd :", size=(105, -1))
        Rf_minus_Rw = wx.StaticText(self, wx.ID_ANY, "Rfree - Rwork difference :", size=(190, -1))

        self.wxc_scale_TextCtrl = wx.TextCtrl(self, wx.ID_ANY, "0.5", style=wx.TE_PROCESS_ENTER, size=(width_TextCtrl, 30))
        self.wxu_scale_TextCtrl = wx.TextCtrl(self, wx.ID_ANY, "1.0", style=wx.TE_PROCESS_ENTER, size=(width_TextCtrl, 30))
        self.bonds_rmsd_TextCtrl = wx.TextCtrl(self, wx.ID_ANY, "None", style=wx.TE_PROCESS_ENTER, size=(width_TextCtrl, 30))
        self.angle_rmsd_TextCtrl = wx.TextCtrl(self, wx.ID_ANY, "None", style=wx.TE_PROCESS_ENTER, size=(width_TextCtrl, 30))
        self.rf_minus_rw = wx.TextCtrl(self, wx.ID_ANY, "None", style=wx.TE_PROCESS_ENTER, size=(width_TextCtrl, 30))

        weights_GridBag.Add(wxc_scale, pos=(0, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=5)
        weights_GridBag.Add(self.wxc_scale_TextCtrl, pos=(0, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=5)
        weights_GridBag.Add((30,-1), pos=(0, 2), flag=wx.ALL, border=5)
        weights_GridBag.Add(wxu_scale, pos=(0, 3), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=5)
        weights_GridBag.Add(self.wxu_scale_TextCtrl, pos=(0, 4), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=5)

        weights_GridBag.Add(bonds_rmsd, pos=(1, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=5)
        weights_GridBag.Add(self.bonds_rmsd_TextCtrl, pos=(1, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=5)
        weights_GridBag.Add((30, -1), pos=(1, 2), flag=wx.ALL, border=5)
        weights_GridBag.Add(angle_rmsd, pos=(1, 3), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=5)
        weights_GridBag.Add(self.angle_rmsd_TextCtrl, pos=(1, 4), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=5)
        weights_GridBag.Add((30, -1), pos=(1, 5), flag=wx.ALL, border=5)
        weights_GridBag.Add(Rf_minus_Rw, pos=(1, 6), span=(1, 2), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=5)
        weights_GridBag.Add(self.rf_minus_rw, pos=(1, 8), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=5)

        self.weights.Add(weights_GridBag, 0, wx.ALL, border=5)

        # Other options
        others = wx.StaticBox(self, 1, "Other options", size=(800, -1))
        self.others = wx.StaticBoxSizer(others, wx.VERTICAL)
        others_GridBag = wx.GridBagSizer(vgap=5, hgap=5)
        self.sim_ann = wx.CheckBox(self, wx.ID_ANY, "Simulated annealing")
        self.ordered_solvent = wx.CheckBox(self, wx.ID_ANY, "Ordered solvent")
        self.density_modification = wx.CheckBox(self, wx.ID_ANY, "Density Modification")
        self.map_sharpening = wx.CheckBox(self, wx.ID_ANY, "Map Sharpening")
        border=10
        others_GridBag.Add(self.ordered_solvent, pos=(0, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        others_GridBag.Add(self.map_sharpening,  pos=(0, 2), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        others_GridBag.Add(self.sim_ann, pos=(0, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        others_GridBag.Add(self.density_modification, pos=(0, 3), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)

        self.others.Add(others_GridBag, 0, wx.ALL, border=10)
        ##

        sim_annealing = wx.StaticBox(self, 1, "Simulated Annealing", size=(500, 250))
        self.simulated_annealing = wx.StaticBoxSizer(sim_annealing, wx.VERTICAL)
        sim_annealing_fgs = wx.FlexGridSizer(2, 4, 10, 10)
        start_T = wx.StaticText(self, wx.ID_ANY, "Start T :")
        final_T = wx.StaticText(self, wx.ID_ANY, "Final T :")
        cool_rate = wx.StaticText(self, wx.ID_ANY, "Cooling Rate:")
        mode = wx.StaticText(self, wx.ID_ANY, "Mode :")

        self.start_T = wx.TextCtrl(self, wx.ID_ANY, "5000", style=wx.TE_PROCESS_ENTER, size=(width_TextCtrl, 30))
        self.final_T = wx.TextCtrl(self, wx.ID_ANY, "300", style=wx.TE_PROCESS_ENTER,size=(width_TextCtrl, 30))
        self.cooling_rate = wx.TextCtrl(self, wx.ID_ANY, "100", style=wx.TE_PROCESS_ENTER,
                                   size=(width_TextCtrl, 30))
        self.mode = wx.Choice(self, wx.ID_ANY, choices=["every_macro_cycle", "second_and_before_last", "once", "first",  "first_half"])
        self.mode.SetSelection(1)
        border=10
        sim_annealing_fgs.AddMany([(start_T, 0, wx.ALIGN_CENTER_VERTICAL, border),
                                   (self.start_T, 0, wx.ALIGN_CENTER_VERTICAL, border),
                                   (final_T, 0, wx.ALIGN_CENTER_VERTICAL, border),
                                   (self.final_T, 0, wx.ALIGN_CENTER_VERTICAL, border),
                                   (cool_rate, 0, wx.ALIGN_CENTER_VERTICAL, border),
                                   (self.cooling_rate, 0, wx.ALIGN_CENTER_VERTICAL, border),
                                   (mode, 0, wx.ALIGN_CENTER_VERTICAL, border),
                                   (self.mode, 0, wx.ALIGN_CENTER_VERTICAL, border)])


        self.simulated_annealing.Add(sim_annealing_fgs, 0, wx.ALL, border)

        DM = wx.StaticBox(self, 1, "Density Modification", size=(280, 250))
        self.DM = wx.StaticBoxSizer(DM, wx.VERTICAL)
        DM_GridBag = wx.GridBagSizer(vgap=5, hgap=5)

        combine = wx.StaticText(self, wx.ID_ANY, "Combine :")
        cycles = wx.StaticText(self, wx.ID_ANY, "Cycles :")
        self.combine = wx.Choice(self, wx.ID_ANY, choices=["PERT", "OMIT"])
        self.cycles = wx.TextCtrl(self, wx.ID_ANY, "10", style=wx.TE_PROCESS_ENTER, size=(width_TextCtrl, 30))
        self.combine.SetSelection(0) # Default

        border=5
        DM_GridBag.Add(combine, pos=(0, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        DM_GridBag.Add(self.combine, pos=(0, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        DM_GridBag.Add(cycles, pos=(1, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        DM_GridBag.Add(self.cycles, pos=(1, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        self.DM.Add(DM_GridBag, 0, wx.ALL, border)

        self.option_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.option_sizer.Add(self.simulated_annealing, 1, wx.EXPAND)
        self.option_sizer.AddSpacer(20)
        self.option_sizer.Add(self.DM, 0, wx.EXPAND)

        self.FinalPhenixSizer = wx.BoxSizer(wx.VERTICAL)
        self.FinalPhenixSizer.AddSpacer(15)
        self.FinalPhenixSizer.Add(self.strategy)
        self.FinalPhenixSizer.AddSpacer(15)
        self.FinalPhenixSizer.Add(self.weights)
        self.FinalPhenixSizer.AddSpacer(15)
        self.FinalPhenixSizer.Add(self.others)
        self.FinalPhenixSizer.AddSpacer(15)
        self.FinalPhenixSizer.Add(self.option_sizer)
        #
        # End Phenix
        #
        #
        # Start Refmac
        #

        refmac_refine = wx.StaticBox(self, 1, "Refinement", size=(800, 800))
        self.refmac_refine = wx.StaticBoxSizer(refmac_refine, wx.VERTICAL)
        refmac_refine_gridbag = wx.GridBagSizer(vgap=10, hgap=10)
        size_TxtCtrl = (50, 30)
        size_Static = (175, -1)

        # Type of Refinement and number of cycle
        ref_type = wx.StaticText(self, wx.ID_ANY, "Type of refinement :", size=size_Static)
        self.REFType = wx.Choice(self, wx.ID_ANY, choices=["RESTrained", "UNREstrained", "RIGId"])
        self.REFType.SetSelection(0)  # Default
        cycles = wx.StaticText(self, wx.ID_ANY, "Number of cycles :", size=size_Static)
        self.refmac_cycles = wx.TextCtrl(self, wx.ID_ANY, "20", style=wx.TE_PROCESS_ENTER, size=size_TxtCtrl)

        # Bfactor refinement and Reset
        Bref = wx.StaticText(self, wx.ID_ANY, "Bfactor refinement :", size=size_Static)
        self.Bref = wx.Choice(self, wx.ID_ANY, choices=["OVERall", "ISOTropic"])
        self.Bref.SetSelection(1) #Default
        bfac_set = wx.StaticText(self, wx.ID_ANY, "Reset Bfactor :", size=size_Static)
        self.bfac_set = wx.TextCtrl(self, wx.ID_ANY, "30", style=wx.TE_PROCESS_ENTER, size=size_TxtCtrl)

        # Tls refinement
        tls = wx.StaticText(self, wx.ID_ANY, "TLS Refinement :", size=size_Static)
        self.TLS = wx.Choice(self, wx.ID_ANY, choices=["True", "False"])
        self.TLS.SetSelection(1) # Default
        tls_cycles = wx.StaticText(self, wx.ID_ANY, "Number of TLS_cycles :", size=size_Static)
        self.tls_cycles = wx.TextCtrl(self, wx.ID_ANY, "20", style=wx.TE_PROCESS_ENTER,size=size_TxtCtrl)

        refmac_refine_gridbag.Add(ref_type, pos=(0, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_refine_gridbag.Add(self.REFType, pos=(0, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_refine_gridbag.Add(cycles, pos=(0, 2), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_refine_gridbag.Add(self.refmac_cycles, pos=(0, 3), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)

        refmac_refine_gridbag.Add(Bref, pos=(1, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_refine_gridbag.Add(self.Bref, pos=(1, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_refine_gridbag.Add(bfac_set, pos=(1, 2), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_refine_gridbag.Add(self.bfac_set, pos=(1, 3), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)

        refmac_refine_gridbag.Add(tls, pos=(2, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_refine_gridbag.Add(self.TLS, pos=(2, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_refine_gridbag.Add(tls_cycles, pos=(2, 2), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_refine_gridbag.Add(self.tls_cycles, pos=(2, 3), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)


        self.refmac_refine.Add(refmac_refine_gridbag, 0, wx.ALL, border)


        #
        # Target Weights
        #
        refmac_weights = wx.StaticBox(self, 1, "Target weights", size=(800, 250))
        self.refmac_weights = wx.StaticBoxSizer(refmac_weights, wx.VERTICAL)
        refmac_weights_gridbag = wx.GridBagSizer(vgap=10, hgap=10)

        # Weight and value
        refmac_weight_static = wx.StaticText(self, wx.ID_ANY, "Weight:", size=size_Static)
        self.AUTO = wx.Choice(self, wx.ID_ANY, choices=["AUTO", "MATRIx"])
        self.AUTO.SetSelection(0)  # Default
        weight_term = wx.StaticText(self, wx.ID_ANY, "Weighting term :", size=size_Static)
        self.weighting_term = wx.TextCtrl(self, wx.ID_ANY, "0.2", style=wx.TE_PROCESS_ENTER,size=size_TxtCtrl)
        exp_sig = wx.StaticText(self, wx.ID_ANY, "Experimental Sigmas:", size=size_Static)
        self.NOEX = wx.Choice(self, wx.ID_ANY, choices = ["NOEX", "EXPE"])
        self.NOEX.SetSelection(0)  # Default

        refmac_weights_gridbag.Add(refmac_weight_static, pos=(0, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_weights_gridbag.Add(self.AUTO, pos=(0, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_weights_gridbag.Add((40, -1), pos=(0, 2), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_weights_gridbag.Add(weight_term, pos=(0, 3), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_weights_gridbag.Add(self.weighting_term, pos=(0, 4), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)

        refmac_weights_gridbag.Add(exp_sig, pos=(1, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_weights_gridbag.Add(self.NOEX, pos=(1, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)

        self.refmac_weights.Add(refmac_weights_gridbag, 0, wx.ALL, border)

        #
        # Restraints
        #
        refmac_restraints = wx.StaticBox(self, 1, "Restraints", size=(800, 200))
        self.refmac_restraints = wx.StaticBoxSizer(refmac_restraints, wx.VERTICAL)
        refmac_restraints_gridbag= wx.GridBagSizer(5, 5)

        jelly_body_ref_static = wx.StaticText(self, wx.ID_ANY, "Jelly body refinement :", size=size_Static)
        self.jelly_body_refinement = wx.Choice(self, wx.ID_ANY, choices=["True", "False"])
        self.jelly_body_refinement.SetSelection(1)  # Default
        jelly_sigmas = wx.StaticText(self, wx.ID_ANY, "Jelly body sigma :", size=size_Static)
        self.jbs_TextCtrl = wx.TextCtrl(self, wx.ID_ANY, "0.3", style=wx.TE_PROCESS_ENTER, size=size_TxtCtrl)

        jbar = wx.StaticText(self, wx.ID_ANY, "Additional restraints :", size=size_Static)
        self.jbar_TextCtrl = wx.TextCtrl(self, wx.ID_ANY, "None", style=wx.TE_PROCESS_ENTER, size=size_TxtCtrl)

        extr = wx.StaticText(self, wx.ID_ANY, "External restraints :", size=size_Static)
        self.extr_TextCtrl = wx.TextCtrl(self, wx.ID_ANY, "None", style=wx.TE_PROCESS_ENTER, size=size_TxtCtrl)

        refmac_restraints_gridbag.Add(jelly_body_ref_static, pos=(0, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_restraints_gridbag.Add(self.jelly_body_refinement, pos=(0, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_restraints_gridbag.Add((40, -1), pos=(0, 2), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_restraints_gridbag.Add(jelly_sigmas, pos=(0, 3), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_restraints_gridbag.Add(self.jbs_TextCtrl, pos=(0, 4), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)

        refmac_restraints_gridbag.Add(jbar, pos=(1, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_restraints_gridbag.Add(self.jbar_TextCtrl, pos=(1, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_restraints_gridbag.Add(extr, pos=(2, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL,
                                      border=border)
        refmac_restraints_gridbag.Add(self.extr_TextCtrl, pos=(2, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL,
                                      border=border)

        self.refmac_restraints.Add(refmac_restraints_gridbag, 0, wx.ALL, border)

        # Other options
        refmac_others = wx.StaticBox(self, 1, "Other options", size=(800, 250))
        self.refmac_others = wx.StaticBoxSizer(refmac_others, wx.VERTICAL)
        refmac_others_gridbag = wx.GridBagSizer(vgap=5, hgap=5)

        refmac_combine = wx.StaticText(self, wx.ID_ANY, "Combine :", size=(80, -1))
        refmac_cycles = wx.StaticText(self, wx.ID_ANY, "Cycles :", size=(50, -1))
        self.refmac_combine = wx.Choice(self, wx.ID_ANY, choices=["PERT", "OMIT"])
        self.refmac_DM_cycles = wx.TextCtrl(self, wx.ID_ANY, "10", style=wx.TE_PROCESS_ENTER, size=size_TxtCtrl)
        self.refmac_combine.SetSelection(0)  # Default
        self.refmac_DM = wx.CheckBox(self, wx.ID_ANY, "Density modification")
        self.twinning = wx.CheckBox(self, wx.ID_ANY, "Twinning Refinement")
        self.REFMAC_map_sharpening = wx.CheckBox(self, wx.ID_ANY, "Map Sharpening")

        border = 5
        refmac_others_gridbag.Add(self.refmac_DM, pos=(0, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_others_gridbag.Add(self.twinning, pos=(1, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_others_gridbag.Add(self.REFMAC_map_sharpening, pos=(2, 0), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_others_gridbag.Add((40, -1), pos=(0, 1), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_others_gridbag.Add(refmac_combine, pos=(0, 2), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_others_gridbag.Add(self.refmac_combine, pos=(0, 3), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_others_gridbag.Add((40, -1), pos=(0, 4), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_others_gridbag.Add(refmac_cycles, pos=(0, 5), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        refmac_others_gridbag.Add(self.refmac_DM_cycles, pos=(0, 6), flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=border)
        self.refmac_others.Add(refmac_others_gridbag, 0, wx.ALL, border)


        self.FinalRefmacSizer = wx.BoxSizer(wx.VERTICAL)
        self.FinalRefmacSizer.Add(self.refmac_refine, 0, wx.ALL, border)
        self.FinalRefmacSizer.Add(self.refmac_weights, 0, wx.ALL, border)
        self.FinalRefmacSizer.Add(self.refmac_restraints, 0, wx.ALL, border)
        self.FinalRefmacSizer.Add(self.refmac_others, 0, wx.ALL, border)

        self.FinalSizer =  wx.BoxSizer(wx.VERTICAL)

        self.FinalSizer.Add(soft_sizer, 0, wx.ALIGN_CENTER | wx.ALL, border=15)
        self.FinalSizer.Add(self.FinalPhenixSizer, 0, wx.GROW | wx.ALL, 5)
        self.FinalSizer.Add(self.FinalRefmacSizer, 0, wx.GROW | wx.ALL, 5)
        self.SetSizer(self.FinalSizer)
        self.Layout()
        self.FinalSizer.Hide(self.FinalRefmacSizer)

        self.SoftChoice.Bind(wx.EVT_CHOICE, self.onSoftChanged)

    def onSoftChanged(self, evt):
        choice = self.SoftChoice.GetSelection()
        if choice == 0:
            self.FinalSizer.Hide(self.FinalRefmacSizer)
            self.FinalSizer.Show(self.FinalPhenixSizer)
        if choice == 1:
            self.FinalSizer.Show(self.FinalRefmacSizer)
            self.FinalSizer.Hide(self.FinalPhenixSizer)
        self.Layout()
        try:
            evt.Skip()
        except AttributeError:
            pass


