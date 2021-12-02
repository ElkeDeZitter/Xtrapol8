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
import wx

class DefaultHorizontalSizer(wx.BoxSizer):


    def __init__(self, parent, label, label_button, multiline=False,padding=10):

        wx.BoxSizer.__init__(self,wx.HORIZONTAL)
        width_TextCtrl = 400
        height = 30
        #defont = wx.Font(11, wx.MODERN, wx.NORMAL, wx.NORMAL, False, 'MS Shell Dlg 2')

        if multiline:
            self.StaticTxt = wx.StaticText(parent, wx.ID_ANY, label, size=(150, 60))
            self.TextCtrl =  wx.TextCtrl(parent, wx.ID_ANY, "",style = wx.TE_READONLY | wx.TE_MULTILINE,size=(width_TextCtrl,60))
        else:
            self.TextCtrl = wx.TextCtrl(parent, wx.ID_ANY, "", style=wx.TE_READONLY,size=(width_TextCtrl,height))
            self.StaticTxt = wx.StaticText(parent, wx.ID_ANY, label, size=(150, -1))

        self.Btn =  wx.Button(parent, wx.ID_ANY, label_button)

        #self.StaticTxt.SetFont(defont)
        #self.TextCtrl.SetFont(defont)
        #self.Btn.SetFont(defont)
        self.Add(self.StaticTxt, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, padding)
        self.Add(self.TextCtrl,  0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, padding)
        self.Add(self.Btn,       0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, padding)
