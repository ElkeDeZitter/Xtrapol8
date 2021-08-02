#----------------------------------------------------------------------
# panelIO.py
#
# Created 10/21/2020
#
# Author: Nicolas Coquelle - nicolas.coquelle@esrf.fr
#
#----------------------------------------------------------------------
import os
import sys
import wx
from functools import partial
from Sizers import DefaultHorizontalSizer
from iotbx.file_reader import any_file

class TabIO(wx.Panel):
    """
    This will be the first tab of the Configure Notebook (which holds all X8 inputs)
    """
    # ----------------------------------------------------------------------
    def __init__(self, parent):

        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
        self.defaultDir = os.getcwd()
        self.refMTZ = False
        self.menu_IDS = [wx.ID_ANY, wx.ID_ANY]
        self.dmin = None
        self.dmax = None
        self.files = {'Reference model': [],
                      'Reference mtz': [],
                      'Triggered mtz': [],
                      'Restraints': []}
        self.SetupUI()


    def SetupUI(self):

        # Setting font
        self.defont = wx.Font(11, wx.MODERN, wx.NORMAL, wx.NORMAL, False, 'MS Shell Dlg 2')
        self.SetFont(self.defont)
        # A majority of the sizers in this panel use a special class defined in Sizers.py

        # Dealing with inputs file
        self.list = wx.ListCtrl(self, wx.ID_ANY, style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        self.list.InsertColumn(0, 'File Path', width=700)
        self.list.InsertColumn(1, 'Data type', width=150)

        HorBox = wx.BoxSizer(wx.HORIZONTAL)
        self.AddButton = wx.Button(self, wx.ID_ANY, label="Add file", size=(100,30))
        self.DelButton = wx.Button(self, wx.ID_ANY, label="Remove file", size=(100, 30))
        highRes = wx.StaticText(self, wx.ID_ANY, label="High resolution :")
        self.highRes = wx.TextCtrl(self, wx.ID_ANY, "", size=(80,30))
        lowRes = wx.StaticText(self, wx.ID_ANY, label="Low resolution :")
        self.lowRes = wx.TextCtrl(self, wx.ID_ANY, "", size=(80,30)) 
        HorBox.Add(self.AddButton, 0, wx.ALIGN_CENTER_VERTICAL)
        HorBox.AddSpacer(10)
        HorBox.Add(self.DelButton, 0, wx.ALIGN_CENTER_VERTICAL)
        HorBox.AddSpacer(30)
        HorBox.Add(lowRes, 0, wx.ALIGN_CENTER_VERTICAL)
        HorBox.AddSpacer(5)
        HorBox.Add(self.lowRes, 0, wx.ALIGN_CENTER_VERTICAL)
        HorBox.AddSpacer(30)
        HorBox.Add(highRes, 0, wx.ALIGN_CENTER_VERTICAL)
        HorBox.AddSpacer(5)
        HorBox.Add(self.highRes, 0, wx.ALIGN_CENTER_VERTICAL)

        # All the input widgets are gathered in a Vertical StaticBoxSizer
        Input = wx.StaticBox(self, 1, "Input Files", size=(820, 300))
        sizerIn = wx.StaticBoxSizer(Input, wx.VERTICAL)
        sizerIn.Add(self.list, 1, wx.ALL | wx.EXPAND, 5)
        sizerIn.Add(HorBox, 0, wx.ALL, 5)

        # Now dealing with the output

        # Output directory and its binding
        self.outdir_sizer = DefaultHorizontalSizer(self, "Output Directory :", "Browse", multiline=False)
        self.outdir_sizer.Btn.Bind(wx.EVT_BUTTON,
                                    partial(self.onBrowseDir, text_static=self.outdir_sizer.TextCtrl, multi=False,
                                            key='outdir'))

        # Horizontal Sizer for the output filenames root
        sizer_hor_out = wx.BoxSizer(wx.HORIZONTAL)
        root_static = wx.StaticText(self, wx.ID_ANY, "Output Prefix :", size=(150, -1))
        self.outname = wx.TextCtrl(self, wx.ID_ANY, "", size=(400,30))
        sizer_hor_out.Add(root_static, 1, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 10)
        sizer_hor_out.Add(self.outname, 0, wx.ALL, 10)

        # Adding the phil only option
        # self.FoFoOnly = wx.CheckBox(self, wx.ID_ANY, "Only Generate FoFo")

        # All outputs are gathered in a vertical StaticBoxSizer
        Output = wx.StaticBox(self, 1, "Output Files", size=(820, 800))
        sizerOut = wx.StaticBoxSizer(Output, wx.VERTICAL)
        sizerOut.Add(self.outdir_sizer, 0, wx.ALL, 5)
        sizerOut.Add(sizer_hor_out, 0, wx.ALL, 5)

        # Inputs and Outputs are gathered in a Vertical Sizer.
        FinalSizer = wx.BoxSizer(wx.VERTICAL)
        FinalSizer.Add(sizerIn, 0, wx.ALL | wx.EXPAND, 5)
        FinalSizer.Add(sizerOut, 0, wx.ALL | wx.EXPAND, 5)
        self.SetSizer(FinalSizer)

        self.AddButton.Bind(wx.EVT_BUTTON, self.onAddFile)
        self.DelButton.Bind(wx.EVT_BUTTON, self.onDelFile)
        self.list.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        self.list.Bind(wx.EVT_LIST_ITEM_SELECTED, self.onSelection)

    def onSelection(self, evt):
        self.index = evt.GetIndex()

    def OnRightClick(self, evt):
        self.list_item_clicked = evt.GetItem()
        self.index = evt.GetIndex()

        _, ext = os.path.splitext(self.list.GetItemText(self.index, col=0))#self.list_item_clicked)

        if ext == '.mtz':
            menu = wx.Menu()
            menu.Append(1, "Reference mtz", kind=wx.ITEM_CHECK)
            menu.Append(2, "Triggered mtz", kind=wx.ITEM_CHECK)
            if 'Triggered' in self.list.GetItemText(self.index, col=1):
                menu.Check(2, True)
            else:
                menu.Check(1, True)
            wx.EVT_MENU(menu, 1, self.MenuSelectionCb)
            wx.EVT_MENU(menu, 2, self.MenuSelectionCb)
            self.PopupMenu(menu, evt.GetPoint())
            menu.Destroy()

    def MenuSelectionCb(self, evt):
        ID = evt.GetId()
        print("********************")
        print("Brefore substituion:")
        print("Reference mtz:")
        print self.files['Reference mtz']
        print("********************")
        print("\n\n")

        fn = self.list.GetItem(self.index, 0).GetText()
        if ID == 1:
            self.files["Triggered mtz"].remove(fn)
            self.files["Reference mtz"].append(fn)
            self.list.SetStringItem(self.index, 1, "Reference mtz")
        if ID == 2:
            self.files["Triggered mtz"].append(fn)
            self.files["Reference mtz"].remove(fn)
            self.list.SetStringItem(self.index, 1, "Triggered mtz")
        print("********************")
        print("After substituion:")
        print("Reference mtz:")
        print self.files['Reference mtz']
        print("********************")
        print("\n\n")

    def onBrowseDir(self, evt, text_static, multi, key):
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
        
        if multi:
            text_static.WriteText(path+'\n')
        else:
            text_static.WriteText(path)

    def onAddFile(self, evt):
        """
        File selection.
        :param evt: wx.EVT_BUTTON (guessed)
        :param text_static: str, modified TextStatic widget
        :param multi: bool, is text_static multiline (This could be determined by text_static attribute)
        :param key: str, used to set the value of this key in self.inputs dict
        :param ext: Could be used to filter file extensions (Not a the moment
        :return: None
        """
        flags = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_MULTIPLE
        dlg = wx.FileDialog(self,
                            'Open File',
                            defaultDir=self.defaultDir,
                            wildcard="All files (*.*)|*.*|Model files (*.pdb)|*.pdb|Reflections files (*.mtz)|*.mtz|Restraints files (*.cif)|*.cif",
                            style=flags)

        if dlg.ShowModal() == wx.ID_CANCEL:
            return

        for path in dlg.GetPaths():
            self.defaultDir = os.path.dirname(path)
            _, ext = os.path.splitext(path)

            if ext not in ['.pdb', '.cif', '.mtz']:
                print("Sorry: this file type is not (yet?) accepted")
                return

            index = self.list.InsertStringItem(sys.maxint, path)
            if ext == '.pdb':
                self.list.SetStringItem(index, 1, 'Reference model')
            elif ext == '.cif':
                try:
                    cif = any_file(path, force_type="pdb")
                    self.list.SetStringItem(index, 1, 'Reference model')
                    self.files['Reference model'].append(path)
                except ValueError:
                    self.list.SetStringItem(index, 1, 'Restraints')
                    self.files["Restraints"].append(path)
            else:
                if self.refMTZ:
                    self.list.SetStringItem(index, 1, 'Triggered mtz')
                    self.files["Triggered mtz"].append(path)
                else:
                    self.list.SetStringItem(index, 1, 'Reference mtz')
                    self.refMTZ = True
                    self.files["Reference mtz"].append(path)
                self.extract_dmin_dmax(path)

    def extract_dmin_dmax(self, path):
        hkl_in = any_file(path, force_type="hkl", raise_sorry_if_errors=True)
        miller_arrays = hkl_in.file_object.as_miller_arrays()
        self._dmax, self._dmin = miller_arrays[0].resolution_range()

        if self.dmax is None:
            self.dmax = self._dmax
        else:
            if self._dmax < self.dmax:
                self.dmax = self._dmax
        self.lowRes.SetValue("%4.2f" % self._dmax)

        if self.dmin is None:
            self.dmin = self._dmin
        else:
            if self._dmin > self.dmin:
                self.dmin = self._dmin
        self.highRes.SetValue("%4.2f" % self._dmin)
        #print self._dmin, self._dmax

    def onDelFile(self, evt):
        if self.list.ItemCount > 0 and hasattr(self,"index"):
            fn = self.list.GetItem(self.index, 0).GetText()
            file_type = self.list.GetItem(self.index, 1).GetText()
            self.files[file_type].remove(fn)
            print self.files[file_type]
            self.list.DeleteItem(self.index)
            evt.Skip()







