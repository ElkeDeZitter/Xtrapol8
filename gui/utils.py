import wx
import string

class CharValidator(wx.PyValidator):
    ''' Validates data as it is entered into the text controls. '''

    def __init__(self, flag):
        wx.PyValidator.__init__(self)
        self.flag = flag
        self.Bind(wx.EVT_CHAR, self.OnChar)

    def Clone(self):
        '''Required Validator method'''
        return CharValidator(self.flag)

    def Validate(self, win):
        return True

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True

    def OnChar(self, event):
        keycode = int(event.GetKeyCode())
        if keycode < 256:
            #print keycode
            key = chr(keycode)
            #print key
            if self.flag == 'no-alpha' and key in string.letters:
                return
            if self.flag == 'no-digit' and key in string.digits:
                return
        event.Skip()



