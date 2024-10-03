# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#  .
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#  .
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.

__authors__ = ["Elke De Zitter", "Nicolas Coquelle", "Paula Oeser", "Thomas Barends", "Jacques Philippe Colletier"]
__license__ = "MIT"
__copyright__ = "Institut de Biologie Structurale - group DYNAMOP - team SNaX"
__date__ = "01/12/2022"
__citation__ = "doi.org/10.1038/s42003-022-03575-7"

major = 1
minor = 2
micro = 7
tag = "Dev" #can be None or "Dev"
# tag = None #can be None or "Dev"

VERSION = ".".join(map(lambda x: str(x) , [major, minor, micro]))
if tag != None:
    VERSION += "-{:s}".format(tag)

if __name__=="__main__":
    import sys
    if len(sys.argv) < 2 :
        print(VERSION)
    elif sys.argv[1] == 'update':
        f = "CITATION.cff"
        print("Updating version number in {:s}".format(f))
        with open(f) as cff:
            c =  cff.readlines()
        o = open(f,"w")
        for line in c:
            if line.startswith("version"):
                new_line = 'version: "{:s}"\n'.format(VERSION)
                o.write(new_line)
            else:
                o.write(line)
        o.close()
    else:
        print("Argument not recognised")
