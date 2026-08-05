# -*- coding: utf-8 -*-
"""
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

import numpy as np


class Pymol_movie:
    """
    Aim is to merge Pymol visualization
    """

    def __init__(
        self,
        occupancies,
        pdblst=None,
        ccp4_maps=None,
        resids_lst=None,
        model_label="",
        ccp4_map_label="",
    ):
        self.resids_lst = resids_lst
        self.model_label = model_label
        self.ccp4_map_label = ccp4_map_label

        if pdblst != None:
            # Sort the pdb files and alphas in order to have alpha from small to large, this is important for fitting
            # This might not work in pyhton3
            alphas = [round(1 / x, 3) for x in occupancies]
            if ccp4_maps != None:
                zipped = zip(alphas, pdblst, ccp4_maps)
                zipped_sorted = sorted(zipped, key=lambda x: x[0])
                alphas, pdblst, ccp4_maps = zip(*zipped_sorted)
                alphas = list(alphas)
                self.pdblst = list(pdblst)
                self.ccp4_maps = list(ccp4_maps)
            else:
                zipped = zip(alphas, pdblst)
                zipped_sorted = sorted(zipped, key=lambda x: x[0])
                alphas, pdblst = zip(*zipped_sorted)
                alphas = list(alphas)
                self.pdblst = list(pdblst)
                self.ccp4_maps = None
        else:
            self.pdblst = None
            self.ccp4_maps = None
            print("No models selected for pymol visualization.")

    def get_residlist(self):
        if self.resids_lst != None:
            with open(self.resids_lst) as rs_lst:
                residlst = rs_lst.read().split("\n")
                residlst = [lne for lne in residlst if len(lne) > 0]
            if len(residlst) == 0:
                self.resids_lst = None
                return
            residlst = np.array(residlst)
            residlst_unique = list(set(residlst))
            residlst_unique.remove(
                [lne for lne in residlst_unique if "Resn Resv Chain Alt" in lne][0]
            )
            self.residlst = np.array(residlst_unique)

    def get_syntax(self):
        self.resids = False
        self.resid_lines = []
        self.waters = False
        self.water_lines = []
        for resid in self.residlst:
            try:
                resn, resv, chain, alt = resid.split("\n")[0].split()[:4]
                if resn == "HOH":
                    self.waters = True
                    self.water_lines.append(
                        "(chain %s and resid %s and altloc %s)" % (chain, resv, alt)
                    )
                else:
                    self.resids = True
                    self.resid_lines.append(
                        "(chain %s and resid %s and altloc %s)" % (chain, resv, alt)
                    )
            except ValueError:
                resn, resv, chain = resid.split("\n")[0].split()[:3]
                if resn == "HOH":
                    self.waters = True
                    self.water_lines.append("(chain %s and resid %s)" % (chain, resv))
                else:
                    self.resids = True
                    self.resid_lines.append("(chain %s and resid %s)" % (chain, resv))

    def load_models(self, outscript):
        for i, pdb in enumerate(self.pdblst):
            outscript.write(
                'cmd.load("%s","models_%s",state=%d)\n' % (pdb, self.model_label, i + 1)
            )

    def show_default(self, outscript):
        outscript.write("cmd.hide()\n")
        outscript.write('cmd.show("lines")\n')
        outscript.write('cmd.show("nonbonded")\n')

    def show_sticks_and_spheres(self, outscript):
        if self.resids_lst == None:
            print("Showing all residues in pymol movie.")
            outscript.write('cmd.show("sticks")\n')
            outscript.write('cmd.show("nb_spheres")\n')
        else:
            if self.resids:
                outscript.write(
                    'cmd.show("sticks", "%s")\n' % (" or ".join(self.resid_lines))
                )
            if self.waters:
                outscript.write(
                    'cmd.show("nb_spheres", "%s")\n' % (" or ".join(self.water_lines))
                )

    def load_maps(self, outscript):
        for i, ccp4 in enumerate(self.ccp4_maps):
            outscript.write(
                'cmd.load("%s","maps_%s",%d,"ccp4")\n'
                % (ccp4, self.ccp4_map_label, i + 1)
            )

    def show_mesh(self, outscript):
        if self.resids_lst == None:
            print("No residues selected for isomesh.")
        else:
            if self.waters and self.resids:
                selection = " or ".join(
                    [" or ".join(self.resid_lines), " or ".join(self.water_lines)]
                )
            elif self.resids:
                selection = " or ".join(self.resid_lines)
            else:
                selection = " or ".join(self.water_lines)
            outscript.write(
                'cmd.isomesh("%s", "maps_%s", 1.0, "%s", carve=1.6)\n'
                % (self.ccp4_map_label, self.ccp4_map_label, selection)
            )

    def additional_settings(self, outscript):
        outscript.write('cmd.set("movie_fps", 5)\n')
        outscript.write('cmd.show("mesh", "F*")\n')
        outscript.write('cmd.show("mesh", "qF*")\n')
        outscript.write('cmd.set("mesh_width", 0.3)\n')

    def write_pymol_script(self, outfile="pymol_movie.py"):
        self.get_residlist()
        if self.resids_lst != None:
            self.get_syntax()
        o = open(outfile, "a")
        self.load_models(o)
        if self.ccp4_maps != None:
            self.load_maps(o)
            self.show_mesh(o)
        o.close()

    def write_pymol_appearance(self, outfile="pymol_movie.py"):
        self.get_residlist()
        if self.resids_lst != None:
            self.get_syntax()
        o = open(outfile, "a")
        self.show_default(o)
        self.show_sticks_and_spheres(o)
        self.additional_settings(o)
        o.close()
