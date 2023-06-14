# -*- coding: utf-8 -*-
"""
Adapted mmtbx.map_tools in order to calculate Fo-Fo difference maps with fom and centric/acentric correction. fobs_in is supposed to be the Miller set containing the Fo-Fo data and needs to be pre-calculated (and q-weighted). The reference fmodel, from which fom will be taken and centric/acentric correction is supposed to have the same reflections (and order), a check is done but no solution if this is not the case. No special cases.

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

Additional license information
-------
cctbx Copyright (c) 2006 - 2018, The Regents of the University of
California, through Lawrence Berkeley National Laboratory (subject to
receipt of any required approvals from the U.S. Dept. of Energy).  All
rights reserved.

see https://github.com/cctbx/cctbx_project/blob/master/LICENSE.txt

-------
"""

from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import miller
from cctbx import maptbx
import mmtbx.map_tools

class combine(object):
  def __init__(self,
               fobs_in, #fobs_in is a Miller set containing the data
               fmodel_2, #fmodel_2 ist the one from whih fom will be taken
               map_type_str,
               fo_scale,
               fc_scale,
               map_calculation_helper_2=None):
    self.mch_2 = map_calculation_helper_2
    self.mnm = mmtbx.map_names(map_name_string = map_type_str)
    self.fobs_in  = fobs_in
    self.fmodel_2 = fmodel_2
    self.fc_scale = fc_scale
    self.fo_scale = fo_scale
    self.f_obs = None
    self.f_model = None
    if(not self.mnm.ml_map):
      self.f_obs = self.fobs_in.data()*fo_scale #calculation of fo_scale*Fo
    else:
      if(self.mch_2 is None): self.mch_2 = self.fmodel_2.map_calculation_helper()
    self.f_obs = self.fobs_in.data()*fo_scale*self.mch_2.fom #calculation of fo_scale*m*Fo

  def map_coefficients(self, f_model=None):
    def compute(fo,miller_set):
      if(type(fo) == flex.double):
        fo = miller.array(miller_set=miller_set, data=fo).phase_transfer(
          phase_source = self.f_model).data()
      return miller.array(miller_set=miller_set, data=fo)
    if(f_model is None): 
      if(not self.mnm.ml_map):
        print("fmodel_2.scale_k1:", self.fmodel_2.scale_k1())
        self.f_model = self.fmodel_2.f_model_scaled_with_k1()
      else:
        self.f_model = self.mch_2.f_model
    else:
      self.f_model = f_model
    result = compute(fo=self.f_obs, miller_set=self.fobs_in)
    return result

class electron_density_map(object):

  def __init__(self,
               fobs_in, #fobs_in is a Miller set containing the data
               fmodel_2, #fmodel_2 ist the one from whih fom will be taken
               map_calculation_helper_2 = None):
    self.fobs_in  = fobs_in
    self.fmodel_2 = fmodel_2
    self.mch_2 = map_calculation_helper_2

  def map_coefficients(self,
                       map_type,
                       acentrics_scale = 2.0,
                       centrics_pre_scale = 1.0,
                       sharp=False,
                       pdb_hierarchy=None): # XXX required for map_type=llg
    map_name_manager = mmtbx.map_names(map_name_string = map_type)
    if(self.mch_2 is None): self.mch_2 = self.fmodel_2.map_calculation_helper()
    assert list(self.fobs_in.indices())==list(self.fmodel_2.f_obs().indices())==list(self.mch_2.f_obs.indices())
    ffs = mmtbx.map_tools.fo_fc_scales(
      fmodel          = self.fmodel_2, #fo_scale and fc_scale can also be calculated with reference fmodel because samen reflections
      map_type_str    = map_type,
      acentrics_scale = acentrics_scale,
      centrics_scale  = centrics_pre_scale)
    fo_scale, fc_scale = ffs.fo_scale, ffs.fc_scale
    coeffs = combine(
      fobs_in                  = self.fobs_in,
      fmodel_2                 = self.fmodel_2,
      map_type_str             = map_type,
      fo_scale                 = fo_scale,
      fc_scale                 = fc_scale,
      map_calculation_helper_2 = self.mch_2,
      ).map_coefficients()
    r_free_flags = None
    if(sharp):
      ss = 1./flex.pow2(coeffs.d_spacings().data()) / 4.
      from cctbx import adptbx
      b = flex.mean(self.fmodel_2.xray_structure.extract_u_iso_or_u_equiv() *
        adptbx.u_as_b(1))/2
      k_sharp = 1./flex.exp(-ss * b)
      coeffs = coeffs.customized_copy(data = coeffs.data()*k_sharp)
    return coeffs
