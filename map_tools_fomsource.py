# -*- coding: utf-8 -*-
"""
Adapted mmtbx.map_tools in order to calculate extrapolated maps according to Genick method in which D-correction shoudl be taken from fmodel with Fextr while fom should be taken from fmodel with Fobs_reference. The two models are supposed to have the same reflections (and order), a check is done but no solution if this is not the case. No special cases.

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
from libtbx.utils import null_out
import mmtbx.map_tools
import libtbx
import random
import boost.python
from six.moves import range
mmtbx_f_model_ext = boost.python.import_ext("mmtbx_f_model_ext")
import mmtbx.masks

class combine(object):
  def __init__(self,
               fmodel_1, #fmodel_1 is the one from which the observations and D will be taken
               fmodel_2, #fmodel_2 ist the one from whih fom will be taken
               map_type_str,
               fo_scale,
               fc_scale,
               #use_shelx_weight, #No Shelx weights
               #shelx_weight_parameter,
               map_calculation_helper_1=None,
               map_calculation_helper_2=None):
    self.mch_1 = map_calculation_helper_1
    self.mch_2 = map_calculation_helper_2
    self.mnm = mmtbx.map_names(map_name_string = map_type_str)
    self.fmodel_1 = fmodel_1
    self.fmodel_2 = fmodel_2
    self.fc_scale = fc_scale
    self.fo_scale = fo_scale
    self.f_obs = None
    self.f_model = None
    if(not self.mnm.ml_map):
      self.f_obs = self.fmodel_1.f_obs().data()*fo_scale #calculation of fo_scale*Fo
    else:
      if(self.mch_1 is None): self.mch_1 = self.fmodel_1.map_calculation_helper()
      if(self.mch_2 is None): self.mch_2 = self.fmodel_2.map_calculation_helper()
      #if(self.fmodel_1.hl_coeffs() is None):
    self.f_obs = self.mch_1.f_obs.data()*fo_scale*self.mch_2.fom #calculation of fo_scale*m*Fo
      #else:
        #cp = fmodel.combine_phases(map_calculation_helper = self.mch_1)
        #self.f_obs = self.mch_1.f_obs.data()*fo_scale*cp.f_obs_phase_and_fom_source()

  def map_coefficients(self, f_model=None):
    def compute(fo,fc,miller_set):
      if(type(fo) == flex.double):
        fo = miller.array(miller_set=miller_set, data=fo).phase_transfer(
          phase_source = self.f_model).data()
      return miller.array(miller_set=miller_set, data=fo+fc)
    if(f_model is None):
      if(not self.mnm.ml_map):
        self.f_model = self.fmodel_1.f_model_scaled_with_k1()
        f_model_data = self.f_model.data()*self.fc_scale #calculation of fc_scale*Fc
      else:
        self.f_model = self.mch_1.f_model
        f_model_data = self.f_model.data()*self.fc_scale*self.mch_1.alpha.data() #calculation of fc_scale*D*Fc
    else:
      self.f_model = f_model
      if(not self.mnm.ml_map):
        f_model_data = self.f_model.data()*self.fc_scale #calculation of fc_scale*Fc
      else:
        f_model_data = self.f_model.data()*self.fc_scale*self.mch_1.alpha.data() #calculation of fc_scale*D*Fc
    # f_model_data may be multiplied by scales like "-1", so it cannot be
    # phase source !
    result = compute(fo=self.f_obs, fc=f_model_data,
      miller_set=self.fmodel_1.f_obs())
    return result

class electron_density_map(object):

  def __init__(self,
               fmodel_1, #fmodel_1 is the one from which the observations and D will be taken
               fmodel_2, #fmodel_2 ist the one from whih fom will be taken
               map_calculation_helper_1 = None,
               map_calculation_helper_2 = None):
    self.fmodel_1 = fmodel_1
    self.fmodel_2 = fmodel_2
    #self.anom_diff = None
    self.mch_1 = map_calculation_helper_1
    self.mch_2 = map_calculation_helper_2
    #if(self.fmodel_1.f_obs().anomalous_flag()):
      #self.anom_diff = self.fmodel_1.f_obs().anomalous_differences()
      #f_model = self.fmodel_1.f_model().as_non_anomalous_array().\
        #merge_equivalents().array()
      #fmodel_match_anom_diff, anom_diff_common = \
        #f_model.common_sets(other =  self.anom_diff)
      #assert anom_diff_common.indices().size()==self.anom_diff.indices().size()
      #self.anom_diff = anom_diff_common.phase_transfer(
        #phase_source = fmodel_match_anom_diff)

  def map_coefficients(self,
                       map_type,
                       acentrics_scale = 2.0,
                       centrics_pre_scale = 1.0,
                       exclude_free_r_reflections=False,
                       fill_missing=False,
                       fill_missing_method="f_model",
                       isotropize=True,
                       sharp=False,
                       pdb_hierarchy=None, # XXX required for map_type=llg
                       #merge_anomalous=None,
                       #use_shelx_weight=False,
                       #shelx_weight_parameter=1.5
                       ):
    map_name_manager = mmtbx.map_names(map_name_string = map_type)
    # No special cases
    # Special case #1: anomalous map
    # Special case #2: anomalous residual map
    # Special case #3: Phaser SAD LLG map
    # Special case #4: Fcalc map
    if(self.mch_1 is None): self.mch_1 = self.fmodel_1.map_calculation_helper()
    if(self.mch_2 is None): self.mch_2 = self.fmodel_2.map_calculation_helper()
    assert list(self.fmodel_1.f_obs().indices())==list(self.fmodel_2.f_obs().indices())==list(self.mch_1.f_obs.indices())==list(self.mch_2.f_obs.indices())
    ffs = mmtbx.map_tools.fo_fc_scales(
      fmodel          = self.fmodel_1,
      map_type_str    = map_type,
      acentrics_scale = acentrics_scale,
      centrics_scale  = centrics_pre_scale)
    fo_scale, fc_scale = ffs.fo_scale, ffs.fc_scale
    coeffs = combine(
      fmodel_1                 = self.fmodel_1,
      fmodel_2                 = self.fmodel_2,
      map_type_str             = map_type,
      fo_scale                 = fo_scale,
      fc_scale                 = fc_scale,
      map_calculation_helper_1 = self.mch_1,
      map_calculation_helper_2 = self.mch_2,
      #use_shelx_weight         = use_shelx_weight,
      #shelx_weight_parameter   = shelx_weight_parameter
      ).map_coefficients()
    r_free_flags = None
    # XXX the default scale array (used for the isotropize option) needs to be
    # calculated and processed now to avoid array size errors
    scale_default = 1. / (self.fmodel_1.k_isotropic()*self.fmodel_1.k_anisotropic())
    scale_array = coeffs.customized_copy(data=scale_default)
    if (exclude_free_r_reflections):
      #if (coeffs.anomalous_flag()):
        #coeffs = coeffs.average_bijvoet_mates()
      r_free_flags = self.fmodel.r_free_flags()
      #if (r_free_flags.anomalous_flag()):
        #r_free_flags = r_free_flags.average_bijvoet_mates()
        #scale_array = scale_array.average_bijvoet_mates()
      coeffs = coeffs.select(~r_free_flags.data())
      scale_array = scale_array.select(~r_free_flags.data())
    scale=None
    if(isotropize):
      if(scale is None):
        if (scale_array.anomalous_flag()) and (not coeffs.anomalous_flag()):
          scale_array = scale_array.average_bijvoet_mates()
        scale = scale_array.data()
      coeffs = coeffs.customized_copy(data = coeffs.data()*scale)
    if(fill_missing):
      #if(coeffs.anomalous_flag()):
        #coeffs = coeffs.average_bijvoet_mates()
      coeffs = mmtbx.map_tools.fill_missing_f_obs(
        coeffs = coeffs,
        fmodel = self.fmodel_1,
        method = fill_missing_method)
    if(sharp):
      ss = 1./flex.pow2(coeffs.d_spacings().data()) / 4.
      from cctbx import adptbx
      b = flex.mean(self.fmodel_1.xray_structure.extract_u_iso_or_u_equiv() *
        adptbx.u_as_b(1))/2
      k_sharp = 1./flex.exp(-ss * b)
      coeffs = coeffs.customized_copy(data = coeffs.data()*k_sharp)
    #if (merge_anomalous) and (coeffs.anomalous_flag()):
      #return coeffs.average_bijvoet_mates()
    return coeffs
