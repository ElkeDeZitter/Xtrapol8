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

import mmtbx.map_tools
from cctbx import adptbx, miller
from cctbx.array_family import flex


class combine(object):
    def __init__(
        self,
        fmodel_1,  # fmodel_1 is the one from which the observations and D will be taken
        fmodel_2,  # fmodel_2 ist the one from whih fom will be taken
        map_type_str,
        fo_scale,
        fc_scale,
        map_calculation_helper_1=None,
        map_calculation_helper_2=None,
    ):
        self.mch_1 = map_calculation_helper_1
        self.mch_2 = map_calculation_helper_2
        self.mnm = mmtbx.map_names(map_name_string=map_type_str)
        self.fmodel_1 = fmodel_1
        self.fmodel_2 = fmodel_2
        self.fc_scale = fc_scale
        self.fo_scale = fo_scale
        self.f_obs = None
        self.f_model = None
        if not self.mnm.ml_map:
            self.f_obs = (
                self.fmodel_1.f_obs().data() * fo_scale
            )  # calculation of fo_scale*Fo
        else:
            if self.mch_1 is None:
                self.mch_1 = self.fmodel_1.map_calculation_helper()
            if self.mch_2 is None:
                self.mch_2 = self.fmodel_2.map_calculation_helper()
        self.f_obs = (
            self.mch_1.f_obs.data() * fo_scale * self.mch_2.fom
        )  # calculation of fo_scale*m*Fo

    def map_coefficients(self, f_model=None):
        def compute(fo, fc, miller_set):
            if type(fo) == flex.double:
                fo = (
                    miller.array(miller_set=miller_set, data=fo)
                    .phase_transfer(phase_source=self.f_model)
                    .data()
                )
            return miller.array(miller_set=miller_set, data=fo + fc)

        if f_model is None:
            if not self.mnm.ml_map:
                self.f_model = self.fmodel_1.f_model_scaled_with_k1()
                f_model_data = (
                    self.f_model.data() * self.fc_scale
                )  # calculation of fc_scale*Fc
            else:
                self.f_model = self.mch_1.f_model
                f_model_data = (
                    self.f_model.data() * self.fc_scale * self.mch_1.alpha.data()
                )  # calculation of fc_scale*D*Fc
        else:
            self.f_model = f_model
            if not self.mnm.ml_map:
                f_model_data = (
                    self.f_model.data() * self.fc_scale
                )  # calculation of fc_scale*Fc
            else:
                f_model_data = (
                    self.f_model.data() * self.fc_scale * self.mch_1.alpha.data()
                )  # calculation of fc_scale*D*Fc
        result = compute(
            fo=self.f_obs, fc=f_model_data, miller_set=self.fmodel_1.f_obs()
        )
        return result


class electron_density_map(object):
    def __init__(
        self,
        fmodel_1,  # fmodel_1 is the one from which the observations and D will be taken
        fmodel_2,  # fmodel_2 ist the one from whih fom will be taken
        map_calculation_helper_1=None,
        map_calculation_helper_2=None,
    ):
        self.fmodel_1 = fmodel_1
        self.fmodel_2 = fmodel_2
        self.mch_1 = map_calculation_helper_1
        self.mch_2 = map_calculation_helper_2

    def map_coefficients(
        self,
        map_type,
        acentrics_scale=2.0,
        centrics_pre_scale=1.0,
        exclude_free_r_reflections=False,
        fill_missing=False,
        fill_missing_method="f_model",
        isotropize=True,
        sharp=False,
        pdb_hierarchy=None,  # XXX required for map_type=llg
    ):
        map_name_manager = mmtbx.map_names(map_name_string=map_type)
        if self.mch_1 is None:
            self.mch_1 = self.fmodel_1.map_calculation_helper()
        if self.mch_2 is None:
            self.mch_2 = self.fmodel_2.map_calculation_helper()
        assert (
            list(self.fmodel_1.f_obs().indices())
            == list(self.fmodel_2.f_obs().indices())
            == list(self.mch_1.f_obs.indices())
            == list(self.mch_2.f_obs.indices())
        )
        ffs = mmtbx.map_tools.fo_fc_scales(
            fmodel=self.fmodel_1,
            map_type_str=map_type,
            acentrics_scale=acentrics_scale,
            centrics_scale=centrics_pre_scale,
        )
        fo_scale, fc_scale = ffs.fo_scale, ffs.fc_scale
        coeffs = combine(
            fmodel_1=self.fmodel_1,
            fmodel_2=self.fmodel_2,
            map_type_str=map_type,
            fo_scale=fo_scale,
            fc_scale=fc_scale,
            map_calculation_helper_1=self.mch_1,
            map_calculation_helper_2=self.mch_2,
        ).map_coefficients()
        r_free_flags = None
        scale_default = 1.0 / (
            self.fmodel_1.k_isotropic() * self.fmodel_1.k_anisotropic()
        )
        scale_array = coeffs.customized_copy(data=scale_default)
        if exclude_free_r_reflections:
            r_free_flags = self.fmodel.r_free_flags()
            coeffs = coeffs.select(~r_free_flags.data())
            scale_array = scale_array.select(~r_free_flags.data())
        scale = None
        if isotropize:
            if scale is None:
                if (scale_array.anomalous_flag()) and (not coeffs.anomalous_flag()):
                    scale_array = scale_array.average_bijvoet_mates()
                scale = scale_array.data()
            coeffs = coeffs.customized_copy(data=coeffs.data() * scale)
        if fill_missing:
            coeffs = mmtbx.map_tools.fill_missing_f_obs(
                coeffs=coeffs, fmodel=self.fmodel_1, method=fill_missing_method
            )
        if sharp:
            ss = 1.0 / flex.pow2(coeffs.d_spacings().data()) / 4.0
            b = (
                flex.mean(
                    self.fmodel_1.xray_structure.extract_u_iso_or_u_equiv()
                    * adptbx.u_as_b(1)
                )
                / 2
            )
            k_sharp = 1.0 / flex.exp(-ss * b)
            coeffs = coeffs.customized_copy(data=coeffs.data() * k_sharp)
        return coeffs
