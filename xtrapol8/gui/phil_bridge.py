"""Bridge between Qt controls and Xtrapol8's authoritative master PHIL.

Scientific imports are lazy so the Qt interface can still start without Phenix.
"""
from pathlib import Path
from PyQt5.QtWidgets import QCheckBox, QComboBox, QLineEdit, QTableWidgetItem


def master():
    try:
        from xtrapol8.master import master_phil
        from iotbx import phil
    except ImportError as exc:
        raise RuntimeError('PHIL support requires Phenix/CCTBX (iotbx.phil). Run the GUI using phenix.python or an environment containing iotbx.') from exc
    return master_phil, phil


def get(obj, path):
    for part in path.split('.'):
        obj = getattr(obj, part)
    return obj


def put(obj, path, value):
    parts = path.split('.')
    for part in parts[:-1]:
        obj = getattr(obj, part)
    setattr(obj, parts[-1], value)


# Qt control name -> exact path in xtrapol8/master.py.
FIELDS = {
    'LowTextCtrl':'occupancies.low_occ', 'HighTextCtrl':'occupancies.high_occ',
    'StepsTextCtrl':'occupancies.steps', 'ListTextCtrl':'occupancies.list_occ',
    'FoChoice':'f_and_maps.fofo_type', 'kscale':'f_and_maps.kweight_scale',
    'peak_detection_thresholdTextCtrl':'map_explorer.peak_detection_threshold',
    'peak_integration_floorTextCtrl':'map_explorer.peak_integration_floor',
    'RadiusTextCtrl':'map_explorer.radius', 'ZscoreTextCtrl':'map_explorer.z_score',
    'OccEstimation':'map_explorer.occupancy_estimation',
    'ScalingChoice':'scaling.b_scaling', 'ScalingHighRes':'scaling.high_resolution',
    'ScalingLowRes':'scaling.low_resolution',
    'RunRef':'refinement.run_refinement',
    'NCyclesReal_TextCtrl':'refinement.phenix_keywords.real_space_refine.cycles',
    'NCyclesReciprocal_TextCtrl':'refinement.phenix_keywords.main.cycles',
    'wxc_scale_TextCtrl':'refinement.phenix_keywords.target_weights.wxc_scale',
    'wxu_scale_TextCtrl':'refinement.phenix_keywords.target_weights.wxu_scale',
    'bonds_rmsd_TextCtrl':'refinement.phenix_keywords.target_weights.weight_selection_criteria.bonds_rmsd',
    'angle_rmsd_TextCtrl':'refinement.phenix_keywords.target_weights.weight_selection_criteria.angles_rmsd',
    'rf_minus_rw':'refinement.phenix_keywords.target_weights.weight_selection_criteria.r_free_minus_r_work',
    'sim_ann':'refinement.phenix_keywords.main.simulated_annealing',
    'ordered_solvent':'refinement.phenix_keywords.main.ordered_solvent',
    'density_modification':'refinement.phenix_keywords.density_modification.density_modification',
    'map_sharpening':'refinement.phenix_keywords.map_sharpening.map_sharpening',
    'start_T':'refinement.phenix_keywords.simulated_annealing.start_temperature',
    'final_T':'refinement.phenix_keywords.simulated_annealing.final_temperature',
    'cooling_rate':'refinement.phenix_keywords.simulated_annealing.cool_rate',
    'mode':'refinement.phenix_keywords.simulated_annealing.mode',
    'combine':'refinement.phenix_keywords.density_modification.combine',
    'cycles':'refinement.phenix_keywords.density_modification.cycles',
    'REFType':'refinement.refmac_keywords.refine.type',
    'refmac_cycles':'refinement.refmac_keywords.refine.cycles',
    'Bref':'refinement.refmac_keywords.refine.Brefinement',
    'bfac_set':'refinement.refmac_keywords.refine.bfac_set',
    'TLS':'refinement.refmac_keywords.refine.TLS',
    'tls_cycles':'refinement.refmac_keywords.refine.TLS_cycles',
    'AUTO':'refinement.refmac_keywords.target_weights.weight',
    'weighting_term':'refinement.refmac_keywords.target_weights.weighting_term',
    'NOEX':'refinement.refmac_keywords.target_weights.experimental_sigmas',
    'jelly_body_refinement':'refinement.refmac_keywords.restraints.jelly_body_refinement',
    'jbs_TextCtrl':'refinement.refmac_keywords.restraints.jelly_body_sigma',
    'jbar_TextCtrl':'refinement.refmac_keywords.restraints.jelly_body_additional_restraints',
    'extr_TextCtrl':'refinement.refmac_keywords.restraints.external_restraints',
    'refmac_combine':'refinement.refmac_keywords.density_modification.combine',
    'refmac_DM_cycles':'refinement.refmac_keywords.density_modification.cycles',
    'refmac_DM':'refinement.refmac_keywords.density_modification.density_modification',
    'twinning':'refinement.refmac_keywords.refine.twinning',
    'REFMAC_map_sharpening':'refinement.refmac_keywords.map_sharpening.map_sharpening',
}
STRATEGIES = ['individual_sites','individual_sites_real_space','rigid_body','individual_adp','group_adp','tls','occupancies','group_anomalous']
FACTORS = ['qfextr','fextr','kfextr','qfgenick','fgenick','kfgenick','qfextr_calc','fextr_calc','kfextr_calc']
CHOICE_ALIASES = {
    'REFType': {'REST':'RESTrained','RIGID':'RIGId','TLS':'UNREstrained'},
    'Bref': {'OVERall':'OVERall','ISOTropic':'ISOTropic'},
    'mode': {'every_macro_cycle':'every_macro_cycle','first':'first','once':'once'},
}

# Fields whose master PHIL default is None, despite their declared numeric type.
# Without explicit conversion, a populated QLineEdit would be saved as str.
NULLABLE_FLOAT_FIELDS = {
    'RadiusTextCtrl',  # map_explorer.radius
    'rf_minus_rw',     # refinement.phenix_keywords.target_weights...
    'ScalingHighRes', 'ScalingLowRes',  # scaling resolution defaults are None
}


def _from_widget(widget, previous, name, path=None):
    """Convert a Qt value using the authoritative PHIL schema, not defaults.

    PHIL defaults may be None, so checking isinstance(previous, float) is
    insufficient. This follows the wx GUI's parse -> master.fetch -> extract
    approach, but applies it to each field for precise error reporting.
    """
    if isinstance(widget, QCheckBox):
        raw = 'True' if widget.isChecked() else 'False'
    elif isinstance(widget, QComboBox):
        value = widget.currentData()
        if value is None:
            value = widget.currentText()
        raw = str(CHOICE_ALIASES.get(name, {}).get(str(value), value))
    else:
        raw = widget.text().strip() or 'None'
        if name == 'ListTextCtrl' and raw != 'None':
            raw = ' '.join(raw.replace(',', ' ').split())
    if path is None:
        return raw
    # 'AUTO' is a literal PHIL choice, but PHIL's special Auto token can
    # be extracted as None when parsed alone. Preserve the declared choice.
    if name == 'AUTO' and raw == 'AUTO':
        return 'AUTO'
    master_phil, phil = master()
    try:
        # PHIL performs type conversion (float/int/bool/floats/path/choice)
        # and validation, even when the master default is None.
        source = phil.parse(f'{path} = {raw}\n')
        parsed = master_phil.fetch(source=source).extract()
        return get(parsed, path)
    except Exception as exc:
        raise ValueError(f'Invalid GUI field {name!r} ({path}): {raw!r}: {exc}') from exc


def _set_widget(widget, value, name):
    if isinstance(widget, QCheckBox):
        widget.setChecked(bool(value))
    elif isinstance(widget, QComboBox):
        target = str(value)
        index = widget.findData(target)
        if index < 0:
            index = widget.findText(target)
        if index < 0:
            for i in range(widget.count()):
                if CHOICE_ALIASES.get(name, {}).get(widget.itemText(i)) == target:
                    index = i; break
        if index >= 0:
            widget.setCurrentIndex(index)
    else:
        if isinstance(value, (list, tuple)):
            value = ' '.join(map(str, value))
        widget.setText('' if value is None else str(value))


def _legacy_refinement_only(filename, phil):
    """Recognize legacy-only selection without mistaking master defaults for input."""
    source = phil.parse(file_name=str(filename), process_includes=True)
    # CCTBX versions differ: all_definitions() may yield object_locator
    # wrappers rather than the PHIL definitions themselves.
    paths = set()
    for locator in source.all_definitions():
        definition = getattr(locator, 'object', locator)
        full_path = getattr(definition, 'full_path', None)
        if callable(full_path):
            paths.add(full_path())
        else:
            # Some libtbx builds expose the path on the locator.
            path = getattr(locator, 'path', None)
            if isinstance(path, str):
                paths.add(path)
    if not paths:
        raise RuntimeError('Cannot inspect explicit PHIL definitions for legacy refinement selectors')
    return ('refinement.use_refmac_instead_of_phenix' in paths
            and 'refinement.reciprocal_space' not in paths
            and 'refinement.real_space' not in paths)


def read_into(window, filename):
    master_phil, phil = master()
    scope = master_phil.fetch(source=phil.parse(file_name=str(filename), process_includes=True))
    params = scope.extract()
    window._phil_params = params
    io = window.tabIO
    io.list.setRowCount(0)
    for key, path in [('Reference model','input.reference_pdb'),('Reference mtz','input.reference_mtz'),('Triggered mtz','input.triggered_mtz')]:
        value = get(params,path)
        if value:
            _add_row(io, str(value), key)
    for value in get(params,'input.additional_files') or []:
        if value:
            _add_row(io, str(value), 'Restraints')
    for widget, path in [(io.lowRes,'input.low_resolution'),(io.highRes,'input.high_resolution'),(io.outdir,'output.outdir'),(io.outname,'output.outname')]:
        _set_widget(widget,get(params,path),'')
    ex = window.tabExtrapolation
    mode = 0 if get(params,'output.generate_fofo_only') else 1 if get(params,'f_and_maps.fast_and_furious') else 2
    ex.mode_radios[mode].setChecked(True)
    ex.set_mode(mode)
    for panel in (ex, window.tabRefinement):
        for name, path in FIELDS.items():
            if name in panel.controls:
                _set_widget(panel.controls[name], get(params,path), name)
    for name in FACTORS:
        ex.controls[name].setChecked(name in (get(params,'f_and_maps.f_extrapolated_and_maps') or []))
    for name in STRATEGIES:
        window.tabRefinement.controls[name].setChecked(name in (get(params,'refinement.phenix_keywords.refine.strategy') or []))
    negative = str(get(params,'f_and_maps.negative_and_missing'))
    # Preserve the exact combined choice, even if the GUI only exposes a simplified pair.
    ex.negChoice.setCurrentText(negative.split('_and_')[0].replace('_no_fill',''))
    ex.missChoice.setCurrentText('no_fill' if negative.endswith('no_fill') else 'fill')
    ref = window.tabRefinement
    # Legacy PHIL files used a single boolean for both refinement programs.
    # Only apply it when neither modern selector is explicitly present in the file.
    reciprocal = get(params, 'refinement.reciprocal_space')
    real = get(params, 'refinement.real_space')
    if _legacy_refinement_only(filename, phil):
        legacy = bool(get(params, 'refinement.use_refmac_instead_of_phenix'))
        reciprocal = 'refmac5' if legacy else 'phenix.refine'
        real = 'coot' if legacy else 'phenix.real_space_refine'
    ref.SoftChoiceReci.setCurrentText(str(reciprocal))
    ref.SoftChoiceReal.setCurrentText(str(real))
    ref.update_reciprocal_sections()
    return params


def _add_row(io, path, category):
    row = io.list.rowCount(); io.list.insertRow(row)
    io.list.setItem(row,0,QTableWidgetItem(path))
    io.list.setItem(row,1,QTableWidgetItem(category))


def collect(window):
    master_phil, _ = master()
    params = getattr(window,'_phil_params',None)
    if params is None:
        params = master_phil.extract()
    io = window.tabIO
    grouped = {'Reference model':[],'Reference mtz':[],'Triggered mtz':[],'Restraints':[]}
    for row in range(io.list.rowCount()):
        file_item, type_item = io.list.item(row,0), io.list.item(row,1)
        if file_item and type_item and type_item.text() in grouped:
            grouped[type_item.text()].append(file_item.text())
    for key,path in [('Reference model','input.reference_pdb'),('Reference mtz','input.reference_mtz'),('Triggered mtz','input.triggered_mtz')]:
        put(params,path,grouped[key][0] if grouped[key] else None)
    put(params,'input.additional_files',grouped['Restraints'])
    for widget,path in [(io.lowRes,'input.low_resolution'),(io.highRes,'input.high_resolution'),(io.outdir,'output.outdir'),(io.outname,'output.outname')]:
        put(params,path,_from_widget(widget,get(params,path),'',path))
    for panel in (window.tabExtrapolation,window.tabRefinement):
        for name,path in FIELDS.items():
            if name in panel.controls:
                put(params,path,_from_widget(panel.controls[name],get(params,path),name,path))
    ex = window.tabExtrapolation
    mode = ex.mode_buttons.checkedId()
    # Coot is launched explicitly from Results, never by the scientific subprocess.
    # This also overrides open_coot=True in an imported PHIL.
    put(params,'output.open_coot',False)
    put(params,'output.generate_fofo_only',mode == 0)
    put(params,'f_and_maps.fast_and_furious',mode == 1)
    put(params,'f_and_maps.f_extrapolated_and_maps',[name for name in FACTORS if ex.controls[name].isChecked()])
    negative = ex.negChoice.currentText()
    missing = ex.missChoice.currentText()
    # Master PHIL has one combined choice, not two independent parameters.
    neg_mapping = {'truncate':'truncate','fref':'fref','fcalc':'fcalc','keep':'keep','reject':'reject','zero':'zero'}
    if negative in neg_mapping:
        put(params,'f_and_maps.negative_and_missing',f'{neg_mapping[negative]}_and_fill' if missing == 'fill' else f'{neg_mapping[negative]}_no_fill')
    ref = window.tabRefinement
    put(params,'refinement.reciprocal_space',ref.SoftChoiceReci.currentText())
    put(params,'refinement.real_space',ref.SoftChoiceReal.currentText())
    # The deprecated combined boolean cannot represent independent selections.
    # Leave it False: modern reciprocal_space / real_space are authoritative.
    put(params,'refinement.use_refmac_instead_of_phenix',False)
    put(params,'refinement.phenix_keywords.refine.strategy',[name for name in STRATEGIES if ref.controls[name].isChecked()])
    return params


def save_from(window, filename):
    master_phil, _ = master()
    params = collect(window)
    try:
        formatted = master_phil.format(python_object=params)
    except (TypeError, ValueError) as exc:
        raise ValueError('PHIL serialization failed after schema-based widget conversion: '
                         + str(exc) + '. Please send the full traceback.') from exc
    with open(filename,'w',encoding='utf-8') as stream:
        formatted.show(out=stream)
    window._phil_params = params
    return Path(filename)
