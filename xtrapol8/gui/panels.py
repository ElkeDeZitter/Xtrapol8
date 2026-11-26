"""Qt configuration panels mirroring the wx GUI groups; no scientific imports."""
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QWidget,QScrollArea,QVBoxLayout,QHBoxLayout,QGridLayout,QGroupBox,
    QLabel,QLineEdit,QComboBox,QCheckBox,QPushButton,QFormLayout,QFrame,QRadioButton,QButtonGroup,QSizePolicy)

class Panel(QScrollArea):
    def __init__(self,parent=None):
        super().__init__(parent); self.setWidgetResizable(True)
        self.body=QWidget(); self.setWidget(self.body)
        self.root=QVBoxLayout(self.body); self.root.setSpacing(14); self.root.setContentsMargins(16,16,16,16)
        self.body.setStyleSheet("QGroupBox { font-size: 13px; font-weight: 600; margin-top: 12px; padding-top: 12px; } QGroupBox QLabel, QGroupBox QRadioButton, QGroupBox QCheckBox, QGroupBox QPushButton, QGroupBox QLineEdit, QGroupBox QComboBox { font-size: 13px; font-weight: 400; }")
        self.controls={}
    def group(self,title):
        g=QGroupBox(title); layout=QVBoxLayout(g); layout.setSpacing(9); layout.setContentsMargins(14,18,14,12); self.root.addWidget(g); return layout
    def row(self,layout,*items):
        r=QHBoxLayout()
        for x in items:
            if isinstance(x,str): x=QLabel(x)
            r.addWidget(x)
        r.addStretch(); layout.addLayout(r)
    def text(self,name,value='',width=95):
        w=QLineEdit(str(value)); w.setFixedWidth(width); self.controls[name]=w; setattr(self,name,w); return w
    def choice(self,name,options,selected=0):
        w=QComboBox(); w.addItems(options); w.setCurrentIndex(selected); w.setMaxVisibleItems(8); self.controls[name]=w; setattr(self,name,w); return w
    def check(self,name,label,checked=False):
        w=QCheckBox(label); w.setChecked(checked); self.controls[name]=w; setattr(self,name,w); return w
    def fields(self,layout,fields):
        f=QFormLayout()
        for label,widget in fields: f.addRow(label,widget)
        layout.addLayout(f)

class ExtrapolationPanel(Panel):
    def __init__(self,parent=None):
        super().__init__(parent)
        modes=self.group('Xtrapol8 Modes')
        self.mode_buttons=QButtonGroup(self)
        self.mode_radios=[]
        r=QHBoxLayout(); modes.addLayout(r)
        for i,label in enumerate(('FoFo only','Fast and Furious','Calm and Curious')):
            btn=QRadioButton(label); self.mode_buttons.addButton(btn,i)
            self.mode_radios.append(btn); r.addWidget(btn)
        r.addStretch()
        g=self.group('Occupancies'); self.occupancies_group=g.parentWidget()

        occupancy_grid=QGridLayout(); occupancy_grid.setHorizontalSpacing(12); occupancy_grid.setVerticalSpacing(10)
        for col,(label,name,value) in enumerate((('Lowest','LowTextCtrl','0.1'),('Highest','HighTextCtrl','0.5'),('Steps','StepsTextCtrl','5'))):
            occupancy_grid.addWidget(QLabel(label),0,col)
            field=self.text(name,value,120); field.setMinimumWidth(95); field.setMaximumWidth(160)
            occupancy_grid.addWidget(field,1,col)
            occupancy_grid.setColumnStretch(col,1)
        occupancy_grid.addWidget(QLabel('Occupancies list'),2,0)
        occupancy_list=self.text('ListTextCtrl','',550)
        occupancy_list.setMinimumWidth(160); occupancy_list.setMaximumWidth(16777215)
        occupancy_grid.addWidget(occupancy_list,3,0,1,3)
        g.addLayout(occupancy_grid)
        g=self.group('Fourier difference map (FoFo)')
        self.row(g,'Type of FoFo Map :',self.choice('FoChoice',['qfofo','fofo','kfofo']),'K weight Scale :',self.text('kscale','0.05'))
        g=self.group('Extrapolated structure factors'); self.extrapolated_group=g.parentWidget()
        names=['qfextr','fextr','kfextr','qfgenick','fgenick','kfgenick','qfextr_calc','fextr_calc','kfextr_calc']
        # Compact, left-anchored controls: no stretching of checkboxes or buttons
        # when the main window becomes wider.
        factors=QHBoxLayout(); factors.setSpacing(26); g.addLayout(factors)
        grid=QGridLayout(); grid.setHorizontalSpacing(12); grid.setVerticalSpacing(8)
        factors.addLayout(grid)
        for i,name in enumerate(names):
            checkbox=self.check(name,name,name=='qfextr')
            checkbox.setMinimumWidth(112)
            grid.addWidget(checkbox,i//3,i%3)
        actions=QGridLayout(); actions.setHorizontalSpacing(8); actions.setVerticalSpacing(6)
        factors.addLayout(actions)
        action_definitions=[
            ('Select all',names),('Clear all',[]),
            ('All q-weighted',[x for x in names if x.startswith('q')]),
            ('All non-weighted',[x for x in names if x.startswith('f')]),
            ('All k-weighted',[x for x in names if x.startswith('k')]),
        ]
        self.selection_buttons=[]
        for i,(label,subset) in enumerate(action_definitions):
            button=QPushButton(label)
            button.setFixedWidth(145)
            button.setMinimumHeight(28)
            actions.addWidget(button,i//2,i%2)
            self.selection_buttons.append(button)
            button.clicked.connect(lambda _,subset=subset: [self.controls[n].setChecked(n in subset) for n in names])
        factors.addStretch(1)
        g=self.group('Map-explorer'); self.map_group=g.parentWidget()
        # Left-aligned, compact grid: labels and inputs stay together on wide windows.
        map_grid=QGridLayout()
        map_grid.setHorizontalSpacing(12)
        map_grid.setVerticalSpacing(9)
        map_fields=[
            ('Peak detection threshold :',self.text('peak_detection_thresholdTextCtrl','4.0')),
            ('Peak integration floor :',self.text('peak_integration_floorTextCtrl','3.5')),
            ('Radius (A) :',self.text('RadiusTextCtrl','1.5')),
            ('z-score :',self.text('ZscoreTextCtrl','2.0')),
            ('Occupancy estimation :',self.choice('OccEstimation',['difference_map_maximization','difference_map_PearsonCC','distance_analysis']))
        ]
        for row,(label,widget) in enumerate(map_fields):
            map_grid.addWidget(QLabel(label),row,0)
            map_grid.addWidget(widget,row,1)
        map_grid.setColumnStretch(0,0)
        map_grid.setColumnStretch(1,0)
        map_row=QHBoxLayout()
        map_row.addLayout(map_grid)
        map_row.addStretch(1)
        g.addLayout(map_row)
        self._occupancy_labels={'difference_map_maximization':'Map maximization','difference_map_PearsonCC':'Map Pearson CC','distance_analysis':'Distance analysis'}
        self.OccEstimation.clear()
        for key, label in self._occupancy_labels.items():
            self.OccEstimation.addItem(label, key)
        self.OccEstimation.setToolTip('Occupancy estimation method (full PHIL identifiers are preserved as item data)')
        self.OccEstimation.setMinimumContentsLength(17)
        self.expert_button=QPushButton('Expert Mode ?'); self.expert_button.setMinimumHeight(34); self.expert_button.setSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed); self.root.addWidget(self.expert_button, alignment=Qt.AlignLeft)
        g=self.group('Negative and Missing reflections'); self.negative_group=g.parentWidget()
        self.row(g,'Negative :',self.choice('negChoice',['truncate','reject','zero','fref','fcalc','keep','absolute']),'Missing :',self.choice('missChoice',['fill','no_fill']))
        g=self.group('Scaling'); self.scaling_group=g.parentWidget()
        self.row(g,'B-factor Scaling :',self.choice('ScalingChoice',['no','isotropic','anisotropic']))
        self.row(g,'Scaling high resolution :',self.text('ScalingHighRes',''),'Scaling low resolution :',self.text('ScalingLowRes',''))
        self.root.addStretch()
        self.expert=False
        self.fextr_names=names
        self.mode_buttons.buttonClicked[int].connect(self.set_mode)
        self.expert_button.clicked.connect(self.toggle_expert)
        self.FoChoice.currentTextChanged.connect(self.update_k_scale)
        self.mode_radios[2].setChecked(True)
        self.set_mode(2)

    def update_k_scale(self, *_):
        self.kscale.setEnabled(self.FoChoice.isEnabled() and self.FoChoice.currentText()=='kfofo')

    def toggle_expert(self):
        self.expert=not self.expert
        self.expert_button.setText('Expert Mode' if self.expert else 'Expert Mode ?')
        self.set_mode(self.mode_buttons.checkedId())

    def set_mode(self, mode):
        # Mirror panelExtrapolation.py's onFoFo/onFastNFurious/onCalmNCurious.
        fofo=(mode==0); fast=(mode==1)
        self.occupancies_group.setVisible(not fofo)
        self.extrapolated_group.setVisible(not fofo)
        self.negative_group.setVisible(self.expert and not fofo)
        self.scaling_group.setVisible(self.expert)
        self.OccEstimation.setVisible(not fofo)
        # QFormLayout labels are separate widgets: hide the entire row.
        form=self.OccEstimation.parentWidget()
        self.FoChoice.setEnabled(not fast)
        if fast:
            self.FoChoice.setCurrentText('qfofo')
            for name in self.fextr_names:
                self.controls[name].setChecked(name=='qfextr')
            self.negChoice.setCurrentText('truncate')
            self.missChoice.setCurrentText('fill')
        for name in self.fextr_names:
            self.controls[name].setEnabled(not fast)
        for btn in self.selection_buttons:
            btn.setEnabled(not fast)
        self.negChoice.setEnabled(not fast)
        self.missChoice.setEnabled(not fast)
        allowed=['difference_map_maximization','difference_map_PearsonCC']
        if not fast: allowed.append('distance_analysis')
        old=self.OccEstimation.currentData()
        self.OccEstimation.blockSignals(True)
        self.OccEstimation.clear()
        for key in allowed:
            self.OccEstimation.addItem(self._occupancy_labels[key], key)
        if old in allowed: self.OccEstimation.setCurrentIndex(allowed.index(old))
        self.OccEstimation.blockSignals(False)
        self.update_k_scale()


class RefinementPanel(Panel):
    def __init__(self,parent=None):
        super().__init__(parent)
        self.root.addWidget(self.check('RunRef','Run refinement',True))
        g=self.group('Real-space refinement')
        self.row(g,'Program :',self.choice('SoftChoiceReal',['phenix.real_space_refine','coot']), 'Cycles (Real Space) :',self.text('NCyclesReal_TextCtrl','3'))
        g=self.group('Reciprocal-space refinement')
        self.row(g,'Program :',self.choice('SoftChoiceReci',['phenix.refine','refmac5']))
        self.phenix_groups=[]
        self.refmac_groups=[]
        g=self.group('Phenix.refine Refinement strategy'); self.phenix_groups.append(g.parentWidget())
        opts=[('individual_sites','XYZ (Reciprocal space)'),('individual_sites_real_space','XYZ (Real space)'),('rigid_body','Rigid Body'),('individual_adp','Individual ADP'),('group_adp','Group ADP'),('tls','TLS'),('occupancies','Occupancies'),('group_anomalous','Group anomalous')]
        grid=QGridLayout(); g.addLayout(grid)
        for i,(name,label) in enumerate(opts): grid.addWidget(self.check(name,label,name in ('individual_sites','individual_adp')),i//2,i%2)
        self.row(g,'Cycles (Reciprocal space) :',self.text('NCyclesReciprocal_TextCtrl','5'))
        g=self.group('Phenix.refine Weights'); self.phenix_groups.append(g.parentWidget())
        self.fields(g,[(label,self.text(name,value)) for label,name,value in [('wxc_scale :','wxc_scale_TextCtrl','0.5'),('wxu_scale :','wxu_scale_TextCtrl','1.0'),('Bonds rmsd :','bonds_rmsd_TextCtrl','0.02'),('Angles rmsd :','angle_rmsd_TextCtrl','2.0'),('Rfree - Rwork :','rf_minus_rw','0.05')]])
        # Place Phenix refinement strategy and weights on one horizontal row.
        # Keep the original widgets so their state and signal connections are preserved.
        phenix_row = QHBoxLayout()
        phenix_row.setSpacing(12)
        for group in self.phenix_groups[:2]:
            self.root.removeWidget(group)
            group.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
            phenix_row.addWidget(group, 1)
            # Avoid centred forms and excessive stretching in the weights panel.
            for item_index in range(group.layout().count()):
                item = group.layout().itemAt(item_index)
                form = item.layout() if item is not None else None
                if isinstance(form, QFormLayout):
                    form.setFormAlignment(Qt.AlignLeft | Qt.AlignTop)
                    form.setLabelAlignment(Qt.AlignRight | Qt.AlignVCenter)
                    form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.root.addLayout(phenix_row)
        g=self.group('Phenix.refine Other options'); self.phenix_groups.append(g.parentWidget())
        self.row(g,*[self.check(n,l) for n,l in [('sim_ann','Simulated annealing'),('ordered_solvent','Ordered solvent'),('density_modification','Density Modification'),('map_sharpening','Map Sharpening')]])
        g=self.group('Simulated annealing'); self.phenix_groups.append(g.parentWidget())
        self.simulated_annealing_group = g.parentWidget()
        self.fields(g,[('Start T :',self.text('start_T','5000')),('Final T :',self.text('final_T','300')),('Cooling Rate :',self.text('cooling_rate','100')),('Mode :',self.choice('mode',['every_macro_cycle','second_and_before_last','once','first','first_half'],selected=1))])
        # Anchor the form to the left, with labels and inputs close together.
        annealing_form = g.itemAt(0).layout()
        annealing_form.setFormAlignment(Qt.AlignLeft | Qt.AlignTop)
        annealing_form.setLabelAlignment(Qt.AlignRight | Qt.AlignVCenter)
        annealing_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        g=self.group('DM Density Modification'); self.phenix_groups.append(g.parentWidget())
        self.row(g,'Combine :',self.choice('combine',['PERT','OMIT']),'Cycles :',self.text('cycles','10'))
        g=self.group('Refmac5 Refinement'); self.refmac_groups.append(g.parentWidget())
        self.fields(g,[('Refinement type :',self.choice('REFType',['REST','RIGID','TLS'])),('Number of cycles :',self.text('refmac_cycles','10')),('Bfactor refinement :',self.choice('Bref',['OVERall','ISOTropic'])),('Reset Bfactor :',self.text('bfac_set','')),('TLS Refinement :',self.choice('TLS',['True','False'])),('TLS cycles :',self.text('tls_cycles','10'))])
        g=self.group('Refmac5 Weights'); self.refmac_groups.append(g.parentWidget())
        self.fields(g,[('Weighting :',self.choice('AUTO',['AUTO','MATRIx'])),('Weight term :',self.text('weighting_term','')),('Experimental sigmas :',self.choice('NOEX',['NOEX','EXPE']))])
        g=self.group('Refmac5 Restraints'); self.refmac_groups.append(g.parentWidget())
        self.fields(g,[('Jelly body refinement :',self.choice('jelly_body_refinement',['True','False'])),('Jelly body sigmas :',self.text('jbs_TextCtrl','')),('Jelly body additional restraints :',self.text('jbar_TextCtrl','')),('External restraints :',self.text('extr_TextCtrl',''))])
        # Keep the three main Refmac5 sections side by side instead of stacked.
        # The group widgets retain their original controls and conditional visibility.
        refmac_columns = QHBoxLayout()
        refmac_columns.setSpacing(12)
        for group in self.refmac_groups[:3]:
            self.root.removeWidget(group)
            group.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
            refmac_columns.addWidget(group, 1)
            # Labels/fields stay close together rather than centred in wide panels.
            form = group.layout().itemAt(0).layout()
            if isinstance(form, QFormLayout):
                form.setFormAlignment(Qt.AlignLeft | Qt.AlignTop)
                form.setLabelAlignment(Qt.AlignRight | Qt.AlignVCenter)
                form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.root.addLayout(refmac_columns)
        g=self.group('Refmac5 Other options'); self.refmac_groups.append(g.parentWidget())
        self.row(g,'Combine :',self.choice('refmac_combine',['PERT','OMIT']),'Cycles :',self.text('refmac_DM_cycles','10'))
        self.row(g,*[self.check(n,l) for n,l in [('refmac_DM','Density modification'),('twinning','Twinning Refinement'),('REFMAC_map_sharpening','Map Sharpening')]])
        self.root.addStretch()
        self.SoftChoiceReal.currentTextChanged.connect(self.update_real_space_cycles)
        self.update_real_space_cycles()
        self.SoftChoiceReci.currentTextChanged.connect(self.update_reciprocal_sections)
        self.sim_ann.toggled.connect(self.update_reciprocal_sections)
        self.update_reciprocal_sections()

    def update_real_space_cycles(self, *_):
        self.NCyclesReal_TextCtrl.setEnabled(self.SoftChoiceReal.currentText() == 'phenix.real_space_refine')

    def update_reciprocal_sections(self, *_):
        """Show only controls belonging to the selected reciprocal refinement engine."""
        program = self.SoftChoiceReci.currentText()
        for group in self.phenix_groups:
            group.setVisible(program == 'phenix.refine')
        self.simulated_annealing_group.setVisible(program == 'phenix.refine' and self.sim_ann.isChecked())
        for group in self.refmac_groups:
            group.setVisible(program == 'refmac5')
