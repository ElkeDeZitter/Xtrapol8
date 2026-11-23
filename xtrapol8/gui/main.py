"""Qt GUI migration milestone 1: window, notebook and Input/Output layout.

Run standalone: python -m xtrapol8.gui.main
No scientific code is imported or modified.
"""
import sys
import os
import tempfile
import uuid
import io
from dataclasses import dataclass
from pathlib import Path
from PyQt5.QtCore import Qt, QSize, QLineF
from .action_icons import action_icon
from .results import ResultsViewer
from .progress import ProgressPanel
from .run_paths import resolve_inputs, unused_output_directory
from .process import FextrProcess
from .run_store import RunStore
from . import phil_bridge
from .panels import ExtrapolationPanel, RefinementPanel
from PyQt5.QtGui import QIcon, QFontDatabase, QColor, QBrush, QPixmap, QPainter, QPen
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
    QHBoxLayout, QTabWidget, QGroupBox, QTableWidget, QTableWidgetItem,
    QPushButton, QLabel, QLineEdit, QFileDialog, QAction, QToolBar,
    QMessageBox, QStyle, QStyleOptionViewItem, QHeaderView, QPlainTextEdit, QProgressBar, QStyledItemDelegate, QComboBox, QSplitter, QMenu, QAbstractItemView, QInputDialog)



class TailElideDelegate(QStyledItemDelegate):
    """Show the end of long identifiers and paths, preserving the full tooltip."""

    def paint(self, painter, option, index):
        opt = QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)
        opt.text = opt.fontMetrics.elidedText(opt.text, Qt.ElideLeft,
                                             max(0, opt.rect.width() - 12))
        opt.textElideMode = Qt.ElideNone
        opt.widget.style().drawControl(QStyle.CE_ItemViewItem, opt, painter, opt.widget)


def locate_icon():
    """Magnifying glass icon drawn in Qt, independent of OS icon themes."""
    pixmap = QPixmap(20, 20)
    pixmap.fill(Qt.transparent)
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.Antialiasing)
    painter.setPen(QPen(QColor('#315b83'), 2.0, Qt.SolidLine, Qt.RoundCap, Qt.RoundJoin))
    painter.drawEllipse(3, 2, 11, 11)
    painter.drawLine(12, 12, 18, 18)
    painter.end()
    return QIcon(pixmap)


def config_icon():
    """Compact file-with-cog glyph, consistent with the other 20px action icons."""
    pixmap = QPixmap(20, 20)
    pixmap.fill(Qt.transparent)
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.Antialiasing)
    pen = QPen(QColor('#315b83'), 1.6, Qt.SolidLine, Qt.RoundCap, Qt.RoundJoin)
    painter.setPen(pen)
    painter.drawLine(3, 2, 11, 2)
    painter.drawLine(11, 2, 15, 6)
    painter.drawLine(15, 6, 15, 9)
    painter.drawLine(3, 2, 3, 18)
    painter.drawLine(3, 18, 10, 18)
    painter.drawLine(11, 2, 11, 6)
    painter.drawLine(11, 6, 15, 6)
    painter.drawEllipse(10, 11, 7, 7)
    painter.drawLine(QLineF(13.5, 9.5, 13.5, 11))
    painter.drawLine(QLineF(13.5, 18, 13.5, 19))
    painter.drawLine(QLineF(8.5, 14.5, 10, 14.5))
    painter.drawLine(QLineF(17, 14.5, 19, 14.5))
    painter.end()
    return QIcon(pixmap)


def trash_icon():
    """Small monochrome bin icon without external assets or platform fonts."""
    pixmap = QPixmap(20, 20)
    pixmap.fill(Qt.transparent)
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.Antialiasing)
    painter.setPen(QPen(QColor('#a22b36'), 1.7, Qt.SolidLine, Qt.RoundCap, Qt.RoundJoin))
    painter.drawLine(4, 6, 16, 6)
    painter.drawLine(8, 4, 12, 4)
    painter.drawLine(6, 8, 7, 17)
    painter.drawLine(14, 8, 13, 17)
    painter.drawLine(7, 17, 13, 17)
    painter.drawLine(9, 9, 9, 15)
    painter.drawLine(11, 9, 11, 15)
    painter.end()
    return QIcon(pixmap)


class InputTypeDelegate(QStyledItemDelegate):
    TYPES = ('Reference model', 'Reference mtz', 'Triggered mtz', 'Restraints')

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        editor.addItems(self.TYPES)
        return editor

    def setEditorData(self, editor, index):
        value = index.data()
        editor.setCurrentIndex(max(0, editor.findText(value)))

    def setModelData(self, editor, model, index):
        model.setData(index, editor.currentText())


class InputOutputPanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.files = {"Reference model": [], "Reference mtz": [],
                      "Triggered mtz": [], "Restraints": []}
        root = QVBoxLayout(self)
        inputs = QGroupBox('Input Files')
        in_layout = QVBoxLayout(inputs)
        self.list = QTableWidget(0, 2)
        self.list.setHorizontalHeaderLabels(['File Path', 'Data type'])
        self.list.setToolTip('Right-click a Data type cell to select its type')
        self.list.setColumnWidth(0, 700)
        self.list.setColumnWidth(1, 150)
        self.list.setItemDelegateForColumn(1, InputTypeDelegate(self.list))
        self.list.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.list.setContextMenuPolicy(Qt.CustomContextMenu)
        self.list.customContextMenuRequested.connect(self.show_input_type_menu)
        self.list.setSelectionBehavior(QTableWidget.SelectRows)
        self.list.setSelectionMode(QTableWidget.SingleSelection)
        in_layout.addWidget(self.list)
        buttons = QHBoxLayout()
        self.AddButton = QPushButton('Add file')
        self.DelButton = QPushButton('Remove file')
        for b in (self.AddButton, self.DelButton):
            b.setFixedSize(100, 30)
        buttons.addWidget(self.AddButton)
        buttons.addSpacing(10)
        buttons.addWidget(self.DelButton)
        buttons.addSpacing(30)
        buttons.addWidget(QLabel('Low resolution :'))
        self.lowRes = QLineEdit()
        self.lowRes.setFixedSize(80, 30)
        buttons.addWidget(self.lowRes)
        buttons.addSpacing(30)
        buttons.addWidget(QLabel('High resolution :'))
        self.highRes = QLineEdit()
        self.highRes.setFixedSize(80, 30)
        buttons.addWidget(self.highRes)
        buttons.addStretch()
        in_layout.addLayout(buttons)
        root.addWidget(inputs)
        outputs = QGroupBox('Output Files')
        out_layout = QVBoxLayout(outputs)
        out_row = QHBoxLayout()
        out_row.addWidget(QLabel('Output Directory :'))
        self.outdir = QLineEdit()
        out_row.addWidget(self.outdir, 1)
        browse = QPushButton('Browse')
        browse.clicked.connect(self.browse_directory)
        out_row.addWidget(browse)
        out_layout.addLayout(out_row)
        prefix_row = QHBoxLayout()
        prefix_label = QLabel('Output Prefix :')
        prefix_label.setFixedWidth(150)
        prefix_row.addWidget(prefix_label)
        self.outname = QLineEdit()
        self.outname.setFixedWidth(400)
        prefix_row.addWidget(self.outname)
        prefix_row.addStretch()
        out_layout.addLayout(prefix_row)
        name_row = QHBoxLayout()
        name_label = QLabel('Run Name :')
        name_label.setFixedWidth(150)
        name_row.addWidget(name_label)
        self.run_name = QLineEdit()
        self.run_name.setFixedWidth(400)
        self.run_name.setPlaceholderText('Optional display name (does not affect output files)')
        name_row.addWidget(self.run_name)
        name_row.addStretch()
        out_layout.addLayout(name_row)
        root.addWidget(outputs)
        root.addStretch()
        self.AddButton.clicked.connect(self.add_file)
        self.DelButton.clicked.connect(self.remove_file)

    def show_input_type_menu(self, position):
        """Open the data-type dropdown on right-click (no double-click required)."""
        index = self.list.indexAt(position)
        if not index.isValid() or index.column() != 1:
            return
        row = index.row()
        item = self.list.item(row, 1)
        current = item.text() if item else 'Unassigned'
        menu = QMenu(self.list)
        for data_type in InputTypeDelegate.TYPES:
            action = menu.addAction(data_type)
            action.setCheckable(True)
            action.setChecked(data_type == current)
        selected = menu.exec_(self.list.viewport().mapToGlobal(position))
        if selected is not None:
            if item is None:
                item = QTableWidgetItem()
                self.list.setItem(row, 1, item)
            item.setText(selected.text())

    def browse_directory(self):
        path = QFileDialog.getExistingDirectory(self, 'Output Directory', self.outdir.text())
        if path:
            self.outdir.setText(path)

    def add_file(self):
        paths, _ = QFileDialog.getOpenFileNames(self, 'Add input files')
        for path in paths:
            row = self.list.rowCount()
            self.list.insertRow(row)
            self.list.setItem(row, 0, QTableWidgetItem(path))
            self.list.setItem(row, 1, QTableWidgetItem('Unassigned'))

    def remove_file(self):
        for row in sorted({i.row() for i in self.list.selectedIndexes()}, reverse=True):
            self.list.removeRow(row)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('XtrapolG8 -- Qt migration (QProcess)')
        self.runner = None  # Legacy alias: the most recently launched run.
        self.runs = []
        self.run_store = RunStore()
        self.selected_run = None
        self.phil_file = None
        self._configuration_snapshot = None
        self._phil_params = None
        self._run_phil_path = None
        self._results_source_directory = None
        self._input_base_directory = None
        self.resize(1100, 1000)
        self.source_phil_label = QLabel('Configuration: not loaded')
        self.source_phil_label.setToolTip('PHIL file currently loaded into the configuration panels')
        self.statusBar().addPermanentWidget(self.source_phil_label, 1)
        file_menu = self.menuBar().addMenu('&File')
        for title in ('Open phil', 'Load Results', 'Save phil'):
            action = QAction(title, self)
            action.triggered.connect(lambda checked=False, name=title: self.open_phil() if name == 'Open phil' else self.save_phil() if name == 'Save phil' else self.load_results())
            file_menu.addAction(action)
        file_menu.addSeparator()
        file_menu.addAction('Quit', self.close)
        toolbar = QToolBar(self)
        self.addToolBar(toolbar)
        toolbar.setIconSize(QSize(32, 32))
        toolbar.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        toolbar.setMinimumHeight(65)
        for title, kind in [('Run', 'run'), ('Stop', 'stop')]:
            action = QAction(action_icon(kind), title, self)
            action.triggered.connect(self.run_fextr if kind == 'run' else self.stop_fextr)
            toolbar.addAction(action)
        self.notebook = QTabWidget()
        self.setCentralWidget(self.notebook)
        configure = QTabWidget()
        self.tabIO = InputOutputPanel()
        configure.addTab(self.tabIO, 'Input / Output')
        self.tabExtrapolation = ExtrapolationPanel()
        self.tabRefinement = RefinementPanel()
        configure.addTab(self.tabExtrapolation, 'FoFo / Extrapolation')
        configure.addTab(self.tabRefinement, 'Refinement')
        self.notebook.addTab(configure, 'Configure')
        self.log_view = None
        self.progress_panel = None
        self.run_results = None
        self.run_count = 0
        self._run_logs = []
        self._run_viewers = []
        self._run_progress = []
        self._build_run_manager()
        self._build_results_tab()
        self._restore_runs()

    def _build_run_manager(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        self.run_table = QTableWidget(0, 5)
        self.run_table.setHorizontalHeaderLabels(['Run Name', 'Status', 'Progress', 'Output directory', 'Actions'])
        self.run_table.setItemDelegateForColumn(0, TailElideDelegate(self.run_table))
        self.run_table.setItemDelegateForColumn(3, TailElideDelegate(self.run_table))
        self.run_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.run_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.Stretch)
        self.run_table.horizontalHeader().setSectionResizeMode(4, QHeaderView.Fixed)
        self.run_table.setColumnWidth(4, 208)
        self.run_table.verticalHeader().setDefaultSectionSize(48)
        self.run_table.verticalHeader().setMinimumSectionSize(48)
        self.run_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.run_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.run_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.run_table.itemSelectionChanged.connect(self._manager_selection_changed)
        layout.addWidget(self.run_table)
        actions = QHBoxLayout()
        for title, callback in [('View Results', self.view_selected_results),
                                ('Rename', self.rename_selected_run),
                                ('Stop all', self.stop_all_runs)]:
            button = QPushButton(title)
            button.clicked.connect(callback)
            actions.addWidget(button)
        actions.addStretch()
        layout.addLayout(actions)
        self.manager_log = QPlainTextEdit()
        self.manager_log.setReadOnly(True)
        self.manager_log.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.manager_log.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
        layout.addWidget(QLabel('Selected run log'))
        layout.addWidget(self.manager_log)
        self.notebook.addTab(page, 'Run Manager')

    def _build_results_tab(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        selector = QHBoxLayout()
        selector.addWidget(QLabel('Run:'))
        self.results_selector = QComboBox()
        self.results_selector.currentIndexChanged.connect(self._results_selection_changed)
        selector.addWidget(self.results_selector, 1)
        back = QPushButton('Back to Run Manager')
        back.clicked.connect(lambda: self.notebook.setCurrentIndex(1))
        selector.addWidget(back)
        layout.addLayout(selector)
        self.results_host = QVBoxLayout()
        layout.addLayout(self.results_host, 1)
        self.empty_results = QLabel('No run selected. Start a calculation or load existing results.')
        self.empty_results.setAlignment(Qt.AlignCenter)
        self.results_host.addWidget(self.empty_results)
        self.notebook.addTab(page, 'Results')

    def _unique_run_name(self, proposed):
        """Keep user-facing run names unique, using 00, 01, ... for collisions."""
        proposed = proposed.strip() or 'Run'
        used = {entry['name'] for entry in self.runs}
        if proposed not in used:
            return proposed
        index = 0
        while f'{proposed}_{index:02d}' in used:
            index += 1
        return f'{proposed}_{index:02d}'

    def _add_run(self, name, directory, runner=None, progress=None, *, restored=None):
        name = self._unique_run_name(name)
        self.run_count += 1
        entry = {'id': restored['id'] if restored else uuid.uuid4().hex,
                 'name': name, 'directory': Path(directory), 'runner': runner,
                 'status': restored['status'] if restored else ('Starting' if runner else 'Loaded'),
                 'progress_value': restored['progress_value'] if restored else (0 if runner else 100),
                 'progress': progress,
                 'lines': [], 'viewer': ResultsViewer(), 'log_widget': None}
        if restored is not None:
            entry['lines'] = self._read_saved_log(entry)
        entry['viewer'].directory = Path(directory)
        if restored is not None and not Path(directory).is_dir():
            entry['status'] = 'Missing'
        if Path(directory).is_dir():
            entry['viewer'].refresh()
        self.runs.append(entry)
        row = self.run_table.rowCount()
        self.run_table.insertRow(row)
        for col, value in enumerate([name, entry['status'], '', str(directory)]):
            item = QTableWidgetItem(value)
            if col in (0, 3):
                item.setToolTip(value)
            self.run_table.setItem(row, col, item)
        progress_cell = QWidget()
        progress_layout = QHBoxLayout(progress_cell)
        progress_layout.setContentsMargins(9, 7, 9, 7)
        bar = QProgressBar()
        bar.setMinimumHeight(25)
        bar.setTextVisible(True)
        progress_layout.addWidget(bar)
        bar.setRange(0, 100)
        bar.setValue(entry['progress_value'])
        bar.setFormat('%p%' if runner else '—')
        entry['table_bar'] = bar
        self.run_table.setCellWidget(row, 2, progress_cell)
        actions = QWidget()
        buttons = QHBoxLayout(actions)
        buttons.setContentsMargins(9, 6, 9, 6)
        buttons.setSpacing(7)
        stop = QPushButton()
        stop.setIcon(self.style().standardIcon(QStyle.SP_MediaStop))
        stop.setIconSize(QSize(16, 16))
        stop.setToolTip('Stop this run')
        stop.setFixedSize(33, 31)
        stop.setStyleSheet('QPushButton { border: 1px solid #b8c0ca; border-radius: 5px; background: #f6f8fb; } QPushButton:hover { background: #fce7e7; border-color: #d48a8a; } QPushButton:disabled { opacity: 0.4; }')
        stop.clicked.connect(lambda _checked=False, e=entry: self._stop_run(e))
        delete = QPushButton()
        delete.setIcon(trash_icon())
        delete.setIconSize(QSize(17, 17))
        delete.setToolTip('Remove this run from the list (files are kept)')
        delete.setFixedSize(33, 31)
        delete.setStyleSheet('QPushButton { border: 1px solid #b8c0ca; border-radius: 5px; background: #f6f8fb; } QPushButton:hover { background: #fce7e7; border-color: #d48a8a; }')
        delete.clicked.connect(lambda _checked=False, e=entry: self._remove_run(e))
        # One contextual action: configuration for a valid run, locate for Missing.
        config_action = QPushButton()
        config_action.setIconSize(QSize(18, 18))
        config_action.setFixedSize(33, 31)
        config_action.setStyleSheet('QPushButton { border: 1px solid #b8c0ca; border-radius: 5px; background: #f6f8fb; } QPushButton:hover { background: #e7f0fc; border-color: #8aa8d4; }')
        config_action.clicked.connect(lambda _checked=False, e=entry: self._run_config_or_locate(e))
        buttons.addWidget(stop)
        buttons.addWidget(config_action)
        buttons.addWidget(delete)
        entry['config_action_button'] = config_action
        entry['stop_button'] = stop
        self.run_table.setCellWidget(row, 4, actions)
        self._update_run_status(entry)
        self.results_selector.addItem(name)
        self.run_table.selectRow(row)
        self._select_run(row)
        if restored is None or entry['status'] == 'Missing':
            self._persist_run(entry)
        return entry

    def _persist_run(self, entry):
        if entry.get('progress') is not None:
            entry['progress_value'] = entry['progress'].percent_complete()
        elif entry['status'] == 'Completed':
            entry['progress_value'] = 100
        self.run_store.save(entry)

    def _restore_runs(self):
        for record in self.run_store.all():
            if record['status'] in ('Running', 'Starting', 'Stopping'):
                record['status'] = 'Interrupted'
                self.run_store.save(record)
            self._add_run(record['name'], record['directory'], restored=record)

    def _select_run(self, index):
        if index < 0 or index >= len(self.runs):
            return
        entry = self.runs[index]
        self._check_run_directory(entry)
        self.selected_run = entry
        if not (entry.get('runner') and entry['runner'].running):
            entry['lines'] = self._read_saved_log(entry)
        self.manager_log.setPlainText('\n'.join(entry['lines']))
        self.results_selector.blockSignals(True)
        self.results_selector.setCurrentIndex(index)
        self.results_selector.blockSignals(False)
        self._show_results_entry(entry)

    def _check_run_directory(self, entry):
        """Detect a removed or detached results directory without deleting history."""
        if entry.get('runner') and entry['runner'].running:
            return
        if not entry['directory'].is_dir() and entry['status'] != 'Missing':
            entry['status'] = 'Missing'
            self._update_run_status(entry)

    def _run_config_or_locate(self, entry):
        if entry['status'] == 'Missing':
            self._locate_run(entry)
        else:
            self._load_run_configuration(entry)

    def _locate_run(self, entry):
        if entry.get('runner') and entry['runner'].running:
            QMessageBox.information(self, 'Run active', 'Wait for this run to finish before relocating it.')
            return
        selected = QFileDialog.getExistingDirectory(
            self, 'Locate existing Xtrapol8 results directory', str(entry['directory'].parent))
        if not selected:
            return
        directory = Path(selected).resolve()
        if not directory.is_dir():
            QMessageBox.warning(self, 'Directory missing', 'The selected directory is not accessible.')
            return
        # Require an Xtrapol8 output marker to avoid accidentally linking unrelated folders.
        # Never associate arbitrary directories with an existing run.
        phil_file = directory / 'Xtrapol8_in.phil'
        if not phil_file.is_file():
            QMessageBox.warning(self, 'Invalid Xtrapol8 results directory',
                'The selected directory has no Xtrapol8_in.phil. The run was not changed.')
            return
        try:
            from libtbx import phil as libtbx_phil
            parsed = libtbx_phil.parse(phil_file.read_text(encoding='utf-8'))
            if not any(scope.name in ('input', 'output', 'f_and_maps') for scope in parsed.objects):
                raise ValueError('Missing Xtrapol8 parameter scopes')
        except Exception as exc:
            QMessageBox.warning(self, 'Invalid Xtrapol8 configuration',
                f'Cannot validate Xtrapol8_in.phil: {exc}\nThe run was not changed.')
            return
        # Interrupted calculations may have no completed maps or log yet.
        entry['directory'] = directory
        entry['viewer'].directory = directory
        entry['viewer'].refresh()
        entry['status'] = 'Loaded' if entry['status'] == 'Missing' else entry['status']
        index = self.runs.index(entry)
        self.run_table.item(index, 3).setText(str(directory))
        self.run_table.item(index, 3).setToolTip(str(directory))
        self._update_run_status(entry)
        if self.selected_run is entry:
            self._show_results_entry(entry)

    def _show_results_entry(self, entry):
        while self.results_host.count():
            item = self.results_host.takeAt(0)
            widget = item.widget()
            if widget:
                widget.setParent(None)
        tabs = entry.get('tabs')
        if tabs is None:
            tabs = QTabWidget()
            log = QPlainTextEdit()
            log.setReadOnly(True)
            log.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
            log.setLineWrapMode(QPlainTextEdit.NoWrap)
            log.setPlainText('\n'.join(entry['lines']))
            entry['log_widget'] = log
            log_page = QWidget()
            log_layout = QVBoxLayout(log_page)
            log_actions = QHBoxLayout()
            refresh_log = QPushButton('Refresh log')
            refresh_log.setToolTip('Reload the Xtrapol8 log from the results directory')
            refresh_log.clicked.connect(lambda _checked=False, e=entry: self._refresh_saved_log(e))
            log_actions.addWidget(refresh_log)
            log_actions.addStretch()
            log_layout.addLayout(log_actions)
            splitter = QSplitter(Qt.Horizontal)
            splitter.addWidget(log)
            if entry['progress'] is not None:
                splitter.addWidget(entry['progress'])
            log_layout.addWidget(splitter)
            tabs.addTab(log_page, 'Log')
            viewer = entry['viewer']
            tabs.addTab(viewer, 'Main')
            tabs.addTab(viewer.occupancies_page, 'Occupancies')
            entry['tabs'] = tabs
        self.results_host.addWidget(tabs)

    def _manager_selection_changed(self):
        selected = self.run_table.selectionModel().selectedRows()
        if selected:
            self._select_run(selected[0].row())

    def _results_selection_changed(self, index):
        if 0 <= index < len(self.runs):
            self._select_run(index)
            self.run_table.selectRow(index)

    def view_selected_results(self):
        if self.selected_run is not None:
            self._check_run_directory(self.selected_run)
            if self.selected_run['status'] != 'Missing':
                self.selected_run['viewer'].refresh()
            self.notebook.setCurrentIndex(2)

    def rename_selected_run(self):
        entry = self.selected_run
        if entry is None:
            return
        name, ok = QInputDialog.getText(self, 'Rename run', 'Run Name:', text=entry['name'])
        if ok and name.strip():
            entry['name'] = (entry['name'] if name.strip() == entry['name'] else self._unique_run_name(name.strip()))
            index = self.runs.index(entry)
            self.run_table.item(index, 0).setText(entry['name'])
            self.run_table.item(index, 0).setToolTip(entry['name'])
            self.results_selector.setItemText(index, entry['name'])
            self._persist_run(entry)

    def _stop_run(self, entry):
        if entry in self.runs and entry['runner'] and entry['runner'].running:
            entry['status'] = 'Stopping'
            self._update_run_status(entry)
            entry['runner'].stop()

    def stop_selected_run(self):
        if self.selected_run is not None:
            self._stop_run(self.selected_run)

    def _remove_run(self, entry):
        if entry not in self.runs:
            return
        if entry['runner'] and entry['runner'].running:
            QMessageBox.information(self, 'Run active', 'Stop this run before removing it from the list.')
            return
        if QMessageBox.question(self, 'Remove run',
                'Remove this run from the lists? Files and results on disk will be preserved.',
                QMessageBox.Yes | QMessageBox.No) != QMessageBox.Yes:
            return
        index = self.runs.index(entry)
        if self.selected_run is entry:
            self.selected_run = None
        self.run_table.removeRow(index)
        self.results_selector.removeItem(index)
        self.runs.pop(index)
        self.run_store.remove(entry['id'])
        if self.runs:
            self.run_table.selectRow(min(index, len(self.runs)-1))
            self._select_run(min(index, len(self.runs)-1))
        else:
            self.manager_log.clear()
            while self.results_host.count():
                item = self.results_host.takeAt(0)
                if item.widget():
                    item.widget().setParent(None)
            self.results_host.addWidget(self.empty_results)
            self.empty_results.show()


    def stop_all_runs(self):
        active = [r for r in self.runs if r['runner'] and r['runner'].running]
        if active and QMessageBox.question(self, 'Stop all runs',
                    f'Stop all {len(active)} running calculations?',
                    QMessageBox.Yes | QMessageBox.No) == QMessageBox.Yes:
            for entry in active:
                entry['status'] = 'Stopping'
                self._update_run_status(entry)
                entry['runner'].stop()

    def _update_run_status(self, entry):
        index = self.runs.index(entry)
        self.run_table.item(index, 1).setText(entry['status'])
        palette = {'Running': ('#e8f2ff', '#1556a0'), 'Starting': ('#e8f2ff', '#1556a0'),
                   'Stopping': ('#fff1da', '#895000'), 'Completed': ('#e3f5e9', '#176c3c'),
                   'Stopped': ('#eeeeee', '#555555'), 'Failed': ('#ffe7e7', '#a02727'),
                   'Loaded': ('#edf0f5', '#46566b'), 'Interrupted': ('#fff1da', '#895000'),
                   'Missing': ('#ffe7e7', '#a02727')}
        bg, fg = palette.get(entry['status'], ('#ffffff', '#333333'))
        for col in (0, 1, 3):
            item = self.run_table.item(index, col)
            if item:
                item.setBackground(QBrush(QColor(bg)))
                item.setForeground(QBrush(QColor(fg)))
        entry['stop_button'].setEnabled(bool(entry['runner'] and entry['runner'].running))
        config_action = entry['config_action_button']
        if entry['status'] == 'Missing':
            config_action.setIcon(locate_icon())
            config_action.setToolTip('Locate missing Xtrapol8 results directory')
            config_action.setEnabled(True)
        else:
            config_action.setIcon(config_icon())
            config_action.setToolTip('Load Xtrapol8_in.phil into Configure (new run)')
            config_action.setEnabled(self._find_results_phil(entry['directory']) is not None)
        self._update_run_progress(entry)
        self._persist_run(entry)

    def _update_run_progress(self, entry):
        bar = entry.get('table_bar')
        if bar is None:
            return
        if entry['progress'] is not None:
            bar.setValue(entry['progress'].percent_complete())
            bar.setFormat('%p%')
        else:
            bar.setValue(entry.get('progress_value', 0))
            bar.setFormat('%p%')

    @staticmethod
    def _scientific_log_path(directory):
        """Find the Xtrapol8 log in the results directory; never create copies."""
        directory = Path(directory)
        if not directory.is_dir():
            return None
        candidates = [p for p in directory.glob('*_Xtrapol8.log') if p.is_file()]
        if not candidates:
            return None
        # Prefer the most recently modified log when a directory contains retries.
        return max(candidates, key=lambda p: (p.stat().st_mtime_ns, p.name))

    def _read_saved_log(self, entry):
        path = self._scientific_log_path(entry['directory'])
        if path is None:
            return []
        try:
            with path.open('r', encoding='utf-8', errors='replace') as stream:
                from collections import deque
                return list(deque((line.rstrip('\n') for line in stream), maxlen=20000))
        except OSError:
            return []

    def _refresh_saved_log(self, entry):
        entry['lines'] = self._read_saved_log(entry)
        if entry.get('log_widget') is not None:
            entry['log_widget'].setPlainText('\n'.join(entry['lines']))
        if self.selected_run is entry:
            self.manager_log.setPlainText('\n'.join(entry['lines']))

    def _load_run_configuration(self, entry):
        """Load the existing scientific PHIL without altering the source run."""
        self._check_run_directory(entry)
        if entry['status'] == 'Missing':
            QMessageBox.information(self, 'Results unavailable', 'Locate the results directory first.')
            return
        phil_path = self._find_results_phil(entry['directory'])
        if phil_path is None:
            QMessageBox.warning(self, 'Configuration unavailable',
                                'Xtrapol8_in.phil was not found in this results directory.')
            return
        if not self._confirm_replace_configuration():
            return
        try:
            phil_bridge.read_into(self, phil_path)
            self._set_loaded_phil(phil_path)
            self._results_source_directory = Path(entry['directory'])
            self.notebook.setCurrentIndex(0)
            self.statusBar().showMessage(f'Configuration loaded from {phil_path}', 12000)
        except Exception as exc:
            QMessageBox.warning(self, 'PHIL load failed', str(exc))

    def _run_log(self, entry, line):
        entry['lines'].append(line)
        if len(entry['lines']) > 20000:
            del entry['lines'][:1000]
        if entry.get('log_widget'):
            entry['log_widget'].appendPlainText(line)
        if self.selected_run is entry:
            self.manager_log.appendPlainText(line)
        if entry['progress'] is not None:
            entry['progress'].ingest(line)
            self._update_run_progress(entry)
            if len(entry['lines']) % 25 == 0:
                self._persist_run(entry)

    def _run_finished(self, entry, code, stopped):
        entry['status'] = 'Stopped' if stopped else ('Completed' if code == 0 else 'Failed')
        self._update_run_status(entry)
        self._run_log(entry, f'[{entry["status"]}: exit code {code}]')
        if entry['progress'] is not None:
            entry['progress'].finish(code, stopped)
            self._update_run_progress(entry)
        entry['viewer'].refresh()
        self.statusBar().showMessage(f'{entry["name"]}: {entry["status"]}')

    def _run_error(self, entry, message):
        self._run_log(entry, f'[Process error] {message}')
        self.statusBar().showMessage(f'{entry["name"]}: {message}')

    def open_phil(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select PHIL configuration', '', 'PHIL (*.phil);;All files (*)')
        if path:
            try:
                phil_bridge.read_into(self, path)
            except (RuntimeError, ValueError, OSError, AttributeError) as exc:
                QMessageBox.warning(self, 'PHIL load failed', str(exc))
                return
            self._set_loaded_phil(path)
            self._results_source_directory = None
            self._input_base_directory = Path(path).resolve().parent
            self.statusBar().showMessage(f'PHIL loaded: {path}')

    def save_phil(self):
        path, _ = QFileDialog.getSaveFileName(self, 'Save PHIL configuration', self.phil_file or 'Xtrapol8_config.phil', 'PHIL (*.phil);;All files (*)')
        if not path:
            return
        try:
            phil_bridge.save_from(self, path)
        except (RuntimeError, ValueError, OSError, AttributeError) as exc:
            QMessageBox.warning(self, 'PHIL save failed', str(exc))
            return
        self._set_loaded_phil(path)
        self.statusBar().showMessage(f'PHIL saved: {path}')

    def _configuration_signature(self):
        """Canonical GUI PHIL; compare to detect unsaved configuration edits."""
        master_phil, _ = phil_bridge.master()
        params = phil_bridge.collect(self)
        output = io.StringIO()
        master_phil.format(python_object=params).show(out=output)
        return output.getvalue()

    def _set_loaded_phil(self, path):
        self.phil_file = str(Path(path).resolve())
        self._configuration_snapshot = self._configuration_signature()
        self.source_phil_label.setText(f'Configuration: {Path(path).name}')
        self.source_phil_label.setToolTip(self.phil_file)

    def _confirm_replace_configuration(self):
        if self._configuration_snapshot is None:
            return True
        try:
            changed = self._configuration_signature() != self._configuration_snapshot
        except (RuntimeError, ValueError, AttributeError) as exc:
            QMessageBox.warning(self, 'Configuration check failed', str(exc))
            return False
        if not changed:
            return True
        answer = QMessageBox.question(
            self, 'Unsaved configuration changes',
            'The configuration panels contain unsaved changes. '
            'Replace them with the PHIL from the selected results directory?',
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        return answer == QMessageBox.Yes

    @staticmethod
    def _find_results_phil(directory):
        # Do not silently choose an unrelated PHIL from a neighbouring run.
        candidate = directory / 'Xtrapol8_in.phil'
        if candidate.is_file():
            return candidate
        return None

    def _input_file_paths(self):
        paths = []
        for row in range(self.tabIO.list.rowCount()):
            item = self.tabIO.list.item(row, 0)
            if item and item.text().strip():
                paths.append(item.text().strip())
        return paths

    def _resolve_run_inputs(self, workdir):
        paths = self._input_file_paths()
        # A loaded results PHIL is commonly inside the output directory,
        # while its relative inputs live in the original data directory.
        bases = [self._input_base_directory, workdir]
        if self._results_source_directory is not None:
            results = self._results_source_directory
            bases.extend([results.parent, results.parent.parent])
        bases = list(dict.fromkeys(Path(p).resolve() for p in bases if p))
        resolved, missing, ambiguous = resolve_inputs(paths, bases)
        if missing or ambiguous:
            details = '\n'.join(missing + [f'{p} (multiple matches)' for p in ambiguous])
            QMessageBox.information(
                self, 'Locate input files',
                'Some relative input files could not be uniquely located:\n'
                + details + '\n\nSelect their original input-data directory.')
            selected = QFileDialog.getExistingDirectory(
                self, 'Select original input-data directory', str(workdir))
            if not selected:
                raise ValueError('Input directory selection cancelled; run not started')
            # Explicitly chosen directory takes precedence over guessed bases.
            resolved, missing, ambiguous = resolve_inputs(paths, [selected])
            if missing or ambiguous:
                raise FileNotFoundError('Input files not found or ambiguous in selected directory: '
                                        + ', '.join(missing + list(ambiguous)))
            self._input_base_directory = Path(selected).resolve()
        # Update GUI only after all inputs have been validated. The original
        # PHIL on disk is never overwritten.
        for row in range(self.tabIO.list.rowCount()):
            item = self.tabIO.list.item(row, 0)
            if item and item.text().strip() in resolved:
                item.setText(resolved[item.text().strip()])

    def run_fextr(self):
        try:
            workdir = (self._input_base_directory or
                       (Path(self.phil_file).resolve().parent if self.phil_file else Path.cwd()))
            workdir = Path(workdir).resolve()
            self._resolve_run_inputs(workdir)
            outdir = self.tabIO.outdir.text().strip()
            if not outdir:
                raise ValueError('Please specify an output directory before running')
            requested = Path(outdir).expanduser()
            if not requested.is_absolute():
                requested = workdir / requested
            target = unused_output_directory(requested)
            self.tabIO.outdir.setText(str(target))
            run_dir = Path(tempfile.mkdtemp(prefix='xtrapol8_qt_'))
            phil_path = phil_bridge.save_from(self, run_dir / 'run.phil')
            name = self.tabIO.run_name.text().strip() or self.tabIO.outname.text().strip() or f'Run #{self.run_count + 1}'
            runner = FextrProcess(self)
            try:
                progress = ProgressPanel(phil_bridge.collect(self))
            except (RuntimeError, ValueError, AttributeError, TypeError):
                progress = None
            entry = self._add_run(name, target, runner, progress)
            runner.log_received.connect(lambda line, e=entry: self._run_log(e, line))
            runner.started.connect(lambda e=entry: (e.update(status='Running'), self._update_run_status(e)))
            runner.finished.connect(lambda code, stopped, e=entry: self._run_finished(e, code, stopped))
            runner.error.connect(lambda message, e=entry: self._run_error(e, message))
            self.runner = runner
            self._run_log(entry, f'[Output directory] {target}')
            runner.start(str(phil_path), working_directory=workdir)
            self.notebook.setCurrentIndex(1)
        except (RuntimeError, FileNotFoundError, ValueError, OSError) as exc:
            QMessageBox.warning(self, 'Xtrapol8', str(exc))

    def load_results(self):
        path = QFileDialog.getExistingDirectory(self, 'Select existing Xtrapol8 results directory', str(Path.home()))
        if not path:
            return
        directory = Path(path).resolve()
        config_path = self._find_results_phil(directory)
        config_message = ''
        if config_path is not None:
            if not self._confirm_replace_configuration():
                return
            try:
                phil_bridge.read_into(self, config_path)
                self._set_loaded_phil(config_path)
                self._results_source_directory = directory
                self._input_base_directory = None
                config_message = f'Configuration restored from {config_path.name}'
            except (RuntimeError, ValueError, OSError, AttributeError) as exc:
                config_message = f'Results loaded; PHIL could not be restored: {exc}'
                QMessageBox.warning(self, 'PHIL restoration failed', str(exc))
        else:
            config_message = 'Results loaded; Xtrapol8_in.phil not found (configuration unchanged)'
        name = self.tabIO.run_name.text().strip() or directory.name
        entry = self._add_run(name, directory)
        self.view_selected_results()
        entry['tabs'].setCurrentIndex(1)
        self.statusBar().showMessage(config_message, 12000)

    def stop_fextr(self):
        self.stop_selected_run()

    def append_log(self, line):
        if self.selected_run:
            self._run_log(self.selected_run, line)

    def process_finished(self, exit_code, stopped):
        if self.selected_run:
            self._run_finished(self.selected_run, exit_code, stopped)

    def process_error(self, message):
        if self.selected_run:
            self._run_error(self.selected_run, message)

    def closeEvent(self, event):
        if any(r['runner'] and r['runner'].running for r in self.runs):
            QMessageBox.warning(self, 'Xtrapol8', 'Stop all running calculations before closing.')
            event.ignore()
        else:
            for entry in self.runs:
                self._persist_run(entry)
            event.accept()

    def not_implemented(self, name):
        QMessageBox.information(self, 'Migration in progress', f'{name}: not yet migrated to Qt.')


def main():
    app = QApplication(sys.argv)
    icon_path = Path(__file__).resolve().parent / "assets" / "Xtrapol8_logo_dock.png"
    if icon_path.is_file():
        icon = QIcon(str(icon_path))
        app.setWindowIcon(icon)
    app.setApplicationName("Xtrapol8")
    app.setApplicationDisplayName("Xtrapol8")
    app.setOrganizationName("Xtrapol8")
    window = MainWindow()
    if icon_path.is_file():
        window.setWindowIcon(icon)
    window.show()
    return app.exec_()


if __name__ == '__main__':
    sys.exit(main())
