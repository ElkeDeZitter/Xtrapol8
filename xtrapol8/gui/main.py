"""Qt GUI migration milestone 1: window, notebook and Input/Output layout.

Run standalone: python -m xtrapol8.gui_qt.main
No scientific code is imported or modified.
"""
import sys
import os
import tempfile
import io
from pathlib import Path
from PyQt5.QtCore import Qt, QSize
from .action_icons import action_icon
from .results import ResultsViewer
from .progress import ProgressPanel
from .run_paths import resolve_inputs, unused_output_directory
from .process import FextrProcess
from . import phil_bridge
from .panels import ExtrapolationPanel, RefinementPanel
from PyQt5.QtGui import QIcon, QFontDatabase
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
    QHBoxLayout, QTabWidget, QGroupBox, QTableWidget, QTableWidgetItem,
    QPushButton, QLabel, QLineEdit, QFileDialog, QAction, QToolBar,
    QMessageBox, QHeaderView, QPlainTextEdit, QStyledItemDelegate, QComboBox, QSplitter, QMenu, QAbstractItemView)



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
        self.runner = FextrProcess(self)
        self.runner.log_received.connect(self.append_log)
        self.runner.started.connect(lambda: self.statusBar().showMessage('Fextr running'))
        self.runner.finished.connect(self.process_finished)
        self.runner.error.connect(self.process_error)
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
        self._run_logs = []
        self._run_viewers = []
        self._run_progress = []
        self.progress_panel = None
        # As in the wx GUI, run results are not shown until a run starts.
        self.run_results = None
        self.run_count = 0

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
        if self.runner.running:
            QMessageBox.information(self, 'Xtrapol8', 'A calculation is already running.')
            return
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
            # The generated PHIL must name the exact directory monitored by
            # the GUI. Never change the source PHIL on disk.
            self.tabIO.outdir.setText(str(target))
            run_dir = Path(tempfile.mkdtemp(prefix='xtrapol8_qt_'))
            self._run_phil_path = phil_bridge.save_from(self, run_dir / 'run.phil')
            self.show_run_results(track_progress=True)
            self._run_viewers[-1].directory = target
            self._run_viewers[-1].refresh()
            self.append_log(f'[Output directory] {target}')
            self.runner.start(str(self._run_phil_path), working_directory=workdir)
        except (RuntimeError, FileNotFoundError, ValueError, OSError) as exc:
            QMessageBox.warning(self, 'Xtrapol8', str(exc))

    def show_run_results(self, track_progress=False):
        # wx creates a Run #N page containing a Log subtab at run time.
        self.run_count += 1
        results = QTabWidget()
        log_view = QPlainTextEdit()
        log_view.setReadOnly(True)
        # Preserve alignment of fixed-width scientific tables in stdout/stderr.
        log_font = QFontDatabase.systemFont(QFontDatabase.FixedFont)
        log_font.setPointSize(11)
        log_view.setFont(log_font)
        log_view.setLineWrapMode(QPlainTextEdit.NoWrap)
        log_view.setTabStopDistance(log_view.fontMetrics().horizontalAdvance(' ') * 4)
        log_view.document().setMaximumBlockCount(20000)
        log_page = QWidget()
        log_layout = QHBoxLayout(log_page)
        log_layout.setContentsMargins(2, 2, 2, 2)
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(log_view)
        panel = None
        if track_progress:
            try:
                panel = ProgressPanel(phil_bridge.collect(self))
            except (RuntimeError, ValueError, AttributeError, TypeError) as exc:
                self.statusBar().showMessage(f'Progress unavailable: {exc}', 10000)
            if panel is not None:
                splitter.addWidget(panel)
                splitter.setSizes([820, 260])
                panel.setMinimumWidth(200)
        log_layout.addWidget(splitter)
        self.progress_panel = panel
        self._run_progress.append(panel)
        results.addTab(log_page, 'Log')
        viewer = ResultsViewer()
        results.addTab(viewer, 'Main')
        results.addTab(viewer.occupancies_page, 'Occupancies')
        self._run_viewers.append(viewer)
        self.run_results = results
        self.log_view = log_view
        self._run_logs.append(log_view)
        self.notebook.addTab(results, f'Run #{self.run_count}')
        self.notebook.setCurrentWidget(results)

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
                # Results remain useful even if the PHIL cannot be loaded.
                config_message = f'Results loaded; PHIL could not be restored: {exc}'
                QMessageBox.warning(self, 'PHIL restoration failed', str(exc))
        else:
            config_message = 'Results loaded; Xtrapol8_in.phil not found (configuration unchanged)'
        self.show_run_results()
        self._run_viewers[-1].directory = directory
        self._run_viewers[-1].refresh()
        self.run_results.setCurrentIndex(1)
        self.statusBar().showMessage(config_message, 12000)

    def stop_fextr(self):
        if self.runner.running:
            self.statusBar().showMessage('Stopping Xtrapol8…')
        self.runner.stop()

    def append_log(self, line):
        if self.log_view is not None:
            self.log_view.appendPlainText(line)
        if self.progress_panel is not None:
            self.progress_panel.ingest(line)

    def process_finished(self, exit_code, stopped):
        status = ('Stopped by user' if stopped else
                  'Completed successfully' if exit_code == 0 else
                  f'Process failed (exit code {exit_code})')
        self.statusBar().showMessage(status)
        self.append_log(f'[{status}]')
        if self.progress_panel is not None:
            self.progress_panel.finish(exit_code, stopped)
        if self.run_results is not None:
            self.notebook.setTabText(self.notebook.indexOf(self.run_results), f'Run #{self.run_count}')

    def process_error(self, message):
        self.append_log(f'[Process error] {message}')
        self.statusBar().showMessage(message)

    def closeEvent(self, event):
        if self.runner.running:
            QMessageBox.warning(self, 'Xtrapol8', 'Stop the running calculation before closing.')
            event.ignore()
        else:
            event.accept()

    def not_implemented(self, name):
        QMessageBox.information(self, 'Migration in progress',
                                f'{name}: not yet migrated to Qt.')


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
