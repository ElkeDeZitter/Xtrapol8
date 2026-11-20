"""Read-only viewer of Xtrapol8 PNG results; no pickle deserialization."""
from pathlib import Path
import io
import pickle
import re
import shutil
from PyQt5.QtCore import Qt, QTimer, QProcess
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QComboBox, QScrollArea, QFileDialog, QTabWidget, QFrame, QMessageBox)


class _DataOnlyUnpickler(pickle.Unpickler):
    """Reject arbitrary Python objects; occupancy recap uses plain dict/list/scalars."""
    def find_class(self, module, name):
        raise pickle.UnpicklingError('Python objects are not allowed in occupancy recap')


def _read_best_occupancies(path):
    if not path.is_file() or path.stat().st_size > 1_000_000:
        return {}, {}
    try:
        result = _DataOnlyUnpickler(io.BytesIO(path.read_bytes())).load()
        if not isinstance(result, dict):
            return {}, {}
        best = {}
        scripts = {}
        for kind, entry in result.items():
            if not isinstance(kind, str) or not isinstance(entry, (list, tuple)) or not entry:
                continue
            value = float(entry[0])
            if 0 < value <= 1:
                key = kind.lower()
                best[key] = value
                if len(entry) > 1 and isinstance(entry[1], (str, bytes)):
                    script = entry[1].decode("utf-8", "replace") if isinstance(entry[1], bytes) else entry[1]
                    if script.endswith(".py"):
                        scripts[key] = script
        return best, scripts
    except (OSError, ValueError, TypeError, pickle.UnpicklingError, OverflowError):
        return {}, {}


def _group_occupancy(name):
    match = re.search(r'occupancy[_-]([0-9]+(?:\.[0-9]+)?)', name, re.I)
    return float(match.group(1)) if match else None


def _group_fextr(name):
    # Only the occupancy directory itself determines the extrapolation type.
    # In particular, plain occupancy_0.200 is unweighted Fextr.
    folder = Path(name).name.lower()
    if re.fullmatch(r'qweight_occupancy_\d+(?:\.\d+)?', folder):
        return 'qfextr'
    if re.fullmatch(r'kweight_occupancy_\d+(?:\.\d+)?', folder):
        return 'kfextr'
    if re.fullmatch(r'occupancy_\d+(?:\.\d+)?', folder):
        return 'fextr'
    return None


class ImageGallery(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.root = QVBoxLayout(self)
        self.note = QLabel('No PNG results found yet.')
        self.root.addWidget(self.note)
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.root.addWidget(self.scroll, 1)
        self.images = []
        self._paths = []
        self.scroll.viewport().installEventFilter(self)

    def eventFilter(self, obj, event):
        if obj is self.scroll.viewport() and event.type() == event.Resize:
            self._resize_images()
        return super().eventFilter(obj, event)

    def _resize_images(self):
        # One full-width graph per row, preserving the source aspect ratio.
        available = max(100, self.scroll.viewport().width() - 64)
        for label, path in self.images:
            pix = QPixmap(str(path))
            if not pix.isNull():
                width = min(available, pix.width())
                label.setPixmap(pix.scaledToWidth(width, Qt.SmoothTransformation))

    def set_files(self, paths):
        paths = sorted(paths, key=lambda p: p.name.lower())
        if paths == self._paths:
            return
        self._paths = paths
        self.images = []
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setAlignment(Qt.AlignTop)
        for i, path in enumerate(paths):
            frame = QFrame()
            frame.setFrameShape(QFrame.StyledPanel)
            box = QVBoxLayout(frame)
            title = QLabel(path.name)
            title.setTextInteractionFlags(Qt.TextSelectableByMouse)
            box.addWidget(title)
            pic = QLabel()
            pic.setAlignment(Qt.AlignCenter)
            pix = QPixmap(str(path))
            if pix.isNull():
                pic.setText('Cannot read image')
            box.addWidget(pic)
            layout.addWidget(frame)
            self.images.append((pic, path))
        self.scroll.takeWidget()
        self.scroll.setWidget(container)
        self._resize_images()
        self.note.setText(f'{len(paths)} graph(s)' if paths else 'No PNG results found yet.')


# The wxPython GUI supports all three families with three weighting choices.
FAMILIES = ('Fextr', 'Fgenick', 'Fextr_calc')
WEIGHTS = ('', 'q', 'k')
METHODS = tuple((weight + family).lower() for family in FAMILIES for weight in WEIGHTS)


def _method_parts(method):
    method = method.lower()
    weight = method[0] if method.startswith(('q', 'k')) else ''
    return method[len(weight):], weight


def _method_folder(method, occupancy):
    _, weight = _method_parts(method)
    return f"{('qweight_' if weight == 'q' else 'kweight_' if weight == 'k' else '')}occupancy_{occupancy:.3f}"


def _method_from_filename(filename):
    # Exact method token, not a prefix: Fextr must never match Fextr_calc.
    match = re.match(r'^(q|k)?(Fextr_calc|Fgenick|Fextr)(?=_|\.)', filename, re.I)
    return match.group(0).lower() if match else None


class ResultsViewer(QWidget):
    """PNG results organized by family, weighting and occupancy.

    A method is selectable only if it has method-specific output files or a
    completed recap entry. An occupancy folder by itself is not evidence.
    """
    def __init__(self, directory=None, parent=None):
        super().__init__(parent)
        self.directory = Path(directory).expanduser().resolve() if directory else None
        self._best, self._scripts, self._groups = {}, {}, {}
        self._selected_occupancy = {}
        self._available = set()
        layout = QVBoxLayout(self)
        top = QHBoxLayout()
        self.location = QLabel()
        self.location.setTextInteractionFlags(Qt.TextSelectableByMouse)
        top.addWidget(self.location, 1)
        choose = QPushButton('Open results folder…')
        choose.clicked.connect(self.choose_directory)
        top.addWidget(choose)
        refresh = QPushButton('Refresh')
        refresh.clicked.connect(self.refresh)
        top.addWidget(refresh)
        layout.addLayout(top)
        self.main = ImageGallery()
        layout.addWidget(self.main, 1)

        occ = QWidget()
        occ_layout = QVBoxLayout(occ)
        actions = QHBoxLayout()
        choose_occ = QPushButton('Open results folder…')
        choose_occ.clicked.connect(self.choose_directory)
        actions.addWidget(choose_occ)
        refresh_occ = QPushButton('Refresh')
        refresh_occ.clicked.connect(self.refresh)
        actions.addWidget(refresh_occ)
        actions.addStretch(1)
        occ_layout.addLayout(actions)
        method_row = QHBoxLayout()
        method_row.addWidget(QLabel('Family:'))
        self.family_selector = QComboBox()
        self.family_selector.currentIndexChanged.connect(self._family_changed)
        method_row.addWidget(self.family_selector)
        method_row.addWidget(QLabel('Weighting:'))
        self.weight_selector = QComboBox()
        self.weight_selector.currentIndexChanged.connect(self._weight_changed)
        method_row.addWidget(self.weight_selector)
        method_row.addStretch(1)
        occ_layout.addLayout(method_row)
        row = QHBoxLayout()
        row.addWidget(QLabel('Occupancy:'))
        self.selector = QComboBox()
        self.selector.currentIndexChanged.connect(self.show_occupation)
        row.addWidget(self.selector)
        self.best_label = QLabel('Best occupancy / α: —')
        row.addWidget(self.best_label)
        self.best_button = QPushButton('Go to best')
        self.best_button.clicked.connect(self.go_to_best)
        row.addWidget(self.best_button)
        self.coot_button = QPushButton('Open in Coot')
        self.coot_button.clicked.connect(self.open_coot)
        row.addWidget(self.coot_button)
        row.addStretch(1)
        occ_layout.addLayout(row)
        self.occupancy = ImageGallery()
        occ_layout.addWidget(self.occupancy, 1)
        self.occupancies_page = occ
        self.timer = QTimer(self)
        self.timer.setInterval(10000)
        self.timer.timeout.connect(self.refresh)
        self.timer.start()
        self.refresh()

    @property
    def current_method(self):
        family = self.family_selector.currentData()
        weight = self.weight_selector.currentData()
        return (weight + family).lower() if family is not None and weight is not None else None

    def _set_options(self, combo, options, selected):
        previous = [(combo.itemText(i), combo.itemData(i)) for i in range(combo.count())]
        if previous == options:
            return
        combo.blockSignals(True)
        combo.clear()
        for label, value in options:
            combo.addItem(label, value)
        idx = combo.findData(selected)
        combo.setCurrentIndex(idx if idx >= 0 else (0 if options else -1))
        combo.blockSignals(False)

    def _family_changed(self, *_):
        self._refresh_weights()
        self._refresh_occupancies()

    def _weight_changed(self, *_):
        self._refresh_occupancies()

    def _refresh_weights(self, preferred=None):
        family = self.family_selector.currentData()
        if preferred is None:
            preferred = self.weight_selector.currentData()
        labels = {'': 'None', 'q': 'q-weighted', 'k': 'k-weighted'}
        options = [(labels[w], w) for w in WEIGHTS if family and (w + family).lower() in self._available]
        self._set_options(self.weight_selector, options, preferred)

    def _refresh_occupancies(self):
        method = self.current_method
        previous = self._selected_occupancy.get(method) if method else None
        options = [(f'{occ:.3f}', occ) for occ in sorted(self._groups.get(method, {}))]
        self._set_options(self.selector, options, previous)
        self.show_occupation()

    def choose_directory(self):
        path = QFileDialog.getExistingDirectory(self, 'Select Xtrapol8 results directory',
                                                str(self.directory or Path.home()))
        if path:
            self.directory = Path(path).resolve()
            self._selected_occupancy.clear()
            self.refresh()

    def refresh(self):
        self.location.setText(f'Results: {self.directory}' if self.directory else 'Select a results directory')
        previous_method = self.current_method
        if previous_method is not None and self.selector.currentData() is not None:
            self._selected_occupancy[previous_method] = self.selector.currentData()
        root_images = []
        groups = {method: {} for method in METHODS}
        self._best, self._scripts = {}, {}
        if self.directory and self.directory.is_dir():
            self._best, self._scripts = _read_best_occupancies(self.directory / 'occupancy_recap.pickle')
            count = 0
            # Do not infer availability from directory names: several methods share them.
            for folder in self.directory.iterdir():
                if not folder.is_dir():
                    continue
                occupancy = _group_occupancy(folder.name)
                if occupancy is None or not re.fullmatch(r'(?:qweight_|kweight_)?occupancy_\d+(?:\.\d+)?', folder.name, re.I):
                    continue
                for path in folder.rglob('*'):
                    if not path.is_file() or len(path.relative_to(folder).parts) > 4:
                        continue
                    method = _method_from_filename(path.name)
                    if method and method in groups and folder.name.lower() == _method_folder(method, occupancy).lower():
                        groups[method].setdefault(occupancy, [])
                        if path.suffix.lower() == '.png':
                            groups[method][occupancy].append(path)
                    count += 1
                    if count >= 3000:
                        break
                if count >= 3000:
                    break
            # A completed recap is independent evidence of computation. Populate
            # only occupancy folders belonging to that specific weighting.
            for method in self._best:
                if method not in groups:
                    continue
                for folder in self.directory.iterdir():
                    if folder.is_dir():
                        occ = _group_occupancy(folder.name)
                        if occ is not None and folder.name.lower() == _method_folder(method, occ).lower():
                            groups[method].setdefault(occ, [])
            root_images = sorted(self.directory.glob('*.png'))[:500]
        self.main.set_files(root_images)
        self._groups = groups
        self._available = {method for method, occs in groups.items() if occs}
        families = [(family, family) for family in FAMILIES if any((w + family).lower() in self._available for w in WEIGHTS)]
        old_family, old_weight = _method_parts(previous_method) if previous_method else (None, None)
        self._set_options(self.family_selector, families, next((f for f in FAMILIES if f.lower() == old_family), None))
        self._refresh_weights(old_weight)
        self._refresh_occupancies()

    def show_occupation(self, *_):
        method = self.current_method
        occupancy = self.selector.currentData()
        if method and occupancy is not None:
            self._selected_occupancy[method] = occupancy
        self.occupancy.set_files(self._groups.get(method, {}).get(occupancy, []))
        best = self._best.get(method)
        self.best_label.setText(f'Best occupancy: {best:.3f}  |  α: {1/best:.2f}'
                                if best is not None and best > 0 else 'Best occupancy / α: —')
        available = best is not None and any(abs(occ - best) < 0.00051 for occ in self._groups.get(method, {}))
        self.best_button.setEnabled(available)
        script = self._script_for_current()
        self.coot_button.setEnabled(bool(available and occupancy is not None and
                                         abs(occupancy - best) < 0.00051 and
                                         script is not None and script.is_file()))

    def _best_for_current(self):
        return self._best.get(self.current_method)

    def _script_for_current(self):
        raw = self._scripts.get(self.current_method)
        if not raw or self.directory is None:
            return None
        path = Path(raw).expanduser()
        return (path if path.is_absolute() else self.directory / path).resolve()

    def go_to_best(self):
        best = self._best_for_current()
        if best is None:
            return
        for i in range(self.selector.count()):
            if abs(self.selector.itemData(i) - best) < 0.00051:
                self.selector.setCurrentIndex(i)
                break

    def open_coot(self):
        if self.coot_button.isEnabled():
            self._launch_coot(self._script_for_current())

    def _launch_coot(self, script):
        if script is None or not script.is_file():
            QMessageBox.warning(self, 'Coot', 'The Coot script is unavailable.')
            return
        executable = next((found for name in ('coot', 'coot1', 'coot-1', 'coot0.9')
                           if (found := shutil.which(name))), None)
        if not executable:
            QMessageBox.warning(self, 'Coot not found',
                                'Could not find Coot on PATH. Launch the Qt GUI from a shell '
                                'where the coot command is available.')
            return
        started, _pid = QProcess.startDetached(executable, ['--script', str(script)],
                                               str(self.directory))
        if not started:
            QMessageBox.warning(self, 'Coot', f'Failed to launch {executable}')
