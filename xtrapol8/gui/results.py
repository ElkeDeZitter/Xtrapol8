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
    # e.g. qweight_occupancy_0.250 -> qFextr; other map types may be present.
    prefix = name.split('occupancy')[0].lower()
    if 'qweight' in prefix or prefix.startswith('q'):
        return 'qfextr'
    if 'kweight' in prefix or prefix.startswith('k'):
        return 'kfextr'
    if 'fextr' in prefix:
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


class ResultsViewer(QWidget):
    """Per-run viewer; scans only PNGs and never evaluates untrusted pickle data."""
    def __init__(self, directory=None, parent=None):
        super().__init__(parent)
        self.directory = Path(directory).expanduser().resolve() if directory else None
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
        occ_row = QHBoxLayout()
        occ_row.addWidget(QLabel('Occupation folder:'))
        self.selector = QComboBox()
        self.selector.currentIndexChanged.connect(self.show_occupation)
        self.selector.setMinimumWidth(190)
        self.selector.setMaximumWidth(340)
        occ_row.addWidget(self.selector, 0)
        self.best_label = QLabel('Best estimate: —')
        occ_row.addWidget(self.best_label)
        self.best_button = QPushButton('Go to best')
        self.best_button.setEnabled(False)
        self.best_button.clicked.connect(self.go_to_best)
        occ_row.addWidget(self.best_button)
        self.coot_button = QPushButton('Open in Coot')
        self.coot_button.setToolTip('Available at the estimated optimal occupancy when a Coot script exists')
        self.coot_button.setEnabled(False)
        self.coot_button.clicked.connect(self.open_coot)
        occ_row.addWidget(self.coot_button)
        occ_row.addStretch(1)
        occ_layout.addLayout(occ_row)
        self.occupancy = ImageGallery()
        occ_layout.addWidget(self.occupancy, 1)
        self.occupancies_page = occ
        # The parent run notebook displays Main and Occupancies as direct tabs.
        occ_actions = QHBoxLayout()
        choose_occ = QPushButton('Open results folder…')
        choose_occ.clicked.connect(self.choose_directory)
        occ_actions.addWidget(choose_occ)
        refresh_occ = QPushButton('Refresh')
        refresh_occ.clicked.connect(self.refresh)
        occ_actions.addWidget(refresh_occ)
        occ_actions.addStretch(1)
        occ_layout.insertLayout(0, occ_actions)
        self._groups = {}
        self._best = {}
        self._scripts = {}
        self.timer = QTimer(self)
        self.timer.setInterval(10000)
        self.timer.timeout.connect(self.refresh)
        self.timer.start()
        self.refresh()

    def choose_directory(self):
        path = QFileDialog.getExistingDirectory(self, 'Select Xtrapol8 results directory',
                                                str(self.directory or Path.home()))
        if path:
            self.directory = Path(path).resolve()
            self.refresh()

    def refresh(self):
        self.location.setText(f'Results: {self.directory}' if self.directory else 'Select a results directory')
        if not self.directory or not self.directory.is_dir():
            self.main.set_files([])
            self._groups = {}
            self._best = {}
            self._scripts = {}
            self.selector.clear()
            self._update_best()
            self.occupancy.set_files([])
            return
        # Limit depth and total files: avoid accidentally walking huge experimental trees.
        root_images = []
        groups = {}
        count = 0
        for path in self.directory.rglob('*.png'):
            if not path.is_file():
                continue
            relative = path.relative_to(self.directory)
            if len(relative.parts) > 5:
                continue
            if len(relative.parts) == 1:
                root_images.append(path)
            else:
                groups.setdefault(str(relative.parent), []).append(path)
            count += 1
            if count >= 500:
                break
        self.main.set_files(root_images)
        self._best, self._scripts = _read_best_occupancies(self.directory / 'occupancy_recap.pickle')
        previous = self.selector.currentText()
        self._groups = groups
        names = sorted(groups, key=str.lower)
        if [self.selector.itemText(i) for i in range(self.selector.count())] != names:
            self.selector.blockSignals(True)
            self.selector.clear()
            self.selector.addItems(names)
            idx = self.selector.findText(previous)
            self.selector.setCurrentIndex(idx if idx >= 0 else (0 if names else -1))
            self.selector.blockSignals(False)
        self.show_occupation()

    def _best_for_current(self):
        kind = _group_fextr(self.selector.currentText())
        if kind in self._best:
            return self._best[kind]
        if len(self._best) == 1 and len({k for k in map(_group_fextr, self._groups) if k}) <= 1:
            return next(iter(self._best.values()))
        return None

    def _script_for_current(self):
        kind = _group_fextr(self.selector.currentText())
        script = self._scripts.get(kind)
        if script is None and len(self._scripts) == 1 and len(
                {k for k in map(_group_fextr, self._groups) if k}) <= 1:
            script = next(iter(self._scripts.values()))
        if not script or self.directory is None:
            return None
        path = Path(script).expanduser()
        if not path.is_absolute():
            path = self.directory / path
        return path.resolve()

    def open_coot(self):
        script = self._script_for_current()
        if script is None or not script.is_file():
            QMessageBox.warning(self, 'Coot', 'The Coot script is unavailable.')
            self._update_best()
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

    def _update_best(self):
        best = self._best_for_current()
        self.best_label.setText(f'Best estimate: {best:.3f}' if best is not None else 'Best estimate: —')
        self.coot_button.setEnabled(best is not None and
            _group_occupancy(self.selector.currentText()) is not None and
            abs(_group_occupancy(self.selector.currentText()) - best) < 0.00051 and
            self._script_for_current() is not None and self._script_for_current().is_file())
        self.best_button.setEnabled(best is not None and any(
            _group_occupancy(name) is not None and
            abs(_group_occupancy(name) - best) < 0.00051 and
            (not _group_fextr(self.selector.currentText()) or
             _group_fextr(name) == _group_fextr(self.selector.currentText()))
            for name in self._groups))

    def go_to_best(self):
        best = self._best_for_current()
        if best is None:
            return
        current_kind = _group_fextr(self.selector.currentText())
        for index in range(self.selector.count()):
            name = self.selector.itemText(index)
            occ = _group_occupancy(name)
            if (occ is not None and abs(occ - best) < 0.00051 and
                    (not current_kind or _group_fextr(name) == current_kind)):
                self.selector.setCurrentIndex(index)
                return

    def show_occupation(self, *_):
        self.occupancy.set_files(self._groups.get(self.selector.currentText(), []))
        self._update_best()
