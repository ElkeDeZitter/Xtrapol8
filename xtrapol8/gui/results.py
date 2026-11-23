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


class _MethodPage(QWidget):
    """Independent occupancy selection, gallery and Coot action for one Fextr type."""

    def __init__(self, viewer, kind):
        super().__init__()
        self.viewer = viewer
        self.kind = kind
        self.folders = {}
        layout = QVBoxLayout(self)
        row = QHBoxLayout()
        row.addWidget(QLabel('Occupancy:'))
        self.selector = QComboBox()
        self.selector.setMinimumWidth(110)
        self.selector.setMaximumWidth(190)
        self.selector.currentIndexChanged.connect(self.show_occupation)
        row.addWidget(self.selector)
        self.best_label = QLabel('Best occupancy / α: —')
        row.addWidget(self.best_label)
        self.best_button = QPushButton('Go to best')
        self.best_button.clicked.connect(self.go_to_best)
        row.addWidget(self.best_button)
        self.coot_button = QPushButton('Open in Coot')
        self.coot_button.setToolTip('Available only at the optimal occupancy with an existing Coot script')
        self.coot_button.clicked.connect(self.open_coot)
        row.addWidget(self.coot_button)
        row.addStretch(1)
        layout.addLayout(row)
        self.gallery = ImageGallery()
        layout.addWidget(self.gallery, 1)
        self.show_occupation()

    def update_folders(self, folders):
        previous = self.selector.currentText()
        self.folders = folders
        names = sorted(folders, key=lambda name: (_group_occupancy(name) or 0, name))
        if [self.selector.itemText(i) for i in range(self.selector.count())] != names:
            self.selector.blockSignals(True)
            self.selector.clear()
            self.selector.addItems(names)
            index = self.selector.findText(previous)
            self.selector.setCurrentIndex(index if index >= 0 else (0 if names else -1))
            self.selector.blockSignals(False)
        self.show_occupation()

    def best(self):
        return self.viewer._best.get(self.kind)

    def script(self):
        raw = self.viewer._scripts.get(self.kind)
        if not raw or self.viewer.directory is None:
            return None
        path = Path(raw).expanduser()
        if not path.is_absolute():
            path = self.viewer.directory / path
        return path.resolve()

    def show_occupation(self, *_):
        self.gallery.set_files(self.folders.get(self.selector.currentText(), []))
        best = self.best()
        self.best_label.setText(
            f'Best occupancy: {best:.3f}  |  α: {1 / best:.2f}'
            if best is not None and best > 0 else 'Best occupancy / α: —'
        )
        available = best is not None and any(
            _group_occupancy(name) is not None and
            abs(_group_occupancy(name) - best) < 0.00051
            for name in self.folders
        )
        self.best_button.setEnabled(available)
        selected = _group_occupancy(self.selector.currentText())
        script = self.script()
        self.coot_button.setEnabled(
            available and selected is not None and
            abs(selected - best) < 0.00051 and
            script is not None and script.is_file()
        )

    def go_to_best(self):
        best = self.best()
        if best is None:
            return
        for index in range(self.selector.count()):
            occ = _group_occupancy(self.selector.itemText(index))
            if occ is not None and abs(occ - best) < 0.00051:
                self.selector.setCurrentIndex(index)
                return

    def open_coot(self):
        if not self.coot_button.isEnabled():
            return
        self.viewer._launch_coot(self.script())


class ResultsViewer(QWidget):
    """Per-run PNG viewer; separate method tabs, each with independent state."""
    METHODS = (('fextr', 'Fextr'), ('qfextr', 'qFextr'), ('kfextr', 'kFextr'))

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
        actions = QHBoxLayout()
        choose_occ = QPushButton('Open results folder…')
        choose_occ.clicked.connect(self.choose_directory)
        actions.addWidget(choose_occ)
        refresh_occ = QPushButton('Refresh')
        refresh_occ.clicked.connect(self.refresh)
        actions.addWidget(refresh_occ)
        actions.addStretch(1)
        occ_layout.addLayout(actions)
        self.method_tabs = QTabWidget()
        occ_layout.addWidget(self.method_tabs, 1)
        # Initialize shared results before _MethodPage.__init__ calls show_occupation().
        self._best = {}
        self._scripts = {}
        self._groups = {}
        self.pages = {kind: _MethodPage(self, kind) for kind, _ in self.METHODS}
        self.occupancies_page = occ
        self.timer = QTimer(self)
        self.timer.setInterval(10000)
        self.timer.timeout.connect(self.refresh)
        self.timer.start()
        self.refresh()

    @property
    def current_page(self):
        return self.method_tabs.currentWidget()

    @property
    def selector(self):
        return self.current_page.selector if self.current_page is not None else None

    @property
    def coot_button(self):
        return self.current_page.coot_button if self.current_page is not None else None

    def _script_for_current(self):
        return self.current_page.script() if self.current_page is not None else None

    def _best_for_current(self):
        return self.current_page.best() if self.current_page is not None else None

    def go_to_best(self):
        if self.current_page is not None:
            self.current_page.go_to_best()

    def open_coot(self):
        if self.current_page is not None:
            self.current_page.open_coot()

    def choose_directory(self):
        path = QFileDialog.getExistingDirectory(self, 'Select Xtrapol8 results directory',
                                                str(self.directory or Path.home()))
        if path:
            self.directory = Path(path).resolve()
            self.refresh()

    def refresh(self):
        self.location.setText(f'Results: {self.directory}' if self.directory else 'Select a results directory')
        root_images = []
        groups = {kind: {} for kind, _ in self.METHODS}
        self._best, self._scripts = {}, {}
        if self.directory and self.directory.is_dir():
            # Discover actual occupancy folders, even if their PNGs are not ready yet.
            for folder in self.directory.iterdir():
                if folder.is_dir():
                    kind = _group_fextr(folder.name)
                    if kind:
                        groups[kind][folder.name] = []
            count = 0
            for path in self.directory.rglob('*.png'):
                if not path.is_file():
                    continue
                relative = path.relative_to(self.directory)
                if len(relative.parts) > 5:
                    continue
                if len(relative.parts) == 1:
                    root_images.append(path)
                elif relative.parts[0] in groups.get(_group_fextr(relative.parts[0]), {}):
                    kind = _group_fextr(relative.parts[0])
                    groups[kind][relative.parts[0]].append(path)
                count += 1
                if count >= 500:
                    break
            self._best, self._scripts = _read_best_occupancies(
                self.directory / 'occupancy_recap.pickle'
            )
        self.main.set_files(root_images)
        self._groups = groups
        active = self.current_page.kind if self.current_page is not None else None
        available = [kind for kind, _ in self.METHODS if groups[kind]]
        current = [self.method_tabs.widget(i) for i in range(self.method_tabs.count())]
        expected = [self.pages[kind] for kind in available]
        if current != expected:
            self.method_tabs.clear()
            for kind, label in self.METHODS:
                if kind in available:
                    self.method_tabs.addTab(self.pages[kind], label)
            if active in available:
                self.method_tabs.setCurrentWidget(self.pages[active])
        for kind, page in self.pages.items():
            page.update_folders(groups[kind])

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
