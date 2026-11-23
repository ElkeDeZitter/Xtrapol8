"""Workflow-aware, log-driven Xtrapol8 progress tracking (no engine changes)."""
import re
from pathlib import Path
import shutil
from PyQt5.QtCore import QProcess
from PyQt5.QtWidgets import QPushButton, QMessageBox
from PyQt5.QtCore import QTimer, Qt
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QProgressBar, QHBoxLayout
from PyQt5.QtGui import QPainter, QPen, QColor

NUMBER = r'(?:[0-9]+(?:\.[0-9]+)?)'
START = re.compile(r'CALCULATING\s+(\S+)\s+TYPE OF ESFAS AND MAPS FOR OCCUPANCY\s+(' + NUMBER + r')', re.I)
RESULT = re.compile(r'--->\s*Results in\s+(.+?)\s*$', re.I)
OCC_DIR = re.compile(r'occupancy_(' + NUMBER + r')', re.I)
REFINE_WITH = re.compile(r'RECIPROCAL\s+SPACE\s+REFINEMENT\s+WITH\s+\S*?occ(' + NUMBER + r')_([a-z][a-z0-9_]*)\.mtz', re.I)
REFINE_START = re.compile(r'\b(?:phenix\.refine|refmac5|refmac|starting\s+(?:reciprocal.space\s+)?refinement|refinement\s+(?:at|for|of)\s+occupancy)\b', re.I)
REFINE_DONE = re.compile(r'\b(?:refinement\s+(?:finished|completed|successful)|refinement\s+done|refinement\s+results\s+in)\b', re.I)
OCC_IN_TEXT = re.compile(r'(?:occupancy[_\s=]+|\bocc(?:upancy)?\s*)([0-9]+(?:\.[0-9]+)?)', re.I)
ALPHA = re.compile(r'\balpha\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)', re.I)
BEST_OCC = re.compile(r'(?:optimal|best|estimated)\s+occupancy\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)', re.I)


def _norm(occ):
    return round(float(occ), 6)


def expected_calculations(params):
    fofo = bool(params.output.generate_fofo_only)
    fast = bool(params.f_and_maps.fast_and_furious)
    if fofo:
        return 0, 'FoFo only'
    occ = params.occupancies.list_occ
    n_occ = len(occ) if occ is not None and len(occ) else int(params.occupancies.steps) + 1
    n_types = 1 if fast else len(params.f_and_maps.f_extrapolated_and_maps or [])
    return n_occ * n_types, "Fast'N'Furious" if fast else "Calm'N'Curious"


class SmoothSpinner(QWidget):
    """Antialiased rotating arc, independent of font glyphs on macOS."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedSize(18, 18)
        self.angle = 0

    def advance(self):
        self.angle = (self.angle + 24) % 360
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        pen = QPen(QColor('#2779ad'), 2.6, Qt.SolidLine, Qt.RoundCap)
        painter.setPen(pen)
        painter.drawArc(3, 3, 12, 12, -self.angle * 16, 260 * 16)
        painter.end()


class ProgressPanel(QWidget):
    def __init__(self, params, parent=None):
        super().__init__(parent)
        self.total, self.mode = expected_calculations(params)
        self.fast = bool(params.f_and_maps.fast_and_furious)
        self.fofo = bool(params.output.generate_fofo_only)
        self.run_refinement = bool(params.refinement.run_refinement)
        self.output_directory = Path(str(getattr(params.output, 'outdir', '.'))).expanduser()
        self._fofo_analysis_done = False
        self.refine_total = 0 if (self.fofo or not self.run_refinement) else (1 if self.fast else self.total)
        self.phase = 'waiting'
        self.active = None
        self.fextr_started = set()
        self.fextr_completed = set()
        self.refine_started = set()
        self.refine_completed = set()
        self.refine_active = None
        self.refine_context = None  # (map type, occupancy) currently eligible for refinement
        self.last_fextr_completed = None
        self.best_occ = None
        self.best_alpha = None
        self.final_status = None
        self.alpha_done = False
        self._spinner_frame = 0
        self._spinner_timer = QTimer(self)
        self._spinner_timer.setInterval(120)
        self._spinner_timer.timeout.connect(self._advance_spinner)
        self._last_detail = ''
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 6, 10, 6)
        heading = QLabel('Processing steps')
        heading.setStyleSheet('font-weight: 600;')
        layout.addWidget(heading)
        layout.addWidget(QLabel(self.mode))
        self.steps = {}
        names = [('preparation', 'Data preparation'), ('fofo', 'FoFo map calculation'),
                 ('analysis', 'FoFo map analysis')]
        if not self.fofo:
            names += [('fextr', 'Fextr calculations')]
            if self.fast:
                names += [('alpha', 'Optimal occupancy / α')]
            if self.run_refinement:
                names += [('refine', 'Refinements')]
        for key, name in names:
            label = QLabel('○  ' + name)
            label.setWordWrap(True)
            layout.addWidget(label)
            self.steps[key] = (label, name)
        self.counter = QLabel()
        self.bar = QProgressBar()
        self.bar.setRange(0, max(1, self.total))
        self.refine_counter = QLabel()
        self.refine_bar = QProgressBar()
        self.refine_bar.setRange(0, max(1, self.refine_total))
        self.detail = QLabel('Waiting to start')
        self.detail.setWordWrap(True)
        self.best_label = QLabel('')
        self.best_label.setWordWrap(True)
        if not self.fofo:
            widgets = [self.counter, self.bar]
            if self.run_refinement:
                widgets += [self.refine_counter, self.refine_bar]
            widgets += [self.best_label, self.detail]
            for widget in widgets:
                layout.addWidget(widget)
        self.fofo_coot_button = None
        if self.fofo:
            self.fofo_coot_button = QPushButton('Open in Coot')
            self.fofo_coot_button.setToolTip('Open the completed FoFo map in Coot')
            self.fofo_coot_button.setEnabled(False)
            self.fofo_coot_button.hide()
            self.fofo_coot_button.clicked.connect(self._open_fofo_coot)
            layout.addWidget(self.fofo_coot_button)
        status_row = QHBoxLayout()
        self.spinner = SmoothSpinner(self)
        self.spinner.hide()
        self.status_icon = QLabel('○')
        self.status_icon.setFixedWidth(20)
        self.status = QLabel('Waiting to start')
        status_row.addWidget(self.spinner)
        status_row.addWidget(self.status_icon)
        status_row.addWidget(self.status)
        status_row.addStretch(1)
        layout.addLayout(status_row)
        layout.addStretch(1)
        self._render()

    def _advance_spinner(self):
        self.spinner.advance()

    def ingest(self, line):
        upper = line.upper()
        if self.final_status is None and not self._spinner_timer.isActive():
            self._spinner_timer.start()
        if 'DATA PREPARATION' in upper:
            self.phase = 'preparation'
        elif 'CALCULATE (WEIGHTED) FO-FO FOURIER DIFFERENCE MAP' in upper:
            self.phase = 'fofo'
        elif 'MAP EXPLORER' in upper and self.phase in ('preparation', 'fofo'):
            self.phase = 'analysis'
        if self.fofo and ('FOFO MAP EXPLORED' in upper or 'FO-FO MAP EXPLORED' in upper or 'FOFO MAP ANALYSIS DONE' in upper):
            self._fofo_analysis_done = True
        elif 'CALCULATE (Q/K-WEIGHTED) FEXTRAPOLATED MAPS' in upper:
            self.phase = 'fextr'

        refine_with = REFINE_WITH.search(line) if not self.fast else None
        match = START.search(line)
        if match:
            # The next calculation is reached only after the preceding
            # refinement has returned. This is a boundary, not a guess based
            # on arbitrary phenix output.
            if not self.fast and self.refine_active is not None:
                self.refine_completed.add(self.refine_active)
                self.refine_active = None
                self.refine_context = None
            key = (match.group(1).lower(), _norm(match.group(2)))
            self.active = key
            self.fextr_started.add(key)
            self.phase = 'fextr'
            if not self.fast:
                self.refine_context = None
                self.refine_active = None
            self._last_detail = f'Fextr running: {match.group(1)} — occupancy {key[1]:g}'

        result = RESULT.search(line)
        if result:
            folder = Path(result.group(1).strip().rstrip('/')).name
            occ_match = OCC_DIR.search(folder)
            if occ_match:
                occ = _norm(occ_match.group(1))
                # The directory name identifies the occupancy; a completion is
                # counted only when the engine explicitly writes 'Results in'.
                candidates = [key for key in self.fextr_started if key[1] == occ]
                if self.active is not None and self.active[1] == occ:
                    key = self.active
                elif len(candidates) == 1:
                    key = candidates[0]
                else:
                    key = None  # never guess which map type completed
                if key is not None:
                    self.fextr_completed.add(key)
                    self.last_fextr_completed = key
                    if not self.fast:
                        self.refine_context = key
                    self._last_detail = f'Fextr completed: {key[0]} — occupancy {occ:g}'
                    if self.active == key:
                        self.active = None

        if refine_with is not None:
            # Xtrapol8's real Calm'N'Curious output does not always print
            # '---> Results in' before entering refinement. The MTZ name
            # explicitly identifies both the map type and occupancy.
            occ = _norm(refine_with.group(1))
            map_type = refine_with.group(2).lower()
            key = (map_type, occ)
            if key in self.fextr_started:
                self.fextr_completed.add(key)
                self.last_fextr_completed = key
                self.refine_context = key
                if self.active == key:
                    self.active = None
                self.refine_active = key
                self.refine_started.add(key)
                self.phase = 'refine'
                self._last_detail = f'Refinement running: {map_type} — occupancy {occ:g}'

        if not self.fofo:
            if self.fast:
                if 'ESTIMATE OPTIMAL OCCUPANCY DONE' in upper:
                    self.alpha_done = True
                if ('OPTIMAL OCCUPANCY' in upper or 'BEST OCCUPANCY' in upper or
                        'ESTIMAT' in upper and 'OCCUPANCY' in upper):
                    self.phase = 'alpha'
                selected = re.search(r'->\s*Optimal occupancy of triggered state\s+(' + NUMBER + r')', line, re.I)
                if selected:
                    self.best_occ = _norm(selected.group(1))
                    self.best_alpha = 1 / self.best_occ if self.best_occ else None
                    self.alpha_done = True
                best = BEST_OCC.search(line)
                if best:
                    self.best_occ = _norm(best.group(1))
                    self.best_alpha = 1 / self.best_occ if self.best_occ else None
                alpha = ALPHA.search(line)
                if alpha and self.phase == 'alpha':
                    self.best_alpha = float(alpha.group(1))
            if self.run_refinement and self.fast and 'FAST AND FURIOUS REFINEMENT' in upper:
                self.alpha_done = True
                self.phase = 'refine'
                key = ('best', self.best_occ)
                self.refine_active = key
                self.refine_started.add(key)
                self._last_detail = 'Refinement running at best occupancy'
            if self.run_refinement and REFINE_START.search(line) and not refine_with and (not self.fast or self.alpha_done):
                # For Calm'N'Curious, refinement is repeated per occupancy;
                # for Fast'N'Furious, it is a single best-occupancy refinement.
                occ_match = OCC_IN_TEXT.search(line)
                occ = _norm(occ_match.group(1)) if occ_match else None
                if self.fast:
                    key = ('best', self.best_occ)
                else:
                    # A refinement belongs to the *completed* map calculation,
                    # not the next occupancy nor another map type sharing its folder.
                    context = self.refine_context
                    if context is not None and (occ is None or occ == context[1]):
                        key = context
                        occ = context[1]
                    else:
                        key = None
                if key is not None:
                    self.refine_active = key
                    self.refine_started.add(key)
                    self._last_detail = ('Refinement running at best occupancy' if self.fast else
                                         f'Refinement running: {key[0]} — occupancy {occ:g}')
                if key is not None:
                    self.phase = 'refine'
            if self.run_refinement and ('FAST AND FURIOUS REFINEMENT DONE' in upper or REFINE_DONE.search(line)) and self.refine_active is not None:
                self.refine_completed.add(self.refine_active)
                self._last_detail = 'Refinement completed'
                self.refine_active = None
                if not self.fast:
                    self.refine_context = None
        self._render()

    def finish(self, exit_code, stopped):
        self._spinner_timer.stop()
        self.final_status = ('Stopped by user' if stopped else 'Completed' if exit_code == 0
                             else f'Failed (exit {exit_code})')
        if exit_code == 0 and not stopped and self.fast and self.best_occ is not None:
            self.alpha_done = True
        # Successful process exit confirms the final refinement, even if the
        # scientific program omitted an explicit refinement completion line.
        if self.run_refinement and exit_code == 0 and not stopped and self.refine_active is not None:
            self.refine_completed.add(self.refine_active)
            self.refine_active = None
        if self.fofo and exit_code == 0 and not stopped:
            self._fofo_analysis_done = True
        self._render()

    def percent_complete(self):
        if self.final_status == 'Completed':
            return 100
        if self.fofo:
            return {'waiting': 0, 'preparation': 15, 'fofo': 45, 'analysis': 75}.get(self.phase, 90)
        initial = {'waiting': 0, 'preparation': 5, 'fofo': 10, 'analysis': 15}.get(self.phase, 20)
        fextr = 65 * min(len(self.fextr_completed), self.total) / max(1, self.total)
        tail = (10 if self.fast and self.alpha_done else 0)
        if self.run_refinement and self.refine_total:
            tail += 5 * min(len(self.refine_completed), self.refine_total) / self.refine_total
        return min(99, round(initial + fextr + tail))

    def _fofo_coot_script(self):
        # Original wxPython panelLog.py checks these in this exact order.
        for prefix in ('q', 'k', ''):
            candidate = self.output_directory / f'coot_all_{prefix}Fo-Fo.py'
            if candidate.is_file():
                return candidate
        return None

    def _open_fofo_coot(self):
        script = self._fofo_coot_script()
        if self.final_status != 'Completed' or script is None:
            QMessageBox.warning(self, 'Coot', 'FoFo Coot script is not available yet.')
            return
        executable = next((path for name in ('coot', 'coot1', 'coot-1', 'coot0.9')
                           if (path := shutil.which(name))), None)
        if executable is None:
            QMessageBox.warning(self, 'Coot', 'Coot executable not found on PATH.')
            return
        started, _pid = QProcess.startDetached(executable, ['--script', str(script)], str(self.output_directory))
        if not started:
            QMessageBox.warning(self, 'Coot', 'Could not launch Coot.')

    def _render(self):
        completed = min(len(self.fextr_completed), self.total)
        self.counter.setText(f'Fextr completed: {completed} / {self.total}')
        self.bar.setValue(completed)
        ref_done = min(len(self.refine_completed), self.refine_total)
        self.refine_counter.setText(f'Refinements completed: {ref_done} / {self.refine_total}')
        self.refine_bar.setValue(ref_done)
        if self.best_occ is not None:
            self.best_label.setText(f'Best occupancy: {self.best_occ:g}   α: {self.best_alpha:.3f}' if self.best_alpha else
                                    f'Best occupancy: {self.best_occ:g}')
        if not self.fofo:
            self.detail.setText(self._last_detail or 'Waiting for Fextr calculations')
        ordering = ['preparation', 'fofo', 'analysis', 'fextr', 'alpha', 'refine']
        if self.fofo_coot_button is not None:
            script = self._fofo_coot_script()
            available = self.final_status == 'Completed' and script is not None
            self.fofo_coot_button.setVisible(available)
            self.fofo_coot_button.setEnabled(available)
        for key, (label, name) in self.steps.items():
            # A successful process exit confirms every FoFo-only stage.
            # The previous phase-based fallback left 'FoFo map calculation'
            # blue when the final phase remained 'fofo' in the log.
            if self.fofo and self.final_status == 'Completed':
                state = 'done'
            elif key == 'fextr':
                state = 'done' if self.total and completed == self.total else 'active' if self.phase == 'fextr' and self.active is not None else 'waiting'
            elif key == 'refine':
                state = 'done' if self.refine_total and ref_done == self.refine_total else 'active' if self.phase == 'refine' and self.refine_active is not None else 'waiting'
            elif key == 'alpha':
                state = 'done' if self.alpha_done else 'active' if self.phase == 'alpha' else 'waiting'
            elif key == 'analysis' and self.fofo:
                state = 'done' if self._fofo_analysis_done else 'active' if self.phase == 'analysis' and self.final_status is None else 'waiting'
            else:
                state = ('done' if self.phase in ordering and ordering.index(self.phase) > ordering.index(key) else
                         'active' if self.phase == key else 'waiting')
            symbol, color = {'done': ('✓', '#28835b'), 'active': ('●', '#2779ad'),
                             'waiting': ('○', 'gray')}[state]
            label.setText(f'{symbol}  {name}')
            label.setStyleSheet(f'color: {color};' + ('font-weight: 600;' if state == 'active' else ''))
        spinning = self._spinner_timer.isActive() and self.final_status is None
        self.spinner.setVisible(spinning)
        self.status_icon.setVisible(not spinning)
        self.status.setText('Done' if self.final_status == 'Completed' else
                            self.final_status or ('Running' if self.phase != 'waiting' else 'Waiting to start'))
        if self.final_status == 'Completed':
            self.status_icon.setText('✓')
            self.status_icon.setStyleSheet('color: #28835b; font-weight: bold;')
            self.status.setStyleSheet('color: #28835b; font-weight: 600;')
        elif self.final_status:
            self.status_icon.setText('■' if self.final_status == 'Stopped by user' else '✕')
            self.status_icon.setStyleSheet('color: gray;' if self.final_status == 'Stopped by user' else 'color: #bb3333;')
        elif self._spinner_timer.isActive():
            self.status_icon.setStyleSheet('color: #2779ad; font-weight: bold;')
        else:
            self.status_icon.setText('○')
