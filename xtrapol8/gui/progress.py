"""Workflow-aware, log-driven Xtrapol8 progress tracking (no engine changes)."""
import re
from pathlib import Path
from PyQt5.QtCore import QTimer, Qt
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QProgressBar, QHBoxLayout
from PyQt5.QtGui import QPainter, QPen, QColor

NUMBER = r'(?:[0-9]+(?:\.[0-9]+)?)'
START = re.compile(r'CALCULATING\s+(\S+)\s+TYPE OF ESFAS AND MAPS FOR OCCUPANCY\s+(' + NUMBER + r')', re.I)
RESULT = re.compile(r'--->\s*Results in\s+(.+?)\s*$', re.I)
OCC_DIR = re.compile(r'occupancy_(' + NUMBER + r')', re.I)
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
        self.refine_total = 0 if self.fofo else (1 if self.fast else self.total)
        self.phase = 'waiting'
        self.active = None
        self.fextr_started = set()
        self.fextr_completed = set()
        self.refine_started = set()
        self.refine_completed = set()
        self.refine_active = None
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
            for widget in (self.counter, self.bar, self.refine_counter, self.refine_bar, self.best_label, self.detail):
                layout.addWidget(widget)
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
        elif 'CALCULATE (Q/K-WEIGHTED) FEXTRAPOLATED MAPS' in upper:
            self.phase = 'fextr'

        match = START.search(line)
        if match:
            key = (match.group(1).lower(), _norm(match.group(2)))
            self.active = key
            self.fextr_started.add(key)
            self.phase = 'fextr'
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
                    self._last_detail = f'Fextr completed: {key[0]} — occupancy {occ:g}'
                    if self.active == key:
                        self.active = None

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
            if self.fast and 'FAST AND FURIOUS REFINEMENT' in upper:
                self.alpha_done = True
                self.phase = 'refine'
                key = ('best', self.best_occ)
                self.refine_active = key
                self.refine_started.add(key)
                self._last_detail = 'Refinement running at best occupancy'
            if REFINE_START.search(line) and (not self.fast or self.alpha_done):
                # For Calm'N'Curious, refinement is repeated per occupancy;
                # for Fast'N'Furious, it is a single best-occupancy refinement.
                occ_match = OCC_IN_TEXT.search(line)
                occ = _norm(occ_match.group(1)) if occ_match else None
                if self.fast:
                    key = ('best', self.best_occ)
                else:
                    if occ is None and self.active:
                        occ = self.active[1]
                    if occ is None and self.last_fextr_completed is not None:
                        occ = self.last_fextr_completed[1]
                    key = ('refine', occ) if occ is not None else None
                if key is not None:
                    self.refine_active = key
                    self.refine_started.add(key)
                    self._last_detail = ('Refinement running at best occupancy' if self.fast else
                                         f'Refinement running: occupancy {occ:g}')
                self.phase = 'refine'
            if ('FAST AND FURIOUS REFINEMENT DONE' in upper or REFINE_DONE.search(line)) and self.refine_active is not None:
                self.refine_completed.add(self.refine_active)
                self._last_detail = 'Refinement completed'
                self.refine_active = None
        self._render()

    def finish(self, exit_code, stopped):
        self._spinner_timer.stop()
        self.final_status = ('Stopped by user' if stopped else 'Completed' if exit_code == 0
                             else f'Failed (exit {exit_code})')
        if exit_code == 0 and not stopped and self.fast and self.best_occ is not None:
            self.alpha_done = True
        # Successful process exit confirms the final refinement, even if the
        # scientific program omitted an explicit refinement completion line.
        if exit_code == 0 and not stopped and self.refine_active is not None:
            self.refine_completed.add(self.refine_active)
            self.refine_active = None
        self._render()

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
        for key, (label, name) in self.steps.items():
            if key == 'fextr':
                state = 'done' if self.total and completed == self.total else 'active' if self.fextr_started else 'waiting'
            elif key == 'refine':
                state = 'done' if self.refine_total and ref_done == self.refine_total else 'active' if self.refine_started else 'waiting'
            elif key == 'alpha':
                state = 'done' if self.alpha_done else 'active' if self.phase == 'alpha' else 'waiting'
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
