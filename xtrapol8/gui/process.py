"""Qt-only process controller for the unchanged Xtrapol8 scientific engine."""
from pathlib import Path
import shutil
from PyQt5.QtCore import QObject, QProcess, QProcessEnvironment, QTimer, pyqtSignal


class FextrProcess(QObject):
    log_received = pyqtSignal(str)
    started = pyqtSignal()
    finished = pyqtSignal(int, bool)  # exit code, was stopped by user
    error = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._read_output)
        self.process.started.connect(self.started.emit)
        self.process.finished.connect(self._finished)
        self.process.errorOccurred.connect(self._error)
        self._buffer = ''
        self._stopping = False
        self._active = False
        self._last_error = None

    @property
    def running(self):
        return self.process.state() != QProcess.NotRunning

    def start(self, phil_file, executable='xtrapol8', working_directory=None):
        if self.running:
            raise RuntimeError('A calculation is already running')
        phil_file = Path(phil_file).resolve()
        if not phil_file.is_file():
            raise FileNotFoundError(f'PHIL file not found: {phil_file}')
        resolved_executable = shutil.which(executable)
        if not resolved_executable:
            raise FileNotFoundError(
                f'Cannot find executable {executable!r} in PATH. '
                'Activate the environment containing xtrapol8 before starting the GUI.')
        self._stopping = False
        self._last_error = None
        self._buffer = ''
        workdir = Path(working_directory).resolve() if working_directory else phil_file.parent
        if not workdir.is_dir():
            raise NotADirectoryError(f'Working directory does not exist: {workdir}')
        self.process.setWorkingDirectory(str(workdir))
        env = QProcessEnvironment.systemEnvironment()
        env.insert("PYTHONUNBUFFERED", "1")
        env.insert("PYTHONIOENCODING", "utf-8")
        self.process.setProcessEnvironment(env)
        self.log_received.emit(f'[Command] {resolved_executable} {phil_file}')
        self.log_received.emit(f'[Working directory] {workdir}')
        self._active = True
        self.process.start(resolved_executable, [str(phil_file)])

    def stop(self):
        if not self.running or self._stopping:
            return
        self._stopping = True
        self.process.terminate()
        QTimer.singleShot(3000, self._kill_if_running)

    def _kill_if_running(self):
        if self._stopping and self.running:
            self.process.kill()

    def _read_output(self):
        self._buffer += bytes(self.process.readAllStandardOutput()).decode('utf-8', errors='replace')
        while '\n' in self._buffer:
            line, self._buffer = self._buffer.split('\n', 1)
            self.log_received.emit(line.rstrip('\r'))

    def _finished(self, exit_code, exit_status):
        self._read_output()
        if self._buffer:
            self.log_received.emit(self._buffer)
            self._buffer = ''
        self._active = False
        self.finished.emit(exit_code, self._stopping)

    def _error(self, error):
        # QProcess.Crashed is also emitted when terminate()/kill() stops a process.
        # That is expected for a user-requested stop, not a scientific crash.
        if self._stopping and error == QProcess.Crashed:
            return
        self._last_error = self.process.errorString()
        self.error.emit(self._last_error)
        if error == QProcess.FailedToStart:
            self._active = False
