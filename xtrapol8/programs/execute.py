import shutil
import subprocess

from cctbx import miller


def _run_process(cmd_path, args, stdin, log_stream=None):
    with subprocess.Popen(
        args=[cmd_path] + args,
        stdin=subprocess.PIPE if stdin else None,
        stdout=log_stream,
        encoding="utf-8",
    ) as process:
        if stdin:
            for line in stdin:
                process.stdin.write(line + "\n")
            process.stdin.close()
        process.wait()
        return process.returncode


def execute(cmd: str, args: list[str], stdin: list[str] = None, log_file: str = None):
    "Execute a command with arguments and input, writing output to a file."
    if not (cmd_path := shutil.which(cmd)):
        raise EnvironmentError(f"Executable '{cmd}' not found")
    if log_file:
        with open(log_file, "w", encoding="utf-8") as log_stream:
            return _run_process(cmd_path, args, stdin, log_stream)
    return _run_process(cmd_path, args, stdin)


def make_obs_mtz(reflections: miller.array, file_name: str):
    "Write the reflections to an mtz-file with default column labels."
    root = "I" if reflections.is_xray_intensity_array() else "F"
    reflections_ms = reflections.as_mtz_dataset(column_root_label=root)
    reflections_ms.mtz_object().write(file_name=file_name)
    labels = reflections_ms.column_labels()[3:]
    return labels
