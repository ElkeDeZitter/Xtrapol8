import os
import shutil
import subprocess

from cctbx import miller

# CCP4 programs (Fortran) truncate file names longer than this.
CCP4_MAX_PATH_LENGTH = 255


def _shorten_path(arg):
    """
    Return a relative path (to the current working directory) for an absolute path
    that is too long for CCP4 programs, which would otherwise truncate it.
    Absolute paths that are short enough, and arguments that are not absolute paths
    to an existing location, are returned as is.
    """
    arg = str(arg)
    if len(arg) <= CCP4_MAX_PATH_LENGTH:
        return arg
    if not os.path.isabs(arg) or not os.path.isdir(os.path.dirname(arg)):
        return arg
    try:
        rel = os.path.relpath(arg)
    except ValueError:
        return arg
    shortest = rel if len(rel) < len(arg) else arg
    if len(shortest) > CCP4_MAX_PATH_LENGTH:
        print(
            "Warning: path longer than %d characters, it may be truncated by the "
            "CCP4 program. Use a shorter output directory.\n  %s"
            % (CCP4_MAX_PATH_LENGTH, shortest)
        )
    return shortest


def shorten_paths_in_keyword(line):
    """
    Shorten the absolute paths that appear as words in a keyword line (e.g. refmac
    'external @/path/to/restraints.txt'). A leading '@' and surrounding quotes are kept.
    """
    words = []
    for word in str(line).split(" "):
        prefix = "@" if word.startswith("@") else ""
        core = word[len(prefix) :]
        quote = ""
        if len(core) > 1 and core[0] == core[-1] and core[0] in "'\"":
            quote, core = core[0], core[1:-1]
        new_core = _shorten_path(core)
        words.append(word if new_core == core else prefix + quote + new_core + quote)
    return " ".join(words)


def _run_process(cmd_path, args, stdin, log_stream=None):
    with subprocess.Popen(
        args=[cmd_path] + [_shorten_path(arg) for arg in args],
        stdin=subprocess.PIPE if stdin else None,
        stdout=log_stream,
        encoding="utf-8",
    ) as process:
        if stdin:
            for line in stdin:
                process.stdin.write(f"{line}\n")
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
