"""Resolve input paths and reserve an unused Xtrapol8 output name.

Pure pathlib helpers; no Qt or scientific imports.
"""
from pathlib import Path
import re


def resolve_inputs(paths, bases):
    """Return (resolved, missing, ambiguous) for the supplied input paths.

    Only accept a unique existing candidate.  Never silently guess between
    distinct files with the same relative name.
    """
    resolved, missing, ambiguous = {}, [], {}
    bases = [Path(base).expanduser().resolve() for base in bases if base]
    for original in paths:
        path = Path(original).expanduser()
        if path.is_absolute():
            if path.is_file():
                resolved[original] = str(path.resolve())
            else:
                missing.append(original)
            continue
        candidates = list(dict.fromkeys(str((base / path).resolve()) for base in bases if (base / path).is_file()))
        if len(candidates) == 1:
            resolved[original] = candidates[0]
        elif len(candidates) > 1:
            ambiguous[original] = candidates
        else:
            missing.append(original)
    return resolved, missing, ambiguous


def unused_output_directory(requested):
    """Use requested if absent; otherwise increment trailing _NN or append _01.

    Does not create a directory: the scientific engine creates it.
    """
    requested = Path(requested).expanduser().resolve()
    if not requested.exists():
        return requested
    match = re.fullmatch(r'(.*)_(\d+)', requested.name)
    if match:
        stem, suffix = match.groups()
        width = max(2, len(suffix))
        start = int(suffix) + 1
    else:
        stem, width, start = requested.name, 2, 1
    for number in range(start, 100000):
        candidate = requested.with_name(f'{stem}_{number:0{width}d}')
        if not candidate.exists():
            return candidate
    raise RuntimeError('No unused output directory suffix available')
