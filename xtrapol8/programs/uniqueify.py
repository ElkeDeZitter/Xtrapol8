from pathlib import Path

from .execute import execute


def uniqueify(mtz_in: str):
    stem = Path(mtz_in).stem
    mtz_out = f"{stem}_rfree.mtz"
    execute("uniqueify", [mtz_in, mtz_out])
    return mtz_out
