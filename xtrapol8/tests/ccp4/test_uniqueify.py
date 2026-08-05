from os import environ
from pathlib import Path

from ...programs.uniqueify import uniqueify
from ..utils import check_and_remove_files


def test_7ocn():
    mtz_in = Path(environ["CCP4"], "examples", "data", "7ocn.mtz")
    mtz_out = uniqueify(mtz_in)
    log_out = Path(mtz_out).with_suffix(".log")
    check_and_remove_files(mtz_out, log_out)
