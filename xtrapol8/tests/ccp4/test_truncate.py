from os import environ
from pathlib import Path

from ...programs.truncate import truncate
from ..utils import check_and_remove_files


def test_insulin():
    truncate(
        hklin=Path(environ["CCP4"], "examples", "data", "insulin.mtz"),
        hklout="truncate.mtz",
        logout="truncate.log",
        labin_line=" IMEAN=IMEAN SIGIMEAN=SIGIMEAN",
        labout_line=" F=F SIGF=SIGF",
    )
    check_and_remove_files("truncate.log", "truncate.mtz")
