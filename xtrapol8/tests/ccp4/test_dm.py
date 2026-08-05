from os import environ
from pathlib import Path

from ...programs.dm import dm
from ..utils import check_and_remove_files


def test_rnase():
    dm(
        mtz_in=Path(environ["CCP4"], "examples", "data", "gere.mtz"),
        mtz_out="dm.mtz",
        solc=0.5,
        combine="OMIT",
        cycles=3,
        f_label="FPHASED",
        log_file="dm.log",
        phi_label="PHIB",
    )
    check_and_remove_files("dm.log", "dm.mtz")
