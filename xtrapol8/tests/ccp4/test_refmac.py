from os import environ
from pathlib import Path

from ...programs.refmac import refmac, refmac_for_dm
from ..utils import check_and_remove_files


def test_refmac_for_dm():
    refmac_for_dm(
        xyzin=Path(environ["CCP4"], "examples", "rnase", "rnase.pdb"),
        hklin=Path(environ["CCP4"], "examples", "rnase", "rnase25.mtz"),
        libins=[],
        f_label="FNAT",
    )
    check_and_remove_files(
        "rnase25_for_dm.mmcif",
        "rnase25_for_dm.mtz",
        "rnase25_for_dm.pdb",
        "rnase25_refmac_for_dm.log",
    )


def test_rnase():
    refmac(
        hklin=Path(environ["CCP4"], "examples", "rnase", "rnase25.mtz"),
        xyzin=Path(environ["CCP4"], "examples", "rnase", "rnase.pdb"),
        hklout="refmac.mtz",
        xyzout="refmac.pdb",
        log_file="refmac.log",
        f_label="FNAT",
    )
    check_and_remove_files("refmac.mmcif", "refmac.mtz", "refmac.pdb", "refmac.log")
