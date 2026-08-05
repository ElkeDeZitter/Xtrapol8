from os import environ
from pathlib import Path

from ...programs.fft import fft
from ..utils import check_and_remove_files


def test_toxd():
    hklin = Path(environ["CCP4"], "examples", "toxd", "toxd_mapcoefs.mtz")
    fft(hklin, "fft.ccp4", "FWT", "PHWT")
    check_and_remove_files("fft.log", "fft.ccp4")
