from os import environ
from pathlib import Path

from iotbx.file_reader import any_file

from ...programs.scaleit import scaleit
from ..utils import check_and_remove_files


def test_rnase():
    path = Path(environ["CCP4"], "examples", "rnase", "rnase25.mtz")
    file = any_file(str(path), force_type="hkl", raise_sorry_if_errors=True)
    arrays = file.file_object.as_miller_arrays()

    fnat = next(a for a in arrays if "FNAT" in a.info().labels)
    fiod = next(a for a in arrays if "FIOD25" in a.info().labels)

    scaled = scaleit(f_obs_ref=fnat, f_obs_2=fiod, b_scaling="no")

    assert scaled.is_xray_amplitude_array()
    assert not scaled.is_empty()
    check_and_remove_files("scaleit.log", "forscaleit.mtz", "fromscaleit.mtz")
