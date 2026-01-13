from os import environ
from pathlib import Path

from iotbx.file_reader import any_file

from ...programs.pointless import pointless
from ..utils import check_and_remove_files


def test_rnase():
    path = Path(environ["CCP4"], "examples", "rnase", "rnase25.mtz")
    file = any_file(str(path), force_type="hkl", raise_sorry_if_errors=True)
    arrays = file.file_object.as_miller_arrays()

    fnat = next(a for a in arrays if "FNAT" in a.info().labels)
    fiod = next(a for a in arrays if "FIOD25" in a.info().labels)

    prefix = "triggered"
    mtz_out, success = pointless(fiod, fnat, prefix)

    assert success
    check_and_remove_files(
        f"{prefix}_forpointless.mtz",
        f"{prefix}_reference_forpointless.mtz",
        f"pointless_{prefix}.log",
        mtz_out,
    )
