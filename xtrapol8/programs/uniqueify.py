from os import environ
from pathlib import Path
from uuid import uuid4

from gemmi import read_mtz_file

from ..Fextr_utils import redirect_stdout_and_stderr
from .execute import execute


def uniqueify(mtz_in: str):
    mtz = read_mtz_file(str(mtz_in))
    a, b, c, alpha, beta, gamma = mtz.cell.parameters
    xname = mtz.datasets[0].crystal_name
    dname = mtz.datasets[0].dataset_name

    temp_dir = Path(environ["CCP4_SCR"]).resolve()
    temp1 = temp_dir / f"uniq1_{uuid4()}.mtz"
    temp2 = temp_dir / f"uniq2_{uuid4()}.mtz"
    temp3 = temp_dir / f"uniq3_{uuid4()}.mtz"
    mtz_out = f"{Path(mtz_in).stem}_rfree.mtz"
    log_out = f"{Path(mtz_in).stem}_rfree.log"

    with redirect_stdout_and_stderr(log_out):
        execute(
            cmd="unique",
            args=["hklout", temp1],
            stdin=[
                f"CELL {a} {b} {c} {alpha} {beta} {gamma}",
                f"RESOLUTION {mtz.resolution_high():.3f} ",
                f"SYMMETRY '{mtz.spacegroup.hm}'",
                "LABOUT F=FUNI SIGF=SIGFUNI",
            ],
        )
        execute(
            cmd="freerflag",
            args=["HKLIN", temp1, "HKLOUT", temp2],
            stdin=["FREERFRAC 0.05", "END"],
        )
        execute(
            cmd="cad",
            args=["HKLIN1", mtz_in, "HKLIN2", temp2, "HKLOUT", temp3],
            stdin=[
                "LABI FILE 1  ALLIN",
                "LABI FILE 2  E1=FreeR_flag",
                f"XNAME FILE 2  E1={xname}",
                f"DNAME FILE 2  E1={dname}",
                "END",
            ],
        )
        execute(
            cmd="freerflag",
            args=["HKLIN", temp3, "HKLOUT", mtz_out],
            stdin=["COMPLETE FREE=FreeR_flag", "END"],
        )

    temp1.unlink(missing_ok=True)
    temp2.unlink(missing_ok=True)
    temp3.unlink(missing_ok=True)

    return mtz_out
