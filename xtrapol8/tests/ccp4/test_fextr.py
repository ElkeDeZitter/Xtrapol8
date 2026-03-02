from os import chdir, environ
from pathlib import Path
from shutil import rmtree
from uuid import uuid4

from ...Fextr import run
from ..utils import download_mtz, download_pdb


def test_fextr():
    original_cwd = Path.cwd()
    try:
        lib_5sq = Path(environ["CLIBD_MON"], "5", "5SQ.cif")
        lib_nfa = Path(environ["CLIBD_MON"], "n", "NFA.cif")
        lib_rc7 = Path(environ["CLIBD_MON"], "r", "RC7.cif")
        outdir = Path(f"Fextr_test_{uuid4().hex}").resolve()
        with download_pdb("6gp0") as pdb_6gp0:
            with download_mtz("6gp0") as mtz_6gp0:
                with download_mtz("6gp1") as mtz_6gp1:
                    args = [
                        f"reference_pdb={pdb_6gp0}",
                        f"reference_mtz={mtz_6gp0}",
                        f"triggered_mtz={mtz_6gp1}",
                        f"additional_files={lib_5sq}",
                        f"additional_files={lib_nfa}",
                        f"additional_files={lib_rc7}",
                        f"outdir={outdir}",
                        "reciprocal_space=refmac5",
                        "real_space=coot",
                        "open_coot=False",
                    ]
                    run(args)
                    rmtree(outdir)
    finally:
        chdir(original_cwd)
