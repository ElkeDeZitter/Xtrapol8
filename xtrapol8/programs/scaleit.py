from cctbx import miller
from iotbx.file_reader import any_file

from ..log import LOG
from .execute import execute


def _make_mtz_for_scaleit(f_obs_ref: miller.array, f_obs_2: miller.array):
    mtz_out = "forscaleit.mtz"
    f_obs_ms = f_obs_ref.as_mtz_dataset(column_root_label="F_obs_ref")
    f_obs_ms.add_miller_array(f_obs_2, column_root_label="F_obs_2")
    f_obs_ms.mtz_object().write(file_name=mtz_out)
    return mtz_out


def _call_scaleit(mtz_in, b_scaling, low_res, high_res):
    mtz_out = "fromscaleit.mtz"
    args = ["HKLIN", mtz_in, "HKLOUT", mtz_out]
    stdin = [
        f"REFINE {b_scaling}",
        f"RESOLUTION {low_res:.2f} {high_res:.2f}",
        "LABIN FP = F_obs_ref SIGFP = SIGF_obs_ref  FPH1 = F_obs_2 SIGFPH1 = SIGF_obs_2",
    ]
    # Optional input for scaleit
    # NOWT
    # converge NCYC 4
    # converge ABS 0.001
    # converge TOLR -7
    LOG.info("Running scaleit, see scaleit.log")
    execute("scaleit", args, stdin, "scaleit.log")
    return mtz_out


def scaleit(
    f_obs_ref: miller.array,
    f_obs_2: miller.array,
    b_scaling: str,
    low_res: float = None,
    high_res: float = None,
):
    b_scaling = b_scaling.upper() if b_scaling.endswith("tropic") else "SCALE"
    dmax, dmin = f_obs_2.d_max_min()
    low_res = low_res or dmax
    high_res = high_res or dmin

    mtz_forscaleit = _make_mtz_for_scaleit(f_obs_ref, f_obs_2)
    mtz_fromscaleit = _call_scaleit(mtz_forscaleit, b_scaling, low_res, high_res)

    file = any_file(mtz_fromscaleit, force_type="hkl", raise_sorry_if_errors=True)
    arrays = file.file_object.as_miller_arrays()
    f_obs_2 = next(a for a in arrays if "F_obs_2" in a.info().labels)

    return f_obs_2
