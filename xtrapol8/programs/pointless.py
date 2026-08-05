from cctbx import miller

from .execute import execute, make_obs_mtz


def pointless(reflections: miller.array, reflections_ref: miller.array, prefix: str):
    """
    Run pointless with a reference data set.
    This is to make sure that the reflections have been processed with the same indexing possibility.
    This should normally be done for the triggered data set.
    Afterwards run truncate as usual to convert the intensities to structure factors
    """
    mtz_in = f"{prefix}_forpointless.mtz"
    mtz_in_ref = f"{prefix}_reference_forpointless.mtz"
    mtz_out = f"{prefix}_frompointless.mtz"
    log_out = f"pointless_{prefix}.log"

    inlabels = make_obs_mtz(reflections, mtz_in)
    inlabels_ref = make_obs_mtz(reflections_ref, mtz_in_ref)

    args = ["HKLIN", mtz_in, "HKLREF", mtz_in_ref, "HKLOUT", mtz_out]
    stdin = [
        "labin ".join(inlabels),
        "labref ".join(inlabels_ref),
    ]
    print("Running pointless to avoid indexing issues")
    return_code = execute("pointless", args, stdin, log_out)
    success = bool(return_code == 0)
    return mtz_out, success
