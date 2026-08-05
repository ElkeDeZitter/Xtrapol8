from .execute import execute


def dm(
    mtz_in: str,
    mtz_out: str,
    solc: float,
    combine: str,
    cycles: int,
    f_label: str,
    log_file: str,
    phi_label: str = "PHIC_ALL",
):
    args = ["hklin", mtz_in, "hklout", mtz_out]
    stdin = [
        f"SOLC {solc:.3f}",
        "MODE SOLV HIST MULTI SAYR",
        f"COMBINE {combine}",
        f"NCYC {cycles}",
        f"LABI FP={f_label} SIGFP=SIG{f_label} PHIO={phi_label} FOMO=FOM",
        "LABO FDM=FDM PHIDM=PHIDM",
    ]
    execute("dm", args, stdin, log_file)
