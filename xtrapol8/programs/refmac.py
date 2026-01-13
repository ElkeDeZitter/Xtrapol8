from pathlib import Path
from .execute import execute


def refmac_for_dm(xyzin: str, hklin: str, libins: list[str], f_label: str):
    stem = Path(hklin).stem
    hklout = f"{stem}_for_dm.mtz"
    xyzout = f"{stem}_for_dm.pdb"
    log_file = f"{stem}_refmac_for_dm.log"

    args = ["XYZIN", xyzin, "HKLIN", hklin, "HKLOUT", hklout, "XYZOUT", xyzout]
    for libin in libins:
        args += ["LIB_IN", libin]

    stdin = [
        f"LABIN FP={f_label} SIGFP=SIG{f_label} FREE=FreeR_flag",
        "REFI TYPE REST RESI MLKF BREF ISOT METH CGMAT",
        "nfree include 1",
        "make check NONE",
        "ncyc 0",
        "MAPC SHAR",
        "NOHARVEST",
        "END",
    ]

    print(
        "Running refmac with zero cycles "
        "to prepare files suitable for running dm afterwards. "
        "Please wait..."
    )
    execute("refmac5", args, stdin, log_file)

    return hklout, xyzout


def refmac(
    hklin: str,
    xyzin: str,
    hklout: str,
    xyzout: str,
    log_file: str,
    libins: list[str] = None,
    f_label: str = "QFEXTR",
    cycles: int = 5,
    weight: str = "AUTO",
    refi_type: str = "RESTrained",
    refi_bref: str = "ISOTropic",
    twinning: bool = False,
    use_tls: bool = False,
    tls_cycles: int = 20,
    bfac_set: str = 30,
    use_jelly_body: bool = False,
    jelly_body_sigma: float = 0.03,
    additional_jelly_body_restraints: list[str] = None,
    external_restraints: list[str] = None,
    map_sharpening: bool = False,
    additional_keywords: list[str] = None,
):
    args = ["HKLIN", hklin, "HKLOUT", hklout, "XYZIN", xyzin, "XYZOUT", xyzout]
    for libin in libins or []:
        args += ["LIB_IN", libin]

    stdin = [
        "MAKE HYDR No",
        "MAKE CHEC NONE",
        "MAKE SS Yes",
        "MAKE SYMM Yes",
        "MAKE SUGA Yes",
        "MAKE CONN No",
        f"LABIN FP={f_label} SIGFP=SIG{f_label} FREE=FreeR_flag",
        "LABO FC=FC PHIC=PHIC FWT=2FOFCWT PHWT=PH2FOFCWT DELFWT=FOFCWT PHDELWT=PHFOFCWT",
        f"WEIGHT {weight}",
        "SCALe TYPE BULK",
        f"NCYC {cycles}",
        "MONI FEW",
        "NOHARVEST",
        f"REFI TYPE {refi_type}",
        "REFI RESI MLKF",
        f"REFI BREF {refi_bref}",
    ]
    if use_tls:
        stdin.append(f"REFI TLSC {tls_cycles}")
        stdin.append(f"BFAC SET {bfac_set}")
    if twinning:
        stdin.append("TWIN")
    if use_jelly_body:
        stdin.append(f"RIDG DIST {jelly_body_sigma}")
        for restraint in additional_jelly_body_restraints or []:
            stdin.append(f"RIDG {restraint}")
    for restraint in external_restraints or []:
        stdin.append(f"external {restraint}")
    if map_sharpening:
        stdin.append("MAPC SHAR")
    if additional_keywords:
        stdin.append(" ".join(additional_keywords))
    stdin.append("END")

    returncode = execute("refmac5", args, stdin, log_file)
    return returncode
