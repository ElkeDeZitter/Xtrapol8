from .execute import execute


def truncate(
    hklin: str,
    hklout: str,
    logout: str,
    labin_line: str,
    labout_line: str,
    ano: bool = False,
    high_res: float = None,
    low_res: float = None,
):
    args = ["HKLIN", hklin, "HKLOUT", hklout]
    stdin = ["truncate YES"]
    stdin.append(f"anomalous {'YES' if ano else 'NO'}")
    if high_res is not None and low_res is not None:
        stdin.append(f"resolution {high_res:.2f} {low_res:.2f}")
    stdin.append("plot OFF")
    stdin.append("header BRIEF BATCH")
    stdin.append(f"labin {labin_line}")
    stdin.append(f"labout {labout_line}")
    execute("truncate", args, stdin, logout)
