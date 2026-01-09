from .execute import execute


def fft(hklin: str, mapout: str, f_label: str, phi_label: str):
    args = ["hklin", hklin, "mapout", mapout]
    stdin = [f"LABI F1={f_label} PHI={phi_label}"]
    execute("fft", args, stdin, "fft.log")
