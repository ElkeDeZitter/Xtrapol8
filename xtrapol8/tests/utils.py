import gzip
from contextlib import contextmanager
from os.path import basename
from pathlib import Path
from shutil import copyfileobj
from subprocess import call
from tempfile import NamedTemporaryFile
from urllib.parse import unquote, urlparse
from urllib.request import urlopen

from iotbx.file_reader import any_file

WWPDB_URL = "https://files.wwpdb.org/pub/pdb/data/structures/divided"


@contextmanager
def download(url: str):
    """
    Downloads a file from the given URL and saves it to a temporary file.
    Yields a string path to the temporary file.
    Use in a with statement to ensure the file is deleted afterwards.
    """
    url_name = unquote(basename(urlparse(url).path))
    with urlopen(url, timeout=30) as response:
        name = response.headers.get_filename() or url_name
        name = name.strip().replace(" ", "_")
        name = "".join(c for c in name if c.isalnum() or c in "-_.")
        with NamedTemporaryFile(suffix=f"_{name}", delete=False) as temp:
            while chunk := response.read(1_000_000):
                temp.write(chunk)
        path = Path(temp.name).resolve()
        try:
            yield path
        finally:
            path.unlink(missing_ok=True)


@contextmanager
def download_mtz(pdb_id: str):
    pdb_id = pdb_id.lower()
    url = f"{WWPDB_URL}/structure_factors/{pdb_id[1:3]}/r{pdb_id}sf.ent.gz"
    with download(url) as cif_path:
        mtz_path = cif_path.with_suffix(".mtz")
        call(["gemmi", "cif2mtz", str(cif_path), str(mtz_path)])
        try:
            yield mtz_path
        finally:
            mtz_path.unlink(missing_ok=True)


@contextmanager
def download_pdb(pdb_id: str):
    pdb_id = pdb_id.lower()
    url = f"{WWPDB_URL}/pdb/{pdb_id[1:3]}/pdb{pdb_id}.ent.gz"
    with download(url) as gz_path:
        pdb_path = gz_path.with_suffix(".pdb")
        with gzip.open(gz_path, "rb") as f_in:
            with pdb_path.open("wb") as f_out:
                copyfileobj(f_in, f_out)
                try:
                    yield pdb_path
                finally:
                    pdb_path.unlink(missing_ok=True)


def check_and_remove_files(*names):
    """
    Checks if files exist and are non-empty, then removes them.
    Ensures MTZ files are recognised by iotbx.
    """
    for name in names:
        path = Path(name)
        assert path.is_file(), f"Expected file {path} not found."
        assert path.stat().st_size > 0, f"File {path} is empty."
        if path.suffix == ".mtz":
            any_file(str(path), force_type="hkl", raise_sorry_if_errors=True)
        path.unlink(missing_ok=True)
