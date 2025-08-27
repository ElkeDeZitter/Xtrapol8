"Xtrapol8"

MAJOR = 1
MINOR = 2
PATCH = 9
TAG = "Dev"  # None or "Dev"

__version_date__ = "28 Jun 2025"
if TAG is None:
    __version_info__ = (MAJOR, MINOR, PATCH)
else:
    __version_info__ = (MAJOR, MINOR, PATCH, TAG)
__version__ = ".".join(map(str, __version_info__))

__authors__ = [
    "Elke De Zitter",
    "Nicolas Coquelle",
    "Paula Oeser",
    "Thomas Barends",
    "Jacques Philippe Colletier",
]
__license__ = "MIT"
__copyright__ = "Institut de Biologie Structurale - group DYNAMOP - team SNaX"
__citation__ = "doi.org/10.1038/s42003-022-03575-7"
