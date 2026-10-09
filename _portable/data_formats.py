"""
Which window opens a data file, and in which format it is.

    detect(path) -> (window, kind)
        window: "lcms", "hrms", "deconvolute" or None (not a data file)
        kind:   "shimadzu_lcd", "agilent_chemstation", "agilent_masshunter",
                "bruker_d", "waters_raw", "thermo_raw", "mzml", "mzxml",
                "andi_cdf", "sciex_wiff" or None
    data_path(path) -> the file or folder to open (a file inside a .D/.d or
        .raw folder: that folder)

Only names, folder listings and a few header bytes are looked at (no data is
read and nothing heavy is imported), so this can run for every entry of a
folder listing.
"""
import os

KIND_NAMES = {
    "shimadzu_lcd": "LabSolutions .lcd",
    "agilent_chemstation": "Agilent .D",
    "agilent_masshunter": "Agilent MassHunter .d",
    "bruker_d": "Bruker .d",
    "waters_raw": "Waters .raw",
    "thermo_raw": "Thermo .raw",
    "mzml": "mzML",
    "mzxml": "mzXML",
    "andi_cdf": "ANDI .cdf",
    "sciex_wiff": "Sciex .wiff",
}

# files that open in LCMS Analysis (Sciex: only to say how to convert it)
LCMS_KINDS = ("shimadzu_lcd", "agilent_chemstation", "waters_raw", "thermo_raw", "mzml", "mzxml", "andi_cdf",
              "sciex_wiff")
# the file dialog of LCMS Analysis (folders are chosen as a file inside them)
LCMS_WILDCARD = ("LC-MS and HPLC data|*.lcd;*.uv;*.ch;*.ms;_FUNC*.DAT;_HEADER.TXT;*.raw;*.mzML;*.mzXML;*.cdf;"
                 "*.dx;*.wiff;*.wiff2|"
                 "Shimadzu LabSolutions (*.lcd)|*.lcd|"
                 "Agilent .D (a file inside), OpenLab .dx|*.uv;*.ch;*.ms;*.dx|"
                 "Waters .raw (a file inside)|_FUNC*.DAT;_HEADER.TXT|"
                 "Thermo (*.raw)|*.raw|"
                 "mzML, mzXML|*.mzML;*.mzXML|"
                 "ANDI netCDF (*.cdf)|*.cdf|"
                 "All files (*.*)|*.*")

_HR_THERMO = ("FTMS", "Orbitrap", "Exploris", "Q Exactive", "Astral", "LTQ FT")


def _names(folder):
    try:
        return os.listdir(folder)
    except OSError:
        return []


def _ext(path):
    return os.path.splitext(path.rstrip("\\/"))[1].lower()


def _head(path, n):
    try:
        with open(path, "rb") as fh:
            return fh.read(n)
    except OSError:
        return b""


def _folder_kind(d):
    """(window, kind) of a .D/.d or .raw folder, or (None, None)."""
    ext = _ext(d)
    names = _names(d)
    low = set(n.lower() for n in names)
    if ext == ".raw":
        if any(n.startswith("_func") and n.endswith(".dat") for n in low) or "_header.txt" in low:
            return "lcms", "waters_raw"
        return None, None
    if ext != ".d":
        return None, None
    if low & {"analysis.baf", "analysis.tdf", "analysis.tsf"}:
        return "hrms", "bruker_d"
    chem_ms = any(n.endswith(".ms") for n in low)
    chem_uv = any(n.endswith((".uv", ".ch")) for n in low)
    acq = os.path.join(d, next((n for n in names if n.lower() == "acqdata"), "AcqData"))
    acq_low = set(n.lower() for n in _names(acq)) if "acqdata" in low else set()
    if chem_ms:
        return "lcms", "agilent_chemstation"
    if acq_low & {"msprofile.bin", "msscan.bin", "mspeak.bin"}:
        return "hrms", "agilent_masshunter"
    if chem_uv or any(n.endswith((".sp", ".cg")) for n in acq_low):
        return "lcms", "agilent_chemstation"
    return None, None


def _thermo_high_res(path):
    """True when the start of a Thermo .raw file names an Orbitrap or FT
    analyzer (FTMS scan filters, instrument model)."""
    b = _head(path, 1 << 18)
    return any(k.encode("utf-16-le") in b for k in _HR_THERMO)


def data_path(path):
    """The data set a path stands for: a file inside a .D/.d or Waters .raw
    folder stands for that folder; anything else for itself."""
    p = path.rstrip("\\/")
    if os.path.isfile(p):
        up = os.path.dirname(p)
        if os.path.basename(up).lower() == "acqdata":
            up = os.path.dirname(up)
        if _ext(up) in (".d", ".raw") and os.path.isdir(up):
            return up
    return p


def detect(path):
    """(window, kind) of a data file or folder; (None, None) if unknown."""
    if not path:
        return None, None
    p = data_path(str(path))
    if os.path.isdir(p):
        return _folder_kind(p)
    if not os.path.isfile(p):
        return None, None
    ext = _ext(p)
    if ext == ".lcd":
        return "lcms", "shimadzu_lcd"
    if ext == ".raw":
        if _head(p, 18) != b"\x01\xa1" + "Finnigan".encode("utf-16-le"):
            return None, None
        return ("hrms" if _thermo_high_res(p) else "lcms"), "thermo_raw"
    if ext == ".mzml" or p.lower().endswith(".mzml.gz"):
        return "hrms", "mzml"
    if ext == ".mzxml":
        return "hrms", "mzxml"
    if ext == ".dx" and _head(p, 2) == b"PK":  # OpenLab CDS archive (a .dx text file is JCAMP-DX)
        return "lcms", "agilent_chemstation"
    if ext == ".cdf":
        return ("lcms", "andi_cdf") if _head(p, 3) == b"CDF" else (None, None)
    if ext in (".wiff", ".wiff2"):
        return "lcms", "sciex_wiff"
    return None, None


def kind_name(kind):
    return KIND_NAMES.get(kind, "")
