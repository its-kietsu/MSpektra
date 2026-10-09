"""
Libraries shipped as wheels in _portable\\wheels (reportlab, python-docx,
rainbow, ...): unpacked the first time they are needed into _portable\\pylibs
(or, when that folder is not writable, a folder in the user profile) and put
on sys.path. No installation, no admin rights.

    import pylibs_loader
    pylibs_loader.ensure_pylibs("rainbow")      # rainbow importable, or RuntimeError
"""
import os
import sys
import zipfile
import importlib
import threading

HERE = os.path.dirname(os.path.abspath(__file__))
WHEELS = os.path.join(HERE, "wheels")
_LOCK = threading.Lock()


def lib_dirs():
    """Folders tried for the unpacked libraries, first writable one used."""
    out = [os.path.join(HERE, "pylibs")]
    la = os.environ.get("LOCALAPPDATA")
    if la:
        out.append(os.path.join(la, "MS Analysis", "pylibs"))
    out.append(os.path.join(os.path.expanduser("~"), ".msanalysis", "pylibs"))
    return out


def _wheels(prefixes=None):
    """Wheel file names in _portable\\wheels; prefixes: only the wheels whose
    distribution name (the part before the first '-', case and '_'/'-'
    ignored) is one of them."""
    try:
        names = sorted(f for f in os.listdir(WHEELS) if f.endswith(".whl"))
    except OSError:
        return []
    if not prefixes:
        return names

    def norm(s):
        return s.lower().replace("-", "_")
    want = set(norm(p) for p in prefixes)
    return [f for f in names if norm(f.split("-", 1)[0]) in want]


def _import_all(modules):
    for m in modules:
        importlib.import_module(m)


def ensure_pylibs(*modules, wheels=None):
    """Make the modules importable (e.g. "rainbow", "reportlab", "docx").
    Already importable: returns "installed". Else unpacks the wheels (all of
    them, or the distributions named in wheels, e.g. ["rainbow_api"]) into the
    first writable folder of lib_dirs(), puts it first on sys.path and imports
    the modules; returns that folder. RuntimeError when it cannot."""
    try:
        _import_all(modules)
        if modules:
            return "installed"
    except ImportError:
        pass
    names = _wheels(wheels)
    if not names:
        raise RuntimeError("the libraries are missing (_portable\\wheels)")
    last = None
    with _LOCK:
        for d in lib_dirs():
            try:
                os.makedirs(d, exist_ok=True)
                for w in names:
                    mark = os.path.join(d, "." + w + ".ok")
                    if not os.path.exists(mark):
                        with zipfile.ZipFile(os.path.join(WHEELS, w)) as z:
                            z.extractall(d)
                        open(mark, "w").close()
                if d not in sys.path:
                    sys.path.insert(0, d)
                importlib.invalidate_caches()
                _import_all(modules)
                return d
            except (OSError, ImportError) as ex:
                last = ex
                continue
    raise RuntimeError("the libraries could not be unpacked (%s)" % (last or "no writable folder"))
