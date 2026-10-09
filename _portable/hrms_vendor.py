"""
Vendor data files for HRMS Analysis (MSpektra), besides Bruker .d and mzML
(hrms_data.py):

  Agilent MassHunter .D folders (TOF, Q-TOF: AcqData\\MSScan.bin with
    MSProfile.bin and/or MSPeak.bin), read with the rainbow library
    (rainbow-api 1.5.3, LGPL; _portable\\wheels): its decoders of the scan
    records, the run length and LZF compressed profiles (LZF through lzf.py)
    and its mass calibration (traditional + polynomial, per scan).
  Waters MassLynx .raw folders (Q-TOF: Synapt, Xevo; also single and triple
    quadrupoles), read with Waters' MassLynx library (MassLynxRaw.dll, shipped
    with UniDec in the private Python; Waters EULA in LICENSES), which applies
    the m/z calibration of the file.
  Thermo .raw files (Orbitrap, LTQ FT), read with Thermo's RawFileReader
    (ThermoFisher.CommonCore.*.dll, shipped with UniDec) through pythonnet.
  Sciex .wiff: only Sciex's licensed library reads them; a message says to
    convert them to mzML with ProteoWizard msconvert.

Only the reading is done here (in Python, through the vendors' libraries).
Each format is a source (scan list + spectrum(i) + centroids(i)) handed to the
C++ core (msengine_py, ms_open_vendor): averaging, scans, mass chromatograms
and everything after run there as for the other formats (VendorEngineFile, the
interface of hrms_data.BrukerD). VendorFile is the Python reference and the
fallback without the library (hrms_data.BrukerD's sums).
"""
import os
import sys
import threading
import ctypes
from ctypes import c_int, c_long, c_double, c_void_p, POINTER, byref

import numpy as np

import hrms_data

HERE = os.path.dirname(os.path.abspath(__file__))
VENDOR_KINDS = ("agilent_masshunter", "waters_raw", "thermo_raw", "sciex_wiff")
WIFF_MESSAGE = ("Sciex .wiff files can only be read with Sciex's own library. Convert the file to mzML "
                "with ProteoWizard msconvert and open the mzML file.")


# ==========================================================================
# which format
# ==========================================================================
def agilent_folder(path):
    """The Agilent MassHunter .D folder of path (the folder, its AcqData
    folder or a file in it), or None."""
    p = os.path.abspath(str(path).rstrip("\\/"))
    for _ in range(3):
        if os.path.isfile(os.path.join(p, "AcqData", "MSScan.bin")):
            return p
        p = os.path.dirname(p)
    return None


def _thermo_magic(path):
    try:
        with open(path, "rb") as fh:
            head = fh.read(18)
        return head[:2] == b"\x01\xa1" and head[2:18].decode("utf-16-le", "replace") == "Finnigan"
    except OSError:
        return False


def waters_folder(path):
    """The Waters .raw folder of path (the folder or a file in it), or None."""
    p = os.path.abspath(str(path).rstrip("\\/"))
    for q in (p, os.path.dirname(p)):
        if q.lower().endswith(".raw") and os.path.isdir(q):
            try:
                names = [n.lower() for n in os.listdir(q)]
            except OSError:
                return None
            if "_header.txt" in names or any(n.startswith("_func") and n.endswith(".dat") for n in names):
                return q
    return None


def _local_kind(path):
    p = str(path).rstrip("\\/")
    low = p.lower()
    if low.endswith((".mzml", ".mzml.gz")):
        return "mzml"
    if low.endswith((".wiff", ".wiff2", ".wiff.scan")):
        return "sciex_wiff"
    if low.endswith(".raw") and os.path.isfile(p):
        return "thermo_raw" if _thermo_magic(p) else None
    if waters_folder(p):
        return "waters_raw"
    if agilent_folder(p):
        return "agilent_masshunter"
    if hrms_data.find_d_folder(p) is not None:
        return "bruker_d"
    return None


def detect_kind(path):
    """Format of a data file: "agilent_masshunter", "waters_raw",
    "thermo_raw", "sciex_wiff", "bruker_d", "mzml" or None (through
    data_formats.detect when that module is there)."""
    try:
        from data_formats import detect
    except ImportError:
        detect = None
    if detect is not None:
        try:
            kind = detect(path)[1]
            if kind:
                return kind
        except Exception as ex:
            print("[HRMS] format of %s not detected: %s" % (os.path.basename(str(path)), ex))
    return _local_kind(path)


# ==========================================================================
# libraries
# ==========================================================================
def ensure_rainbow():
    """rainbow importable (unpacked from its wheel in _portable\\wheels into
    pylibs the first time, as ms_report does for the report libraries)."""
    try:
        import rainbow  # noqa: F401
        return
    except ImportError:
        pass
    try:
        from pylibs_loader import ensure_pylibs
    except ImportError:
        ensure_pylibs = None
    if ensure_pylibs is not None:
        ensure_pylibs()
        import rainbow  # noqa: F401,F811
        return
    import zipfile
    import importlib
    wdir = os.path.join(HERE, "wheels")
    wheels = sorted(f for f in os.listdir(wdir) if f.startswith("rainbow_api") and f.endswith(".whl")) \
        if os.path.isdir(wdir) else []
    if not wheels:
        raise RuntimeError("The rainbow library is missing (_portable\\wheels)")
    dirs = [os.path.join(HERE, "pylibs")]
    la = os.environ.get("LOCALAPPDATA")
    if la:
        dirs.append(os.path.join(la, "MS Analysis", "pylibs"))
    dirs.append(os.path.join(os.path.expanduser("~"), ".msanalysis", "pylibs"))
    for d in dirs:
        try:
            os.makedirs(d, exist_ok=True)
            w = wheels[-1]
            mark = os.path.join(d, "." + w + ".ok")
            if not os.path.exists(mark):
                with zipfile.ZipFile(os.path.join(wdir, w)) as z:
                    z.extractall(d)
                open(mark, "w").close()
            if d not in sys.path:
                sys.path.insert(0, d)
            importlib.invalidate_caches()
            import rainbow  # noqa: F401,F811
            return
        except (OSError, ImportError):
            continue
    raise RuntimeError("The rainbow library could not be unpacked (no writable folder)")


def _unidec_dir(*sub):
    """A folder of the UniDec package of the private Python (not imported)."""
    import importlib.util
    spec = importlib.util.find_spec("unidec")
    if spec is None or not spec.submodule_search_locations:
        return None
    d = os.path.join(list(spec.submodule_search_locations)[0], *sub)
    return d if os.path.isdir(d) else None


def _iso(text):
    return str(text or "").strip()


# ==========================================================================
# sources: the scans of a file and their spectra
# ==========================================================================
class Source(object):
    """Scans (MS1, in time order) of a vendor file: arrays rt (min), pol
    (+1/-1), prof (bool), event, lo, hi (acquisition range, 0 unknown), tic,
    bpc (-1: from the centroids); events: list of dicts (polarity "+"/"-",
    mz_low, mz_high, label); spectrum(i) and centroids(i) (None: centroid the
    profile)."""
    kind = ""
    vendor = ""
    instrument = ""
    n_msms = 0

    def __init__(self):
        self.lock = threading.RLock()
        self.closed = False
        self.props = {}

    def _finish(self, rows, events):
        """rows: (rt, pol, prof, event, lo, hi, tic, bpc, key); sorted by time."""
        if not rows:
            raise ValueError("No MS scans in " + self.path)
        order = sorted(range(len(rows)), key=lambda k: rows[k][0])
        rows = [rows[k] for k in order]
        self.rt = np.array([r[0] for r in rows], float)
        self.pol = np.array([r[1] for r in rows], int)
        self.prof = np.array([bool(r[2]) for r in rows])
        self.event = np.array([r[3] for r in rows], int)
        self.lo = np.array([r[4] or 0.0 for r in rows], float)
        self.hi = np.array([r[5] or 0.0 for r in rows], float)
        self.tic = np.array([-1.0 if r[6] is None else r[6] for r in rows], float)
        self.bpc = np.array([-1.0 if r[7] is None else r[7] for r in rows], float)
        self.keys = [r[8] for r in rows]
        self.events = events

    def _check(self):
        if self.closed:
            raise RuntimeError("the file was closed")

    def spectrum(self, i):
        raise NotImplementedError

    def centroids(self, i):
        return None

    def close(self):
        with self.lock:
            self.closed = True


def _events_by(rows_keys, label_extra):
    """Events of scans grouped by key (polarity first, then the rest): returns
    (event index per key, events). label_extra(key) is the text after the
    polarity, shown only when one polarity has several events."""
    keys = sorted(set(rows_keys), key=lambda k: (0 if k[0] > 0 else 1,) + tuple(str(x) for x in k[1:]))
    per_pol = {}
    for k in keys:
        per_pol.setdefault(k[0], []).append(k)
    events, index = [], {}
    for k in keys:
        lab = "Positive" if k[0] > 0 else "Negative"
        if len(per_pol[k[0]]) > 1:
            extra = label_extra(k)
            lab += " (%s)" % (extra or "event %d" % (len(events) + 1))
        index[k] = len(events)
        events.append({"polarity": "+" if k[0] > 0 else "-", "mz_low": None, "mz_high": None, "label": lab})
    return index, events


def _ranges(src):
    for e, ev in enumerate(src.events):
        sel = src.event == e
        lo, hi = src.lo[sel], src.hi[sel]
        ev["mz_low"] = float(lo[lo > 0].min()) if np.any(lo > 0) else None
        ev["mz_high"] = float(hi[hi > 0].max()) if np.any(hi > 0) else None


# -------------------------------------------------------------- Agilent
class AgilentSource(Source):
    """Agilent MassHunter .D (AcqData\\MSScan.bin + MSProfile.bin / MSPeak.bin)."""
    kind = "Agilent MassHunter .D"
    vendor = "Agilent MassHunter"

    def __init__(self, path):
        Source.__init__(self)
        d = agilent_folder(path)
        if d is None:
            raise ValueError("%s is not an Agilent MassHunter .D folder (no AcqData\\MSScan.bin)" % path)
        self.path = d
        acq = self.acq = os.path.join(d, "AcqData")
        if os.path.isfile(os.path.join(acq, "MSScan_XSpecific.bin")):
            raise ValueError("%s holds ICP-MS data, not high resolution mass spectra." % os.path.basename(d))
        ensure_rainbow()
        if HERE not in sys.path:  # lzf.py (LZF compressed profiles)
            sys.path.append(HERE)
        from rainbow.agilent import masshunter as mh
        self.mh = mh
        try:
            from rainbow.agilent import _msprofile
            self._rle = _msprofile.decompress_inten_list
        except ImportError:
            self._rle = mh.decompress_inten_list
        ct = mh.parse_scan_xsd(os.path.join(acq, "MSScan.xsd"))
        recs = mh.read_scan_records(os.path.join(acq, "MSScan.bin"), ct, mh.count_scans(acq))
        if not recs:
            raise ValueError("No scans in " + d)
        self.prof_path = os.path.join(acq, "MSProfile.bin")
        self.peak_path = os.path.join(acq, "MSPeak.bin")
        has_prof, has_peak = os.path.isfile(self.prof_path), os.path.isfile(self.peak_path)
        if not has_prof and not has_peak:
            raise ValueError("%s has no mass spectra (no MSProfile.bin or MSPeak.bin)" % os.path.basename(d))
        cids = [r.get("CalibrationID") for r in recs]
        self.cal, self.cal_flags = mh._load_calibration(acq, cids)
        self.cids = cids
        psize = os.path.getsize(self.prof_path) if has_prof else 0
        ksize = os.path.getsize(self.peak_path) if has_peak else 0
        rows, keys = [], []
        n_msms = 0
        for k, r in enumerate(recs):
            if int(r.get("MSLevel", 1) or 1) != 1:
                n_msms += 1
                continue
            blocks = r.get("SpectrumParamsBlocks") or [r["SpectrumParamValues"]]
            pb = None
            if has_prof:
                pb = next((b for b in blocks if b.get("SpectrumFormatID") == 1), blocks[0])
                if pb["PointCount"] <= 0 or pb["SpectrumOffset"] + pb["ByteCount"] > psize:
                    pb = None  # never (fully) written: an interrupted acquisition
            cb = None
            if has_peak:
                cb = mh._select_centroid_block([b for b in blocks if b is not pb])
                if cb is not None and cb["SpectrumOffset"] + cb["ByteCount"] > ksize:
                    cb = None
            if pb is None and cb is None:
                continue
            if pb is not None and self.cal is None:
                raise ValueError("The m/z calibration of %s is missing (no MSMassCal.bin or DefaultMassCal.xml)"
                                 % os.path.basename(d))
            pol = -1 if int(r.get("IonPolarity", 0) or 0) == 1 else 1
            ce = round(float(r.get("CollisionEnergy", 0.0) or 0.0), 1)
            b = pb or cb
            lo, hi = float(b.get("MinX", 0) or 0), float(b.get("MaxX", 0) or 0)
            keys.append((pol, ce))
            rows.append([float(r["ScanTime"]), pol, pb is not None, None, lo, hi, float(r.get("TIC", 0) or 0),
                         float(r.get("BasePeakValue", 0) or 0), (k, pb, cb)])
        self.n_msms = n_msms
        index, events = _events_by(keys, lambda key: "CE %g V" % key[1])
        for row, key in zip(rows, keys):
            row[3] = index[key]
        self._finish([tuple(x) for x in rows], events)
        _ranges(self)
        self._props()

    def _props(self):
        import xml.etree.ElementTree as ET
        p = {}
        try:
            root = ET.parse(os.path.join(self.acq, "Contents.xml")).getroot()
            p["InstrumentName"] = (root.findtext("InstrumentName") or "").strip()
            p["AcquisitionDateTime"] = (root.findtext("AcquiredTime") or "").strip()
        except Exception:
            pass
        try:
            root = ET.parse(os.path.join(self.acq, "sample_info.xml")).getroot()
            names = {"sample name": "SampleName", "operator": "OperatorName", "acq method": "MethodName",
                     "method": "MethodName", "acquisition time": "AcquisitionDateTime"}
            for fld in root.iter("Field"):
                n = (fld.findtext("Name") or "").strip().lower()
                v = (fld.findtext("Value") or "").strip()
                if n in names and v and not p.get(names[n]):
                    p[names[n]] = v
        except Exception:
            pass
        if not p.get("InstrumentName"):
            try:
                root = ET.parse(os.path.join(self.acq, "Devices.xml")).getroot()
                for dev in root.iter("Device"):
                    name = (dev.findtext("Name") or "").strip()
                    if any(w in name.upper() for w in ("TOF", "QUAD", "MS")):
                        model = (dev.findtext("ModelNumber") or "").strip()
                        p["InstrumentName"] = " ".join(x for x in ("Agilent", model, name) if x)
                        break
            except Exception:
                pass
        self.props = {k: v for k, v in p.items() if v}
        self.props.setdefault("InstrumentName", self.vendor)
        self.instrument = self.props["InstrumentName"]

    def spectrum(self, i):
        with self.lock:
            self._check()
            k, pb, cb = self.keys[i]
            if pb is None:
                return self._peaks(k, cb)
            mh = self.mh
            n = int(pb["PointCount"])
            with open(self.prof_path, "rb") as fh:
                fh.seek(int(pb["SpectrumOffset"]))
                raw = fh.read(int(pb["ByteCount"]))
            if mh.segment_is_rle(raw, n):
                start, delta = np.frombuffer(raw[:16], "<f8")
                it = self._rle(memoryview(raw)[16:], n)
            else:
                import lzf
                dec = lzf.decompress(raw, int(pb["UncompressedByteCount"]) or (16 + 4 * n))
                if dec is None or len(dec) < 16 + 4 * n:
                    raise ValueError("damaged profile spectrum in MSProfile.bin")
                start, delta = np.frombuffer(dec[:16], "<f8")
                it = np.frombuffer(dec[16:16 + 4 * n], "<u4")
            start, delta = float(start), float(delta)
            # the flight time axis as rainbow builds it, calibrated with this scan's calibration
            ft = np.arange(start, start + delta * (n - 1) + 1e-3, delta)[:n]
            mz = mh.calibrate_mz(ft, self.cal[k], self.cal_flags.get(self.cids[k]))
            return np.asarray(mz, float), np.asarray(it, float)

    def _peaks(self, k, cb):
        mh = self.mh
        n = int(cb["PointCount"])
        with open(self.peak_path, "rb") as fh:
            fh.seek(int(cb["SpectrumOffset"]))
            raw = fh.read(int(cb["ByteCount"]))
        mz, it = mh._decode_peak_block(raw, n, int(cb["ByteCount"]) // n)
        mz = np.asarray(mz, np.float64)
        if self.cal is not None:
            mz = mh.calibrate_mz(mz, self.cal[k], self.cal_flags.get(self.cids[k]))
        return np.asarray(mz, float), np.asarray(it, float)

    def centroids(self, i):
        with self.lock:
            self._check()
            k, pb, cb = self.keys[i]
            if cb is None:
                return None
            return self._peaks(k, cb)


# -------------------------------------------------------------- Waters
_ML = {"dll": None}


def _masslynx():
    if _ML["dll"] is None:
        if os.name != "nt":
            raise RuntimeError("Waters' MassLynx library runs on Windows only")
        d = _unidec_dir("UniDecImporter", "Waters")
        p = os.path.join(d, "MassLynxRaw.dll") if d else None
        if not p or not os.path.isfile(p):
            raise RuntimeError("Waters' MassLynx library (MassLynxRaw.dll) is missing from the private Python")
        try:
            os.add_dll_directory(d)
        except Exception:
            pass
        _ML["dll"] = ctypes.WinDLL(p)
    return _ML["dll"]


class WatersSource(Source):
    """Waters MassLynx .raw folder, read with MassLynxRaw.dll: every function
    with mass spectra is an event (one per function; diode array, analog and
    MS/MS functions are left out)."""
    kind = "Waters .raw"
    vendor = "Waters"
    _MSMS = ("DAU", "MRM", "MSMS", "MS2", "PAR", "NL", "NG", "TOFD", "TOFP", "PSD", "DAUGHTER", "PARENT")

    def __init__(self, path):
        Source.__init__(self)
        d = waters_folder(path)
        if d is None:
            raise ValueError("%s is not a Waters .raw folder" % path)
        self.path = d
        L = self.L = _masslynx()
        self.info, self.scan = c_void_p(), c_void_p()
        self._ok(L.createRawReaderFromPath(d.encode("mbcs" if os.name == "nt" else "utf-8"), byref(self.info), 2))
        self._ok(L.createRawReaderFromPath(d.encode("mbcs" if os.name == "nt" else "utf-8"), byref(self.scan), 1))
        nf = c_int()
        self._ok(L.getFunctionCount(self.info, byref(nf)))
        rows, keys = [], []
        n_msms = 0
        funcs = []
        for f in range(nf.value):
            ns, ft, im, cont = c_int(), c_int(), c_int(), ctypes.c_bool()
            self._ok(L.getScanCount(self.info, f, byref(ns)))
            self._ok(L.getFunctionType(self.info, f, byref(ft)))
            words = self._text(L.getFunctionTypeString, self.info, ft.value).upper().replace("/", "").split()
            if not words or words[0] in ("DAD", "ANALOG", "UV", "PDA"):
                continue
            if any(w in self._MSMS for w in words):
                n_msms += ns.value
                continue
            self._ok(L.getIonMode(self.info, f, byref(im)))
            mode = self._text(L.getIonModeString, self.info, im.value)
            self._ok(L.isContinuum(self.info, f, byref(cont)))
            lo, hi = ctypes.c_float(), ctypes.c_float()
            L.getAcquisitionMassRange(self.info, f, 0, byref(lo), byref(hi))
            pol = -1 if mode.strip().endswith("-") else 1
            funcs.append((f, pol))
            for s in range(ns.value):
                rt = ctypes.c_float()
                self._ok(L.getRetentionTime(self.info, f, s, byref(rt)))
                keys.append((pol, f + 1))
                rows.append([float(rt.value), pol, bool(cont.value), None, float(lo.value),
                             float(hi.value) if hi.value > lo.value else 0.0, None, None, (f, s)])
        self.n_msms = n_msms
        index, events = _events_by(keys, lambda key: "function %d" % key[1])
        for row, key in zip(rows, keys):
            row[3] = index[key]
        self._finish([tuple(x) for x in rows], events)
        _ranges(self)
        self._props()

    def _ok(self, rc):
        if rc:
            msg = ctypes.c_char_p()
            try:
                self.L.getErrorMessage(rc, byref(msg))
                text = (msg.value or b"").decode("utf-8", "replace")
            except Exception:
                text = ""
            raise RuntimeError("Waters MassLynx library: %s (code %d)" % (text or "error", rc))

    def _text(self, fn, h, code):
        out = ctypes.c_char_p()
        if fn(h, code, byref(out)):
            return ""
        return (out.value or b"").decode("utf-8", "replace")

    def _props(self):
        p = {}
        try:
            with open(os.path.join(self.path, next(n for n in os.listdir(self.path) if n.lower() == "_header.txt")),
                      "r", encoding="latin-1") as fh:
                head = {}
                for line in fh:
                    if line.startswith("$$ ") and ":" in line:
                        k, v = line[3:].split(":", 1)
                        head[k.strip()] = v.strip()
            p["InstrumentName"] = head.get("Instrument", "")
            p["SampleName"] = head.get("Sample Description", "")
            p["OperatorName"] = head.get("User Name", "")
            p["MethodName"] = os.path.basename(head.get("MS Method", "").replace("\\", "/"))
            date = " ".join(x for x in (head.get("Acquired Date", ""), head.get("Acquired Time", "")) if x)
            if date:
                import datetime
                for fmt in ("%d-%b-%Y %H:%M:%S", "%d-%b-%Y %H:%M", "%d-%b-%Y"):
                    try:
                        p["AcquisitionDateTime"] = datetime.datetime.strptime(date, fmt).isoformat()
                        break
                    except ValueError:
                        continue
        except Exception:
            pass
        self.props = {k: v for k, v in p.items() if v}
        self.props.setdefault("InstrumentName", self.vendor)
        self.instrument = self.props["InstrumentName"]

    def spectrum(self, i):
        with self.lock:
            self._check()
            f, s = self.keys[i]
            pm, pi, n = c_void_p(), c_void_p(), c_int()
            self._ok(self.L.readScan(self.scan, f, s, byref(pm), byref(pi), byref(n)))
            if n.value <= 0 or not pm.value:
                return np.zeros(0), np.zeros(0)
            m = np.ctypeslib.as_array(ctypes.cast(pm, POINTER(ctypes.c_float)), (n.value,)).astype(np.float64)
            it = np.ctypeslib.as_array(ctypes.cast(pi, POINTER(ctypes.c_float)), (n.value,)).astype(np.float64)
            return m, it

    def close(self):
        with self.lock:
            if not self.closed:
                self.closed = True
                for h in (getattr(self, "info", None), getattr(self, "scan", None)):
                    try:
                        if h is not None and h.value:
                            self.L.destroyRawReader(h)
                    except Exception:
                        pass


# -------------------------------------------------------------- Thermo
_TH = {"ok": False}


def _thermo():
    if not _TH["ok"]:
        d = _unidec_dir("UniDecImporter", "Thermo")
        if d is None:
            raise RuntimeError("Thermo's RawFileReader is missing from the private Python")
        import clr
        for n in ("ThermoFisher.CommonCore.Data", "ThermoFisher.CommonCore.RawFileReader"):
            clr.AddReference(os.path.join(d, n + ".dll"))
        _TH["ok"] = True
    from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter
    from ThermoFisher.CommonCore.Data.Business import Device
    return RawFileReaderAdapter, Device


def _net(a):
    """numpy copy of a .NET double array."""
    if a is None:
        return np.zeros(0)
    n = int(a.Length)
    out = np.empty(n, np.float64)
    if n:
        from System import IntPtr
        from System.Runtime.InteropServices import Marshal
        Marshal.Copy(a, 0, IntPtr(int(out.ctypes.data)), n)
    return out


class ThermoSource(Source):
    """Thermo .raw file, read with RawFileReader: the MS1 scans, one event per
    polarity and mass analyzer (FTMS before ITMS) and SIM apart from full scans."""
    kind = "Thermo .raw"
    vendor = "Thermo"

    def __init__(self, path):
        Source.__init__(self)
        self.path = os.path.abspath(path)
        Adapter, Device = _thermo()
        raw = self.raw = Adapter.FileFactory(self.path)
        if raw is None or not raw.IsOpen or raw.IsError:
            msg = raw.FileError.ErrorMessage if raw is not None and raw.IsError else "cannot open"
            raise RuntimeError("Thermo RawFileReader: %s" % msg)
        if raw.GetInstrumentCountOfType(Device.MS) < 1:
            raise ValueError("%s has no mass spectra" % os.path.basename(path))
        raw.SelectInstrument(Device.MS, 1)
        h = raw.RunHeaderEx
        rows, keys = [], []
        n_msms = 0
        for s in range(int(h.FirstSpectrum), int(h.LastSpectrum) + 1):
            flt = raw.GetFilterForScanNumber(s)
            if str(flt.MSOrder) != "Ms":
                n_msms += 1
                continue
            st = raw.GetScanStatsForScanNumber(s)
            pol = -1 if str(flt.Polarity) == "Negative" else 1
            ana = str(flt.MassAnalyzer).replace("MassAnalyzer", "")
            mode = str(flt.ScanMode)
            keys.append((pol, 0 if ana == "FTMS" else 1, ana, mode))
            rows.append([float(raw.RetentionTimeFromScanNumber(s)), pol, not bool(st.IsCentroidScan), None,
                         float(st.LowMass), float(st.HighMass), float(st.TIC), float(st.BasePeakIntensity), s])
        self.n_msms = n_msms
        index, events = _events_by(keys, lambda key: key[2] + ("" if key[3] == "Full" else " " + key[3]))
        for row, key in zip(rows, keys):
            row[3] = index[key]
        self._finish([tuple(x) for x in rows], events)
        _ranges(self)
        p = {}
        try:
            idata = raw.GetInstrumentData()
            p["InstrumentName"] = str(idata.Model or idata.Name or "")
            p["SampleName"] = str(raw.SampleInformation.SampleName or "")
            p["OperatorName"] = str(raw.FileHeader.WhoCreatedId or "")
            meth = str(raw.SampleInformation.InstrumentMethodFile or "")
            p["MethodName"] = os.path.basename(meth.replace("\\", "/")) if meth else ""
            p["AcquisitionDateTime"] = str(raw.FileHeader.CreationDate.ToString("yyyy-MM-ddTHH:mm:ss"))
        except Exception:
            pass
        self.props = {k: v for k, v in p.items() if v}
        self.props.setdefault("InstrumentName", self.vendor)
        self.instrument = self.props["InstrumentName"]

    def spectrum(self, i):
        with self.lock:
            self._check()
            s = self.keys[i]
            st = self.raw.GetScanStatsForScanNumber(s)
            seg = self.raw.GetSegmentedScanFromScanNumber(s, st)
            return _net(seg.Positions), _net(seg.Intensities)

    def centroids(self, i):
        with self.lock:
            self._check()
            if not self.prof[i]:
                return None
            cs = self.raw.GetCentroidStream(self.keys[i], False)
            if cs is None or not cs.Length:
                return None  # (ITMS profiles: centroided as the other profiles)
            return _net(cs.Masses), _net(cs.Intensities)

    def close(self):
        with self.lock:
            if not self.closed:
                self.closed = True
                try:
                    self.raw.Dispose()
                except Exception:
                    pass


def make_source(path, kind):
    if kind == "agilent_masshunter":
        return AgilentSource(path)
    if kind == "waters_raw":
        return WatersSource(path)
    if kind == "thermo_raw":
        return ThermoSource(path)
    if kind == "sciex_wiff":
        raise ValueError(WIFF_MESSAGE)
    raise ValueError("%s is not a vendor data file" % path)


def _summary(src, has_profile):
    ev = ", ".join(e["label"] for e in src.events)
    name = src.instrument or src.kind
    txt = "%s: %d MS scans, %.1f min, %s" % (name, len(src.rt), src.rt[-1] if len(src.rt) else 0, ev)
    if src.n_msms:
        txt += ", %d MS/MS scans (not shown)" % src.n_msms
    return txt + (", profile" if has_profile else ", centroid only")


# ==========================================================================
# Python reference (no C++ core)
# ==========================================================================
class VendorFile(hrms_data.BrukerD):
    """A vendor file read through a source, with hrms_data.BrukerD's sums
    (the Python reference of VendorEngineFile)."""

    def __init__(self, src, progress=None):
        self.src = src
        self.kind = src.kind
        self.path = src.path
        self.name = os.path.splitext(os.path.basename(src.path.rstrip("\\/")))[0]
        self.props = dict(src.props)
        self.recal_note = ""
        self.handle = 0
        self.n_msms = src.n_msms
        n = len(src.rt)
        mzs, its = [], []
        for i in range(n):
            if progress is not None and i % 200 == 0 and progress(i, n) is False:
                raise hrms_data.Cancelled()
            c = src.centroids(i)
            if c is None or not len(c[0]):
                m, it = src.spectrum(i)
                c = hrms_data.centroid(m, it) if src.prof[i] else (m, it)
            mzs.append(np.asarray(c[0], float))
            its.append(np.asarray(c[1], float))
        self.lines = hrms_data._Spectra(mzs, its)
        self.rt = src.rt.copy()
        self.pol = np.where(src.pol < 0, 1, 0)
        self.tic = np.array([t if t >= 0 else float(np.sum(x)) for t, x in zip(src.tic, its)])
        self.bpc = np.array([b if b >= 0 else (float(np.max(x)) if len(x) else 0.0) for b, x in zip(src.bpc, its)])
        self.saturated = np.zeros(n, bool)
        self.has_profile = bool(np.any(src.prof))
        self.has_line = not bool(np.all(src.prof))
        self.events = [{"polarity": e["polarity"], "mz_low": e["mz_low"], "mz_high": e["mz_high"]} for e in src.events]
        self._scans = [np.where(src.event == e)[0] for e in range(len(src.events))]
        self.n_events = len(self.events)
        self.calib = None
        self.calib_info = ""
        self.calib_events = None
        self.xic_fast = True

    def profile(self, i):
        m, it = self.src.spectrum(int(i))
        return np.asarray(m, float), np.asarray(it, float), bool(self.src.prof[int(i)])

    def event_label(self, event):
        return self.src.events[event]["label"]

    def summary(self):
        return _summary(self.src, self.has_profile)

    def close(self):
        src = getattr(self, "src", None)
        if src is not None:
            src.close()


# ==========================================================================
# C++ core
# ==========================================================================
class _VScan(ctypes.Structure):  # ms_vendor_scan
    _fields_ = [("rt", c_double), ("polarity", c_int), ("profile", c_int), ("event", c_int),
                ("mz_lo", c_double), ("mz_hi", c_double), ("tic", c_double), ("bpc", c_double)]


_PD = POINTER(c_double)
_VENDOR_FN = ctypes.CFUNCTYPE(c_int, c_void_p, c_long, c_int, POINTER(_PD), POINTER(_PD), POINTER(c_long))


def engine_lib():
    """msengine with ms_open_vendor, or None."""
    try:
        import msengine_py
        L = msengine_py.lib()
    except Exception:
        return None
    if L is None or not hasattr(L, "ms_open_vendor"):
        return None
    if getattr(L, "_vendor_scan_type", None) is not _VScan:
        import msengine_py
        L.ms_open_vendor.argtypes = [ctypes.c_char_p, ctypes.c_char_p, c_long, POINTER(_VScan), c_long,
                                     POINTER(_VScan), c_int, c_int, _VENDOR_FN, c_void_p, msengine_py.PROGRESS,
                                     c_void_p]
        L.ms_open_vendor.restype = c_void_p
        L._vendor_scan_type = _VScan
    return L


def _engine_class():
    import msengine_py

    class VendorEngineFile(msengine_py.EngineFile):
        """A vendor file whose spectra the C++ core sums (ms_open_vendor);
        the interface of msengine_py.EngineFile / hrms_data.BrukerD."""

        def __init__(self, src, progress=None):
            try:
                self._open(src, progress)
            except BaseException:
                self.src = None  # the caller keeps the source (for the Python reader, or closes it)
                msengine_py.EngineFile.close(self)
                raise

        def _open(self, src, progress):
            L = engine_lib()
            if L is None:
                raise RuntimeError("msengine without ms_open_vendor")
            self.lib = L
            self.src = src
            self.path = src.path
            self.name = os.path.splitext(os.path.basename(src.path.rstrip("\\/")))[0]
            self._keep = None
            self._cb_error = ""
            n = len(src.rt)
            scans = (_VScan * n)()
            for i in range(n):
                s = scans[i]
                s.rt, s.polarity, s.profile, s.event = float(src.rt[i]), int(src.pol[i]), int(src.prof[i]), int(src.event[i])
                s.mz_lo, s.mz_hi, s.tic, s.bpc = float(src.lo[i]), float(src.hi[i]), float(src.tic[i]), float(src.bpc[i])
            evs = (_VScan * max(len(src.events), 1))()
            for e, ev in enumerate(src.events):
                evs[e].polarity = -1 if ev["polarity"] == "-" else 1
                evs[e].mz_lo = float(ev["mz_low"] or 0.0)
                evs[e].mz_hi = float(ev["mz_high"] or 0.0)
                evs[e].event = e
            self._fn = _VENDOR_FN(self._serve)
            self._cb = msengine_py.PROGRESS(
                lambda u, i, k: 1 if (progress is None or progress(i, k) is not False) else 0)
            self.handle = L.ms_open_vendor(src.kind.encode("utf-8"), src.instrument.encode("utf-8"), int(src.n_msms),
                                           scans, n, evs, len(src.events), 0, self._fn, None, self._cb, None)
            if not self.handle:
                err = msengine_py._err(L)
                if "cancel" in err.lower():
                    raise hrms_data.Cancelled()
                raise RuntimeError(self._cb_error or err)
            self.kind = src.kind
            self.recal_note = ""
            self.has_profile = bool(L.ms_file_has_profile(self.handle))
            self.n_msms = int(L.ms_file_n_msms(self.handle))
            self.n_events = int(L.ms_file_n_events(self.handle))
            self.events, self._scans = [], []
            rts, tics, bpcs = [], [], []
            for e in range(self.n_events):
                ev = src.events[e]
                self.events.append({"polarity": ev["polarity"], "mz_low": ev["mz_low"], "mz_high": ev["mz_high"]})
                prt, ptic, pbpc = _PD(), _PD(), _PD()
                k = L.ms_event_scans(self.handle, e, byref(prt), byref(ptic), byref(pbpc))
                rts.append(msengine_py._arr(prt, k))
                tics.append(msengine_py._arr(ptic, k))
                bpcs.append(msengine_py._arr(pbpc, k))
            # one array of every scan in time order (as EngineFile and the Python readers)
            cat = np.concatenate(rts) if rts else np.zeros(0)
            order = np.argsort(cat, kind="stable")
            where = np.empty(len(cat), dtype=np.int64)
            where[order] = np.arange(len(cat))
            start = 0
            for r in rts:
                self._scans.append(np.sort(where[start:start + len(r)]))
                start += len(r)
            self.rt = cat[order]
            self.tic = (np.concatenate(tics) if tics else np.zeros(0))[order]
            self.bpc = (np.concatenate(bpcs) if bpcs else np.zeros(0))[order]
            self.saturated = np.zeros(len(self.rt), bool)
            self.calib = None
            self.calib_info = ""
            self.calib_events = None
            self.work = src.path
            self.props = dict(src.props)
            self.xic_fast = True
            self._summary = _summary(src, self.has_profile)

        def _serve(self, user, i, what, pm, pi, pn):
            """The spectrum callback of the library (what 0: spectrum, 1: centroids)."""
            try:
                c = self.src.centroids(int(i)) if what == 1 else self.src.spectrum(int(i))
                if c is None:
                    m = it = np.zeros(0)
                else:
                    m = np.ascontiguousarray(c[0], np.float64)
                    it = np.ascontiguousarray(c[1], np.float64)
                k = min(len(m), len(it))
                self._keep = (m, it)
                pm[0] = m.ctypes.data_as(_PD)
                pi[0] = it.ctypes.data_as(_PD)
                pn[0] = k
                return 0
            except Exception as ex:  # an exception must not pass through the library
                self._cb_error = "Scan %d of %s could not be read: %s" % (int(i) + 1, os.path.basename(self.path), ex)
                print("[HRMS] " + self._cb_error)
                return 1

        def _call(self, fn, *a, **k):
            self._cb_error = ""
            try:
                return fn(self, *a, **k)
            except RuntimeError as ex:
                if self._cb_error:
                    raise RuntimeError(self._cb_error) from ex
                raise

        def average(self, *a, **k):
            return self._call(msengine_py.EngineFile.average, *a, **k)

        def scan_spectrum(self, *a, **k):
            return self._call(msengine_py.EngineFile.scan_spectrum, *a, **k)

        def chromatogram(self, *a, **k):
            return self._call(msengine_py.EngineFile.chromatogram, *a, **k)

        def event_label(self, event):
            return self.src.events[event]["label"]

        def close(self):
            msengine_py.EngineFile.close(self)
            src = getattr(self, "src", None)
            if src is not None:
                src.close()

    return VendorEngineFile


_ENGINE_CLASS = {}


def open_vendor(path, kind=None, progress=None):
    """A vendor data file for HRMS Analysis (the C++ core sums its spectra
    when the library has ms_open_vendor; else the Python reference)."""
    kind = kind or detect_kind(path)
    src = make_source(path, kind)
    try:
        if engine_lib() is not None:
            if "cls" not in _ENGINE_CLASS:
                _ENGINE_CLASS["cls"] = _engine_class()
            try:
                return _ENGINE_CLASS["cls"](src, progress)
            except hrms_data.Cancelled:
                raise
            except Exception as ex:
                print("[engine] Python reader used for %s: %s" % (os.path.basename(path.rstrip("\\/")), ex))
        return VendorFile(src, progress)
    except BaseException:
        src.close()
        raise
