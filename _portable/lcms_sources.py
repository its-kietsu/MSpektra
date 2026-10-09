"""
Data sources of LCMS Analysis: the LC-MS and HPLC files of other vendors,
read into the objects the MS, PDA and Compare views use.

    ms   scans by scan event with the interface of lcms_data.LCDFile /
         msengine_py.EngineFile (rt, tic, bpc, events, event_scans,
         event_label, chromatogram, average, scan_spectrum, summary,
         saturated, props)
    pda  lcms_pda.PDAData (times in min, wavelengths in nm, A[time, wl] in mAU)

Formats (kinds of data_formats.detect):
    agilent_chemstation  Agilent ChemStation / OpenLab .D folders: DAD spectra
                         (*.uv), single wavelength signals (*.ch), single quad
                         MS (*.ms), through rainbow (LGPL, _portable\\wheels)
    waters_raw           Waters .raw folders: PDA and MS functions, through rainbow
    thermo_raw           Thermo .raw files: MS and PDA, through Thermo's
                         RawFileReader (.NET, shipped with the UniDec package)
    mzml                 msengine (C++) or hrms_data, as in HRMS Analysis
    mzxml                pyteomics
    andi_cdf             ANDI / AIA netCDF exports (MS or chromatograms)
Shimadzu .lcd files are read by unilcms as before (msengine, lcms_data).

Only the reading of the files is done here in Python. The scans go to the C++
core (msengine ms_open_arrays: sums on m/z bins and mass chromatograms as for
Shimadzu files); ArrayLCDFile (lcms_data.LCDFile on the same arrays) is the
Python reference and the fallback.
"""
import os
import re
import datetime
import threading

import numpy as np

import data_formats
import lcms_data
import lcms_pda

WIFF_TEXT = "Sciex .wiff files cannot be read. Convert them to mzML (ProteoWizard MSConvert)."
_RB_LOCK = threading.Lock()
_NET_LOCK = threading.Lock()
_NET = {}


# ==========================================================================
# MS data: the scans of a file as arrays
# ==========================================================================
class Scans(object):
    """MS1 scans collected by a reader, any order; build() makes the data object."""

    def __init__(self):
        self.blocks = []  # (event, rt, counts, mz, it, tic or None, bpc or None)
        self.events = []  # dicts: polarity ("+", "-" or None), mz_low, mz_high, label
        self.n_msms = 0

    @property
    def n_scans(self):
        return sum(len(b[1]) for b in self.blocks)

    def add_event(self, polarity=None, label="", mz_low=None, mz_high=None):
        self.events.append({"polarity": polarity, "mz_low": mz_low, "mz_high": mz_high, "label": label})
        return len(self.events) - 1

    def add(self, event, rt_min, mz, it, tic=None, bpc=None):
        """One scan (times in minutes)."""
        mz = np.asarray(mz, dtype=np.float64).ravel()
        it = np.asarray(it, dtype=np.float64).ravel()
        n = min(len(mz), len(it))
        self.add_block(event, [rt_min], [n], mz[:n], it[:n], None if tic is None else [tic],
                       None if bpc is None else [bpc])

    def add_block(self, event, rt_min, counts, mz, it, tic=None, bpc=None):
        """Scans of one event: their m/z intensity pairs one after the other
        (counts per scan). Zero and not finite points are left out."""
        counts = np.asarray(counts, dtype=np.int64)
        mz = np.asarray(mz, dtype=np.float64).ravel()
        it = np.asarray(it, dtype=np.float64).ravel()
        keep = (it != 0) & np.isfinite(mz) & np.isfinite(it)
        if not keep.all():
            row = np.repeat(np.arange(len(counts)), counts)
            counts = np.bincount(row[keep], minlength=len(counts)).astype(np.int64)
            mz, it = mz[keep], it[keep]
        self.blocks.append((int(event), np.asarray(rt_min, dtype=np.float64), counts, mz.astype(np.float32),
                            it.astype(np.float32), None if tic is None else np.asarray(tic, dtype=np.float64),
                            None if bpc is None else np.asarray(bpc, dtype=np.float64)))

    def arrays(self):
        """Scans in time order: rt, event, offsets, mz, it (float32), tic, bpc
        (the sum and the largest of the intensities when the file has none)."""
        rt = np.concatenate([b[1] for b in self.blocks])
        event = np.concatenate([np.full(len(b[1]), b[0], dtype=np.int64) for b in self.blocks])
        counts = np.concatenate([b[2] for b in self.blocks])
        mz0 = np.concatenate([b[3] for b in self.blocks])
        it0 = np.concatenate([b[4] for b in self.blocks])
        tic_f = np.concatenate([b[5] if b[5] is not None else np.full(len(b[1]), np.nan) for b in self.blocks])
        bpc_f = np.concatenate([b[6] if b[6] is not None else np.full(len(b[1]), np.nan) for b in self.blocks])
        off0 = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
        order = np.argsort(rt, kind="stable")
        rt, event, counts = rt[order], event[order], counts[order]
        offsets = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
        src = np.repeat(off0[order] - offsets[:-1], counts) + np.arange(int(offsets[-1]), dtype=np.int64)
        mz, it = mz0[src], it0[src]
        cs = np.concatenate([[0.0], np.cumsum(it.astype(np.float64))])
        tic = cs[offsets[1:]] - cs[offsets[:-1]]
        bpc = np.zeros(len(rt))
        ne = counts > 0
        if ne.any():
            bpc[ne] = np.maximum.reduceat(it, offsets[:-1][ne]).astype(np.float64)
        tf, bf = tic_f[order], bpc_f[order]
        tic = np.where(np.isfinite(tf), tf, tic)
        bpc = np.where(np.isfinite(bf), bf, bpc)
        return rt, event, offsets, mz, it, tic, bpc


class ArrayLCDFile(lcms_data.LCDFile):
    """lcms_data.LCDFile on scans given as arrays: the Python reference of
    msengine ms_open_arrays (the same float32 spectra, bins and windows)."""

    def __init__(self, path, rt, event, offsets, mz, it, tic, bpc, events):
        self.path = path
        self.name = os.path.splitext(os.path.basename(path.rstrip("\\/")))[0]
        self.reader = None
        self.rt = np.asarray(rt, dtype=np.float64)
        self.offsets = np.asarray(offsets, dtype=np.int64)
        self.mz = np.asarray(mz, dtype=np.float32)
        self.inten = np.asarray(it, dtype=np.float32)
        self.tic = np.asarray(tic, dtype=np.float64)
        self.bpc = np.asarray(bpc, dtype=np.float64)
        self.bpm = np.zeros(len(self.rt))
        self.n_events = len(events)
        self.scan_event = np.asarray(event, dtype=np.int64)
        self.events = [{"polarity": e["polarity"], "mz_low": e["mz_low"], "mz_high": e["mz_high"]} for e in events]
        self.saturated = np.zeros(len(self.rt), bool)
        self.mz_scale = 1.0


def _event_labels(events):
    """The polarity, or polarity and name when that is not unique."""
    pol = [{"+": "Positive", "-": "Negative"}.get(e["polarity"], "") for e in events]
    out = []
    for e, p in zip(events, pol):
        name = e.get("label") or ""
        if p and pol.count(p) == 1:
            out.append(p)
        else:
            out.append(" ".join(x for x in (name, p) if x) or "Scans")
    return out


def build_ms(path, kind, scans, props=None, engine=True):
    """The MS data object of the scans: msengine (C++) when it has
    ms_open_arrays, else ArrayLCDFile. None when there are no scans."""
    if not scans.n_scans:
        return None
    rt, event, offsets, mz, it, tic, bpc = scans.arrays()
    used = sorted(set(event.tolist()))
    if len(used) < len(scans.events):  # events without scans are left out
        remap = {e: k for k, e in enumerate(used)}
        event = np.array([remap[e] for e in event.tolist()], dtype=np.int64)
        scans.events = [scans.events[e] for e in used]
    events = scans.events
    labels = _event_labels(events)
    data = None
    if engine:
        try:
            import msengine_py
            pol = [-1 if e["polarity"] == "-" else (1 if e["polarity"] == "+" else 0) for e in events]
            data = msengine_py.open_arrays(path, data_formats.kind_name(kind) or kind, (props or {}).get(
                "InstrumentName", ""), rt, event, offsets, mz, it, tic, bpc, pol,
                [e["mz_low"] or 0 for e in events], [e["mz_high"] or 0 for e in events], scans.n_msms)
        except Exception as ex:
            print("[engine] Python reference used for %s: %s" % (os.path.basename(path.rstrip("\\/")), ex))
            data = None
        if data is not None:
            data.events = [{"polarity": e["polarity"], "mz_low": e["mz_low"], "mz_high": e["mz_high"]}
                           for e in events]
    if data is None:
        data = ArrayLCDFile(path, rt, event, offsets, mz, it, tic, bpc, events)
    data.labels = labels
    data.kind = data_formats.kind_name(kind) or kind
    data.props = dict(props or {})
    data.n_msms = scans.n_msms
    txt = "%d scans, %.1f min, %d scan event(s): %s" % (len(rt), rt[-1] if len(rt) else 0, len(events),
                                                        ", ".join(labels))
    if scans.n_msms:
        txt += ", %d MS/MS scans (not shown)" % scans.n_msms
    data._summary = txt
    data.summary = lambda: txt
    data.event_label = lambda e: labels[e]
    return data


# ==========================================================================
# PDA data
# ==========================================================================
def build_pda(times, wavelengths, A, units="mAU"):
    """lcms_pda.PDAData of an absorbance matrix (times x wavelengths); the
    values are converted to mAU (AU, mAU or uAU given)."""
    t = np.asarray(times, dtype=np.float64)
    w = np.asarray(wavelengths, dtype=np.float64)
    a = np.asarray(A, dtype=np.float64).reshape(len(t), len(w))
    u = (units or "mAU").strip().replace("µ", "u").replace("μ", "u").lower()
    scale = {"au": 1000.0, "mau": 1.0, "uau": 0.001}.get(u, 1.0)
    if len(t) > 1 and np.any(np.diff(t) < 0):
        o = np.argsort(t, kind="stable")
        t, a = t[o], a[o]
    o = np.argsort(w, kind="stable")
    obj = lcms_pda.PDAData.__new__(lcms_pda.PDAData)
    obj.wavelengths = w[o]
    obj.A = np.ascontiguousarray(a[:, o] * scale, dtype=np.float32)
    obj.times = t
    obj.interval_s = float(np.median(np.diff(t)) * 60.0) if len(t) > 1 else 0.16
    obj.units = "mAU"
    return obj


def pda_from_channels(channels):
    """PDAData from single wavelength chromatograms [(wavelength, times, y,
    units)]: on the time axis of the longest one (the others interpolated),
    one column per wavelength (the first of a wavelength recorded twice)."""
    chans, seen = [], set()
    for wl, t, y, u in channels:
        if wl in seen or len(t) < 2:
            continue
        seen.add(wl)
        chans.append((float(wl), np.asarray(t, float), np.asarray(y, float), u))
    if not chans:
        return None
    ref = max(chans, key=lambda c: len(c[1]))[1]
    cols, wls = [], []
    for wl, t, y, u in chans:
        scale = {"au": 1000.0, "uau": 0.001}.get((u or "mAU").replace("µ", "u").lower(), 1.0)
        yy = y if (len(t) == len(ref) and np.allclose(t, ref)) else np.interp(ref, t, y, left=0.0, right=0.0)
        cols.append(yy * scale)
        wls.append(wl)
    return build_pda(ref, wls, np.column_stack(cols), "mAU")


# ==========================================================================
# rainbow: Agilent .D and Waters .raw
# ==========================================================================
_SENTINEL_BIN = 0.0123456789  # bin width asked from rainbow: its MS pairs are taken unbinned


def _rainbow():
    import pylibs_loader
    pylibs_loader.ensure_pylibs("rainbow", wheels=["rainbow_api"])
    import rainbow
    return rainbow


def _rainbow_read(path):
    """rainbow.read of a .D or .raw folder with the m/z intensity pairs of
    every MS channel kept as recorded (rainbow sums them on an m/z grid,
    as wide as the scan range at the finest bin width; here they are kept per
    scan instead). Returns (DataDirectory, {placeholder: (mz, it, counts)})."""
    rb = _rainbow()
    import sys
    import warnings
    raw = {}
    try:
        import rainbow._binning as B
        import rainbow.agilent.chemstation  # noqa: F401  (the parsers that bin MS data)
        import rainbow.agilent.openlab  # noqa: F401
        import rainbow.waters.masslynx  # noqa: F401
        orig = B.bin_datapairs
    except (ImportError, AttributeError):  # another rainbow version: its own grid (nominal mass)
        return rb.read(path), raw

    def keep(keys, values, pair_counts, bin_width, display_precision=None, data_dtype=np.int64,
             labels_only=False):
        if bin_width != _SENTINEL_BIN or labels_only:
            return orig(keys, values, pair_counts, bin_width, display_precision=display_precision,
                        data_dtype=data_dtype, labels_only=labels_only)
        tag = -float(len(raw) + 1)
        raw[tag] = (np.asarray(keys, dtype=np.float64), np.asarray(values, dtype=np.float64),
                    np.asarray(pair_counts, dtype=np.int64))
        n = np.asarray(pair_counts).size
        return np.array([tag]), np.zeros((n, 1), dtype=data_dtype)
    with _RB_LOCK:
        mods = [m for n, m in list(sys.modules.items()) if n.startswith("rainbow.") and m is not None and
                getattr(m, "bin_datapairs", None) is orig and m is not B]
        for m in mods:
            m.bin_datapairs = keep
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                dd = rb.read(path, bin_width=_SENTINEL_BIN)
        finally:
            for m in mods:
                m.bin_datapairs = orig
    return dd, raw


def _ms_pairs(df, raw):
    """(mz, it, counts) of an MS DataFile: the kept pairs, or the matrix."""
    y = np.asarray(df.ylabels)
    if y.size == 1 and y.dtype.kind == "f" and float(y[0]) in raw:
        return raw[float(y[0])]
    d = np.asarray(df.data, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    nz = d != 0
    counts = nz.sum(axis=1)
    rows, cols = np.nonzero(nz)
    return y[cols], d[rows, cols], counts


def _agilent_signal_polarities(path):
    """{signal number: '+'/'-'} from the MS acquisition method text of a
    ChemStation .D (ACQ.M\\ACQ.TXT, acq.txt): '[Signal n]' ... 'Polarity : Positive'."""
    out = {}
    cands = []
    for root in (path, os.path.join(path, "ACQ.M")):
        try:
            cands += [os.path.join(root, n) for n in os.listdir(root) if n.lower() == "acq.txt"]
        except OSError:
            pass
    for p in cands:
        try:
            with open(p, "rb") as fh:
                b = fh.read(4 << 20)
        except OSError:
            continue
        txt = b.decode("utf-16", "ignore") if b[:2] in (b"\xff\xfe", b"\xfe\xff") or b[1:2] == b"\x00" \
            else b.decode("latin-1", "ignore")
        sig = None
        for line in txt.splitlines():
            m = re.match(r"\s*\[Signal\s+(\d+)\]", line)
            if m:
                sig = int(m.group(1))
                continue
            m = re.match(r"\s*Polarity\s*:\s*(Positive|Negative)", line, re.I)
            if m and sig is not None and sig not in out:
                out[sig] = "+" if m.group(1).lower() == "positive" else "-"
        if out:
            break
    return out


def _uv_wavelength(df):
    md = df.metadata or {}
    if md.get("wavelength"):
        return float(md["wavelength"])
    m = re.search(r"(?:Sig\s*=\s*|\b)(\d{3}(?:\.\d+)?)\s*(?:nm|,)", str(md.get("signal", "")))
    return float(m.group(1)) if m else None


def _pda_from_rainbow(dd, scale_units=None):
    """PDAData from the UV channels of a DataDirectory (the spectra with the
    most wavelengths, else the single wavelength signals); names of the other
    detector channels."""
    spectra, singles, other = [], [], []
    for df in dd.datafiles:
        det = df.detector
        if det == "MS":
            continue
        if det == "UV":
            y = np.asarray(df.ylabels)
            if y.size > 1:
                try:
                    w = y.astype(float)
                except (TypeError, ValueError):
                    other.append(df.name)
                    continue
                if np.all(np.diff(w) > 0) and 150 <= w[0] and w[-1] <= 1100:
                    spectra.append((df, w))
                else:
                    other.append(df.name)
                continue
            wl = _uv_wavelength(df)
            if wl is not None:
                singles.append((wl, df))
                continue
        other.append("%s (%s)" % (df.name, det) if det else df.name)
    pda = None
    if spectra:
        df, w = max(spectra, key=lambda s: len(s[1]))
        pda = build_pda(df.xlabels, w, df.data, scale_units or (df.metadata or {}).get("unit", "mAU"))
    elif singles:
        pda = pda_from_channels([(wl, df.xlabels, np.asarray(df.data, float)[:, 0],
                                  scale_units or (df.metadata or {}).get("unit", "mAU")) for wl, df in singles])
    return pda, other


def _base(path):
    """File name of a path written on another computer (any separator)."""
    return str(path).replace("/", "\\").rstrip("\\").split("\\")[-1]


def _date(text):
    """datetime of a vendor's date text, or None."""
    if not text:
        return None
    if isinstance(text, datetime.datetime):
        return text
    s = re.sub(r"\s+", " ", str(text).strip())
    s = re.sub(r"\s*[+-]\d{4}$", "", s)  # time zone
    s = re.sub(r"\s*[+-]\d{2}:\d{2}$", "", s)
    for f in ("%d %b %y %I:%M %p", "%d-%b-%y, %H:%M:%S", "%d-%b-%Y %H:%M:%S", "%d-%b-%Y, %H:%M:%S",
              "%d %b %Y %H:%M:%S", "%a %b %d %H:%M:%S %Y", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S",
              "%Y%m%d%H%M%S", "%d.%m.%Y %H:%M:%S", "%m/%d/%Y %I:%M:%S %p"):
        try:
            return datetime.datetime.strptime(s, f)
        except ValueError:
            continue
    m = re.match(r"(\d{4})-(\d{2})-(\d{2})[T ](\d{2}):(\d{2})(?::(\d{2}))?", s)
    if m:
        try:
            return datetime.datetime(*[int(v) for v in m.groups(default="0")])
        except ValueError:
            return None
    return None


def _waters_header(path):
    """Sample information of a Waters .raw folder (_HEADER.TXT)."""
    p = next((os.path.join(path, n) for n in os.listdir(path) if n.lower() == "_header.txt"), None)
    if p is None:
        return {}
    vals = {}
    with open(p, "r", encoding="latin-1", errors="replace") as fh:
        for line in fh:
            m = re.match(r"\$\$ ([^:]+):\s?(.*)", line.rstrip("\r\n"))
            if m:
                vals[m.group(1).strip()] = m.group(2).strip()
    info = {}
    for key, name in (("sample_name", "Acquired Name"), ("comment", "Sample Description"),
                      ("sample_id", "SampleID"), ("operator", "User Name"), ("vial", "Bottle Number"),
                      ("instrument", "Instrument"), ("method_file", "Inlet Method")):
        if vals.get(name):
            info[key] = vals[name]
    if info.get("method_file"):
        info["method_file"] = _base(info["method_file"])
    d = _date(("%s %s" % (vals.get("Acquired Date", ""), vals.get("Acquired Time", ""))).strip())
    if d:
        info["acquired"] = d
    return info


def _agilent_info(path, md=None):
    """Sample information of an Agilent .D folder (rainbow's metadata)."""
    if md is None:
        try:
            md = _rainbow().read_metadata(path)
        except Exception:
            md = {}
        md = md.get("metadata", md) if isinstance(md, dict) else {}
    info = {}
    for key, names in (("sample_name", ("sample", "notebook")), ("operator", ("operator",)),
                       ("vial", ("vialpos",)), ("method_file", ("acq_method", "method"))):
        for n in names:
            if md.get(n):
                info[key] = _base(md[n]) if key == "method_file" else str(md[n])
                break
    d = _date(md.get("date"))
    if d:
        info["acquired"] = d
    mods = md.get("modules") or []
    if mods:
        info["instrument"] = "Agilent " + ", ".join(str(m.get("model") or m.get("name")) for m in mods)
    else:
        try:
            ms = any(n.lower().endswith(".ms") for n in os.listdir(path)) if os.path.isdir(path) else False
        except OSError:
            ms = False
        tech = "GC" if str(md.get("technique", "")).upper() == "GC" else "LC"
        info["instrument"] = "Agilent " + ("%s-MS" % tech if ms else "HPLC")
    return info


def _props(info, instrument):
    """props (as Bruker files have them: the reports read these) of sample information."""
    p = {}
    for key, name in (("sample_name", "SampleName"), ("operator", "OperatorName"), ("method_file", "MethodName")):
        if info.get(key):
            p[name] = info[key]
    p["InstrumentName"] = info.get("instrument") or instrument
    if info.get("acquired"):
        p["AcquisitionDateTime"] = info["acquired"].strftime("%Y-%m-%dT%H:%M:%S")
    return p


def read_rainbow(path, kind, progress=None):
    dd, raw = _rainbow_read(path)
    waters = kind == "waters_raw"
    md = dict(dd.metadata or {})
    for df in dd.datafiles:  # sample name and method in the files' own headers
        for k in ("notebook", "method"):
            if (df.metadata or {}).get(k) and k not in md:
                md[k] = df.metadata[k]
    info = _waters_header(path) if waters else _agilent_info(path, md)
    pols = {} if waters else _agilent_signal_polarities(path)
    sc = Scans()
    ms_files = [df for df in dd.datafiles if df.detector == "MS"]
    for n, df in enumerate(ms_files):
        mz, it, counts = _ms_pairs(df, raw)
        md_f = df.metadata or {}
        if waters:
            pol = md_f.get("polarity")
            m = re.search(r"(\d+)", df.name)
            name = "Function %d" % int(m.group(1)) if m else df.name
        else:
            m = re.search(r"MSD(\d+)", df.name, re.I)
            pol = pols.get(int(m.group(1))) if m else None
            name = os.path.splitext(df.name)[0].upper()
        ions = np.unique(np.round(mz, 1))
        sim = ions.size <= 20  # (rainbow tags the channels it did not bin as SIM)
        if sim and ions.size <= 3:
            name += " SIM " + ", ".join("%g" % v for v in ions)
        ev = sc.add_event(pol if pol in ("+", "-") else None, name,
                          float(mz.min()) if len(mz) and not sim else None,
                          float(mz.max()) if len(mz) and not sim else None)
        sc.add_block(ev, df.xlabels, counts, mz, it)
        if progress is not None:
            progress(n + 1, len(ms_files) + 1)
    tech = "GC" if str(md.get("technique", "")).upper() == "GC" else "LC"
    inst = ("Waters" if waters else "Agilent") + (" %s-MS" % tech if ms_files else " HPLC")
    info.setdefault("instrument", inst)
    ms = build_ms(path, kind, sc, _props(info, inst))
    pda, other = _pda_from_rainbow(dd, "uAU" if waters else None)
    notes = ["not shown: " + ", ".join(other)] if other else []
    return ms, pda, notes, info


# ==========================================================================
# Thermo .raw (RawFileReader)
# ==========================================================================
def _thermo_dir():
    import importlib.util
    spec = importlib.util.find_spec("unidec")
    if spec is None or not spec.submodule_search_locations:
        raise RuntimeError("the Thermo reader (RawFileReader) is missing")
    return os.path.join(list(spec.submodule_search_locations)[0], "UniDecImporter", "Thermo")


def _thermo():
    """The RawFileReader classes (loaded once through pythonnet)."""
    with _NET_LOCK:
        if not _NET:
            import clr
            d = _thermo_dir()
            for n in ("ThermoFisher.CommonCore.Data", "ThermoFisher.CommonCore.RawFileReader"):
                clr.AddReference(os.path.join(d, n + ".dll"))
            from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter
            from ThermoFisher.CommonCore.Data.Business import Device
            _NET.update(adapter=RawFileReaderAdapter, Device=Device)
    return _NET


def _np(arr):
    """numpy copy of a .NET double array."""
    if arr is None:
        return np.zeros(0)
    try:
        return np.array(memoryview(arr), dtype=np.float64)
    except (TypeError, ValueError, BufferError):
        return np.fromiter(arr, dtype=np.float64, count=arr.Length)


def _net_date(d):
    try:
        return datetime.datetime(d.Year, d.Month, d.Day, d.Hour, d.Minute, d.Second)
    except Exception:
        return None


def _thermo_open(path):
    N = _thermo()
    r = N["adapter"].FileFactory(path)
    if r is None or r.IsError or not r.IsOpen:
        msg = r.FileError.ErrorMessage if r is not None and r.IsError else "cannot open"
        if r is not None:
            r.Dispose()
        raise RuntimeError("Thermo .raw: %s" % msg)
    return r


def _thermo_info(r):
    info = {}
    try:
        si = r.SampleInformation
        for key, v in (("sample_name", si.SampleName), ("sample_id", si.SampleId), ("vial", si.Vial),
                       ("comment", si.Comment), ("method_file", si.InstrumentMethodFile)):
            if v:
                info[key] = str(v)
        if info.get("method_file"):
            info["method_file"] = _base(info["method_file"])
        if si.InjectionVolume:
            info["inj_vol"] = "%g" % si.InjectionVolume
    except Exception:
        pass
    try:
        who = str(r.FileHeader.WhoCreatedId or "")
        if who:
            info["operator"] = who
    except Exception:
        pass
    d = _net_date(r.CreationDate)
    if d:
        info["acquired"] = d
    try:
        N = _thermo()
        if r.GetInstrumentCountOfType(N["Device"].MS):
            r.SelectInstrument(N["Device"].MS, 1)
            idata = r.GetInstrumentData()
            info["instrument"] = ("Thermo " + str(idata.Model or idata.Name or "")).strip()
    except Exception:
        pass
    return info


def read_thermo(path, progress=None):
    N = _thermo()
    Dev = N["Device"]
    r = _thermo_open(path)
    try:
        info = _thermo_info(r)
        sc = Scans()
        keys = {}
        inst = info.setdefault("instrument", "Thermo")
        if r.GetInstrumentCountOfType(Dev.MS):
            r.SelectInstrument(Dev.MS, 1)
            h = r.RunHeaderEx
            first, last = int(h.FirstSpectrum), int(h.LastSpectrum)
            rows = []
            for s in range(first, last + 1):
                f = r.GetFilterForScanNumber(s)
                if int(f.MSOrder) != 1:
                    sc.n_msms += 1
                    continue
                pol = "+" if int(f.Polarity) == 1 else ("-" if int(f.Polarity) == 0 else None)
                ana = str(f.MassAnalyzer).replace("MassAnalyzer", "")
                st = r.GetScanStatsForScanNumber(s)
                mz = it = None
                if not st.IsCentroidScan:
                    cs = r.GetCentroidStream(s, False)
                    if cs is not None and cs.Length > 0:
                        mz, it = _np(cs.Masses), _np(cs.Intensities)
                if mz is None:
                    seg = r.GetSegmentedScanFromScanNumber(s, st)
                    mz, it = _np(seg.Positions), _np(seg.Intensities)
                rows.append((pol, ana, r.RetentionTimeFromScanNumber(s), mz, it, float(st.TIC),
                             float(st.BasePeakIntensity), float(st.LowMass), float(st.HighMass)))
                if progress is not None and s % 200 == 0:
                    progress(s - first, last - first + 1)
            anas = set(r_[1] for r_ in rows)
            for pol, ana, rt, mz, it, tic, bpc, lo, hi in rows:
                key = (pol, ana)
                if key not in keys:
                    keys[key] = sc.add_event(pol, ana if len(anas) > 1 else "", lo, hi)
                e = sc.events[keys[key]]
                e["mz_low"], e["mz_high"] = min(e["mz_low"], lo), max(e["mz_high"], hi)
                sc.add(keys[key], rt, mz, it, tic, bpc)
        ms = build_ms(path, "thermo_raw", sc, _props(info, inst))
        pda, notes = _thermo_pda(r, Dev)
    finally:
        r.Dispose()
    return ms, pda, notes, info


def _thermo_units(r):
    try:
        u = str(r.GetInstrumentData().Units)
    except Exception:
        return "mAU"
    return {"AbsorbanceUnits": "AU", "MilliAbsorbanceUnits": "mAU", "MicroAbsorbanceUnits": "uAU"}.get(u, "mAU")


def _thermo_pda(r, Dev):
    """PDA spectra (or UV detector channels) of a Thermo file; names of the
    channels not shown."""
    notes = []
    try:
        if r.GetInstrumentCountOfType(Dev.Pda):
            r.SelectInstrument(Dev.Pda, 1)
            h = r.RunHeaderEx
            times, rows, wl = [], [], None
            for s in range(int(h.FirstSpectrum), int(h.LastSpectrum) + 1):
                seg = r.GetSegmentedScanFromScanNumber(s, None)
                w, a = _np(seg.Positions), _np(seg.Intensities)
                if wl is None:
                    wl = w
                elif len(w) != len(wl) or not np.allclose(w, wl):
                    a = np.interp(wl, w, a)
                times.append(r.RetentionTimeFromScanNumber(s))
                rows.append(a)
            if rows and wl is not None and len(wl):
                return build_pda(times, wl, np.vstack(rows), _thermo_units(r)), notes
        n_uv = r.GetInstrumentCountOfType(Dev.UV)
        chans = []
        for k in range(1, n_uv + 1):
            r.SelectInstrument(Dev.UV, k)
            h = r.RunHeaderEx
            try:
                labels = [str(x) for x in r.GetInstrumentData().ChannelLabels]
            except Exception:
                labels = []
            times, rows = [], []
            for s in range(int(h.FirstSpectrum), int(h.LastSpectrum) + 1):
                seg = r.GetSegmentedScanFromScanNumber(s, None)
                times.append(r.RetentionTimeFromScanNumber(s))
                rows.append(_np(seg.Intensities))
            if not rows:
                continue
            n = min(len(x) for x in rows)
            M = np.vstack([x[:n] for x in rows])
            for c in range(n):
                lab = labels[c] if c < len(labels) else "UV %d" % (c + 1)
                m = re.search(r"(\d{3}(?:\.\d+)?)\s*nm", lab)
                if m:
                    chans.append((float(m.group(1)), times, M[:, c], _thermo_units(r)))
                else:
                    notes.append(lab)
        if chans:
            return pda_from_channels(chans), (["not shown: " + ", ".join(notes)] if notes else [])
    except Exception as ex:
        notes.append("UV data: %s" % ex)
    return None, (["not shown: " + ", ".join(notes)] if notes else [])


# ==========================================================================
# mzXML (pyteomics) and ANDI netCDF (scipy)
# ==========================================================================
def read_mzxml(path, progress=None):
    from pyteomics import mzxml
    sc = Scans()
    keys = {}
    with mzxml.read(path, use_index=False, decode_binary=True) as rd:
        for k, s in enumerate(rd):
            if int(s.get("msLevel", 1)) != 1:
                sc.n_msms += 1
                continue
            pol = {"+": "+", "-": "-"}.get(str(s.get("polarity", "")).strip())
            if pol not in keys:
                keys[pol] = sc.add_event(pol, "")
            e = sc.events[keys[pol]]
            for key, f in (("mz_low", min), ("mz_high", max)):
                v = s.get("lowMz" if key == "mz_low" else "highMz")
                if v is not None:
                    e[key] = float(v) if e[key] is None else f(e[key], float(v))
            rt = float(s.get("retentionTime", 0.0))  # pyteomics: minutes
            tic = s.get("totIonCurrent")
            sc.add(keys[pol], rt, s.get("m/z array", []), s.get("intensity array", []),
                   float(tic) if tic else None)
            if progress is not None and k % 500 == 0:
                progress(k, 0)
    return build_ms(path, "mzxml", sc, {"InstrumentName": ""}), None, [], {}


def _cdf_text(v):
    if v is None:
        return ""
    if isinstance(v, bytes):
        v = v.decode("latin-1", "replace")
    return str(v).strip("\x00 ").strip()


def _andi_info(A):
    """Sample information of the global attributes of an ANDI file."""
    info = {}
    for key, names in (("sample_name", ("sample_name", "sample_id")), ("operator", ("operator_name",)),
                       ("method_file", ("experiment_title", "detection_method_name")),
                       ("instrument", ("instrument_name", "instrument_mfr")), ("comment", ("sample_comments",))):
        for n in names:
            if A.get(n):
                info[key] = A[n]
                break
    d = _date((A.get("injection_date_time_stamp") or A.get("experiment_date_time_stamp") or "")[:14])
    if d:
        info["acquired"] = d
    return info


def read_andi(path, progress=None):
    from scipy.io import netcdf_file
    f = netcdf_file(path, "r", mmap=False)
    try:
        A = {k: _cdf_text(v) for k, v in f._attributes.items()}
        V = f.variables
        info = _andi_info(A)
        if "mass_values" in V:  # ANDI/MS
            t = np.asarray(V["scan_acquisition_time"][:], float) / 60.0
            idx = np.asarray(V["scan_index"][:], np.int64)
            cnt = np.asarray(V["point_count"][:], np.int64)
            mz = np.asarray(V["mass_values"][:], float)
            it = np.asarray(V["intensity_values"][:], float)
            sf = V["mass_values"]
            if getattr(sf, "scale_factor", None):
                mz = mz * float(sf.scale_factor)
            sfi = V["intensity_values"]
            if getattr(sfi, "scale_factor", None):
                it = it * float(sfi.scale_factor)
            tot = np.asarray(V["total_intensity"][:], float) if "total_intensity" in V else None
            ptxt = A.get("test_ionization_polarity", "").lower()
            pol = "+" if ptxt.startswith("pos") else ("-" if ptxt.startswith("neg") else None)
            sc = Scans()
            ev = sc.add_event(pol, "")
            if len(idx) and np.array_equal(idx, np.concatenate([[0], np.cumsum(cnt)[:-1]])):
                n = int(cnt.sum())
                sc.add_block(ev, t, cnt, mz[:n], it[:n], tot)
            else:
                for k in range(len(t)):
                    a = int(idx[k])
                    sc.add(ev, t[k], mz[a:a + cnt[k]], it[a:a + cnt[k]], None if tot is None else tot[k])
            return build_ms(path, "andi_cdf", sc, _props(info, "")), None, [], info
        if "ordinate_values" in V:  # ANDI/CHROM
            y = np.asarray(V["ordinate_values"][:], float)
            dt = float(np.asarray(V["actual_sampling_interval"].data).ravel()[0]) if \
                "actual_sampling_interval" in V else 1.0
            t0 = float(np.asarray(V["actual_delay_time"].data).ravel()[0]) if "actual_delay_time" in V else 0.0
            unit = A.get("retention_unit", "seconds").lower()
            div = 1.0 if unit.startswith("min") else 60.0
            t = (t0 + np.arange(len(y)) * dt) / div
            if "raw_data_retention" in V and len(V["raw_data_retention"][:]) == len(y):
                t = np.asarray(V["raw_data_retention"][:], float) / div
            u = A.get("detector_unit", "mAU")
            names = " ".join(A.get(n, "") for n in ("detector_name", "detection_method_name", "detector_maximum_value"))
            m = re.search(r"(\d{3}(?:\.\d+)?)\s*nm", names) or re.search(r"Sig\s*=\s*(\d{3}(?:\.\d+)?)", names)
            ul = u.replace("µ", "u").lower()
            if not m and ul not in ("au", "mau", "uau"):
                return None, None, ["not shown: %s (%s)" % (A.get("detector_name") or "detector", u)], info
            wl = float(m.group(1)) if m else 0.0
            return None, pda_from_channels([(wl, t, y, u)]), [], info
        raise ValueError("no MS or chromatogram data in this netCDF file")
    finally:
        f.close()


# ==========================================================================
# entry points
# ==========================================================================
def read_mzml(path, progress=None):
    import hrms_data
    ms = hrms_data.open_hrms(path, progress)
    return ms, None, [], {}


READERS = {"agilent_chemstation": lambda p, pr: read_rainbow(p, "agilent_chemstation", pr),
           "waters_raw": lambda p, pr: read_rainbow(p, "waters_raw", pr),
           "thermo_raw": read_thermo, "mzml": read_mzml, "mzxml": read_mzxml, "andi_cdf": read_andi}


def read(path, progress=None):
    """{"ms": ..., "pda": ..., "info": sample information}, [notes and errors]
    for a file of another vendor (not .lcd)."""
    path = data_formats.data_path(path)
    window, kind = data_formats.detect(path)
    if kind == "sciex_wiff":
        raise RuntimeError(WIFF_TEXT)
    fn = READERS.get(kind)
    if fn is None:
        raise RuntimeError("not a data file LCMS Analysis can read")
    ms, pda, notes, info = fn(path, progress)
    return {"ms": ms, "pda": pda, "info": info or {}}, list(notes or [])


def read_sample_info(path):
    """Sample name, operator, date ... of a data file as lcms_data.read_sample_info
    gives them for .lcd files; {} when there is none. Reads only the headers."""
    try:
        window, kind = data_formats.detect(path)
        p = data_formats.data_path(path)
        if kind == "shimadzu_lcd":
            return lcms_data.read_sample_info(p)
        if kind == "waters_raw":
            return _waters_header(p)
        if kind == "agilent_chemstation":
            return _agilent_info(p)
        if kind == "andi_cdf":
            from scipy.io import netcdf_file
            f = netcdf_file(p, "r", mmap=True)
            try:
                return _andi_info({k: _cdf_text(v) for k, v in f._attributes.items()})
            finally:
                f.close()
        if kind == "thermo_raw":
            r = _thermo_open(p)
            try:
                return _thermo_info(r)
            finally:
                r.Dispose()
    except Exception:
        return {}
    return {}
