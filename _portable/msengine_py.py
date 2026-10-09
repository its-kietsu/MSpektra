"""
Bridge to msengine (the C++ core of MSpektra, _portable\\msengine\\msengine.dll;
see the cpp folder of the project for the sources).

EngineFile has the interface of hrms_data.BrukerD, so the HRMS window and
the reports use it unchanged; maxent() runs the maximum entropy
deconvolution of ms_deconv in C++. When the library is missing or fails to
load, available() is False and the Python code is used (hrms_data.open_hrms
and ms_deconv decide). The setting "engine": "python" in
config\\portable_settings.json turns the bridge off.
"""
import os
import ctypes
from ctypes import c_char_p, c_int, c_long, c_double, c_void_p, POINTER, byref

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_LIB = {"lib": None, "error": None}
PROGRESS = ctypes.CFUNCTYPE(c_int, c_void_p, c_long, c_long)


def _path():
    name = "msengine.dll" if os.name == "nt" else "msengine.so"
    for d in (os.path.join(HERE, "msengine"), HERE):
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return p
    return None


def _enabled():
    if os.environ.get("MSENGINE_DISABLE"):  # tests of the Python reference
        return False
    try:
        import unidec_theme as T
        return str(T._load().get("engine", "cpp")).lower() != "python"
    except Exception:
        return True


def lib():
    """The loaded library or None (the reason in _LIB["error"])."""
    if _LIB["lib"] is not None or _LIB["error"] is not None:
        return _LIB["lib"]
    try:
        if not _enabled():
            raise RuntimeError("turned off in the settings")
        p = _path()
        if p is None:
            raise RuntimeError("msengine library not found")
        # OpenMP workers must not spin wait: on a laptop with few cores the
        # interface thread and the SDK thread would be starved
        os.environ.setdefault("OMP_WAIT_POLICY", "PASSIVE")
        os.environ.setdefault("GOMP_SPINCOUNT", "0")
        L = ctypes.CDLL(p)
        L.ms_version.restype = c_char_p
        L.ms_last_error.restype = c_char_p
        L.ms_set_threads.argtypes = [c_int]
        L.ms_get_threads.restype = c_int
        if hasattr(L, "ms_memory_budget"):  # since 2.7
            L.ms_memory_budget.restype = c_double
        L.ms_open.argtypes = [c_char_p, c_char_p, PROGRESS, c_void_p]
        L.ms_open.restype = c_void_p
        L.ms_close.argtypes = [c_void_p]
        for f in ("ms_file_kind", "ms_file_summary", "ms_file_recalibration"):
            getattr(L, f).argtypes = [c_void_p]
            getattr(L, f).restype = c_char_p
        for f in ("ms_file_has_profile", "ms_file_n_events"):
            getattr(L, f).argtypes = [c_void_p]
            getattr(L, f).restype = c_int
        for f in ("ms_file_n_scans", "ms_file_n_msms"):
            getattr(L, f).argtypes = [c_void_p]
            getattr(L, f).restype = c_long
        L.ms_file_event_polarity.argtypes = [c_void_p, c_int]
        L.ms_file_event_polarity.restype = c_int
        L.ms_file_event_mz_range.argtypes = [c_void_p, c_int, POINTER(c_double), POINTER(c_double)]
        L.ms_file_event_mz_range.restype = c_int
        PD = POINTER(c_double)
        L.ms_event_scans.argtypes = [c_void_p, c_int, POINTER(PD), POINTER(PD), POINTER(PD)]
        L.ms_event_scans.restype = c_long
        L.ms_xic.argtypes = [c_void_p, c_int, c_double, c_double, POINTER(PD), POINTER(c_long)]
        L.ms_xic.restype = c_int
        L.ms_average.argtypes = [c_void_p, c_int, c_double, c_double, PD, c_int, POINTER(PD), POINTER(PD),
                                 POINTER(c_long), POINTER(c_long), PROGRESS, c_void_p]
        L.ms_average.restype = c_int
        L.ms_scan.argtypes = [c_void_p, c_int, c_double, POINTER(PD), POINTER(PD), POINTER(c_long), POINTER(c_long)]
        L.ms_scan.restype = c_int
        L.ms_maxent.argtypes = [PD, PD, c_long, POINTER(MaxentParams), POINTER(c_void_p), PROGRESS, c_void_p]
        L.ms_maxent.restype = c_int
        L.ms_maxent_free.argtypes = [c_void_p]
        for f in ("ms_maxent_mass", "ms_maxent_fit", "ms_maxent_zdist"):
            getattr(L, f).argtypes = [c_void_p, POINTER(PD), POINTER(PD)]
            getattr(L, f).restype = c_long
        L.ms_maxent_peaks.argtypes = [c_void_p, POINTER(PD), POINTER(PD), POINTER(PD), POINTER(PD), POINTER(POINTER(c_int))]
        L.ms_maxent_peaks.restype = c_long
        for f in ("ms_maxent_r2", "ms_maxent_chi2"):
            getattr(L, f).argtypes = [c_void_p]
            getattr(L, f).restype = c_double
        L.ms_maxent_rounds.argtypes = [c_void_p]
        L.ms_maxent_rounds.restype = c_int
        L.ms_maxent_notes.argtypes = [c_void_p]
        L.ms_maxent_notes.restype = c_char_p
        # phase 2
        L.ms_pda.argtypes = [c_void_p, POINTER(PD), POINTER(c_long), POINTER(PD), POINTER(c_long), POINTER(PD)]
        L.ms_pda.restype = c_int
        L.ms_sample_info.argtypes = [c_void_p]
        L.ms_sample_info.restype = c_char_p
        L.ms_file_set_bin_width.argtypes = [c_void_p, c_double]
        L.ms_file_set_bin_width.restype = c_int
        L.ms_unidec.argtypes = [PD, PD, c_long, c_char_p, c_char_p, c_char_p, c_char_p, POINTER(c_void_p), PROGRESS, c_void_p]
        L.ms_unidec.restype = c_int
        L.ms_unidec_free.argtypes = [c_void_p]
        for f in ("ms_unidec_mass", "ms_unidec_fit", "ms_unidec_zdist"):
            getattr(L, f).argtypes = [c_void_p, POINTER(PD), POINTER(PD)]
            getattr(L, f).restype = c_long
        L.ms_unidec_peaks.argtypes = [c_void_p, POINTER(PD), POINTER(PD), POINTER(PD), POINTER(PD), POINTER(PD),
                                      POINTER(POINTER(c_int))]
        L.ms_unidec_peaks.restype = c_long
        for f in ("ms_unidec_r2", "ms_unidec_uniscore"):
            getattr(L, f).argtypes = [c_void_p]
            getattr(L, f).restype = c_double
        if hasattr(L, "ms_unidec_quality"):  # since 2.8
            L.ms_unidec_quality.argtypes = [c_void_p, POINTER(c_int), POINTER(c_char_p)]
            L.ms_unidec_quality.restype = c_double
        for f in ("ms_unidec_notes", "ms_unidec_folder", "ms_isodec_notes"):
            getattr(L, f).argtypes = [c_void_p]
            getattr(L, f).restype = c_char_p
        L.ms_isodec.argtypes = [PD, PD, c_long, c_char_p, c_char_p, POINTER(c_void_p)]
        L.ms_isodec.restype = c_int
        L.ms_isodec_free.argtypes = [c_void_p]
        L.ms_isodec_peaks.argtypes = [c_void_p, POINTER(PD), POINTER(PD), POINTER(PD), POINTER(POINTER(c_int)),
                                      POINTER(POINTER(c_char_p))]
        L.ms_isodec_peaks.restype = c_long
        L.ms_isodec_mass.argtypes = [c_void_p, POINTER(PD), POINTER(PD)]
        L.ms_isodec_mass.restype = c_long
        if hasattr(L, "ms_event_flags"):  # since 2.8: saturated scans (Shimadzu)
            L.ms_event_flags.argtypes = [c_void_p, c_int, POINTER(POINTER(ctypes.c_ubyte))]
            L.ms_event_flags.restype = c_long
        if hasattr(L, "ms_open_arrays"):  # scans read in Python (lcms_sources)
            PI, PLL = POINTER(c_int), POINTER(ctypes.c_longlong)
            L.ms_open_arrays.argtypes = [c_char_p, c_char_p, c_long, PD, PI, PLL, PD, PD, PD, PD, c_int, PI, PD, PD,
                                         c_long]
            L.ms_open_arrays.restype = c_void_p
        _LIB["lib"] = L
        print("[engine] %s loaded from %s" % (L.ms_version().decode(), p))
    except Exception as ex:
        _LIB["error"] = str(ex)
        print("[engine] C++ core not used: %s" % ex)
    return _LIB["lib"]


def available():
    return lib() is not None


def _err(L):
    return (L.ms_last_error() or b"").decode("utf-8", "replace")


def _arr(ptr, n):
    """Copy of a library array (the library owns and reuses its buffers)."""
    if n <= 0 or not ptr:
        return np.zeros(0)
    return np.ctypeslib.as_array(ptr, shape=(int(n),)).astype(np.float64)


class MaxentParams(ctypes.Structure):
    _fields_ = [("z_lo", c_int), ("z_hi", c_int), ("mass_lo", c_double), ("mass_hi", c_double),
                ("mass_step", c_double), ("adduct_mass", c_double), ("mz_lo", c_double), ("mz_hi", c_double),
                ("min_intensity", c_double), ("min_intensity_pct", c_int), ("peak_width", c_double),
                ("resolution", c_double), ("resolved_isotopes", c_int), ("rounds", c_int), ("min_rounds", c_int),
                ("iterations", c_int), ("peak_window", c_double), ("peak_threshold", c_double),
                ("baseline", c_int), ("baseline_width", c_double),
                # since 2.8 (0 = automatic): envelope nodes across the mass range, noise model, chi squared target
                ("envelope_nodes", c_int), ("noise_model", c_int), ("chi2_target", c_double)]


# ============================================================================
# data files
# ============================================================================
class Cancelled(Exception):
    pass


class EngineFile(object):
    """A data file read by the C++ core; same interface as hrms_data.BrukerD."""

    def __init__(self, path, progress=None, handle=None):
        L = lib()
        if L is None:
            raise RuntimeError(_LIB["error"] or "msengine not available")
        self.lib = L
        self.path = path
        self.name = os.path.splitext(os.path.basename(path.rstrip("\\/")))[0]
        self._cb = PROGRESS(lambda u, i, n: 1 if (progress is None or progress(i, n) is not False) else 0)
        # handle: a file opened already (open_arrays)
        self.handle = handle or L.ms_open(path.encode("utf-8"), HERE.encode("utf-8"), self._cb, None)
        if not self.handle:
            raise RuntimeError(_err(L))
        self.kind = L.ms_file_kind(self.handle).decode("utf-8", "replace")
        self.binned = self.kind == "Shimadzu .lcd" or handle is not None  # sums on bins (bin width setting)
        self._summary = L.ms_file_summary(self.handle).decode("utf-8", "replace")
        self.recal_note = L.ms_file_recalibration(self.handle).decode("utf-8", "replace")
        self.has_profile = bool(L.ms_file_has_profile(self.handle))
        self.n_msms = int(L.ms_file_n_msms(self.handle))
        self.n_events = int(L.ms_file_n_events(self.handle))
        self.events, self._scans = [], []
        rts, tics, bpcs, sats = [], [], [], []
        PD = POINTER(c_double)
        for e in range(self.n_events):
            pol = L.ms_file_event_polarity(self.handle, e)
            lo, hi = c_double(), c_double()
            ok = L.ms_file_event_mz_range(self.handle, e, byref(lo), byref(hi)) == 0
            self.events.append({"polarity": "-" if pol < 0 else "+", "mz_low": lo.value if ok else None,
                                "mz_high": hi.value if ok else None})
            prt, ptic, pbpc = PD(), PD(), PD()
            n = L.ms_event_scans(self.handle, e, byref(prt), byref(ptic), byref(pbpc))
            rts.append(_arr(prt, n))
            tics.append(_arr(ptic, n))
            bpcs.append(_arr(pbpc, n))
            sat = np.zeros(max(int(n), 0), bool)  # scans with saturated readings left out (Shimadzu)
            if hasattr(L, "ms_event_flags") and n > 0:
                pf = POINTER(ctypes.c_ubyte)()
                nf = L.ms_event_flags(self.handle, e, byref(pf))
                if nf == n and pf:
                    sat = (np.ctypeslib.as_array(pf, shape=(int(n),)) & 1).astype(bool)
            sats.append(sat)
        # the scans of all events in one array in time order, as the Python readers
        # (hrms_data) have them: event_scans(e) indexes rt, tic, bpc, saturated, and
        # the scan index returned by scan_spectrum is a position in these arrays
        # (with polarity switching the events alternate; event by event in one array
        # put the time of another scan into the spectrum description)
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
        self.saturated = (np.concatenate(sats) if sats else np.zeros(0, bool))[order]
        self.calib = None
        self.calib_info = ""
        self.calib_events = None
        self.work = path
        # instrument, sample, operator, method and date (the reports use them), as the Python readers
        try:
            import hrms_data
            self.props = hrms_data.file_properties(path, self.kind) if handle is None else {}
        except Exception as ex:
            print("[engine] file properties of %s not read: %s" % (os.path.basename(path), ex))
            self.props = {}
        # mass chromatograms are quick (centroids cached in the library), except timsTOF
        # data (extracted by the SDK from the whole file each time)
        self.xic_fast = self.kind != "Bruker .d (timsTOF)"

    # ------------------------------------------------------------------ info
    def close(self):
        try:
            if getattr(self, "handle", 0):
                self.lib.ms_close(self.handle)
                self.handle = 0
        except Exception:
            pass

    def __del__(self):
        self.close()

    def event_scans(self, event):
        return self._scans[event]

    def event_label(self, event):
        labels = getattr(self, "labels", None)  # set by the reader of the file (lcms_sources)
        if labels:
            return labels[event]
        name = {"+": "Positive", "-": "Negative"}.get(self.events[event]["polarity"], "Scans")
        if self.kind == "Shimadzu .lcd" and self.n_events > 1:
            name += " (event %d)" % (event + 1)
        return name

    def _bins(self, binw):
        """Bin width of the Shimadzu sums (the LCMS window's setting)."""
        if not getattr(self, "binned", False):
            return
        w = float(binw) if binw else 0.05
        if abs(w - getattr(self, "_binw", 0.05)) < 1e-12:
            return
        if self.lib.ms_file_set_bin_width(self.handle, w):
            raise ValueError(_err(self.lib))
        self._binw = w

    def pda(self):
        """lcms_pda.PDAData of a Shimadzu file, or None."""
        L = self.lib
        PD = POINTER(c_double)
        wl, t, a, nw, nt = PD(), PD(), PD(), c_long(), c_long()
        rc = L.ms_pda(self.handle, byref(wl), byref(nw), byref(t), byref(nt), byref(a))
        if rc == -4:
            return None
        if rc:
            raise RuntimeError(_err(L))
        import lcms_pda
        obj = lcms_pda.PDAData.__new__(lcms_pda.PDAData)
        obj.wavelengths = _arr(wl, nw.value)
        obj.times = _arr(t, nt.value)
        obj.A = np.ascontiguousarray(_arr(a, nw.value * nt.value).reshape(nt.value, nw.value), dtype=np.float32)
        obj.interval_s = float(np.median(np.diff(obj.times)) * 60.0) if nt.value > 1 else 0.16
        obj.units = "mAU"
        return obj

    def sample_info(self):
        txt = (self.lib.ms_sample_info(self.handle) or b"").decode("utf-8", "replace")
        out = {}
        for line in txt.split("\n"):
            if "\t" in line:
                k, v = line.split("\t", 1)
                out[k] = v
        return out

    def adduct_sign(self, event):
        return -1 if self.events[event]["polarity"] == "-" else 1

    def set_calibration(self, func, info="", events=None):
        """m/z calibration (callable: measured m/z -> calibrated m/z) for the
        scan events given (None: every event). A calibration made on the
        ions of one polarity is not used for the other one, and it replaces
        only the calibration of its own events: with polarity switching each
        polarity keeps its own (the second one replaced the first). func None
        removes the calibration of the events given (None: of every event).
        calib, calib_info and calib_events describe the calibration set last
        (calibs: every one, oldest first, as (func, info, events or None))."""
        evs = None if events is None else set(int(e) for e in events)
        cals = []
        if evs is not None:  # the calibrations of the other events stay
            for f, i, ev in getattr(self, "calibs", None) or (
                    [(self.calib, self.calib_info, getattr(self, "calib_events", None))] if self.calib else []):
                left = (set(range(self.n_events)) if ev is None else set(ev)) - evs
                if left:
                    cals.append((f, i, left))
        if func is not None:
            cals.append((func, info, evs))
        self.calibs = cals
        f, i, ev = cals[-1] if cals else (None, "", None)
        self.calib, self.calib_info, self.calib_events = f, i, ev

    def calibration_for(self, event):
        """The calibration used for the spectra of this event, or None (also
        while calib is set to None: the m/z as recorded)."""
        if self.calib is None:
            return None
        cals = getattr(self, "calibs", None)
        if not cals:
            ev = getattr(self, "calib_events", None)
            return self.calib if (ev is None or int(event) in ev) else None
        for f, i, ev in reversed(cals):
            if ev is None or int(event) in ev:
                return f
        return None

    def _cal(self, mz, event=None):
        cal = self.calib if event is None else self.calibration_for(event)
        return cal(mz) if cal is not None else mz

    def _uncal(self, mz, cal=None):
        """Measured m/z that the calibration maps onto mz."""
        cal = cal or self.calib
        r = float(mz)
        for _ in range(4):
            r -= float(np.atleast_1d(cal(np.array([r])))[0]) - float(mz)
        return r

    def summary(self):
        return self._summary

    # ----------------------------------------------------------- chromatograms
    def chromatogram(self, event=0, kind="tic", mz=None, tol=0.01):
        idx = self.event_scans(event)
        t = self.rt[idx]
        if kind == "tic":
            return t, self.tic[idx]
        if kind == "bpc":
            return t, self.bpc[idx]
        if kind == "xic":
            # the window is in calibrated m/z: the m/z as recorded at its ends
            lo, hi = float(mz) - float(tol), float(mz) + float(tol)
            cal = self.calibration_for(event)
            if cal is not None:
                lo, hi = self._uncal(lo, cal), self._uncal(hi, cal)
            c, w = 0.5 * (lo + hi), 0.5 * abs(hi - lo)
            # timsTOF: every mass chromatogram is extracted from the whole file (seconds for a
            # long run), and the window draws them again for each change of the view: kept
            key = (int(event), round(c, 9), round(w, 9))
            cache = None if self.xic_fast else self.__dict__.setdefault("_xic_cache", {})
            if cache is not None and key in cache:
                return t, cache[key].copy()
            py, n = POINTER(c_double)(), c_long()
            if self.lib.ms_xic(self.handle, event, c, w, byref(py), byref(n)):
                raise RuntimeError(_err(self.lib))
            y = _arr(py, n.value)
            if cache is not None:
                if len(cache) >= 64:
                    cache.pop(next(iter(cache)))
                cache[key] = y.copy()
            return t, y
        raise ValueError(kind)

    # ---------------------------------------------------------------- spectra
    def average(self, event, t0, t1, bg=None, binw=None, fill=True, progress=None):
        L = self.lib
        self._bins(binw)
        ranges = []
        if bg:
            ranges = [bg] if isinstance(bg[0], (int, float)) else list(bg)
        flat = np.array([v for a, b in ranges for v in (a, b)], dtype=np.float64)
        pbg = flat.ctypes.data_as(POINTER(c_double)) if len(flat) else POINTER(c_double)()
        pm, pi, n, ns = POINTER(c_double)(), POINTER(c_double)(), c_long(), c_long()
        state = {"cancel": False}

        def cb(u, i, k):
            try:
                if progress is not None and progress(i, k) is False:
                    state["cancel"] = True
                    return 0
            except Exception as ex:  # an exception in a ctypes callback would be lost
                print("[engine] progress callback failed:", ex)
                return 1
            return 1
        rc = L.ms_average(self.handle, event, float(min(t0, t1)), float(max(t0, t1)), pbg, len(ranges),
                          byref(pm), byref(pi), byref(n), byref(ns), PROGRESS(cb), None)
        if rc == -2 or state["cancel"]:
            import hrms_data
            raise hrms_data.Cancelled()
        if rc:
            raise RuntimeError(_err(L))
        if n.value <= 0:
            return np.zeros((0, 2)), 0
        s = np.column_stack([_arr(pm, n.value), _arr(pi, n.value)])
        s[:, 0] = self._cal(s[:, 0], event)
        if not fill:
            s = s[s[:, 1] > 0]
        return s, int(ns.value)

    def scan_spectrum(self, event, t, binw=None, fill=True):
        """The scan of the event nearest to t; returns (spectrum, k) with k the
        position of that scan in rt (the library's own scan index counts the
        scans of every event in the file's order)."""
        self._bins(binw)
        idx = self.event_scans(event)
        k = int(idx[int(np.argmin(np.abs(self.rt[idx] - float(t))))]) if len(idx) else 0
        pm, pi, n, kk = POINTER(c_double)(), POINTER(c_double)(), c_long(), c_long()
        if self.lib.ms_scan(self.handle, event, float(t), byref(pm), byref(pi), byref(n), byref(kk)):
            raise RuntimeError(_err(self.lib))
        if n.value <= 0:
            return np.zeros((0, 2)), k
        s = np.column_stack([_arr(pm, n.value), _arr(pi, n.value)])
        s[:, 0] = self._cal(s[:, 0], event)
        if not fill:
            s = s[s[:, 1] > 0]
        return s, k


def open_file(path, progress=None):
    """EngineFile or None when the C++ core is not available."""
    if not available():
        return None
    return EngineFile(path, progress)


def open_arrays(path, kind, instrument, rt, event, offsets, mz, it, tic, bpc, ev_pol, ev_lo=None, ev_hi=None,
                n_msms=0):
    """EngineFile of scans read in Python (see ms_open_arrays in msengine.h), or
    None when the library has no ms_open_arrays."""
    L = lib()
    if L is None or not hasattr(L, "ms_open_arrays"):
        return None
    PD = POINTER(c_double)
    keep = []

    def d(a):
        a = np.ascontiguousarray(a, dtype=np.float64)
        keep.append(a)
        return a.ctypes.data_as(PD)

    def i32(a):
        a = np.ascontiguousarray(a, dtype=np.int32)
        keep.append(a)
        return a.ctypes.data_as(POINTER(c_int))
    off = np.ascontiguousarray(offsets, dtype=np.int64)
    n_ev = len(ev_pol)
    n = len(rt)
    if len(off) != n + 1 or not (len(event) == len(tic) == len(bpc) == n) or \
            (n and (min(len(mz), len(it)) < off[-1] or off[0] != 0)):
        raise ValueError("open_arrays: arrays of the wrong length")
    h = L.ms_open_arrays(kind.encode("utf-8"), (instrument or "").encode("utf-8"), len(rt), d(rt), i32(event),
                         off.ctypes.data_as(POINTER(ctypes.c_longlong)), d(mz), d(it), d(tic), d(bpc), n_ev,
                         i32(ev_pol), d(np.nan_to_num(np.asarray(ev_lo if ev_lo is not None else [0] * n_ev, float))),
                         d(np.nan_to_num(np.asarray(ev_hi if ev_hi is not None else [0] * n_ev, float))),
                         int(n_msms))
    if not h:
        raise RuntimeError(_err(L))
    return EngineFile(path, None, handle=h)


# ============================================================================
# maximum entropy
# ============================================================================
def maxent(spec, p, progress=None):
    """Native-only MaxEnt pipeline, including display-resolution processing."""
    from ms_maxent_native import run
    return run(spec, p, progress)


# ============================================================================
# UniDec engine and IsoDec through the C++ core
# ============================================================================
def _param_text(p, extra=()):
    lines = {"z_lo": p["z_range"][0], "z_hi": p["z_range"][1], "mass_lo": p["mass_range"][0],
             "mass_hi": p["mass_range"][1], "mass_step": p.get("mass_step", 1.0), "peak_width": p.get("peak_width", 0) or 0,
             "peak_window": p.get("peak_window", 10), "peak_threshold": p.get("peak_thresh", 0.1),
             "adduct_mass": p.get("adduct_mass", 0), "sign": p.get("sign", 1),
             "mz_lo": (p.get("mz_range") or (0, 0))[0], "mz_hi": (p.get("mz_range") or (0, 0))[1],
             "min_intensity": p.get("min_int", 0) or 0, "min_intensity_pct": 1 if p.get("min_int_unit") == "pct" else 0,
             "baseline": 1 if p.get("baseline") else 0, "baseline_width": p.get("baseline_width", 15),
             "isotopemode": p.get("isotopemode", 0), "psfun": p.get("psfun", 0), "poolflag": p.get("poolflag", 2),
             "binning": p.get("binning", "auto")}
    for k in ("numit", "zzsig", "psig", "beta", "msig", "smooth", "mzbins", "subbuff", "subtype", "peaknorm",
              "phaseres", "threads") + tuple(extra) + tuple(k for k in p if k.startswith("iso_")):
        if p.get(k) not in (None, "") and k not in lines:
            lines[k] = p[k]
    return "".join("%s %s\n" % (k, v) for k, v in lines.items())


def _unidec_dir():
    """The unidec package folder, without importing the package (slow)."""
    import sys
    cands = [os.path.join(os.path.dirname(sys.executable), "Lib", "site-packages", "unidec")]
    try:
        import importlib.util
        spec = importlib.util.find_spec("unidec")
        if spec and spec.origin:
            cands.append(os.path.dirname(spec.origin))
    except Exception:
        pass
    for d in cands:
        if os.path.isdir(d):
            return d
    return None


def _engine_exe():
    """unidec.exe of the UniDec package: only a library older than 3.1 (no
    ms_unidec_engine) starts it; since 3.1 the engine runs inside msengine
    and this path is passed for compatibility only ("" when the library has
    the engine)."""
    if unidec_in_library():
        return ""
    d = _unidec_dir()
    return os.path.join(d, "bin", "unidec.exe") if d else ""


def unidec_in_library():
    """True when the library runs UniDec's engine itself (3.1 and later):
    no unidec.exe is needed for ms_unidec."""
    L = lib()
    if L is None or not hasattr(L, "ms_unidec_engine"):
        return False
    L.ms_unidec_engine.restype = c_char_p
    L.ms_unidec_engine.argtypes = []
    return True


def unidec(spec, folder, name, p, progress=None):
    """ms_deconv.unidec_deconvolute through the C++ core (None if the core
    is not available). Since 3.1 the core runs UniDec's engine itself (all
    processor cores); an older core needs unidec.exe."""
    L = lib()
    if L is None:
        return None
    exe = _engine_exe()
    if not unidec_in_library() and not os.path.isfile(exe):
        return None
    spec = np.ascontiguousarray(np.asarray(spec, np.float64))
    mz = np.ascontiguousarray(spec[:, 0])
    it = np.ascontiguousarray(spec[:, 1])
    import ms_deconv
    p = dict(p)
    p["adduct_mass"] = ms_deconv.carrier_mass(p)
    txt = _param_text(p)
    stage = {"last": 0}

    in_lib = unidec_in_library()

    def cb(u, i, n):
        if progress is not None and i != stage["last"]:
            stage["last"] = i
            try:
                if in_lib:   # 0 to 100: 10 to 90 while the engine iterates
                    k = int(i)
                    text = ("preparing the data" if k < 5 else "processing the data" if k < 10 else
                            "deconvoluting (UniDec engine in the library, all processor cores)" if k < 90 else
                            "picking peaks" if k < 100 else "done")
                else:
                    text = ["preparing the data", "processing the data",
                            "deconvoluting (UniDec engine, all processor cores)", "picking peaks",
                            "done"][min(max(int(i) - 1, 0), 4)]
                progress(int(i), int(n), text)
            except Exception:
                pass
        return 1
    out = c_void_p()
    rc = L.ms_unidec(mz.ctypes.data_as(POINTER(c_double)), it.ctypes.data_as(POINTER(c_double)), len(mz),
                     exe.encode("utf-8"), folder.encode("utf-8"), name.encode("utf-8"), txt.encode("utf-8"),
                     byref(out), PROGRESS(cb), None)
    if rc == -3:
        raise ValueError(_err(L))
    if rc:
        # "The UniDec engine ..." (it stopped or crashed): ms_deconv.run_method then
        # reports it instead of running the same engine again in Python
        raise RuntimeError(_err(L))
    try:
        PD = POINTER(c_double)
        a, b = PD(), PD()
        n = L.ms_unidec_mass(out, byref(a), byref(b))
        mass = np.column_stack([_arr(a, n), _arr(b, n)]) if n else np.zeros((0, 2))
        n = L.ms_unidec_fit(out, byref(a), byref(b))
        fit = np.column_stack([_arr(a, n), _arr(b, n)]) if n else None
        n = L.ms_unidec_zdist(out, byref(a), byref(b))
        zdist = np.column_stack([_arr(a, n), _arr(b, n)]) if n else None
        pm, ph, pa, ps, px, pn = PD(), PD(), PD(), PD(), PD(), POINTER(c_int)()
        n = L.ms_unidec_peaks(out, byref(pm), byref(ph), byref(pa), byref(ps), byref(px), byref(pn))
        peaks = []
        for k in range(int(n)):
            d = {"mass": float(pm[k]), "height": float(ph[k]), "area": float(pa[k]), "score": float(ps[k])}
            if pn and int(pn[k]) > 0:  # grouped isotopes: average mass, most abundant isotope
                d["apex"] = float(px[k])
                d["n_iso"] = int(pn[k])
                d["z"] = "most abundant %.3f" % d["apex"]
            peaks.append(d)
        ms_deconv._rel_frac(peaks)
        step = float(p.get("mass_step", 1.0))
        notes = (L.ms_unidec_notes(out) or b"").decode("utf-8", "replace") + "; C++ core"
        if not in_lib:
            notes += " with unidec.exe"
        iso_pk = ms_deconv.isotope_peaks(mass) if step < 0.1 and len(mass) else None
        res = {"method": "UniDec", "mass": ms_deconv._thin_mass_axis(mass), "fit": fit, "data": spec, "peaks": peaks,
               "r2": float(L.ms_unidec_r2(out)), "uniscore": float(L.ms_unidec_uniscore(out)), "zdist": zdist,
               "notes": notes, "folder": (L.ms_unidec_folder(out) or b"").decode("utf-8", "replace"),
               "isotope_peaks": iso_pk, "engine": "C++ (engine in the library)" if in_lib else "C++"}
        if hasattr(L, "ms_unidec_quality"):
            lv, tx = c_int(2), c_char_p()
            sc = float(L.ms_unidec_quality(out, byref(lv), byref(tx)))
            res["quality"] = {"score": sc, "level": ("good", "fair", "poor")[min(max(lv.value, 0), 2)],
                              "text": (tx.value or b"").decode("utf-8", "replace")}
        return res
    finally:
        L.ms_unidec_free(out)


def _isodec_dir():
    d = _unidec_dir()
    if d:
        d = os.path.join(d, "IsoDec")
        if os.path.isfile(os.path.join(d, "isodeclib.dll")) or os.path.isfile(os.path.join(d, "isodeclib.so")):
            return d
    return None


def isodec(spec, p):
    """ms_deconv.isodec_deconvolute through the C++ core (None if not available)."""
    L = lib()
    d = _isodec_dir()
    if L is None or d is None:
        return None
    spec = np.ascontiguousarray(np.asarray(spec, np.float64))
    mz = np.ascontiguousarray(spec[:, 0])
    it = np.ascontiguousarray(spec[:, 1])
    import ms_deconv
    q = dict(p)
    q["adduct_mass"] = ms_deconv.carrier_mass(q)
    txt = _param_text(q)
    out = c_void_p()
    rc = L.ms_isodec(mz.ctypes.data_as(POINTER(c_double)), it.ctypes.data_as(POINTER(c_double)), len(mz),
                     d.encode("utf-8"), txt.encode("utf-8"), byref(out))
    if rc == -3:
        raise ValueError(_err(L))
    if rc:
        raise RuntimeError(_err(L))
    try:
        PD = POINTER(c_double)
        pm, ph, px, pz, pc = PD(), PD(), PD(), POINTER(c_int)(), POINTER(c_char_p)()
        n = L.ms_isodec_peaks(out, byref(pm), byref(ph), byref(px), byref(pz), byref(pc))
        peaks = []
        for k in range(int(n)):
            zs = (pc[k] or b"").decode("utf-8", "replace") if pc else str(int(pz[k]))
            peaks.append({"mass": float(pm[k]), "height": float(ph[k]), "area": float(ph[k]), "z": zs,
                          "mz": float(px[k]), "avg": 0.0})
        ms_deconv._rel_frac(peaks)
        a, b = PD(), PD()
        n = L.ms_isodec_mass(out, byref(a), byref(b))
        mass = np.column_stack([_arr(a, n), _arr(b, n)]) if n else np.zeros((0, 2))
        notes = (L.ms_isodec_notes(out) or b"").decode("utf-8", "replace") + "; C++ core"
        return {"method": "IsoDec", "mass": mass, "fit": None, "data": spec, "peaks": peaks, "r2": None,
                "zdist": None, "notes": notes, "engine": "C++"}
    finally:
        L.ms_isodec_free(out)


# ---- 3.1 additions (agent P): the remaining computations of the windows ----------------------
# Each wrapper below runs one function of the Python modules in the library (hrms_calib, ms_formula,
# lcms_integrate, lcms_pda, ms_deconv: see include/msengine.h). It returns None when the library or
# that function (an older DLL) is not there, or when the library refuses the input: the module then
# runs its own Python code, the reference, which also raises the same errors as before.
_P31 = {"lib": None, "done": False}
_PD = POINTER(c_double)
# sum() of floats compensates the rounding since Python 3.12 (ms_formula's sums: the library does as the caller)
_PYSUM = 1 if __import__("sys").version_info >= (3, 12) else 0


class LcPeak2(ctypes.Structure):
    _fields_ = [(k, c_double) for k in ("t0", "t1", "rt", "height", "area", "base0", "base1", "apex_y")] + \
               [("i0", c_long), ("i1", c_long)]


class CmpOpts(ctypes.Structure):
    _fields_ = [("smooth", c_int), ("shift", c_double), ("t0", c_double), ("t1", c_double), ("baseline", c_int),
                ("rolling_min", c_double), ("scale", c_int), ("ref_t", c_double), ("ref_win", c_double)]


def _declare31(L):
    PD, PL, PI = _PD, POINTER(c_long), POINTER(c_int)
    PPD, PF = POINTER(PD), POINTER(ctypes.c_float)
    sig = {
        "ms_calibration_apply": ([c_int, c_int, PD, c_int, c_double, c_double, PD, PD, c_long], c_int),
        "ms_calibration_fit2": ([c_int, c_int, PD, PD, c_long, PD, PI, PD, PD, PD, PD], c_int),
        "ms_peak_position": ([PD, PD, c_long, c_long, c_double, c_double, c_int, PD, PD], c_int),
        "ms_peak_near": ([PD, PD, c_long, c_double, c_double, c_int, c_double, PD, PD], c_int),
        "ms_profile_apex": ([PD, PD, c_long, c_long], c_double),
        "ms_find_calibrants2": ([PD, PD, c_long, c_long, c_int, PD, c_long, c_double, c_double, PD, PD, PD, PI], c_int),
        "ms_isotope_pattern2": ([c_char_p, c_double, c_double, c_int, PPD, PPD, PL], c_int),
        "ms_find_formulas2": ([c_double, c_double, c_double, c_int, c_double, c_char_p, PL, PL, c_char_p, c_char_p,
                               c_int, c_int, POINTER(PL), PPD, PPD, PPD, PL, PI], c_int),
        "ms_isotope_match": ([PD, PD, c_long, PD, PD, c_long], c_double),
        "ms_measured_pattern": ([PD, PD, c_long, c_long, c_double, c_int, c_int, c_double, PD, PD], c_int),
        "ms_smooth": ([PD, c_long, c_int, PD], c_int),
        "ms_lc_integrate": ([PD, PD, c_long, c_double, c_double, c_int, c_int, c_double, c_double, c_double,
                             POINTER(POINTER(LcPeak2)), PL], c_int),
        "ms_lc_manual": ([PD, PD, c_long, c_double, c_double, POINTER(LcPeak2)], c_int),
        "ms_lc_peak_at": ([PD, PD, c_long, c_double, c_double, c_double, POINTER(LcPeak2)], c_int),
        "ms_lc_split": ([PD, PD, c_long, c_long, c_long, c_double, c_double, c_double, POINTER(LcPeak2)], c_int),
        "ms_pda_chromatogram": ([PF, c_long, c_long, PD, c_double, c_double, PF], c_int),
        "ms_pda_max_plot": ([PF, c_long, c_long, PF], c_int),
        "ms_pda_spectrum": ([PF, c_long, c_long, PD, c_double, c_double, c_double, c_int, c_double, c_double, PD], c_int),
        "ms_isotope_peaks": ([PD, PD, c_long, c_double, c_long, PPD, PPD, PL], c_int),
        "ms_refine_peak_mass": ([PD, PD, c_long, c_double], c_double),
        "ms_refine_species": ([PD, PD, c_long, PD, PD, c_long, c_double, c_int, c_int, c_double, c_double, PD], c_int),
        "ms_group_isotopes2": ([PD, PD, c_long, c_double, c_double, PPD, PPD, PPD, PPD, POINTER(PI), PPD, PPD, PL], c_int),
        "ms_restrict": ([PD, PD, c_long, c_double, c_double, PPD, PPD, PL], c_int),
        "ms_peak_mask": ([PD, PD, c_long, c_double, POINTER(ctypes.c_ubyte), PL], c_int),
        "ms_subtract_baseline": ([PD, PD, c_long, c_double, c_double, PD], c_int),
        "ms_label_maxima": ([PD, PD, c_long, c_int, c_double, c_double, PL, PL], c_int),
        "ms_xic_multi": ([c_void_p, c_int, PD, PD, c_int, PD, c_long], c_int),
        # 3.31: the Compare view (lcms_compare_core)
        "ms_cmp_subtract_blank": ([PD, PD, c_long, PD, PD, c_long, PD], c_int),
        "ms_cmp_find_apex": ([PD, PD, c_long, c_double, c_double, c_int, PD], c_int),
        "ms_cmp_baseline": ([PD, PD, c_long, c_int, c_double, PD], c_int),
        "ms_cmp_process": ([PD, PD, c_long, POINTER(CmpOpts), PD, PD, PL, PD, PD, PI], c_int),
        "ms_cmp_stack": ([PD, PD, PI, c_long, c_double, c_int, c_double, c_int, PD, PD], c_int),
        "ms_cmp_peaks": ([PD, PD, c_long, c_double, c_double, PD, c_long, POINTER(POINTER(LcPeak2)), PL, PL, PD,
                          PD], c_int),
        "ms_cmp_lambda_max": ([PF, c_long, c_long, PD, PD, c_double, c_double, c_double, c_int, c_double, c_double,
                               c_double, PD, PD, PI], c_int),
        # 3.65: area of a region (lcms_compare_core._region_peak)
        "ms_cmp_region": ([PD, PD, c_long, c_double, c_double, PD, PD, PD, c_long, PL, PD], c_int),
        # 3.33: suggested neutral mass from the ESI+ and ESI- spectra (ms_adducts)
        "ms_neutral_masses": ([PD, PD, c_long, PD, PD, c_long, c_double, c_double, c_int, c_int,
                               POINTER(POINTER(MassIon)), PL, POINTER(POINTER(MassCand)), PL], c_int),
        "ms_adduct_name": ([c_int], c_char_p),
        "ms_neutral_masses2": ([PD, PD, c_long, PD, PD, c_long, c_double, c_double, c_int, c_int,
                                POINTER(POINTER(MassIon)), PL, POINTER(POINTER(MassCand2)), PL, PD, PI], c_int),
        "ms_additive_name": ([c_int], c_char_p),
        # 3.35: the mass shift finder (ms_shifts)
        "ms_mass_shifts": ([PD, PD, PD, c_long, PD, c_long, PD, c_long, c_double, c_double, c_long, c_int, c_int,
                            c_int, c_double, POINTER(POINTER(ShiftPair)), PL, POINTER(POINTER(ShiftSpecies)),
                            POINTER(PI), PPD, PL], c_int),
    }
    for name, (args, res) in sig.items():
        if hasattr(L, name):
            f = getattr(L, name)
            f.argtypes = args
            f.restype = res


def _lib31(name):
    """The library when it has the function `name` (an older DLL may lack it), else None."""
    if not _P31["done"]:
        _P31["done"] = True
        L = lib()
        if L is not None:
            try:
                _declare31(L)
                _P31["lib"] = L
            except Exception as ex:
                print("[engine] 3.1 functions not used: %s" % ex)
    L = _P31["lib"]
    return L if L is not None and hasattr(L, name) else None


def _d(a):
    """Contiguous float64 copy (or the array itself) and its pointer."""
    a = np.ascontiguousarray(a, dtype=np.float64)
    return a, a.ctypes.data_as(_PD)


def _xy(x, y):
    """Two columns for the library: the columns of one n x 2 float64 array as they are (stride 2, no
    copy), else contiguous copies (stride 1). Returns the arrays (keep them while the call runs), their
    pointers and the stride."""
    x = np.asarray(x)
    y = np.asarray(y)
    if (x.dtype == np.float64 and y.dtype == np.float64 and x.ndim == 1 and y.ndim == 1 and len(x) == len(y)
            and len(x) > 0 and x.strides == (16,) and y.strides == (16,) and y.ctypes.data == x.ctypes.data + 8):
        return x, y, x.ctypes.data_as(_PD), y.ctypes.data_as(_PD), 2
    a, pa = _d(x)
    b, pb = _d(y)
    return a, b, pa, pb, 1


def _out(n):
    a = np.empty(int(n), dtype=np.float64)
    return a, a.ctypes.data_as(_PD)


# ------------------------------------------------------------------ hrms_calib
def calibration_fit(model_c, order, measured, reference):
    """ms_calibration_fit2: (coef (x0, s, model coefficients), err_ppm, rms, cv_err or None, cv_rms or None)."""
    L = _lib31("ms_calibration_fit2")
    if L is None:
        return None
    m, pm = _d(measured)
    r, pr = _d(reference)
    n = len(m)
    coef, pc = _out(16)
    nc = c_int(0)
    err, pe = _out(n)
    cv, pcv = _out(n)
    rms, cvr = c_double(), c_double()
    if L.ms_calibration_fit2(int(model_c), int(order), pm, pr, n, pc, byref(nc), pe, byref(rms), pcv, byref(cvr)):
        return None
    cv_rms = float(cvr.value)
    if not np.isfinite(cv_rms):
        return coef[:nc.value].copy(), err, float(rms.value), None, None
    return coef[:nc.value].copy(), err, float(rms.value), cv, cv_rms


def calibration_apply(model_c, coef, lo, hi, mz):
    """ms_calibration_apply on an array of m/z (None: not available)."""
    L = _lib31("ms_calibration_apply")
    if L is None:
        return None
    mz = np.asarray(mz, dtype=np.float64)
    shape = mz.shape
    x, px = _d(mz.ravel())
    c, pc = _d(coef)
    out, po = _out(len(x))
    order = len(c) - 6 if int(model_c) == 4 else 0
    if L.ms_calibration_apply(int(model_c), order, pc, len(c), float(lo), float(hi), px, po, len(x)):
        return None
    out = out.reshape(shape)
    return out[()] if out.ndim == 0 else out   # a numpy scalar for a scalar, as numpy gives


def peak_position(mz, it, lo, hi, profile):
    L = _lib31("ms_peak_position")
    if L is None:
        return None
    x, y, px, py, st = _xy(mz, it)
    pos, h = c_double(), c_double()
    rc = L.ms_peak_position(px, py, len(x), st, float(lo), float(hi), 1 if profile else 0, byref(pos), byref(h))
    if rc == -4:
        return None, 0.0
    return (float(pos.value), float(h.value)) if rc == 0 else None


def peak_near(mz, it, target, d, profile, min_rel):
    L = _lib31("ms_peak_near")
    if L is None:
        return None
    x, px = _d(mz)
    y, py = _d(it)
    pos, h = c_double(), c_double()
    rc = L.ms_peak_near(px, py, len(x), float(target), float(d), 1 if profile else 0, float(min_rel), byref(pos), byref(h))
    if rc == -4:
        return None, 0.0
    return (float(pos.value), float(h.value)) if rc == 0 else None


def profile_apex(mz, it, i):
    L = _lib31("ms_profile_apex")
    if L is None:
        return None
    x, px = _d(mz)
    y, py = _d(it)
    return float(L.ms_profile_apex(px, py, len(x), int(i)))


def find_calibrants(mz, it, refs, tol_ppm, min_rel, profile):
    """ms_find_calibrants2: (measured, height, rel %, found) arrays for the references."""
    L = _lib31("ms_find_calibrants2")
    if L is None:
        return None
    x, y, px, py, st = _xy(mz, it)
    r, pr = _d(refs)
    n = len(r)
    m, pm = _out(n)
    h, ph = _out(n)
    rel, prel = _out(n)
    found = np.zeros(n, dtype=np.int32)
    if L.ms_find_calibrants2(px, py, len(x), st, 1 if profile else 0, pr, n, float(tol_ppm), float(min_rel), pm, ph, prel,
                             found.ctypes.data_as(POINTER(c_int))):
        return None
    return m, h, rel, found


# ------------------------------------------------------------------ ms_formula
def isotope_pattern(formula_text, merge, min_rel):
    """ms_isotope_pattern2: (neutral masses, relative intensities %) or None."""
    L = _lib31("ms_isotope_pattern2")
    if L is None:
        return None
    a, b, n = _PD(), _PD(), c_long()
    if L.ms_isotope_pattern2(formula_text.encode("ascii"), float(merge), float(min_rel), _PYSUM, byref(a), byref(b),
                             byref(n)):
        return None
    return _arr(a, n.value), _arr(b, n.value)


def find_formulas(mz, neutral, tol, zabs, electrons, elements, lo, hi, add, rem, ion_itself):
    """ms_find_formulas2: (counts (n x ncols, int), ion m/z, error ppm, RDB) or None."""
    L = _lib31("ms_find_formulas2")
    if L is None:
        return None
    lo_a = (c_long * len(lo))(*[int(v) for v in lo])
    hi_a = (c_long * len(hi))(*[int(v) for v in hi])
    pc, pt, pe, pr = POINTER(c_long)(), _PD(), _PD(), _PD()
    n, nc = c_long(), c_int()
    rc = L.ms_find_formulas2(float(mz), float(neutral), float(tol), int(zabs), float(electrons),
                             " ".join(elements).encode("ascii"), lo_a, hi_a, (add or "").encode("ascii"),
                             (rem or "").encode("ascii"), 1 if ion_itself else 0, _PYSUM, byref(pc), byref(pt), byref(pe),
                             byref(pr), byref(n), byref(nc))
    if rc:
        return None
    k, w = int(n.value), int(nc.value)
    counts = np.ctypeslib.as_array(pc, shape=(k * w,)).astype(np.int64).reshape(k, w) if k else np.zeros((0, w), np.int64)
    return counts, _arr(pt, k), _arr(pe, k), _arr(pr, k)


def isotope_match(pattern, observed):
    L = _lib31("ms_isotope_match")
    if L is None:
        return None
    pat = np.asarray(pattern, dtype=np.float64)
    obs = np.asarray(observed, dtype=np.float64)
    if pat.ndim != 2 or obs.ndim != 2 or pat.shape[1] < 2 or obs.shape[1] < 2:
        return None
    a, pa = _d(pat[:, 0])
    b, pb = _d(pat[:, 1])
    c, pc = _d(obs[:, 0])
    e, pe = _d(obs[:, 1])
    return float(L.ms_isotope_match(pa, pb, len(a), pc, pe, len(c)))


def measured_pattern(spec, mono_mz, z, n, ppm):
    L = _lib31("ms_measured_pattern")
    if L is None:
        return None
    spec = np.asarray(spec, dtype=np.float64)
    x, y, px, py, st = _xy(spec[:, 0], spec[:, 1])
    om, pom = _out(n)
    oh, poh = _out(n)
    if L.ms_measured_pattern(px, py, len(x), st, float(mono_mz), int(z), int(n), float(ppm), pom, poh):
        return None
    return [(float(a), float(b)) for a, b in zip(om, oh)]


# ------------------------------------------------------------------ lcms_integrate
def _peak_dict(p):
    return {"rt": float(p.rt), "t0": float(p.t0), "t1": float(p.t1), "i0": int(p.i0), "i1": int(p.i1),
            "height": float(p.height), "area": float(p.area), "b0": float(p.base0), "b1": float(p.base1),
            "apex_y": float(p.apex_y)}


def smooth(y, points):
    L = _lib31("ms_smooth")
    if L is None:
        return None
    a, pa = _d(y)
    out, po = _out(len(a))
    if L.ms_smooth(pa, len(a), int(points), po):
        return None
    return out


def lc_integrate(t, y, threshold_pct, min_width_s, smooth_points, t_range, baseline_window_min):
    L = _lib31("ms_lc_integrate")
    if L is None:
        return None
    a, pa = _d(t)
    b, pb = _d(y)
    pk, n = POINTER(LcPeak2)(), c_long()
    lo, hi = (min(t_range), max(t_range)) if t_range is not None else (0.0, 0.0)
    if L.ms_lc_integrate(pa, pb, len(a), float(threshold_pct), float(min_width_s), int(smooth_points or 0),
                         1 if t_range is not None else 0, float(lo), float(hi), float(baseline_window_min),
                         byref(pk), byref(n)):
        return None
    return [_peak_dict(pk[k]) for k in range(int(n.value))]


def lc_manual(t, y, t0, t1):
    L = _lib31("ms_lc_manual")
    if L is None:
        return None
    a, pa = _d(t)
    b, pb = _d(y)
    p = LcPeak2()
    rc = L.ms_lc_manual(pa, pb, len(a), float(t0), float(t1), byref(p))
    if rc == -4:
        return False
    return _peak_dict(p) if rc == 0 else None


def lc_peak_at(t, y, x, search_s, baseline_window_min):
    L = _lib31("ms_lc_peak_at")
    if L is None:
        return None
    a, pa = _d(t)
    b, pb = _d(y)
    p = LcPeak2()
    rc = L.ms_lc_peak_at(pa, pb, len(a), float(x), float(search_s), float(baseline_window_min), byref(p))
    if rc == -4:
        return False
    return _peak_dict(p) if rc == 0 else None


def lc_split(t, y, i0, i1, b0, b1, x):
    L = _lib31("ms_lc_split")
    if L is None:
        return None
    a, pa = _d(t)
    b, pb = _d(y)
    out = (LcPeak2 * 2)()
    rc = L.ms_lc_split(pa, pb, len(a), int(i0), int(i1), float(b0), float(b1), float(x), out)
    if rc == -4:
        return False
    return (_peak_dict(out[0]), _peak_dict(out[1])) if rc == 0 else None


# ------------------------------------------------------------------ lcms_pda
def _f32(A):
    A = np.ascontiguousarray(A, dtype=np.float32)
    return A, A.ctypes.data_as(POINTER(ctypes.c_float))


def pda_chromatogram(A, wavelengths, wl, bw):
    L = _lib31("ms_pda_chromatogram")
    if L is None or np.ndim(A) != 2 or not np.size(A):
        return None
    A, pA = _f32(A)
    w, pw = _d(wavelengths)
    if len(w) != A.shape[1]:
        return None
    out = np.empty(A.shape[0], dtype=np.float32)
    if L.ms_pda_chromatogram(pA, A.shape[0], A.shape[1], pw, float(wl), float(bw),
                             out.ctypes.data_as(POINTER(ctypes.c_float))):
        return None
    return out


def pda_max_plot(A):
    L = _lib31("ms_pda_max_plot")
    if L is None or np.ndim(A) != 2 or not np.size(A):
        return None
    A, pA = _f32(A)
    out = np.empty(A.shape[0], dtype=np.float32)
    if L.ms_pda_max_plot(pA, A.shape[0], A.shape[1], out.ctypes.data_as(POINTER(ctypes.c_float))):
        return None
    return out


def pda_spectrum(A, times, interval_s, t0, t1, bg):
    """bg: None, a time, or (b0, b1)."""
    L = _lib31("ms_pda_spectrum")
    if L is None or np.ndim(A) != 2 or not np.size(A):
        return None
    A, pA = _f32(A)
    t, pt = _d(times)
    if len(t) != A.shape[0]:
        return None
    if bg is None:
        mode, b0, b1 = 0, 0.0, 0.0
    elif isinstance(bg, (int, float)):
        mode, b0, b1 = 1, float(bg), 0.0
    else:
        mode, b0, b1 = 2, float(bg[0]), float(bg[1])
    out, po = _out(A.shape[1])
    if L.ms_pda_spectrum(pA, A.shape[0], A.shape[1], pt, float(interval_s), float(t0),
                         float("nan") if t1 is None else float(t1), mode, b0, b1, po):
        return None
    return out


# ------------------------------------------------------------------ ms_deconv helpers
def isotope_peaks(mass, threshold, max_n):
    L = _lib31("ms_isotope_peaks")
    if L is None:
        return None
    mass = np.asarray(mass, dtype=np.float64)
    x, px = _d(mass[:, 0])
    y, py = _d(mass[:, 1])
    a, b, n = _PD(), _PD(), c_long()
    if L.ms_isotope_peaks(px, py, len(x), float(threshold), int(max_n), byref(a), byref(b), byref(n)):
        return None
    k = int(n.value)
    return [(float(a[i]), float(b[i])) for i in range(k)]


def refine_peak_mass(mass, m0):
    L = _lib31("ms_refine_peak_mass")
    if L is None:
        return None
    mass = np.asarray(mass, dtype=np.float64)
    x, px = _d(mass[:, 0])
    y, py = _d(mass[:, 1])
    return float(L.ms_refine_peak_mass(px, py, len(x), float(m0)))


def refine_species_one(x, y, data_x, data_y, carrier, z0, z1, spacing, apex):
    """One species: (checked, apex, apex2 or None, apex2_rel, source, best) or None."""
    L = _lib31("ms_refine_species")
    if L is None:
        return None
    out, po = _out(6)
    nd = 0 if data_x is None else len(data_x)
    rc = L.ms_refine_species(x.ctypes.data_as(_PD), y.ctypes.data_as(_PD), len(x),
                             None if data_x is None else data_x.ctypes.data_as(_PD),
                             None if data_y is None else data_y.ctypes.data_as(_PD), nd, float(carrier), int(z0),
                             int(z1), float(spacing), float(apex), po)
    if rc:
        return None
    a2 = float(out[2])
    return (bool(out[0]), float(out[1]), None if not np.isfinite(a2) else a2, float(out[3]),
            "m/z data" if out[4] else "deconvolution", int(out[5]))


def group_isotopes(mass, threshold, spacing):
    """ms_group_isotopes2: list of (avg, apex, height, area, n_iso, first, last) or None."""
    L = _lib31("ms_group_isotopes2")
    if L is None:
        return None
    mass = np.asarray(mass, dtype=np.float64)
    x, px = _d(mass[:, 0])
    y, py = _d(mass[:, 1])
    p = [_PD() for _ in range(4)]
    ni, f0, f1, n = POINTER(c_int)(), _PD(), _PD(), c_long()
    if L.ms_group_isotopes2(px, py, len(x), float(threshold), float(spacing), byref(p[0]), byref(p[1]), byref(p[2]),
                            byref(p[3]), byref(ni), byref(f0), byref(f1), byref(n)):
        return None
    return [(float(p[0][k]), float(p[1][k]), float(p[2][k]), float(p[3][k]), int(ni[k]), float(f0[k]), float(f1[k]))
            for k in range(int(n.value))]


def restrict(spec, lo, hi):
    L = _lib31("ms_restrict")
    if L is None:
        return None
    x, px = _d(spec[:, 0])
    y, py = _d(spec[:, 1])
    a, b, n = _PD(), _PD(), c_long()
    if L.ms_restrict(px, py, len(x), float(lo), float(hi), byref(a), byref(b), byref(n)):
        return None
    k = int(n.value)
    out = np.empty((k, 2), dtype=np.float64)
    if k:
        out[:, 0] = np.ctypeslib.as_array(a, shape=(k,))
        out[:, 1] = np.ctypeslib.as_array(b, shape=(k,))
    return out


def peak_mask(spec, thr):
    L = _lib31("ms_peak_mask")
    if L is None:
        return None
    x, px = _d(spec[:, 0])
    y, py = _d(spec[:, 1])
    mask = np.zeros(len(x), dtype=np.uint8)
    npk = c_long()
    if L.ms_peak_mask(px, py, len(x), float(thr), mask.ctypes.data_as(POINTER(ctypes.c_ubyte)), byref(npk)):
        return None
    return mask.view(bool), int(npk.value)


def subtract_baseline(spec, width, factor):
    L = _lib31("ms_subtract_baseline")
    if L is None:
        return None
    x, px = _d(spec[:, 0])
    y, py = _d(spec[:, 1])
    out, po = _out(len(x))
    if L.ms_subtract_baseline(px, py, len(x), float(width or 0.0), float(factor), po):
        return None
    return out


def label_maxima(x, y, n, sep, rel):
    L = _lib31("ms_label_maxima")
    if L is None:
        return None
    a, pa = _d(x)
    b, pb = _d(y)
    idx = np.zeros(max(int(n), 1), dtype=np.int64 if ctypes.sizeof(c_long) == 8 else np.int32)
    k = c_long()
    if L.ms_label_maxima(pa, pb, len(a), int(n), float("nan") if sep is None else float(sep), float(rel),
                         idx.ctypes.data_as(POINTER(c_long)), byref(k)):
        return None
    return [int(v) for v in idx[:k.value]]


def xic_multi(data, event, mzs, tols):
    """Mass chromatograms of the windows mzs[k] +- tols[k] (m/z as recorded, no calibration in use)
    in one pass over the scans of an EngineFile (ms_xic_multi, in libraries that have it): the values
    of data.chromatogram(event, "xic", mz, tol) for each window, as rows; None: not available."""
    L = _lib31("ms_xic_multi")
    if L is None or not isinstance(data, EngineFile) or data.calibration_for(event) is not None or not data.xic_fast:
        return None
    c, w = [], []
    for m, t in zip(mzs, tols):
        lo, hi = float(m) - float(t), float(m) + float(t)   # as EngineFile.chromatogram
        c.append(0.5 * (lo + hi))
        w.append(0.5 * abs(hi - lo))
    n = len(data.event_scans(event))
    ca, pc = _d(c)
    wa, pw = _d(w)
    out, po = _out(len(c) * n)
    if L.ms_xic_multi(data.handle, int(event), pc, pw, len(c), po, n):
        return None
    return out.reshape(len(c), n)


# ------------------------------------------------------------------ lcms_compare_core (3.31)
_NAN = float("nan")
_BASELINE_KIND = {"none": 0, "offset": 1, "drift": 2, "rolling": 3}
_SCALE_KIND = {"abs": 0, "max": 1, "ref": 2}


def _opt(v):
    return _NAN if v is None else float(v)


def cmp_subtract_blank(t, y, bt, by):
    L = _lib31("ms_cmp_subtract_blank")
    if L is None:
        return None
    a, pa = _d(t)
    b, pb = _d(y)
    c, pc = _d(bt)
    d, pd = _d(by)
    if len(a) != len(b) or len(c) != len(d):
        return None
    out, po = _out(len(a))
    if L.ms_cmp_subtract_blank(pa, pb, len(a), pc, pd, len(c), po):
        return None
    return out


def cmp_find_apex(t, y, tc, win, smooth_points):
    """The apex time, None without a peak there; False when the library cannot (the caller computes it)."""
    L = _lib31("ms_cmp_find_apex")
    if L is None:
        return False
    a, pa = _d(t)
    b, pb = _d(y)
    if len(a) != len(b):
        return False
    r = c_double()
    if L.ms_cmp_find_apex(pa, pb, len(a), float(tc), float(win), int(smooth_points or 0), byref(r)):
        return False
    return None if r.value != r.value else float(r.value)


def cmp_baseline(t, y, kind, rolling_min):
    L = _lib31("ms_cmp_baseline")
    if L is None or kind not in _BASELINE_KIND:
        return None
    a, pa = _d(t)
    b, pb = _d(y)
    if len(a) != len(b):
        return None
    out, po = _out(len(a))
    if L.ms_cmp_baseline(pa, pb, len(a), _BASELINE_KIND[kind], float(rolling_min), po):
        return None
    return out


def cmp_process(t, y, s, shift):
    """process() of lcms_compare_core in the library: its dict, or None."""
    L = _lib31("ms_cmp_process")
    if L is None:
        return None
    kind, sc = s.get("baseline", "none"), s.get("scale", "abs")
    if kind not in _BASELINE_KIND or sc not in _SCALE_KIND:
        return None
    a, pa = _d(t)
    b, pb = _d(y)
    if len(a) != len(b):
        return None
    o = CmpOpts(int(s.get("smooth", 0) or 0), float(shift or 0.0), _opt(s.get("t0")), _opt(s.get("t1")),
                _BASELINE_KIND[kind], float(s.get("rolling_min", 1.0)), _SCALE_KIND[sc],
                _opt(s.get("ref_t")) if sc == "ref" else _NAN, float(s.get("align_win", 0.3)))
    to, pto = _out(len(a))
    yo, pyo = _out(len(a))
    n, f, top, miss = c_long(), c_double(), c_double(), c_int()
    if L.ms_cmp_process(pa, pb, len(a), byref(o), pto, pyo, byref(n), byref(f), byref(top), byref(miss)):
        return None
    k = int(n.value)
    return {"t": to[:k].copy(), "y": yo[:k].copy(), "factor": float(f.value), "top": float(top.value), "ok": k > 1,
            "ref_missing": bool(miss.value)}


def cmp_stack(tops, spans, ok, spacing, reverse, skew, stacked):
    L = _lib31("ms_cmp_stack")
    if L is None:
        return None
    a, pa = _d(tops)
    b, pb = _d(spans)
    k = np.ascontiguousarray([1 if v else 0 for v in ok], dtype=np.intc)
    n = len(a)
    if len(b) != n or len(k) != n:
        return None
    dx, pdx = _out(n)
    dy, pdy = _out(n)
    if L.ms_cmp_stack(pa, pb, k.ctypes.data_as(POINTER(c_int)), n, float(spacing), 1 if reverse else 0, float(skew),
                      1 if stacked else 0, pdx, pdy):
        return None
    return [(float(x), float(v)) for x, v in zip(dx, dy)]


def cmp_peaks(t, y, factor, thr_pct, guides):
    """(peak dicts sorted by time, index of the tallest or -1, its area % (NaN), area % per guide (NaN)), or
    None."""
    L = _lib31("ms_cmp_peaks")
    if L is None:
        return None
    a, pa = _d(t)
    b, pb = _d(y)
    if len(a) != len(b):
        return None
    g, pg = _d(list(guides) if len(guides) else [0.0])
    ng = len(guides)
    gp, pgp = _out(max(ng, 1))
    pk, n, main, mp = POINTER(LcPeak2)(), c_long(), c_long(), c_double()
    if L.ms_cmp_peaks(pa, pb, len(a), float(factor or 0.0), float(thr_pct), pg, ng, byref(pk), byref(n),
                      byref(main), byref(mp), pgp):
        return None
    peaks = [_peak_dict(pk[k]) for k in range(int(n.value))]
    return peaks, int(main.value), float(mp.value), [float(v) for v in gp[:ng]]


def cmp_lambda_max(A, wavelengths, times, interval_s, wl0, bw, smooth, shift, t0, t1):
    """(lambda max, time of the peak, background subtracted) or None (no such peak); False when the
    library cannot (the caller computes it)."""
    L = _lib31("ms_cmp_lambda_max")
    if L is None or np.ndim(A) != 2 or not np.size(A):
        return False
    A, pA = _f32(A)
    w, pw = _d(wavelengths)
    t, pt = _d(times)
    if len(w) != A.shape[1] or len(t) != A.shape[0]:
        return False
    lm, ta, bg = c_double(), c_double(), c_int()
    rc = L.ms_cmp_lambda_max(pA, A.shape[0], A.shape[1], pw, pt, float(interval_s), _opt(wl0), float(bw),
                             int(smooth or 0), float(shift or 0.0), _opt(t0), _opt(t1), byref(lm), byref(ta), byref(bg))
    if rc == -4:  # MS_NOT_FOUND
        return None
    if rc:
        return False
    return float(lm.value), float(ta.value), bool(bg.value)


def cmp_region(t, y, a, b):
    """_region_peak of lcms_compare_core in the library: its dict, or None (no data in the region); False when
    the library cannot (the caller computes it)."""
    L = _lib31("ms_cmp_region")
    if L is None:
        return False
    ta, pt = _d(t)
    ya, py = _d(y)
    if len(ta) != len(ya):
        return None
    cap = 2 * len(ta) + 2
    to, pto = _out(cap)
    yo, pyo = _out(cap)
    bo, pbo = _out(cap)
    res, pres = _out(7)
    n = c_long()
    rc = L.ms_cmp_region(pt, py, len(ta), float(a), float(b), pto, pyo, pbo, cap, byref(n), pres)
    if rc == -4:  # MS_NOT_FOUND
        return None
    if rc:
        return False
    k = int(n.value)
    return {"t0": float(res[0]), "t1": float(res[1]), "rt": float(res[2]), "height": float(res[3]),
            "area": float(res[4]), "b0": float(res[5]), "b1": float(res[6]),
            "t": to[:k].copy(), "y": yo[:k].copy(), "baseline": bo[:k].copy()}


# ------------------------------------------------------------------ ms_adducts (3.33)
class MassIon(ctypes.Structure):   # ms_mass_ion (3.33)
    _fields_ = [("polarity", c_int), ("adduct", c_int), ("mz", c_double), ("rel", c_double), ("mass", c_double)]


class MassCand(ctypes.Structure):  # ms_mass_cand (3.33)
    _fields_ = [("mass", c_double), ("score", c_double), ("both", c_int), ("assumed", c_int), ("first", c_long),
                ("n", c_long)]


class MassCand2(ctypes.Structure):  # ms_mass_cand2 (3.33)
    _fields_ = MassCand._fields_ + [("additive", c_int), ("m2", c_double), ("m2_ratio", c_double),
                                    ("m2_first", c_long), ("m2_n", c_long)]


_ADDUCT_NAMES = {}


def adduct_name(code):
    """Name of an adduct code of the library ("" for an unknown code or without the library)."""
    if code not in _ADDUCT_NAMES:
        L = _lib31("ms_adduct_name")
        _ADDUCT_NAMES[code] = (L.ms_adduct_name(int(code)) or b"").decode("ascii", "replace") if L is not None else ""
    return _ADDUCT_NAMES[code]


def additive_name(code):
    """Name of a mobile phase additive code of the library ("" for an unknown code or without the library)."""
    L = _lib31("ms_additive_name")
    return (L.ms_additive_name(int(code)) or b"").decode("ascii", "replace") if L is not None else ""


def neutral_masses_info(pmz, pit, nmz, nit, tol, min_rel, max_ions, top):
    """ms_neutral_masses2: the dict of ms_adducts.neutral_masses_info (candidates with their ions, base
    peaks), or None when the library lacks it or refuses the arguments (the caller computes it)."""
    L = _lib31("ms_neutral_masses2")
    if L is None or _lib31("ms_adduct_name") is None or _lib31("ms_additive_name") is None:
        return None
    a, pa = _d(pmz)
    b, pb = _d(pit)
    c, pc = _d(nmz)
    d, pd = _d(nit)
    if len(a) != len(b) or len(c) != len(d):
        return None
    ions, ni, cands, nc = POINTER(MassIon)(), c_long(), POINTER(MassCand2)(), c_long()
    base, pbase = _out(2)
    expl = (c_int * 2)()
    if L.ms_neutral_masses2(pa, pb, len(a), pc, pd, len(c), float(tol), float(min_rel), int(max_ions), int(top),
                            byref(ions), byref(ni), byref(cands), byref(nc), pbase, expl):
        return None

    def ion_list(first, n):
        out = []
        for j in range(int(first), int(first + n)):
            q = ions[j]
            out.append({"polarity": "+" if q.polarity > 0 else "-", "adduct": adduct_name(q.adduct),
                        "mz": float(q.mz), "rel": float(q.rel), "mass": float(q.mass)})
        return out
    masses = []
    for k in range(int(nc.value)):
        m = cands[k]
        has2 = m.m2_n > 0
        masses.append({"mass": float(m.mass), "score": float(m.score), "both": bool(m.both),
                       "assumed": bool(m.assumed), "additive": additive_name(m.additive) if m.additive >= 0 else None,
                       "m2": float(m.m2) if has2 else None, "m2_ratio": float(m.m2_ratio) if has2 else None,
                       "ions": ion_list(m.first, m.n), "m2_ions": ion_list(m.m2_first, m.m2_n)})
    nb = [None if v != v else float(v) for v in base]
    return {"masses": masses, "base": {"+": nb[0], "-": nb[1]},
            "explained": {"+": bool(expl[0]), "-": bool(expl[1])}}


def neutral_masses(pmz, pit, nmz, nit, tol, min_rel, max_ions, top):
    """The candidate list of neutral_masses_info, or None (the caller computes it)."""
    r = neutral_masses_info(pmz, pit, nmz, nit, tol, min_rel, max_ions, top)
    return None if r is None else r["masses"]


# ------------------------------------------------------------------ ms_shifts (3.35)
class ShiftPair(ctypes.Structure):  # ms_shift_pair (3.35)
    _fields_ = [("from_", c_long), ("to", c_long), ("delta", c_double), ("value", c_double), ("error", c_double),
                ("kind", c_int), ("sign", c_int), ("tag", c_int), ("k", c_int), ("s1", c_int), ("s2", c_int),
                ("iso", c_int), ("n_alt", c_int), ("alt_tag", c_int), ("alt_k", c_int), ("alt_s1", c_int),
                ("alt_s2", c_int), ("alt_sign", c_int), ("alt_iso", c_int), ("alt_error", c_double),
                ("ref", c_int), ("link", c_int)]


class ShiftSpecies(ctypes.Structure):  # ms_shift_species (3.35)
    _fields_ = [("status", c_int), ("parent", c_long), ("link", c_long), ("value", c_double), ("error", c_double),
                ("height", c_double), ("area", c_double)]


def mass_shifts(m, h, a, s, t, tol_da, tol_ppm, ref, kmax, max_extra, flags, min_rel=0.0):
    """ms_mass_shifts on checked arguments (ms_shifts.check_args): the dict of ms_shifts.find_shifts, or None
    when the library lacks the function or refuses the arguments (the caller computes it)."""
    L = _lib31("ms_mass_shifts")
    if L is None:
        return None
    pm, ph, pa = m.ctypes.data_as(_PD), h.ctypes.data_as(_PD), a.ctypes.data_as(_PD)
    ps, pt = s.ctypes.data_as(_PD), t.ctypes.data_as(_PD)
    n, ns, nt = len(m), len(s), len(t)
    pairs, npairs = POINTER(ShiftPair)(), c_long()
    species, comp, conj, ref_used = POINTER(ShiftSpecies)(), POINTER(c_int)(), _PD(), c_long()
    if L.ms_mass_shifts(pm, ph, pa, n, ps, ns, pt, nt, float(tol_da), float(tol_ppm), int(ref), int(kmax),
                        int(max_extra), int(flags), float(min_rel), byref(pairs), byref(npairs), byref(species),
                        byref(comp),
                        byref(conj), byref(ref_used)):
        return None
    out_pairs = []
    for i in range(int(npairs.value)):
        q = pairs[i]
        alt = None
        if q.n_alt > 0 and q.alt_sign != 0:
            alt = {"tag": q.alt_tag, "k": q.alt_k, "s1": q.alt_s1, "s2": q.alt_s2, "sign": q.alt_sign,
                   "iso": q.alt_iso, "error": float(q.alt_error)}
        out_pairs.append({"from": int(q.from_), "to": int(q.to), "delta": float(q.delta), "value": float(q.value),
                          "error": float(q.error), "kind": q.kind, "sign": q.sign, "tag": q.tag, "k": q.k,
                          "s1": q.s1, "s2": q.s2, "iso": q.iso, "n_alt": q.n_alt, "alt": alt, "ref": bool(q.ref),
                          "link": bool(q.link)})
    nc = nt + ns + 1
    out_species = []
    for i in range(n):
        q = species[i]
        row = [comp[i * nc + j] for j in range(nc)]
        out_species.append({"status": q.status, "parent": int(q.parent), "link": int(q.link), "value": float(q.value),
                            "error": float(q.error), "comp": (tuple(row[:nt]), tuple(row[nt:nt + ns]), row[-1]),
                            "height": float(q.height), "area": float(q.area)})
    nb = 2 * int(kmax) + 5
    out_conj = []
    for ti in range(nt):
        v = [float(conj[ti * nb + j]) for j in range(nb)]
        k1 = int(kmax) + 1
        out_conj.append({"height": v[:k1], "area": v[k1:2 * k1], "avg_height": v[2 * k1], "avg_area": v[2 * k1 + 1],
                         "n": int(v[-1]) if v[-1] == v[-1] else 0})
    return {"ref": int(ref_used.value), "pairs": out_pairs, "species": out_species, "conj": out_conj}
