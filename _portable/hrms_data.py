"""
Bruker .d reader for HRMS Postrun (MS Analysis), for the maXis / micrOTOF /
compact / impact QTOF series and other instruments that write analysis.baf,
and for the newer TSF format (analysis.tsf + analysis.tsf_bin, written by
otofControl 6 and later on maXis II / impact II, without ion mobility),
and timsTOF data (analysis.tdf + analysis.tdf_bin): the MS1 frames are read
with the ion mobility dimension summed, MS/MS (PASEF, DIA) frames are left out.

Reading analysis.baf uses Bruker's Baf2Sql library (baf2sql_c.dll, bundled in
_portable\\baf2sql; see the licence notes there). It creates an SQLite index
(analysis.sqlite) inside the .d folder the first time a file is opened; when
the folder is write-protected, the .d folder is copied to a temporary folder
first. Reading analysis.tsf and analysis.tdf uses Bruker's TDF SDK library (timsdata.dll,
bundled in _portable\\timsdata with its licence files); it only reads.

The class has the same interface as lcms_data.LCDFile, so the mass
spectrometry tab of LCMS Postrun works with it unchanged:
  n_events, events, event_scans(e), event_label(e), adduct_sign(e), rt,
  saturated, chromatogram(e, kind, mz, tol), average(e, t0, t1, bg, binw),
  scan_spectrum(e, t, binw), summary(), name, path
plus an m/z calibration that can be replaced (internal calibration), one
per polarity: set_calibration(func or None, info, events), calibration_for(e).

Also reads mzML files (converted with ProteoWizard msconvert), as a fallback
when a .d folder cannot be read.
"""
import itertools
import os
import sys
import shutil
import sqlite3
import stat
import tempfile
import ctypes
from ctypes import c_uint64, c_uint32, c_char_p, c_int, c_double, POINTER

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_LIB = None


def _lib_path():
    names = ["baf2sql_c.dll"] if os.name == "nt" else ["libbaf2sql_c.so"]
    dirs = [os.path.join(HERE, "baf2sql"), HERE, "/tmp/pybaf2sql/Baf2Sql"]
    for d in dirs:
        for n in names:
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return p
    return None


def _lib():
    global _LIB
    if _LIB is not None:
        return _LIB
    path = _lib_path()
    if path is None:
        raise RuntimeError("Bruker's Baf2Sql library (baf2sql_c.dll) is missing from _portable\\baf2sql")
    if os.name == "nt":
        try:
            os.add_dll_directory(os.path.dirname(path))
        except Exception:
            pass
    lib = ctypes.cdll.LoadLibrary(path)
    lib.baf2sql_array_close_storage.argtypes = [c_uint64]
    lib.baf2sql_array_close_storage.restype = None
    lib.baf2sql_array_get_num_elements.argtypes = [c_uint64, c_uint64, POINTER(c_uint64)]
    lib.baf2sql_array_get_num_elements.restype = c_int
    lib.baf2sql_array_open_storage.argtypes = [c_int, c_char_p]
    lib.baf2sql_array_open_storage.restype = c_uint64
    lib.baf2sql_array_read_double.argtypes = [c_uint64, c_uint64, POINTER(c_double)]
    lib.baf2sql_array_read_double.restype = c_int
    lib.baf2sql_get_last_error_string.argtypes = [c_char_p, c_uint32]
    lib.baf2sql_get_last_error_string.restype = c_uint32
    lib.baf2sql_get_sqlite_cache_filename_v2.argtypes = [c_char_p, c_uint32, c_char_p, c_int]
    lib.baf2sql_get_sqlite_cache_filename_v2.restype = c_uint32
    _LIB = lib
    return lib


def _last_error(lib):
    n = lib.baf2sql_get_last_error_string(None, 0)
    buf = ctypes.create_string_buffer(max(int(n), 1))
    lib.baf2sql_get_last_error_string(buf, n)
    return buf.value.decode("utf-8", "replace") or "unknown Baf2Sql error"


_TIMS = None


def _tims_path():
    names = ["timsdata.dll"] if os.name == "nt" else ["libtimsdata.so"]
    for d in (os.path.join(HERE, "timsdata"), HERE):
        for n in names:
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return p
    return None


def _tims():
    """Bruker TDF SDK (timsdata.dll): tsf_ functions for analysis.tsf,
    tims_ functions for analysis.tdf."""
    global _TIMS
    if _TIMS is not None:
        return _TIMS
    path = _tims_path()
    if path is None:
        raise RuntimeError("Bruker's TDF SDK library (timsdata.dll) is missing from _portable\\timsdata")
    if os.name == "nt":
        try:
            os.add_dll_directory(os.path.dirname(path))
        except Exception:
            pass
    from ctypes import c_int64, c_int32, c_float
    lib = ctypes.cdll.LoadLibrary(path)
    lib.tsf_open.argtypes = [c_char_p, c_uint32]
    lib.tsf_open.restype = c_uint64
    lib.tsf_close.argtypes = [c_uint64]
    lib.tsf_close.restype = None
    lib.tsf_get_last_error_string.argtypes = [c_char_p, c_uint32]
    lib.tsf_get_last_error_string.restype = c_uint32
    lib.tsf_has_recalibrated_state.argtypes = [c_uint64]
    lib.tsf_has_recalibrated_state.restype = c_uint32
    lib.tsf_read_line_spectrum_v2.argtypes = [c_uint64, c_int64, POINTER(c_double), POINTER(c_float), c_int32]
    lib.tsf_read_line_spectrum_v2.restype = c_int32
    lib.tsf_read_profile_spectrum_v2.argtypes = [c_uint64, c_int64, POINTER(c_uint32), c_int32]
    lib.tsf_read_profile_spectrum_v2.restype = c_int32
    lib.tsf_index_to_mz.argtypes = [c_uint64, c_int64, POINTER(c_double), POINTER(c_double), c_uint32]
    lib.tsf_index_to_mz.restype = c_uint32
    lib.tims_open.argtypes = [c_char_p, c_uint32]
    lib.tims_open.restype = c_uint64
    lib.tims_close.argtypes = [c_uint64]
    lib.tims_close.restype = None
    lib.tims_get_last_error_string.argtypes = [c_char_p, c_uint32]
    lib.tims_get_last_error_string.restype = c_uint32
    lib.tims_has_recalibrated_state.argtypes = [c_uint64]
    lib.tims_has_recalibrated_state.restype = c_uint32
    lib.tims_index_to_mz.argtypes = [c_uint64, c_int64, POINTER(c_double), POINTER(c_double), c_uint32]
    lib.tims_index_to_mz.restype = c_uint32
    lib.tims_extract_profile_for_frame.argtypes = [c_uint64, c_int64, c_uint32, c_uint32, _PROFILE_CB, ctypes.c_void_p]
    lib.tims_extract_profile_for_frame.restype = c_uint32
    lib.tims_extract_chromatograms.argtypes = [c_uint64, _CHROM_JOBS, _CHROM_SINK, ctypes.c_void_p]
    lib.tims_extract_chromatograms.restype = c_uint32
    _TIMS = lib
    return lib


class _ChromJob(ctypes.Structure):
    _fields_ = [("id", ctypes.c_int64), ("time_begin", c_double), ("time_end", c_double),
                ("mz_min", c_double), ("mz_max", c_double), ("ook0_min", c_double), ("ook0_max", c_double)]


_PROFILE_CB = ctypes.CFUNCTYPE(None, ctypes.c_int64, c_uint32, POINTER(ctypes.c_int32), ctypes.c_void_p)
_CHROM_JOBS = ctypes.CFUNCTYPE(c_uint32, POINTER(_ChromJob), ctypes.c_void_p)
_CHROM_SINK = ctypes.CFUNCTYPE(c_uint32, ctypes.c_int64, c_uint32, POINTER(ctypes.c_int64),
                               POINTER(ctypes.c_uint64), ctypes.c_void_p)


def _tims_error(lib):
    n = lib.tims_get_last_error_string(None, 0)
    buf = ctypes.create_string_buffer(max(int(n), 1))
    lib.tims_get_last_error_string(buf, n)
    return buf.value.decode("utf-8", "replace") or "unknown TDF SDK error"


def _sqlite_ro(path):
    """Opens an SQLite file without writing to it (read only URI)."""
    try:
        import pathlib
        return sqlite3.connect(pathlib.Path(os.path.abspath(path)).as_uri() + "?mode=ro", uri=True)
    except (sqlite3.Error, ValueError):
        return sqlite3.connect(path)


def _recal_note(d):
    """Who recalibrated the file and when (calibration.sqlite, written by
    DataAnalysis), e.g. "DataAnalysis, 15.09.2026". Of several calibration
    states the latest is the one the Bruker libraries use."""
    try:
        con = _sqlite_ro(os.path.join(d, "calibration.sqlite"))
        con.text_factory = bytes
        try:
            info = {}
            for key, val in con.execute("SELECT KeyName, Value FROM CalibrationInfo WHERE KeyName IN "
                                        "('CalibrationDateTime', 'CalibrationSoftware') "
                                        "ORDER BY CalibrationState DESC"):
                info.setdefault(key.decode("utf-8", "replace"), val.decode("utf-8", "replace"))
        finally:
            con.close()
        date = info.get("CalibrationDateTime", "")[:10]
        if len(date) == 10 and date[4] == "-":
            date = "%s.%s.%s" % (date[8:10], date[5:7], date[:4])
        return ", ".join(x for x in (info.get("CalibrationSoftware", "DataAnalysis"), date) if x)
    except Exception:
        return "DataAnalysis"


def _baf_copy_dir(d):
    """Temporary copy of a write protected .d folder with analysis.baf, named
    after the folder and its analysis.baf (path, size, time), so that two files
    of the same name never share a copy (the C++ reader names it the same way)."""
    import hashlib
    name = os.path.splitext(os.path.basename(d.rstrip("\\/")))[0]
    try:
        st = os.stat(os.path.join(d, "analysis.baf"))
        tag = hashlib.md5(("%s|%d|%d" % (os.path.abspath(d), st.st_size, int(st.st_mtime))).encode("utf-8"))
        tag = tag.hexdigest()[:8]
    except OSError:
        tag = "copy"
    return os.path.join(tempfile.gettempdir(), "MSAnalysis_d", "%s_%s.d" % (name, tag))


_COPY_COUNT = itertools.count()


def _remove_tree(p):
    """Deletes a folder (best effort; copies of read only files are read only)."""
    def again(func, path, exc):
        try:
            os.chmod(path, stat.S_IWRITE)
            func(path)
        except OSError:
            pass
    if os.path.isdir(p):
        if sys.version_info >= (3, 12):
            shutil.rmtree(p, onexc=again)
        else:
            shutil.rmtree(p, onerror=again)


def _copy_tree_complete(src, dst):
    """Copies the folder src to dst through a temporary folder renamed to dst
    when the copy is complete, so that dst never exists half copied (an
    interrupted copy would else be taken for the copy forever)."""
    part = "%s.part%d_%d" % (dst, os.getpid(), next(_COPY_COUNT))
    _remove_tree(part)
    try:
        shutil.copytree(src, part)
        os.rename(part, dst)
    except BaseException:
        _remove_tree(part)
        if not os.path.isdir(dst):  # (else another reader made the copy meanwhile)
            raise


def file_properties(path, kind=None):
    """Instrument, sample, operator, method and acquisition date of a Bruker .d
    folder as {key: value}: GlobalMetadata of analysis.tsf / analysis.tdf, or
    the Properties of the Baf2Sql index of analysis.baf (written when the file
    is opened); {} for other files. kind: the reader's kind, when known."""
    if str(path).lower().endswith((".mzml", ".mzml.gz")):
        return {}
    d = find_d_folder(path)
    if d is None:
        return {}
    kind = kind or ""
    has_baf = os.path.isfile(os.path.join(d, "analysis.baf"))
    if "TSF" in kind or (not kind and not has_baf and is_tsf(d)):
        cands = [(os.path.join(d, "analysis.tsf"), "GlobalMetadata")]
    elif "timsTOF" in kind or (not kind and not has_baf and is_tdf(d)):
        cands = [(os.path.join(d, "analysis.tdf"), "GlobalMetadata")]
    else:
        cands = [(os.path.join(d, "analysis.sqlite"), "Properties"),
                 (os.path.join(_baf_copy_dir(d), "analysis.sqlite"), "Properties")]
    for p, table in cands:
        if not os.path.isfile(p):
            continue
        try:
            con = _sqlite_ro(p)
            try:
                return {str(k): v for k, v in con.execute("SELECT Key, Value FROM %s" % table).fetchall()
                        if k is not None}
            finally:
                con.close()
        except sqlite3.Error:
            continue
    return {}


def _tsf_error(lib):
    n = lib.tsf_get_last_error_string(None, 0)
    buf = ctypes.create_string_buffer(max(int(n), 1))
    lib.tsf_get_last_error_string(buf, n)
    return buf.value.decode("utf-8", "replace") or "unknown TDF SDK error"


def is_tsf(path):
    """Newer Bruker format without ion mobility (analysis.tsf + .tsf_bin)."""
    return os.path.isdir(path) and os.path.isfile(os.path.join(path, "analysis.tsf")) and \
        os.path.isfile(os.path.join(path, "analysis.tsf_bin"))


def is_tdf(path):
    """timsTOF format with ion mobility (analysis.tdf + .tdf_bin)."""
    return os.path.isdir(path) and os.path.isfile(os.path.join(path, "analysis.tdf")) and \
        os.path.isfile(os.path.join(path, "analysis.tdf_bin"))


def is_bruker_d(path):
    return os.path.isdir(path) and (os.path.isfile(os.path.join(path, "analysis.baf")) or is_tsf(path)
                                    or is_tdf(path))


def find_d_folder(path):
    """Accepts the .d folder, a file inside it (e.g. analysis.baf) or a
    folder containing exactly one .d folder."""
    path = os.path.abspath(path)
    if os.path.isfile(path) and (os.path.basename(path).lower().startswith("analysis.")
                                 or is_bruker_d(os.path.dirname(path))):
        path = os.path.dirname(path)
    if is_bruker_d(path):
        return path
    if os.path.isdir(path):
        subs = [os.path.join(path, n) for n in os.listdir(path) if n.lower().endswith(".d")]
        subs = [s for s in subs if is_bruker_d(s)]
        if len(subs) == 1:
            return subs[0]
    return None


class Cancelled(Exception):
    """The user cancelled a long averaging."""


def _tick(progress, i, n):
    """Progress callback of the long loops: progress(i, n) returning False
    cancels."""
    if progress is not None and progress(i, n) is False:
        raise Cancelled()


class _LogBins(object):
    """Streaming sum of profile spectra whose m/z axes differ (zero trimmed
    profiles): every spectrum is added into bins as wide as its point
    spacing (relative, as TOF spacing grows with m/z) as it comes, so the
    memory stays that of one spectrum however many scans are averaged
    (collecting all points first took gigabytes for 1000 scans)."""

    def __init__(self, r):
        self.r = float(r)
        self.k0 = None
        self.s_i = self.s_w = self.s_m = None

    def add(self, m, it):
        ok = m > 0
        if not ok.all():
            m, it = m[ok], it[ok]
        if not len(m):
            return
        it = np.asarray(it, float)
        keys = np.round(np.log(m) / self.r).astype(np.int64)
        lo, hi = int(keys.min()), int(keys.max())
        if self.k0 is None:
            self.k0 = lo
            n = hi - lo + 1
            self.s_i, self.s_w, self.s_m = np.zeros(n), np.zeros(n), np.zeros(n)
        elif lo < self.k0 or hi > self.k0 + len(self.s_i) - 1:
            new0 = min(lo, self.k0)
            n = max(hi, self.k0 + len(self.s_i) - 1) - new0 + 1
            off = self.k0 - new0
            for name in ("s_i", "s_w", "s_m"):
                a = np.zeros(n)
                old = getattr(self, name)
                a[off:off + len(old)] = old
                setattr(self, name, a)
            self.k0 = new0
        idx = keys - lo
        n = hi - lo + 1
        w = np.maximum(it, 1e-30)
        sl = slice(lo - self.k0, lo - self.k0 + n)
        self.s_i[sl] += np.bincount(idx, weights=it, minlength=n)
        self.s_w[sl] += np.bincount(idx, weights=w, minlength=n)
        self.s_m[sl] += np.bincount(idx, weights=m * w, minlength=n)

    def result(self):
        if self.s_w is None:
            return np.zeros((0, 2))
        used = self.s_w > 0
        return np.column_stack([self.s_m[used] / self.s_w[used], self.s_i[used]])


def _clean_axis(m, it):
    """Points with a finite m/z above 0, sorted by m/z (stable)."""
    m = np.asarray(m, float)
    it = np.asarray(it, float)
    ok = (m > 0) & np.isfinite(m)
    if not ok.all():
        m, it = m[ok], it[ok]
    if len(m) > 1 and np.any(m[1:] < m[:-1]):
        o = np.argsort(m, kind="stable")
        m, it = m[o], it[o]
    return m, it


def _interp_gap(x, xp, fp):
    """Linear interpolation of the profile (xp, fp) at x (sorted) that reads
    a gap of a zero trimmed profile as zero: 0 outside xp, and 0 between two
    points more than 1.5 x the spacing next to them apart (the points left
    out there were zeros). At a point of xp its own value."""
    n = len(xp)
    out = np.zeros(len(x))
    if n == 0 or not len(x):
        return out
    j = np.searchsorted(xp, x, "right") - 1
    ok = (j >= 0) & (x <= xp[-1])
    jc = np.clip(j, 0, n - 1)
    exact = ok & (xp[jc] == x)
    out[exact] = fp[jc[exact]]
    inner = ok & ~exact & (j < n - 1)
    if np.any(inner):
        jj = j[inner]
        h = xp[jj + 1] - xp[jj]
        d = np.diff(xp)
        nd = len(d)
        left = np.where(jj >= 1, d[np.clip(jj - 1, 0, nd - 1)], np.inf)
        right = np.where(jj + 1 < nd, d[np.clip(jj + 1, 0, nd - 1)], np.inf)
        ref = np.minimum(left, right)
        v = fp[jj] + (fp[jj + 1] - fp[jj]) * (x[inner] - xp[jj]) / h
        out[inner] = np.where(h > 1.5 * ref, 0.0, v)
    return out


def _resample_add(grid, acc, m, it):
    """Adds the profile (m, it) to the sum acc on the m/z axis grid (both
    cleaned and sorted) and returns the new (grid, acc). The spectrum is read
    on the grid by linear interpolation (_interp_gap), so that spectra whose
    m/z axes differ add up as the same spectrum would (heights and apex as
    in one scan when the scans are equal); its points farther than 0.75 x
    its own point spacing from every grid point (signal where the grid has
    none, as with zero trimmed profiles) are added to the grid first, with
    the sum read there from the grid."""
    m, it = _clean_axis(m, it)
    if not len(m):
        return grid, acc
    if len(m) == len(grid) and np.array_equal(m, grid):
        acc += it
        return grid, acc
    n, g = len(m), len(grid)
    if n > 1:
        d = np.diff(m)
        sp = np.minimum(np.concatenate([[np.inf], d]), np.concatenate([d, [np.inf]]))
    else:
        sp = np.zeros(1)
    j = np.searchsorted(grid, m, "left")
    dl = np.where(j > 0, m - grid[np.clip(j - 1, 0, max(g - 1, 0))], np.inf) if g else np.full(n, np.inf)
    dr = np.where(j < g, grid[np.clip(j, 0, max(g - 1, 0))] - m, np.inf) if g else np.full(n, np.inf)
    new = np.minimum(dl, dr) > 0.75 * sp
    if np.any(new):
        nm = m[new]
        grid, acc = np.insert(grid, j[new], nm), np.insert(acc, j[new], _interp_gap(nm, grid, acc))
    a, b = np.searchsorted(grid, m[0], "left"), np.searchsorted(grid, m[-1], "right")
    acc[a:b] += _interp_gap(grid[a:b], m, it)
    return grid, acc


MERGE_STICKS_PPM = 3.0


def _merge_sticks(mz, it, ppm=MERGE_STICKS_PPM):
    """Summed centroids (sorted): the sticks less than ppm apart from the one
    before are one ion, whose centroid scatters a little from scan to scan;
    each such group becomes one stick at the intensity weighted m/z with the
    summed intensity (the sums in order, as in the C++ core)."""
    mz = np.asarray(mz, float)
    it = np.asarray(it, float)
    n = len(mz)
    if n < 2:
        return mz, it
    brk = np.diff(mz) > mz[:-1] * (ppm * 1e-6)
    if brk.all():
        return mz, it
    st = np.concatenate([[0], np.nonzero(brk)[0] + 1])
    ln = np.diff(np.concatenate([st, [n]]))
    w = np.maximum(it, 1e-30)
    si, sw, sm = np.zeros(len(st)), np.zeros(len(st)), np.zeros(len(st))
    for k in range(int(ln.max())):
        g = ln > k
        idx = st[g] + k
        si[g] += it[idx]
        sw[g] += w[idx]
        sm[g] += mz[idx] * w[idx]
    out = sm / sw
    one = ln == 1
    out[one] = mz[st[one]]  # a single stick as it is
    return out, si


class _Spectra(object):
    """Concatenated centroid (line) spectra of many scans, for chromatograms."""

    def __init__(self, mzs, ints):
        lens = np.array([len(m) for m in mzs], dtype=np.int64)
        self.offsets = np.concatenate([[0], np.cumsum(lens)])
        self.mz = np.concatenate(mzs) if len(mzs) else np.zeros(0)
        self.inten = np.concatenate(ints) if len(ints) else np.zeros(0)

    def xic(self, mz, tol, mzfun=None):
        m = self.mz if mzfun is None else mzfun(self.mz)
        w = np.where(np.abs(m - float(mz)) <= float(tol), self.inten, 0.0)
        cs = np.concatenate([[0.0], np.cumsum(w)])
        return cs[self.offsets[1:]] - cs[self.offsets[:-1]]


def centroid(mz, it, rel=0.002):
    """Simple centroiding of a profile spectrum: local maxima above rel x the
    largest point, position from the three top points (Gaussian fit)."""
    mz = np.asarray(mz, float)
    it = np.asarray(it, float)
    if len(it) < 3 or it.max() <= 0:
        return np.zeros(0), np.zeros(0)
    thr = rel * it.max()
    k = np.where((it[1:-1] > it[:-2]) & (it[1:-1] >= it[2:]) & (it[1:-1] > thr))[0] + 1
    a, b, c = it[k - 1], it[k], it[k + 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        la, lb, lc = np.log(np.maximum(a, 1e-12)), np.log(b), np.log(np.maximum(c, 1e-12))
        den = la - 2 * lb + lc
        off = np.where((a > 0) & (c > 0) & (den < 0), 0.5 * (la - lc) / den, 0.0)
    off = np.clip(off, -0.5, 0.5)
    step = np.where(off >= 0, mz[np.minimum(k + 1, len(mz) - 1)] - mz[k], mz[k] - mz[k - 1])
    return mz[k] + off * step, b


class BrukerD(object):
    kind = "Bruker .d"

    def __init__(self, path, progress=None, raw_calibration=False):
        # an error or a cancel while opening closes the SDK handle at once (the
        # files would stay locked until the half made object is collected)
        try:
            self._open(path, progress, raw_calibration)
        except BaseException:
            self.close()
            raise

    def _open(self, path, progress, raw_calibration):
        d = find_d_folder(path)
        if d is None:
            raise ValueError("%s is not a Bruker .d folder with analysis.baf" % path)
        self.path = d
        self.name = os.path.splitext(os.path.basename(d.rstrip("\\/")))[0]
        self.lib = _lib()
        self.work = d
        try:
            cache = self._cache(d)
        except Exception:
            # write-protected folder: work on a temporary copy
            tmp = _baf_copy_dir(d)
            if not os.path.isdir(tmp):
                _copy_tree_complete(d, tmp)
            self.work = tmp
            cache = self._cache(tmp)
        self.handle = self.lib.baf2sql_array_open_storage(1 if raw_calibration else 0,
                                                          self.work.encode("utf-8"))
        if self.handle == 0:
            raise RuntimeError(_last_error(self.lib))
        self.recal_note = ""
        con = sqlite3.connect(cache)
        self.props = dict(con.execute("SELECT Key, Value FROM Properties").fetchall())
        rows = con.execute(
            "SELECT s.Id, s.Rt, s.SumIntensity, s.MaxIntensity, s.ProfileMzId, s.ProfileIntensityId, s.LineMzId, "
            "s.LineIntensityId, s.MzAcqRangeLower, s.MzAcqRangeUpper, ak.Polarity, ak.MsLevel, ak.ScanMode "
            "FROM Spectra s LEFT JOIN AcquisitionKeys ak ON s.AcquisitionKey = ak.Id ORDER BY s.Rt, s.Id").fetchall()
        con.close()
        if not rows:
            raise ValueError("No spectra in " + d)
        lvl = np.array([r[11] if r[11] is not None else 0 for r in rows])
        ms1 = [r for r, l in zip(rows, lvl) if l == lvl.min()]
        self.n_msms = len(rows) - len(ms1)
        self.rows = ms1
        n = len(ms1)
        self.rt = np.array([r[1] for r in ms1], float) / 60.0
        self.tic = np.array([r[2] or 0.0 for r in ms1], float)
        self.bpc = np.array([r[3] or 0.0 for r in ms1], float)
        self.saturated = np.zeros(n, bool)
        self.pol = np.array([1 if r[10] == 1 else 0 for r in ms1])
        self.has_profile = any(r[4] for r in ms1)
        self.has_line = any(r[6] for r in ms1)
        self.events = []
        self._scans = []
        for p in sorted(set(self.pol.tolist())):
            idx = np.where(self.pol == p)[0]
            lo = min(r[8] for r in (ms1[i] for i in idx) if r[8] is not None) if any(ms1[i][8] for i in idx) else None
            hi = max(r[9] for r in (ms1[i] for i in idx) if r[9] is not None) if any(ms1[i][9] for i in idx) else None
            self.events.append({"polarity": "-" if p == 1 else "+", "mz_low": lo, "mz_high": hi})
            self._scans.append(idx)
        self.n_events = len(self.events)
        self.calib = None  # callable m/z -> corrected m/z
        self.calib_info = ""
        self.calib_events = None
        self.xic_fast = True
        # centroid spectra of every scan, for mass chromatograms
        mzs, its = [], []
        for i, r in enumerate(ms1):
            if r[6] and r[7]:
                m, it = self._read(r[6]), self._read(r[7])
            elif r[4] and r[5]:
                m, it = centroid(self._read(r[4]), self._read(r[5]))
            else:
                m, it = np.zeros(0), np.zeros(0)
            mzs.append(m)
            its.append(it)
            if progress is not None and i % 200 == 0:
                progress(i, n)
        self.lines = _Spectra(mzs, its)
        if not np.any(self.tic):
            self.tic = np.array([float(x.sum()) for x in its])
            self.bpc = np.array([float(x.max()) if len(x) else 0.0 for x in its])
        if not raw_calibration:
            self.recal_note = self._baf_recal(ms1)

    def _baf_recal(self, rows):
        """A recalibration saved by DataAnalysis (calibration.sqlite or
        Calibrator.ami in the .d folder) is applied by Baf2Sql when it reads
        the spectra; it is recognised by comparing the m/z of the first scan
        with those of the calibration stored in analysis.baf."""
        r = next((r for r in rows if (r[6] and r[7]) or (r[4] and r[5])), None)
        if r is None:
            return ""
        ident = r[6] if (r[6] and r[7]) else r[4]
        raw = self.lib.baf2sql_array_open_storage(1, self.work.encode("utf-8"))
        if not raw:
            return ""
        try:
            n = c_uint64(0)
            if not self.lib.baf2sql_array_get_num_elements(raw, int(ident), ctypes.byref(n)) or not n.value:
                return ""
            a = np.empty(n.value, dtype=np.float64)
            if not self.lib.baf2sql_array_read_double(raw, int(ident), a.ctypes.data_as(POINTER(c_double))):
                return ""
        finally:
            self.lib.baf2sql_array_close_storage(raw)
        b = self._read(ident)
        if len(a) != len(b) or not np.any(np.abs(a - b) > 1e-9 * np.maximum(np.abs(b), 1.0)):
            return ""
        if os.path.isfile(os.path.join(self.path, "calibration.sqlite")):
            return _recal_note(self.path)
        return "DataAnalysis"

    # ---------------------------------------------------------------- low level
    def _cache(self, d):
        u8 = d.encode("utf-8")
        n = self.lib.baf2sql_get_sqlite_cache_filename_v2(None, 0, u8, 0)
        if n == 0:
            raise RuntimeError(_last_error(self.lib))
        buf = ctypes.create_string_buffer(int(n))
        self.lib.baf2sql_get_sqlite_cache_filename_v2(buf, n, u8, 0)
        return buf.value.decode("utf-8")

    def _read(self, ident):
        n = c_uint64(0)
        if not self.lib.baf2sql_array_get_num_elements(self.handle, int(ident), ctypes.byref(n)):
            raise RuntimeError(_last_error(self.lib))
        buf = np.empty(n.value, dtype=np.float64)
        if n.value and not self.lib.baf2sql_array_read_double(self.handle, int(ident),
                                                              buf.ctypes.data_as(POINTER(c_double))):
            raise RuntimeError(_last_error(self.lib))
        return buf

    def close(self):
        try:
            if getattr(self, "handle", 0):
                self.lib.baf2sql_array_close_storage(self.handle)
                self.handle = 0
        except Exception:
            pass

    def __del__(self):
        self.close()

    # ------------------------------------------------------------------ info
    def event_scans(self, event):
        return self._scans[event]

    def event_label(self, event):
        return {"+": "Positive", "-": "Negative"}.get(self.events[event]["polarity"], "Scans")

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

    def chromatogram(self, event=0, kind="tic", mz=None, tol=0.01):
        idx = self.event_scans(event)
        t = self.rt[idx]
        if kind == "tic":
            return t, self.tic[idx]
        if kind == "bpc":
            return t, self.bpc[idx]
        if kind == "xic":
            y = self.lines.xic(mz, tol, self.calibration_for(event))
            return t, y[idx]
        raise ValueError(kind)

    # ---------------------------------------------------------------- spectra
    def _scans_in(self, event, t0, t1):
        idx = self.event_scans(event)
        lo, hi = min(t0, t1), max(t0, t1)
        return idx[(self.rt[idx] >= lo) & (self.rt[idx] <= hi)]

    def profile(self, i):
        r = self.rows[i]
        if r[4] and r[5]:
            return self._read(r[4]), self._read(r[5]), True
        if r[6] and r[7]:
            return self._read(r[6]), self._read(r[7]), False
        return np.zeros(0), np.zeros(0), False

    def _sum(self, scans, progress=None):
        """Summed spectrum of the scans. Profile spectra on one common m/z
        axis are added point by point; otherwise (zero-trimmed profiles,
        m/z axes that differ) every spectrum is read on the m/z axis of the
        sum by linear interpolation, the axis taking in the points of the
        spectra where it has none (_resample_add). progress(i, n) is called
        every few scans; returning False cancels."""
        grid = acc = None
        clean = False
        n = len(scans)
        for k, i in enumerate(scans):
            if k % 10 == 0:
                _tick(progress, k, n)
            m, it, prof = self.profile(int(i))
            if not len(m):
                continue
            if not prof:
                return self._sum_lines(scans, progress)
            if grid is None:
                grid, acc = m, np.asarray(it, float).copy()
            elif len(m) == len(grid) and np.array_equal(m, grid):
                acc += it
            else:
                if not clean:
                    grid, acc = _clean_axis(grid, acc)
                    clean = True
                grid, acc = _resample_add(grid, acc, m, it)
        if grid is None:
            return np.zeros((0, 2))
        return np.column_stack([grid, acc])

    def _sum_lines(self, scans, progress=None):
        """Summed centroid spectrum: the centroids in 1 ppm bins; with several
        scans the sticks of one ion (its centroid differs a little from scan
        to scan) are then merged into one (_merge_sticks)."""
        bins = _LogBins(1e-6)  # 1 ppm bins
        n = len(scans)
        used = 0
        for k, i in enumerate(scans):
            if k % 50 == 0:
                _tick(progress, k, n)
            a, b = self.lines.offsets[i], self.lines.offsets[i + 1]
            if b > a:
                bins.add(self.lines.mz[a:b], self.lines.inten[a:b])
                used += 1
        s = bins.result()
        if used > 1 and len(s) > 1:
            m, it = _merge_sticks(s[:, 0], s[:, 1])
            s = np.column_stack([m, it])
        return s

    def average(self, event, t0, t1, bg=None, binw=None, fill=True, progress=None):
        """Averaged spectrum (N x 2) of the scans of one polarity between t0
        and t1 (min); bg = (b0, b1) or a list of ranges to subtract.
        progress(i, n), called during the summing, cancels when it returns
        False (hrms_data.Cancelled is raised)."""
        scans = self._scans_in(event, t0, t1)
        if not len(scans):
            return np.zeros((0, 2)), 0
        s = self._sum(scans, progress)
        if not len(s):
            return s, 0
        s[:, 1] /= float(len(scans))
        if bg:
            ranges = [bg] if isinstance(bg[0], (int, float)) else bg
            bscans = np.unique(np.concatenate([self._scans_in(event, a, b) for a, b in ranges]))
            if len(bscans):
                b = self._sum(bscans, progress)
                if len(b):
                    s[:, 1] = np.clip(s[:, 1] - np.interp(s[:, 0], b[:, 0], b[:, 1] / len(bscans), left=0, right=0),
                                      0, None)
        s[:, 0] = self._cal(s[:, 0], event)
        return s, len(scans)

    def scan_spectrum(self, event, t, binw=None, fill=True):
        idx = self.event_scans(event)
        k = int(idx[np.argmin(np.abs(self.rt[idx] - t))])
        s = self._sum([k])
        if len(s):
            s[:, 0] = self._cal(s[:, 0], event)
        return s, k

    def summary(self):
        ev = ", ".join(self.event_label(e) for e in range(self.n_events))
        inst = self.props.get("InstrumentName", "")
        txt = "%s%d MS scans, %.1f min, %s" % ((inst + ": ") if inst else "", len(self.rt),
                                             self.rt[-1] if len(self.rt) else 0, ev)
        if self.n_msms:
            txt += ", %d MS/MS scans (not shown)" % self.n_msms
        txt += ", profile" if self.has_profile else ", centroid only"
        if getattr(self, "recal_note", ""):
            txt += ", recalibrated (%s)" % self.recal_note
        return txt


class BrukerTSF(BrukerD):
    """Bruker .d folder in the TSF format (analysis.tsf: SQLite with the
    frames; analysis.tsf_bin: the spectra), read with Bruker's TDF SDK.
    Each frame is one spectrum (profile on the digitizer sample axis, plus
    line spectra); MS/MS frames are left out, as for analysis.baf."""
    kind = "Bruker .d (TSF)"

    def _open(self, path, progress, raw_calibration):
        d = find_d_folder(path)
        if d is None or not is_tsf(d):
            raise ValueError("%s is not a Bruker .d folder with analysis.tsf" % path)
        self.path = self.work = d
        self.name = os.path.splitext(os.path.basename(d.rstrip("\\/")))[0]
        self.lib = _tims()
        # the recalibration of DataAnalysis is used when there is one (as for analysis.baf)
        self.handle = self.lib.tsf_open(d.encode("utf-8"), 0 if raw_calibration else 1)
        if not self.handle:
            raise RuntimeError("The file could not be opened: " + _tsf_error(self.lib))
        # a recalibration saved by DataAnalysis (calibration.sqlite) is used, as in DataAnalysis
        self.recal_note = _recal_note(d) if self.lib.tsf_has_recalibrated_state(self.handle) else ""
        try:
            import pathlib
            uri = pathlib.Path(os.path.abspath(os.path.join(d, "analysis.tsf"))).as_uri() + "?mode=ro"
            con = sqlite3.connect(uri, uri=True)
        except (sqlite3.Error, ValueError):
            con = sqlite3.connect(os.path.join(d, "analysis.tsf"))
        try:
            self.props = dict(con.execute("SELECT Key, Value FROM GlobalMetadata").fetchall())
            rows = con.execute("SELECT Id, Time, Polarity, MsMsType, SummedIntensities, MaxIntensity FROM Frames "
                               "ORDER BY Time, Id").fetchall()
        finally:
            con.close()
        if not rows:
            raise ValueError("No spectra in " + d)
        lvl = np.array([int(r[3] or 0) for r in rows])
        ms1 = [r for r, l in zip(rows, lvl) if l == lvl.min()]
        self.n_msms = len(rows) - len(ms1)
        self.ids = np.array([int(r[0]) for r in ms1], dtype=np.int64)
        n = len(ms1)
        self.rt = np.array([r[1] for r in ms1], float) / 60.0
        self.tic = np.array([r[4] or 0.0 for r in ms1], float)
        self.bpc = np.array([r[5] or 0.0 for r in ms1], float)
        self.saturated = np.zeros(n, bool)
        self.pol = np.array([1 if r[2] == "-" else 0 for r in ms1])
        self.has_profile = str(self.props.get("HasProfileSpectra", "1")) == "1"
        self.has_line = str(self.props.get("HasLineSpectra", "0")) == "1"
        try:
            self.n_samples = int(float(self.props.get("DigitizerNumSamples", 0)))
        except ValueError:
            self.n_samples = 0
        try:
            lo, hi = float(self.props.get("MzAcqRangeLower")), float(self.props.get("MzAcqRangeUpper"))
        except (TypeError, ValueError):
            lo = hi = None
        self.events, self._scans = [], []
        for p in sorted(set(self.pol.tolist())):
            self.events.append({"polarity": "-" if p == 1 else "+", "mz_low": lo, "mz_high": hi})
            self._scans.append(np.where(self.pol == p)[0])
        self.n_events = len(self.events)
        self.calib = None
        self.calib_info = ""
        self.calib_events = None
        self.xic_fast = True
        self._prof_buf = max(self.n_samples, 1024)
        self._line_buf = 65536
        # line spectra of every frame (or centroided profiles), for mass chromatograms
        mzs, its = [], []
        for i in range(n):
            m, it = self._line(i) if self.has_line else (np.zeros(0), np.zeros(0))
            if not len(m) and self.has_profile:
                pm, pit = self._profile_frame(i)
                m, it = centroid(pm, pit)
            mzs.append(m)
            its.append(it)
            if progress is not None and i % 200 == 0:
                progress(i, n)
        self.lines = _Spectra(mzs, its)
        if not np.any(self.tic):
            self.tic = np.array([float(x.sum()) for x in its])
            self.bpc = np.array([float(x.max()) if len(x) else 0.0 for x in its])

    # ---------------------------------------------------------------- low level
    def _to_mz(self, i, index):
        index = np.ascontiguousarray(index, dtype=np.float64)
        out = np.empty(len(index), dtype=np.float64)
        if len(index) and not self.lib.tsf_index_to_mz(self.handle, int(self.ids[i]),
                                                       index.ctypes.data_as(POINTER(c_double)),
                                                       out.ctypes.data_as(POINTER(c_double)), len(index)):
            raise RuntimeError(_tsf_error(self.lib))
        return out

    def _profile_raw(self, i):
        """Profile intensities of frame i on the digitizer sample axis."""
        while True:
            buf = np.empty(self._prof_buf, dtype=np.uint32)
            n = self.lib.tsf_read_profile_spectrum_v2(self.handle, int(self.ids[i]),
                                                      buf.ctypes.data_as(POINTER(c_uint32)), self._prof_buf)
            if n < 0:
                raise RuntimeError(_tsf_error(self.lib))
            if n > self._prof_buf:
                self._prof_buf = int(n)
                continue
            return buf[:n]

    def _line(self, i):
        from ctypes import c_float
        while True:
            ib = np.empty(self._line_buf, dtype=np.float64)
            vb = np.empty(self._line_buf, dtype=np.float32)
            n = self.lib.tsf_read_line_spectrum_v2(self.handle, int(self.ids[i]), ib.ctypes.data_as(POINTER(c_double)),
                                                   vb.ctypes.data_as(POINTER(c_float)), self._line_buf)
            if n < 0:
                raise RuntimeError(_tsf_error(self.lib))
            if n > self._line_buf:
                self._line_buf = int(n)
                continue
            return self._to_mz(i, ib[:n]), vb[:n].astype(np.float64)

    def _profile_frame(self, i):
        it = self._profile_raw(i).astype(np.float64)
        return self._to_mz(i, np.arange(len(it), dtype=np.float64)), it

    def close(self):
        try:
            if getattr(self, "handle", 0):
                self.lib.tsf_close(self.handle)
                self.handle = 0
        except Exception:
            pass

    # ---------------------------------------------------------------- spectra
    def profile(self, i):
        if self.has_profile:
            m, it = self._profile_frame(i)
            return m, it, True
        a, b = self.lines.offsets[i], self.lines.offsets[i + 1]
        return self.lines.mz[a:b], self.lines.inten[a:b], False

    def _sum(self, scans, progress=None):
        """Summed profile spectrum: the frames are added on the digitizer
        sample axis (as the instrument records them), then converted to m/z
        with the calibration of the middle frame (the frames' calibrations
        differ by well below 1 ppm)."""
        scans = [int(i) for i in scans]
        if not scans:
            return np.zeros((0, 2))
        if not self.has_profile:
            return self._sum_lines(scans, progress)
        acc = None
        n = len(scans)
        for k, i in enumerate(scans):
            if k % 10 == 0:
                _tick(progress, k, n)
            it = self._profile_raw(i)
            if not len(it):
                continue
            if acc is None:
                acc = it.astype(np.float64)
            elif len(it) == len(acc):
                acc += it
            else:
                if len(it) > len(acc):
                    acc = np.concatenate([acc, np.zeros(len(it) - len(acc))])
                acc[:len(it)] += it
        if acc is None or not len(acc):
            return np.zeros((0, 2))
        mid = scans[len(scans) // 2]
        return np.column_stack([self._to_mz(mid, np.arange(len(acc), dtype=np.float64)), acc])

    def summary(self):
        ev = ", ".join(self.event_label(e) for e in range(self.n_events))
        inst = self.props.get("InstrumentName", "")
        txt = "%s%d MS spectra, %.1f min, %s" % ((inst + ": ") if inst else "", len(self.rt),
                                               self.rt[-1] if len(self.rt) else 0, ev)
        if self.n_msms:
            txt += ", %d MS/MS spectra (not shown)" % self.n_msms
        txt += ", TSF format, " + ("profile" if self.has_profile else "centroid only")
        if self.recal_note:
            txt += ", recalibrated (%s)" % self.recal_note
        return txt


class BrukerTDF(BrukerTSF):
    """timsTOF .d folder (analysis.tdf: SQLite with the frames;
    analysis.tdf_bin: the spectra), read with Bruker's TDF SDK. Each MS1 frame
    is one spectrum with the ion mobility scans summed (the SDK's profile of
    the frame on the digitizer sample axis); MS/MS frames (PASEF, DIA) are
    left out. Mass chromatograms come from the SDK's chromatogram extraction,
    so nothing has to be read in advance."""
    kind = "Bruker .d (timsTOF)"

    def _open(self, path, progress, raw_calibration):
        d = find_d_folder(path)
        if d is None or not is_tdf(d):
            raise ValueError("%s is not a timsTOF .d folder with analysis.tdf" % path)
        self.path = self.work = d
        self.name = os.path.splitext(os.path.basename(d.rstrip("\\/")))[0]
        self.lib = _tims()
        self.handle = self.lib.tims_open(d.encode("utf-8"), 0 if raw_calibration else 1)
        if not self.handle:
            raise RuntimeError("The file could not be opened: " + _tims_error(self.lib))
        self.recal_note = _recal_note(d) if self.lib.tims_has_recalibrated_state(self.handle) else ""
        con = _sqlite_ro(os.path.join(d, "analysis.tdf"))
        try:
            self.props = dict(con.execute("SELECT Key, Value FROM GlobalMetadata").fetchall())
            rows = con.execute("SELECT Id, Time, Polarity, MsMsType, SummedIntensities, MaxIntensity, NumScans "
                               "FROM Frames ORDER BY Time, Id").fetchall()
        finally:
            con.close()
        if not rows:
            raise ValueError("No spectra in " + d)
        lvl = np.array([int(r[3] or 0) for r in rows])
        ms1 = [r for r, l in zip(rows, lvl) if l == lvl.min()]
        self.n_msms = len(rows) - len(ms1)
        n = len(ms1)
        self.ids = np.array([int(r[0]) for r in ms1], dtype=np.int64)
        self.nscans = np.array([int(r[6] or 0) for r in ms1], dtype=np.int64)
        self.time_s = np.array([r[1] for r in ms1], float)
        self.rt = self.time_s / 60.0
        self.tic = np.array([r[4] or 0.0 for r in ms1], float)
        self.bpc = np.array([r[5] or 0.0 for r in ms1], float)
        self.saturated = np.zeros(n, bool)
        self.pol = np.array([1 if r[2] == "-" else 0 for r in ms1])
        self.has_profile = True
        self.has_line = False
        self.lines = None
        try:
            lo, hi = float(self.props.get("MzAcqRangeLower")), float(self.props.get("MzAcqRangeUpper"))
        except (TypeError, ValueError):
            lo = hi = None
        self.events, self._scans = [], []
        for p in sorted(set(self.pol.tolist())):
            self.events.append({"polarity": "-" if p == 1 else "+", "mz_low": lo, "mz_high": hi})
            self._scans.append(np.where(self.pol == p)[0])
        self.n_events = len(self.events)
        self.calib = None
        self.calib_info = ""
        self.calib_events = None
        self.xic_fast = False  # each mass chromatogram is extracted by the SDK from the whole file
        if progress is not None:
            progress(n, n)

    # ---------------------------------------------------------------- low level
    def _to_mz(self, i, index):
        index = np.ascontiguousarray(index, dtype=np.float64)
        out = np.empty(len(index), dtype=np.float64)
        if len(index) and not self.lib.tims_index_to_mz(self.handle, int(self.ids[i]),
                                                        index.ctypes.data_as(POINTER(c_double)),
                                                        out.ctypes.data_as(POINTER(c_double)), len(index)):
            raise RuntimeError(_tims_error(self.lib))
        return out

    def _profile_raw(self, i):
        """Intensities of frame i (all mobility scans) on the digitizer sample axis."""
        res = []

        def cb(ident, n, values, user):
            res.append(np.ctypeslib.as_array(values, shape=(int(n),)).copy() if n else np.zeros(0, np.int32))

        fn = _PROFILE_CB(cb)
        if not self.lib.tims_extract_profile_for_frame(self.handle, int(self.ids[i]), 0,
                                                       max(int(self.nscans[i]), 1), fn, None):
            raise RuntimeError(_tims_error(self.lib))
        return res[0] if res else np.zeros(0, np.int32)

    def close(self):
        try:
            if getattr(self, "handle", 0):
                self.lib.tims_close(self.handle)
                self.handle = 0
        except Exception:
            pass

    def _uncal(self, mz, cal=None):
        """Measured m/z that the replaced calibration maps onto mz."""
        cal = cal or self.calib
        r = float(mz)
        for _ in range(4):
            r -= float(np.atleast_1d(cal(np.array([r])))[0]) - float(mz)
        return r

    def _xic_raw(self, lo, hi):
        """Summed intensity between m/z lo and hi for every MS1 frame."""
        job = [_ChromJob(1, -1.0, float(self.time_s[-1]) + 1e6, float(lo), float(hi), 0.0, 1e3)]
        out = {}

        def gen(ptr, user):
            if not job:
                return 2
            ptr[0] = job.pop()
            return 1

        def sink(ident, n, frames, values, user):
            n = int(n)
            if n:
                out["f"] = np.ctypeslib.as_array(frames, shape=(n,)).copy()
                out["v"] = np.ctypeslib.as_array(values, shape=(n,)).astype(np.float64)
            return 1

        g, s = _CHROM_JOBS(gen), _CHROM_SINK(sink)
        if not self.lib.tims_extract_chromatograms(self.handle, g, s, None):
            raise RuntimeError(_tims_error(self.lib))
        y = np.zeros(len(self.ids))
        if "f" in out:
            order = np.argsort(self.ids)
            pos = np.searchsorted(self.ids[order], out["f"])
            ok = (pos < len(order))
            ok[ok] &= self.ids[order][pos[ok]] == out["f"][ok]
            np.add.at(y, order[pos[ok]], out["v"][ok])
        return y

    def chromatogram(self, event=0, kind="tic", mz=None, tol=0.01):
        if kind != "xic":
            return BrukerTSF.chromatogram(self, event, kind, mz, tol)
        idx = self.event_scans(event)
        lo, hi = float(mz) - float(tol), float(mz) + float(tol)
        cal = self.calibration_for(event)
        if cal is not None:
            lo, hi = self._uncal(lo, cal), self._uncal(hi, cal)
        return self.rt[idx], self._xic_raw(min(lo, hi), max(lo, hi))[idx]

    def summary(self):
        ev = ", ".join(self.event_label(e) for e in range(self.n_events))
        inst = self.props.get("InstrumentName", "")
        txt = "%s%d MS spectra, %.1f min, %s" % ((inst + ": ") if inst else "", len(self.rt),
                                               self.rt[-1] if len(self.rt) else 0, ev)
        if self.n_msms:
            txt += ", %d MS/MS frames (not shown)" % self.n_msms
        txt += ", timsTOF format (ion mobility summed), profile"
        if self.recal_note:
            txt += ", recalibrated (%s)" % self.recal_note
        return txt


# ============================================================================
# mzML (for instruments or files the Baf2Sql library cannot read)
# ============================================================================
# binary data array types of mzML (little endian, as the standard prescribes)
_MZML_TYPES = {"MS:1000521": np.float32, "MS:1000523": np.float64, "MS:1000519": np.int32, "MS:1000522": np.int64}
# compressions other than zlib: MS-Numpress (msconvert --numpressLinear, --numpressPic,
# --numpressSlof) and the truncation / prediction schemes of msconvert --mzTruncation
_MZML_OTHER_COMPRESSION = {
    "MS:1002312": "MS-Numpress linear prediction compression", "MS:1002313": "MS-Numpress positive integer compression",
    "MS:1002314": "MS-Numpress short logged float compression",
    "MS:1002746": "MS-Numpress linear prediction compression followed by zlib compression",
    "MS:1002747": "MS-Numpress positive integer compression followed by zlib compression",
    "MS:1002748": "MS-Numpress short logged float compression followed by zlib compression",
    "MS:1003088": "truncation, delta prediction and zlib compression",
    "MS:1003089": "truncation, linear prediction and zlib compression",
    "MS:1003090": "delta prediction and zlib compression", "MS:1003091": "linear prediction and zlib compression"}


def _mzml_compression_error(name):
    return ("This mzML file stores its spectra with %s, which MS Analysis cannot read. Convert the data again with "
            "ProteoWizard msconvert without that option (zlib compression is fine)." % name)


def _read_mzml(path, progress=None, want_offsets=False):
    """Minimal mzML reader (no extra packages): yields one dict per spectrum
    with ms level, time (min), polarity, profile flag, scan window, m/z and
    intensity arrays. Supports 32/64 bit floats and integers, zlib or no
    compression (MS-Numpress and other compressions are reported as an error).
    path may also be an open binary file (e.g. positioned at one <spectrum>
    of an indexed mzML, wrapped by _one_spectrum)."""
    import base64
    import zlib
    import xml.etree.ElementTree as ET
    cur = None
    arrays = None
    k = 0
    for ev, el in ET.iterparse(path, events=("start", "end")):
        tag = el.tag.rsplit("}", 1)[-1]
        if ev == "start":
            if tag == "spectrum":
                cur = {"level": 1, "rt": 0.0, "pol": 0, "prof": False, "lo": None, "hi": None, "mz": None, "it": None}
            elif tag == "binaryDataArray" and cur is not None:
                arrays = {"kind": None, "dtype": np.float64, "zlib": False, "other": None}
            continue
        if cur is None:
            el.clear()
            continue
        if tag == "cvParam":
            acc = el.get("accession")
            val = el.get("value")
            if arrays is not None:
                if acc == "MS:1000514":
                    arrays["kind"] = "mz"
                elif acc == "MS:1000515":
                    arrays["kind"] = "it"
                elif acc in _MZML_TYPES:
                    arrays["dtype"] = _MZML_TYPES[acc]
                elif acc == "MS:1000574":
                    arrays["zlib"] = True
                elif acc in _MZML_OTHER_COMPRESSION:
                    arrays["other"] = el.get("name") or _MZML_OTHER_COMPRESSION[acc]
            elif acc == "MS:1000511":
                cur["level"] = int(float(val))
            elif acc == "MS:1000129":
                cur["pol"] = 1
            elif acc == "MS:1000128":
                cur["prof"] = True
            elif acc == "MS:1000016":
                t = float(val)
                # unit by name or accession (UO:0000010 second, UO:0000031 minute, UO:0000028 millisecond)
                unit = (el.get("unitName") or "").lower()
                ua = el.get("unitAccession") or ""
                if unit.startswith("second") or ua == "UO:0000010":
                    t /= 60.0
                elif unit.startswith("millisecond") or ua == "UO:0000028":
                    t /= 60000.0
                cur["rt"] = t
            elif acc == "MS:1000501":
                cur["lo"] = float(val)
            elif acc == "MS:1000500":
                cur["hi"] = float(val)
            elif acc == "MS:1000285":
                cur["tic"] = float(val)
        elif tag == "binary" and arrays is not None:
            if arrays["other"] and arrays["kind"]:
                raise ValueError(_mzml_compression_error(arrays["other"]))
            raw = base64.b64decode(el.text or "")
            if arrays["zlib"] and raw:
                raw = zlib.decompress(raw)
            itemsize = np.dtype(arrays["dtype"]).itemsize
            if len(raw) % itemsize:
                raw = raw[:len(raw) - len(raw) % itemsize]
            arrays["data"] = np.frombuffer(raw, dtype=arrays["dtype"]).astype(np.float64)
        elif tag == "binaryDataArray" and arrays is not None:
            if arrays["kind"] and "data" in arrays:
                cur[arrays["kind"]] = arrays["data"]
            arrays = None
        elif tag == "spectrum":
            if cur["mz"] is None:
                cur["mz"] = np.zeros(0)
            if cur["it"] is None:
                cur["it"] = np.zeros(len(cur["mz"]))
            yield cur
            k += 1
            if progress is not None and k % 200 == 0:
                progress(k, 0)
            cur = None
            el.clear()
        elif tag in ("chromatogram",):
            el.clear()


def _mzml_offsets(path):
    """Byte offsets of the spectra of an indexed mzML (the index at the end
    of the file), or None."""
    import re
    try:
        size = os.path.getsize(path)
        with open(path, "rb") as fh:
            fh.seek(max(0, size - 4096))
            tail = fh.read()
            m = re.search(rb"<indexListOffset>\s*(\d+)\s*</indexListOffset>", tail)
            if not m:
                return None
            fh.seek(int(m.group(1)))
            idx = fh.read(size - int(m.group(1)))
        start = idx.find(b'<index name="spectrum">')
        if start < 0:
            return None
        end = idx.find(b"</index>", start)
        out = [int(x) for x in re.findall(rb"<offset[^>]*>\s*(\d+)\s*</offset>", idx[start:end if end > 0 else None])]
        return out or None
    except (OSError, ValueError):
        return None


class _OneSpectrum(object):
    """File object that serves one <spectrum> element at the given offset
    (wrapped so the XML parser sees a complete document)."""

    def __init__(self, path, offset):
        self.fh = open(path, "rb")
        self.fh.seek(offset)
        self.parts = [b"<mzML>"]
        self.done = False
        self.buf = b""

    def read(self, n=-1):
        if self.parts:
            return self.parts.pop(0)
        if self.done:
            return b""
        chunk = self.fh.read(65536)
        if not chunk:
            self.done = True
            self.fh.close()
            return b"</mzML>"
        self.buf += chunk
        k = self.buf.find(b"</spectrum>")
        if k >= 0:
            out = self.buf[:k + len(b"</spectrum>")] + b"</mzML>"
            self.done = True
            self.fh.close()
            return out
        out, self.buf = self.buf[:-12], self.buf[-12:]  # keep a possible partial end tag
        return out


class MzMLFile(BrukerD):
    kind = "mzML"
    CACHE_BYTES = 400 * 1024 * 1024  # profile scans kept in memory

    def __init__(self, path, progress=None):
        self.path = path
        self.name = os.path.splitext(os.path.basename(path))[0]
        self.props = {}
        offsets = _mzml_offsets(path)
        rows, specs, mzs, its, tics = [], [], [], [], []
        self._offsets = []
        kept = 0
        n_msms = 0
        shared = None  # the same m/z axis is stored once when scans share it
        for k, sp in enumerate(_read_mzml(path, progress)):
            if sp["level"] != 1:
                n_msms += 1
                continue
            rows.append((sp["rt"], sp["pol"], sp["prof"], sp["lo"], sp["hi"]))
            self._offsets.append(offsets[k] if offsets is not None and k < len(offsets) else None)
            m, it = sp["mz"], sp["it"]
            if len(m) != len(it):  # a damaged spectrum: the points both arrays have (as the C++ reader)
                n = min(len(m), len(it))
                m, it = m[:n], it[:n]
            # TIC is measured before centroiding; apex heights discard profile signal.
            tic = sp.get("tic")
            if tic is None or not np.isfinite(tic) or tic < 0:
                tic = float(np.sum(it, dtype=np.float64))
            tics.append(float(tic))
            if sp["prof"]:
                cm, ci = centroid(m, it)
            else:
                cm, ci = m, it
            mzs.append(cm)
            its.append(ci)
            if shared is not None and len(shared) == len(m) and np.array_equal(shared, m):
                m = shared
            else:
                shared = m
                kept += m.nbytes
            it32 = np.asarray(it, np.float32)
            kept += it32.nbytes
            # profiles stay in memory up to CACHE_BYTES when they can be read
            # again from the file by offset; otherwise all stay (float32)
            if self._offsets[-1] is not None and kept > self.CACHE_BYTES:
                specs.append(None)
            else:
                specs.append((m, it32))
        if not rows:
            raise ValueError("No MS1 spectra in " + path)
        order = np.argsort([r[0] for r in rows], kind="stable")
        rows = [rows[i] for i in order]
        specs = [specs[i] for i in order]
        mzs = [mzs[i] for i in order]
        its = [its[i] for i in order]
        self._offsets = [self._offsets[i] for i in order]
        self._specs = specs
        self._prof = [r[2] for r in rows]
        self._tic = [tics[i] for i in order]
        self.n_msms = n_msms
        self.rt = np.array([r[0] for r in rows])
        self.pol = np.array([r[1] for r in rows])
        self.has_profile = any(self._prof)
        self.has_line = not all(self._prof)
        self.saturated = np.zeros(len(rows), bool)
        self.lines = _Spectra(mzs, its)
        self.tic = np.array(self._tic)
        self.bpc = np.array([float(it.max()) if len(it) else 0.0 for it in its])
        self.events, self._scans = [], []
        for p in sorted(set(self.pol.tolist())):
            idx = np.where(self.pol == p)[0]
            los = [rows[i][3] for i in idx if rows[i][3] is not None]
            his = [rows[i][4] for i in idx if rows[i][4] is not None]
            self.events.append({"polarity": "-" if p == 1 else "+", "mz_low": min(los) if los else None,
                                "mz_high": max(his) if his else None})
            self._scans.append(idx)
        self.n_events = len(self.events)
        self.calib = None
        self.calib_info = ""
        self.calib_events = None
        self.xic_fast = True
        self.handle = 0

    def profile(self, i):
        sp = self._specs[i]
        if sp is None:  # not kept in memory: read again from the file
            try:
                for one in _read_mzml(_OneSpectrum(self.path, self._offsets[i])):
                    m, it = one["mz"], one["it"]
                    n = min(len(m), len(it))
                    return m[:n], it[:n], self._prof[i]
                raise ValueError("no spectrum at the position given by the index of the file")
            except Exception as ex:
                # no silent switch to the centroids: a spectrum summed partly from
                # profiles and partly from centroids would be wrong
                raise RuntimeError("Scan %d of %s could not be read again from the file (moved, changed or deleted "
                                   "while open?): %s" % (i + 1, os.path.basename(self.path), ex))
        m, it = sp
        return m, np.asarray(it, np.float64), self._prof[i]

    def close(self):
        pass


def open_hrms(path, progress=None):
    # the C++ core (msengine) reads every format when it is available; the
    # Python readers below are the fallback
    try:
        import msengine_py
        if msengine_py.available():
            target = path
            if not path.lower().endswith((".mzml", ".mzml.gz")):
                target = find_d_folder(path) or path
            return msengine_py.EngineFile(target, progress)
    except Exception as ex:
        print("[engine] Python reader used for %s: %s" % (os.path.basename(path), ex))
    if path.lower().endswith((".mzml", ".mzml.gz")):
        return MzMLFile(path, progress)
    d = find_d_folder(path)
    if d is not None and not os.path.isfile(os.path.join(d, "analysis.baf")):
        if is_tsf(d):
            return BrukerTSF(d, progress)
        if is_tdf(d):
            return BrukerTDF(d, progress)
    return BrukerD(path, progress)


def write_spectrum_txt(path, data):
    """m/z (6 decimals: 0.001 ppm at m/z 1000) and intensity (10 significant
    digits: averaged intensities below 0.01 are kept), tab separated."""
    np.savetxt(path, np.asarray(data, float).reshape(-1, 2), fmt="%.6f\t%.10g", delimiter="\t")


def write_spectrum_jdx(path, data, title, polarity="+", t0=None, t1=None):
    """JCAMP-DX mass spectrum (XY pairs) with the m/z to 6 decimals, the
    polarity and the time range (the export of HRMS Postrun)."""
    data = np.asarray(data, float).reshape(-1, 2)
    lines = ["##TITLE= %s" % title, "##JCAMP-DX= 4.24", "##DATA TYPE= MASS SPECTRUM",
             "##ORIGIN= MS Analysis, HRMS Postrun", "##OWNER= ", "##IONIZATION MODE= ESI%s" % (polarity or "+"),
             "##XUNITS= M/Z", "##YUNITS= RELATIVE ABUNDANCE"]
    if t0 is not None:
        lines.append("##RETENTION TIME= %.3f - %.3f" % (t0 * 60, (t1 if t1 is not None else t0) * 60))
        lines.append("##SCAN TIME UNITS= SECONDS")
    lines += ["##NPOINTS= %d" % len(data), "##XYDATA= (XY..XY)"]
    lines += ["%.6f,%.10g" % (a, b) for a, b in data]
    lines.append("##END=")
    with open(path, "w", encoding="ascii", errors="replace") as fh:
        fh.write("\n".join(lines) + "\n")
