"""
Data layer for LCMS Analysis (MSpektra): Shimadzu LabSolutions LC-MS files (.lcd).

MS scans are read with OpenSZRaw (Apache-2.0, clean-room reader, bundled).
On top of it this module
  * splits alternating scan events (e.g. positive/negative switching), which
    OpenSZRaw returns as one interleaved list, by the event number in the
    header of each scan (single quad; else by the pattern of the time steps),
  * reads each event's polarity and m/z range from the MS method stored in
    the file,
  * corrects the m/z scale: for LCMS-2020 files OpenSZRaw 0.2.0 reports m/z
    values twice too high (checked against a LabSolutions export); the scan
    range in the method is used to detect and undo this,
  * leaves out saturated detector readings, as LabSolutions does,
  * builds TIC, base peak and mass (extracted ion) chromatograms per event,
  * averages spectra over a time range, optionally minus a background range.

Checked against LabSolutions (KUA525-POZi_027, positive, 2 to 10 min):
same m/z values, intensities agree to r = 0.9995 when summed per peak.

PDA (UV/Vis) data are decoded by lcms_pda.py.
"""
import os
import re
import struct

import numpy as np

SATURATION = 2.1e9  # detector ceiling seen in LabSolutions files (int32 limit)
MAX_BINS = 20000000  # bins of one summed spectrum (a tiny bin width would need gigabytes)


def _method_events(path):
    """[(polarity '+'/'-'/None, mz_low, mz_high), ...] from the MS method."""
    try:
        import olefile
    except Exception:
        return []
    try:
        ole = olefile.OleFileIO(path)
    except Exception:
        return []
    try:
        streams = ["/".join(e) for e in ole.listdir()
                   if e[-1] == "GUC.1.METHOD" and len(e) >= 2 and "LCMS" in e[-2].upper()]
        if not streams:
            return []
        text = ole.openstream(streams[0]).read().decode("utf-16-le", "ignore")
    finally:
        ole.close()
    i = text.find('UPD ID="Mass"')
    if i < 0:
        return []
    m = re.search(r"([0-9A-Fa-f]{64,})", text[i:i + 20000])
    if not m:
        return []
    b = bytes.fromhex(m.group(1))
    events = []
    # an event block holds: end time (ms), event time (ms), 0, m/z start*10000, m/z end*10000
    off = 0
    while off < len(b) - 20:
        t_end, t_ev, zero, a, c = struct.unpack_from("<5I", b, off)
        off += 1
        if zero != 0 or not (0 < t_ev < 100000) or not (1000 <= t_end <= 100000000):
            continue
        if not (0 < a < c):
            continue
        lo, hi = a / 10000.0, c / 10000.0
        if not (5 <= lo < hi <= 10000):
            continue
        start = off - 1
        pol = None
        j = b.find(b"\xa0\xa5", start + 20, start + 260)
        if j > 0:
            pol = "-" if b[j - 1] & 0x04 else "+"
        events.append((pol, lo, hi))
        off = start + 20
    return events


def _scan_headers(path, n):
    """Event number (u16 at 0x08: 0, 1, ... in the order of the method) and
    polarity (bit 0x04 of the byte at 0x2C, set for negative ions) from the
    64 byte header of every single quadrupole scan ("Mass Raw Data"), for
    the scans OpenSZRaw decodes (it leaves out the others, so the same tests
    are made here); None for the other variants or when the counts differ."""
    try:
        import olefile
        ole = olefile.OleFileIO(path)
    except Exception:
        return None
    try:
        if ole.exists("TTFL Raw Data") or ole.exists("QTFL RawData/Centroid Index") or \
                not ole.exists("Mass Raw Data/MS Raw Data") or not ole.exists("Mass Raw Data/Spectrum Index"):
            return None
        si = ole.openstream("Mass Raw Data/Spectrum Index").read()
        raw = ole.openstream("Mass Raw Data/MS Raw Data").read()
    except Exception:
        return None
    finally:
        ole.close()
    off = np.frombuffer(si[:len(si) // 4 * 4], dtype="<u4").astype(np.int64)
    ends = np.r_[off[1:], len(raw)] if len(off) else off
    ev, pol = [], []
    for a, b in zip(off.tolist(), ends.tolist()):
        if a > b or b > len(raw) or b - a < 64:
            continue
        npk = raw[a + 0x36] | (raw[a + 0x37] << 8)
        pay = b - a - 64
        if npk == 0 or pay % npk or not (2 < pay // npk <= 6):
            continue
        ev.append(raw[a + 8] | (raw[a + 9] << 8))
        pol.append(-1 if raw[a + 0x2C] & 0x04 else 1)
    if len(ev) != n:
        return None
    return np.array(ev, dtype=np.int64), np.array(pol, dtype=np.int64)


def _event_period(rts, max_period=6):
    """Number of interleaved scan events, from the repeating pattern of
    retention time steps (e.g. 0.141 s, 0.159 s, 0.141 s ... = 2 events)."""
    d = np.diff(np.asarray(rts, dtype=float))
    if len(d) < 20:
        return 1
    d = d[:2000]
    for p in range(1, max_period + 1):
        a, b = d[:-p], d[p:]
        if np.all(np.abs(a - b) <= 0.002 + 0.02 * np.abs(a)):
            # p steps repeat; the pattern is only real if within-period steps differ
            if p == 1 or np.std(d[:p]) > 1e-4:
                return p
    return 1


class LCDFile(object):
    def __init__(self, path, progress=None):
        import openszraw
        self.path = path
        self.name = os.path.splitext(os.path.basename(path))[0]
        self.reader = openszraw.RawReader(path)
        n = self.reader.scan_count
        rts = np.zeros(n)
        tic = np.zeros(n)
        bpc = np.zeros(n)
        bpm = np.zeros(n)
        mzs, ints, lens = [], [], np.zeros(n, dtype=np.int64)
        for i in range(n):
            s = self.reader.read_spectrum(i)
            mz = np.asarray(s.mz, dtype=np.float64)
            it = np.asarray(s.intensity, dtype=np.float64)
            rts[i] = s.retention_time_sec
            tic[i] = s.total_ion_current if s.total_ion_current is not None else it.sum()
            if len(it):
                k = int(np.argmax(it))
                bpc[i], bpm[i] = it[k], mz[k]
            mzs.append(mz.astype(np.float32))
            ints.append(it.astype(np.float32))
            lens[i] = len(mz)
            if progress is not None and i % 500 == 0:
                progress(i, n)
        self.rt = rts / 60.0  # minutes
        self.offsets = np.concatenate([[0], np.cumsum(lens)])
        self.mz = np.concatenate(mzs) if mzs else np.zeros(0, np.float32)
        self.inten = np.concatenate(ints) if ints else np.zeros(0, np.float32)
        # scan events: the event number of the scan headers (single quad) when
        # they name 2 to 8 events, else the repeating pattern of the time steps
        # (scan k in event k % n). OpenSZRaw leaves out a scan it cannot decode:
        # with the pattern alone, one such scan made the whole file a single
        # event of mixed positive and negative scans
        hdr = _scan_headers(path, n)
        ids = np.unique(hdr[0]) if hdr is not None else np.zeros(0, np.int64)
        if hdr is not None and 2 <= len(ids) <= 8:
            self.n_events = len(ids)
            self.scan_event = np.searchsorted(ids, hdr[0])
        else:
            self.n_events = _event_period(rts)
            self.scan_event = np.arange(n) % self.n_events
        meth = _method_events(path)

        # m/z scale check against the scan range in the method
        self.mz_scale = 1.0
        if meth and len(self.mz):
            hi = max(m[2] for m in meth)
            k = float(self.mz.max()) / hi if hi else 1.0
            if round(k) >= 2 and abs(k - round(k)) < 0.05:
                self.mz_scale = 1.0 / round(k)
                self.mz *= self.mz_scale
                bpm *= self.mz_scale
        # saturated readings: left out (LabSolutions ignores them too)
        sat = self.inten >= SATURATION
        self.saturated = np.zeros(n, bool)
        if np.any(sat):
            scan_of = np.repeat(np.arange(n), lens)
            self.saturated[np.unique(scan_of[sat])] = True
            self.inten[sat] = 0
            for i in np.where(self.saturated)[0]:
                a, b = self.offsets[i], self.offsets[i + 1]
                seg = self.inten[a:b]
                tic[i] = float(seg.sum())
                if b > a:
                    k = int(np.argmax(seg))
                    bpc[i], bpm[i] = seg[k], self.mz[a + k]
        self.tic, self.bpc, self.bpm = tic, bpc, bpm
        self.events = []
        for e in range(self.n_events):
            pol, lo, hi = (meth[e] if e < len(meth) and len(meth) == self.n_events else (None, None, None))
            if pol is None and hdr is not None:
                # not in the method: the polarity of the event's scan headers (majority)
                p = hdr[1][self.scan_event == e]
                if len(p):
                    pol = "-" if int((p < 0).sum()) > int((p > 0).sum()) else "+"
            self.events.append({"polarity": pol, "mz_low": lo, "mz_high": hi})
        if self.n_events == 1 and self.events[0]["polarity"] is None and len(meth) == 1:
            self.events[0]["polarity"] = meth[0][0]

    # ------------------------------------------------------------------ info
    def event_scans(self, event):
        return np.flatnonzero(self.scan_event == event)

    def event_label(self, event):
        pol = self.events[event]["polarity"]
        name = {"+": "Positive", "-": "Negative"}.get(pol, "Scan event %d" % (event + 1))
        if self.n_events > 1 and pol is not None:
            name += " (event %d)" % (event + 1)
        return name

    def adduct_sign(self, event):
        return -1 if self.events[event]["polarity"] == "-" else 1

    # ------------------------------------------------------------ chromatograms
    def chromatogram(self, event=0, kind="tic", mz=None, tol=0.5):
        idx = self.event_scans(event)
        t = self.rt[idx]
        if kind == "tic":
            return t, self.tic[idx]
        if kind == "bpc":
            return t, self.bpc[idx]
        if kind == "xic":
            w = np.where(np.abs(self.mz - float(mz)) <= float(tol), self.inten, 0).astype(np.float64)
            cs = np.concatenate([[0.0], np.cumsum(w)])
            y = cs[self.offsets[1:]] - cs[self.offsets[:-1]]
            return t, y[idx]
        raise ValueError(kind)

    # ---------------------------------------------------------------- spectra
    def _scans_in(self, event, t0, t1):
        idx = self.event_scans(event)
        lo, hi = min(t0, t1), max(t0, t1)
        return idx[(self.rt[idx] >= lo) & (self.rt[idx] <= hi)]

    def _binned(self, scans, binw):
        if len(scans) == 0:
            return None, None
        parts_m = [self.mz[self.offsets[s]:self.offsets[s + 1]] for s in scans]
        parts_i = [self.inten[self.offsets[s]:self.offsets[s + 1]] for s in scans]
        m = np.concatenate(parts_m).astype(np.float64)
        it = np.concatenate(parts_i).astype(np.float64)
        if len(m) == 0:
            return None, None
        keys = np.round(m / binw).astype(np.int64)
        nb = int(keys.max()) - int(keys.min()) + 1
        if nb > MAX_BINS:
            raise ValueError("the bin width %g m/z is too small for m/z %g to %g (%d bins, at most %d)"
                             % (binw, keys.min() * binw, keys.max() * binw, nb, MAX_BINS))
        u, inv = np.unique(keys, return_inverse=True)
        s_i = np.bincount(inv, weights=it)
        s_m = np.bincount(inv, weights=m * it)
        with np.errstate(invalid="ignore", divide="ignore"):
            cm = np.where(s_i > 0, s_m / s_i, u * binw)
        return u, (cm, s_i / float(len(scans)))

    @staticmethod
    def _fill(keys, mz, it, binw):
        """Put the binned points on a complete m/z grid with zeros in empty
        bins, so a line plot drops to the baseline between peaks."""
        if len(keys) == 0:
            return np.zeros((0, 2))
        lo, hi = int(keys.min()), int(keys.max())
        grid = np.arange(lo - 1, hi + 2)
        out_mz = grid * float(binw)
        out_it = np.zeros(len(grid))
        pos = keys - (lo - 1)
        out_mz[pos] = mz
        out_it[pos] = it
        return np.column_stack([out_mz, out_it])

    def average(self, event, t0, t1, bg=None, binw=0.05, fill=True):
        """Averaged spectrum (N x 2: m/z, intensity) of one scan event over
        t0..t1 (minutes); bg = (b0, b1) or list of ranges to subtract.
        fill=True returns every bin of the m/z grid (zeros included)."""
        scans = self._scans_in(event, t0, t1)
        keys, res = self._binned(scans, binw)
        if res is None:
            return np.zeros((0, 2)), 0
        mz, it = res
        if bg:
            ranges = [bg] if np.ndim(bg[0]) == 0 else bg  # one (b0, b1) range or a list of them
            bscans = np.unique(np.concatenate([self._scans_in(event, a, b) for a, b in ranges]))
            bkeys, bres = self._binned(bscans, binw)
            if bres is not None:
                sub = np.zeros_like(it)
                pos = np.searchsorted(bkeys, keys)
                ok = (pos < len(bkeys)) & (bkeys[np.minimum(pos, len(bkeys) - 1)] == keys)
                sub[ok] = bres[1][pos[ok]]
                it = np.clip(it - sub, 0, None)
        if fill:
            return self._fill(keys, mz, it, binw), len(scans)
        keep = it > 0
        return np.column_stack([mz[keep], it[keep]]), len(scans)

    def scan_spectrum(self, event, t, binw=0.05, fill=True):
        """Spectrum of the single scan of this event nearest to time t:
        (N x 2, scan index)."""
        idx = self.event_scans(event)
        k = int(idx[np.argmin(np.abs(self.rt[idx] - t))])
        keys, res = self._binned(np.array([k]), binw)
        if res is None:
            return np.zeros((0, 2)), k
        mz, it = res
        if fill:
            return self._fill(keys, mz, it, binw), k
        keep = it > 0
        return np.column_stack([mz[keep], it[keep]]), k

    def summary(self):
        ev = ", ".join(self.event_label(e) for e in range(self.n_events))
        return "%d scans, %.1f min, %d scan event(s): %s" % (len(self.rt), self.rt[-1] if len(self.rt) else 0,
                                                            self.n_events, ev)


def write_spectrum_txt(path, data):
    np.savetxt(path, data, fmt="%.4f\t%.2f", delimiter="\t")


def write_spectrum_jdx(path, data, title, polarity="+", t0=None, t1=None):
    lines = ["##TITLE= %s" % title, "##JCAMP-DX= 4.24", "##DATA_TYPE= MASS SPECTRUM",
             "##ORIGIN= MSpektra, LCMS Analysis", "##IONIZATION_MODE= ESI%s" % (polarity or "+"),
             "##XUNITS= m/z", "##YUNITS= RELATIVE ABUNDANCE"]
    if t0 is not None:
        lines.append("##RETENTION_TIME= %.3f - %.3f" % (t0 * 60, t1 * 60))
        lines.append("##SCAN_TIME_UNITS= Seconds")
    lines += ["##NPOINTS= %d" % len(data), "##XYDATA= (XY..XY)"]
    lines += ["%.4f,%.2f" % (a, b) for a, b in data]
    lines.append("##END=")
    with open(path, "w", encoding="ascii", errors="replace") as fh:
        fh.write("\n".join(lines) + "\n")


# --------------------------------------------------------------------------
# sample information (LabSolutions "File Property" stream of the .lcd file)
# --------------------------------------------------------------------------
def _stox(v):
    """LabSolutions hex text: @StoX@<hex of the text>, @FtoX@<hex of a float>."""
    import struct
    if v is None:
        return ""
    v = v.strip()
    try:
        if v.startswith("@StoX@"):
            return bytes.fromhex(v[6:]).decode("cp1252", "replace").strip()
        if v.startswith("@FtoX@"):
            h = v[6:]
            if len(h) == 8:
                return "%g" % struct.unpack(">f", bytes.fromhex(h))[0]
            return h
    except (ValueError, struct.error):
        return ""
    return v


def read_sample_info(path):
    """Sample name, ID, type, operator, acquisition date, vial, injection
    volume, method and batch files, as LabSolutions shows them; {} if the
    file has none. Reads only a small part of the file."""
    import datetime
    try:
        import olefile
        f = olefile.OleFileIO(path)
    except Exception:
        return {}
    info = {}
    try:
        if not f.exists("File Property"):
            return {}
        txt = f.openstream("File Property").read().decode("latin-1", "replace")

        def tag(block, name):
            m = re.search(r"<%s>(.*?)</%s>" % (re.escape(name), re.escape(name)), block, re.S)
            return m.group(1) if m else None

        m = re.search(r"<SampleInfo>.*?</SampleInfo>", txt, re.S)
        if m:
            b = m.group(0)
            for key, name in (("sample_name", "smpl_name"), ("sample_id", "smpl_id"), ("sample_type", "smpl_type"),
                              ("operator", "operator_name"), ("vial", "szVialNum"), ("tray", "tray_name"),
                              ("inj_vol", "inj_vol")):
                v = _stox(tag(b, name))
                if v:
                    info[key] = v
            try:
                lo, hi = int(tag(b, "dwLowDateTime")), int(tag(b, "dwHighDateTime"))
                ft = (hi << 32) | (lo & 0xFFFFFFFF)
                if ft > 0:
                    info["acquired"] = datetime.datetime.fromtimestamp((ft - 116444736000000000) / 1e7)
            except (TypeError, ValueError, OSError, OverflowError):
                pass
        m = re.search(r"<SampleInfoFile>.*?</SampleInfoFile>", txt, re.S)
        if m:
            for key, name in (("method_file", "methodfile"), ("batch_file", "batchfile"), ("tune_file", "tunefile")):
                v = _stox(tag(m.group(0), name))
                if v:
                    info[key] = v.replace("/", "\\").split("\\")[-1]
        if f.exists("File Comment"):
            c = f.openstream("File Comment").read()
            try:
                c = c.decode("utf-16-le") if len(c) > 1 and c[1:2] == b"\x00" else c.decode("cp1252")
            except UnicodeDecodeError:
                c = ""
            c = c.strip("\x00 \r\n")
            if c:
                info["comment"] = c
    except Exception:
        pass
    finally:
        try:
            f.close()
        except Exception:
            pass
    return info
