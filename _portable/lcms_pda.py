"""
PDA (photodiode array, UV/Vis) data from Shimadzu LabSolutions .lcd files.

Clean-room decoder, written for LCMS Postrun (MS Analysis) from the structure of the file
itself (no LabSolutions code or third-party decoder was used).

Layout of the OLE stream "PDA 3D Raw Data/3D Raw Data", one record per
spectrum:
  16 byte header  b"RC\\0\\0", uint32 1, uint32 n_wavelengths, uint32 length
   8 zero bytes
  blocks of up to 256 wavelengths: uint16 n, n code bytes, uint16 n
Each value is a signed integer with a variable length code; the top three
bits of the first byte give the length (0: 1 byte/5 bit, 1: 2 bytes/13 bit,
2: 3 bytes/21 bit, 3: 4 bytes/29 bit). Inside a block the values are
differences along the wavelength axis; the first value of a block is
absolute. Units are micro-AU.

Checked against LabSolutions for KUA525-POZi_027.lcd: the 280 nm
chromatogram (correlation 0.9998) and the spectrum at 9.398 min
(differences within a few mAU) match the LabSolutions plots.
"""
import struct

import numpy as np


def _lib():
    """msengine_py when the C++ library is in use (it also decodes the PDA data when it reads
    the file: msengine_py.EngineFile.pda); the chromatograms and spectra below then run there
    and the Python code is the fallback (the reference of tests/check_p.py)."""
    try:
        import msengine_py
        return msengine_py if msengine_py.available() else None
    except Exception:
        return None


def _decode_py(raw, nrec, nwl):
    out = np.zeros((nrec, nwl), dtype=np.float32)
    off = 0
    n_raw = len(raw)
    for r in range(nrec):
        if off + 24 > n_raw or raw[off] != 0x52 or raw[off + 1] != 0x43:
            return out[:r]
        reclen = raw[off + 12] | (raw[off + 13] << 8) | (raw[off + 14] << 16) | (raw[off + 15] << 24)
        pos = off + 24
        end = off + reclen
        k = 0
        row = out[r]
        vals = []
        while pos + 2 <= end and k < nwl:
            n = raw[pos] | (raw[pos + 1] << 8)
            i = pos + 2
            stop = i + n
            acc = 0
            first = True
            while i < stop:
                b0 = raw[i]
                pre = b0 >> 5
                if pre == 0:
                    v = b0 & 0x1F
                    if v & 0x10:
                        v -= 0x20
                    i += 1
                elif pre == 1:
                    v = ((b0 & 0x1F) << 8) | raw[i + 1]
                    if v & 0x1000:
                        v -= 0x2000
                    i += 2
                elif pre == 2:
                    v = ((b0 & 0x1F) << 16) | (raw[i + 1] << 8) | raw[i + 2]
                    if v & 0x100000:
                        v -= 0x200000
                    i += 3
                else:
                    v = ((b0 & 0x1F) << 24) | (raw[i + 1] << 16) | (raw[i + 2] << 8) | raw[i + 3]
                    if v & 0x10000000:
                        v -= 0x20000000
                    i += 4
                acc = v if first else acc + v
                first = False
                vals.append(acc)
            k = len(vals)
            pos = stop + 2
        m = min(len(vals), nwl)
        row[:m] = vals[:m]
        off += reclen
    return out


def has_pda(path):
    try:
        import olefile
        with olefile.OleFileIO(path) as ole:
            return ole.exists("PDA 3D Raw Data/3D Raw Data") and ole.get_size("PDA 3D Raw Data/3D Raw Data") > 0
    except Exception:
        return False


class PDAData(object):
    """times (min), wavelengths (nm), absorbance matrix A[time, wavelength]
    in mAU."""

    def __init__(self, path):
        import olefile
        with olefile.OleFileIO(path) as ole:
            base = "PDA 3D Raw Data/"
            raw = ole.openstream(base + "3D Raw Data").read()
            wt = ole.openstream(base + "Wavelength Table").read()
            st = ole.openstream(base + "Status").read() if ole.exists(base + "Status") else b""
        nwl = struct.unpack_from("<i", wt, 0)[0]
        wl = np.array(struct.unpack_from("<%di" % nwl, wt, 4), dtype=float) / 100.0
        nrec, wl_lo, wl_hi, t0_ms, interval = 0, None, None, 0, None
        if len(st) >= 28:
            vals = struct.unpack_from("<7i", st, 0)
            nrec, wl_lo, wl_hi, t0_ms, interval = vals[1], vals[2] / 100.0, vals[3] / 100.0, vals[4], vals[6]
        if not nrec:
            # count records by walking the headers
            off = 0
            while off + 16 <= len(raw) and raw[off:off + 2] == b"RC":
                nrec += 1
                off += struct.unpack_from("<I", raw, off + 12)[0]
        if not interval:
            interval = 160
        A = _decode_py(raw, nrec, nwl) / 1000.0  # micro-AU -> mAU
        # drop wavelengths outside the acquired range (the first table entry
        # lies just below it and holds no data)
        keep = np.ones(nwl, bool)
        if wl_lo is not None and wl_hi is not None and wl_hi > wl_lo:
            keep = (wl >= wl_lo - 0.01) & (wl <= wl_hi + 0.61)
        self.wavelengths = wl[keep]
        self.A = np.ascontiguousarray(A[:, keep], dtype=np.float32)
        self.times = (t0_ms + np.arange(len(self.A)) * interval) / 60000.0
        self.interval_s = interval / 1000.0
        self.units = "mAU"

    # ----------------------------------------------------------------- access
    def _wl_index(self, wl, bw=0.0):
        wl = float(wl)
        half = max(float(bw), 0.0) / 2.0
        sel = np.where(np.abs(self.wavelengths - wl) <= half + 1e-9)[0]
        if len(sel) == 0:
            sel = np.array([int(np.argmin(np.abs(self.wavelengths - wl)))])
        return sel

    def chromatogram(self, wl, bw=4.0):
        """Absorbance at wl (nm), averaged over wl +- bw/2 (as 'wl, bw nm' in
        LabSolutions)."""
        E = _lib()
        if E is not None:
            y = E.pda_chromatogram(self.A, self.wavelengths, wl, bw)
            if y is not None:
                return self.times, y
        sel = self._wl_index(wl, bw)
        return self.times, self.A[:, sel].mean(axis=1)

    def max_plot(self):
        E = _lib()
        if E is not None:
            y = E.pda_max_plot(self.A)
            if y is not None:
                return self.times, y
        return self.times, self.A.max(axis=1)

    def spectrum(self, t0, t1=None, bg=None):
        """Spectrum at time t0 (min) or averaged over t0..t1; bg is a time or
        (b0, b1) range whose spectrum is subtracted."""
        E = _lib()
        if E is not None:
            s = E.pda_spectrum(self.A, self.times, self.interval_s, t0, t1, bg)
            if s is not None:
                return self.wavelengths, s
        s = self._mean_spectrum(t0, t1)
        if bg is not None:
            if isinstance(bg, (int, float)):
                s = s - self._mean_spectrum(bg, None)
            else:
                s = s - self._mean_spectrum(bg[0], bg[1])
        return self.wavelengths, s

    def _mean_spectrum(self, t0, t1=None):
        if t1 is None or abs(t1 - t0) < self.interval_s / 60.0:
            i = int(np.argmin(np.abs(self.times - t0)))
            return self.A[i].astype(float)
        lo, hi = min(t0, t1), max(t0, t1)
        m = (self.times >= lo) & (self.times <= hi)
        if not np.any(m):
            i = int(np.argmin(np.abs(self.times - (lo + hi) / 2)))
            return self.A[i].astype(float)
        return self.A[m].mean(axis=0).astype(float)

    def summary(self):
        return "%d spectra, %.0f to %.0f nm, every %.2f s" % (
            len(self.times), self.wavelengths[0], self.wavelengths[-1], self.interval_s)


def write_xy(path, x, y, header):
    with open(path, "w", encoding="ascii", errors="replace") as fh:
        fh.write(header + "\n")
        for a, b in zip(x, y):
            fh.write("%.4f\t%.4f\n" % (a, b))
