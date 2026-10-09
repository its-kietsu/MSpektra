"""
JCAMP-DX (.jdx / .dx / .jcamp) support for UniDec (portable add-on).

* Lets every part of UniDec that opens data files (UniDec, MetaUniDec,
  UniChrom, batch processing, IsoDec) read JCAMP-DX mass spectra directly,
  e.g. the averaged spectra exported by Shimadzu LabSolutions.
  On opening, UniDec itself saves the imported spectrum as
  <name>_unidecfiles\\<name>_rawdata.txt.
* Can also be run as a converter:
      python unidec_jcamp.py file1.jdx file2.jdx folder ...
  writes <name>_unidec.txt (tab separated m/z and intensity) next to each file.

Supported JCAMP-DX content
  ##XYDATA= / ##PEAK TABLE= / ##XYPOINTS= / ##DATA TABLE=
  forms (XY..XY), (XYW..XYW), (XY), (X++(Y..Y)) with AFFN or compressed
  ASDF (SQZ/DIF/DUP) values, XFACTOR/YFACTOR, several blocks per file and
  NTUPLES pages (each page/block becomes one scan with its retention time).

UniDec itself is unchanged; this module registers an extra importer at start.
"""
import os
import re
import sys

import numpy as np

JCAMP_EXTENSIONS = (".jdx", ".dx", ".jcamp", ".jcm")

_DATA_LABELS = ("XYDATA", "PEAKTABLE", "XYPOINTS", "DATATABLE")

_SQZ = {"@": 0, "A": 1, "B": 2, "C": 3, "D": 4, "E": 5, "F": 6, "G": 7, "H": 8, "I": 9,
        "a": -1, "b": -2, "c": -3, "d": -4, "e": -5, "f": -6, "g": -7, "h": -8, "i": -9}
_DIF = {"%": 0, "J": 1, "K": 2, "L": 3, "M": 4, "N": 5, "O": 6, "P": 7, "Q": 8, "R": 9,
        "j": -1, "k": -2, "l": -3, "m": -4, "n": -5, "o": -6, "p": -7, "q": -8, "r": -9}
_DUP = {"S": 1, "T": 2, "U": 3, "V": 4, "W": 5, "X": 6, "Y": 7, "Z": 8, "s": 9}

_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
_AFFN_LINE = re.compile(r"^[\s,;]*(?:%s[\s,;]*)+$" % _NUM)
_NUM_RE = re.compile(_NUM)


# --------------------------------------------------------------------------
# Low level parsing
# --------------------------------------------------------------------------
def _norm_label(label):
    """JCAMP labels ignore case, spaces, '-', '/', '_'."""
    return re.sub(r"[\s\-/_]", "", label).upper()


def _to_float(text, default=None):
    try:
        m = _NUM_RE.search(text)
        return float(m.group(0)) if m else default
    except (TypeError, ValueError):
        return default


def _all_floats(text):
    return [float(x) for x in _NUM_RE.findall(text)]


def _decode_asdf_line(line):
    """Decode one line of (X++(Y..Y)) data. Returns (values, ended_with_dif)."""
    line = line.strip()
    if not line:
        return [], False
    if _AFFN_LINE.match(line):
        return _all_floats(line), False

    tokens = []  # (kind, text)  kind in {'val', 'dif', 'dup'}
    i, n = 0, len(line)
    while i < n:
        ch = line[i]
        if ch in " \t,;":
            i += 1
            continue
        if ch in _SQZ or ch in _DIF or ch in _DUP:
            kind = "val" if ch in _SQZ else ("dif" if ch in _DIF else "dup")
            first = _SQZ.get(ch, _DIF.get(ch, _DUP.get(ch)))
            j = i + 1
            while j < n and (line[j].isdigit() or line[j] == "."):
                j += 1
            digits = line[i + 1:j]
            if kind == "dup":
                text = str(first) + digits
            else:
                sign = "-" if first < 0 else ""
                text = sign + str(abs(first)) + digits
            tokens.append((kind, text))
            i = j
            continue
        if ch in "+-" or ch.isdigit() or ch == ".":
            m = _NUM_RE.match(line, i)
            if m:
                # an 'E'/'e' after a number is an SQZ character, not an exponent,
                # unless the whole line was AFFN (handled above)
                txt = re.match(r"[-+]?(?:\d+\.?\d*|\.\d+)", line[i:]).group(0)
                tokens.append(("val", txt))
                i += len(txt)
                continue
        i += 1  # unknown character: skip

    values = []
    last_kind = None
    last_step = None
    for kind, text in tokens:
        if kind == "val":
            v = float(text)
            values.append(v)
            last_kind, last_step = "val", v
        elif kind == "dif":
            d = float(text)
            values.append((values[-1] if values else 0.0) + d)
            last_kind, last_step = "dif", d
        else:  # dup: repeat the previous token (count includes the original)
            count = int(float(text))
            for _ in range(count - 1):
                if last_kind == "dif":
                    values.append(values[-1] + last_step)
                else:
                    values.append(last_step)
    ended_with_dif = bool(tokens) and (
        tokens[-1][0] == "dif" or (tokens[-1][0] == "dup" and last_kind == "dif"))
    return values, ended_with_dif


def _parse_xy_pairs(lines, width=2):
    nums = []
    for line in lines:
        nums.extend(_all_floats(line))
    usable = len(nums) - (len(nums) % width)
    arr = np.array(nums[:usable], dtype=float).reshape(-1, width)
    return arr[:, :2]


def _ntuple_field(meta, key, sym):
    """One field of an NTUPLES variable list (##FIRST=, ##LAST=, ##FACTOR=, ##UNITS=): the field of the
    variable with symbol sym (order of ##SYMBOL=, else X, Y, ...). None if missing or empty."""
    if key not in meta:
        return None
    fields = [f.strip() for f in meta[key].split(",")]
    syms = [s.strip().upper() for s in meta.get("SYMBOL", "").split(",")] if meta.get("SYMBOL") else []
    i = syms.index(sym) if sym in syms else {"X": 0, "Y": 1}.get(sym)
    if i is None or i >= len(fields) or not fields[i]:
        return None
    return fields[i]


def _factors(meta):
    """(XFACTOR, YFACTOR); in NTUPLES one ##FACTOR= field per variable."""
    xfac = _to_float(meta.get("XFACTOR", ""))
    yfac = _to_float(meta.get("YFACTOR", ""))
    if xfac is None:
        xfac = _to_float(_ntuple_field(meta, "FACTOR", "X") or "")
    if yfac is None:
        yfac = _to_float(_ntuple_field(meta, "FACTOR", "Y") or "")
    return (xfac or 1.0), (yfac or 1.0)


def _first_last(meta):
    """(FIRSTX, LASTX) in real units; in NTUPLES ##FIRST= / ##LAST= fields of X."""
    firstx = _to_float(meta.get("FIRSTX", ""))
    lastx = _to_float(meta.get("LASTX", ""))
    if firstx is None:
        firstx = _to_float(_ntuple_field(meta, "FIRST", "X") or "")
    if lastx is None:
        lastx = _to_float(_ntuple_field(meta, "LAST", "X") or "")
    return firstx, lastx


def _parse_xppyy(lines, meta):
    """(X++(Y..Y)) data; x axis rebuilt from FIRSTX/LASTX/NPOINTS if present.
    Each line starts with the x of its first y value; after a line ending in DIF form that first y is
    the y check (the last point of the line before), so the line's x belongs to that earlier point."""
    ys = []
    xs_line = []  # (x of the line as written, index of the point it belongs to)
    prev_dif = False
    for line in lines:
        vals, ended_dif = _decode_asdf_line(line)
        if len(vals) < 1:
            continue
        x0, yv = vals[0], vals[1:]
        if not yv:
            prev_dif = ended_dif
            continue
        if prev_dif and ys:
            xs_line.append((x0, len(ys) - 1))
            yv = yv[1:]  # DIF y-check value repeats the last y of previous line
        else:
            xs_line.append((x0, len(ys)))
        ys.extend(yv)
        prev_dif = ended_dif
    ys = np.array(ys, dtype=float)
    npts = len(ys)
    if npts == 0:
        return np.zeros((0, 2))
    xfac = _factors(meta)[0]
    firstx, lastx = _first_last(meta)
    nexp = _to_float(meta.get("NPOINTS", ""))
    nexp = int(round(nexp)) if nexp is not None and nexp > 0 else None
    if nexp is not None and nexp != npts:
        print("JCAMP-DX: %d points decoded, %d expected (NPOINTS); x positions taken from the lines" % (npts, nexp))
    if firstx is not None and lastx is not None and npts > 1 and nexp in (None, npts):
        xs = np.linspace(firstx, lastx, npts)
    else:
        deltax = _to_float(meta.get("DELTAX", ""))
        if deltax is None and firstx is not None and lastx is not None and nexp and nexp > 1:
            deltax = (lastx - firstx) / (nexp - 1)
        if deltax is None and len(xs_line) > 1:
            (xa, ia), (xb, ib) = xs_line[0], xs_line[-1]
            deltax = (xb - xa) * xfac / max(ib - ia, 1)
        deltax = deltax or 1.0
        if nexp is not None and nexp != npts:
            # point count off (a missing or extra y check): each line placed at its own x
            xs = np.empty(npts)
            bounds = [i for _, i in xs_line] + [npts]
            for (x0, i0), i1 in zip(xs_line, bounds[1:]):
                xs[i0:max(i1, i0 + 1)] = x0 * xfac + np.arange(max(i1, i0 + 1) - i0) * deltax
        else:
            start = firstx if firstx is not None else xs_line[0][0] * xfac
            xs = start + np.arange(npts) * deltax
    return np.column_stack([xs, ys])


def _data_form(value):
    v = value.replace(" ", "").upper()
    if "X++" in v:
        return "xppyy"
    if "XYW" in v:
        return "xyw"
    return "xy"


def _finish_block(meta, data_label, data_form, data_lines, page_meta=None):
    """Convert collected lines of one spectrum into an Nx2 array and metadata."""
    if data_label is None:
        return None
    m = dict(meta)
    if page_meta:
        m.update(page_meta)
    xfac, yfac = _factors(m)  # NTUPLES: one ##FACTOR= field per variable
    if data_form == "xppyy":
        arr = _parse_xppyy(data_lines, m)  # x already scaled (FIRSTX/LASTX are in real units)
        arr[:, 1] *= yfac
    else:
        arr = _parse_xy_pairs(data_lines, 3 if data_form == "xyw" else 2)
        arr[:, 0] *= xfac
        arr[:, 1] *= yfac
    if len(arr) == 0:
        return None
    arr = arr[np.all(np.isfinite(arr), axis=1)]
    arr = arr[np.argsort(arr[:, 0], kind="stable")]
    return {"meta": m, "data": arr}


def _retention_minutes(meta):
    """Midpoint retention time in minutes (None if unknown)."""
    for key in ("RETENTIONTIME", "PAGETIME", "SCANTIME", "TIME_RT"):
        if key in meta:
            vals = [abs(v) for v in _all_floats(meta[key])]
            if vals:
                t = float(np.mean(vals[:2]))
                units = (meta.get("SCANTIMEUNITS", "") + " " + meta.get("RETENTIONTIMEUNITS", "")).upper()
                if key == "PAGETIME" and not units.strip():
                    # NTUPLES: the unit of the page variable (##PAGE=T=60 -> the T field of ##UNITS=)
                    sym = meta.get("PAGE", "").partition("=")[0].strip().upper() if "=" in meta.get("PAGE", "") else ""
                    units = (_ntuple_field(meta, "UNITS", sym) if sym else None) or meta.get("UNITS", "")
                    units = units.upper()
                words = re.findall(r"[A-Z]+", units)
                # milliseconds first: "MILLISECONDS" also contains "SEC"
                if any(w.startswith("MILLISEC") or w in ("MS", "MSEC") for w in words):
                    t /= 60000.0
                elif any(w.startswith("SEC") or w in ("S", "SECONDS") for w in words):
                    t /= 60.0
                return t
    return None


def _polarity(meta):
    text = " ".join(meta.get(k, "") for k in ("IONIZATIONMODE", "POLARITY", "IONMODE", "MSPOLARITY")).upper()
    if "NEG" in text or text.strip().endswith("-") or "ESI-" in text or "APCI-" in text:
        return "Negative"
    return "Positive"


def parse_jcamp(path):
    """Return a list of spectra: [{'meta': {...}, 'data': Nx2 array}, ...]."""
    with open(path, "r", encoding="latin-1", errors="replace") as fh:
        raw_lines = fh.read().splitlines()

    spectra = []
    meta = {}
    data_label = None
    data_form = None
    data_lines = []
    in_ntuples = False
    page_meta = None

    def flush():
        nonlocal data_label, data_lines, data_form
        if data_label is not None:
            blk = _finish_block(meta, data_label, data_form, data_lines, page_meta)
            if blk is not None:
                spectra.append(blk)
        data_label, data_form, data_lines = None, None, []

    for raw in raw_lines:
        line = raw.split("$$", 1)[0].rstrip()
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("##"):
            label, _, value = stripped[2:].partition("=")
            key = _norm_label(label)
            value = value.strip()
            if key in _DATA_LABELS:
                flush()
                data_label, data_form, data_lines = key, _data_form(value), []
                continue
            if key == "TITLE":
                flush()
                if not in_ntuples:
                    meta = {}
                    page_meta = None
            elif key == "NTUPLES":
                flush()
                in_ntuples = True
            elif key == "PAGE":
                flush()
                page_meta = {"PAGE": value}
                t = _to_float(value)
                if t is not None:
                    page_meta["PAGETIME"] = str(t)
                continue
            elif key == "ENDNTUPLES":
                flush()
                in_ntuples = False
                page_meta = None
                continue
            elif key == "END":
                flush()
                continue
            if data_label is not None:
                flush()
            if page_meta is not None and in_ntuples:
                page_meta[key] = value
            else:
                meta[key] = value
        elif data_label is not None:
            data_lines.append(stripped)
    flush()

    for s in spectra:
        s["time"] = _retention_minutes(s["meta"])
        s["polarity"] = _polarity(s["meta"])
    return spectra


# --------------------------------------------------------------------------
# UniDec importer
# --------------------------------------------------------------------------
def _make_importer_class():
    from unidec.UniDecImporter.Importer import Importer

    class JCAMPImporter(Importer):
        """UniDec importer for JCAMP-DX mass spectra."""

        def __init__(self, file_path, **kwargs):
            super().__init__(file_path, **kwargs)
            spectra = parse_jcamp(file_path)
            if not spectra:
                raise IOError("No spectrum found in JCAMP-DX file: %s" % file_path)
            self.spectra = spectra
            self.data = [s["data"] for s in spectra]
            n = len(self.data)
            self.scans = np.arange(1, n + 1)
            times = [s["time"] for s in spectra]
            if any(t is None for t in times):
                times = [float(i) for i in range(n)] if n > 1 else [times[0] or 0.0]
            self.times = np.array(times, dtype=float)
            self.levels = np.ones(n, dtype=int)
            self.scan_range = [1, n]
            self.polarity = spectra[0]["polarity"]
            self.cdms_support = False
            self.imms_support = False
            self.chrom_support = n > 1
            m = spectra[0]["meta"]
            print("JCAMP-DX import: %d spectrum(s), %d points, %s, %s" % (
                n, len(self.data[0]), m.get("SPECTROMETERSYSTEM", "unknown system"), self.polarity))
            if "SCANNUMBER" in m or "RETENTIONTIME" in m:
                print("  scans %s, retention time %s %s" % (
                    m.get("SCANNUMBER", "?"), m.get("RETENTIONTIME", "?"), m.get("SCANTIMEUNITS", "")))

        def __len__(self):
            return len(self.data[0])

        def get_all_scans(self, *args, **kwargs):
            return self.data

        def get_single_scan(self, scan=None):
            if scan is None:
                return self.data[0]
            return self.data[self.get_scan_index(scan)]

        def get_avg_scan(self, scan_range=None, time_range=None):
            if len(self.data) == 1:
                return self.data[0].copy()
            return self.avg_fast(scan_range, time_range)

        def get_tic(self):
            tic = [float(np.sum(d[:, 1])) for d in self.data]
            return np.transpose([self.times, tic])

        def get_imms_avg_scan(self, *args, **kwargs):
            raise IOError("Ion mobility data is not available in JCAMP-DX files")

        def close(self):
            pass

    return JCAMPImporter


_INSTALLED = False


def install():
    """Register the JCAMP-DX importer with UniDec's ImporterFactory."""
    global _INSTALLED
    if _INSTALLED:
        return
    import importlib
    # (the package re-exports the class under the same name, so import the module explicitly)
    factory_module = importlib.import_module("unidec.UniDecImporter.ImporterFactory")

    JCAMPImporter = _make_importer_class()
    factory_module.JCAMPImporter = JCAMPImporter

    for ext in JCAMP_EXTENSIONS:
        if ext not in factory_module.recognized_types:
            factory_module.recognized_types.append(ext)

    original = factory_module.ImporterFactory.create_importer

    def create_importer(file_path, **kwargs):
        if os.path.splitext(str(file_path))[1].lower() in JCAMP_EXTENSIONS:
            return JCAMPImporter(file_path, **kwargs)
        return original(file_path, **kwargs)

    factory_module.ImporterFactory.create_importer = staticmethod(create_importer)
    _INSTALLED = True
    print("JCAMP-DX (.jdx/.dx/.jcamp) import enabled")


# --------------------------------------------------------------------------
# Stand-alone converter
# --------------------------------------------------------------------------
def convert_file(path, out_path=None):
    spectra = parse_jcamp(path)
    if not spectra:
        raise IOError("no spectrum found")
    if len(spectra) == 1:
        data = spectra[0]["data"]
    else:  # several spectra: merge exactly as UniDec does when opening the file
        from unidec.UniDecImporter.ImportTools import merge_spectra
        data = merge_spectra([s["data"] for s in spectra])
    if out_path is None:
        out_path = os.path.splitext(path)[0] + "_unidec.txt"
    np.savetxt(out_path, data, fmt="%.10g", delimiter="\t")
    return out_path, len(data), len(spectra)


def _collect(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            for name in sorted(os.listdir(p)):
                if os.path.splitext(name)[1].lower() in JCAMP_EXTENSIONS:
                    files.append(os.path.join(p, name))
        elif os.path.isfile(p):
            files.append(p)
    return files


def main(argv):
    files = _collect(argv)
    if not files:
        print("Drag .jdx files (or a folder containing them) onto Convert_JDX_to_TXT.bat")
        return 1
    ok = 0
    for f in files:
        try:
            out, npts, nspec = convert_file(f)
            print("OK   %s -> %s (%d points%s)" % (
                os.path.basename(f), os.path.basename(out), npts,
                "" if nspec == 1 else ", merged from %d spectra" % nspec))
            ok += 1
        except Exception as e:
            print("FAIL %s: %s" % (os.path.basename(f), e))
    print("\n%d of %d file(s) converted." % (ok, len(files)))
    return 0 if ok == len(files) else 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
