"""Transport for the native MaxEnt pipeline; no Python numerical fallback.

msmaxent.dll performs preprocessing, calls the preserved msengine.dll solver,
and calculates display-resolution broadening, peak tables and relative values.
Python only transports settings, arrays, progress and result dictionaries.
"""
import ctypes as C
import os
import numpy as np

_LIB = None
PD = C.POINTER(C.c_double)


def _lib(E):
    global _LIB
    if _LIB is not None:
        return _LIB
    path = os.path.join(os.path.dirname(__file__), "msengine", "msmaxent.dll")
    if not os.path.isfile(path):
        path = os.path.join(os.path.dirname(__file__), "msmaxent.dll")
    if not os.path.isfile(path):
        raise RuntimeError("Native MaxEnt component is missing: msmaxent.dll")
    L = C.CDLL(path)
    L.mx_abi.restype = C.c_uint
    L.mx_params_size.restype = C.c_uint
    if L.mx_abi() != 1 or L.mx_params_size() != C.sizeof(E.MaxentParams):
        raise RuntimeError("Native MaxEnt component has an incompatible parameter interface")
    L.mx_error.restype = C.c_char_p
    L.mx_run.argtypes = [C.c_wchar_p, PD, PD, C.c_long, C.POINTER(E.MaxentParams), C.c_double,
                         E.PROGRESS, C.c_void_p, C.POINTER(C.c_void_p)]
    L.mx_run.restype = C.c_int
    L.mx_free.argtypes = [C.c_void_p]
    L.mx_array.argtypes = [C.c_void_p, C.c_int, C.POINTER(PD), C.POINTER(PD)]
    L.mx_array.restype = C.c_long
    L.mx_peaks.argtypes = [C.c_void_p, C.POINTER(PD)]
    L.mx_peaks.restype = C.c_long
    for name in ("mx_summary", "mx_metrics"):
        getattr(L, name).argtypes = [C.c_void_p]
        getattr(L, name).restype = PD
    L.mx_notes.argtypes = [C.c_void_p]
    L.mx_notes.restype = C.c_char_p
    _LIB = L
    return L


def _integer(value, name):
    n = int(value)
    if n != value:
        raise ValueError("%s must be a whole number" % name)
    return n


def run(spec, p, progress=None):
    import msengine_py as E
    if not E.available():
        raise RuntimeError("MaxEnt requires the native engine. Enable the C++ engine or restore msengine.dll.")
    L = _lib(E)
    spec = np.ascontiguousarray(spec, dtype=np.float64)
    if spec.ndim != 2 or spec.shape[1] != 2:
        raise ValueError("A mass spectrum must have two columns: m/z and intensity")
    x, y = np.ascontiguousarray(spec[:, 0]), np.ascontiguousarray(spec[:, 1])
    q = E.MaxentParams()
    q.z_lo, q.z_hi = (_integer(v, "Charge") for v in p["z_range"])
    q.mass_lo, q.mass_hi = p["mass_range"]
    q.mass_step = float(p.get("mass_step", 1.0))
    # Explicit carrier settings and the original signed-proton convention.
    carrier = p.get("adduct_mass")
    try:
        carrier = float(carrier) if carrier not in (None, "") else 0.0
    except (TypeError, ValueError):
        carrier = 0.0
    if carrier == 0:
        carrier = -1.007276467 if (p.get("sign", 1) or 1) < 0 else 1.007276467
    q.adduct_mass = float(carrier)
    rng = p.get("mz_range")
    q.mz_lo, q.mz_hi = rng if rng else (0, 0)
    q.min_intensity = float(p.get("min_int", 0) or 0)
    q.min_intensity_pct = int(p.get("min_int_unit") == "pct")
    q.peak_width = float(p.get("peak_width", 0) or 0)
    res = float(p.get("resolution", 0) or 0)
    q.resolution = res if p.get("tof", False) or res > 0 else -1.0
    q.resolved_isotopes = int(p.get("isotopes", "envelope") == "resolved")
    for field, default in (("rounds", 10), ("min_rounds", 4), ("iterations", 50)):
        setattr(q, field, _integer(p.get(field, default), field))
    q.peak_window = float(p.get("peak_window", 10.0))
    q.peak_threshold = float(p.get("peak_thresh", .1))
    q.baseline = int(bool(p.get("baseline")))
    q.baseline_width = float(p.get("baseline_width", 15.0))
    q.envelope_nodes = _integer(p.get("envelope_nodes", 0) or 0, "Envelope nodes")
    q.noise_model = _integer(p.get("noise_model", 0) or 0, "Noise model")
    q.chi2_target = float(p.get("chi2_target", 0) or 0)
    callback_error = []

    def callback(user, i, n):
        if progress is not None:
            try:
                progress(i, n, float("nan"))
            except Exception as ex:
                callback_error.append(ex)
                return 0
        return 1

    out = C.c_void_p()
    cb = E.PROGRESS(callback)
    rc = L.mx_run(E._path(), x.ctypes.data_as(PD), y.ctypes.data_as(PD), len(x), C.byref(q),
                  float(p.get("display_resolution", 0) or 0), cb, None, C.byref(out))
    try:
        if callback_error:
            raise callback_error[0]
        if rc:
            raise ValueError((L.mx_error() or b"Native MaxEnt failed").decode("utf-8", "replace"))

        def array(kind):
            a, b = PD(), PD()
            n = L.mx_array(out, kind, C.byref(a), C.byref(b))
            if n < 0:
                raise RuntimeError("Invalid native MaxEnt output")
            return np.column_stack((np.ctypeslib.as_array(a, (n,)), np.ctypeslib.as_array(b, (n,)))) if n else np.empty((0, 2))

        rows = PD()
        n = L.mx_peaks(out, C.byref(rows))
        peaks = []
        for i in range(n):
            v = rows[i * 12:(i + 1) * 12]
            peak = dict(mass=v[0], height=v[1], area=v[2], apex=v[3], rel=v[5], frac=v[6])
            if v[4]:
                peak["n_iso"] = int(v[4])
                peak["z"] = "most abundant %.3f" % v[3]
            if not np.isnan(v[7]):
                peak["apex_dec"] = v[7]
                peak["apex_source"] = "m/z data" if v[10] else "deconvolution"
            if not np.isnan(v[8]):
                peak["apex2"], peak["apex2_rel"] = v[8], v[9]
            peaks.append(peak)
        s, m = L.mx_summary(out), L.mx_metrics(out)
        info = dict(base=s[0], base_mz=s[1], threshold=s[2], points=int(s[3]), points_used=int(s[4]),
                    peaks_used=int(s[5]), mz=(s[6], s[7]), value=s[8], unit=p.get("min_int_unit", "abs"))
        input_data = array(3)
        isotope = array(4)
        return dict(method="Maximum entropy", mass=array(0), fit=array(1), zdist=array(2), peaks=peaks,
                    data=input_data, input_data=input_data, data_info=info, r2=m[0],
                    history=[(0.0, m[1])] * int(m[2]), notes=L.mx_notes(out).decode("utf-8", "replace"),
                    isotope_peaks=[(float(mass), float(height)) for mass, height in isotope] if len(isotope) else None,
                    engine="C++ (native MaxEnt pipeline)")
    finally:
        if out.value:
            L.mx_free(out)
