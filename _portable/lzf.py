"""
LZF decompression with the interface of python-lzf (lzf.decompress), for the
rainbow reader of Agilent MassHunter data (MSProfile.bin stores the profile
spectra of TOF and Q-TOF instruments LZF compressed; python-lzf has no wheel
for Windows and Python 3.12). Decompresses in the C++ core (msengine,
ms_lzf_decompress); decompress_py is the Python reference and the fallback.
"""
import ctypes

_FN = {"fn": None, "tried": False}


def _engine():
    if not _FN["tried"]:
        _FN["tried"] = True
        try:
            import msengine_py
            L = msengine_py.lib()
            if L is not None and hasattr(L, "ms_lzf_decompress"):
                f = L.ms_lzf_decompress
                f.argtypes = [ctypes.c_char_p, ctypes.c_long, ctypes.c_void_p, ctypes.c_long,
                              ctypes.POINTER(ctypes.c_long)]
                f.restype = ctypes.c_int
                _FN["fn"] = (f, L)
        except Exception:
            _FN["fn"] = None
    return _FN["fn"]


def decompress_py(data, max_len):
    """Python reference of ms_lzf_decompress: the bytes, or None when they
    would be longer than max_len (as python-lzf); ValueError for damaged data."""
    data = bytes(data)
    n_in = len(data)
    out = bytearray()
    ip = 0
    while ip < n_in:
        ctrl = data[ip]
        ip += 1
        if ctrl < 32:
            ln = ctrl + 1
            if ip + ln > n_in:
                raise ValueError("damaged LZF data (a literal run goes past the end)")
            if len(out) + ln > max_len:
                return None
            out += data[ip:ip + ln]
            ip += ln
        else:
            ln = ctrl >> 5
            ref = len(out) - ((ctrl & 0x1F) << 8) - 1
            if ln == 7:
                if ip >= n_in:
                    raise ValueError("damaged LZF data (truncated)")
                ln += data[ip]
                ip += 1
            if ip >= n_in:
                raise ValueError("damaged LZF data (truncated)")
            ref -= data[ip]
            ip += 1
            ln += 2
            if ref < 0:
                raise ValueError("damaged LZF data (a copy starts before the output)")
            if len(out) + ln > max_len:
                return None
            for k in range(ln):  # the copy may overlap itself
                out.append(out[ref + k])
    return bytes(out)


def decompress(data, max_len):
    """The decompressed bytes of data, or None when they would be longer than
    max_len (python-lzf's lzf.decompress)."""
    max_len = int(max_len)
    e = _engine()
    if e is None:
        return decompress_py(data, max_len)
    f, L = e
    data = bytes(data)
    buf = ctypes.create_string_buffer(max(max_len, 1))
    n = ctypes.c_long(0)
    rc = f(data, len(data), buf, max_len, ctypes.byref(n))
    if rc == -4:  # MS_NOT_FOUND: more than max_len bytes
        return None
    if rc:
        raise ValueError((L.ms_last_error() or b"").decode("utf-8", "replace") or "damaged LZF data")
    return buf.raw[:n.value]
