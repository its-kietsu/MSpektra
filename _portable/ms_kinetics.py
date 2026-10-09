"""Array transport to native kinetic least squares. No Python fitting fallback."""
import ctypes as C
from pathlib import Path
import numpy as np

MODELS = ('Linear trend', 'First-order decay', 'First-order growth')
EQUATIONS = ('y = C + slope × (X − X₀)', 'y = C + A × exp[−k × (X − X₀)]',
             'y = C + A × (1 − exp[−k × (X − X₀)])')
PD = C.POINTER(C.c_double)
_LIB = None

class Result(C.Structure):
    _fields_ = [(k, C.c_int) for k in ('model', 'count', 'parameters', 'dof', 'warnings')] + [
        ('origin', C.c_double), ('values', C.c_double * 3), ('errors', C.c_double * 3)] + [
        (k, C.c_double) for k in ('sse', 'rmse', 'r2', 'half_life', 'half_life_error')]

def library():
    global _LIB
    if _LIB is None:
        root = Path(__file__).resolve().parent
        paths = [root/'msengine'/'mskinetics.dll', root/'native'/'mskinetics.dll']
        path = next((p for p in paths if p.is_file()), None)
        if path is None:
            raise RuntimeError('mskinetics.dll is missing')
        dll = C.CDLL(str(path))
        dll.kin_abi.restype = dll.kin_result_size.restype = C.c_uint
        if dll.kin_abi() != 1 or dll.kin_result_size() != C.sizeof(Result):
            raise RuntimeError('mskinetics.dll does not match this version')
        dll.kin_error.restype = C.c_char_p
        dll.kin_fit.argtypes = [PD, PD, C.c_int, C.c_int, C.c_int, C.POINTER(Result), PD, PD]
        dll.kin_fit.restype = C.c_int
        dll.kin_predict.argtypes = [C.POINTER(Result), PD, C.c_int, PD]
        dll.kin_predict.restype = C.c_int
        _LIB = dll
    return _LIB

def array(values):
    a = np.ascontiguousarray(values, dtype=np.float64)
    if a.ndim != 1:
        raise ValueError('Fit inputs must be one-dimensional')
    return a

def pointer(a):
    return a.ctypes.data_as(PD)

def fit(x, y, model=1, offset=True):
    x, y = array(x), array(y)
    if len(x) != len(y):
        raise ValueError('X and area counts must match')
    if model not in range(len(MODELS)):
        raise ValueError('Choose a valid fitting model')
    dll = library()
    result = Result()
    fitted, residuals = np.empty_like(y), np.empty_like(y)
    if dll.kin_fit(pointer(x), pointer(y), len(x), model, int(bool(offset)), C.byref(result),
                   pointer(fitted), pointer(residuals)):
        raise ValueError(dll.kin_error().decode('utf-8', 'replace'))
    return result, fitted, residuals

def predict(result, x):
    x = array(x)
    if not len(x):
        return x.copy()
    out = np.empty_like(x)
    dll = library()
    if dll.kin_predict(C.byref(result), pointer(x), len(x), pointer(out)):
        raise ValueError(dll.kin_error().decode('utf-8', 'replace'))
    return out

def warnings(result):
    messages = []
    if result.warnings & 1:
        messages.append('Rate at the search limit.')
    if result.warnings & 2:
        messages.append('Uncertainties unavailable.')
    if result.warnings & 4:
        messages.append('Rate uncertainty as large as the rate.')
    return messages
