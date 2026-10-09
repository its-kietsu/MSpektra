"""ctypes transport to native polymer calculations; no Python numerical fallback."""
import ctypes as C
from pathlib import Path
import numpy as np
PD=C.POINTER(C.c_double)
PI=C.POINTER(C.c_int)
_LIB=None
class Summary(C.Structure):
    _fields_=[('used',C.c_int)]+[(k,C.c_double) for k in ('mn','mw','dispersity')]
class Candidate(C.Structure):
    _fields_=[('mass',C.c_double),('score',C.c_double),('support',C.c_int),('single',C.c_int)]
def library():
    global _LIB
    if _LIB is None:
        root=Path(__file__).resolve().parent
        path=next((p for p in (root/'msengine'/'mspolymer.dll',root/'native'/'mspolymer.dll') if p.is_file()),None)
        if path is None:raise RuntimeError('Native polymer component is missing: mspolymer.dll')
        dll=C.CDLL(str(path))
        for name in ('poly_abi','poly_summary_size','poly_candidate_size'):getattr(dll,name).restype=C.c_uint
        if dll.poly_abi()!=1 or dll.poly_summary_size()!=C.sizeof(Summary) or dll.poly_candidate_size()!=C.sizeof(Candidate):
            raise RuntimeError('Native polymer component is incompatible')
        dll.poly_error.restype=C.c_char_p
        dll.poly_distribution.argtypes=[PD,PD,PD,C.c_int,C.c_int,C.POINTER(Summary),PD,PD]
        dll.poly_assign.argtypes=[PD,C.c_int,C.c_double,C.c_double,C.c_double,PD,PD,PD,PD,PI]
        dll.poly_offsets.argtypes=[PD,C.c_int,C.c_double,C.c_double,C.POINTER(Candidate),C.c_int]
        dll.poly_repeats.argtypes=[PD,C.c_int,C.c_double,C.c_double,C.c_double,C.POINTER(Candidate),C.c_int]
        for name in ('poly_distribution','poly_assign','poly_offsets','poly_repeats'):getattr(dll,name).restype=C.c_int
        _LIB=dll
    return _LIB
def array(x):
    a=np.ascontiguousarray(x,dtype=np.float64)
    if a.ndim!=1:raise ValueError('Polymer inputs must be one dimensional')
    return a
def ptr(x):return x.ctypes.data_as(PD)
def checked(code):
    if code<0:raise ValueError(library().poly_error().decode('utf-8','replace'))
    return code
def distribution(mass,amount,response=None,basis=0):
    m,a=array(mass),array(amount);r=array(response if response is not None else np.ones(len(m)))
    if len(m)!=len(a) or len(m)!=len(r):raise ValueError('Mass, signal and response counts must match')
    if basis not in (0,1):raise ValueError('Choose number or mass weighting')
    summary=Summary();nf=np.empty_like(m);wf=np.empty_like(m)
    checked(library().poly_distribution(ptr(m),ptr(a),ptr(r),len(m),basis,C.byref(summary),ptr(nf),ptr(wf)))
    return summary,nf,wf
def assign(mass,repeat,ends,tolerance=.05):
    m=array(mass);outputs=[np.empty_like(m) for i in range(4)];ok=np.empty(len(m),dtype=np.int32)
    checked(library().poly_assign(ptr(m),len(m),float(repeat),float(ends),float(tolerance),
                                 *(ptr(a) for a in outputs),ok.ctypes.data_as(PI)))
    return dict(zip(('dp','predicted','delta','ppm','matched'),outputs+[ok]))
def candidates(mass,repeat=None,lo=10,hi=1000,tolerance=.05):
    m=array(mass);out=(Candidate*20)()
    if repeat is None:count=checked(library().poly_repeats(ptr(m),len(m),float(lo),float(hi),float(tolerance),out,20))
    else:count=checked(library().poly_offsets(ptr(m),len(m),float(repeat),float(tolerance),out,20))
    return [dict(mass=q.mass,score=q.score,support=q.support,single=q.single) for q in out[:count]]
