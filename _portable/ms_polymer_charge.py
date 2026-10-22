"""Transport for native polymer charge assignments and unique mass groups."""
import ctypes as C
from pathlib import Path
import numpy as np
import ms_polymer as N

class Params(C.Structure):
    _fields_=[(k,C.c_double) for k in ('carrier','mz_min','mz_max','threshold_percent',
        'isotope_tolerance','merge_da','merge_ppm')]+[(k,C.c_int) for k in
        ('profile','charge_min','charge_max','sign','min_isotopes','mass_mode')]
class Peak(C.Structure):
    _fields_=[(k,C.c_double) for k in ('mass','mz','area','height','first_mz','last_mz','spacing_error','response')]+[
        (k,C.c_int) for k in ('charge','isotopes','state','alternate')]
class Neutral(C.Structure):
    _fields_=[(k,C.c_double) for k in ('mass','area','height')]+[(k,C.c_int) for k in ('count','representative')]
_LIB=None
DLL='mspolymer_charge5.dll'
class Unavailable(RuntimeError):
    """The native polymer library cannot be loaded (missing, quarantined or of another version)."""
STATES={0:'Assigned',1:'Uncertain',2:'Unassigned',3:'Manual charge',4:'Deconvoluted',5:'Outside range'}
DEFAULTS=dict(charge_min=1,charge_max=3,carrier=1.007276467,carrier_label='H+ (proton)',
              threshold_percent=5.,isotope_tolerance=.01,min_isotopes=3,mass_mode=0,
              merge_da=.02,merge_ppm=5.,profile=True)

def library():
    global _LIB
    if _LIB is None:
        path=Path(__file__).resolve().parent/'msengine'/DLL
        try:dll=C.CDLL(str(path))
        except OSError:
            raise Unavailable('Polymer tools not available: the file %s %s'%(
                DLL,'is missing' if not path.is_file() else 'cannot be opened'))
        for name in ('abi','params_size','peak_size','neutral_size'):
            getattr(dll,'poly_charge_'+name).restype=C.c_uint
        if (dll.poly_charge_abi()!=5 or dll.poly_charge_params_size()!=C.sizeof(Params) or
            dll.poly_charge_peak_size()!=C.sizeof(Peak) or dll.poly_charge_neutral_size()!=C.sizeof(Neutral)):
            raise Unavailable('Polymer tools not available: the file %s is of another version'%DLL)
        dll.poly_charge_error.restype=C.c_char_p
        dll.poly_charge_find.argtypes=[N.PD,N.PD,C.c_int,C.POINTER(Params),C.POINTER(Peak),C.c_int,
                                      C.POINTER(N.DirectMember),C.c_int,N.PI]
        dll.poly_charge_merge.argtypes=[C.POINTER(Peak),N.PI,C.c_int,C.c_double,C.c_double,
                                       C.POINTER(Neutral),C.c_int,N.PI]
        dll.poly_charge_resolve.argtypes=[C.POINTER(Peak),C.POINTER(N.DirectMember),C.c_int,C.c_int,
                                         C.c_double,C.c_int,C.c_double,N.PI]
        dll.poly_charge_pick.argtypes=[N.PD,N.PD,C.c_int,C.POINTER(Params),C.c_double,C.c_int,C.c_int,
                                      C.POINTER(Peak),C.POINTER(N.DirectMember),C.c_int,N.PI]
        for name in ('find','merge','resolve','pick'):getattr(dll,'poly_charge_'+name).restype=C.c_int
        _LIB=dll
    return _LIB

def available():
    """'' when the native library loads, else the message to show."""
    try:library();return ''
    except Unavailable as ex:return str(ex)

def whole(value,what='Charge and isotope settings must be whole numbers'):
    """An int from an int, a NumPy integer or a float without fraction (settings read back from files)."""
    if isinstance(value,bool):raise ValueError(what)
    if isinstance(value,(int,np.integer)):return int(value)
    if isinstance(value,(float,np.floating)) and float(value).is_integer():return int(value)
    raise ValueError(what)

def params(options,sign):
    ints={k for k,t in Params._fields_ if t is C.c_int}
    return Params(**{k:(int(bool(options[k])) if k=='profile' else int(sign) if k=='sign' else
                        whole(options[k]) if k in ints else float(options[k])) for k,_ in Params._fields_})

def checked(value):
    if value<0:raise ValueError(library().poly_charge_error().decode('utf-8','replace'))
    return value

def packed(row):
    values={k:row.get(k,1. if k=='response' else 0) for k,_ in Peak._fields_}
    return Peak(**{k:0 if v is None else v for k,v in values.items()})

def unpack(q):
    return {k:getattr(q,k) for k,_ in Peak._fields_}

def find(spectrum,options,sign=1):
    options={**DEFAULTS,**options}
    spec=np.asarray(spectrum,dtype=np.float64)
    if spec.ndim!=2 or spec.shape[1]!=2:raise ValueError('Choose an m/z and intensity spectrum')
    p=params(options,sign)
    x,y=N.array(spec[:,0]),N.array(spec[:,1])
    out=(Peak*5000)();capacity=max(1,min(len(x),200000));links=(N.DirectMember*capacity)();count=C.c_int()
    number=checked(library().poly_charge_find(N.ptr(x),N.ptr(y),len(x),C.byref(p),out,5000,
                                             links,capacity,C.byref(count)))
    result=[dict(unpack(q),members=[],include=q.state in (0,3),response=1.) for q in out[:number]]
    for m in links[:count.value]:
        result[m.species]['members'].append({k:getattr(m,k) for k in ('isotope','mz','area','height')})
    return result

def merge(rows,da=.02,ppm=5.,charge=None):
    if not rows:return [],[]
    if len(rows)>5000:raise ValueError('Use at most 5000 charge assignments')
    peaks=(Peak*len(rows))(*(packed(r) for r in rows))
    included=np.asarray([int(bool(r.get('include',True)) and (charge is None or r['charge']==charge))
                         for r in rows],dtype=np.int32)
    mapping=np.empty(len(rows),dtype=np.int32);output=(Neutral*len(rows))()
    count=checked(library().poly_charge_merge(peaks,included.ctypes.data_as(N.PI),len(rows),float(da),float(ppm),
                                              output,len(rows),mapping.ctypes.data_as(N.PI)))
    groups=[dict((k,getattr(q,k)) for k,_ in Neutral._fields_) for q in output[:count]]
    for group in groups:group['rows']=[]
    for i,j in enumerate(mapping):
        if j>=0:groups[j]['rows'].append(i)
    return groups,mapping.tolist()

def resolve(row,charge,options):
    charge=whole(charge,'Choose a whole charge state')
    members=row.get('members',[])
    if not members:raise ValueError('No measured peaks for this assignment')
    native=(N.DirectMember*len(members))(*(N.DirectMember(0,i,m['mz'],m['area'],m['height']) for i,m in enumerate(members)))
    q=packed(row);mask=(C.c_int*len(members))()
    checked(library().poly_charge_resolve(C.byref(q),native,len(members),charge,float(options['carrier']),
        whole(options.get('mass_mode',0)),float(options.get('isotope_tolerance',.01)),mask))
    selected=dict(row,**unpack(q),include=True,members=[]);remaining=[]
    for i,m in enumerate(members):
        if mask[i]:selected['members'].append(dict(m,isotope=len(selected['members'])))
        else:remaining.append(dict(mass=0.,mz=m['mz'],area=m['area'],height=m['height'],first_mz=m['mz'],
            last_mz=m['mz'],spacing_error=0.,response=row.get('response',1.),charge=0,isotopes=1,state=2,alternate=0,
            members=[dict(m,isotope=0)],include=False))
    return selected,remaining

def pick(spectrum,target,charge,options,isotopes=False):
    charge=whole(charge,'Choose a whole charge state')
    options={**DEFAULTS,**options}
    spec=np.asarray(spectrum,dtype=np.float64)
    if spec.ndim!=2 or spec.shape[1]!=2:raise ValueError('Choose a measured spectrum')
    p=params(options,options.get('sign',1))
    x,y=N.array(spec[:,0]),N.array(spec[:,1]);out=Peak();members=(N.DirectMember*5000)();count=C.c_int()
    checked(library().poly_charge_pick(N.ptr(x),N.ptr(y),len(x),C.byref(p),float(target),int(charge),int(bool(isotopes)),
                                      C.byref(out),members,5000,C.byref(count)))
    return dict(unpack(out),include=True,members=[{k:getattr(m,k) for k in ('isotope','mz','area','height')} for m in members[:count.value]])


