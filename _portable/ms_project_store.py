"""Portable project archives: typed JSON and non-pickled NumPy arrays, atomic writes."""
import hashlib,io,json,math,os,shutil,tempfile,time,zipfile
import numpy as np

FORMAT='MS Analysis project'
VERSION=1
NAME='session.msanalysis'
EXT='.msanalysis'
# Limits of the reader; the writer refuses a project that would break them.
MAX_SOURCES=100
MAX_ARRAYS=4094
MAX_HEADER=16*1024**2
MAX_TOTAL=2*1024**3
MAX_NDIM=8
STALE_TMP_S=600  # temporary files of a save older than this were left by a crash

def project_path(path):
    if os.path.isdir(path):path=os.path.join(path,NAME)
    return os.path.abspath(path)

def is_project(path):
    return str(path).lower().endswith(EXT) or (os.path.isdir(path) and os.path.isfile(os.path.join(path,NAME)))

class Encoder:
    def __init__(self,docs):self.docs={id(d):i for i,d in enumerate(docs)};self.arrays=[];self.ids={}
    def encode(self,v):
        if id(v) in self.docs:return ['doc',self.docs[id(v)]]
        if isinstance(v,np.ndarray):
            if v.dtype.hasobject:raise ValueError('Object arrays cannot be saved in a project')
            if id(v) not in self.ids:self.ids[id(v)]=len(self.arrays);self.arrays.append(v)
            return ['array',self.ids[id(v)]]
        if isinstance(v,np.generic):return self.encode(v.item())
        if v is None or isinstance(v,(str,bool,int)):return ['value',v]
        if isinstance(v,float):return ['float',v if math.isfinite(v) else str(v)]
        if isinstance(v,dict):return ['dict',[[self.encode(k),self.encode(x)] for k,x in v.items()]]
        if isinstance(v,(list,tuple,set,frozenset)):
            return [type(v).__name__,[self.encode(x) for x in v]]
        if type(v).__name__=='Calibration' and type(v).__module__=='hrms_calib':
            return ['calibration',self.encode(vars(v))]
        raise ValueError('Cannot save project value: '+type(v).__name__)
    def signature(self,state):
        tree=self.encode(state)
        identity=[(id(a),a.shape,str(a.dtype)) for a in self.arrays]
        return hashlib.sha256(json.dumps([tree,identity],sort_keys=True).encode()).hexdigest(),tree

def decode(tree,arrays,docs,depth=0):
    if depth>100:raise ValueError('Project nesting is excessive')
    if not isinstance(tree,list) or len(tree)!=2:raise ValueError('Invalid project value')
    tag,v=tree
    dec=lambda x:decode(x,arrays,docs,depth+1)
    if tag=='value':
        if v is not None and not isinstance(v,(str,bool,int)):raise ValueError('Invalid scalar')
        return v
    if tag=='float':return float(v)
    if tag in ('array','doc'):
        seq=arrays if tag=='array' else docs
        if type(v) is not int or not 0<=v<len(seq):raise ValueError('Invalid project reference')
        return seq[v]
    if tag=='dict':return {dec(k):dec(x) for k,x in v}
    if tag in ('list','tuple','set','frozenset'):
        return {'list':list,'tuple':tuple,'set':set,'frozenset':frozenset}[tag](dec(x) for x in v)
    if tag=='calibration':
        import hrms_calib
        attrs=dec(v)
        if not isinstance(attrs,dict) or attrs.get('model') not in range(6):raise ValueError('Invalid saved calibration')
        cal=hrms_calib.Calibration.__new__(hrms_calib.Calibration);cal.__dict__.update(attrs);return cal
    raise ValueError('Unknown project value type: '+str(tag))

def check(header,arrays):
    # The limits _header() and read() enforce, checked before anything is written.
    if header.get('kind') not in ('lcms','hrms'):raise ValueError('Unknown kind of analysis')
    n=len(header.get('sources') or [])
    if not 1<=n<=MAX_SOURCES:raise ValueError('An analysis can hold at most %d files (this one: %d)'%(MAX_SOURCES,n))
    if header.get('array_count')!=len(arrays):raise ValueError('Array count does not match the arrays')
    if len(arrays)>MAX_ARRAYS:raise ValueError('The analysis is too large to save: %d data arrays (at most %d)'%(len(arrays),MAX_ARRAYS))
    if not 0<=header.get('anchor',0)<n:raise ValueError('Invalid project anchor')
    for a in arrays:
        if a.dtype.hasobject or a.ndim>MAX_NDIM:raise ValueError('An array of the analysis cannot be saved')
    # .npy header of each array: at most a few hundred bytes
    if sum(a.nbytes for a in arrays)+512*len(arrays)>MAX_TOTAL-MAX_HEADER:
        raise ValueError('The analysis is too large to save (more than 2 GB of data)')

def _clean(folder):
    # Temporary files left by a save that crashed (other saves in progress are younger).
    try:names=os.listdir(folder)
    except OSError:return
    now=time.time()
    for n in names:
        if n.endswith('.tmp') and (n.startswith('.session-') or EXT+'.' in n.lower()):
            p=os.path.join(folder,n)
            try:
                if now-os.path.getmtime(p)>STALE_TMP_S:os.remove(p)
            except OSError:pass

def write(path,header,arrays,keep=None):
    """keep: None (the file replaced goes to .previous), or a suffix: the project replaced and its .previous are
    kept as <name>.<keep> and <name>.previous-<keep>, which later saves never rotate or replace."""
    check(header,arrays)
    meta=json.dumps(header,allow_nan=False).encode('utf-8')
    if len(meta)>MAX_HEADER:raise ValueError('The analysis is too large to save (its description exceeds 16 MB)')
    folder=os.path.dirname(os.path.abspath(path));os.makedirs(folder,exist_ok=True)
    _clean(folder)
    fd,tmp=tempfile.mkstemp(prefix='.session-',suffix='.tmp',dir=folder);os.close(fd)
    try:
        with zipfile.ZipFile(tmp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1) as z:
            z.writestr('project.json',meta)
            for i,a in enumerate(arrays):
                with z.open('arrays/%d.npy'%i,'w',force_zip64=True) as out:np.lib.format.write_array(out,a,allow_pickle=False)
        with open(tmp,'r+b') as f:f.flush();os.fsync(f.fileno())
        if keep:
            if os.path.isfile(path):_keep(path,path+'.'+keep)
            if os.path.isfile(path+'.previous'):_keep(path+'.previous',path+'.previous-'+keep)
        elif os.path.isfile(path):
            previous=path+'.previous'
            pending=previous+'.tmp';shutil.copyfile(path,pending);os.replace(pending,previous)
        os.replace(tmp,path)
    finally:
        if os.path.isfile(tmp):os.remove(tmp)

def _keep(src,name):
    # a copy of src under name (name-2, name-3, ... when an earlier one is there)
    dst=name;i=1
    while os.path.exists(dst) and i<100:i+=1;dst='%s-%d'%(name,i)
    pending=dst+'.tmp';shutil.copyfile(src,pending);os.replace(pending,dst)

def _header(z):
    if z.getinfo('project.json').file_size>MAX_HEADER:raise ValueError('Project metadata is too large')
    h=json.loads(z.read('project.json'))
    if not isinstance(h,dict) or h.get('format')!=FORMAT:raise ValueError('This is not an MS Analysis project')
    if h.get('version')!=VERSION:
        app=h.get('app_version')
        raise ValueError('Saved by a newer MS Analysis%s; update to open it'%(' (%s)'%app if isinstance(app,str) else '')
                         if isinstance(h.get('version'),int) and h['version']>VERSION else 'Unsupported project version')
    if h.get('kind') not in ('lcms','hrms') or not isinstance(h.get('sources'),list) or not 1<=len(h['sources'])<=MAX_SOURCES:
        raise ValueError('Invalid project source list')
    if type(h.get('array_count')) is not int or not 0<=h['array_count']<=MAX_ARRAYS:raise ValueError('Invalid project array count')
    if type(h.get('anchor',0)) is not int or not 0<=h.get('anchor',0)<len(h['sources']):raise ValueError('Invalid project anchor')
    return h

def peek(path):
    with zipfile.ZipFile(project_path(path)) as z:return _header(z)

def read(path):
    path=project_path(path)
    with zipfile.ZipFile(path) as z:
        infos=z.infolist()
        if len(infos)>MAX_ARRAYS+2 or sum(i.file_size for i in infos)>MAX_TOTAL:raise ValueError('Project archive is too large')
        h=_header(z)
        arrays=[]
        for i in range(h.get('array_count',0)):
            data=z.read('arrays/%d.npy'%i);stream=io.BytesIO(data)
            version=np.lib.format.read_magic(stream)
            if version==(1,0):shape,order,dtype=np.lib.format.read_array_header_1_0(stream)
            elif version==(2,0):shape,order,dtype=np.lib.format.read_array_header_2_0(stream)
            else:raise ValueError('Unsupported project array format')
            size=math.prod(shape)*dtype.itemsize
            if dtype.hasobject or size>len(data)-stream.tell() or len(shape)>MAX_NDIM:raise ValueError('Invalid project array')
            arrays.append(np.load(io.BytesIO(data),allow_pickle=False))
        return h,arrays

def _candidates(src,folder):
    rel=src.get('relative')  # None when the raw data was on another drive or share
    return ([os.path.abspath(os.path.join(folder,rel))] if isinstance(rel,str) else [])+[src['absolute']]

def resolve_sources(header,path):
    folder=os.path.dirname(project_path(path));paths=[]
    for src in header['sources']:
        found=next((p for p in _candidates(src,folder) if os.path.exists(p)),None)
        if found is None:raise FileNotFoundError('Raw data is missing: '+src['absolute'])
        paths.append(found)
    return paths

def contains(header,path,raw):
    """True when one of the sources of the project at path is the raw data file raw."""
    folder=os.path.dirname(project_path(path))
    key=os.path.normcase(os.path.abspath(raw.rstrip('\\/')))
    for src in header.get('sources') or []:
        try:
            if any(os.path.normcase(os.path.abspath(str(p).rstrip('\\/')))==key for p in _candidates(src,folder)):return True
        except (TypeError,ValueError,KeyError,AttributeError):continue
    return False

# ---------------------------------------------------------------- lock
# <project>.lock holds the process (id and start time) and the window that writes the project. A second
# window (of this or another process) does not write it while that window is open; a lock of a process
# that ended (a crash) is taken over.
_HELD={}

def _started(pid):
    if os.name!='nt':
        try:
            with open('/proc/%d/stat'%pid) as f:return int(f.read().rsplit(')',1)[1].split()[19])
        except Exception:return None
    import ctypes
    from ctypes import wintypes
    k=ctypes.WinDLL('kernel32',use_last_error=True)
    k.OpenProcess.restype=wintypes.HANDLE
    h=k.OpenProcess(0x1000,False,int(pid))  # PROCESS_QUERY_LIMITED_INFORMATION
    if not h:return -1 if ctypes.get_last_error()==5 else None  # access denied: it runs
    try:
        code=wintypes.DWORD()
        if k.GetExitCodeProcess(h,ctypes.byref(code)) and code.value!=259:return None  # 259: still active
        t=[wintypes.FILETIME() for _ in range(4)]
        if not k.GetProcessTimes(h,*[ctypes.byref(x) for x in t]):return -1
        return (t[0].dwHighDateTime<<32)|t[0].dwLowDateTime
    finally:k.CloseHandle(h)

def _holder(lp):
    """The owner of a lock file when it is alive, else None."""
    try:
        with open(lp,encoding='utf-8') as f:info=json.load(f)
        pid,started,owner=int(info['pid']),info.get('started'),str(info['owner'])
    except (OSError,ValueError,KeyError,TypeError):
        try:return 'unknown' if time.time()-os.path.getmtime(lp)<10 else None  # being written, or damaged
        except OSError:return None
    if pid==os.getpid():return owner if _HELD.get(os.path.normcase(lp))==owner else None
    now=_started(pid)
    if now is None:return None
    return owner if now==-1 or started in (None,now) else None

def lock(path,owner):
    """True when the window owner may write the project at path (the lock is taken, or it has it)."""
    lp=project_path(path)+'.lock';key=os.path.normcase(lp)
    if _HELD.get(key)==owner:return True
    for _ in range(2):
        try:fd=os.open(lp,os.O_WRONLY|os.O_CREAT|os.O_EXCL)
        except FileExistsError:
            if _holder(lp) is not None:return False
            try:os.remove(lp)
            except OSError:return False
            continue
        except OSError:return True  # no lock possible here (read only); the save reports its own error
        with os.fdopen(fd,'w',encoding='utf-8') as f:json.dump(dict(pid=os.getpid(),started=_started(os.getpid()),owner=owner),f)
        _HELD[key]=owner
        return True
    return False

def locked(path,owner):
    """True when another window holds the lock of the project at path."""
    lp=project_path(path)+'.lock';key=os.path.normcase(lp)
    if _HELD.get(key)==owner or not os.path.exists(lp):return False
    return _holder(lp) not in (None,owner)

def unlock(path,owner):
    lp=project_path(path)+'.lock';key=os.path.normcase(lp)
    if _HELD.get(key)!=owner:return
    del _HELD[key]
    try:os.remove(lp)
    except OSError:pass
