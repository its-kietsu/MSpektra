"""The polymer analysis of HRMS Analysis: a tile below the mass spectrum (the column of tiles fits the window
while it is open). Numerical processing is in C++ libraries."""
import copy
import math
import threading
import uuid
import wx
import numpy as np
import unilcms as U
import ms_polymer as N
import ms_polymer_charge as Q
import ms_polymer_widgets as PW
from unidec_theme import C,ui_font
from ms_polymer_results import PolymerResults,capture_source,file_key,same_file,rebind,FOCUS

COLOURS=(C['accent'],'#168579','#9D58B1','#9B6A20')  # charges +1, +2, +3, +4 (and on)
SERIES_COLOURS=('#2A62C4','#D97706','#7C3AED','#0E9F6E','#DB2777','#0891B2','#65A30D','#9B6A20')
GREY='#AEB6C2'
RB_H=24  # (DIP: a row of small ribbon buttons)
MASS_NAMES=['First isotope peak','Tallest isotope peak','Mean of the isotope peaks']

def charge_colour(z):return COLOURS[(abs(int(z))-1)%len(COLOURS)] if z else GREY
def number(value):
    return '%.4f'%value if isinstance(value,(float,np.floating)) and math.isfinite(value) else '—' if value is None else str(value)
def signed(z):return ('+%d'%z if z>0 else '−%d'%-z) if z else '—'
def height_text(v):return '%.1e'%v if v else '0'

def carriers(sign):
    """Charge carriers as in the deconvolution settings (Na⁺, ...): [(label, mass per charge)], then Custom mass."""
    import deconv_tab as D
    return [(label,mass) for label,mass in D.carriers(sign)]

# carrier names of 4.56 (saved analyses and presets): the same carrier by its mass
def carrier_index(choices,label,mass=None):
    i=next((i for i,(l,m) in enumerate(choices) if l==label),None)
    if i is None and mass is not None:
        i=next((i for i,(l,m) in enumerate(choices) if m is not None and abs(m-float(mass))<1e-6),None)
    return i

# Known repeat units and end groups (sum of both ends, Da) for the names of series.
REPEATS=[(44.0262147,'PEG'),(58.0418648,'PPG'),(74.0187569,'PDMS'),(72.0211294,'PLA'),(114.0680795,'PCL'),
         (100.0524295,'PMMA'),(104.0626,'PS'),(86.0367794,'PMA'),(56.0626,'PB'),(42.04695,'PP'),(28.0313,'PE')]
ENDS=[(18.0105647,'H / OH'),(32.0262147,'CH₃ / OH'),(46.0418648,'C₂H₅ / OH'),(2.0156501,'H / H'),
      (16.0313001,'CH₃ / H'),(30.0469502,'CH₃ / CH₃'),(60.0211294,'H / OCOCH₃'),
      (62.0367794,'CH₃ / OCH₂OH')]

def series_name(repeat,ends,taken=()):
    """'PEG · H / OH' for a known repeat unit and end groups (ends counted modulo the repeat), otherwise the
    masses."""
    base=next((name for m,name in REPEATS if abs(m-repeat)<.006),None)
    end=None
    if repeat>0:
        r=(ends%repeat)
        end=next((name for m,name in ENDS if min(abs(r-m%repeat),repeat-abs(r-m%repeat))<.01),None)
    if base=='PEG' and end=='CH₃ / OH':base='mPEG'
    name=('%s · %s'%(base,end) if base and end else '%s · ends %.2f Da'%(base,ends) if base else
          '%.2f Da · ends %.2f Da'%(repeat,ends))
    out,k=name,2
    while out in taken:out='%s (%d)'%(name,k);k+=1
    return out

FIELD_DEFAULTS=dict(repeat='',ends='',tolerance='0.05',search_min='10',search_max='1000')
FIELD_KEYS=('repeat','ends','tolerance','search_min','search_max')

def _method_dict(options,carrier_chosen,profile_explicit,charge_min,charge_max,fields,weight,basis):
    """The method values (presets, undo) of the workspace state given."""
    out={k:v for k,v in options.items() if k in Q.DEFAULTS and k!='profile'}
    out.update(carrier_label=options['carrier_label'] if carrier_chosen else '',
               spectrum_type=('profile' if options['profile'] else 'centroid') if profile_explicit else 'auto',
               charge_min=int(charge_min),charge_max=int(charge_max),weight=int(weight),basis=int(basis))
    out.update({k:fields[k] for k in FIELD_KEYS})
    return out

def _method_parts(values,options=None,profile_explicit=False):
    """(options, values, profile explicit) after method values were applied."""
    values=dict(values);options=dict(Q.DEFAULTS if options is None else options)
    if 'spectrum_type' in values:
        profile_explicit=values['spectrum_type']!='auto'
        if profile_explicit:values['profile']=values['spectrum_type']=='profile'
        else:values.pop('profile',None)
    if 'profile' in values:profile_explicit=True
    options.update({k:v for k,v in values.items() if k in Q.DEFAULTS})
    return options,values,profile_explicit

def default_method(tab):
    """The method of a workspace as it opens on this tab: undo compares it with the state before the workspace
    existed, so opening it is no step."""
    saved=dict(getattr(tab,'_polymer_method',{}))
    options,values,explicit=_method_parts(saved,profile_explicit='profile' in saved)
    fields={k:str(values[k]) if k in values else FIELD_DEFAULTS[k] for k in FIELD_KEYS}
    choice=lambda k:values[k] if values.get(k) in (0,1) else 0
    zmin=min(100,max(1,int(options['charge_min'])));zmax=min(100,max(1,int(options['charge_max'])))
    return _method_dict(options,bool(saved.get('carrier_label')),explicit,zmin,zmax,fields,choice('weight'),choice('basis'))

def method_get(tab):
    dock=getattr(tab,'polymer_results',None)
    if isinstance(dock,PolymerWorkspace):return dock.method_values()
    return dict(getattr(tab,'_polymer_method',{}))
def method_set(tab,values):
    tab._polymer_method={**getattr(tab,'_polymer_method',{}),**values}
    dock=getattr(tab,'polymer_results',None)
    if isinstance(dock,PolymerWorkspace):dock.set_method(tab._polymer_method)

def ensure_workspace(tab):
    old=getattr(tab,'polymer_results',None)
    if isinstance(old,PolymerWorkspace):return old
    legacy=copy.deepcopy(old.snapshot()) if old is not None and old.series else None
    if old is not None:
        for values in (tab.panels,getattr(tab,'_polymer_frames',[])):
            if old in values:values.remove(old)
        old.Hide();old.Destroy()
    dock=PolymerWorkspace(tab);tab.polymer_results=dock
    tab.extras['polymer']=(dock,None,3.5);tab.panels.append(dock)
    if legacy:dock.import_legacy(legacy)
    else:
        try:  # (its state as it opens is the reference: no undo step)
            import undo
            m=undo.manager(tab)
            if m is not None:m.rebase([k for k,d in m.domains.items() if k.startswith('polymer:') and d.tab and d.tab() is tab])
        except Exception:pass
    return dock

def open_workspace(tab,page=0,use_deconvolution=False):
    try:
        dock=ensure_workspace(tab)
        dock.visible=True;dock.refresh_inputs();dock.set_page(page)
        dock._watch();dock.fit_layout(True);tab._layout_rows();dock.Show()
        dock.refill();dock.redraw()
        if callable(getattr(tab,'keep_in_view',None)):tab.keep_in_view([dock])  # (a small window: scrolled to it)
        if use_deconvolution:dock.use_deconvolution()
        return dock
    except (ValueError,RuntimeError,OSError) as ex:
        tab.status(str(ex));return None


ROW_FLOATS=('mass','mz','area','height','first_mz','last_mz','spacing_error','response')
ROW_INTS=('charge','isotopes','state','alternate','include','deconv_index')

def pack_rows(rows):
    """Rows of a run as arrays (a project holds these as data files, not in its description: a run of 5000
    rows with their isotope peaks fits easily)."""
    n=len(rows);f=np.full((n,len(ROW_FLOATS)),np.nan);g=np.full((n,len(ROW_INTS)),-1,dtype=np.int64)
    members=[];extra=[];none_mz=[]
    for i,q in enumerate(rows):
        for j,k in enumerate(ROW_FLOATS):
            v=q.get(k)
            if v is None:
                if k=='mz':none_mz.append(i)
                continue
            f[i,j]=float(v)
        for j,k in enumerate(ROW_INTS):
            if k in q and q[k] is not None:g[i,j]=int(bool(q[k])) if k=='include' else int(q[k])
        for m in q.get('members') or []:
            members.append((i,m.get('isotope',0),m['mz'],m['area'],m['height']))
        more={k:v for k,v in q.items() if k not in ROW_FLOATS and k not in ROW_INTS and k!='members'}
        if more:extra.append([i,more])
    return dict(kind='packed_rows',floats=f,ints=g,members=np.asarray(members,dtype=np.float64).reshape(-1,5),
                extra=extra,none_mz=none_mz)

def unpack_rows(packed):
    if not isinstance(packed,dict) or packed.get('kind')!='packed_rows':return packed
    f=np.asarray(packed['floats'],dtype=np.float64);g=np.asarray(packed['ints'],dtype=np.int64);rows=[]
    for i in range(len(f)):
        q={k:float(f[i,j]) for j,k in enumerate(ROW_FLOATS) if not np.isnan(f[i,j])}
        for j,k in enumerate(ROW_INTS):
            if g[i,j]>=0 or k in ('charge','alternate'):
                q[k]=bool(g[i,j]) if k=='include' else int(g[i,j])
        q['members']=[]
        rows.append(q)
    for i in packed.get('none_mz',[]):rows[int(i)]['mz']=None
    for m in np.asarray(packed['members'],dtype=np.float64).reshape(-1,5):
        rows[int(m[0])]['members'].append(dict(isotope=int(m[1]),mz=float(m[2]),area=float(m[3]),height=float(m[4])))
    for i,more in packed.get('extra',[]):rows[int(i)].update(more)
    return rows

def pack_matches(matches):
    keys=sorted(int(k) for k in matches)
    get=lambda k:matches[k] if k in matches else matches[str(k)]
    return dict(kind='packed_matches',rows=np.asarray(keys,dtype=np.int64),
                n=np.asarray([int(get(k)['n']) for k in keys],dtype=np.int64),
                values=np.asarray([(get(k)['delta'],get(k)['ppm']) for k in keys],dtype=np.float64).reshape(-1,2))

def unpack_matches(m):
    if isinstance(m,dict) and m.get('kind')=='packed_matches':
        return {int(i):dict(n=int(n),delta=float(d),ppm=float(p)) for i,n,(d,p) in zip(m['rows'],m['n'],np.asarray(m['values']).reshape(-1,2))}
    return {int(i):v for i,v in (m or {}).items()}

def pack_groups(groups):
    rows=[r for g in groups for r in g.get('rows',[])];ends=np.cumsum([len(g.get('rows',[])) for g in groups]) if groups else []
    return dict(kind='packed_groups',values=np.asarray([(g['mass'],g['area'],g['height']) for g in groups],dtype=np.float64).reshape(-1,3),
                counts=np.asarray([(g.get('count',0),g.get('representative',-1)) for g in groups],dtype=np.int64).reshape(-1,2),
                rows=np.asarray(rows,dtype=np.int64),ends=np.asarray(ends,dtype=np.int64))

def unpack_groups(groups):
    if not(isinstance(groups,dict) and groups.get('kind')=='packed_groups'):return groups
    out=[];start=0
    for (m,a,h),(c,r),end in zip(np.asarray(groups['values']).reshape(-1,3),np.asarray(groups['counts']).reshape(-1,2),groups['ends']):
        out.append(dict(mass=float(m),area=float(a),height=float(h),count=int(c),representative=int(r),
                        rows=[int(i) for i in groups['rows'][start:int(end)]]));start=int(end)
    return out

def pack_run(run,shared):
    """A run for a project: arrays instead of long lists; a spectrum kept once (shared: the spectra already
    stored, the file's own included)."""
    out=dict(run);out.pop('chains',None)
    out['rows']=pack_rows(run['rows'])
    out['families']=[dict(f,matches=pack_matches(f['matches'])) for f in run.get('families',[])]
    d=run.get('distribution')
    if d:out['distribution']=dict(d,groups=pack_groups(d['groups']))
    source=dict(run['source']);columns=[]
    for c in source.get('columns',[]):
        spec=c.get('spec')
        if isinstance(spec,np.ndarray):
            same=next((a for a in shared if a.shape==spec.shape and a.dtype==spec.dtype and np.array_equal(a,spec)),None)
            if same is None:shared.append(spec)
            else:c=dict(c,spec=same)
        columns.append(c)
    source['columns']=columns;out['source']=source
    return out

def unpack_run(run):
    run=dict(run);run['rows']=unpack_rows(run['rows']);run['chains']=run['rows']
    run['families']=[dict(f,matches=unpack_matches(f.get('matches'))) for f in run.get('families',[])]
    d=run.get('distribution')
    if d:run['distribution']=dict(d,groups=unpack_groups(d['groups']))
    return run


class Settings(wx.Dialog):
    """Settings of Find species. inputs: the spectra shown (a choice when there are several)."""
    def __init__(self,parent,options,inputs=(),index=0,search=('10','1000')):
        super().__init__(parent,title='Polymer settings',style=wx.DEFAULT_DIALOG_STYLE)
        self.options=dict(options);self.controls={}
        self.SetFont(ui_font(9));self.SetBackgroundColour(C['bg'])
        grid=wx.FlexGridSizer(cols=2,hgap=self.FromDIP(12),vgap=self.FromDIP(9));grid.AddGrowableCol(1)
        self.inputs=list(inputs);self.spectrum=None
        if len(self.inputs)>1:
            self.spectrum=wx.Choice(self,choices=[d['label'] for d in self.inputs])
            self.spectrum.SetSelection(min(max(0,index),len(self.inputs)-1))
            grid.Add(wx.StaticText(self,label='Spectrum'),0,wx.ALIGN_CENTER_VERTICAL);grid.Add(self.spectrum,1,wx.EXPAND)
        fields=[('mz_min','m/z from'),('mz_max','m/z to'),('threshold_percent','Smallest peak (% of the tallest)'),
                ('min_isotopes','Isotope peaks needed'),('isotope_tolerance','Peak position tolerance (m/z)'),
                ('merge_da','Same mass within (Da)'),('merge_ppm','Same mass within (ppm)')]
        for key,label in fields:
            c=wx.TextCtrl(self,value='%.15g'%options[key],size=self.FromDIP(wx.Size(120,-1)));self.controls[key]=c
            grid.Add(wx.StaticText(self,label=label),0,wx.ALIGN_CENTER_VERTICAL);grid.Add(c,1,wx.EXPAND)
        self.search=[wx.TextCtrl(self,value=str(v),size=self.FromDIP(wx.Size(120,-1))) for v in search]
        for label,c in zip(('Repeat units from (Da)','Repeat units to (Da)'),self.search):
            grid.Add(wx.StaticText(self,label=label),0,wx.ALIGN_CENTER_VERTICAL);grid.Add(c,1,wx.EXPAND)
        self.profile=wx.Choice(self,choices=['Profile','Centroid']);self.profile.SetSelection(0 if options['profile'] else 1)
        self.mass=wx.Choice(self,choices=MASS_NAMES);self.mass.SetSelection(int(options['mass_mode']))
        for label,c in [('Data',self.profile),('Mass from',self.mass)]:
            grid.Add(wx.StaticText(self,label=label),0,wx.ALIGN_CENTER_VERTICAL);grid.Add(c,1,wx.EXPAND)
        box=wx.BoxSizer(wx.VERTICAL);box.Add(grid,0,wx.EXPAND|wx.ALL,self.FromDIP(16))
        box.Add(self.CreateButtonSizer(wx.OK|wx.CANCEL),0,wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM,self.FromDIP(14))
        self.SetSizerAndFit(box);self.CentreOnParent()
    def values(self):
        out=dict(self.options)
        for key,c in self.controls.items():
            try:v=float(c.GetValue().replace(',','.'))
            except ValueError:raise ValueError('Enter numbers in the settings')
            if not math.isfinite(v):raise ValueError('Enter finite settings')
            if key=='min_isotopes':
                if not v.is_integer():raise ValueError('Isotope peaks needed must be a whole number')
                v=int(v)
            out[key]=v
        out['profile']=self.profile.GetSelection()==0;out['mass_mode']=self.mass.GetSelection()
        if out['mz_min']<0 or out['mz_max']<=out['mz_min']:raise ValueError('Check the m/z range')
        if not 0<out['threshold_percent']<=100:raise ValueError('The smallest peak must be above 0 and at most 100 %')
        if not 0<out['isotope_tolerance']<.25:raise ValueError('The peak position tolerance must be above 0 and below 0.25')
        if not 2<=out['min_isotopes']<=20:raise ValueError('Choose 2 to 20 isotope peaks')
        if not 0<=out['merge_da']<=1 or not 0<=out['merge_ppm']<=1000:raise ValueError('Check the same mass tolerance')
        self.search_range()
        return out
    def search_range(self):
        try:lo,hi=[float(c.GetValue().replace(',','.')) for c in self.search]
        except ValueError:raise ValueError('Check the range of the repeat units')
        if not 0<lo<hi:raise ValueError('Check the range of the repeat units')
        return '%g'%lo,'%g'%hi
    def input_index(self):
        return self.spectrum.GetSelection() if self.spectrum is not None else 0

class PeakCharge(wx.Dialog):
    """Assign a charge to a measured peak (Other… of the charge menu)."""
    def __init__(self,parent,mz,sign,spectrum=None,charge=None,carrier=None):
        super().__init__(parent,title='Assign charge')
        self.SetFont(ui_font(9));self.SetBackgroundColour(C['bg'])
        self.sign=sign;self.choices=carriers(sign)
        grid=wx.FlexGridSizer(cols=2,hgap=self.FromDIP(12),vgap=self.FromDIP(10));grid.AddGrowableCol(1)
        self.charge=wx.SpinCtrl(self,min=1,max=100,initial=max(1,abs(int(charge or parent.z0.GetValue()))))
        self.carrier=wx.Choice(self,choices=[q[0] for q in self.choices])
        label,mass=carrier if carrier else (parent.options.get('carrier_label'),parent.options.get('carrier'))
        i=carrier_index(self.choices,label,mass)
        self.carrier.SetSelection(i if i is not None else 0)
        self.mass=wx.TextCtrl(self,value='%.15g'%(mass if mass is not None else parent.options['carrier']))
        for text,c in [('Charge (%s)'%('−' if sign<0 else '+'),self.charge),('Carrier',self.carrier),('Custom mass (Da)',self.mass)]:
            grid.Add(wx.StaticText(self,label=text),0,wx.ALIGN_CENTER_VERTICAL);grid.Add(c,1,wx.EXPAND)
        self.target=mz;self.spectrum=spectrum;self.options=dict(parent.current()['options'] if parent.current() else parent.options,sign=sign)
        self.isotopes=wx.CheckBox(self,label='With its isotope peaks');self.isotopes.SetValue(False)
        self.preview=wx.StaticText(self,label='')
        self.members=wx.ListCtrl(self,style=wx.LC_REPORT|wx.LC_SINGLE_SEL|wx.BORDER_SIMPLE,size=self.FromDIP(wx.Size(330,140)))
        self.members.InsertColumn(0,'m/z',wx.LIST_FORMAT_RIGHT,width=self.FromDIP(150))
        self.members.InsertColumn(1,'Height',wx.LIST_FORMAT_RIGHT,width=self.FromDIP(150))
        self.preview_result=None;self.preview_error=None
        box=wx.BoxSizer(wx.VERTICAL)
        box.Add(wx.StaticText(self,label='Peak at m/z %.4f'%mz),0,wx.ALL,self.FromDIP(14))
        box.Add(grid,0,wx.EXPAND|wx.LEFT|wx.RIGHT,self.FromDIP(14))
        box.Add(self.isotopes,0,wx.ALL,self.FromDIP(14))
        box.Add(self.preview,0,wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM,self.FromDIP(14))
        box.Add(self.members,1,wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM,self.FromDIP(14))
        box.Add(self.CreateButtonSizer(wx.OK|wx.CANCEL),0,wx.EXPAND|wx.ALL,self.FromDIP(14))
        self.SetSizer(box)
        self.carrier.Bind(wx.EVT_CHOICE,self.carrier_changed)
        self.isotopes.Bind(wx.EVT_CHECKBOX,self.update_preview)
        self.charge.Bind(wx.EVT_SPINCTRL,self.update_preview);self.mass.Bind(wx.EVT_TEXT,self.update_preview)
        self.carrier_changed()
    def carrier_changed(self,e=None):
        self.mass.Enable(self.choices[self.carrier.GetSelection()][1] is None)
        self.update_preview()
    def update_preview(self,e=None):
        self.preview_result=None;self.preview_error=None;expanded=self.isotopes.GetValue()
        self.members.DeleteAllItems();self.members.Show(expanded)
        if not expanded:self.preview.SetLabel('This peak only')
        else:
            try:
                charge,label,carrier,_=self.values()
                if self.spectrum is None:raise ValueError('No measured spectrum')
                q=Q.pick(self.spectrum,self.target,charge,dict(self.options,carrier=carrier),isotopes=True)
                self.preview_result=q
                self.preview.SetLabel('%d peaks, m/z %.4f to %.4f'%(q['isotopes'],q['first_mz'],q['last_mz'])
                    if q['isotopes']>1 else 'No isotope peaks at this charge: this peak only')
                for m in q['members']:
                    i=self.members.InsertItem(self.members.GetItemCount(),'%.4f'%m['mz'])
                    self.members.SetItem(i,1,height_text(m['height']))
            except (ValueError,RuntimeError,OSError) as ex:
                self.preview_error=str(ex);self.preview.SetLabel(self.preview_error)
        self.GetSizer().Layout();self.Fit();self.SetMinSize(self.FromDIP(wx.Size(380,-1)))
    def values(self):
        label,mass=self.choices[self.carrier.GetSelection()]
        if mass is None:
            try:mass=float(self.mass.GetValue().replace(',','.'))
            except ValueError:raise ValueError('Enter the carrier mass')
        if not math.isfinite(mass) or mass==0:raise ValueError('Enter a finite carrier mass')
        return self.sign*self.charge.GetValue(),label,mass,self.isotopes.GetValue()

class FlatPlot(U.PlotCard):
    """A plot inside the polymer tile: no edge and no title (the legend in its first row; the toggles above name
    it). The title stays empty when a saved view is put back (4.57 named these plots)."""
    title=property(lambda self:'',lambda self,value:None)
    def _paint(self,e):
        if self.pair is not None:return U.PlotCard._paint(self,e)
        U.paint_card(self,self._hh,U.inner_frame,self._draw_header,self._header_key())

class PolymerWorkspace(PolymerResults):
    show_empty=True
    BODY_DIP=300  # (height of the plot and the table: the plots with their axis titles)
    SIDE_DIP=12  # (the space inside the edge of the tile)
    LIST_DIP=184  # (a list of suggestions: seven rows and more by scrolling)
    def __init__(self,tab):
        wx.Panel.__init__(self,tab.rows,style=wx.BORDER_NONE)
        self.tab=tab;self.SetBackgroundColour(C['bg']);self.SetFont(ui_font(9))
        self.series=[];self.active=None;self.chain=None;self.visible=True;self.outputs=None
        self._bindings={};self._restoring_source=False;self._filling=False;self._busy=False
        self.inputs=[];self.input_index=0;self.shown=[];self.groups=[];self.mapping=[];self._source_signature=None
        self.options=dict(Q.DEFAULTS);self._method=dict(getattr(tab,'_polymer_method',{}))
        self._profile_explicit='profile' in self._method
        self._carrier_chosen=bool(self._method.get('carrier_label'))
        self._work_token=0;self._exact_fields={};self.unavailable='';self._layout=None
        self._said=None;self.charge_values=[None,'uncertain'];self._plot_groups=[];self._assigned_rows=[]
        self._plot=0;self._table=1  # plot: 0 neutral masses, 1 mass distribution; table: 0 series, 1 charge assignments
        self._base=None;self._base_key=None
        self.carriers=carriers(1)
        card=self.header=self.card=PW.Card(self);card.SetFont(ui_font(9));card.no_tile_menu=True
        self._build_header();self._build_finder();self._build_body()
        side,gap=self.FromDIP(self.SIDE_DIP),self.FromDIP(8);box=wx.BoxSizer(wx.VERTICAL)
        for w,top in [(self.flow,side),(PW.Line(card),gap),(self.finder,gap),(PW.Line(card),gap),(self.split,gap)]:
            box.Add(w,1 if w is self.split else 0,wx.EXPAND|wx.TOP,top)
        box.AddSpacer(side)
        frame=wx.BoxSizer(wx.HORIZONTAL);frame.Add(box,1,wx.EXPAND|wx.LEFT|wx.RIGHT,side)
        card.SetSizer(frame)
        outer=wx.BoxSizer(wx.VERTICAL);outer.Add(card,1,wx.EXPAND);self.SetSizer(outer)
        self.SetMinSize(self.FromDIP(wx.Size(100,200)))
        # the wheel over the tile scrolls the column of tiles
        for w in [card,self.flow,self.finder,self.split,self.plot_column,self.table_column]+self.groups_panels:
            w.Bind(wx.EVT_MOUSEWHEEL,lambda e,w=w:U.forward_wheel(w,e))
        for c in (self.mass_card,self.dist_card):c.Bind(wx.EVT_SHOW,lambda e,c=c:(e.Skip(),wx.CallAfter(self._shown_plot,c)))
        self.set_method(self._method);self.refresh_inputs();self._watch();self.set_views(0,1);self.Hide()

    # ------------------------------------------------------------------ building
    def _panel(self,parent):
        p=wx.Panel(parent,style=wx.BORDER_NONE);p.SetBackgroundColour(PW.WHITE);p.SetFont(ui_font(9))
        self.groups_panels.append(p);return p
    def _add(self,row,parent,items):
        for w,space in items:
            if isinstance(w,str):w=PW.label(parent,w)
            row.Add(w,0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(space))

    def _build_header(self):
        """The tool row: the run and Find species, then what is shown (charge, series); settings at the right. A
        narrow window puts the second part on a line of its own."""
        self.groups_panels=[]
        self.flow=PW.Flow(self.card,gap=9,vgap=8,lines=False);self.flow.on_height=self._height_changed
        g1=self._panel(self.flow);g2=self._panel(self.flow);g3=self._panel(self.flow)
        title=PW.label(g1,'Polymer analysis',bold=True)
        self.runs=wx.Choice(g1);self.runs.SetMinSize(self.FromDIP(wx.Size(170,-1)))
        self.z0=wx.SpinCtrl(g1,min=1,max=100,initial=1,size=self.FromDIP(wx.Size(54,-1)))
        self.z1=wx.SpinCtrl(g1,min=1,max=100,initial=3,size=self.FromDIP(wx.Size(54,-1)))
        self.carrier=wx.Choice(g1);self.carrier.SetMinSize(self.FromDIP(wx.Size(96,-1)))
        self.find_btn=U.flat(g1,'Find species','primary',handler=self.find_species)
        row=wx.BoxSizer(wx.HORIZONTAL)
        self._add(row,g1,[(title,14),(self.runs,18),('Charges',6),(self.z0,6),('to',6),(self.z1,18),('Carrier',6),
                          (self.carrier,14),(self.find_btn,0)])
        g1.SetSizer(row)
        # what is shown: the charge (everything below shows its masses) and the series
        self.charge=PW.KeyChoice(g2,self.filter_changed);self.charge.SetItems(['All charges'],['All · 0 masses'])
        self.family=wx.Choice(g2,choices=['All series']);self.family.SetSelection(0)
        row=wx.BoxSizer(wx.HORIZONTAL)
        self._add(row,g2,[(PW.label(g2,'Show',muted=True),6),(self.charge.choice,6),(self.family,0)])
        g2.SetSizer(row)
        self.settings_btn=U.flat(g3,'⚙',handler=self.settings,tooltip='Settings',padx=10)
        self.more_btn=U.flat(g3,'•••',handler=lambda e:self.more_menu(),tooltip='More',padx=10)
        row=wx.BoxSizer(wx.HORIZONTAL);self._add(row,g3,[(self.settings_btn,6),(self.more_btn,0)]);g3.SetSizer(row)
        self.tool_groups=[g1,g2];self.flow.right=g3;self.flow.set_items(self.tool_groups)
        self.carrier.Bind(wx.EVT_CHOICE,self.carrier_changed)
        self.runs.Bind(wx.EVT_CHOICE,lambda e:self.select_series(self.runs.GetSelection()))
        for c in (self.z0,self.z1):c.Bind(wx.EVT_SPINCTRL,self.charges_changed)
        self.family.Bind(wx.EVT_CHOICE,self.filter_changed)

    def _suggestions(self,parent,columns,kind):
        t=PW.Table(parent,'',columns,empty='Press Find',flat=True,note=True)
        t.SetMinSize(wx.Size(-1,self.FromDIP(self.LIST_DIP)))
        t.on_select=lambda i:self.use_candidate(i,kind);t.on_activate=lambda i:self.use_candidate(i,kind)
        t.extra_menu=lambda:[('Use this mass',lambda:self.use_candidate(t.selected()[0],kind) if t.selected() else None)]
        return t

    def _build_finder(self):
        """The series finder: repeat unit and end groups, each with its suggestions; the tolerance and Assign series,
        the weighting of the mass distribution and its result."""
        f=self.finder=wx.Panel(self.card,style=wx.BORDER_NONE);f.SetBackgroundColour(PW.WHITE);f.SetFont(ui_font(9))
        F=FIELD_DEFAULTS;size=self.FromDIP(wx.Size(76,-1));C_=PW.Column
        # repeat unit (the range searched is in the settings)
        left=self._panel(f)
        self.repeat=wx.TextCtrl(left,value=F['repeat'],size=size);self.repeat.SetToolTip('Da')
        self.lo=wx.TextCtrl(left,value=F['search_min'],size=self.FromDIP(wx.Size(60,-1)))
        self.hi=wx.TextCtrl(left,value=F['search_max'],size=self.FromDIP(wx.Size(60,-1)))
        for c in (self.lo,self.hi):c.Hide()
        self.find_repeat_btn=U.flat(left,'Find',height=26,padx=10,handler=self.find_repeats,tooltip='Find repeat units')
        self.repeat_list=self._suggestions(left,[C_('mass','Mass (Da)','right',80,1.),C_('support','Gaps','right',50,.5),
                                                 C_('single','Single gaps','right',70,.5)],'repeat')
        # end groups
        mid=self._panel(f)
        self.ends=wx.TextCtrl(mid,value=F['ends'],size=size);self.ends.SetToolTip('Da')
        self.find_ends_btn=U.flat(mid,'Find',height=26,padx=10,handler=self.find_offsets,tooltip='Find end groups')
        self.ends_list=self._suggestions(mid,[C_('mass','Mass (Da)','right',80,1.),C_('support','Masses','right',60,.5)],'offset')
        for p,label,field,btn,lst in ((left,'Repeat unit',self.repeat,self.find_repeat_btn,self.repeat_list),
                                      (mid,'End groups',self.ends,self.find_ends_btn,self.ends_list)):
            row=wx.BoxSizer(wx.HORIZONTAL);name=PW.label(p,label);name.SetFont(ui_font(9,600))
            self._add(row,p,[(name,8),(field,6),(btn,0)])
            box=wx.BoxSizer(wx.VERTICAL);box.Add(row,0);box.Add(lst,1,wx.EXPAND|wx.TOP,self.FromDIP(6));p.SetSizer(box)
        # tolerance and Assign series; the mass distribution
        right=self._panel(f)
        self.tol=wx.TextCtrl(right,value=F['tolerance'],size=self.FromDIP(wx.Size(56,-1)));self.tol.SetToolTip('Tolerance (Da)')
        self.assign_btn=U.flat(right,'Assign series','primary',height=26,padx=10,handler=self.assign_series)
        self.weight=wx.Choice(right,choices=['Area','Height']);self.weight.SetSelection(0)
        self.basis=wx.Choice(right,choices=['Per molecule','Per mass']);self.basis.SetSelection(0)
        self.weight.SetToolTip('Abundance from');self.basis.SetToolTip('Count')
        self.calc_btn=U.flat(right,'Calculate',height=26,padx=10,handler=self.calculate_distribution)
        self.summary=PW.Summary(right);self.summary.SetFont(ui_font(9,600))
        r1=wx.BoxSizer(wx.HORIZONTAL);name=PW.label(right,'±');name.SetFont(ui_font(9,600))
        self._add(r1,right,[(name,6),(self.tol,10),(self.assign_btn,0)])
        r2=wx.BoxSizer(wx.HORIZONTAL);name=PW.label(right,'Distribution');name.SetFont(ui_font(9,600))
        self._add(r2,right,[(name,8),(self.weight,6),(self.basis,8),(self.calc_btn,0)])
        box=wx.BoxSizer(wx.VERTICAL);box.Add(r1,0);box.Add(r2,0,wx.TOP,self.FromDIP(14))
        box.Add(self.summary,0,wx.TOP|wx.LEFT,self.FromDIP(10));right.SetSizer(box)
        row=wx.BoxSizer(wx.HORIZONTAL)
        for k,(p,share) in enumerate(((left,1.),(mid,.8),(right,1.))):
            if k:row.Add(PW.Line(f,wx.VERTICAL),0,wx.EXPAND|wx.LEFT|wx.RIGHT,self.FromDIP(14))
            row.Add(p,int(share*10),wx.EXPAND)
        f.SetSizer(row)
        for c in (self.repeat,self.ends,self.tol,self.lo,self.hi):c.Bind(wx.EVT_TEXT,self.fields_changed)
        for c in (self.repeat,self.ends):c.Bind(wx.EVT_KILL_FOCUS,self.format_field)
        self.weight.Bind(wx.EVT_CHOICE,self.invalidate_distribution);self.basis.Bind(wx.EVT_CHOICE,self.invalidate_distribution)

    def _column(self,parent):
        p=wx.Panel(parent,style=wx.BORDER_NONE);p.SetBackgroundColour(PW.WHITE);p.SetFont(ui_font(9))
        top=wx.BoxSizer(wx.HORIZONTAL);box=wx.BoxSizer(wx.VERTICAL)
        box.Add(top,0,wx.EXPAND|wx.BOTTOM,self.FromDIP(4));p.SetSizer(box)
        return p,top,box

    def _build_body(self):
        split=self.split=U.SplitBox(self.card,wx.HORIZONTAL,gap=22,min_size=150);split.SetBackgroundColour(PW.WHITE)
        split.SetMinSize(wx.Size(-1,self.FromDIP(self.BODY_DIP)))
        # the plot: neutral masses of the charge chosen, or the mass distribution
        col,top,box=self._column(split);self.plot_column=col
        self.plot_buttons=[PW.Toggle(col,'Masses',lambda:self.set_views(plot=0,user=True),icon='ms'),
                           PW.Toggle(col,'Distribution',lambda:self.set_views(plot=1,user=True),icon='pda')]
        for k,b in enumerate(self.plot_buttons):top.Add(b,0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(4))
        self.mass_card=FlatPlot(col,'','',mode='zoom');self.mass_card.no_tile_menu=True
        self.mass_card.on_pick=self.pick_mass;self.mass_card.on_view=self.draw_mass
        self.mass_card.on_menu=lambda x,y:[('Show all masses',self.mass_card.reset_view)]
        self.mass_card.image_name=lambda:self.image_name('neutral_masses')
        card=self.dist_card=FlatPlot(col,'','',mode='zoom');card.no_tile_menu=True
        card.on_menu=lambda x,y:[('Show all masses',card.reset_view)]
        card.image_name=lambda:self.image_name('mass_distribution')
        for c in (self.mass_card,self.dist_card):box.Add(c,1,wx.EXPAND)
        # the series assigned, or the charge assignments of every peak
        C_=PW.Column
        col,top,box=self._column(split);self.table_column=col
        self.table_buttons=[PW.Toggle(col,'Series',lambda:self.set_views(table=0,user=True),icon='list'),
                            PW.Toggle(col,'Charge assignments',lambda:self.set_views(table=1,user=True),icon='table')]
        for b in self.table_buttons:top.Add(b,0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(4))
        self.table_info=self._info_label(col)
        top.Add(self.table_info,1,wx.ALIGN_CENTER_VERTICAL|wx.LEFT,self.FromDIP(6))
        self.assign_table=PW.Table(col,'Series',[C_('series','Series','left',90,2.,swatch=True),C_('n','Units (n)','right',56,.4),
            C_('mass','Mass (Da)','right',80,1.),C_('error','Error (Da)','right',66,.6)],empty='Assign a series to list its chains',flat=True)
        self.assign_table.on_select=self.select_assignment
        self.assign_table.extra_menu=lambda:[('Remove this series',self.remove_series)] if self.family.GetSelection()>0 else []
        self.chain_table=PW.Table(col,'Charge assignments',[
            C_('use','Use','center',30,0,check=True),C_('mass','Mass (Da)','right',70,.6),C_('z','z','center',34,.2),
            C_('mz','m/z','right',70,.6),C_('height','Height','right',56,.3),C_('iso','Isotopes','right',62,.2),
            C_('series','Series','left',110,3.,swatch=True),C_('n','Units (n)','right',56,.4,shown=False,optional=True),
            C_('area','Area','right',60,.6,shown=False,optional=True),
            C_('response','Response factor','right',70,.5,shown=False,optional=True),
            C_('state','Assignment','left',80,.8,shown=False,optional=True)],
            multi=True,empty='Find species to list their charge states',flat=True)
        t=self.chain_table
        t.on_select=self.select_table;t.on_check=self.check_row;t.on_cell=self.table_cell;t.on_delete=self.delete_assignments
        t.on_activate=lambda n:self.open_charge_menu_for_row(self.shown[n]) if 0<=n<len(self.shown) else None
        t.extra_menu=lambda:[('Set charge',self.charge_items(self.chain)) if self.chain is not None else ('Set charge',None),
            ('Include or exclude',self.toggle_selected),('Response factor…',self.edit_response),(None,None),
            ('Delete selected peaks',self.delete_assignments),('Export Excel…',self.export)]
        for table in (self.assign_table,self.chain_table):
            table.on_info=lambda text,table=table:self._table_info(table,text);box.Add(table,1,wx.EXPAND)
        split.set_panes([self.plot_column,self.table_column],[1.1,1.])

    def _shown_plot(self,card):
        """A plot shown again: drawn at its size now."""
        if not self or self.IsBeingDeleted() or not card.IsShown():return
        card.place()
        if card is self.dist_card:self.draw_distribution()
        else:self.draw_mass()
    def steps_changed(self):
        """A control of the tool row or the finder changed its width (charge and series names, the result)."""
        for p in getattr(self,'groups_panels',[]):p.InvalidateBestSize();p.Layout()
        self.finder.Layout();self.flow.place()
    def _info_label(self,parent):
        t=wx.StaticText(parent,label='',style=wx.ST_ELLIPSIZE_END)
        t.SetBackgroundColour(PW.WHITE);t.SetForegroundColour(wx.Colour(C['muted']));t.SetMinSize(wx.Size(self.FromDIP(20),-1))
        return t
    def _set_info(self,label,text):
        if label.GetLabel()!=text:label.SetLabel(text);label.GetContainingSizer().Layout()
    def _table_info(self,table,text):
        if table is (self.assign_table,self.chain_table)[self._table]:self._set_info(self.table_info,text)

    def image_name(self,what):
        s=self.current()
        base=self.tab._spec_name(s['source']['event']) if s and self.tab.data else 'polymer'
        return '%s_%s'%(base,what)

    # ------------------------------------------------------------------ layout in HRMS Analysis
    def fit_layout(self,on):
        """While the analysis is open this tile has its own height (it does not fit the window): the chromatogram
        and the spectrum keep the heights they have without it (taller when the window has room) and the column
        scrolls when the window is small. Closing it gives back the tile heights of before."""
        rows=self.tab.rows
        keys=('frac','fit_share','fit_first','fit_floor','fit_label','mode_key','default_mode')
        if on and self._layout is None:
            self._layout={k:copy.copy(getattr(rows,k)) for k in keys}
            rows.frac={};rows.fit_share=self._fit_share;rows.fit_floor=0
            rows.fit_label='Plots fill the window'
            rows.mode_key,rows.default_mode='polymer_tile_mode','two'
        elif not on and self._layout is not None:
            for k,v in self._layout.items():setattr(rows,k,v)
            self._layout=None
        elif on:
            rows.frac={}  # (tile heights of another layout, such as those of a saved analysis)
        self._base_key=None
        try:rows.relayout();self.tab.scroller.refit();self.tab.scroller.Scroll(0,0)
        except (RuntimeError,AttributeError):pass
    def _base_heights(self):
        """{tile: height} of the other tiles as they are without this tile (the layout of before)."""
        rows=self.tab.rows;lay=self._layout
        if lay is None:return {}
        panes=[p for p in rows.panes if p is not self];weights=[w for p,w in zip(rows.panes,rows.weights) if p is not self]
        st=U.settings()
        key=(rows.view_h,rows.GetClientSize()[0],tuple(id(p) for p in panes),st.get(lay['mode_key']),st.get('tile_scale'),
             tuple(sorted((id(p),f) for p,f in lay['frac'].items())))
        if key==self._base_key:return self._base
        saved={k:getattr(rows,k) for k in tuple(lay)+('panes','weights')}
        try:
            for k,v in lay.items():setattr(rows,k,v)
            rows.panes,rows.weights=panes,weights
            hs=rows.heights()
        finally:
            for k,v in saved.items():setattr(rows,k,v)
        self._base,self._base_key=dict(zip(panes,hs)),key
        return self._base
    def _fit_share(self,pane):
        """(share, pixels): this tile its own height; the chromatogram its height of before; the spectrum at least a
        readable height and the room the window has left (never lower to make room for this tile: the column
        scrolls); the other tiles their heights of before."""
        tab=self.tab;rows=tab.rows
        if pane is self:return 0,self.tile_height()
        n=len(getattr(pane,'cards',[]) or [1])
        if pane in tab._chrom_tiles():  # (its height of before: never a thin strip)
            return 1e-6,self._base_heights().get(pane) or max(rows.min,self.FromDIP(150))*n
        if pane in tab._spec_tiles():
            return 1.0*n,max(rows.min,self.FromDIP(230 if n==1 else 170))*n
        h=self._base_heights().get(pane)
        return None if h is None else (0,h)
    def tile_height(self):
        """The height of this tile: the tool row (on two lines when the window is narrow), the series finder, the
        plot and the table."""
        w=self.GetSize()[0] if self.GetSize()[0]>50 else self.tab.rows.GetClientSize()[0]
        inner=max(1,w-2*self.FromDIP(self.SIDE_DIP))
        self.flow.SetMinSize(wx.Size(-1,self.flow.height_for(inner)))
        return self.card.GetSizer().CalcMin()[1]
    def _height_changed(self):
        """The tool row needs another height (a narrower window): the column is laid out again."""
        if not self or self.IsBeingDeleted() or not self.visible or self._layout is None:return
        try:self.tab.rows.relayout();self.tab.scroller.refit()
        except (RuntimeError,AttributeError):pass

    # ------------------------------------------------------------------ the plot and the table shown
    def page(self):
        """The section of 4.57 (saved analyses, the polymer menu): 0 charge states, 1 repeat units, 2 mass
        distribution."""
        return 2 if self._plot==1 else 1 if self._table==0 else 0
    def set_page(self,i):
        i=min(2,max(0,int(i)))
        self.set_views(*{0:(0,1),1:(0,0),2:(1,0)}[i])
    def views(self):return [self._plot,self._table]
    def set_views(self,plot=None,table=None,user=False):
        if plot is not None:self._plot=min(1,max(0,int(plot)))
        if table is not None:self._table=min(1,max(0,int(table)))
        for k,b in enumerate(self.plot_buttons):b.set_active(k==self._plot)
        for k,b in enumerate(self.table_buttons):b.set_active(k==self._table)
        self.mass_card.Show(self._plot==0);self.dist_card.Show(self._plot==1)
        self.assign_table.Show(self._table==0);self.chain_table.Show(self._table==1)
        for c in (self.plot_column,self.table_column):c.Layout()
        table=(self.assign_table,self.chain_table)[self._table];self._set_info(self.table_info,table.info_text())
        if self._plot==1:self.draw_distribution()
        else:self.draw_mass()
        if user:
            if plot==1 and self.current() and not self.current().get('distribution') and self.scoped_groups():
                self.calculate_distribution()
            self._changed()
    def section_changed(self,i):
        self.set_page(i)
        if i==2 and self.current() and not self.current().get('distribution') and self.scoped_groups():
            self.calculate_distribution()
        self._changed()

    # ------------------------------------------------------------------ state
    def current(self):
        return self.series[self.active] if self.active is not None and 0<=self.active<len(self.series) else None
    def say(self,text):
        """A short message in the status bar, cleared after a while (no old texts stay)."""
        self.tab.status(text);self._said=text
        def clear():
            try:
                bar=self.tab.frame.GetStatusBar()
                if bar and self._said==text and bar.GetStatusText(0)==text:self.tab.status('')
            except (RuntimeError,AttributeError):pass
        try:wx.CallLater(9000,clear)
        except Exception:pass
    def fail(self,ex):
        self.say(str(ex));return None
    def _changed(self):
        super()._changed()
        import undo
        undo.schedule(self.tab)
    def method_values(self):
        fields=dict(zip(FIELD_KEYS,(self.field_text(self.repeat),self.field_text(self.ends),self.tol.GetValue(),
                                    self.lo.GetValue(),self.hi.GetValue())))
        return _method_dict(self.options,self._carrier_chosen,self._profile_explicit,self.z0.GetValue(),
                            self.z1.GetValue(),fields,self.weight.GetSelection(),self.basis.GetSelection())
    def set_method(self,values):
        if 'carrier_label' in values:self._carrier_chosen=bool(values['carrier_label'])
        self.options,values,self._profile_explicit=_method_parts(values,self.options,self._profile_explicit)
        self.z0.SetValue(int(self.options['charge_min']));self.z1.SetValue(int(self.options['charge_max']))
        for key,c in [('repeat',self.repeat),('ends',self.ends),('tolerance',self.tol),('search_min',self.lo),('search_max',self.hi)]:
            if key in values:self.set_field(c,values[key])
        for key,c in [('weight',self.weight),('basis',self.basis)]:
            if values.get(key) in (0,1):c.SetSelection(values[key])
        self._method=self.method_values();self.tab._polymer_method=dict(self._method)
        if self.inputs:self.set_carriers()
    def charges_changed(self,e=None):
        self.options.update(charge_min=self.z0.GetValue(),charge_max=self.z1.GetValue())
        self.tab._polymer_method=self.method_values()
        if e:e.Skip()
    def fields_changed(self,e=None):
        if e:self._exact_fields.pop(e.GetEventObject(),None)
        self.tab._polymer_method=self.method_values()
        if self.current():self.current()['fields']=self.field_values()
        self._changed()
        if e:e.Skip()
    def set_field(self,control,value):
        raw=str(value);display=raw
        if control in (self.repeat,self.ends) and raw.strip():
            try:
                v=float(raw)
                if math.isfinite(v):display='%.4f'%v
            except ValueError:pass
        control.ChangeValue(display);self._exact_fields[control]=(display,raw)
    def field_text(self,control):
        value=control.GetValue();stored=self._exact_fields.get(control)
        return stored[1] if stored and value==stored[0] else value
    def format_field(self,e):
        c=e.GetEventObject();self.set_field(c,self.field_text(c));e.Skip()
    def field_values(self):return [self.field_text(c) for c in (self.repeat,self.ends,self.tol,self.lo,self.hi)]

    # ------------------------------------------------------------------ input spectrum and carrier
    def refresh_inputs(self):
        old=self.inputs[self.input_index]['label'] if 0<=self.input_index<len(self.inputs) else None
        self.inputs=list(self.tab.deconv_inputs())
        self.input_index=next((i for i,d in enumerate(self.inputs) if d['label']==old),0)
        self.source_changed()
    def source_changed(self,e=None):
        if self.inputs:
            d=self.inputs[min(self.input_index,len(self.inputs)-1)];spec=d['spec']
            self.options.update(mz_min=float(spec[0,0]),mz_max=float(spec[-1,0]))
            if self.options['mz_max']==self.options['mz_min']:
                self.options['mz_min']=max(0.,self.options['mz_min']-.01);self.options['mz_max']+=.01
            profile=getattr(self.tab.data,'has_profile',True)
            if callable(profile):profile=profile()
            if not self._profile_explicit:self.options['profile']=bool(profile)
            self.set_carriers()
        self.find_btn.Enable(bool(self.inputs) and not self._busy)
    def input(self):
        return self.inputs[min(self.input_index,len(self.inputs)-1)] if self.inputs else None
    def deconvolution_carrier(self,sign):
        """The carrier of the deconvolution settings (the one chosen for these data), or None."""
        try:
            cfg=self.tab.dec.cfg;lst=carriers(sign)
            k=int(cfg.get('adduct_neg' if sign<0 else 'adduct',0) or 0)
            return lst[k] if 0<=k<len(lst) and lst[k][1] is not None else None
        except Exception:return None
    def set_carriers(self):
        sign=(self.input() or {}).get('sign',1)
        self.carriers=carriers(sign)
        self.carrier.SetItems([q[0] if q[1] is not None else 'Custom mass…' for q in self.carriers])
        if self._carrier_chosen:
            i=carrier_index(self.carriers,self._method.get('carrier_label',self.options.get('carrier_label')),self.options.get('carrier'))
            if i is None or self.carriers[i][1] is None:
                i=len(self.carriers)-1  # (a custom mass)
        else:
            q=self.deconvolution_carrier(sign) or self.carriers[0]
            i=carrier_index(self.carriers,q[0],q[1]) or 0
            self.options.update(carrier=self.carriers[i][1],carrier_label=self.carriers[i][0])
        self.carrier.SetSelection(i)
        self.carrier.SetMinSize(wx.Size(max(self.FromDIP(96),self.carrier.GetBestSize()[0]),-1));self.header.Layout()
        self.carrier.SetToolTip('%.6f Da per charge'%self.options['carrier'] if self.carriers[i][1] is None else '')
    def carrier_changed(self,e=None):
        i=self.carrier.GetSelection()
        if i<0:return
        label,value=self.carriers[i]
        if value is None:
            with wx.TextEntryDialog(self,'Carrier mass per charge (Da, with its sign)','Charge carrier','%.6f'%self.options['carrier']) as dlg:
                if dlg.ShowModal()!=wx.ID_OK:self.set_carriers();return
                try:value=float(dlg.GetValue().replace(',','.'));assert math.isfinite(value) and value!=0
                except (ValueError,AssertionError):self.set_carriers();self.fail('Enter a finite carrier mass');return
        self._carrier_chosen=True
        self.options.update(carrier=value,carrier_label=label);self._method.update(carrier=value,carrier_label=label)
        self.tab._polymer_method=self.method_values()
    def settings(self,e=None):
        if not self.inputs:return self.fail('Show a spectrum first')
        dlg=Settings(self,self.options,self.inputs,self.input_index,(self.field_text(self.lo),self.field_text(self.hi)))
        try:
            while dlg.ShowModal()==wx.ID_OK:
                try:
                    self.options=dlg.values();self.input_index=dlg.input_index();self._profile_explicit=True
                    lo,hi=dlg.search_range();self.set_field(self.lo,lo);self.set_field(self.hi,hi);self.fields_changed()
                    self._method=self.method_values();self.tab._polymer_method=dict(self._method);self.set_carriers();self._changed();break
                except ValueError as ex:wx.MessageBox(str(ex),'Polymer settings',wx.OK|wx.ICON_INFORMATION,dlg)
        finally:dlg.Destroy()
    def more_items(self):
        s=self.current()
        return [('Use deconvoluted masses',self.use_deconvolution),(None,None),
                ('Export Excel…',self.export if s else None),('Remove analysis',self.remove_analysis if s else None),
                (None,None),('Close',self.hide_results)]
    def more_menu(self):
        m=U.build_menu(self,self.more_items());U.popup_menu(self.more_btn,m,wx.Point(0,self.more_btn.GetSize()[1]));m.Destroy()
    def analysis_options(self):
        d=self.input()
        if d is None:raise ValueError('Show a spectrum first')
        options=dict(self.options,charge_min=self.z0.GetValue(),charge_max=self.z1.GetValue(),
                     source_label=d['label'],source_event=int(d.get('e',0)),sign=d.get('sign',1))
        if options['charge_max']<options['charge_min']:raise ValueError('Check the charges')
        return d,options

    # ------------------------------------------------------------------ runs
    def find_species(self,e=None,synchronous=False):
        if self._busy:return None
        if len(self.series)>=100:return self.fail('Keep at most 100 polymer analyses in one file')
        if callable(getattr(self.tab,'busy_averaging',None)) and self.tab.busy_averaging():return self.fail('Wait for the spectrum to load')
        message=Q.available()
        if message:return self.fail(message)
        try:
            self.refresh_source_only();d,options=self.analysis_options()
            source=capture_source(self.tab,options);spec=np.array(d['spec'],copy=True);sign=d.get('sign',1)
        except (ValueError,RuntimeError) as ex:return self.fail(ex)
        self._work_token+=1;token=self._work_token
        self._set_busy(True);self.tab.status('Finding polymer species…')
        def work():
            try:result=Q.find(spec,options,sign);error=None
            except Exception as ex:result=None;error=str(ex)
            if synchronous:return self._finished(token,source,options,result,error)
            wx.CallAfter(self._finished,token,source,options,result,error)
        if synchronous:return work()
        threading.Thread(target=work,name='Polymer charges',daemon=True).start()
    def refresh_source_only(self):
        current=list(self.tab.deconv_inputs())
        if [d['label'] for d in current]!=[d['label'] for d in self.inputs]:
            self.refresh_inputs();raise ValueError('The spectrum changed. Find species again')
        self.inputs=current
    def _set_busy(self,busy):
        self._busy=busy
        for c in (self.z0,self.z1,self.carrier,self.find_btn,self.settings_btn):c.Enable(not busy)
    def _finished(self,token,source,options,rows,error):
        if not self or self.IsBeingDeleted() or token!=self._work_token:return None
        self._set_busy(False)
        if error:return self.fail(error)
        if not same_file(self.tab,source['file']):return self.fail('The data file changed. Find species again')
        col=next((c for c in self.tab.cols if c['e']==source['event']),None)
        saved=next((c['spec'] for c in source['columns'] if c['e']==source['event']),None)
        if col is None or col.get('spec') is None or not np.array_equal(col['spec'],saved):return self.fail('The spectrum changed. Find species again')
        record=dict(id=uuid.uuid4().hex,source=source,options=options,rows=rows,chains=rows,
                    families=[],fields=self.field_values(),weight=self.weight.GetSelection(),basis=self.basis.GetSelection())
        # Running the same settings on the same spectrum again replaces that run.
        old=next((i for i,s in enumerate(self.series) if s['options']==options and
                  s['source']['kind']=='direct' and s['source']['event']==source['event'] and
                  all(np.array_equal(a['spec'],b['spec']) for a,b in zip(s['source']['columns'],source['columns']))),None)
        if old is None:self.series.append(record);self.active=len(self.series)-1
        else:record['id']=self.series[old]['id'];self.series[old]=record;self.active=old
        self.chain=None
        self.visible=True;self._bindings.clear();self.charge.SetSelection(0);self.family.SetSelection(0)
        self._changed();self.select_series(self.active);self.tab._layout_rows()
        n=sum(q['state'] in (0,3) for q in rows)
        self.say('%d species with a charge, %d peaks without'%(n,len(rows)-n));return rows
    def use_deconvolution(self,e=None):
        if self._busy:return None
        if len(self.series)>=100:return self.fail('Keep at most 100 polymer analyses in one file')
        try:
            from ms_polymer_ui import source_rows
            rows,meta=source_rows(self.tab);source=capture_source(self.tab)
            options=dict(self.options);params=source['result'].get('params',{})
            if 'adduct_mass' in params:
                options['carrier']=float(params['adduct_mass'])
                lst=carriers(params.get('sign',1));i=carrier_index(lst,None,options['carrier'])
                options['carrier_label']=lst[i][0] if i is not None else 'Custom mass'
            if params.get('z_range'):
                a,b=params['z_range'];options.update(charge_min=min(abs(int(a)),abs(int(b))),charge_max=max(abs(int(a)),abs(int(b))))
            options['sign']=params.get('sign',1)
            for q in rows:q.update(charge=0,state=4,alternate=0,isotopes=0,members=[],spacing_error=0,response=1.)
            record=dict(id=uuid.uuid4().hex,source=source,options=options,rows=rows,chains=rows,
                        families=[],fields=self.field_values(),metadata=meta,weight=0,basis=0)
            self.series.append(record);self.active=len(self.series)-1;self.visible=True;self.chain=None
            self._changed();self.select_series(self.active);self.tab._layout_rows();return rows
        except (ValueError,RuntimeError) as ex:return self.fail(ex)
    def run_label(self,i,run):
        label=run['source'].get('label','') or ''
        what=label.split(': ',1)[-1]
        if len(self.inputs)>1 or len(getattr(self.tab,'cols',[]))>1:
            what='(%s) %s'%('−' if run['options'].get('sign',run['source'].get('sign',1))<0 else '+',what)
        if run['source'].get('kind')!='direct':what+=' · deconvoluted'
        return 'Run %d · %s'%(i+1,what)
    def select_series(self,index):
        if self._filling or not 0<=index<len(self.series):return
        self.active=index;self.chain=None;self._bindings.clear();s=self.current()
        self.options.update({k:v for k,v in s['options'].items() if k not in ('source_label','source_event')})
        self._profile_explicit=True;self._carrier_chosen=True
        self.z0.SetValue(int(self.options.get('charge_min',1)));self.z1.SetValue(int(self.options.get('charge_max',3)))
        self._method.update({k:v for k,v in self.options.items() if k in Q.DEFAULTS})
        for c,v in zip((self.repeat,self.ends,self.tol,self.lo,self.hi),s.get('fields',[])):self.set_field(c,v)
        self.weight.SetSelection(s.get('weight',0));self.basis.SetSelection(s.get('basis',0))
        self.show_source();self.refresh_inputs()
        self.input_index=next((i for i,d in enumerate(self.inputs) if d.get('e',0)==s['source'].get('event',0)),self.input_index)
        if self.inputs:self.set_carriers()
        self.rebuild_filters();self.refill();self.draw_distribution();self.redraw()
    def family_name(self,s,i):
        f=s['families'][i];names=[]
        for k,g in enumerate(s['families'][:i]):names.append(series_name(g['repeat'],g['ends'],names))
        return series_name(f['repeat'],f['ends'],names)
    def family_colour(self,i):return SERIES_COLOURS[i%len(SERIES_COLOURS)]
    def rebuild_filters(self):
        s=self.current();old=self.charge.GetStringSelection();family=self.family.GetSelection()
        rows=s['rows'] if s else []
        charges=sorted({q['charge'] for q in rows if q['charge'] and q['state'] in (0,3)},key=abs)
        self.charge_values=[None]+charges+['uncertain']
        keys=['All charges']+[signed(z) for z in charges]+['Uncertain']
        uncertain=sum(q['state'] in (1,2,5) for q in rows)
        labels=(['All · %d masses'%len(rows)]+['%s · %d masses'%(signed(z),sum(q['charge']==z for q in rows)) for z in charges]+
                ['Uncertain · %d'%uncertain])
        order=list(range(1,len(charges)+1))+[0]+([len(keys)-1] if uncertain else [])  # (the charges, All, Uncertain)
        self.charge.SetItems(keys,labels,order)
        self.charge.SetSelection(keys.index(old) if old in keys and keys.index(old) in order else 0)
        self.family.SetItems(['All series']+[self.family_name(s,i) for i in range(len(s['families']))] if s else ['All series'])
        self.family.InvalidateBestSize();self.steps_changed()
        self.family.SetSelection(min(max(0,family),self.family.GetCount()-1))
    def charge_value(self):
        i=self.charge.GetSelection()
        return self.charge_values[i] if 0<=i<len(self.charge_values) else None
    def family_rows(self):
        s=self.current();i=self.family.GetSelection()-1
        return set(s['families'][i]['matches']) if s and 0<=i<len(s['families']) else None
    def filter_changed(self,e=None):
        self.refill();self.redraw();self._changed()
    def merge(self,rows,options,charge=None):
        """Q.merge, or no masses (and the message) when the native library cannot be loaded: the saved rows
        are shown all the same."""
        try:
            groups,mapping=Q.merge(rows,options.get('merge_da',.02),options.get('merge_ppm',5.),charge=charge)
            self.unavailable='';return groups,mapping
        except Q.Unavailable as ex:
            self.unavailable=str(ex);return [],[-1]*len(rows)
    def groups_for(self,s,charge=None,family=0):
        """Neutral masses of run s: included rows of that charge (None all, 'uncertain' none) and series
        (0 all, else the series number)."""
        if not s or charge=='uncertain':return []
        f=s['families'][family-1] if 0<family<=len(s['families']) else None
        matches=set(f['matches']) if f else None
        rows=[dict(q,include=bool(q.get('include')) and (matches is None or i in matches)) for i,q in enumerate(s['rows'])]
        return self.merge(rows,s['options'],charge)[0]
    def scoped_groups(self,include_family=True):
        return self.groups_for(self.current(),self.charge_value(),self.family.GetSelection() if include_family else 0)
    def mass_text(self,s,mass):
        if not mass or mass<=0:return '—'
        if s['source'].get('kind')!='direct':
            import deconv_tab as D
            p=(s['source'].get('result') or {}).get('params',{})
            return D._mass_fmt(mass,float(p.get('mass_step',1.) or 1.))
        return '%.4f'%mass
    def row_series(self,s,i):
        for k,f in enumerate(s['families']):
            if i in f['matches']:return k,f['matches'][i]
        return None,None
    def table_row(self,s,i):
        q=s['rows'][i];k,match=self.row_series(s,i)
        uncertain=q['state'] in (1,2,5)
        cells=['',self.mass_text(s,q['mass']),signed(q['charge']) if q['state'] in (0,3) else '—',
               number(q.get('mz')) if q.get('mz') else '—',height_text(q['height']),
               q['isotopes'] if q.get('isotopes') else '—',
               self.family_name(s,k) if k is not None else ('uncertain' if q['state']==1 else ''),
               match['n'] if match else '—','%.3g'%q.get('area',0),'%.3g'%q.get('response',1.),Q.STATES.get(q['state'],'')]
        colours={2:charge_colour(q['charge'])} if q['state'] in (0,3) else {}
        if k is not None:colours[6]=self.family_colour(k)
        return dict(cells=cells,colours=colours,muted=uncertain,check=bool(q.get('include')),check_disabled=uncertain)
    def refill(self):
        self._filling=True
        try:
            s=self.current()
            self.runs.SetItems([self.run_label(i,r) for i,r in enumerate(self.series)] or ['No run yet'])
            self.runs.SetSelection(self.active if self.active is not None else 0)
            self.runs.SetMinSize(wx.Size(max(self.FromDIP(170),self.runs.GetBestSize()[0]),-1))
            self.shown=[];rows=[];self.groups=[];self.mapping=[]
            if s:
                self.groups,self.mapping=self.merge(s['rows'],s['options'])
                family=self.family_rows();charge=self.charge_value()
                for i,q in enumerate(s['rows']):
                    if charge=='uncertain' and q['state'] not in (1,2,5):continue
                    if isinstance(charge,int) and q['charge']!=charge:continue
                    if family is not None and i not in family:continue
                    self.shown.append(i);rows.append(self.table_row(s,i))
                n=len(self.scoped_groups())
                info=self.unavailable or '%d assignments · %d mass%s'%(len(self.shown),n,'' if n==1 else 'es')
            else:info=self.unavailable
            self.chain_table.set_rows(rows,info)
            if self.chain is not None and self.chain not in self.shown:self.chain=None
            if self.chain is not None:self.chain_table.select(self.shown.index(self.chain))
            self.header.Layout();self.draw_mass();self.draw_assignments();self.Layout()
        finally:self._filling=False
    def select_table(self,index):
        if not self._filling and 0<=index<len(self.shown):self.select_chain(self.shown[index],preserve_table=True)
    def select_chain(self,index,preserve_table=False):
        s=self.current()
        if s is None or not 0<=index<len(s['rows']):return
        self.chain=index;self._filling=True
        try:
            if index in self.shown:self.chain_table.select(self.shown.index(index),only=not preserve_table)
        finally:self._filling=False
        self.show_source();self.draw_mass();self.redraw()
    def draw_mass(self):
        """Neutral masses as sticks in the colour of their series (grey: in no series, orange: selected);
        labels on the tallest peaks that stand apart."""
        card=self.mass_card;groups=self.scoped_groups();self._plot_groups=groups;s=self.current()
        scope=(s['id'] if s else None,self.charge_value(),self.family.GetSelection(),tuple((g['mass'],g['height']) for g in groups))
        view=(card.ax.get_xlim(),card.ax.get_ylim()) if scope==getattr(self,'_mass_scope',None) and card.full else None
        self._mass_scope=scope;card.reset();ax=card.ax
        if not groups:
            ax.set_axis_off()
            ax.text(.5,.5,self.unavailable or ('No masses with these choices' if s else 'Find species to see the neutral masses'),
                    transform=ax.transAxes,ha='center',va='center',fontsize=9,color=C['muted'])
            card.full=None;card.set_legend([]);card.draw();return
        families=s['families'];selected=None
        if self.chain is not None:
            gid=self.mapping[self.chain] if self.chain<len(self.mapping) else -1
            selected=next((g for g in groups if self.chain in g['rows'] or gid>=0 and any(self.mapping[i]==gid for i in g['rows'])),None)
        colours=[]
        for g in groups:
            k=next((k for k,f in enumerate(families) if any(i in f['matches'] for i in g['rows'])),None)
            colours.append(self.family_colour(k) if k is not None else (GREY if families else C['accent']))
        x=np.array([g['mass'] for g in groups]);y=np.array([g['height'] for g in groups])
        ax.vlines(x,0,y,colors=colours,lw=1.)
        if selected is not None:ax.vlines([selected['mass']],0,[selected['height']],color=FOCUS,lw=1.6,zorder=5)
        span=x.max()-x.min();pad=max(1.,span*.05)
        card.full=(x.min()-pad,x.max()+pad);ax.set_xlim(card.full);ax.set_ylim(0,y.max()*1.18 if y.max()>0 else 1)
        if view:ax.set_xlim(view[0]);ax.set_ylim(view[1])
        U.style_axes(ax,'Mass (Da)','Intensity')
        self._label_masses(ax,x,y,selected)
        legend=[(self.family_colour(k),self.family_name(s,k)) for k in range(len(families))]
        if families:legend.append((GREY,'no series'))
        legend.append((FOCUS,'selected'))
        card.set_legend(legend);card.draw()
    def _label_masses(self,ax,x,y,selected):
        """Labels (the size of the spectrum's peak labels) on the tallest peaks, each at least a label width
        from the others and with no taller peak under it: no label covers another at any width."""
        lo,hi=ax.get_xlim();top=ax.get_ylim()[1]
        bbox=ax.get_window_extent();px=bbox.width/max(hi-lo,1e-12)
        if px<=0:return
        dpi=ax.figure.dpi;char=7*dpi/72.*.6;s=self.current();fmt='%.2f'
        if s and s['source'].get('kind')!='direct':
            step=float(((s['source'].get('result') or {}).get('params') or {}).get('mass_step',1.) or 1.)
            fmt='%.4f' if step<.05 else '%.2f' if step<.5 else '%.1f'  # (the decimals of the deconvolution)
        inside=np.where((x>=lo)&(x<=hi)&(y<=top))[0]
        order=sorted(inside,key=lambda i:-y[i])
        if selected is not None:order=[i for i in order if x[i]==selected['mass']]+[i for i in order if x[i]!=selected['mass']]
        placed=[];limit=max(1,int(bbox.width/(len(fmt%x.max())*char*1.6)))
        for i in order:
            if len(placed)>=limit:break
            w=len(fmt%x[i])*char/px  # (label width in Da)
            if x[i]-w*.55<lo or x[i]+w*.55>hi:continue  # (inside the plot)
            if any(abs(x[i]-x[j])<w*1.15 for j in placed):continue
            near=(x>x[i]-w*.6)&(x<x[i]+w*.6)
            if np.any(y[near]>y[i]*1.02):continue
            if y[i]<.06*top and i not in [k for k in placed]:continue
            sel=selected is not None and x[i]==selected['mass']
            t=ax.annotate(fmt%x[i],(x[i],y[i]),xytext=(0,3),textcoords='offset points',ha='center',va='bottom',
                          fontsize=7,color=FOCUS if sel else C['text'],fontweight='bold' if sel else 'normal')
            t._no_scale=True;placed.append(i)
    def pick_mass(self,x,y=None):
        groups=self._plot_groups
        if not groups:return False
        group=min(groups,key=lambda g:abs(g['mass']-x))
        limits=self.mass_card.ax.get_xlim();pixels=max(1,self.mass_card.ax.bbox.width)
        if abs(group['mass']-x)>max(.002,12*abs(limits[1]-limits[0])/pixels):return False
        self.select_chain(group['rows'][0]);return True

    # ------------------------------------------------------------------ charges
    def charge_items(self,index=None,col=None,mz=None):
        """The charge menu of a row (or of a peak of the spectrum without a row): +1 to +4 (or more), not
        assigned, other; the present one ticked."""
        s=self.current();sign=(s['options'].get('sign',1) if s else (self.input() or {}).get('sign',1))
        sign=-1 if sign<0 else 1
        q=s['rows'][index] if s and index is not None and 0<=index<len(s['rows']) else None
        now=q['charge'] if q and q['state'] in (0,3) else 0
        top=max(4,self.z1.GetValue(),abs(now))
        items=[(signed(sign*z),(lambda z=z:self.set_charge(index,sign*z,col,mz)),now==sign*z) for z in range(1,top+1)]
        items+=[('Not assigned',lambda:self.set_charge(index,0,col,mz),now==0),(None,None),
                ('Other…',lambda:self.other_charge(index,col,mz))]
        return items
    def spectrum_charge_items(self,col,mz):
        """The charge menu of the peak at mz of the spectrum (right click on it): of its row, or of the peak."""
        s=self.current()
        if s and (s['source']['kind']!='direct' or not self.source_matches(col,s)):return [('Polymer charge',None)]
        tol=max(.004,abs(mz)*10e-6)
        row=next((i for i,q in enumerate(s['rows']) for m in q.get('members',[]) if abs(m['mz']-mz)<=tol),None) if s else None
        return self.charge_items(row,col,mz)
    def popup_charge_menu(self,win,items,pos=None):
        m=U.build_menu(win,items)
        try:U.popup_menu(win,m,pos)
        finally:m.Destroy()
    def table_cell(self,n,key,pos):
        if key!='z' or not 0<=n<len(self.shown):return False
        i=self.shown[n];self.select_chain(i)
        wx.CallAfter(self.popup_charge_menu,self.chain_table.body,self.charge_items(i),pos);return True
    def open_charge_menu_for_row(self,i):
        self.select_chain(i);self.popup_charge_menu(self.chain_table.body,self.charge_items(i))
    def source_column(self,s):
        col=next((c for c in self.tab.cols if s and c['e']==s['source']['event']),None)
        return col if col is not None and col.get('spec') is not None and self.source_matches(col,s) else None
    def set_charge(self,index,z,col=None,mz=None):
        """The charge of row index (or of the peak at mz without a row): z, or 0 for not assigned."""
        s=self.current()
        if index is None:
            if not z:return None
            return self.assign_from_spectrum(col,mz,value=z,isotopes=True)
        if s is None or not 0<=index<len(s['rows']):return None
        q=s['rows'][index]
        if not z:
            q.update(charge=0,state=2,mass=0.,include=False,alternate=0)
            self.rows_changed(changed=[index]);self.rebuild_filters();self.refill();self.redraw();self._changed()
            return q
        if s['source']['kind']!='direct':return self.fail('Deconvoluted masses have no charge to set')
        col=self.source_column(s)
        peak=max(q.get('members') or [dict(mz=q.get('mz'),height=1)],key=lambda m:m['height'])['mz']
        if col is not None and peak:
            return self.assign_from_spectrum(col,peak,value=z,isotopes=True,run=s)
        self.chain=index;return self.manual_charge(value=z)
    def other_charge(self,index=None,col=None,mz=None):
        s=self.current()
        if index is not None and s:
            q=s['rows'][index];col=self.source_column(s)
            mz=max(q.get('members') or [dict(mz=q.get('mz'),height=1)],key=lambda m:m['height'])['mz']
        if col is None or mz is None:return self.fail('Show the spectrum of this run first')
        return self.assign_from_spectrum(col,mz,run=s if index is not None else None)
    def manual_charge(self,e=None,value=None):
        s=self.current()
        if s is None or self.chain is None:return
        q=s['rows'][self.chain]
        if s['source']['kind']!='direct':return self.fail('Deconvoluted masses have no charge to set')
        if value is None:return self.other_charge(self.chain)
        try:
            sign=1 if s['source'].get('sign',s['options'].get('sign',1))>=0 else -1
            value=Q.whole(value,'Choose a whole charge')
            if sign<0 and value>0:value=-value  # (negative ions: 2 is -2)
            if value*sign<=0:raise ValueError('The charge must have the sign of the spectrum')
            replacement,remaining=Q.resolve(q,value,s['options'])
            s['rows'][self.chain]=replacement;s['rows'].extend(remaining);s['chains']=s['rows']
            self.rows_changed(changed=[self.chain]);self._changed();self.rebuild_filters();self.refill();self.redraw();return True
        except (ValueError,RuntimeError,OSError) as ex:return self.fail(ex)
    def assign_from_spectrum(self,col,mz,value=None,isotopes=False,run=None):
        """A charge for the peak at mz of the spectrum col (with its isotope peaks when isotopes): into run (its
        settings), the present run when its settings and spectrum fit, or a new run."""
        if self._busy or callable(getattr(self.tab,'busy_averaging',None)) and self.tab.busy_averaging():
            return self.fail('Wait for the spectrum to load')
        message=Q.available()
        if message:return self.fail(message)
        try:
            if col is None or col.get('spec') is None:raise ValueError('Show the spectrum first')
            self.refresh_inputs()
            index=next((i for i,d in enumerate(self.inputs) if d.get('e',0)==col['e']),None)
            if index is None:raise ValueError('Show the spectrum first')
            self.input_index=index;self.source_changed()
            sign=self.inputs[index].get('sign',1)
            if value is None:
                carrier=(run['options'].get('carrier_label'),run['options'].get('carrier')) if run else None
                with PeakCharge(self,mz,sign,spectrum=col['spec'],carrier=carrier) as dlg:
                    while True:
                        if dlg.ShowModal()!=wx.ID_OK:return None
                        try:value,label,carrier,isotopes=dlg.values();break
                        except ValueError as ex:wx.MessageBox(str(ex),'Assign charge',wx.OK|wx.ICON_INFORMATION,dlg)
                if run is not None and abs(carrier-run['options'].get('carrier',carrier))>1e-9:run=None
                self._carrier_chosen=True
                self.options.update(carrier=carrier,carrier_label=label);self._method.update(carrier=carrier,carrier_label=label)
                self.set_carriers()
            if run is not None:options=dict(run['options'])
            else:_,options=self.analysis_options()
            options['sign']=sign
            source=capture_source(self.tab,dict(options,source_label=self.inputs[index]['label'],source_event=int(col['e'])))
            replacement=Q.pick(np.array(col['spec'],copy=True),mz,value,options,isotopes)
            s=run if run is not None else self.current()
            saved=next((q['spec'] for q in s['source']['columns'] if q['e']==col['e']),None) if s else None
            same=(s is not None and s['source']['kind']=='direct' and same_file(self.tab,s['source']['file']) and
                  s['source']['event']==source['event'] and np.array_equal(saved,col['spec']) and
                  all(s['options'].get(k,Q.DEFAULTS.get(k))==options.get(k,Q.DEFAULTS.get(k))
                      for k in ('carrier','mass_mode','profile','isotope_tolerance')))
            if not same:
                if len(self.series)>=100:raise ValueError('Keep at most 100 polymer analyses in one file')
                rows=[]
                s=dict(id=uuid.uuid4().hex,source=source,options=options,rows=rows,chains=rows,
                       families=[],fields=self.field_values(),weight=self.weight.GetSelection(),basis=self.basis.GetSelection())
                self.series.append(s)
            self.active=self.series.index(s)
            # The measured peaks go to the new assignment; peaks of an old one that it does not take stay.
            claimed=[m['mz'] for m in replacement['members']];rows=[];moved={}
            for old,q in enumerate(s['rows']):
                members=q.get('members',[])
                remaining=[m for m in members if not any(abs(m['mz']-x)<1e-7 for x in claimed)]
                if len(remaining)==len(members):moved[old]=len(rows);rows.append(q);continue
                for m in remaining:
                    rows.append(dict(mass=0.,mz=m['mz'],area=m['area'],height=m['height'],first_mz=m['mz'],last_mz=m['mz'],
                        spacing_error=0.,response=q.get('response',1.),charge=0,isotopes=1,state=2,alternate=0,
                        members=[dict(m,isotope=0)],include=False))
            rows.append(replacement);s['rows']=s['chains']=rows;self.chain=len(rows)-1
            self.rows_changed(moved,changed=[self.chain])
            if self.chain not in self.shown:self.charge.SetSelection(0);self.family.SetSelection(0)
            self._changed();self.refill();self.redraw()
            self.say('%d peak%s set to %s'%(replacement['isotopes'],'' if replacement['isotopes']==1 else 's',signed(value)))
            return replacement
        except (ValueError,RuntimeError,OSError) as ex:return self.fail(ex)

    # ------------------------------------------------------------------ rows
    def delete_assignments(self,indices=None):
        s=self.current()
        if s is None or self._busy:return 0
        if not isinstance(indices,(list,tuple)):indices=self.chain_table.selected()
        removed={self.shown[i] for i in indices if isinstance(i,int) and 0<=i<len(self.shown)}
        if not removed:return 0
        old=self.chain;kept=[i for i in range(len(s['rows'])) if i not in removed]
        s['rows']=[s['rows'][i] for i in kept];s['chains']=s['rows']
        self.chain=(min(old-sum(i<old for i in removed),len(s['rows'])-1) if old is not None and old not in removed and s['rows'] else None)
        self.rows_changed({o:n for n,o in enumerate(kept)})
        self.refill();self.redraw();self._changed()
        self.say('Deleted %d peak%s'%(len(removed),'' if len(removed)==1 else 's'))
        return len(removed)
    def check_row(self,n):
        if 0<=n<len(self.shown):self.chain=self.shown[n];self.toggle_selected()
    def toggle_table(self,index):self.check_row(index)
    def toggle_selected(self,e=None):
        s=self.current()
        if s is None or self.chain is None:return
        q=s['rows'][self.chain]
        if q['state'] in (1,2,5):return self.fail('Set a charge before using this peak')
        q['include']=not q.get('include',True);self.rows_changed();self._changed();self.refill();self.redraw()
    def edit_response(self,e=None):
        s=self.current()
        if s is None or self.chain is None:return
        q=s['rows'][self.chain]
        with wx.TextEntryDialog(self,'Response factor','Response factor','%g'%q.get('response',1.)) as dlg:
            if dlg.ShowModal()!=wx.ID_OK:return
            try:
                value=float(dlg.GetValue().replace(',','.'))
                if not math.isfinite(value) or value<=0:raise ValueError('The response factor must be above zero')
                q['response']=value;self.rows_changed();self._changed();self.refill()
            except ValueError as ex:self.fail(ex)
    def rows_changed(self,moved=None,changed=()):
        """Rows of the present run edited. moved: {old index: new index} after rows were removed or rebuilt
        (rows left out are gone); changed: indices of rows whose mass changed. Series matches follow the rows
        (changed rows are matched again), a calculated mass distribution is calculated again."""
        s=self.current()
        if not s:return
        for f in s['families']:
            matches=unpack_matches(f['matches'])
            if moved is not None:matches={moved[i]:m for i,m in matches.items() if i in moved}
            for i in changed:
                matches.pop(i,None)
                q=s['rows'][i] if 0<=i<len(s['rows']) else None
                if q is None or not q.get('mass') or q['mass']<=0 or q['state'] not in (0,3,4):continue
                try:r=N.assign([q['mass']],f['repeat'],f['ends'],f['tolerance'])
                except (ValueError,RuntimeError,OSError):continue
                if r['matched'][0]:matches[i]=dict(n=int(r['dp'][0]),delta=float(r['delta'][0]),ppm=float(r['ppm'][0]))
            f['matches']=matches
        if s.get('distribution'):self.recalculate_distribution(s,s['distribution'])
        self.rebuild_filters();self.draw_assignments();self.draw_distribution()

    # ------------------------------------------------------------------ series
    def find_repeats(self,e=None):return self._candidates(False)
    def find_offsets(self,e=None):return self._candidates(True)
    def _number(self,control,what):
        try:v=float(self.field_text(control).replace(',','.'))
        except ValueError:raise ValueError('Enter the %s'%what)
        if not math.isfinite(v):raise ValueError('Enter the %s'%what)
        return v
    def _candidates(self,offset):
        try:
            groups=self.scoped_groups()
            if not groups:raise ValueError('Find species first')
            result=N.candidates([g['mass'] for g in groups],repeat=self._number(self.repeat,'repeat unit') if offset else None,
                                lo=self._number(self.lo,'search range'),hi=self._number(self.hi,'search range'),
                                tolerance=self._number(self.tol,'tolerance'))
            s=self.current();found=dict(kind='offset' if offset else 'repeat',values=result)
            if offset:
                s['offsets']=found
                if (s.get('candidates') or {}).get('kind')=='offset':s.pop('candidates')
            else:
                if (s.get('candidates') or {}).get('kind')=='offset' and 'offsets' not in s:s['offsets']=s['candidates']
                s['candidates']=found
            self.draw_candidates();self._changed()
            self.say('%d %s'%(len(result),'end group masses' if offset else 'repeat units'));return result
        except (ValueError,RuntimeError,OSError) as ex:return self.fail(ex)
    def suggestions(self,s,kind):
        """The suggestions of run s found by Find: kind 'repeat' (repeat units) or 'offset' (end groups)."""
        if not s:return []
        c=s.get('offsets') if kind=='offset' else None
        if c is None:
            c=s.get('candidates')
            if not c or c.get('kind','repeat')!=kind:return []
        return c.get('values') or []
    def draw_candidates(self):
        s=self.current()
        self.repeat_list.set_rows([dict(cells=['%.4f'%v['mass'],v['support'],v['single']]) for v in self.suggestions(s,'repeat')])
        self.ends_list.set_rows([dict(cells=['%.4f'%v['mass'],v['support']]) for v in self.suggestions(s,'offset')])
        for lst,kind,field in ((self.repeat_list,'repeat',self.repeat),(self.ends_list,'offset',self.ends)):
            try:v=float(self.field_text(field))  # (the suggestion in the field is marked)
            except ValueError:continue
            k=next((k for k,c in enumerate(self.suggestions(s,kind)) if abs(c['mass']-v)<1e-9),None)
            if k is not None:lst.select(k)
    def use_candidate(self,index,kind='repeat'):
        """A suggestion clicked: its mass in the field."""
        values=self.suggestions(self.current(),kind)
        if 0<=index<len(values):
            control=self.ends if kind=='offset' else self.repeat
            self.set_field(control,repr(float(values[index]['mass'])));self.fields_changed()
            (self.ends_list if kind=='offset' else self.repeat_list).select(index,show=False)
    def assign_series(self,e=None):
        try:
            # A shown series must not keep a revised series from taking rows it did not match before; excluded
            # rows and the charge shown still count.
            s=self.current();groups=self.scoped_groups(include_family=False)
            if not s or not groups:raise ValueError('Find species first')
            repeat,ends,tol=self._number(self.repeat,'repeat unit'),self._number(self.ends,'end groups'),self._number(self.tol,'tolerance')
            result=N.assign([g['mass'] for g in groups],repeat,ends,tol);matches={}
            for i,g in enumerate(groups):
                if result['matched'][i]:
                    for row in g['rows']:matches[row]=dict(n=int(result['dp'][i]),delta=float(result['delta'][i]),ppm=float(result['ppm'][i]))
            if not matches:raise ValueError('No masses fit this series')
            family=dict(repeat=repeat,ends=ends,tolerance=tol,matches=matches)
            old=next((i for i,f in enumerate(s['families']) if f['repeat']==repeat and f['ends']==ends),None)
            if old is None:s['families'].append(family);old=len(s['families'])-1
            else:s['families'][old]=family
            self.rebuild_filters();self.family.SetSelection(old+1);self.set_views(table=0);self.refill();self.redraw();self._changed()
            self.say('%s: %d masses'%(self.family_name(s,old),len({g['mass'] for g in groups if any(r in matches for r in g['rows'])})))
            return family
        except (ValueError,RuntimeError,OSError) as ex:return self.fail(ex)
    def remove_series(self):
        s=self.current();i=self.family.GetSelection()-1
        if not s or not 0<=i<len(s['families']):return
        del s['families'][i]
        d=s.get('distribution')
        if d and d.get('family',0)>i:
            if d['family']==i+1:s['distribution']=None
            else:d['family']-=1
        self.family.SetSelection(0);self.rebuild_filters();self.refill();self.redraw();self.draw_distribution();self._changed()
    def draw_assignments(self):
        s=self.current();index=self.family.GetSelection()-1;self._assigned_rows=[];rows=[]
        families=range(len(s['families'])) if s and index<0 else [index] if s and index<len(s['families']) else []
        for k in families:
            f=s['families'][k]
            for g in self.groups_for(s,self.charge_value(),k+1):
                i=next((i for i in g['rows'] if i in f['matches']),None)
                if i is None:continue
                m=f['matches'][i]
                rows.append(dict(cells=[self.family_name(s,k),m['n'],self.mass_text(s,g['mass']),'%+.4f'%m['delta']],
                                 colours={0:self.family_colour(k)}));self._assigned_rows.append(i)
        n=len(s['families']) if s else 0
        self.assign_table.set_rows(rows,'%d series'%n if n else '')
        self.draw_candidates()
    def select_assignment(self,index):
        if not self._filling and 0<=index<len(self._assigned_rows):self.select_chain(self._assigned_rows[index])

    # ------------------------------------------------------------------ mass distribution
    def invalidate_distribution(self,e=None):
        """The signal or count choice changed: a calculated distribution is calculated again with it."""
        s=self.current()
        if s:
            s['weight']=self.weight.GetSelection();s['basis']=self.basis.GetSelection()
            d=s.get('distribution')
            if d:self.recalculate_distribution(s,dict(d,signal=s['weight'],basis=s['basis']))
        self.draw_distribution()
        if e:self._changed();e.Skip()
    def recalculate_distribution(self,s,d):
        """The distribution of run s for the choices of d (charge, series, signal, count), or none."""
        try:
            groups=self.groups_for(s,d.get('charge'),d.get('family',0))
            if not groups:raise ValueError('No masses left for the mass distribution')
            weights=[g['area' if d.get('signal',0)==0 else 'height'] for g in groups]
            summary,nf,wf=N.distribution([g['mass'] for g in groups],weights,basis=d.get('basis',0))
            s['distribution']=dict(mn=summary.mn,mw=summary.mw,dispersity=summary.dispersity,used=summary.used,
                groups=copy.deepcopy(groups),number=nf.copy(),weight=wf.copy(),charge=d.get('charge'),
                family=d.get('family',0),signal=d.get('signal',0),basis=d.get('basis',0))
            return s['distribution']
        except (ValueError,RuntimeError,OSError) as ex:
            s['distribution']=None;self.say(str(ex));return None
    def calculate_distribution(self,e=None):
        try:
            s=self.current()
            if not s or not self.scoped_groups():raise ValueError('Find species first' if not s else 'No masses with these choices')
            d=self.recalculate_distribution(s,dict(charge=self.charge_value(),family=max(0,self.family.GetSelection()),
                                                   signal=self.weight.GetSelection(),basis=self.basis.GetSelection()))
            if d is None:return None
            self.set_views(plot=1);self._changed()
            from types import SimpleNamespace
            return SimpleNamespace(mn=d['mn'],mw=d['mw'],dispersity=d['dispersity'],used=d['used'])
        except (ValueError,RuntimeError,OSError) as ex:return self.fail(ex)
    def distribution_scope(self,s,d,masses=True):
        k=d.get('family',0)
        what=self.family_name(s,k-1) if 0<k<=len(s['families']) else 'All series'
        if d.get('charge') not in (None,'uncertain'):what+=' · %s'%signed(d['charge'])
        return '%s · %d masses'%(what,len(d['groups'])) if masses else what
    def draw_distribution(self):
        card=self.dist_card;card.reset();card.set_legend([]);s=self.current();d=s.get('distribution') if s else None;ax=card.ax
        if d:
            x=np.array([g['mass'] for g in d['groups']]);k=d.get('family',0)
            colour=self.family_colour(k-1) if 0<k<=len(s['families']) else C['accent']
            step=np.median(np.diff(np.sort(x))) if len(x)>1 else 1.
            w=max(step*.36,1e-3)
            ax.bar(x-w/2,d['number'],width=w,color=colour,label='Number fraction')
            ax.bar(x+w/2,d['weight'],width=w,color='#E39A3B',label='Mass fraction')
            card.set_legend([(colour,'Number fraction'),('#E39A3B','Mass fraction')])
            card.full=(x.min()-step,x.max()+step);ax.set_xlim(card.full)
            ax.set_ylim(0,max(np.max(d['number']),np.max(d['weight']))*1.1)
            U.style_axes(ax,'Mass (Da)','Fraction')
            self.summary.set_values([('Mn','%.1f Da'%d['mn']),('Mw','%.1f Da'%d['mw']),('Đ','%.4f'%d['dispersity']),
                                     ('Masses','%d'%d['used'])],self.distribution_scope(s,d,masses=False))  # (as in the report)
        else:
            ax.set_axis_off()
            ax.text(.5,.5,'Find species first' if not s else 'Calculate to see the mass distribution',
                    transform=ax.transAxes,ha='center',va='center',fontsize=9,color=C['muted'])
            card.full=None;card.set_legend([])
            self.summary.set_values([],'',note='')
        card.draw();self.steps_changed()

    # ------------------------------------------------------------------ the mass spectrum
    def overlay(self,col,ax):
        card=col['spec_card']
        if not self.marks_for(col):
            if col.get('_polymer_key'):col['_polymer_key']=False;card.set_legend([])
            return
        s=self.current()
        if s['source']['kind']!='direct':
            # Deconvoluted masses: the charge states of the selected mass only.
            if self.chain is not None and 0<=self.chain<len(s['rows']):
                from deconv_tab import _charges_of
                for z,mz in _charges_of(s['rows'][self.chain]['mass'],s['source']['result']):
                    mark=ax.axvline(mz,color=FOCUS,lw=1.,ls='--',alpha=.8);mark._no_scale=True
            col['_polymer_key']=True;card.set_legend([(FOCUS,'charge states of the selected mass')])
            return
        self._thin_labels(card,ax)
        for line in ax.lines:
            if type(line).__name__=='EnvelopeLine':line.set_color(C['faint']);line.set_alpha(.45);break
        for a in ax.collections:
            if type(a).__name__=='SpectrumSticks':a.set_color(C['faint']);a.set_alpha(.45);break
        selected_group=self.mapping[self.chain] if self.chain is not None and self.chain<len(self.mapping) else -1
        selected=s['rows'][self.chain] if self.chain is not None and self.chain<len(s['rows']) else None
        single=bool(selected and selected['state']==3 and selected['isotopes']==1)
        buckets={};present=set()
        for i in self.shown:
            q=s['rows'][i]
            focus=i==self.chain or not single and selected_group>=0 and self.mapping[i]==selected_group
            if focus:key='selected'
            elif q['state'] in (1,2,5):key='none'
            elif not q.get('include'):continue
            else:key=q['charge'];present.add(abs(q['charge']))
            buckets.setdefault(key,[]).extend((m['mz'],m['height']) for m in q.get('members',[]))
        for key,points in sorted(buckets.items(),key=lambda kv:{'none':0,'selected':2}.get(kv[0],1) if isinstance(kv[0],str) else 1):
            if not points:continue
            x,y=zip(*points)
            colour=FOCUS if key=='selected' else GREY if key=='none' else charge_colour(key)
            marks=ax.vlines(x,0,y,color=colour,lw=1.8 if key=='selected' else 1.,alpha=.95,zorder=6 if key=='selected' else 5)
            marks._no_scale=True;marks.set_gid('polymer-charge-'+str(key))
        sign='−' if s['options'].get('sign',1)<0 else '+'
        legend=[(charge_colour(z),'%s%d'%(sign,z)) for z in sorted(present)]+[(GREY,'unassigned'),(FOCUS,'selected')]
        col['_polymer_key']=True;card.set_legend(legend)
    def _thin_labels(self,card,ax):
        """The peak labels of the spectrum that would cover one another (a low spectrum tile) are left out."""
        try:
            renderer=card.canvas.get_renderer();kept=[]
            for t in sorted([t for t in ax.texts if t.get_visible()],key=lambda t:-(t.xy[1] if hasattr(t,'xy') else 0)):
                box=t.get_window_extent(renderer).expanded(1.04,1.0)
                if any(box.overlaps(b) for b in kept):t.set_visible(False)
                else:kept.append(box)
        except Exception:pass
    def report_data(self):
        """What the HRMS report shows of the present run: neutral masses with their colours, the charge states
        as listed, the mass distribution."""
        s=self.current();groups=self.scoped_groups();masses=[]
        for g in groups:
            k=next((k for k,f in enumerate(s['families']) if any(i in f['matches'] for i in g['rows'])),None)
            masses.append((g['mass'],g['height'],self.family_colour(k) if k is not None else (GREY if s['families'] else C['accent'])))
        rows=[]
        for i in self.shown:
            q=s['rows'][i]
            if q['state'] not in (0,3,4):continue
            k,_=self.row_series(s,i)
            rows.append([self.mass_text(s,q['mass']),signed(q['charge']) if q['state'] in (0,3) else '—',
                         number(q.get('mz')) if q.get('mz') else '—',height_text(q['height']),
                         q['isotopes'] or '—',self.family_name(s,k) if k is not None else ''])
        d=s.get('distribution');dist=None
        if d:
            k=d.get('family',0)
            dist=dict(mass=[g['mass'] for g in d['groups']],number=list(d['number']),weight=list(d['weight']),
                      mn=d['mn'],mw=d['mw'],dispersity=d['dispersity'],used=d['used'],scope=self.distribution_scope(s,d),
                      colour=self.family_colour(k-1) if 0<k<=len(s['families']) else C['accent'])
        legend=[(self.family_colour(k),self.family_name(s,k)) for k in range(len(s['families']))]
        return dict(masses=masses,rows=rows,distribution=dist,legend=legend)

    # ------------------------------------------------------------------ saving, undo, legacy
    def snapshot(self):
        shared=[c['spec'] for c in getattr(self.tab,'cols',[]) if isinstance(c.get('spec'),np.ndarray)]
        series=[pack_run(run,shared) for run in self.series]
        return dict(kind='polymer_workspace',series=series,active=self.active,chain=self.chain,visible=self.visible,
                    charge=self.charge.GetStringSelection(),family=self.family.GetSelection(),page=self.page(),
                    views=self.views(),method=self.method_values(),view=__import__('ms_project').views(self.tab))
    def undo_state(self):
        import undo
        if getattr(self,'_undo_token',None) is not self.outputs:
            self._undo_data=undo._copy(self.series);self._undo_token=self.outputs
        return dict(series=getattr(self,'_undo_data',[]),method=self.method_values(),active=self.active,
                    chain=self.chain,visible=self.visible,charge=self.charge.GetStringSelection(),
                    family=self.family.GetSelection(),page=self.page(),views=self.views(),kind='polymer_workspace')
    def restore(self,state):
        if state.get('kind')=='polymer_dock':
            self.reset();self.import_legacy(state)
            active=state.get('active')
            if isinstance(active,int) and 0<=active<len(self.series):self.select_series(active)
            old_chain=state.get('chain');legacy=state.get('series',[])
            if self.active is not None and self.active<len(legacy) and isinstance(old_chain,int) and 0<=old_chain<len(legacy[self.active]['chains']):
                mass=legacy[self.active]['chains'][old_chain]['mass']
                index=next((i for i,q in enumerate(self.current()['rows']) if q['mass']==mass),None)
                if index is not None:self.select_chain(index)
            self.visible=bool(state.get('visible',True));self.fit_layout(self.visible);self.tab._layout_rows()
            if state.get('view'):__import__('ms_project').set_views(self.tab,state['view'])
            return
        self._work_token+=1;self._set_busy(False)
        self.series=[unpack_run(run) for run in copy.deepcopy(state.get('series',[]))];self.active=state.get('active')
        rebind(self.tab,self.series)
        self.active=self.active if isinstance(self.active,int) and 0<=self.active<len(self.series) else 0 if self.series else None
        self.chain=state.get('chain');chain=self.chain
        self.visible=bool(state.get('visible',True));self._bindings.clear();self.set_method(state.get('method',{}))
        if self.active is not None:self.select_series(self.active)
        else:self.rebuild_filters();self.refill()
        if state.get('charge') in self.charge.GetStrings():self.charge.SetStringSelection(state['charge'])
        self.family.SetSelection(min(max(0,state.get('family',0)),self.family.GetCount()-1))
        views=state.get('views')
        if isinstance(views,(list,tuple)) and len(views)==2:self.set_views(*views)
        else:self.set_page(state.get('page',0))
        s=self.current();self.chain=chain if s and isinstance(chain,int) and 0<=chain<len(s['rows']) else None
        self.fit_layout(self.visible)
        self.refill();self.draw_distribution();self.redraw();self._watch();self.outputs=object();self.tab._layout_rows()
        if state.get('view'):__import__('ms_project').set_views(self.tab,state['view'])
        if self.visible:
            self.fit_layout(True)
            wx.CallAfter(lambda:self and self.visible and self.fit_layout(True))
    def import_legacy(self,state):
        start=len(self.series)
        for s in state.get('series',[]):
            rows=copy.deepcopy(s.get('rows') or s['chains'])
            options=s.get('options') or dict(Q.DEFAULTS)
            kind=s['source']['kind'];sign=next((d.get('sign',1) for d in self.inputs if d.get('e',0)==s['source']['event']),1)
            for row_index,q in enumerate(rows):
                grid=s.get('grid',[])
                if row_index<len(grid):
                    saved=grid[row_index];q.update(include=saved[0]=='1',mass=float(saved[1]),area=float(saved[2]),
                                                  height=float(saved[3]),response=float(saved[4]))
                q.update(state=4 if kind=='deconv' else 3,charge=0 if kind=='deconv' else sign,
                         alternate=0,spacing_error=0,response=q.get('response',1.),isotopes=q.get('isotopes',0),members=q.get('members',[]))
            matches={}
            for i,q in enumerate(rows):
                old=next((p for p in s['chains'] if p['mass']==q['mass']),None)
                if old:matches[i]={k:old[k] for k in ('n','delta','ppm')}
            self.series.append(dict(id=s['id'],source=s['source'],rows=rows,chains=rows,options=options,
                families=[dict(repeat=s['repeat'],ends=s['ends'],tolerance=s['tolerance'],matches=matches)],
                fields=[str(s['repeat']),str(s['ends']),str(s['tolerance']),'10','1000'],weight=s.get('weight',0),basis=s.get('basis',0),
                metadata=s.get('metadata',[])))
        rebind(self.tab,self.series[start:])
        if len(self.series)>start:self.active=start;self.select_series(start)
        self.visible=state.get('visible',True);self._changed()
    def add(self,record):
        if not same_file(self.tab,record['source']['file']):raise ValueError('These results belong to another file')
        if len(self.series)>=100:raise ValueError('Keep at most 100 polymer analyses in one file')
        self.import_legacy(dict(series=[copy.deepcopy(record)],visible=True));self.fit_layout(True);self.tab._layout_rows()
    def show_results(self):
        self.visible=True;self._watch();self.refresh_inputs();self.fit_layout(True);self.tab._layout_rows();self.refill();self.redraw()
    def hide_results(self,e=None):
        self.visible=False;self.Hide();self.fit_layout(False);self._changed();self.tab._layout_rows();self.redraw()
    def remove_analysis(self,e=None):
        if self._busy or self.active is None:return
        del self.series[self.active]
        self.active=min(self.active,len(self.series)-1) if self.series else None;self.chain=None
        if self.active is not None:self.select_series(self.active)
        else:self.rebuild_filters();self.refill();self.draw_distribution();self.redraw()
        self._changed()
    def reset(self):
        self._work_token+=1;self._set_busy(False);self.series=[];self.active=self.chain=None;self.visible=False
        self._bindings.clear();self.outputs=None;self.Hide();self.fit_layout(False);self.rebuild_filters();self.refill()
        if self in getattr(self.tab,'_polymer_frames',[]):self.tab._polymer_frames.remove(self)
    def export(self,e=None,path=None):
        s=self.current()
        if not s:return None
        if path is None:path=self.tab.frame.ask_save('Export polymer results','polymer_results.xlsx','Excel workbook (*.xlsx)|*.xlsx')
        if not path:return None
        if not path.lower().endswith('.xlsx'):path+='.xlsx'
        headers=['Analysis','Included','Neutral mass (Da)','Charge','m/z','Area / signal','Height','Response factor',
                 'Isotope peaks','Assignment','Series','Units (n)','First isotope m/z','Last isotope m/z','Spectrum']
        rows=[];meta=[('Tool','Polymer analysis'),('Mass equation','M = abs(charge) * (m/z - carrier per charge)'),
                      ('Precision','Full calculation precision'),('Charge merge','One neutral mass entry; corrected signals summed')]
        for ai,run in enumerate(self.series,1):
            meta.append(('Analysis %d input'%ai,'Measured spectrum' if run['source']['kind']=='direct' else 'Deconvoluted neutral masses'))
            meta.extend(('Analysis %d %s'%(ai,k),v) for k,v in run.get('metadata',[]))
            for i,q in enumerate(run['rows']):
                k,match=self.row_series(run,i)
                rows.append([ai,bool(q.get('include')),q['mass'] or '',q['charge'] or '',q.get('mz') or '',q['area'],q['height'],
                    q.get('response',1),q['isotopes'],Q.STATES[q['state']],self.family_name(run,k) if k is not None else '',
                    match['n'] if match else '',q.get('first_mz',''),q.get('last_mz',''),run['source']['label']])
            meta.extend(('Analysis %d %s'%(ai,k),v) for k,v in run['options'].items())
            for i,q in enumerate(run['rows'],1):
                for j,m in enumerate(q.get('members',[]),1):
                    meta.append(('Analysis %d assignment %d isotope %d'%(ai,i,j),'m/z %.17g; height %.17g; area %.17g'%(m['mz'],m['height'],m['area'])))
            for fi,f in enumerate(run['families'],1):
                meta.extend([('Analysis %d series %d name'%(ai,fi),self.family_name(run,fi-1)),
                             ('Analysis %d series %d repeat (Da)'%(ai,fi),f['repeat']),('Analysis %d series %d ends (Da)'%(ai,fi),f['ends'])])
                for index,m in f['matches'].items():meta.append(('Analysis %d series %d assignment %s'%(ai,fi,index),'n %d; error %.17g Da; %.17g ppm'%(m['n'],m['delta'],m['ppm'])))
            d=run.get('distribution')
            if d:
                for k in ('mn','mw','dispersity','used','charge','family','signal','basis'):meta.append(('Analysis %d distribution %s'%(ai,k),d[k]))
                for i,g in enumerate(d['groups']):meta.append(('Analysis %d neutral mass %d'%(ai,i+1),'mass %.17g; number fraction %.17g; mass fraction %.17g'%(g['mass'],d['number'][i],d['weight'][i])))
        try:__import__('lcms_compare_core').write_area_xlsx(path,headers,rows,meta)
        except (ValueError,OSError) as ex:return self.fail(ex)
        self.say('Polymer results exported');return path
