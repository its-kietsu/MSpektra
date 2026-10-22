"""One native HRMS polymer workspace. Numerical processing is in C++ libraries."""
import copy
import math
import threading
import uuid
import wx
import numpy as np
import unilcms as U
import ms_polymer as N
import ms_polymer_charge as Q
from unidec_theme import C,ui_font
from deconv_tab import TableCard
from ms_polymer_results import PolymerResults,capture_source,file_key,FOCUS
from ms_polymer_direct import POSITIVE,NEGATIVE,MASS_NAMES

COLOURS=(C['accent'],'#168579','#9D58B1','#9B6A20')
def number(value):
    return '%.4f'%value if isinstance(value,(float,np.floating)) and math.isfinite(value) else '—' if value is None else str(value)
def signed(z):return '%+d'%z if z else '—'
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
    return dock

def open_workspace(tab,page=0,use_deconvolution=False):
    try:
        dock=ensure_workspace(tab)
        dock.visible=True;dock.refresh_inputs();dock.book.SetSelection(page)
        dock._watch();tab._layout_rows();dock.Show();dock._keep_visible()
        if use_deconvolution:dock.use_deconvolution()
        return dock
    except (ValueError,RuntimeError,OSError) as ex:
        tab.status(str(ex));return None

class Settings(wx.Dialog):
    def __init__(self,parent,options):
        super().__init__(parent,title='Polymer settings',style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        self.options=dict(options);self.controls={}
        self.SetFont(ui_font(10));self.SetBackgroundColour(C['bg'])
        panel=wx.ScrolledWindow(self,style=wx.VSCROLL);panel.SetScrollRate(0,self.FromDIP(15))
        grid=wx.FlexGridSizer(cols=2,hgap=self.FromDIP(12),vgap=self.FromDIP(10));grid.AddGrowableCol(1)
        fields=[('mz_min','m/z from'),('mz_max','m/z to'),('threshold_percent','Minimum height (%)'),
                ('isotope_tolerance','Spacing tolerance (m/z)'),('min_isotopes','Minimum isotope peaks'),
                ('merge_da','Mass matching (Da)'),('merge_ppm','Mass matching (ppm)')]
        for key,label in fields:
            c=wx.TextCtrl(panel,value='%.15g'%options[key]);self.controls[key]=c
            grid.Add(wx.StaticText(panel,label=label),0,wx.ALIGN_CENTER_VERTICAL);grid.Add(c,1,wx.EXPAND)
        self.profile=wx.Choice(panel,choices=['Profile','Centroid']);self.profile.SetSelection(0 if options['profile'] else 1)
        self.mass=wx.Choice(panel,choices=MASS_NAMES);self.mass.SetSelection(options['mass_mode'])
        for label,c in [('Spectrum',self.profile),('Mass basis',self.mass)]:
            grid.Add(wx.StaticText(panel,label=label),0,wx.ALIGN_CENTER_VERTICAL);grid.Add(c,1,wx.EXPAND)
        box=wx.BoxSizer(wx.VERTICAL);box.Add(grid,1,wx.EXPAND|wx.ALL,self.FromDIP(14));panel.SetSizer(box)
        layout=wx.BoxSizer(wx.VERTICAL);layout.Add(panel,1,wx.EXPAND)
        layout.Add(self.CreateButtonSizer(wx.OK|wx.CANCEL),0,wx.EXPAND|wx.ALL,self.FromDIP(12));self.SetSizer(layout)
        self.SetSize(self.FromDIP(wx.Size(480,520)));self.SetMinSize(self.FromDIP(wx.Size(350,300)));self.CentreOnParent()
    def values(self):
        out=dict(self.options)
        for key,c in self.controls.items():
            v=float(c.GetValue())
            if not math.isfinite(v):raise ValueError('Enter finite settings')
            if key=='min_isotopes':
                if not v.is_integer():raise ValueError('Isotope count must be a whole number')
                v=int(v)
            out[key]=v
        out['profile']=self.profile.GetSelection()==0;out['mass_mode']=self.mass.GetSelection()
        if out['mz_min']<0 or out['mz_max']<=out['mz_min']:raise ValueError('Check the m/z range')
        if not 0<out['threshold_percent']<=100:raise ValueError('Minimum height must be above 0 and at most 100')
        if not 0<out['isotope_tolerance']<.25:raise ValueError('Spacing tolerance must be above 0 and below 0.25')
        if not 2<=out['min_isotopes']<=20:raise ValueError('Choose 2 to 20 isotope peaks')
        if not 0<=out['merge_da']<=1 or not 0<=out['merge_ppm']<=1000:raise ValueError('Check the mass matching tolerance')
        return out

class PeakCharge(wx.Dialog):
    def __init__(self,parent,mz,sign,spectrum=None):
        super().__init__(parent,title='Assign polymer charge')
        self.SetFont(ui_font(10));self.SetBackgroundColour(C['bg'])
        self.sign=sign;self.choices=(NEGATIVE if sign<0 else POSITIVE)+[('Custom mass',None)]
        grid=wx.FlexGridSizer(cols=2,hgap=self.FromDIP(12),vgap=self.FromDIP(10));grid.AddGrowableCol(1)
        self.charge=wx.SpinCtrl(self,min=1,max=100,initial=max(1,parent.z0.GetValue()))
        self.carrier=wx.Choice(self,choices=['Choose carrier…']+[q[0] for q in self.choices])
        label=parent.options.get('carrier_label') if parent._carrier_chosen else None
        self.carrier.SetSelection(next((i+1 for i,q in enumerate(self.choices) if q[0]==label),0))
        self.mass=wx.TextCtrl(self,value='%.15g'%parent.options['carrier'])
        for text,c in [('Charge (%s)'%('−' if sign<0 else '+'),self.charge),('Carrier',self.carrier),('Custom mass (Da)',self.mass)]:
            grid.Add(wx.StaticText(self,label=text),0,wx.ALIGN_CENTER_VERTICAL);grid.Add(c,1,wx.EXPAND)
        self.target=mz;self.spectrum=spectrum;self.options=dict(parent.options,sign=sign)
        self.isotopes=wx.CheckBox(self,label='Include isotope group');self.isotopes.SetValue(False)
        self.preview=wx.StaticText(self,label='This peak only')
        self.members=wx.ListCtrl(self,style=wx.LC_REPORT|wx.LC_SINGLE_SEL|wx.BORDER_SIMPLE,size=self.FromDIP(wx.Size(350,150)))
        self.members.InsertColumn(0,'m/z',width=self.FromDIP(150));self.members.InsertColumn(1,'Intensity',width=self.FromDIP(150))
        self.members.Hide();self.preview_result=None;self.preview_error=None
        self.carrier.Bind(wx.EVT_CHOICE,self.carrier_changed);self.carrier_changed()
        box=wx.BoxSizer(wx.VERTICAL)
        box.Add(wx.StaticText(self,label='m/z %.4f'%mz),0,wx.ALL,self.FromDIP(14))
        box.Add(grid,0,wx.EXPAND|wx.LEFT|wx.RIGHT,self.FromDIP(14))
        box.Add(self.isotopes,0,wx.ALL,self.FromDIP(14))
        box.Add(self.preview,0,wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM,self.FromDIP(14))
        box.Add(self.members,1,wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM,self.FromDIP(14))
        box.Add(self.CreateButtonSizer(wx.OK|wx.CANCEL),0,wx.EXPAND|wx.ALL,self.FromDIP(14))
        self.SetSizerAndFit(box);self.SetMinSize(self.FromDIP(wx.Size(390,-1)));self.fit_visible()
        self.isotopes.Bind(wx.EVT_CHECKBOX,self.update_preview)
        self.charge.Bind(wx.EVT_SPINCTRL,self.update_preview);self.mass.Bind(wx.EVT_TEXT,self.update_preview)
    def carrier_changed(self,e=None):
        self.mass.Enable(self.carrier.GetSelection()==len(self.choices))
        if hasattr(self,'preview'):self.update_preview()
    def update_preview(self,e=None):
        self.preview_result=None;self.preview_error=None;expanded=self.isotopes.GetValue()
        self.members.DeleteAllItems();self.members.Show(expanded)
        if not expanded:self.preview.SetLabel('This peak only')
        else:
            try:
                charge,label,carrier,_=self.values()
                if self.spectrum is None:raise ValueError('Choose a measured spectrum')
                options=dict(self.options,carrier=carrier)
                q=Q.pick(self.spectrum,self.target,charge,options,isotopes=True)
                self.preview_result=q
                self.preview.SetLabel('%d peaks · m/z %.4f to %.4f'%(q['isotopes'],q['first_mz'],q['last_mz'])
                    if q['isotopes']>1 else 'No matching group; this peak only')
                for m in q['members']:
                    i=self.members.InsertItem(self.members.GetItemCount(),number(m['mz']))
                    self.members.SetItem(i,1,number(m['height']))
            except (ValueError,RuntimeError,OSError) as ex:
                self.preview_error=str(ex);self.preview.SetLabel(self.preview_error)
        if self.GetSizer():self.fit_visible()
    def fit_visible(self):
        self.GetSizer().Layout();self.Fit()
        owner=wx.GetTopLevelParent(self.GetParent());rect=owner.GetScreenRect();w,h=self.GetSize()
        display=wx.Display.GetFromWindow(owner);area=wx.Display(max(0,display)).GetClientArea()
        x=rect.x+(rect.width-w)//2;y=rect.y+(rect.height-h)//2
        self.SetPosition((max(area.x,min(x,area.GetRight()+1-w)),max(area.y,min(y,area.GetBottom()+1-h))))
    def values(self):
        i=self.carrier.GetSelection()-1
        if i<0:raise ValueError('Choose a charge carrier')
        label,mass=self.choices[i]
        if mass is None:mass=float(self.mass.GetValue())
        if not math.isfinite(mass) or mass==0:raise ValueError('Enter a finite carrier mass')
        return self.sign*self.charge.GetValue(),label,mass,self.isotopes.GetValue()

class PolymerWorkspace(PolymerResults):
    show_empty=True
    def tile_min_height(self):
        return max(self.FromDIP(460),self.book.GetPosition().y+self.FromDIP(300))
    def __init__(self,tab):
        wx.Panel.__init__(self,tab.rows,style=wx.BORDER_NONE)
        self.tab=tab;self.SetBackgroundColour(C['bg']);self.SetFont(ui_font(9))
        self.series=[];self.active=None;self.chain=None;self.visible=True;self.outputs=None
        self._bindings={};self._restoring_source=False;self._filling=False;self._busy=False
        self.inputs=[];self.shown=[];self.groups=[];self.mapping=[];self._source_signature=None
        self.options=dict(Q.DEFAULTS);self._method=dict(getattr(tab,'_polymer_method',{}))
        self._profile_explicit='profile' in self._method
        self._carrier_chosen=bool(self._method.get('carrier_label'))
        self._work_token=0;self._exact_fields={}
        outer=wx.BoxSizer(wx.VERTICAL)
        heading=wx.WrapSizer(wx.HORIZONTAL)
        title=wx.StaticText(self,label='Polymer analysis');title.SetFont(ui_font(10,600))
        heading.Add(title,0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(14))
        self.runs=wx.Choice(self);self.runs.SetMinSize(self.FromDIP(wx.Size(300,-1)))
        heading.Add(self.runs,0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(10))
        self.export_btn=U.flat(self,'Export Excel',handler=self.export)
        for c in (U.flat(self,'Use deconvoluted masses',handler=self.use_deconvolution),self.export_btn,
                  U.flat(self,'Remove analysis',handler=self.remove_analysis),
                  U.flat(self,'Close',handler=self.hide_results)):
            heading.Add(c,0,wx.RIGHT,self.FromDIP(6))
        outer.Add(heading,0,wx.EXPAND|wx.ALL,self.FromDIP(7))
        self.source=wx.Choice(self);self.source.SetMinSize(self.FromDIP(wx.Size(245,-1)))
        self.z0=wx.SpinCtrl(self,min=1,max=100,initial=1,size=self.FromDIP(wx.Size(64,-1)))
        self.z1=wx.SpinCtrl(self,min=1,max=100,initial=3,size=self.FromDIP(wx.Size(64,-1)))
        self.carrier=wx.Choice(self);self.carrier.SetMinSize(self.FromDIP(wx.Size(150,-1)))
        self.find_btn=U.flat(self,'Find species','primary',handler=self.find_species)
        self.settings_btn=U.flat(self,'Settings',handler=self.settings)
        toolbar=wx.WrapSizer(wx.HORIZONTAL)
        for label,c in [('Spectrum',self.source),('Charge from',self.z0),('to',self.z1),('Carrier',self.carrier)]:
            toolbar.Add(wx.StaticText(self,label=label),0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(5))
            toolbar.Add(c,0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(9))
        for c in (self.find_btn,self.settings_btn):toolbar.Add(c,0,wx.RIGHT,self.FromDIP(6))
        outer.Add(toolbar,0,wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM,self.FromDIP(7))
        filters=wx.WrapSizer(wx.HORIZONTAL)
        self.charge=wx.Choice(self,choices=['All charges','Uncertain']);self.charge.SetSelection(0)
        self.family=wx.Choice(self,choices=['All species']);self.family.SetSelection(0)
        for label,c in [('Show charge',self.charge),('Polymer series',self.family)]:
            filters.Add(wx.StaticText(self,label=label),0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(5))
            filters.Add(c,0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(10))
        self.note=wx.StaticText(self,label='');self.note.SetForegroundColour(C['muted'])
        filters.Add(self.note,0,wx.ALIGN_CENTER_VERTICAL)
        outer.Add(filters,0,wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM,self.FromDIP(7))
        self.book=wx.Notebook(self)
        self.peak_page=wx.Panel(self.book);self.repeat_page=wx.Panel(self.book);self.dist_page=wx.Panel(self.book)
        for page,label in [(self.peak_page,'Charge states'),(self.repeat_page,'Repeat units'),(self.dist_page,'Mass distribution')]:
            page.SetBackgroundColour(C['bg']);self.book.AddPage(page,label)
        self.split=U.SplitBox(self.peak_page,wx.HORIZONTAL,min_size=100)
        self.mass_card=U.PlotCard(self.split,'Neutral mass spectrum','',mode='zoom');self.mass_card.on_pick=self.pick_mass
        self.mass_card.on_view=self.draw_mass
        self.chain_table=TableCard(self.split,'Charge assignments',[
            ('Use',40),('Units (n)',65),('Mass (Da)',110),('Charge',65),('m/z',110),
            ('Height',85),('Isotopes',65),('Assignment',100)],on_select=self.select_table,on_activate=self.toggle_table,multi=True)
        self.chain_table.on_delete=self.delete_assignments
        self.chain_table.list.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK,self.assignment_context)
        self.chain_table.extra_menu=lambda:[('Delete selected peaks',self.delete_assignments),
            ('Include or exclude',self.toggle_selected),('Assign charge…',self.manual_charge),
            ('Response factor…',self.edit_response),('Export Excel…',self.export)]
        self.split.set_panes([self.mass_card,self.chain_table],[1.6,1.5])
        pbox=wx.BoxSizer(wx.VERTICAL);pbox.Add(self.split,1,wx.EXPAND)
        buttons=wx.WrapSizer(wx.HORIZONTAL)
        for text,handler in [('Include or exclude',self.toggle_selected),('Assign charge…',self.manual_charge),('Delete selected',self.delete_assignments)]:
            buttons.Add(U.flat(self.peak_page,text,handler=handler),0,wx.RIGHT,self.FromDIP(6))
        pbox.Add(buttons,0,wx.ALL,self.FromDIP(6));self.peak_page.SetSizer(pbox)
        # Repeat finder, end group offsets and assignments share one page.
        self.repeat=wx.TextCtrl(self.repeat_page,value='',size=self.FromDIP(wx.Size(95,-1)))
        self.ends=wx.TextCtrl(self.repeat_page,value='',size=self.FromDIP(wx.Size(95,-1)))
        self.tol=wx.TextCtrl(self.repeat_page,value='.05',size=self.FromDIP(wx.Size(65,-1)))
        self.lo=wx.TextCtrl(self.repeat_page,value='10',size=self.FromDIP(wx.Size(65,-1)))
        self.hi=wx.TextCtrl(self.repeat_page,value='1000',size=self.FromDIP(wx.Size(65,-1)))
        repeatbar=wx.WrapSizer(wx.HORIZONTAL)
        for label,c in [('Repeat (Da)',self.repeat),('End groups (Da)',self.ends),('Tolerance (Da)',self.tol)]:
            repeatbar.Add(wx.StaticText(self.repeat_page,label=label),0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(5))
            repeatbar.Add(c,0,wx.RIGHT,self.FromDIP(8))
        searchbar=wx.WrapSizer(wx.HORIZONTAL)
        for label,c in [('Search from',self.lo),('to (Da)',self.hi)]:
            searchbar.Add(wx.StaticText(self.repeat_page,label=label),0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(5));searchbar.Add(c,0,wx.RIGHT,self.FromDIP(8))
        for text,handler in [('Find repeat units',self.find_repeats),('Find end groups',self.find_offsets)]:
            searchbar.Add(U.flat(self.repeat_page,text,handler=handler),0,wx.RIGHT,self.FromDIP(6))
        self.assign_btn=U.flat(self.repeat_page,'Assign series','primary',handler=self.assign_series)
        searchbar.Add(self.assign_btn,0,wx.RIGHT,self.FromDIP(6))
        self.repeat_split=U.SplitBox(self.repeat_page,wx.HORIZONTAL,min_size=100)
        self.candidates=TableCard(self.repeat_split,'Candidates',[('Mass (Da)',110),('Supported gaps',95),('Single gaps',85)],on_activate=self.use_candidate)
        self.assign_table=TableCard(self.repeat_split,'Chain assignments',[('Units (n)',75),('Mass (Da)',110),('Error (Da)',100)],on_select=self.select_assignment)
        self.candidates.extra_menu=lambda:[('Use candidate',lambda:self.use_candidate(self.candidates.selected()[0]) if self.candidates.selected() else None)]
        self.repeat_split.set_panes([self.candidates,self.assign_table],[1.,1.2])
        rbox=wx.BoxSizer(wx.VERTICAL)
        for row in (repeatbar,searchbar):rbox.Add(row,0,wx.EXPAND|wx.ALL,self.FromDIP(6))
        rbox.Add(self.repeat_split,1,wx.EXPAND);self.repeat_page.SetSizer(rbox)
        self.weight=wx.Choice(self.dist_page,choices=['Peak area','Peak height']);self.weight.SetSelection(0)
        self.basis=wx.Choice(self.dist_page,choices=['Signal per molecule','Signal per unit mass']);self.basis.SetSelection(0)
        distbar=wx.WrapSizer(wx.HORIZONTAL)
        for c in (self.weight,self.basis,U.flat(self.dist_page,'Calculate','primary',handler=self.calculate_distribution)):
            distbar.Add(c,0,wx.RIGHT,self.FromDIP(8))
        self.dist_split=U.SplitBox(self.dist_page,wx.HORIZONTAL,min_size=100)
        self.dist_card=U.PlotCard(self.dist_split,'Mass distribution','',mode='zoom')
        self.summary=TableCard(self.dist_split,'Results',[('Measure',100),('Value',125)])
        self.dist_split.set_panes([self.dist_card,self.summary],[3.,1.])
        dbox=wx.BoxSizer(wx.VERTICAL);dbox.Add(distbar,0,wx.EXPAND|wx.ALL,self.FromDIP(6));dbox.Add(self.dist_split,1,wx.EXPAND);self.dist_page.SetSizer(dbox)
        outer.Add(self.book,1,wx.EXPAND)
        self.detail=wx.StaticText(self,label='');outer.Add(self.detail,0,wx.EXPAND|wx.ALL,self.FromDIP(6))
        self.SetSizer(outer);self.SetMinSize(self.FromDIP(wx.Size(100,300)))
        self.source.Bind(wx.EVT_CHOICE,self.source_changed)
        self.carrier.Bind(wx.EVT_CHOICE,self.carrier_changed)
        self.runs.Bind(wx.EVT_CHOICE,lambda e:self.select_series(self.runs.GetSelection()))
        self.charge.Bind(wx.EVT_CHOICE,self.filter_changed);self.family.Bind(wx.EVT_CHOICE,self.filter_changed)
        self.weight.Bind(wx.EVT_CHOICE,self.invalidate_distribution);self.basis.Bind(wx.EVT_CHOICE,self.invalidate_distribution)
        for c in (self.repeat,self.ends,self.tol,self.lo,self.hi):c.Bind(wx.EVT_TEXT,self.fields_changed)
        for c in (self.repeat,self.ends):c.Bind(wx.EVT_KILL_FOCUS,self.format_field)
        self.book.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED,lambda e:(self.Layout(),e.Skip()))
        self.set_method(self._method);self.refresh_inputs();self._watch();self.Hide()

    def current(self):
        return self.series[self.active] if self.active is not None and 0<=self.active<len(self.series) else None
    def fail(self,ex):
        self.tab.status(str(ex));return None
    def _changed(self):
        super()._changed()
        import undo
        undo.schedule(self.tab)
    def method_values(self):
        return {**{k:v for k,v in self.options.items() if k not in ('mz_min','mz_max','source_label','source_event')},
                'carrier_label':self.options['carrier_label'] if self._carrier_chosen else '',
                'spectrum_type':('profile' if self.options['profile'] else 'centroid') if self._profile_explicit else 'auto',
                'charge_min':self.z0.GetValue(),'charge_max':self.z1.GetValue(),
                'repeat':self.field_text(self.repeat),'ends':self.field_text(self.ends),'tolerance':self.tol.GetValue(),
                'search_min':self.lo.GetValue(),'search_max':self.hi.GetValue(),
                'weight':self.weight.GetSelection(),'basis':self.basis.GetSelection()}
    def set_method(self,values):
        values=dict(values)
        if 'spectrum_type' in values:
            self._profile_explicit=values['spectrum_type']!='auto'
            if self._profile_explicit:values['profile']=values['spectrum_type']=='profile'
            else:values.pop('profile',None)
        if 'profile' in values:self._profile_explicit=True
        if 'carrier_label' in values:self._carrier_chosen=bool(values['carrier_label'])
        self.options.update({k:v for k,v in values.items() if k in Q.DEFAULTS})
        self.z0.SetValue(int(self.options['charge_min']));self.z1.SetValue(int(self.options['charge_max']))
        for key,c in [('repeat',self.repeat),('ends',self.ends),('tolerance',self.tol),('search_min',self.lo),('search_max',self.hi)]:
            if key in values:self.set_field(c,values[key])
        for key,c in [('weight',self.weight),('basis',self.basis)]:
            if values.get(key) in (0,1):c.SetSelection(values[key])
        self._method=self.method_values();self.tab._polymer_method=dict(self._method)
        if self.inputs:self.set_carriers()
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
    def refresh_inputs(self):
        selected=self.source.GetStringSelection()
        self.inputs=list(self.tab.deconv_inputs());self.source.SetItems([d['label'] for d in self.inputs])
        if self.inputs:self.source.SetSelection(next((i for i,d in enumerate(self.inputs) if d['label']==selected),0))
        self.source_changed()
    def source_changed(self,e=None):
        if self.inputs and self.source.GetSelection()>=0:
            d=self.inputs[self.source.GetSelection()];spec=d['spec']
            self.options.update(mz_min=float(spec[0,0]),mz_max=float(spec[-1,0]))
            if self.options['mz_max']==self.options['mz_min']:
                self.options['mz_min']=max(0.,self.options['mz_min']-.01);self.options['mz_max']+=.01
            profile=getattr(self.tab.data,'has_profile',True)
            if callable(profile):profile=profile()
            if not self._profile_explicit:self.options['profile']=bool(profile)
            self.set_carriers()
        self.find_btn.Enable(bool(self.inputs) and not self._busy)
    def set_carriers(self):
        sign=self.inputs[self.source.GetSelection()].get('sign',1)
        self.carriers=(NEGATIVE if sign<0 else POSITIVE)+[('Custom mass',None)]
        self.carrier.SetItems(['Choose carrier…']+[q[0] for q in self.carriers])
        label=self._method.get('carrier_label') if self._carrier_chosen else None
        self.carrier.SetSelection(next((i+1 for i,q in enumerate(self.carriers) if q[0]==label),0))
    def carrier_changed(self,e=None):
        i=self.carrier.GetSelection()-1
        if i<0:return
        label,value=self.carriers[i]
        if value is None:
            with wx.TextEntryDialog(self,'Carrier mass per charge (Da, signed)','Charge carrier',str(self.options['carrier'])) as dlg:
                if dlg.ShowModal()!=wx.ID_OK:self.carrier.SetSelection(0);return
                try:value=float(dlg.GetValue());assert math.isfinite(value) and value!=0
                except (ValueError,AssertionError):self.carrier.SetSelection(0);self.fail('Enter a finite carrier mass');return
        self._carrier_chosen=True
        self.options.update(carrier=value,carrier_label=label);self._method.update(carrier=value,carrier_label=label)
        self.tab._polymer_method=self.method_values()
    def settings(self,e=None):
        if not self.inputs:return self.fail('Choose a spectrum first')
        dlg=Settings(self,self.options)
        try:
            while dlg.ShowModal()==wx.ID_OK:
                try:self.options=dlg.values();self._profile_explicit=True;self._method=self.method_values();self.tab._polymer_method=dict(self._method);self._changed();break
                except ValueError as ex:wx.MessageBox(str(ex),'Polymer settings',wx.OK|wx.ICON_INFORMATION,dlg)
        finally:dlg.Destroy()
    def analysis_options(self):
        if not self.inputs:raise ValueError('Choose a spectrum first')
        if self.carrier.GetSelection()<1:raise ValueError('Choose a charge carrier')
        d=self.inputs[self.source.GetSelection()]
        options=dict(self.options,charge_min=self.z0.GetValue(),charge_max=self.z1.GetValue(),
                     source_label=d['label'],source_event=int(d.get('e',0)),sign=d.get('sign',1))
        if options['charge_max']<options['charge_min']:raise ValueError('Check the charge range')
        return d,options
    def find_species(self,e=None,synchronous=False):
        if self._busy:return None
        if len(self.series)>=100:return self.fail('Keep at most 100 polymer analyses in one file')
        if callable(getattr(self.tab,'busy_averaging',None)) and self.tab.busy_averaging():return self.fail('Wait for the spectrum to load')
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
        current=list(self.tab.deconv_inputs());label=self.source.GetStringSelection()
        if [d['label'] for d in current]!=[d['label'] for d in self.inputs]:
            self.refresh_inputs();raise ValueError('Spectrum changed. Choose it again')
        self.inputs=current
    def _set_busy(self,busy):
        self._busy=busy
        for c in (self.source,self.z0,self.z1,self.carrier,self.find_btn,self.settings_btn):c.Enable(not busy)
    def _finished(self,token,source,options,rows,error):
        if not self or self.IsBeingDeleted() or token!=self._work_token:return None
        self._set_busy(False)
        if error:return self.fail(error)
        if source['file']!=file_key(self.tab):return self.fail('The source file changed. Run again')
        col=next((c for c in self.tab.cols if c['e']==source['event']),None)
        saved=next((c['spec'] for c in source['columns'] if c['e']==source['event']),None)
        if col is None or col.get('spec') is None or not np.array_equal(col['spec'],saved):return self.fail('The spectrum changed. Run again')
        record=dict(id=uuid.uuid4().hex,source=source,options=options,rows=rows,chains=rows,
                    families=[],fields=self.field_values(),weight=self.weight.GetSelection(),basis=self.basis.GetSelection())
        # Re-running identical settings on the identical source replaces that run.
        old=next((i for i,s in enumerate(self.series) if s['options']==options and
                  s['source']['kind']=='direct' and s['source']['event']==source['event'] and
                  all(np.array_equal(a['spec'],b['spec']) for a,b in zip(s['source']['columns'],source['columns']))),None)
        if old is None:self.series.append(record);self.active=len(self.series)-1
        else:record['id']=self.series[old]['id'];self.series[old]=record;self.active=old
        self.chain=None
        self.visible=True;self._bindings.clear();self.charge.SetSelection(0);self.family.SetSelection(0)
        self._changed();self.select_series(self.active);self.tab._layout_rows()
        self.tab.status('%d charge assignments'%len(rows));return rows
    def use_deconvolution(self,e=None):
        if self._busy:return None
        if len(self.series)>=100:return self.fail('Keep at most 100 polymer analyses in one file')
        try:
            from ms_polymer_ui import source_rows
            rows,meta=source_rows(self.tab);source=capture_source(self.tab)
            options=dict(self.options);params=source['result'].get('params',{})
            if 'adduct_mass' in params:
                options['carrier']=float(params['adduct_mass'])
                choices=NEGATIVE if params.get('sign',1)<0 else POSITIVE
                options['carrier_label']=next((label for label,value in choices if abs(value-options['carrier'])<1e-8),'Custom mass')
            if params.get('z_range'):
                a,b=params['z_range'];options.update(charge_min=min(abs(int(a)),abs(int(b))),charge_max=max(abs(int(a)),abs(int(b))))
            options['sign']=params.get('sign',1)
            for q in rows:q.update(charge=0,state=4,alternate=0,isotopes=0,members=[],spacing_error=0,response=1.)
            record=dict(id=uuid.uuid4().hex,source=source,options=options,rows=rows,chains=rows,
                        families=[],fields=self.field_values(),metadata=meta,weight=0,basis=0)
            self.series.append(record);self.active=len(self.series)-1;self.visible=True;self.chain=None
            self._changed();self.select_series(self.active);self.tab._layout_rows();return rows
        except (ValueError,RuntimeError) as ex:return self.fail(ex)
    def select_series(self,index):
        if self._filling or not 0<=index<len(self.series):return
        self.active=index;self.chain=None;self._bindings.clear();s=self.current()
        self.options.update(s['options'])
        self._profile_explicit=True;self._carrier_chosen=True
        self.z0.SetValue(self.options.get('charge_min',1));self.z1.SetValue(self.options.get('charge_max',3))
        self._method.update({k:v for k,v in self.options.items() if k in Q.DEFAULTS})
        for c,v in zip((self.repeat,self.ends,self.tol,self.lo,self.hi),s.get('fields',[])):self.set_field(c,v)
        self.weight.SetSelection(s.get('weight',0));self.basis.SetSelection(s.get('basis',0))
        self.show_source();self.refresh_inputs()
        label=s['options'].get('source_label',s['source']['label'])
        self.source.SetSelection(next((i for i,d in enumerate(self.inputs) if d['label']==label),0) if self.inputs else wx.NOT_FOUND)
        self.options.update(s['options'])
        if self.inputs:self.set_carriers()
        self.rebuild_filters();self.refill();self.draw_distribution();self.redraw()
    def rebuild_filters(self):
        s=self.current();old=self.charge.GetStringSelection();family=self.family.GetSelection()
        charges=sorted({q['charge'] for q in s['rows'] if q['charge']}) if s else []
        self.charge_values=[None]+charges+['uncertain']
        labels=['All charges']+[signed(z) for z in charges]+['Uncertain']
        self.charge.SetItems(labels);self.charge.SetSelection(labels.index(old) if old in labels else 0)
        self.family.SetItems(['All species']+['Series %d · %.4f Da'%(i+1,f['repeat']) for i,f in enumerate(s['families'] if s else [])])
        self.family.SetSelection(min(max(0,family),self.family.GetCount()-1))
    def charge_value(self):return self.charge_values[self.charge.GetSelection()] if hasattr(self,'charge_values') else None
    def family_rows(self):
        s=self.current();i=self.family.GetSelection()-1
        return set(s['families'][i]['matches']) if s and 0<=i<len(s['families']) else None
    def filter_changed(self,e=None):
        self.refill();self.redraw();self.invalidate_distribution();self._changed()
    def scoped_groups(self,include_family=True):
        s=self.current()
        if not s or self.charge_value()=='uncertain':return []
        family=self.family_rows() if include_family else None
        rows=[dict(q,include=bool(q.get('include')) and (family is None or i in family)) for i,q in enumerate(s['rows'])]
        return Q.merge(rows,s['options'].get('merge_da',.02),s['options'].get('merge_ppm',5.),charge=self.charge_value())[0]
    def refill(self):
        self._filling=True
        try:
            s=self.current();self.runs.SetItems(['%d · %s'%(i+1,r['source']['label']) for i,r in enumerate(self.series)])
            if self.active is not None:self.runs.SetSelection(self.active)
            self.shown=[];display=[];self.groups=[];self.mapping=[]
            if s:
                self.groups,self.mapping=Q.merge(s['rows'],s['options'].get('merge_da',.02),s['options'].get('merge_ppm',5.))
                family=self.family_rows();charge=self.charge_value();active_family=self.family.GetSelection()-1
                matches=s['families'][active_family]['matches'] if active_family>=0 else {}
                for i,q in enumerate(s['rows']):
                    if charge=='uncertain' and q['state'] not in (1,2,5):continue
                    if isinstance(charge,int) and q['charge']!=charge:continue
                    if family is not None and i not in family:continue
                    self.shown.append(i);match=matches.get(i,matches.get(str(i),{}))
                    display.append(['Yes' if q.get('include') else 'No',match.get('n','—'),
                        number(q['mass']) if q['mass']>0 else '—',signed(q['charge']),number(q.get('mz')),
                        number(q['height']),q['isotopes'],Q.STATES.get(q['state'],'')])
                self.note.SetLabel('%d assignments · %d masses'%(len(self.shown),len(self.scoped_groups())))
            else:self.note.SetLabel('')
            self.chain_table.set_rows(display)
            if self.chain in self.shown:self.chain_table.select(self.shown.index(self.chain))
            self.export_btn.Enable(bool(s));self.draw_mass();self.draw_assignments();self.update_detail();self.Layout()
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
        self.show_source();self.draw_mass();self.redraw();self.update_detail();self._keep_visible()
    def update_detail(self):
        s=self.current()
        if s is None or self.chain is None:self.detail.SetLabel('');return
        q=s['rows'][self.chain]
        self.detail.SetLabel('%s%s · m/z %s · %s'%(number(q['mass'])+' Da · ' if q['mass'] else '',
            signed(q['charge']),number(q.get('mz')),Q.STATES.get(q['state'],'')))
    def draw_mass(self):
        groups=self.scoped_groups();self._plot_groups=groups
        scope=(self.current()['id'] if self.current() else None,self.charge_value(),self.family.GetSelection(),
               tuple((g['mass'],g['height']) for g in groups))
        view=(self.mass_card.ax.get_xlim(),self.mass_card.ax.get_ylim()) if scope==getattr(self,'_mass_scope',None) else None
        self._mass_scope=scope;self.mass_card.reset()
        if groups:
            ax=self.mass_card.ax;x=[q['mass'] for q in groups];y=[q['height'] for q in groups]
            ax.vlines(x,0,y,color=C['accent'],lw=.75)
            span=max(x)-min(x);pad=max(1.,span*.07)
            self.mass_card.full=(min(x)-pad,max(x)+pad);ax.set_xlim(self.mass_card.full)
            ax.set_ylim(0,max(y)*1.32 if max(y)>0 else 1)
            if view:ax.set_xlim(view[0]);ax.set_ylim(view[1])
            # Horizontal labels in measured pixel bounds, with a few vertical
            # lanes. Dense views omit collisions; zooming reconsiders labels.
            from matplotlib.font_manager import FontProperties
            lo,hi=ax.get_xlim();top=ax.get_ylim()[1];placed=[]
            renderer=self.mass_card.canvas.get_renderer();bounds=ax.bbox
            scale=self.mass_card.figure.dpi/72. if hasattr(self.mass_card,'figure') else ax.figure.dpi/72.
            candidates=sorted(groups,key=lambda g:(self.chain not in g['rows'],-g['height']))
            for g in candidates:
                if not lo<=g['mass']<=hi or g['height']>top:continue
                selected=self.chain in g['rows']
                text='%.4f'%g['mass'];font=FontProperties(size=8.5,weight='bold' if selected else 'normal')
                width,height,_=renderer.get_text_width_height_descent(text,font,ismath=False)
                px,py=ax.transData.transform((g['mass'],g['height']));half=width/2
                center=max(bounds.x0+half+3,min(bounds.x1-half-3,px))
                for dy in (7,20,33,46):
                    bottom=py+dy*scale;rect=(center-half-2,bottom-2,center+half+2,bottom+height+3)
                    if rect[3]>bounds.y1-2 or any(rect[0]<r[2] and rect[2]>r[0] and rect[1]<r[3] and rect[3]>r[1] for r in placed):continue
                    ax.annotate(text,(g['mass'],g['height']),xytext=((center-px)/scale,dy),textcoords='offset points',ha='center',va='bottom',
                                fontsize=8.5,rotation=0,color=FOCUS if selected else C['text'],fontweight='bold' if selected else 'normal')
                    placed.append(rect);break
            if self.chain is not None:
                group_id=self.mapping[self.chain] if self.chain<len(self.mapping) else -1
                g=next((g for g in groups if self.chain in g['rows'] or group_id>=0 and
                        any(self.mapping[i]==group_id for i in g['rows'])),None)
                if g:ax.vlines([g['mass']],0,[g['height']],color=FOCUS,lw=1.2)
        U.style_axes(self.mass_card.ax,'Mass (Da)','Intensity');self.mass_card.draw()
    def assignment_context(self,e):
        index=e.GetIndex()
        if 0<=index<len(self.shown) and index not in self.chain_table.selected():
            self.chain_table.select(index)
        e.Skip()
    def delete_assignments(self,indices=None):
        s=self.current()
        if s is None or self._busy:return 0
        if not isinstance(indices,(list,tuple)):indices=self.chain_table.selected()
        removed={self.shown[i] for i in indices if isinstance(i,int) and 0<=i<len(self.shown)}
        if not removed:return 0
        old=self.chain
        s['rows']=[q for i,q in enumerate(s['rows']) if i not in removed];s['chains']=s['rows']
        self.chain=(min(old-sum(i<old for i in removed),len(s['rows'])-1) if old is not None and s['rows'] else None)
        self.invalidate_tools();self.rebuild_filters();self.refill();self.redraw();self._changed()
        self.tab.status('Deleted %d peak%s'%(len(removed),'' if len(removed)==1 else 's'))
        return len(removed)
    def pick_mass(self,x,y=None):
        groups=getattr(self,'_plot_groups',[])
        if not groups:return False
        group=min(groups,key=lambda g:abs(g['mass']-x))
        limits=self.mass_card.ax.get_xlim();pixels=max(1,self.mass_card.ax.bbox.width)
        if abs(group['mass']-x)>max(.002,12*abs(limits[1]-limits[0])/pixels):return False
        self.select_chain(group['rows'][0]);return True
    def toggle_table(self,index):
        if 0<=index<len(self.shown):self.chain=self.shown[index];self.toggle_selected()
    def toggle_selected(self,e=None):
        s=self.current()
        if s is None or self.chain is None:return
        q=s['rows'][self.chain]
        if q['state'] in (1,2,5):return self.fail('Assign a charge before including this peak')
        q['include']=not q.get('include',True);self.invalidate_tools();self._changed();self.refill();self.redraw()
    def manual_charge(self,e=None,value=None):
        s=self.current()
        if s is None or self.chain is None:return
        q=s['rows'][self.chain]
        if s['source']['kind']!='direct':return self.fail('Choose measured spectrum input to assign a charge')
        if value is None:
            with wx.TextEntryDialog(self,'Signed charge state','Assign charge',str(q['charge'] or 1)) as dlg:
                if dlg.ShowModal()!=wx.ID_OK:return
                try:value=int(dlg.GetValue())
                except ValueError:return self.fail('Choose a whole charge state')
        try:
            sign=1 if s['source'].get('sign',s['options'].get('sign',1))>=0 else -1
            if value*sign<=0:raise ValueError('Charge sign must match the source spectrum')
            replacement,remaining=Q.resolve(q,value,s['options'])
            s['rows'][self.chain]=replacement;s['rows'].extend(remaining);s['chains']=s['rows']
            self.invalidate_tools();self._changed();self.rebuild_filters();self.refill();self.redraw();return True
        except (ValueError,RuntimeError,OSError) as ex:return self.fail(ex)
    def assign_from_spectrum(self,col,mz,value=None,isotopes=False):
        if self._busy or callable(getattr(self.tab,'busy_averaging',None)) and self.tab.busy_averaging():
            return self.fail('Wait for the spectrum to finish')
        try:
            self.refresh_inputs()
            index=next((i for i,d in enumerate(self.inputs) if d.get('e',0)==col['e']),None)
            if index is None or col.get('spec') is None:raise ValueError('Choose a measured spectrum')
            self.source.SetSelection(index);self.source_changed()
            sign=self.inputs[index].get('sign',1)
            if value is None:
                with PeakCharge(self,mz,sign,spectrum=col['spec']) as dlg:
                    while True:
                        if dlg.ShowModal()!=wx.ID_OK:return None
                        try:value,label,carrier,isotopes=dlg.values();break
                        except ValueError as ex:wx.MessageBox(str(ex),'Assign polymer charge',wx.OK|wx.ICON_INFORMATION,dlg)
                self._carrier_chosen=True
                self.options.update(carrier=carrier,carrier_label=label);self._method.update(carrier=carrier,carrier_label=label)
                self.set_carriers()
            _,options=self.analysis_options();options['sign']=sign
            source=capture_source(self.tab,options)
            replacement=Q.pick(np.array(col['spec'],copy=True),mz,value,options,isotopes)
            s=self.current()
            saved=next((q['spec'] for q in s['source']['columns'] if q['e']==col['e']),None) if s else None
            same=(s is not None and s['source']['kind']=='direct' and s['source']['file']==source['file'] and
                  s['source']['event']==source['event'] and np.array_equal(saved,col['spec']) and
                  all(s['options'].get(k,Q.DEFAULTS.get(k))==options.get(k,Q.DEFAULTS.get(k))
                      for k in ('carrier','mass_mode','profile','isotope_tolerance')))
            if not same:
                if len(self.series)>=100:raise ValueError('Keep at most 100 polymer analyses in one file')
                rows=[]
                s=dict(id=uuid.uuid4().hex,source=source,options=options,rows=rows,chains=rows,
                       families=[],fields=self.field_values(),weight=self.weight.GetSelection(),basis=self.basis.GetSelection())
                self.series.append(s);self.active=len(self.series)-1
            # Transfer ownership of measured members, leaving unmatched old members available for review.
            claimed=[m['mz'] for m in replacement['members']];rows=[]
            for q in s['rows']:
                members=q.get('members',[])
                remaining=[m for m in members if not any(abs(m['mz']-x)<1e-7 for x in claimed)]
                if len(remaining)==len(members):rows.append(q);continue
                for m in remaining:
                    rows.append(dict(mass=0.,mz=m['mz'],area=m['area'],height=m['height'],first_mz=m['mz'],last_mz=m['mz'],
                        spacing_error=0.,response=q.get('response',1.),charge=0,isotopes=1,state=2,alternate=0,
                        members=[dict(m,isotope=0)],include=False))
            rows.append(replacement);s['rows']=s['chains']=rows;self.chain=len(rows)-1
            self.invalidate_tools();self.charge.SetSelection(0);self.family.SetSelection(0)
            self._changed();self.refill();self.redraw();self._keep_visible()
            self.tab.status('Assigned %d peak%s to %s'%(replacement['isotopes'],'' if replacement['isotopes']==1 else 's',signed(value)))
            return replacement
        except (ValueError,RuntimeError,OSError) as ex:return self.fail(ex)
    def edit_response(self,e=None):
        s=self.current()
        if s is None or self.chain is None:return
        q=s['rows'][self.chain]
        with wx.TextEntryDialog(self,'Response factor','Response factor',str(q.get('response',1.))) as dlg:
            if dlg.ShowModal()!=wx.ID_OK:return
            try:
                value=float(dlg.GetValue())
                if not math.isfinite(value) or value<=0:raise ValueError('Response factor must be above zero')
                q['response']=value;self.invalidate_tools();self._changed();self.refill()
            except ValueError as ex:self.fail(ex)
    def invalidate_distribution(self,e=None):
        s=self.current()
        if s:s['distribution']=None;s['weight']=self.weight.GetSelection();s['basis']=self.basis.GetSelection()
        self.summary.set_rows([]);self.dist_card.reset();self.dist_card.draw()
        if e:self._changed();e.Skip()
    def invalidate_tools(self):
        s=self.current()
        if s:s['families']=[];s['candidates']=None
        self.rebuild_filters();self.candidates.set_rows([]);self.assign_table.set_rows([]);self.invalidate_distribution()
    def find_repeats(self,e=None):
        return self._candidates(False)
    def find_offsets(self,e=None):return self._candidates(True)
    def _candidates(self,offset):
        try:
            groups=self.scoped_groups()
            if not groups:raise ValueError('Include charge assignments first')
            result=N.candidates([g['mass'] for g in groups],repeat=float(self.field_text(self.repeat)) if offset else None,
                                lo=float(self.lo.GetValue()),hi=float(self.hi.GetValue()),tolerance=float(self.tol.GetValue()))
            self.current()['candidates']=dict(kind='offset' if offset else 'repeat',values=result)
            self.candidates.title='End group offsets' if offset else 'Repeat candidates'
            self.candidates.set_rows([[number(c['mass']),c['support'],c['single']] for c in result]);self._changed()
            self.tab.status('%d candidates'%len(result));return result
        except (ValueError,RuntimeError,OSError) as ex:return self.fail(ex)
    def use_candidate(self,index):
        s=self.current();c=s.get('candidates') if s else None
        if c and 0<=index<len(c['values']):
            control=self.ends if c['kind']=='offset' else self.repeat
            self.set_field(control,repr(float(c['values'][index]['mass'])));self.fields_changed()
    def assign_series(self,e=None):
        try:
            # A display filter must not prevent a revised series from acquiring
            # previously unmatched rows. Explicit exclusions and charge still apply.
            s=self.current();groups=self.scoped_groups(include_family=False)
            if not s or not groups:raise ValueError('Include charge assignments first')
            repeat,ends,tol=[float(self.field_text(c)) for c in (self.repeat,self.ends,self.tol)]
            result=N.assign([g['mass'] for g in groups],repeat,ends,tol);matches={}
            for i,g in enumerate(groups):
                if result['matched'][i]:
                    for row in g['rows']:matches[row]=dict(n=int(result['dp'][i]),delta=float(result['delta'][i]),ppm=float(result['ppm'][i]))
            if not matches:raise ValueError('No chains match these settings')
            family=dict(repeat=repeat,ends=ends,tolerance=tol,matches=matches)
            old=next((i for i,f in enumerate(s['families']) if f['repeat']==repeat and f['ends']==ends),None)
            if old is None:s['families'].append(family);old=len(s['families'])-1
            else:s['families'][old]=family
            self.rebuild_filters();self.family.SetSelection(old+1);self.refill();self.redraw();self._changed();return family
        except (ValueError,RuntimeError,OSError) as ex:return self.fail(ex)
    def draw_assignments(self):
        s=self.current();index=self.family.GetSelection()-1;self._assigned_rows=[];rows=[]
        if s and 0<=index<len(s['families']):
            f=s['families'][index]
            for g in self.scoped_groups():
                i=next((i for i in g['rows'] if i in f['matches'] or str(i) in f['matches']),None)
                if i is not None:
                    q=f['matches'].get(i,f['matches'].get(str(i)));rows.append([q['n'],number(g['mass']),number(q['delta'])]);self._assigned_rows.append(i)
        self.assign_table.set_rows(rows)
    def select_assignment(self,index):
        if not self._filling and 0<=index<len(self._assigned_rows):self.select_chain(self._assigned_rows[index])
    def calculate_distribution(self,e=None):
        try:
            groups=self.scoped_groups()
            if not groups:raise ValueError('Include charge assignments first')
            weights=[g['area' if self.weight.GetSelection()==0 else 'height'] for g in groups]
            summary,nf,wf=N.distribution([g['mass'] for g in groups],weights,basis=self.basis.GetSelection())
            self.current()['distribution']=dict(mn=summary.mn,mw=summary.mw,dispersity=summary.dispersity,
                used=summary.used,groups=copy.deepcopy(groups),number=nf.copy(),weight=wf.copy(),
                charge=self.charge_value(),family=self.family.GetSelection(),signal=self.weight.GetSelection(),basis=self.basis.GetSelection())
            self.draw_distribution();self._changed();return summary
        except (ValueError,RuntimeError,OSError) as ex:return self.fail(ex)
    def draw_distribution(self):
        self.dist_card.reset();s=self.current();d=s.get('distribution') if s else None
        if d:
            x=[g['mass'] for g in d['groups']]
            self.dist_card.ax.plot(x,d['number'],'o-',color=C['accent'],label='Number fraction')
            self.dist_card.ax.plot(x,d['weight'],'s-',color=FOCUS,label='Mass fraction')
            self.dist_card.ax.legend(fontsize=8);self.dist_card.full=(min(x)-.5,max(x)+.5)
            self.summary.set_rows([['Mₙ (Da)',number(d['mn'])],['M𝑤 (Da)',number(d['mw'])],
                ['Dispersity',number(d['dispersity'])],['Unique masses',d['used']]],'Signal based')
        else:self.summary.set_rows([])
        U.style_axes(self.dist_card.ax,'Mass (Da)','Fraction');self.dist_card.draw()
    def overlay(self,col,ax):
        if not self.marks_for(col):return
        s=self.current()
        if s['source']['kind']!='direct':
            from deconv_tab import _charges_of
            for i in self.shown:
                q=s['rows'][i]
                if not q.get('include'):continue
                for z,mz in _charges_of(q['mass'],s['source']['result']):
                    mark=ax.axvline(mz,color=FOCUS if i==self.chain else C['accent'],lw=1.,ls='--',alpha=.75)
                    mark._no_scale=True
            return
        for line in ax.lines:
            if type(line).__name__=='EnvelopeLine':line.set_color(C['faint']);line.set_alpha(.45);break
        for a in ax.collections:
            if type(a).__name__=='SpectrumSticks':a.set_color(C['faint']);a.set_alpha(.45);break
        selected_group=self.mapping[self.chain] if self.chain is not None and self.chain<len(self.mapping) else -1
        selected=s['rows'][self.chain] if self.chain is not None and self.chain<len(s['rows']) else None
        single=bool(selected and selected['state']==3 and selected['isotopes']==1)
        buckets={}
        for i in self.shown:
            q=s['rows'][i]
            if not q.get('include') and i!=self.chain:continue
            focus=i==self.chain or not single and selected_group>=0 and self.mapping[i]==selected_group
            key=('selected' if focus else q['charge'])
            buckets.setdefault(key,[]).extend((m['mz'],m['height']) for m in q.get('members',[]))
        for key,points in buckets.items():
            if not points:continue
            x,y=zip(*points);colour=FOCUS if key=='selected' else COLOURS[(abs(key)-1)%len(COLOURS)]
            marks=ax.vlines(x,0,y,color=colour,lw=1.8 if key=='selected' else 1.,alpha=.9,zorder=5)
            marks._no_scale=True;marks.set_gid('polymer-charge-'+str(key))
    def pick(self,col,mz,y=None):
        if not self.marks_for(col):return False
        s=self.current();targets=[]
        if s['source']['kind']=='direct':
            targets=[(m['mz'],i) for i in self.shown for m in s['rows'][i].get('members',[])]
        else:
            from deconv_tab import _charges_of
            targets=[(m,i) for i in self.shown for _,m in _charges_of(s['rows'][i]['mass'],s['source']['result'])]
        if not targets:return False
        peak=self.tab._nearest_peak(col,mz)
        if peak is None:return False
        target,index=min(targets,key=lambda q:abs(q[0]-peak))
        if abs(target-peak)>max(.004,abs(peak)*10e-6):return False
        self.select_chain(index);return True
    def snapshot(self):
        return dict(kind='polymer_workspace',series=self.series,active=self.active,chain=self.chain,visible=self.visible,
                    charge=self.charge.GetStringSelection(),family=self.family.GetSelection(),page=self.book.GetSelection(),
                    method=self.method_values(),view=__import__('ms_project').views(self.tab))
    def undo_state(self):
        import undo
        if getattr(self,'_undo_token',None) is not self.outputs:
            self._undo_data=undo._copy(self.series);self._undo_token=self.outputs
        return dict(series=getattr(self,'_undo_data',[]),method=self.method_values(),active=self.active,
                    chain=self.chain,visible=self.visible,charge=self.charge.GetStringSelection(),
                    family=self.family.GetSelection(),page=self.book.GetSelection(),kind='polymer_workspace')
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
            self.visible=bool(state.get('visible',True));self.tab._layout_rows()
            if state.get('view'):__import__('ms_project').set_views(self.tab,state['view'])
            return
        self._work_token+=1;self._set_busy(False)
        self.series=copy.deepcopy(state.get('series',[]));self.active=state.get('active')
        for run in self.series:
            for family in run['families']:family['matches']={int(i):m for i,m in family['matches'].items()}
        self.active=self.active if isinstance(self.active,int) and 0<=self.active<len(self.series) else 0 if self.series else None
        self.chain=state.get('chain');chain=self.chain
        self.visible=bool(state.get('visible',True));self._bindings.clear();self.set_method(state.get('method',{}))
        if self.active is not None:self.select_series(self.active)
        else:self.rebuild_filters();self.refill()
        if state.get('charge') in self.charge.GetStrings():self.charge.SetStringSelection(state['charge'])
        self.family.SetSelection(min(max(0,state.get('family',0)),self.family.GetCount()-1))
        self.book.SetSelection(min(2,max(0,state.get('page',0))))
        s=self.current();self.chain=chain if s and isinstance(chain,int) and 0<=chain<len(s['rows']) else None
        self.refill();self.draw_distribution();self.redraw();self._watch();self.outputs=object();self.tab._layout_rows()
        if state.get('view'):__import__('ms_project').set_views(self.tab,state['view'])
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
        if len(self.series)>start:self.active=start;self.select_series(start)
        self.visible=state.get('visible',True);self._changed()
    def add(self,record):
        if record['source']['file']!=file_key(self.tab):raise ValueError('These results belong to another file')
        if len(self.series)>=100:raise ValueError('Keep at most 100 polymer analyses in one file')
        self.import_legacy(dict(series=[copy.deepcopy(record)],visible=True));self.tab._layout_rows()
    def show_results(self):
        self.visible=True;self._watch();self.refresh_inputs();self.tab._layout_rows();self.refill();self.redraw()
    def remove_analysis(self,e=None):
        if self._busy or self.active is None:return
        del self.series[self.active]
        self.active=min(self.active,len(self.series)-1) if self.series else None;self.chain=None
        if self.active is not None:self.select_series(self.active)
        else:self.rebuild_filters();self.refill();self.draw_distribution();self.redraw()
        self._changed()
    def reset(self):
        self._work_token+=1;self._set_busy(False);self.series=[];self.active=self.chain=None;self.visible=False
        self._bindings.clear();self.outputs=None;self.Hide();self.rebuild_filters();self.refill()
        if self in getattr(self.tab,'_polymer_frames',[]):self.tab._polymer_frames.remove(self)
    def export(self,e=None,path=None):
        s=self.current()
        if not s:return None
        if path is None:path=self.tab.frame.ask_save('Export polymer results','polymer_results.xlsx','Excel workbook (*.xlsx)|*.xlsx')
        if not path:return None
        if not path.lower().endswith('.xlsx'):path+='.xlsx'
        headers=['Analysis','Included','Neutral mass (Da)','Charge','m/z','Area / signal','Height','Response factor',
                 'Isotope peaks','Assignment','First isotope m/z','Last isotope m/z','Spectrum']
        rows=[];meta=[('Tool','Polymer analysis'),('Mass equation','M = abs(charge) * (m/z - carrier per charge)'),
                      ('Precision','Full calculation precision'),('Charge merge','One neutral mass entry; corrected signals summed')]
        for ai,run in enumerate(self.series,1):
            meta.append(('Analysis %d input'%ai,'Measured spectrum' if run['source']['kind']=='direct' else 'Deconvoluted neutral masses'))
            meta.extend(('Analysis %d %s'%(ai,k),v) for k,v in run.get('metadata',[]))
            for q in run['rows']:
                rows.append([ai,bool(q.get('include')),q['mass'] or '',q['charge'] or '',q.get('mz') or '',q['area'],q['height'],
                    q.get('response',1),q['isotopes'],Q.STATES[q['state']],q.get('first_mz',''),q.get('last_mz',''),run['source']['label']])
            meta.extend(('Analysis %d %s'%(ai,k),v) for k,v in run['options'].items())
            for i,q in enumerate(run['rows'],1):
                for j,m in enumerate(q.get('members',[]),1):
                    meta.append(('Analysis %d assignment %d isotope %d'%(ai,i,j),'m/z %.17g; height %.17g; area %.17g'%(m['mz'],m['height'],m['area'])))
            for fi,f in enumerate(run['families'],1):
                meta.extend([('Analysis %d series %d repeat (Da)'%(ai,fi),f['repeat']),('Analysis %d series %d ends (Da)'%(ai,fi),f['ends'])])
                for index,m in f['matches'].items():meta.append(('Analysis %d series %d assignment %s'%(ai,fi,index),'n %d; error %.17g Da; %.17g ppm'%(m['n'],m['delta'],m['ppm'])))
            d=run.get('distribution')
            if d:
                for k in ('mn','mw','dispersity','used','charge','family','signal','basis'):meta.append(('Analysis %d distribution %s'%(ai,k),d[k]))
                for i,g in enumerate(d['groups']):meta.append(('Analysis %d neutral mass %d'%(ai,i+1),'mass %.17g; number fraction %.17g; mass fraction %.17g'%(g['mass'],d['number'][i],d['weight'][i])))
        try:__import__('lcms_compare_core').write_area_xlsx(path,headers,rows,meta)
        except (ValueError,OSError) as ex:return self.fail(ex)
        self.tab.status('Polymer results exported');return path
