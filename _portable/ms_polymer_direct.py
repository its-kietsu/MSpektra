"""UI for direct singly charged polymer input; numerical work is in mspolymer_direct.dll."""
import wx
import numpy as np
import ms_polymer as N
import unilcms as U
from unidec_theme import C, ui_font

POSITIVE=[('H+ (proton)',1.007276467),('Na+ (sodium)',22.989218),
          ('K+ (potassium)',38.963158),('NH4+ (ammonium)',18.033823)]
NEGATIVE=[('H+ loss',-1.007276467),('Cl- (chloride)',34.969401),
          ('Formate',44.998203),('Acetate',59.013853)]
MASS_NAMES=['First observed isotope','Most abundant isotope','Signal-weighted mean']

class DirectDialog(wx.Dialog):
    def __init__(self,tab,options=None):
        super().__init__(wx.GetTopLevelParent(tab),title='Singly charged polymer spectrum',
                         style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        self.tab=tab;self.options=dict(options or {});self.result=None;self._ready=False
        self.inputs=list(tab.deconv_inputs())
        if not self.inputs:raise ValueError('Select a scan or a time range in HRMS first.')
        self.SetFont(ui_font(10));self.SetBackgroundColour(C['bg'])
        panel=wx.ScrolledWindow(self,style=wx.VSCROLL)
        panel.SetScrollRate(0,self.FromDIP(15))
        self.source=wx.Choice(panel,choices=[d['label'] for d in self.inputs])
        selected=next((i for i,d in enumerate(self.inputs) if d['label']==self.options.get('source_label')),0)
        self.source.SetSelection(selected)
        self.carrier=wx.Choice(panel)
        self.custom=wx.TextCtrl(panel,value=str(self.options.get('custom','')))
        self.custom.Enable(False)
        self.low=wx.TextCtrl(panel);self.high=wx.TextCtrl(panel)
        self.threshold=wx.TextCtrl(panel,value=str(self.options.get('threshold_percent',5)))
        self.tolerance=wx.TextCtrl(panel,value=str(self.options.get('isotope_tolerance',.01)))
        self.format=wx.Choice(panel,choices=['Profile spectrum','Centroid peaks'])
        profile=getattr(getattr(tab,'data',None),'has_profile',True)
        if callable(profile):profile=profile()
        self.format.SetSelection(0 if self.options.get('profile',bool(profile)) else 1)
        self.group=wx.CheckBox(panel,label='Combine adjacent isotope peaks')
        self.group.SetValue(self.options.get('group_isotopes',True))
        self.mass=wx.Choice(panel,choices=MASS_NAMES)
        self.mass.SetSelection(int(self.options.get('mass_mode',0)))
        grid=wx.FlexGridSizer(cols=2,hgap=self.FromDIP(12),vgap=self.FromDIP(8));grid.AddGrowableCol(1)
        for label,control in [('Spectrum',self.source),('Charge carrier',self.carrier),
                ('Custom carrier mass (Da, signed)',self.custom),('m/z from',self.low),('m/z to',self.high),
                ('Minimum peak above background (%)',self.threshold),('Spectrum type',self.format),
                ('Isotope grouping',self.group),('Isotope spacing tolerance (Da)',self.tolerance),
                ('Mass reported for each group',self.mass)]:
            grid.Add(wx.StaticText(panel,label=label),0,wx.ALIGN_CENTER_VERTICAL)
            grid.Add(control,1,wx.EXPAND)
        self.mass_note=wx.StaticText(panel,label='First observed is not necessarily the monoisotopic peak.')
        self.preview=wx.ListCtrl(panel,style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
        for i,(label,width) in enumerate([('Mass (Da)',135),('m/z',125),('Isotope peaks',105),
                                         ('Area / signal',125),('Height',110)]):
            self.preview.InsertColumn(i,label,width=self.FromDIP(width))
        self.preview.SetMinSize(self.FromDIP(wx.Size(650,155)))
        preview_btn=U.flat(panel,'Preview peaks',handler=self.calculate)
        self.notice=wx.StaticText(panel,label='One charge per ion. Select one adduct series, then preview the peaks.')
        layout=wx.BoxSizer(wx.VERTICAL)
        for control,flags in [(self.notice,wx.EXPAND),(grid,wx.EXPAND),(self.mass_note,wx.EXPAND),
                               (preview_btn,0),(self.preview,wx.EXPAND)]:
            layout.Add(control,0,flags|wx.ALL,self.FromDIP(9))
        panel.SetSizer(layout)
        self.ok=wx.Button(self,wx.ID_OK,'Use these masses');self.ok.Enable(False)
        cancel=wx.Button(self,wx.ID_CANCEL,'Cancel')
        buttons=wx.BoxSizer(wx.HORIZONTAL);buttons.AddStretchSpacer()
        for b in (cancel,self.ok):buttons.Add(b,0,wx.ALL,self.FromDIP(8))
        outer=wx.BoxSizer(wx.VERTICAL);outer.Add(panel,1,wx.EXPAND);outer.Add(buttons,0,wx.EXPAND)
        self.SetSizer(outer);self.SetSize(self.FromDIP(wx.Size(820,740)))
        self.SetMinSize(self.FromDIP(wx.Size(650,440)))
        self.CentreOnParent()
        self.source.Bind(wx.EVT_CHOICE,self.change_source)
        self.carrier.Bind(wx.EVT_CHOICE,self.change_carrier)
        for c in (self.custom,self.low,self.high,self.threshold,self.tolerance):c.Bind(wx.EVT_TEXT,self.invalidate)
        for c in (self.format,self.mass):c.Bind(wx.EVT_CHOICE,self.invalidate)
        self.group.Bind(wx.EVT_CHECKBOX,self.invalidate)
        self.Bind(wx.EVT_BUTTON,self.accept,id=wx.ID_OK)
        self.change_source(initial=True);self._ready=True
        self.mass.Enable(self.group.GetValue());self.tolerance.Enable(self.group.GetValue())
        panel.FitInside()
        # Keep the dialog usable on a small display or high-DPI monitor.
        display=wx.Display.GetFromWindow(self)
        if display!=wx.NOT_FOUND:
            bounds=wx.Display(display).GetClientArea()
            self.SetSize(wx.Size(min(self.GetSize().width,bounds.width-20),min(self.GetSize().height,bounds.height-30)))
            self.CentreOnParent()

    def invalidate(self,e=None):
        if not self._ready:return
        self.result=None;self.ok.Enable(False);self.preview.DeleteAllItems()
        mode=self.mass.GetSelection()
        self.mass_note.SetLabel(['First observed is not necessarily the monoisotopic peak.',
            'The strongest isotope can change with chain length.',
            'The mean includes only isotope peaks detected above the threshold.'][mode])
        self.mass.Enable(self.group.GetValue());self.tolerance.Enable(self.group.GetValue())
        self.notice.SetLabel('Settings changed. Preview the peaks again.')
        if e:e.Skip()

    def change_source(self,e=None,initial=False):
        d=self.inputs[self.source.GetSelection()]
        self.carriers=(NEGATIVE if d.get('sign',1)<0 else POSITIVE)+[('Custom mass',None)]
        self.carrier.SetItems(['Choose a carrier…']+[name for name,_ in self.carriers])
        choice=self.options.get('carrier_label') if initial else None
        labels=[name for name,_ in self.carriers]
        self.carrier.SetSelection(labels.index(choice)+1 if choice in labels else 0)
        spec=np.asarray(d['spec'])
        same=initial and d['label']==self.options.get('source_label')
        self.low.ChangeValue(str(self.options['mz_min'] if same else float(spec[0,0])))
        self.high.ChangeValue(str(self.options['mz_max'] if same else float(spec[-1,0])))
        self.change_carrier()

    def change_carrier(self,e=None):
        self.custom.Enable(self.carrier.GetSelection()==len(self.carriers))
        self.invalidate(e)

    def calculate(self,e=None):
        self.invalidate()
        try:
            index=self.carrier.GetSelection()-1
            if index<0:raise ValueError('Choose the charge carrier for your polymer.')
            label,carrier=self.carriers[index]
            if carrier is None:
                carrier=float(self.custom.GetValue())
                if not np.isfinite(carrier) or carrier==0:raise ValueError('Enter a finite, nonzero carrier mass.')
            d=self.inputs[self.source.GetSelection()]
            options=dict(source_label=d['label'],source_event=int(d.get('e',0)),carrier_label=label,custom=self.custom.GetValue(),
                carrier=carrier,mz_min=float(self.low.GetValue()),mz_max=float(self.high.GetValue()),
                threshold_percent=float(self.threshold.GetValue()),
                isotope_tolerance=float(self.tolerance.GetValue()) if self.group.GetValue() else .01,
                profile=self.format.GetSelection()==0,group_isotopes=self.group.GetValue(),
                mass_mode=self.mass.GetSelection() if self.group.GetValue() else 0)
            kwargs={k:options[k] for k in ('carrier','mz_min','mz_max','threshold_percent','isotope_tolerance',
                                           'profile','group_isotopes','mass_mode')}
            with wx.BusyCursor():peaks=N.singly_charged(d['spec'],with_members=True,**kwargs)
            rows=[dict(q,include=True) for q in peaks]
            meta=[('Source file',getattr(self.tab.frame,'path','') or getattr(self.tab.data,'path','')),
                  ('Spectrum selection',d['label']),('Input','Singly charged spectrum; no deconvolution'),
                  ('Charge carrier',label),('Carrier mass added per charge (Da)',carrier),
                  ('Ion charge',-1 if d.get('sign',1)<0 else 1),
                  ('Mass source',MASS_NAMES[options['mass_mode']] if options['group_isotopes'] else 'Individual detected peaks'),
                  ('Mass equation','Neutral mass = m/z - signed carrier mass'),
                  ('m/z minimum',options['mz_min']),('m/z maximum',options['mz_max']),
                  ('Minimum peak above background (%)',options['threshold_percent']),('Group isotope peaks',options['group_isotopes']),
                  ('Isotope tolerance (Da)',options['isotope_tolerance']),
                  ('Signal basis','Baseline-subtracted profile peak areas' if options['profile'] else 'Summed centroid intensities'),
                  ('Isotope grouping','Adjacent peaks at 1.00335483507 m/z; no protein composition model')]
            self.preview.Freeze()
            try:
                for q in rows:
                    i=self.preview.InsertItem(self.preview.GetItemCount(),'%.6f'%q['mass'])
                    for j,v in enumerate(('%.6f'%q['mz'],str(q['isotopes']),'%.8g'%q['area'],'%.8g'%q['height']),1):
                        self.preview.SetItem(i,j,v)
            finally:self.preview.Thaw()
            self.result=(rows,meta,options);self.ok.Enable(True)
            self.notice.SetLabel('%d masses ready. Review the list; exclude unrelated species in the next window.'%len(rows))
            return self.result
        except (ValueError,RuntimeError,OSError,OverflowError) as ex:
            self.notice.SetLabel('Cannot prepare peaks: '+str(ex));return None

    def accept(self,e=None):
        if self.result is not None:self.EndModal(wx.ID_OK)

def prepare(tab,options=None):
    if callable(getattr(tab,'busy_averaging',None)) and tab.busy_averaging():
        raise ValueError('Wait for the spectrum to finish loading, then try again.')
    dlg=DirectDialog(tab,options)
    try:return dlg.result if dlg.ShowModal()==wx.ID_OK else None
    finally:dlg.Destroy()
