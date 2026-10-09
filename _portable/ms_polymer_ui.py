"""HRMS polymer workspace. All distribution, series and DP calculations are native."""
import math
import wx
import wx.grid as G
import numpy as np
from unidec_theme import C,ui_font
import unilcms as U
import lcms_compare_core as K
import ms_polymer as N

def source_rows(tab):
    import deconv_tab as D
    res=tab.dec.result
    if not res or not res.get('peaks'):
        raise ValueError('No deconvolution result. Deconvolute a spectrum, then click Reload active result.')
    if res.get('tag')!='isodec' and res.get('params',{}).get('isotopes')=='resolved' and not D._grouped(res):
        raise ValueError('Group the isotope peaks into species first.')
    rows=[dict(mass=q['mass'],area=q.get('area',0),height=q.get('height',0),include=not q.get('added',False))
          for q in res['peaks']]
    if len(rows)>5000:raise ValueError('Too many species (at most 5000).')
    metadata=[('Source file',getattr(tab.frame,'path','') or getattr(getattr(tab,'data',None),'path','')),
                 ('Deconvolution method',res.get('tag','')),('Spectrum selection',str(res.get('input',''))),
                 ('Result name',res.get('name','')),
                 ('Mass source','Neutral species mass (envelope mass for isotope species)')]
    for key in ('mass_range','z_range','mass_step','adduct_mass','isotopes'):
        if key in res.get('params',{}):metadata.append(('Deconvolution '+key,str(res['params'][key])))
    return rows,metadata

class PolymerFrame(wx.Frame):
    def __init__(self,tab,rows,metadata,page=0):
        wx.Frame.__init__(self,wx.GetTopLevelParent(tab),title='Polymer analysis',size=tab.FromDIP(wx.Size(1040,800)))
        self.SetMinSize(self.FromDIP(wx.Size(800,650)))
        self.CentreOnParent()  # (inside the screen when it is shown: unidec_theme.fit_to_screen)
        self.SetFont(ui_font(10));self.SetBackgroundColour(C['bg'])
        self.tab=tab;self.metadata=list(metadata);self.outputs=None;self.mode=None
        self.weight=wx.Choice(self,choices=['Peak area','Peak height'])
        self.weight.SetSelection(0)
        self.basis=wx.Choice(self,choices=['Signal per molecule','Signal per unit mass'])
        self.basis.SetSelection(0)
        reload_btn=U.flat(self,'Reload active result',handler=self.reload)
        self.export_btn=U.flat(self,'Export results to Excel…',icon='export',handler=self.export)
        self.export_btn.Enable(False)
        bar=wx.BoxSizer(wx.HORIZONTAL)
        for control in (self.weight,self.basis,reload_btn,self.export_btn):
            bar.Add(control,0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(10))
        self.input=G.Grid(self)
        self.input.CreateGrid(0,5)
        for i,label in enumerate(('Use','Species mass (Da)','Integrated area','Peak height','Response factor')):
            self.input.SetColLabelValue(i,label);self.input.SetColSize(i,self.FromDIP(70 if i==0 else 150))
        self.input.SetColFormatBool(0);self.input.SetRowLabelSize(self.FromDIP(48))
        self.input.SetMinSize(self.FromDIP(wx.Size(-1,150)))
        self.status=wx.StaticText(self,label='')
        self.book=wx.Notebook(self)
        self.dist=self.page('Mass distribution')
        self.repeat=self.page('Repeat unit finder')
        self.chain=self.page('Chain lengths / end groups')
        self.result_text=wx.TextCtrl(self.dist,style=wx.TE_MULTILINE|wx.TE_READONLY,size=self.FromDIP(wx.Size(-1,65)))
        calc=U.flat(self.dist,'Calculate Mn / Mw / Đ','primary',handler=self.calculate_distribution)
        self.dist_card=U.PlotCard(self.dist,'Mass distribution','',mode='zoom')
        self.dist_card.SetMinSize(self.FromDIP(wx.Size(-1,190)))
        ds=wx.BoxSizer(wx.VERTICAL);ds.Add(calc,0,wx.ALL,self.FromDIP(8));ds.Add(self.result_text,0,wx.EXPAND|wx.LEFT|wx.RIGHT,self.FromDIP(8));ds.Add(self.dist_card,1,wx.EXPAND|wx.ALL,self.FromDIP(8));self.dist.SetSizer(ds)
        self.lo,self.hi,self.gap_tol=[wx.TextCtrl(self.repeat,value=s,size=self.FromDIP(wx.Size(85,-1))) for s in ('10','1000','0.05')]
        self.repeat_table=self.table(self.repeat,('Repeat mass (Da)','Supported gaps','One-repeat gaps','Ranking score'))
        rs=wx.BoxSizer(wx.VERTICAL);rb=wx.BoxSizer(wx.HORIZONTAL)
        for control in (wx.StaticText(self.repeat,label='Search from'),self.lo,wx.StaticText(self.repeat,label='to (Da)'),self.hi,
                        wx.StaticText(self.repeat,label='Gap tolerance (Da)'),self.gap_tol,U.flat(self.repeat,'Find candidates','primary',handler=self.find_repeats)):
            rb.Add(control,0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(8))
        rs.Add(rb,0,wx.ALL,self.FromDIP(8));rs.Add(self.repeat_table,1,wx.EXPAND|wx.ALL,self.FromDIP(8));self.repeat.SetSizer(rs)
        self.repeat_table.Bind(wx.EVT_LIST_ITEM_ACTIVATED,self.use_repeat)
        self.repeat_table.SetToolTip('Double click a row to use its repeat mass')
        self.rmass,self.ends,self.match_tol=[wx.TextCtrl(self.chain,value=s,size=self.FromDIP(wx.Size(100,-1))) for s in ('','','0.05')]
        cs=wx.BoxSizer(wx.VERTICAL);cb=wx.BoxSizer(wx.HORIZONTAL)
        for control in (wx.StaticText(self.chain,label='Repeat mass (Da)'),self.rmass,wx.StaticText(self.chain,label='End groups E (Da)'),self.ends,
                        wx.StaticText(self.chain,label='Tolerance (Da)'),self.match_tol):
            cb.Add(control,0,wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,self.FromDIP(8))
        cs.Add(cb,0,wx.ALL,self.FromDIP(8))
        buttons=wx.BoxSizer(wx.HORIZONTAL)
        for control in (U.flat(self.chain,'Assign chain lengths','primary',handler=self.assign),U.flat(self.chain,'Find end group offsets',handler=self.find_offsets)):
            buttons.Add(control,0,wx.RIGHT,self.FromDIP(8))
        cs.Add(buttons,0,wx.ALL,self.FromDIP(8))
        self.chain_table=self.table(self.chain,('Species','Mass (Da)','n / DP','Predicted (Da)','Error (Da)','Error (ppm)','Matched'))
        cs.Add(self.chain_table,1,wx.EXPAND|wx.ALL,self.FromDIP(8));self.chain.SetSizer(cs)
        layout=wx.BoxSizer(wx.VERTICAL);layout.Add(bar,0,wx.ALL,self.FromDIP(10));layout.Add(self.input,2,wx.EXPAND|wx.LEFT|wx.RIGHT,self.FromDIP(10));layout.Add(self.status,0,wx.EXPAND|wx.ALL,self.FromDIP(10));layout.Add(self.book,4,wx.EXPAND|wx.ALL,self.FromDIP(10));self.SetSizer(layout)
        self.weight.Bind(wx.EVT_CHOICE,self.invalidate);self.basis.Bind(wx.EVT_CHOICE,self.invalidate)
        self.input.Bind(G.EVT_GRID_CELL_CHANGED,self.invalidate)
        for control in (self.lo,self.hi,self.gap_tol,self.rmass,self.ends,self.match_tol):control.Bind(wx.EVT_TEXT,self.invalidate)
        self.load_rows(rows);self.book.SetSelection(page)
        self.dist_card.image_name=lambda:'polymer_distribution'
        self.dist_card.on_menu=lambda x,y:[('Calculate distribution',self.calculate_distribution),('Export distribution to Excel…',self.export)]
        self.book.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED,self.on_page)

    def on_page(self,e=None):
        pages={'Mass distribution':0,'Repeat-unit candidates':1,'Chain lengths':2,'End-group offset candidates':2}
        self.export_btn.Enable(self.outputs is not None and pages.get(self.mode)==self.book.GetSelection())
        if e:e.Skip()

    def page(self,label):
        panel=wx.ScrolledWindow(self.book,style=wx.VSCROLL)
        panel.SetBackgroundColour(C['bg']);panel.SetScrollRate(0,self.FromDIP(15));self.book.AddPage(panel,label)
        panel.Bind(wx.EVT_SIZE,lambda e,p=panel:(p.FitInside(),e.Skip()))
        return panel
    def table(self,parent,headers):
        table=wx.ListCtrl(parent,style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
        for i,h in enumerate(headers):table.InsertColumn(i,h,width=self.FromDIP(145 if i<2 else 115))
        return table
    def fill(self,table,headers,rows):
        table.DeleteAllItems()
        # Column sets differ between assignments and modulo-offset candidates.
        while table.GetColumnCount():table.DeleteColumn(0)
        for i,h in enumerate(headers):table.InsertColumn(i,h,width=self.FromDIP(145 if i<2 else 115))
        for row in rows:
            idx=table.InsertItem(table.GetItemCount(),self.format(row[0]))
            for j,v in enumerate(row[1:],1):table.SetItem(idx,j,self.format(v))
    @staticmethod
    def format(v):
        if isinstance(v,(float,np.floating)):
            return '—' if not math.isfinite(v) else '%.10g'%v
        return str(v)
    def invalidate(self,e=None):
        self.outputs=self.mode=None;self.export_btn.Enable(False)
        self.result_text.ChangeValue('')
        self.repeat_table.DeleteAllItems();self.chain_table.DeleteAllItems();self.dist_card.reset();self.dist_card.draw()
        self.status.SetLabel('Inputs changed: results cleared.' if e else '')
        if e:e.Skip()
    def load_rows(self,rows):
        if self.input.GetNumberRows():self.input.DeleteRows(0,self.input.GetNumberRows())
        if rows:self.input.AppendRows(len(rows))
        for i,r in enumerate(rows):
            for j,value in enumerate(('1' if r['include'] else '',r['mass'],r['area'],r['height'],1)):
                self.input.SetCellValue(i,j,str(value))
            self.input.SetReadOnly(i,2);self.input.SetReadOnly(i,3)
        self.invalidate()
    def reload(self,e=None):
        try:rows,meta=source_rows(self.tab)
        except ValueError as ex:self.invalidate();self.status.SetLabel(str(ex));return None
        self.metadata=meta;self.load_rows(rows);return True
    def data(self):
        if self.input.IsCellEditControlEnabled():self.input.SaveEditControlValue();self.input.HideCellEditControl();self.input.DisableCellEditControl()
        rows=[]
        for i in range(self.input.GetNumberRows()):
            use=self.input.GetCellValue(i,0)
            if use not in ('','0','1'):raise ValueError('Use must be ticked or unticked')
            if use!='1':continue
            values=[float(self.input.GetCellValue(i,j)) for j in (1,2 if self.weight.GetSelection()==0 else 3,4)]
            m,a,r=values
            if not all(math.isfinite(v) for v in values) or m<=0 or a<0 or r<=0:
                raise ValueError('Row %d: mass and response must be above 0, signal 0 or more'%(i+1))
            rows.append((i+1,m,a,r))
        if not rows:raise ValueError('Select at least one species')
        return rows
    def start(self):
        self.invalidate()
        return self.data()
    def finish(self,mode,header,rows,meta,inputs):
        self.mode=mode;self.outputs=(header,rows,self.metadata+[
            ('Tool',mode),('Weight signal',self.weight.GetStringSelection()),('Weight interpretation',self.basis.GetStringSelection()),
            ('Response correction','Signal / response factor'),
            ('Included species',len(inputs)),('Input masses (Da)',', '.join('%.15g'%v[1] for v in inputs))]+meta)
        self.on_page()
    def failure(self,ex):
        self.invalidate();self.status.SetLabel('Cannot calculate: %s'%ex);return None
    def calculate_distribution(self,e=None):
        try:
            rows=self.start();s,nf,wf=N.distribution([r[1] for r in rows],[r[2] for r in rows],[r[3] for r in rows],self.basis.GetSelection())
            self.result_text.ChangeValue('%d species with signal\nMn = %.10g g/mol    Mw = %.10g g/mol    Đ = %.10g'%(s.used,s.mn,s.mw,s.dispersity))
            self.dist_card.ax.plot([r[1] for r in rows],nf,'o',label='Number fraction',color=C['accent'])
            self.dist_card.ax.plot([r[1] for r in rows],wf,'s',label='Mass fraction',color='#E56B32')
            self.dist_card.ax.legend(fontsize=8);U.style_axes(self.dist_card.ax,'Species mass (Da)','Fraction')
            mass=np.array([r[1] for r in rows]);self.dist_card.full=(float(mass.min())-.5,float(mass.max())+.5);self.dist_card.ylock=self.dist_card._auto=None;self.dist_card.draw()
            output=[[r[0],r[1],r[2],r[3],nf[i],wf[i]] for i,r in enumerate(rows)]
            self.finish('Mass distribution',['Species','Mass (Da)','Signal','Response factor','Number fraction','Mass fraction'],output,
                        [('Mn (g/mol, estimate)',s.mn),('Mw (g/mol, estimate)',s.mw),('Dispersity',s.dispersity),('Species with signal',s.used)],rows)
            self.status.SetLabel('Mass distribution calculated.');return s
        except (ValueError,RuntimeError,OSError) as ex:return self.failure(ex)
    def find_repeats(self,e=None):
        try:
            rows=self.start();lo,hi,tol=[float(c.GetValue()) for c in (self.lo,self.hi,self.gap_tol)]
            result=N.candidates([r[1] for r in rows],lo=lo,hi=hi,tolerance=tol)
            output=[[r['mass'],r['support'],r['single'],r['score']] for r in result]
            header=['Repeat mass (Da)','Supported gaps','One-repeat gaps','Ranking score'];self.fill(self.repeat_table,header,output)
            self.finish('Repeat-unit candidates',header,output,[('Search minimum (Da)',lo),('Search maximum (Da)',hi),('Gap tolerance (Da)',tol)],rows)
            self.status.SetLabel('%d candidate spacings'%len(result));return result
        except (ValueError,RuntimeError,OSError) as ex:return self.failure(ex)
    def use_repeat(self,e):
        value=self.repeat_table.GetItemText(e.GetIndex())
        self.rmass.ChangeValue(value);self.book.SetSelection(2);self.invalidate()
    def assign(self,e=None):
        try:
            rows=self.start();repeat,ends,tol=[float(c.GetValue()) for c in (self.rmass,self.ends,self.match_tol)]
            result=N.assign([r[1] for r in rows],repeat,ends,tol)
            output=[[r[0],r[1],result['dp'][i],result['predicted'][i],result['delta'][i],result['ppm'][i],bool(result['matched'][i])] for i,r in enumerate(rows)]
            header=['Species','Mass (Da)','n / DP','Predicted (Da)','Error (Da)','Error (ppm)','Matched'];self.fill(self.chain_table,header,output)
            self.finish('Chain lengths',header,output,[('Repeat mass (Da)',repeat),('End groups E (Da)',ends),('Tolerance (Da)',tol),('Equation','M = n × repeat + E')],rows)
            self.status.SetLabel('%d of %d species match the series'%(sum(result['matched']),len(rows)));return result
        except (ValueError,RuntimeError,OSError) as ex:return self.failure(ex)
    def find_offsets(self,e=None):
        try:
            rows=self.start();repeat=float(self.rmass.GetValue());tol=float(self.match_tol.GetValue())
            result=N.candidates([r[1] for r in rows],repeat=repeat,tolerance=tol)
            output=[[r['mass'],r['support']] for r in result];header=['End group offset mod repeat (Da)','Species supported'];self.fill(self.chain_table,header,output)
            self.finish('End-group offset candidates',header,output,[('Repeat mass (Da)',repeat),('Tolerance (Da)',tol)],rows)
            self.status.SetLabel('%d end group offsets (modulo the repeat mass)'%len(result));return result
        except (ValueError,RuntimeError,OSError) as ex:return self.failure(ex)
    def export(self,e=None,path=None):
        if self.input.IsCellEditControlEnabled():
            self.input.SaveEditControlValue();self.input.HideCellEditControl();self.input.DisableCellEditControl()
        self.on_page()
        if self.outputs is None or not self.export_btn.IsEnabled():return None
        if not path:
            with wx.FileDialog(self,'Export polymer results',defaultDir=self.tab.frame.out_dir(),defaultFile='polymer_results.xlsx',wildcard='Excel workbook (*.xlsx)|*.xlsx',style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT) as dlg:
                if dlg.ShowModal()!=wx.ID_OK:return None
                path=dlg.GetPath()
            if not path.lower().endswith('.xlsx'):path+='.xlsx'
        header,rows,meta=self.outputs
        rows=[[v if not isinstance(v,(float,np.floating)) or math.isfinite(v) else '' for v in row] for row in rows]
        try:K.write_area_xlsx(path,header,rows,meta)
        except (ValueError,OSError) as ex:self.status.SetLabel('Export failed: %s'%ex);return None
        return path

def open_tools(tab,page=0):
    notice=None
    try:
        try:rows,metadata=source_rows(tab)
        except ValueError as ex:
            rows,metadata=[],[];notice=str(ex)
        frame=PolymerFrame(tab,rows,metadata,page)
    except Exception as ex:
        import traceback
        traceback.print_exc()
        wx.MessageBox('Polymer analysis could not open.\n\n'+str(ex),'Polymer analysis',wx.OK|wx.ICON_ERROR,wx.GetTopLevelParent(tab))
        return None
    if notice:
        frame.status.SetLabel(notice)
        frame.result_text.ChangeValue(notice)
        tab.status(notice)
    if not hasattr(tab,'_polymer_frames'):tab._polymer_frames=[]
    tab._polymer_frames.append(frame)
    frame.Bind(wx.EVT_CLOSE,lambda e:(tab._polymer_frames.remove(frame) if frame in tab._polymer_frames else None,e.Skip()))
    frame.Show();return frame

def menu_items(tab):
    return [('Mass distribution (Mn / Mw / Đ)…',lambda:open_tools(tab,0)),
            ('Repeat unit finder…',lambda:open_tools(tab,1)),
            ('Chain lengths / end groups…',lambda:open_tools(tab,2))]
