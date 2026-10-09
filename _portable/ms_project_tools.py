"""Snapshots of open kinetics and polymer tools; numerical engines stay unchanged."""
import wx
import numpy as np

def capture_polymer(tab):
    result=[]
    for w in getattr(tab,'_polymer_frames',[]):
        if not w or not w.IsShown():continue
        grid=[[w.input.GetCellValue(i,j) for j in range(5)] for i in range(w.input.GetNumberRows())]
        tables=[]
        for t in (w.repeat_table,w.chain_table):
            tables.append(([t.GetColumn(j).GetText() for j in range(t.GetColumnCount())],
                           [[t.GetItemText(i,j) for j in range(t.GetColumnCount())] for i in range(t.GetItemCount())]))
        lines=[dict(x=np.asarray(l.get_xdata()),y=np.asarray(l.get_ydata()),colour=l.get_color(),
                    marker=l.get_marker(),style=l.get_linestyle(),label=l.get_label()) for l in w.dist_card.ax.lines]
        result.append(dict(metadata=w.metadata,grid=grid,weight=w.weight.GetSelection(),basis=w.basis.GetSelection(),
                           fields=[c.GetValue() for c in (w.lo,w.hi,w.gap_tol,w.rmass,w.ends,w.match_tol)],
                           page=w.book.GetSelection(),outputs=w.outputs,mode=w.mode,tables=tables,lines=lines,
                           text=w.result_text.GetValue(),status=w.status.GetLabel(),view=__import__('ms_project').views(w)))
    return result

def restore_polymer(tab,states):
    from ms_polymer_ui import PolymerFrame
    for st in states:
        # The editable table, including excluded rows and response factors, is retained verbatim.
        rows=[dict(include=r[0]=='1',mass=0,area=0,height=0) for r in st['grid']]
        w=PolymerFrame(tab,rows,st['metadata'],st['page'])
        for i,row in enumerate(st['grid']):
            for j,value in enumerate(row):w.input.SetCellValue(i,j,value)
        w.weight.SetSelection(st['weight']);w.basis.SetSelection(st['basis'])
        for c,v in zip((w.lo,w.hi,w.gap_tol,w.rmass,w.ends,w.match_tol),st['fields']):c.ChangeValue(v)
        for table,(headers,rows) in zip((w.repeat_table,w.chain_table),st['tables']):w.fill(table,headers,rows)
        w.dist_card.reset()
        for line in st['lines']:
            w.dist_card.ax.plot(line['x'],line['y'],color=line['colour'],marker=line['marker'],linestyle=line['style'],label=line['label'])
        if st['lines']:
            import unilcms as U
            U.style_axes(w.dist_card.ax,'Species mass (Da)','Fraction');w.dist_card.ax.legend(fontsize=8)
        w.outputs=st['outputs'];w.mode=st['mode'];w.on_page()
        w.result_text.ChangeValue(st['text']);w.status.SetLabel(st['status']);w.dist_card.draw()
        if not hasattr(tab,'_polymer_frames'):tab._polymer_frames=[]
        tab._polymer_frames.append(w)
        w.Bind(wx.EVT_CLOSE,lambda e,f=w:(tab._polymer_frames.remove(f) if f in tab._polymer_frames else None,e.Skip()))
        w.Show();__import__('ms_project').set_views(w,st['view'])

def _fit_fields(fit):
    return None if fit is None else {k:(list(getattr(fit,k)) if k in ('values','errors') else getattr(fit,k)) for k,_ in fit._fields_}

def _fit_from(fields):
    import ms_kinetics as N
    if fields is None:return None
    fit=N.Result()
    for k,v in fields.items():
        if not hasattr(fit,k):continue  # (a field of another version)
        if k in ('values','errors'):
            for i,x in enumerate(v):getattr(fit,k)[i]=x
        else:setattr(fit,k,v)
    if 'aicc' not in fields:  # saved before 4.1: the criteria were not kept
        fit.aic=fit.aicc=fit.k2=fit.k2_error=float('nan');fit.order=(0,1,1,2)[fit.model] if 0<=fit.model<4 else 0
    return fit

def capture_kinetics(frame):
    tab=frame.__dict__.get('_compare_tab');result=[]
    for w in getattr(tab,'_area_frames',[]):
        if not w or type(w).__name__!='KineticsFrame' or not w.IsShown():continue
        fit=w.fit_result
        fits=[(m,_fit_fields(r),f,res,fx,fy) for m,r,f,res,fx,fy in getattr(w,'fits',[])]
        best=next((i for i,(m,r,*_) in enumerate(getattr(w,'fits',[])) if r is fit and r is not None),None)
        result.append(dict(rows=w.rows,metadata=w.metadata,relative=w.relative,xlabel=w.xlabel,ylabel=w.ylabel,
                           model=w.choice_key(),offset=w.offset.GetValue(),unit=w.unit.GetValue(),
                           fit=_fit_fields(fit),data=w.fit_data,fits=fits,best=best,
                           summary=w.summary.GetValue(),view=__import__('ms_project').views(w)))
    return result

def restore_kinetics(frame,states):
    if not states:return
    import ms_kinetics as N,unilcms as U
    from ms_kinetics_ui import KineticsFrame
    tab=frame.compare_tab()
    for st in states:
        w=KineticsFrame(tab,st['rows'],st['relative'],st['metadata']);w.xlabel=st['xlabel'];w.ylabel=st['ylabel']
        w.set_choice(st['model']);w.offset.SetValue(st['offset']);w.unit.ChangeValue(st['unit']);w.on_settings()
        if st['fit'] is not None:
            fits=[(m,_fit_from(r),f,res,fx,fy) for m,r,f,res,fx,fy in st.get('fits') or []]
            best=st.get('best')
            fit=fits[best][1] if best is not None and 0<=best<len(fits) else _fit_from(st['fit'])
            w.fits=fits;w.fit_result=fit;w.fit_data=st['data']
            w.draw_fit();w.show_summary();w.export_button.Enable(True)
        else:w.summary.ChangeValue(st['summary'])
        tab._area_frames.append(w)
        w.Bind(wx.EVT_CLOSE,lambda e,f=w:(tab._area_frames.remove(f) if f in tab._area_frames else None,e.Skip()))
        w.Show();__import__('ms_project').set_views(w,st['view'])
        w.card.set_title('Kinetic fit', '')
        w.residual_card.set_title('Fit difference', '')
