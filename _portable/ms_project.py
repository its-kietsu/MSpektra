"""Automatic Postrun workspace recovery and explicit Save, using safe portable archives.

What is written where
* Nothing is written for a file that was only opened. After the first edit (an undo step, or a change of an
  open polymer or kinetics window) the analysis of the window is saved as <file>_analysis\\session.msanalysis
  next to each file that was edited (a change of the Compare view counts for the file shown), and into the
  project that was opened. The project holds every open file.
* Saves follow edits; view changes (zoom, scroll, the view shown) alone never save, they go with the next save
  and with the save when a file or the window is closed.
* A project whose recovery failed, or that another window has open (<project>.lock), is not written; Save
  asks before replacing a project that failed to open and keeps the old one beside it.
"""
import datetime,os,traceback
from concurrent.futures import ThreadPoolExecutor
import wx
import ms_project_store as S

DISPLAY_KEYS=('show_charge_states','tile_mode','hrms_tile_mode','tile_scale')
OLD_SHORTCUT='Resume analysis.cmd'  # written by 3.4 to 3.55; removed when the folder is saved again

def _walk(root):
    yield root
    for child in root.GetChildren():yield from _walk(child)

def views(page):
    import unilcms as U
    cards=[];layouts=[]
    for w in _walk(page):
        if hasattr(w,'ax') and hasattr(w,'canvas'):
            cards.append(dict(x=list(w.ax.get_xlim()),y=list(w.ax.get_ylim()),lock=getattr(w,'ylock',None),
                              fixed=getattr(w,'y_fixed',None),full=getattr(w,'full',None),title=getattr(w,'title',''),subtitle=getattr(w,'subtitle','')))
        if isinstance(w,(U.SplitBox,U.StackBox)):
            if hasattr(w,'_sync'):w._sync()
            layouts.append(dict(weights=list(w.all_weights),closed=[i for i,p in enumerate(w.all_panes) if p in w.closed],
                                frac={i:v for i,p in enumerate(w.all_panes) if (v:=getattr(w,'frac',{}).get(p)) is not None}))
    sc=getattr(page,'scroller',None)
    return dict(cards=cards,layouts=layouts,scroll=list(sc.GetViewStart()) if sc else None)

def set_views(page,state):
    import unilcms as U
    cards=[w for w in _walk(page) if hasattr(w,'ax') and hasattr(w,'canvas')]
    layouts=[w for w in _walk(page) if isinstance(w,(U.SplitBox,U.StackBox))]
    for w,st in zip(layouts,state.get('layouts',[])):
        if len(st['weights'])!=len(w.all_panes):continue
        w.all_weights=list(st['weights']);w.closed={w.all_panes[i] for i in st['closed'] if 0<=i<len(w.all_panes)}
        if hasattr(w,'frac'):w.frac={w.all_panes[i]:v for i,v in st['frac'].items() if 0<=i<len(w.all_panes)}
        (w._visible if isinstance(w,U.SplitBox) else w._update_visible)();w.relayout()
    for w,st in zip(cards,state.get('cards',[])):
        w.ylock=st['lock'];w.y_fixed=st['fixed'];w.full=st['full']
        w.title=st['title'];w.subtitle=st['subtitle']
        w.ax.set_xlim(st['x']);w.ax.set_ylim(st['y']);w.canvas.draw_idle()
    sc=getattr(page,'scroller',None)
    if sc and state.get('scroll') is not None:sc.refit();sc.Scroll(*state['scroll'])

def capture(frame,docs):
    import undo as Z
    states=[]
    for d in docs:
        ms=d.attrs.get('ms');pda=d.attrs.get('pda');pn=getattr(ms,'dec',None)
        dec=None
        if pn:
            dec=dict(cfg=Z._dec_cfg(pn),active=pn.results.index(pn.active) if pn.active in pn.results else None,
                     results=[dict(res=Z._copy(e.res),visible=e.visible,sel=e.sel,view=e.view,z_weights=e.z_weights,
                                   charges=pn.charge_tile_shown(e)) for e in pn.results],shifts=pn.shift_state())
        st=dict(ms=Z.ms_get(ms) if ms else None,pda=Z.pda_get(pda) if pda else None,
                cal=Z.cal_get(ms.cal) if hasattr(ms,'cal') else None,dec=dec,
                spectra=[{k:c.get(k) for k in ('spec','range','desc')} for c in getattr(ms,'cols',[])],
                page=d.book.GetSelection(),views=[views(pg) for _,pg in d.pages],
                tool=getattr(ms,'tool','select'),finder=ms.finder_get() if hasattr(ms,'finder_get') else None,
                polymer=__import__('ms_project_tools').capture_polymer(ms))
        states.append(st)
    cmp=Z.cmp_get(frame)
    compare_view=views(frame.__dict__['_compare_tab']) if frame.__dict__.get('_compare_tab') else None
    return dict(documents=states,compare=cmp,compare_view=compare_view,side=frame._side_shown,
                active=docs.index(frame.active) if frame.active in docs else 0,
                link=Z._link_get(),shifts=Z._shiftset_get(),
                kinetics=__import__('ms_project_tools').capture_kinetics(frame),
                display={k:v for k,v in __import__('unidec_theme')._load().items() if k in DISPLAY_KEYS})

def _display(saved):
    # Only the display choices capture() writes, with values of their type.
    out={}
    for k,v in (saved.items() if isinstance(saved,dict) else ()):
        if k=='show_charge_states' and isinstance(v,bool):out[k]=v
        elif k in ('tile_mode','hrms_tile_mode') and isinstance(v,str) and len(v)<40:out[k]=v
        elif k=='tile_scale' and isinstance(v,(int,float)) and not isinstance(v,bool) and 0<v<100:out[k]=v
    return out

def _merge_shifts(saved):
    """The mass shifts and tags of a project that the user's list lacks are added to it (same masses: the same
    entry); the user's entries, their on/off choices and the finder settings stay. Returns (added, not added)."""
    import deconv_shifts as DS,ms_shifts as MS
    if not isinstance(saved,dict):return 0,0
    cur=DS.settings();new=MS.clean_settings(saved);added=missed=0
    for key,limit in (('shifts',MS.MAX_SHIFTS),('tags',MS.MAX_TAGS)):
        have={(round(e['avg'],4),round(e['mono'],4)) for e in cur[key]}
        for e in new[key]:
            k=(round(e['avg'],4),round(e['mono'],4))
            if k in have:continue
            if len(cur[key])>=limit:missed+=1;continue
            cur[key].append(dict(e));have.add(k);added+=1
    if added:DS.save_settings(cur)
    return added,missed

def _restore_doc(d,st):
    import undo as Z
    import deconv_tab as D
    ms=d.attrs.get('ms');pda=d.attrs.get('pda')
    if st.get('cal') is not None:
        avg,pick=ms.compute_average,ms.chrom_pick
        ms.compute_average=lambda *a,**k:None;ms.chrom_pick=lambda *a,**k:None
        # "Calibrate automatically on opening" is the user's choice for files opened later: it stays as it is
        try:Z.cal_set(ms.cal,dict(st['cal'],auto=bool(ms.cal.auto.GetValue())),save=False)
        finally:ms.compute_average=avg;ms.chrom_pick=pick
    if st.get('ms') is not None:
        # Preserve the processed spectrum itself instead of recomputing it or
        # launching averaging threads while restoring the saved workspace.
        avg,pick=ms.compute_average,ms.chrom_pick
        ms.compute_average=lambda *a,**k:None;ms.chrom_pick=lambda *a,**k:None
        try:Z.ms_set(ms,st['ms'])
        finally:ms.compute_average=avg;ms.chrom_pick=pick
        for col,spec in zip(ms.cols,st.get('spectra',[])):col.update(spec)
        ms.spectra_changed()
        for col in ms.cols:ms.plot_spec(col)
        if st.get('finder') is not None:ms.finder_set(st['finder'])
        ms.tool=st.get('tool','select');ms.tools.set_mode(ms.tool)
    if st.get('pda') is not None:Z.pda_set(pda,st['pda'])
    dec=st.get('dec')
    if dec is not None:
        pn=ms.dec;entries=[]
        for row in dec['results']:
            res=row['res']
            from ms_brand import display
            for key in ('method','notes','engine'):
                if isinstance(res.get(key),str):res[key]=display(res[key])
            if res.get('saved'):res['saved']=os.path.join(d.out_dir(),os.path.basename(res['saved']))
            e=D.DeconvResult(res);e.visible=row['visible'];e.sel=row['sel'];e.view=row['view'];e.z_weights=row['z_weights']
            entries.append(e)
        active=dec.get('active')
        Z.dec_set(pn,dict(results=[(e,e.visible,bool(e.res.get('label_all')),e.res.get('peaks',[])) for e in entries],
                         active=entries[active] if active is not None and 0<=active<len(entries) else None,
                         cfg=dec['cfg'],sel=[e.sel for e in entries]),save=False)
        # the per result choices only: the saved list of shifts was merged into the user's (_merge_shifts)
        sh=dec['shifts']
        pn.set_shift_state({k:v for k,v in sh.items() if k!='settings'} if isinstance(sh,dict) else sh)
        for e,row in zip(entries,dec['results']):pn._apply_z(e,row['charges'])
    d.book.SetSelection(min(st['page'],len(d.pages)-1))

def _relative(path,folder):
    # No relative path exists across drives or to a UNC share; the reader then uses the absolute path.
    try:return os.path.relpath(path,folder)
    except ValueError:return None

def _key(path):
    return os.path.normcase(os.path.abspath(path))

def _candidates(path):
    """Where the analysis of a raw data file is saved (Postrun out_dir: next to it, else in the user folder)."""
    path=path.rstrip('\\/');base=os.path.splitext(os.path.basename(path))[0]
    return [os.path.join(os.path.dirname(path),base+'_analysis',S.NAME),
            os.path.join(os.path.expanduser('~'),'MS Analysis',base,S.NAME)]

def _remove_old_shortcut(folder):
    p=os.path.join(folder,OLD_SHORTCUT)
    try:
        with open(p,encoding='utf-8',errors='replace') as f:ours=f.read(64).startswith('@echo off\nchcp 65001')
        if ours:os.remove(p)
    except OSError:pass

def ask_resume(frame,path,project,header):
    """True: resume the saved analysis of the raw data file path; False: open the file without it."""
    n=len(header.get('sources') or [])
    name=os.path.basename(path.rstrip('\\/'))
    text='%s has a saved analysis%s.'%(name,' with %d files'%n if n>1 else '')
    dlg=wx.MessageDialog(frame,text,'Saved analysis',wx.YES_NO|wx.YES_DEFAULT|wx.ICON_QUESTION)
    try:
        dlg.SetYesNoLabels('Resume the saved analysis','Open the file without it')
        return dlg.ShowModal()==wx.ID_YES
    finally:dlg.Destroy()

class Controller:
    def __init__(self,frame):
        self.frame=frame;self.executor=ThreadPoolExecutor(max_workers=1,thread_name_prefix='analysis-save')
        self.future=None;self.signature=None;self.restoring=False;self.pending=None;self.closed=False
        self.owner='%d-%d'%(os.getpid(),id(self))  # this window, in project locks
        self.blocked=set()  # projects whose recovery failed: kept unchanged until Save replaces them
        self.keep={}  # project -> suffix under which the next save keeps the file it replaces
        self.targets={}  # doc -> project it is saved into (None until the first save): the files edited
        self.locks=set()  # projects this window holds
        self.notes={}  # project -> why it is not written (status bar, once)
        self.dirty=False;self._mark=None;self._tools={}
        self._reading_project=False;self._deferred_loads=[]
        self.timer=wx.Timer(frame);frame.Bind(wx.EVT_TIMER,self.tick,self.timer);self.timer.Start(5000)
        frame.Bind(wx.EVT_CLOSE,self.close)
        self.old_load=frame.load;self.old_loaded=frame._loaded;self.old_close=frame.close_doc
        frame.load=self.load;frame._loaded=self.loaded;frame.close_doc=self.close_doc
        accepts,browse=frame.accepts,frame.browse_accepts
        frame.accepts=lambda p:S.is_project(p) or accepts(p)
        frame.browse_accepts=lambda p:S.is_project(p) or browse(p)
        import undo
        m=undo.manager(frame)
        if m is not None:m.listeners.append(self.undo_changed);self._mark=self._undo_mark(m)

    # ------------------------------------------------------------ what changed
    @staticmethod
    def _undo_mark(m):
        top=m.undo_stack[-1] if m.undo_stack else None
        return len(m.undo_stack),len(m.redo_stack),id(top) if top else None,top.t if top else None

    def undo_changed(self,m):
        """Undo listener: a new step, Undo or Redo is an edit of the files of that step."""
        mark,old=self._undo_mark(m),self._mark;self._mark=mark
        if self.closed or self.restoring or m.restoring or old is None:return
        if mark[0]>old[0] or (mark[0] and mark[2]==old[2] and mark[3]!=old[3]):step=m.undo_stack[-1]
        elif mark[1]>old[1]:step=m.redo_stack[-1]
        else:return  # a file closed (its steps removed) or the restored state became the reference
        docs=[]
        for k in step.changes:
            dom=m.domains.get(k)
            if dom is not None:docs.append(dom.doc if dom.doc is not None else self.frame.active)
        self.edited(docs)

    def edited(self,docs):
        ok=self.eligible()
        for d in docs:
            if d in ok and d not in self.targets:self.targets[d]=None
        self.dirty=True
        wx.CallLater(500,self.tick)

    def _tool_windows(self):
        """Open polymer and kinetics windows: {id: (file or None, what they show)}; cheap, no hashing."""
        out={}
        for d in self.frame.docs:
            for w in getattr(d.attrs.get('ms'),'_polymer_frames',None) or []:
                if w:out[id(w)]=(d,(getattr(w,'mode',None),id(getattr(w,'outputs',None)),w.IsShown()))
        tab=self.frame.__dict__.get('_compare_tab')
        for w in getattr(tab,'_area_frames',None) or []:
            if w and type(w).__name__=='KineticsFrame':out[id(w)]=(None,(id(getattr(w,'fit_result',None)),w.IsShown()))
        return out

    def _tools_changed(self):
        new,old=self._tool_windows(),self._tools;self._tools=new
        changed=[v[0] for k,v in new.items() if old.get(k)!=v]+[v[0] for k,v in old.items() if k not in new]
        docs=[d if d is not None else self.frame.active for d in changed]
        if docs:self.edited(docs)
        return bool(docs)

    # ------------------------------------------------------------ state
    def eligible(self):
        fr=self.frame
        return [d for d in fr.docs if d.path and not d.loading and
                (d.attrs.get('ms_data') is not None or d.attrs.get('pda_data') is not None)]

    def busy(self):
        return bool(getattr(self.frame,'_load_queue',None)) or bool(getattr(self.frame,'_load_scheduled',False)) or any(d.loading or any(getattr(p,'_averaging',False) or getattr(getattr(p,'dec',None),'busy',False)
                                   for _,p in d.pages) for d in self.frame.docs)

    def project_of(self,doc):
        """The project the file is saved into, without creating its folder."""
        p=self.targets.get(doc)
        if p:return p
        folder=doc.attrs.get('analysis_folder')
        return os.path.join(folder,S.NAME) if folder else _candidates(doc.path)[0]

    def is_blocked(self,doc):
        if self.targets.get(doc) or doc.attrs.get('analysis_folder'):return _key(self.project_of(doc)) in self.blocked
        return any(_key(p) in self.blocked for p in _candidates(doc.path))

    def tick(self,e=None):
        if self.closed or self.restoring or self.busy():return
        if self.future and not self.future.done():return
        if not self.dirty and not self._tools_changed():return
        try:self.save()
        except Exception as ex:self.dirty=True;self.error('Automatic save failed',ex,dialog=False)

    def save(self,manual=False,synchronous=False,closing=False):
        """manual: Save (also into the file shown); closing: a file or the window closes (view changes are
        saved too). Otherwise an automatic save after an edit."""
        if self.restoring:raise ValueError('Wait until the saved analysis has finished opening')
        if self.busy():
            if manual:raise ValueError('Wait for loading, averaging or deconvolution to finish')
            return None
        docs=self.eligible()
        if not docs:
            if manual:raise ValueError('Open a data file first')
            return None
        for d in [d for d in self.targets if d not in docs]:self.targets.pop(d,None)
        if manual:
            shown=self.frame.active if self.frame.active in docs else docs[0]
            self.targets.setdefault(shown,None)
        if not self.targets:self.dirty=False;return None
        self.dirty=False  # edits made from now on make the next save
        state=capture(self.frame,docs);enc=S.Encoder(docs);signature,tree=enc.signature(state)
        if signature==self.signature and not manual:return None
        import unidec_theme as T
        header=dict(format=S.FORMAT,version=S.VERSION,app_version=T.APP_VERSION,
                    kind='hrms' if type(self.frame).__name__=='HRMSFrame' else 'lcms',
                    sources=[dict(absolute=d.path) for d in docs],state=tree,array_count=len(enc.arrays),
                    saved=datetime.datetime.now(datetime.timezone.utc).isoformat())
        S.check(header,enc.arrays)
        jobs=[];skipped=[];seen=set()
        for d in docs:
            if d not in self.targets:continue
            path=self.targets[d] or os.path.join(d.out_dir(),S.NAME)
            self.targets[d]=path
            if _key(path) in seen:continue
            seen.add(_key(path))
            if _key(path) in self.blocked:
                skipped.append((path,'kept unchanged, it could not be opened'));continue
            if not S.lock(path,self.owner):
                if manual:raise ValueError('%s is open in another window. Close it there first.'%path)
                skipped.append((path,'open in another window'));continue
            self.locks.add(_key(path))
            jobs.append((docs.index(d),path,self.keep.get(_key(path))))
        for path,why in skipped:
            if self.notes.get(_key(path))!=why:
                self.notes[_key(path)]=why;self.frame.SetStatusText('Not saved to %s: %s'%(path,why),0)
        if not jobs:
            if manual:raise ValueError('Nothing could be saved')
            return None
        # One immutable snapshot is written on the worker; no wx access there.
        arrays=[a.copy() for a in enc.arrays]
        sources=header['sources']
        def write_all():
            paths=[]
            for anchor,path,keep in jobs:
                folder=os.path.dirname(path)
                h=dict(header,anchor=anchor,sources=[dict(s,relative=_relative(s['absolute'],folder)) for s in sources])
                S.write(path,h,arrays,keep=keep);paths.append(path)
                _remove_old_shortcut(folder)
            return paths
        if self.future:
            try:self.future.result()
            except Exception:pass  # a failed earlier checkpoint must not prevent a retry
        if synchronous:
            paths=write_all();self.saved(signature,paths,manual,notify=False);return paths
        self.future=self.executor.submit(write_all)
        def done(f):
            try:paths=f.result();wx.CallAfter(self.saved,signature,paths,manual)
            except Exception as ex:wx.CallAfter(self.failed_save,str(ex),manual)
        self.future.add_done_callback(done)
        if manual:self.frame.SetStatusText('Saving analysis…',1)
        return self.future

    def saved(self,signature,paths,manual,notify=True):
        for p in paths:
            if self.keep.pop(_key(p),None):self.notes.pop(_key(p),None)
        if self.closed or not notify:
            self.signature=signature;return
        self.signature=signature;self.frame.SetStatusText('Analysis saved' if manual else 'Analysis saved automatically',1)
        self.frame.files.refresh_folder(force=True)

    def failed_save(self,ex,manual):
        self.dirty=True  # tried again at the next tick
        self.error('Analysis save failed',ex,manual)

    def error(self,title,ex,dialog=True):
        print(title+': '+str(ex))
        if self.closed:return
        self.frame.SetStatusText(title+': '+str(ex),0)
        self.frame.SetStatusText('Analysis NOT saved',1)
        if dialog:wx.MessageBox(title+'\n\n'+str(ex),S.FORMAT,wx.OK|wx.ICON_ERROR,self.frame)

    def close_unsaved(self,ex):
        # A save that cannot succeed (no write access, a vanished folder) must not trap the user.
        self.error('Analysis could not be saved',ex,dialog=False)
        return wx.MessageBox('The analysis could not be saved: %s\n\nClose without saving?'%ex,
                             S.FORMAT,wx.YES_NO|wx.NO_DEFAULT|wx.ICON_WARNING,self.frame)==wx.YES

    def release(self,path):
        if _key(path) in self.locks:
            self.locks.discard(_key(path));S.unlock(path,self.owner)

    # ------------------------------------------------------------ opening
    def saved_project(self,path):
        """The saved analysis of the raw data file path: next to it or in the user folder, and only when the
        file is one of its sources."""
        for p in _candidates(path):
            if not os.path.isfile(p):continue
            try:h=S.peek(p)
            except Exception as ex:
                print('Saved analysis %s not read: %s'%(p,ex))
                self.blocked.add(_key(p))  # damaged or newer: kept unchanged, never saved over
                continue
            if S.contains(h,p,path):return p,h
        return None,None

    def load(self,path):
        path=os.path.abspath(path)
        if self.restoring:
            if path not in self._deferred_loads:self._deferred_loads.append(path)
            return
        if S.is_project(path):return self.open_project(path)
        for d in self.frame.docs:
            if d.path and os.path.normcase(d.path)==os.path.normcase(path):return self.old_load(path)
        project,h=self.saved_project(path)
        if project and _key(project) not in self.blocked:
            if ask_resume(self.frame,path,project,h):return self.open_project(project)
            self.keep[_key(project)]='replaced'  # opened without it: the first save keeps it beside the new one
        return self.old_load(path)

    def open_project(self,path):
        if self.restoring:
            if path not in self._deferred_loads:self._deferred_loads.append(path)
            return
        self.restoring=True;self._reading_project=True
        self.frame.SetStatusText('Opening saved analysis…',0)
        def read():
            h,arrays=S.read(path)
            return h,arrays,S.resolve_sources(h,path)
        def completed(future):
            try:
                data=future.result()
                wx.CallAfter(self._project_read,path,data,None)
            except Exception as ex:wx.CallAfter(self._project_read,path,None,str(ex))
        self.executor.submit(read).add_done_callback(completed)

    def _resume_deferred(self):
        paths,self._deferred_loads=self._deferred_loads,[]
        for path in paths:wx.CallAfter(self.load,path)

    def _project_read(self,path,data,error):
        if self.closed or not self.frame:return
        self._reading_project=False
        try:
            if error:raise ValueError(error)
            h,arrays,paths=data
            kind='hrms' if type(self.frame).__name__=='HRMSFrame' else 'lcms'
            existing=self.eligible()
            if existing and all(any(os.path.normcase(d.path)==os.path.normcase(p) for d in existing) for p in paths):
                target=paths[h.get('anchor',0)];d=next(d for d in existing if os.path.normcase(d.path)==os.path.normcase(target))
                self.frame.activate(d);self.restoring=False;self._resume_deferred();return
            if kind!=h['kind'] or existing:
                import hrms,unilcms
                f=(hrms.HRMSFrame if h['kind']=='hrms' else unilcms.LCMSFrame)(S.project_path(path));f.Show()
                self.restoring=False;self._resume_deferred();return f
            self.restoring=True;self.pending=(h,arrays,paths,S.project_path(path))
            for p in paths:self.old_load(p)
        except Exception as ex:
            self.restoring=False;self.pending=None;self._deferred_loads=[]
            self.blocked.add(_key(S.project_path(path)))
            self.error('Could not open the saved analysis (kept unchanged)',ex)

    def loaded(self,doc,path,data,errors):
        cal=getattr(doc.attrs.get('ms'),'cal',None) if self.restoring else None
        auto=None
        if cal is not None and hasattr(cal,'auto'):
            # no automatic calibration while the saved one is restored; the box keeps the user's choice
            auto=cal.auto.GetValue();cal.auto.SetValue(False)
        try:self.old_loaded(doc,path,data,errors)
        finally:
            if auto is not None:
                try:cal.auto.SetValue(auto)
                except RuntimeError:pass
        if self.pending:wx.CallAfter(self.finish_restore)

    def finish_restore(self):
        if not self.pending or getattr(self.frame,'_load_queue',None) or getattr(self.frame,'_load_scheduled',False) or any(d.loading for d in self.frame.docs):return
        h,arrays,paths,project=self.pending
        self.pending=None  # multiple file-loaded callbacks must not restore twice
        import undo as Z
        um=Z.manager(self.frame)
        if um is not None:um.restoring=True  # the restored analysis is the starting point, never an undo step
        try:
            docs=[]
            for p in paths:
                doc=next((d for d in self.frame.docs if d.path and os.path.normcase(d.path)==os.path.normcase(p)),None)
                if doc is None:raise ValueError('Raw data could not be loaded: '+p)
                docs.append(doc)
            anchor=h.get('anchor',0)
            docs[anchor].attrs['analysis_folder']=os.path.dirname(project)
            state=S.decode(h['state'],arrays,docs)
            if len(state['documents'])!=len(docs):raise ValueError('The saved analysis does not match its raw data files')
            display=_display(state.get('display'))
            if display:__import__('unidec_theme')._save(display)
            # "Link MS and PDA times" is a setting of the program, not of an analysis: the user's choice stays.
            added,missed=_merge_shifts(state.get('shifts'))
            for d,st in zip(docs,state['documents']):_restore_doc(d,st)
            if state.get('compare') is not None:
                self.frame.compare_tab();Z.cmp_set(self.frame,state['compare'])
            self.frame.set_side(state['side'],remember=False)
            self.frame.activate(docs[anchor])
            def final():
                try:
                    for d,st in zip(docs,state['documents']):
                        for (_,pg),v in zip(d.pages,st['views']):set_views(pg,v)
                    if state.get('compare_view') is not None:set_views(self.frame.__dict__['_compare_tab'],state['compare_view'])
                    import ms_project_tools as tools
                    for d,st in zip(docs,state['documents']):tools.restore_polymer(d.attrs.get('ms'),st.get('polymer',[]))
                    tools.restore_kinetics(self.frame,state.get('kinetics',[]))
                    Z.project_restored(self.frame,docs)
                    self.restoring=False;self.pending=None;self.dirty=False
                    self.targets[docs[anchor]]=project  # updated after the next edit
                    self._tools=self._tool_windows()
                    try:
                        live=self.eligible();self.signature=S.Encoder(live).signature(capture(self.frame,live))[0]
                    except Exception:self.signature=None
                    note=''
                    if S.locked(project,self.owner):note=' (open in another window: changes here are not saved into it)'
                    if added:note+='; %d mass shift%s added to your list'%(added,'s' if added>1 else '')
                    if missed:note+='; %d mass shift%s not added (list full)'%(missed,'s' if missed>1 else '')
                    self.frame.SetStatusText('Opened '+project+note,0)
                    self.frame.SetStatusText('Saved analysis opened',1)
                    self._resume_deferred()
                except Exception as ex:self.failed_restore(ex,project)
            wx.CallAfter(final)
        except Exception as ex:self.failed_restore(ex,project)

    def failed_restore(self,ex,project):
        self.pending=None;self.restoring=False;self._deferred_loads=[]
        self.blocked.add(_key(project))
        try:
            import undo as Z
            um=Z.manager(self.frame)
            if um is not None:um.restoring=False
        except Exception:pass
        traceback.print_exc();self.error('Could not restore the saved analysis (kept unchanged)',ex)

    # ------------------------------------------------------------ closing
    def close_doc(self,doc):
        # A failed reader uses this sentinel to discard an empty loading slot.
        # Other files in the import queue must not prevent that internal cleanup.
        if doc and doc.path=='x' and doc.attrs.get('ms_data') is None and doc.attrs.get('pda_data') is None:
            return self.old_close(doc)
        if self.restoring or self.busy():
            self.frame.SetStatusText('Wait for processing to finish before closing a file',0);return
        if doc in self.targets:
            if self.is_blocked(doc):
                if wx.MessageBox('Close %s without saving? Its saved analysis could not be opened and stays '
                                 'unchanged.'%doc.file_name,S.FORMAT,wx.YES_NO|wx.NO_DEFAULT|wx.ICON_WARNING,
                                 self.frame)!=wx.YES:return
            else:
                try:self.save(synchronous=True,closing=True)
                except Exception as ex:
                    if not self.close_unsaved(ex):return
            path=self.targets.pop(doc,None)
            if path and not any(p and _key(p)==_key(path) for p in self.targets.values()):self.release(path)
        self.old_close(doc)

    def close(self,event):
        if self.closed:event.Skip();return
        if self.restoring or self.busy():
            if event.CanVeto():
                event.Veto();self.frame.SetStatusText('Wait for processing to finish, or cancel it',0);return
        try:
            if self.targets and self.eligible():self.save(synchronous=True,closing=True)
        except Exception as ex:
            if event.CanVeto() and not self.close_unsaved(ex):event.Veto();return
        self.closed=True;self.timer.Stop();self.executor.shutdown(wait=True)
        for p in list(self.locks):
            self.locks.discard(p);S.unlock(p,self.owner)
        event.Skip()

def install(frame):
    frame.__dict__['_project']=Controller(frame)

def save(frame,e=None):
    c=frame.__dict__['_project']
    try:
        docs=c.eligible()
        shown=frame.active if frame.active in docs else (docs[0] if docs else None)
        blocked=[d for d in docs if (d in c.targets or d is shown) and c.is_blocked(d)]
        if blocked:
            if wx.MessageBox('Replace the saved analysis that could not be opened? The old file is kept beside it (.before-recovery).',S.FORMAT,wx.YES_NO|wx.NO_DEFAULT|wx.ICON_WARNING,frame)!=wx.YES:return
            for d in blocked:
                for p in [c.project_of(d)]+_candidates(d.path):
                    if _key(p) in c.blocked:c.blocked.discard(_key(p));c.keep[_key(p)]='before-recovery';c.notes.pop(_key(p),None)
        return c.save(manual=True)
    except Exception as ex:c.error('Analysis could not be saved',ex)

def open_file(frame):
    dlg=wx.FileDialog(frame,'Open a saved analysis',defaultDir=frame.folder(),wildcard='MS Analysis project (*.msanalysis)|*.msanalysis',style=wx.FD_OPEN|wx.FD_FILE_MUST_EXIST)
    try:
        if dlg.ShowModal()==wx.ID_OK:frame.load(dlg.GetPath())
    finally:dlg.Destroy()

def open_folder(frame):
    dlg=wx.DirDialog(frame,'Open an analysis folder',defaultPath=frame.folder(),style=wx.DD_DIR_MUST_EXIST)
    try:
        if dlg.ShowModal()==wx.ID_OK:
            path=dlg.GetPath()
            if not S.is_project(path):raise ValueError('This folder has no saved analysis (session.msanalysis)')
            frame.load(path)
    except Exception as ex:frame.__dict__['_project'].error('Analysis folder could not be opened',ex)
    finally:dlg.Destroy()
