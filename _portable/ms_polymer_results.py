"""Linked HRMS polymer tables. Calculations remain in the native polymer engines."""
import copy
import math
import os
import uuid
import numpy as np
import wx
import unilcms as U
from unidec_theme import C, ui_font
from deconv_tab import TableCard

FOCUS = '#E0602F'  # existing deconvolution fit/selection colour


def file_key(tab):
    path = getattr(getattr(tab, 'data', None), 'path', None) or getattr(tab.frame, 'path', '') or '.'
    return os.path.normcase(os.path.abspath(path))


def same_file(tab, key):
    """Results of this file: the same path, or the same file name (a data file moved with its analysis to
    another folder or drive). Whether they belong to the spectrum shown is checked on the spectra themselves."""
    here = file_key(tab)
    if not isinstance(key, str):
        return False
    key = os.path.normcase(key)
    return key == here or os.path.basename(key.rstrip(os.sep + '/')) == os.path.basename(here.rstrip(os.sep + '/'))


def rebind(tab, series):
    """Saved results opened with their analysis belong to the file of that analysis, wherever it is now (the
    analysis resolves its data file relative to itself): their file key is set to it."""
    key = file_key(tab)
    for run in series or []:
        source = run.get('source') if isinstance(run, dict) else None
        if isinstance(source, dict):
            source['file'] = key


def capture_source(tab, options=None):
    """Keep the actual plotted spectra, not a later scan with a similar label."""
    inputs = list(tab.deconv_inputs())
    if options is not None:
        d = next((d for d in inputs if d['label'] == options.get('source_label') and
                  d.get('e', 0) == options.get('source_event', d.get('e', 0))), None)
        if d is None:
            raise ValueError('The selected spectrum changed. Choose the spectrum again.')
        event = int(d.get('e', 0))
        label = d['label']
        kind, result = 'direct', None
    else:
        result = getattr(tab.dec, 'result', None)
        if not result:
            raise ValueError('No deconvolution result is selected.')
        event, label, kind = int(result.get('e', 0)), result.get('input', ''), 'deconv'
        from deconv_tab import _spec_key
        d = next((d for d in inputs if d.get('e', 0) == event and
                  (result.get('spec_key') is None or _spec_key(event, d['spec']) == result['spec_key'])), None)
        if d is None:
            raise ValueError('Show the spectrum used for this deconvolution before assigning polymer chains.')
    columns = []
    for col in getattr(tab, 'cols', []):
        spec = col.get('spec')
        columns.append(dict(e=col['e'], spec=None if spec is None else np.array(spec, copy=True),
                            range=col.get('range'), desc=col.get('desc', '')))
    return dict(kind=kind, event=event, label=label,
                file=file_key(tab),
                columns=columns, avg=getattr(tab, 'avg', None), bg=getattr(tab, 'bg', None),
                pick_t=getattr(tab, 'pick_t', None),
                use_bg=tab.use_bg.GetValue() if hasattr(tab, 'use_bg') else False,
                result=copy.deepcopy(result) if result is not None else None)


def make_assignment(frame, inputs, assigned, repeat, ends, tolerance):
    """Join native assignments to immutable source-peak records by input row ID."""
    if frame.link_source is None:
        raise ValueError('Choose the source spectrum again before linking these results.')
    chains = []
    for j, values in enumerate(inputs):
        if not assigned['matched'][j]:
            continue
        row_id, mass, signal, response = values
        original = frame.link_rows[row_id-1]
        q = copy.deepcopy(original)
        q.update(row=row_id, mass=mass, signal=signal, response=response,
                 n=int(assigned['dp'][j]), delta=float(assigned['delta'][j]),
                 ppm=float(assigned['ppm'][j]))
        if frame.direct_options is not None and not q.get('members'):
            raise ValueError('Isotope links are missing. Choose the spectrum again to refresh them.')
        chains.append(q)
    if not chains:
        raise ValueError('No chains match these settings. Check the repeat and end-group masses.')
    source = copy.deepcopy(frame.link_source)
    return dict(id=uuid.uuid4().hex, repeat=float(repeat), ends=float(ends), tolerance=float(tolerance),
                chains=chains, rows=copy.deepcopy(frame.link_rows), metadata=list(frame.metadata),
                grid=[[frame.input.GetCellValue(i,j) for j in range(frame.input.GetNumberCols())]
                      for i in range(frame.input.GetNumberRows())],
                options=copy.deepcopy(frame.direct_options), source=source,
                weight=frame.weight.GetSelection(), basis=frame.basis.GetSelection())


def ensure_results(tab):
    dock = getattr(tab, 'polymer_results', None)
    if dock is None:
        dock = tab.polymer_results = PolymerResults(tab)
        # Keep this row alive across native plot-column rebuilds. HRMS lays it out.
        tab.extras['polymer'] = (dock, None, 2.8)
        tab.panels.append(dock)
    return dock


class PolymerResults(wx.Panel):
    own_refresh = True
    polymer_dock = True

    def __init__(self, tab):
        super().__init__(tab.rows, style=wx.BORDER_NONE)
        self.tab = tab
        self.SetBackgroundColour(C['bg'])
        self.SetFont(ui_font(9))
        self.series, self.active, self.chain = [], None, None
        self.visible = False
        self.mode, self.outputs = 'Polymer series', None
        self._filling, self._restoring_source = False, False
        self._bindings = {}
        self.caption = wx.StaticText(self, label='Polymer analysis')
        self.caption.SetFont(ui_font(10, 700))
        self.note = wx.StaticText(self, label='', style=wx.ST_ELLIPSIZE_END)
        self.note.SetForegroundColour(wx.Colour(C['muted']))
        self.note.SetMinSize(self.FromDIP(wx.Size(180, -1)))
        bar = wx.WrapSizer(wx.HORIZONTAL)
        bar.Add(self.caption, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, self.FromDIP(14))
        bar.Add(self.note, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, self.FromDIP(14))
        for label, handler in [('Add series…', self.add_series), ('Export Excel…', self.export),
                               ('Close results', self.hide_results)]:
            bar.Add(U.flat(self, label, handler=handler), 0, wx.RIGHT, self.FromDIP(6))
        self.split = U.SplitBox(self, wx.HORIZONTAL, min_size=120)
        self.series_table = TableCard(self.split, 'Polymer series',
            [('Repeat (Da)', 100), ('Ends (Da)', 95), ('Chains', 58)], on_select=self.select_series,
            empty='Assign a polymer series')
        self.chain_table = TableCard(self.split, 'Chain peaks',
            [('Units (n)', 68), ('Mass (Da)', 116), ('Source m/z', 118), ('Intensity %', 88),
             ('Area / signal', 112), ('Error (Da)', 92)], on_select=self.select_chain)
        self.series_table.extra_menu = lambda: [('Remove selected series', self.remove_active)] if self.series else []
        self.chain_table.extra_menu = lambda: [('Export polymer results…', self.export)] if self.series else []
        self.split.set_panes([self.series_table, self.chain_table], [1., 1.5])
        layout = wx.BoxSizer(wx.VERTICAL)
        layout.Add(bar, 0, wx.EXPAND | wx.BOTTOM, self.FromDIP(5))
        layout.Add(self.split, 1, wx.EXPAND)
        self.SetSizer(layout)
        self.SetMinSize(self.FromDIP(wx.Size(100, 190)))
        self.Hide()

    def _watch(self):
        if not hasattr(self.tab, '_polymer_frames'):
            self.tab._polymer_frames = []
        if self not in self.tab._polymer_frames:
            self.tab._polymer_frames.append(self)

    def _changed(self):
        # Existing autosave watches the identity of this token, not selections/zoom.
        self.outputs = object()
        self._watch()

    def current(self):
        return self.series[self.active] if self.active is not None and 0 <= self.active < len(self.series) else None

    def add(self, record):
        if not same_file(self.tab, record['source']['file']):
            raise ValueError('These polymer results belong to another file. Reopen the source file first.')
        if callable(getattr(self.tab, 'busy_averaging', None)) and self.tab.busy_averaging():
            raise ValueError('Wait for averaging to finish before showing polymer results.')
        # Reassigning the same source and settings updates its row, without duplicates.
        from deconv_tab import _spec_key
        def key(s):
            src = s['source']
            spec = next((c['spec'] for c in src['columns'] if c['e'] == src['event']), None)
            return (src['kind'], src['file'], src['event'], _spec_key(src['event'], spec),
                    s['repeat'], s['ends'], s['tolerance'], repr(s['options']))
        index = next((i for i, old in enumerate(self.series) if key(old) == key(record)), None)
        if index is None:
            if len(self.series) >= 100:
                raise ValueError('Keep at most 100 polymer series in one file.')
            self.series.append(record)
            index = len(self.series)-1
        else:
            record['id'] = self.series[index]['id']
            self.series[index] = record
        self.visible = True
        self.active, self.chain = index, None
        self._bindings.clear()
        self._changed()
        self.refill()
        self.tab._layout_rows()
        self.select_series(index)

    def refill(self):
        self._filling = True
        try:
            self.series_table.set_rows([['%.6f' % s['repeat'], '%.6f' % s['ends'], len(s['chains'])]
                                        for s in self.series], '%d assigned' % len(self.series) if self.series else '')
            if self.active is not None:
                self.series_table.select(self.active)
            current = self.current()
            rows = current['chains'] if current else []
            top = max((q.get('height', 0.) for q in rows), default=0.) or 1.
            direct = current and current['source']['kind'] == 'direct'
            self.chain_table.set_columns([('Units (n)', 68), ('Mass (Da)', 116),
                ('Source m/z' if direct else 'Charge-state m/z', 118), ('Intensity %', 88),
                ('Area / signal', 112), ('Error (Da)', 92)])
            self.chain_table.set_rows([[q['n'], '%.6f' % q['mass'],
                ('%.6f' % q['mz']) if q.get('mz') is not None else '—',
                '%.2f' % (100*q.get('height', 0.)/top), '%.8g' % q.get('area', 0.), '%+.6f' % q['delta']]
                for q in rows], 'Repeat %.6f Da' % current['repeat'] if current else '')
            if self.chain is not None:
                self.chain_table.select(self.chain)
            label = current['source']['label'] if current else ''
            self.note.SetLabel(label[:90])
            self.Layout()
        finally:
            self._filling = False

    def select_series(self, index):
        if self._filling or not 0 <= index < len(self.series):
            return
        self.active, self.chain = index, None
        for doc in getattr(self.tab.frame, 'docs', []):
            if doc.attrs.get('ms') is self.tab:
                self.tab.frame.activate(doc)
                break
        self.refill()
        self.show_source()
        self.redraw()
        s = self.current()
        self.tab.status('Repeat %.6f Da; %d matching chains' % (s['repeat'], len(s['chains'])))
        self._keep_visible()

    def select_chain(self, index):
        s = self.current()
        if self._filling or s is None or not 0 <= index < len(s['chains']):
            return
        self.chain = index
        self._filling = True
        try:
            self.chain_table.select(index)
        finally:
            self._filling = False
        self.show_source()
        self.redraw()
        q = s['chains'][index]
        self.tab.status('n = %d; mass %.6f Da%s' % (q['n'], q['mass'],
                        '; m/z %.6f' % q['mz'] if q.get('mz') is not None else ''))
        self._keep_visible()

    def _keep_visible(self):
        current = self.current()
        col = next((c for c in self.tab.cols if current and c['e'] == current['source']['event']), None)
        if col and hasattr(self.tab, 'keep_in_view'):
            card = col['spec_card']
            self.tab.keep_in_view([card.pair or card, self])

    def source_matches(self, col, record=None):
        s = self.current() if record is None else record
        if s is None or not same_file(self.tab, s['source']['file']) or col['e'] != s['source']['event'] or col.get('spec') is None:
            return False
        key = (s['id'], col['e'])
        spec = col['spec']
        cached = self._bindings.get(key)
        if cached is not None and cached[0] is spec:
            return cached[1]
        original = next((c['spec'] for c in s['source']['columns'] if c['e'] == col['e']), None)
        same = original is not None and np.array_equal(spec, original)
        self._bindings[key] = (spec, same)
        return same

    def show_source(self):
        s = self.current()
        if s is None or self._restoring_source:
            return
        if callable(getattr(self.tab, 'busy_averaging', None)) and self.tab.busy_averaging():
            self.tab.status('Wait for averaging to finish before showing this polymer spectrum.')
            return
        source = s['source']
        if not same_file(self.tab, source.get('file')):
            return  # (results of another file: the spectra shown stay)
        col = next((c for c in self.tab.cols if c['e'] == source['event']), None)
        if col is None or self.source_matches(col):
            return
        self._restoring_source = True
        try:
            for saved in source['columns']:
                target = next((c for c in self.tab.cols if c['e'] == saved['e']), None)
                if target is None:
                    continue
                target.update(spec=None if saved['spec'] is None else np.array(saved['spec'], copy=True),
                              range=saved['range'], desc=saved['desc'], measure=None, m1=None)
                self.tab.spectrum_replaced(target)
            self.tab.avg, self.tab.bg, self.tab.pick_t = source['avg'], source['bg'], source['pick_t']
            if hasattr(self.tab, 'use_bg'):
                self.tab.use_bg.SetValue(source['use_bg'])
            if hasattr(self.tab, '_remember_spectrum_selection'):
                self.tab._remember_spectrum_selection()
                self.tab._restore_spectrum_selection()
            self.tab.spectra_changed()
            self._bindings.clear()
        finally:
            self._restoring_source = False

    def spectrum_replaced(self, col):
        self._bindings.clear()

    def redraw(self):
        for col in self.tab.cols:
            self.tab.plot_spec(col, keep_view=True)

    def marks_for(self, col):
        return bool(self.visible and self.current() and self.source_matches(col))

    def overlay(self, col, ax):
        if not self.marks_for(col):
            return
        s = self.current()
        if s['source']['kind'] == 'direct':
            # Tone down only the native raw trace; fits, labels and calibrants retain their style.
            for line in ax.lines:
                if type(line).__name__ == 'EnvelopeLine':
                    line.set_color(C['faint'])
                    line.set_alpha(.5)
                    break
            for artist in ax.collections:
                if type(artist).__name__ == 'SpectrumSticks':
                    artist.set_color(C['faint'])
                    artist.set_alpha(.5)
                    break
        x0, x1 = sorted(ax.get_xlim())
        focus = None
        for i, q in enumerate(s['chains']):
            selected = i == self.chain
            colour = FOCUS if selected else C['accent']
            if s['source']['kind'] == 'direct':
                members = q.get('members', [])
                points = [(m['mz'], m['height']) for m in members if x0 <= m['mz'] <= x1]
                if points:
                    x, y = zip(*points)
                    marks = ax.vlines(x, 0, y, color=colour, lw=1.5 if selected else 1., alpha=.9, zorder=5)
                    marks._no_scale = True
                    marks.set_gid('polymer-chain-%d' % q['row'])
                    dots, = ax.plot(x, y, ls='none', marker='o' if selected else 's', ms=4 if selected else 2.5,
                                    color=colour, zorder=6)
                    dots._no_scale = True
                    if selected:
                        focus = max(points, key=lambda p: p[1])
            else:
                # Preserve the existing deconvolution charge-state mapping and sign conventions.
                from deconv_tab import _charges_of
                for z, mz in _charges_of(q['mass'], s['source']['result']):
                    if x0 <= mz <= x1:
                        mark = ax.axvline(mz, color=colour, lw=1.3 if selected else .8,
                                         ls=(0, (3, 2)), alpha=.8, zorder=4)
                        mark._no_scale = True
                        mark.set_gid('polymer-charge-%d' % q['row'])
        if focus is not None:
            q = s['chains'][self.chain]
            label = ax.annotate('n %d' % q['n'], focus, xytext=(10, 14),
                                textcoords='offset points', ha='left', fontsize=8, color=C['text'],
                                bbox=dict(fc=C['panel'], ec='none', alpha=.85, pad=1), zorder=8)
            label._no_scale = True

    def pick(self, col, mz, y=None):
        if not self.marks_for(col):
            return False
        targets = []
        for series_index, s in enumerate(self.series):
            if not self.source_matches(col, s):
                continue
            if s['source']['kind'] == 'direct':
                targets += [(m['mz'], series_index, i) for i, q in enumerate(s['chains']) for m in q.get('members', [])]
            else:
                from deconv_tab import _charges_of
                targets += [(x, series_index, i) for i, q in enumerate(s['chains']) for _, x in _charges_of(q['mass'], s['source']['result'])]
        if not targets:
            return False
        # Hit tolerance follows screen pixels, but snap to an actual peak first for direct data.
        peak = self.tab._nearest_peak(col, mz)
        if peak is None:
            return False
        target, series_index, index = min(targets, key=lambda p: (abs(p[0]-peak), p[1] != self.active))
        tolerance = max(.004, abs(peak)*10e-6)
        if abs(target-peak) > tolerance:
            return False
        if series_index != self.active:
            self.select_series(series_index)
        self.select_chain(index)
        return True

    def add_series(self, e=None):
        s = self.current()
        if s is None:
            return
        from ms_polymer_ui import PolymerFrame, watch_frame
        frame = PolymerFrame(self.tab, copy.deepcopy(s['rows']), s['metadata'], page=2,
                             direct_options=copy.deepcopy(s['options']), link_source=copy.deepcopy(s['source']))
        for i, row in enumerate(s.get('grid', [])):
            for j, value in enumerate(row):
                frame.input.SetCellValue(i, j, value)
        frame.rmass.ChangeValue('%.15g' % s['repeat'])
        frame.ends.ChangeValue('%.15g' % s['ends'])
        frame.match_tol.ChangeValue('%.15g' % s['tolerance'])
        frame.weight.SetSelection(s['weight'])
        frame.basis.SetSelection(s['basis'])
        watch_frame(self.tab, frame)
        frame.Show()

    def show_results(self):
        if not self.series:
            return
        self.visible = True
        self.split.closed.clear()
        self.split._visible()
        self._changed()
        self.tab._layout_rows()
        chain = self.chain
        self.select_series(self.active if self.active is not None else 0)
        if chain is not None and self.current() and 0 <= chain < len(self.current()['chains']):
            self.select_chain(chain)

    def hide_results(self, e=None):
        self.visible = False
        self.Hide()
        self._changed()
        self.tab._layout_rows()
        self.redraw()

    def remove_active(self):
        if self.active is None:
            return
        del self.series[self.active]
        self.active = min(self.active, len(self.series)-1) if self.series else None
        self.chain = None
        if not self.series:
            self.visible = False
            self.Hide()
        self._changed()
        self.refill()
        self.tab._layout_rows()
        if self.current():
            self.show_source()
        self.redraw()

    def reset(self):
        self.series, self.active, self.chain = [], None, None
        self.visible = False
        self.outputs = None
        self._bindings.clear()
        self.Hide()
        self.refill()
        if self in getattr(self.tab, '_polymer_frames', []):
            self.tab._polymer_frames.remove(self)

    def snapshot(self):
        return dict(kind='polymer_dock', series=self.series, active=self.active, chain=self.chain,
                    visible=self.visible, view=__import__('ms_project').views(self.tab))

    def restore(self, state):
        self.series = state.get('series', [])
        rebind(self.tab, self.series)
        self.active = state.get('active')
        if self.active is None or not 0 <= self.active < len(self.series):
            self.active = 0 if self.series else None
        self.chain = state.get('chain')
        s = self.current()
        if not s or self.chain is None or not 0 <= self.chain < len(s['chains']):
            self.chain = None
        self.visible = bool(state.get('visible', True) and self.series)
        self._watch()
        self.outputs = object()
        self.refill()
        self.tab._layout_rows()
        self.redraw()
        if state.get('view'):
            __import__('ms_project').set_views(self.tab, state['view'])

    def export(self, e=None, path=None):
        if not self.series:
            return None
        if not path:
            path = self.tab.frame.ask_save('Export polymer results', 'polymer_series.xlsx',
                                            'Excel workbook (*.xlsx)|*.xlsx')
        if not path:
            return None
        if not path.lower().endswith('.xlsx'):
            path += '.xlsx'
        header = ['Series', 'Repeat (Da)', 'End groups (Da)', 'Units (n)', 'Mass (Da)', 'Source m/z',
                  'Peak height', 'Area / signal', 'Isotope peaks', 'Error (Da)', 'Error (ppm)', 'Spectrum']
        rows, metadata = [], [('Tool', 'Linked polymer series')]
        for i, s in enumerate(self.series):
            for q in s['chains']:
                rows.append([i+1, s['repeat'], s['ends'], q['n'], q['mass'], q.get('mz', ''),
                             q.get('height', 0), q.get('area', 0), q.get('isotopes', ''),
                             q['delta'], q['ppm'], s['source']['label']])
            metadata += [('Series %d tolerance (Da)' % (i+1), s['tolerance'])]
            metadata += [('Series %d %s' % (i+1, key), value) for key, value in s['metadata']]
            for q in s['chains']:
                for member in q.get('members', []):
                    metadata.append(('Series %d input %d isotope %d' % (i+1, q['row'], member['isotope']+1),
                                     'm/z %.15g; height %.15g; area %.15g' %
                                     (member['mz'], member['height'], member['area'])))
        try:
            __import__('lcms_compare_core').write_area_xlsx(path, header, rows, metadata)
        except (OSError, ValueError) as ex:
            self.tab.status('Export failed: ' + str(ex))
            return None
        self.tab.status('Polymer results exported to ' + path)
        return path
