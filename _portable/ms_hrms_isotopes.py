"""HRMS detail tile. Peak finding/grouping stays in the existing native engine."""
import csv
import numpy as np
import wx
import unilcms as U
import deconv_tab as D


def species_isotopes(ent, index):
    """Detected isotope peaks within the engine's species boundaries.

    This is a display list, not a predicted pattern or a monoisotopic assignment.
    Cache it on the result view, without changing or saving numerical results.
    """
    res = ent.res
    peaks = res.get('peaks') or []
    if index is None or not 0 <= index < len(peaks) or not D._grouped(res):
        return []
    q = peaks[index]
    if q.get('added') or not q.get('n_iso'):
        return []
    mass = res.get('mass')
    if mass is None or len(mass) < 5:
        return []
    threshold = float(res.get('params', {}).get('peak_thresh', .1))
    key = (id(mass), len(mass), threshold)
    cache = getattr(ent, '_isotope_details', None)
    if cache is None or cache[0] != key:
        ranges = []
        groups = D.D.group_isotopes(mass, threshold, ranges=ranges)
        cache = (key, groups, ranges, {})
        ent._isotope_details = cache
    _, groups, ranges, rows = cache
    signature = (float(q['mass']), int(q['n_iso']))
    if signature in rows:
        return rows[signature]
    if not groups:
        return []
    j = min(range(len(groups)), key=lambda k: abs(groups[k]['mass'] - q['mass']))
    # Do not attribute a neighbouring species to a hand-edited/legacy result.
    if abs(groups[j]['mass'] - q['mass']) > max(.05, float(res.get('params', {}).get('mass_step', .01)) * 2):
        return []
    lo, hi = ranges[j]
    a, b = np.searchsorted(mass[:, 0], [lo - .45, hi + .45])
    segment = mass[a:b]
    if len(segment) < 5 or np.max(segment[:, 1]) <= 0:
        return []
    cutoff = threshold * .2 * float(np.max(mass[:, 1])) / float(np.max(segment[:, 1]))
    found = D.D.isotope_peaks(segment, threshold=cutoff, max_n=max(400, int(q['n_iso']) * 4))
    # Use the same apex masses as the spectrum's labels when available.
    labels = res.get('isotope_peaks') or []
    out = []
    for x, height in sorted(found):
        if not lo - .3 <= x <= hi + .3:
            continue
        k = D._nearest_iso(labels, x, .3)
        if k is not None:
            x, height = labels[k]
        if np.isfinite(x) and np.isfinite(height) and height > 0:
            out.append((float(x), float(height)))
    out = sorted(set(out))
    rows[signature] = out
    return out


class DetailsCard(wx.Panel):
    COLS = [('#', 34), ('Isotope mass (Da)', 145), ('Intensity', 105),
            ('Relative %', 90), ('Spacing (Da)', 100)]

    def __init__(self, tab):
        super().__init__(tab.bottom, style=wx.BORDER_NONE)
        self.tab = tab
        self.SetBackgroundColour(U.C['bg'])
        self.mode = wx.Choice(self, choices=['Isotope peaks', 'Chromatogram peaks'])
        self.mode.SetMinSize(wx.Size(10, -1))
        self.mode.SetSelection(1)
        self.peaks = tab.table
        self.peaks.Reparent(self)
        self.isotopes = D.TableCard(self, 'Isotope peaks', self.COLS,
                                  on_select=self.select_isotope,
                                  empty='Select a mass in the Masses table')
        self.isotopes.extra_menu = self.menu_items
        self.isotopes.Bind(wx.EVT_MOUSEWHEEL, lambda e: U.forward_wheel(self, e))
        self.values = []
        self._updating = False
        self._ent = None
        self._last = None
        s = wx.BoxSizer(wx.VERTICAL)
        s.Add(self.mode, 0, wx.EXPAND | wx.BOTTOM, self.FromDIP(4))
        s.Add(self.peaks, 1, wx.EXPAND)
        s.Add(self.isotopes, 1, wx.EXPAND)
        self.SetSizer(s)
        self.mode.Bind(wx.EVT_CHOICE, lambda e: self.show_mode())
        self.show_mode()

    def show_mode(self):
        isotope = self.mode.GetSelection() == 0
        self.peaks.Show(not isotope)
        self.isotopes.Show(isotope)
        self.Layout()

    def refresh_result(self, ent):
        if ent is not self._ent:
            self._ent = ent
            self.mode.SetSelection(0 if ent is not None else 1)
            self.show_mode()
        if ent is None:
            self.values = []
            self.isotopes.set_rows([])
            self._last = None
            return
        sel = ent.sel
        peaks = ent.res.get('peaks') or []
        q = peaks[sel] if sel is not None and 0 <= sel < len(peaks) else None
        key = (id(ent), sel, D._display_apex(q) if q else None,
               tuple((p['mass'], p.get('n_iso')) for p in peaks))
        if key == self._last:
            return
        self._last = key
        self._updating = True
        try:
            self.values = species_isotopes(ent, sel)
            top = max((h for _, h in self.values), default=1.)
            rows = [[str(i + 1), '%.4f' % x, '%.8g' % h, '%.2f' % (100 * h / top),
                     '' if i == 0 else '%.4f' % (x - self.values[i - 1][0])]
                    for i, (x, h) in enumerate(self.values)]
            self.isotopes.empty = ('No resolved isotope peaks available' if q else
                                   'Select a mass in the Masses table')
            self.isotopes.set_rows(rows, ('Selected mass %.4f Da' % D._apex_mass(ent.res, q))
                                   if q and rows else '')
            if q and rows:
                k = D._nearest_iso(self.values, D._apex_mass(ent.res, q))
                if k is not None:
                    self.isotopes.select(k)
        finally:
            self._updating = False

    def select_isotope(self, index):
        dec = self.tab.dec
        if self._updating or self._ent is not dec.active or not 0 <= index < len(self.values):
            return
        if dec.sel is None or not 0 <= dec.sel < len(dec.result['peaks']):
            return
        q = dec.result['peaks'][dec.sel]
        x, h = self.values[index]
        if abs(D._display_apex(q) - x) < 1e-9:
            return
        q['selected_apex'] = x
        dec.plot_results()
        dec._replot_spec()
        self.tab.status('Isotope %.4f Da; intensity %.8g' % (x, h))

    def menu_items(self):
        return [('Export isotope table…', self.export)] if self.values else []

    def export(self):
        if not self.values:
            return
        path = self.tab.frame.ask_save('Export isotope table',
                                      self.tab.frame.base_name() + '_isotopes.csv',
                                      'CSV (*.csv)|*.csv')
        if not path:
            return
        try:
            top = max(h for _, h in self.values)
            with open(path, 'w', newline='', encoding='utf-8-sig') as fh:
                writer = csv.writer(fh)
                writer.writerow([c[0] for c in self.COLS])
                for i, (x, h) in enumerate(self.values):
                    writer.writerow([i + 1, x, h, 100 * h / top,
                                     '' if i == 0 else x - self.values[i - 1][0]])
            self.tab.status('Isotope table exported to ' + path)
        except OSError as ex:
            wx.MessageBox('The isotope table could not be saved:\n' + str(ex),
                          'Export isotope table', wx.OK | wx.ICON_ERROR)
