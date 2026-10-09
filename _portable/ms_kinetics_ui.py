"""Kinetics UI. Numerical fitting, prediction, residuals and uncertainty run in mskinetics.dll."""
import math
import numpy as np
import wx
import unilcms as U
from unidec_theme import C, ui_font
import lcms_compare_core as K
import ms_kinetics as N

def number(value):
    return 'unavailable' if not math.isfinite(value) else '%.8g' % value

class KineticsFrame(wx.Frame):
    """Owns a snapshot of the selected region; changing Compare cannot alter a saved fit."""
    def __init__(self, tab, rows, relative, metadata):
        wx.Frame.__init__(self, tab.frame, title='Kinetic fitting', size=tab.FromDIP(wx.Size(1000, 780)))
        self.SetMinSize(self.FromDIP(wx.Size(880, 620)))
        self.SetBackgroundColour(C['bg'])
        self.SetFont(ui_font(10))
        self.rows = [dict(r) for r in rows]
        self.metadata = list(metadata)
        self.relative = relative
        self.xlabel = tab.area_xlabel.GetValue().strip() or 'X'
        self.ylabel = 'Relative area (%)' if relative else 'Actual area (%s × s)' % (rows[0].get('units') or 'signal')
        self.x = np.asarray([r['x'] for r in rows], dtype=float)
        self.y = np.asarray([r['relative'] if relative else r['area'] for r in rows], dtype=float)
        self.fit_result = None
        self.fit_data = None
        self.model = wx.Choice(self, choices=['Straight line', 'Exponential decay', 'Rise to a final level'])
        self.model.SetSelection(1)
        self.offset = wx.CheckBox(self, label='Fit baseline')
        self.offset.SetValue(True)
        self.unit = wx.TextCtrl(self, value='', size=self.FromDIP(wx.Size(75, -1)))
        self.unit.SetToolTip('e.g. min, h or s')
        self.fit_button = U.flat(self, 'Fit', 'primary', handler=self.on_fit)
        self.export_button = U.flat(self, 'Export to Excel…', icon='export', handler=self.on_export)
        self.export_button.Enable(False)
        bar = wx.BoxSizer(wx.VERTICAL)
        for controls in ((wx.StaticText(self, label='Model'), self.model, self.offset, self.fit_button),
                         (wx.StaticText(self, label='X unit'), self.unit, self.export_button)):
            row = wx.BoxSizer(wx.HORIZONTAL)
            for w in controls:
                row.Add(w, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, self.FromDIP(10))
            bar.Add(row, 0, wx.EXPAND | wx.BOTTOM, self.FromDIP(6))
        # Retain the compact text snapshot for existing saved projects.
        self.summary = wx.TextCtrl(self, style=wx.TE_MULTILINE | wx.TE_READONLY)
        self.summary.Hide()
        self.body = wx.Panel(self)
        self.body.SetBackgroundColour(C['bg'])
        self.charts = wx.ScrolledWindow(self.body, style=wx.VSCROLL)
        self.charts.SetMinSize(self.FromDIP(wx.Size(350, -1)))
        self.charts.SetScrollRate(0, self.FromDIP(15))
        self.charts.SetBackgroundColour(C['bg'])
        self.card = U.PlotCard(self.charts, 'Kinetic fit', '', mode='zoom')
        self.residual_card = U.PlotCard(self.charts, 'Fit difference', '', mode='zoom')
        self.card.SetMinSize(self.FromDIP(wx.Size(-1, 240)))
        self.residual_card.SetMinSize(self.FromDIP(wx.Size(-1, 180)))
        self.card.image_name = lambda: tab.image_name() + '_kinetic_fit'
        self.residual_card.image_name = lambda: tab.image_name() + '_kinetic_residuals'
        self.card.on_menu = self.residual_card.on_menu = lambda x, y: [('Export to Excel…', self.on_export)]
        from deconv_tab import TableCard
        self.result_side = wx.Panel(self.body)
        self.result_side.SetBackgroundColour(C['bg'])
        self.result_side.SetMinSize(self.FromDIP(wx.Size(310, -1)))
        self.results_card = TableCard(self.result_side, 'Fit results', [('Result', 170), ('Value', 100)],
                                      empty='Press Fit')
        self.equation = wx.StaticText(self.result_side)
        self.equation.SetFont(ui_font(9))
        side = wx.BoxSizer(wx.VERTICAL)
        side.Add(self.results_card, 1, wx.EXPAND)
        side.Add(self.equation, 0, wx.EXPAND | wx.ALL, self.FromDIP(8))
        self.result_side.SetSizer(side)
        self.result_side.Bind(wx.EVT_SIZE, self.on_result_size)
        charts_layout = wx.BoxSizer(wx.VERTICAL)
        charts_layout.Add(self.card, 3, wx.EXPAND | wx.LEFT | wx.RIGHT, self.FromDIP(10))
        charts_layout.Add(self.residual_card, 2, wx.EXPAND | wx.ALL, self.FromDIP(10))
        self.charts.SetSizer(charts_layout)
        self.charts.SetAutoLayout(False)
        self.charts.Bind(wx.EVT_SIZE, self.on_chart_size)
        content = wx.BoxSizer(wx.HORIZONTAL)
        content.Add(self.charts, 1, wx.EXPAND)
        content.Add(self.result_side, 0, wx.EXPAND | wx.RIGHT | wx.BOTTOM, self.FromDIP(10))
        self.body.SetSizer(content)
        layout = wx.BoxSizer(wx.VERTICAL)
        layout.Add(bar, 0, wx.EXPAND | wx.ALL, self.FromDIP(10))
        layout.Add(self.body, 1, wx.EXPAND)
        self.SetSizer(layout)
        self.model.Bind(wx.EVT_CHOICE, self.on_settings)
        self.offset.Bind(wx.EVT_CHECKBOX, self.on_settings)
        self.unit.Bind(wx.EVT_TEXT, self.on_unit)
        self.on_settings()

    def on_chart_size(self, event):
        w, h = self.charts.GetClientSize()
        height = max(h, self.FromDIP(450))
        self.charts.SetVirtualSize(w, height)
        self.charts.GetSizer().SetDimension(0, 0, w, height)
        event.Skip()

    def on_settings(self, e=None):
        model = self.model.GetSelection()
        self.offset.Enable(model != 0)
        self.offset.SetLabel('Fit final level' if model == 1 else 'Fit starting level')
        self.equation.SetLabel(N.EQUATIONS[model])
        self.equation.Wrap(self.FromDIP(290))
        self.clear_result('Press Fit')

    def on_result_size(self, event):
        self.equation.SetLabel(N.EQUATIONS[self.model.GetSelection()])
        self.equation.Wrap(max(10, self.result_side.GetClientSize().width - self.FromDIP(16)))
        self.result_side.Layout()
        event.Skip()

    def clear_result(self, message):
        self.fit_result = self.fit_data = None
        self.export_button.Enable(False)
        self.summary.ChangeValue(message)
        self.results_card.set_rows([], message)
        self.card.reset()
        self.card.ylock = self.card._auto = None
        self.card.extra_scale = [(self.x, self.y)]
        self.card.ax.scatter(self.x, self.y, color=C['accent'], s=25, label='Measured')
        U.style_axes(self.card.ax, self.axis_label(), self.ylabel)
        self.residual_card.reset()
        self.residual_card.ylock = self.residual_card._auto = None
        U.style_axes(self.residual_card.ax, self.axis_label(), 'Difference')
        self.card.full = self.residual_card.full = (float(min(self.x)), float(max(self.x)))
        self.card.draw(); self.residual_card.draw()

    def on_unit(self, e=None):
        self.card.ax.set_xlabel(self.axis_label())
        self.residual_card.ax.set_xlabel(self.axis_label())
        self.card.draw(); self.residual_card.draw()
        if self.fit_result is not None:
            self.show_summary()

    def axis_label(self):
        unit = self.unit.GetValue().strip()
        return self.xlabel + (' (%s)' % unit if unit else '')

    def on_fit(self, e=None):
        self.clear_result('Fitting…')
        try:
            result, fitted, residuals = N.fit(self.x, self.y, self.model.GetSelection(), self.offset.GetValue())
            xx = np.linspace(float(min(self.x)), float(max(self.x)), 401)
            yy = N.predict(result, xx)
        except (ValueError, RuntimeError, OSError) as ex:
            self.summary.ChangeValue('Cannot fit: %s' % ex)
            self.results_card.set_rows([], 'Fit failed')
            wx.MessageBox(str(ex), 'Kinetic fitting', wx.OK | wx.ICON_INFORMATION, parent=self)
            return None
        self.fit_result = result
        self.fit_data = (fitted, residuals, xx, yy)
        self.card.ax.plot(xx, yy, color='#E56B32', lw=1.6, label='Fit')
        self.card.ax.legend(loc='best', fontsize=8)
        self.residual_card.ax.scatter(self.x, residuals, color=C['accent'], s=25)
        self.residual_card.extra_scale = [(self.x, residuals)]
        self.residual_card.ax.axhline(0, color=U.INK, lw=.7)
        self.card.ax.relim()
        self.card.ax.update_datalim(np.column_stack((self.x, self.y)))
        self.card.ax.autoscale_view()
        self.residual_card.ax.update_datalim(np.column_stack((self.x, residuals)))
        self.residual_card.ax.autoscale_view()
        self.show_summary()
        self.export_button.Enable(True)
        self.card.draw();self.residual_card.draw()
        return result

    def parameters(self):
        r = self.fit_result
        unit = self.unit.GetValue().strip() or 'X unit'
        name = ('Starting value', 'Final value', 'Starting value')[r.model]
        items = [(name, r.values[0], r.errors[0], self.ylabel)]
        if r.model == 0:
            items.append(('Slope', r.values[1], r.errors[1], self.ylabel + ' / ' + unit))
        else:
            items += [('Signal change', r.values[1], r.errors[1], self.ylabel),
                      ('Rate', r.values[2], r.errors[2], '1 / ' + unit),
                      ('Half-life' if r.model == 1 else 'Half-rise time', r.half_life, r.half_life_error, unit)]
        return items

    def show_summary(self):
        r = self.fit_result
        if r is None:
            return
        unit = self.unit.GetValue().strip()
        time_suffix = (' (%s)' % unit) if unit else ''
        rows = [('Model', self.model.GetStringSelection()), ('Data points', str(r.count)),
                ('First X value' + time_suffix, number(r.origin)),
                ('Fit quality (R²)', number(r.r2)), ('Fit error (RMSE)', number(r.rmse))]
        for name, value, error, _ in self.parameters():
            suffix = time_suffix if name in ('Half-life', 'Half-rise time') else (
                (' (1/%s)' % unit) if name == 'Rate' and unit else '')
            rows.append((name + suffix, number(value)))
            if math.isfinite(error):
                rows.append((name + ' uncertainty' + suffix, number(error)))
        if r.warnings & 1:
            rows.append(('Rate reliability', 'Poorly resolved'))
        if r.warnings & 2:
            rows.append(('Uncertainty', 'Unavailable'))
        if r.warnings & 4:
            rows.append(('Time points', 'More needed'))
        self.results_card.set_rows(rows)
        self.summary.ChangeValue('\n'.join('%s: %s' % row for row in rows))

    def on_export(self, e=None, path=None):
        if self.fit_result is None:
            self.summary.ChangeValue('Press Fit before exporting.')
            self.results_card.set_rows([], 'Press Fit')
            return None
        if not path:
            with wx.FileDialog(self, 'Export kinetic fit', defaultFile='kinetic_fit.xlsx',
                               wildcard='Excel workbook (*.xlsx)|*.xlsx', style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as dlg:
                if dlg.ShowModal() != wx.ID_OK:
                    return None
                path = dlg.GetPath()
                if not path.lower().endswith('.xlsx'):
                    path += '.xlsx'
        r = self.fit_result
        meta = self.metadata + [('Model', N.MODELS[r.model]), ('Equation', N.EQUATIONS[r.model]),
            ('X axis', self.xlabel), ('X unit', self.unit.GetValue().strip() or 'unspecified X unit'),
            ('Y', self.ylabel), ('First X value', r.origin), ('Data points', r.count),
            ('Fitted parameters', r.parameters), ('Residual degrees of freedom', r.dof),
            ('Fit baseline', r.model == 0 or self.offset.GetValue()),
            ('R squared', r.r2), ('RMSE', r.rmse), ('Residual sum of squares', r.sse)]
        for name, v, se, unit in self.parameters():
            meta += [(name, v), (name + ' uncertainty', se if math.isfinite(se) else 'unavailable'), (name + ' unit', unit)]
        meta += [('Fit warnings', ' '.join(N.warnings(r)) or 'none'),
                 ('Method', 'Least squares')]
        fitted, residuals, _, _ = self.fit_data
        rows = [[q['label'], q.get('file', ''), q['x'], self.y[i], fitted[i], residuals[i],
                 q['area'], q.get('relative') if q.get('relative') is not None else '', q.get('reference', False)]
                for i, q in enumerate(self.rows)]
        try:
            K.write_area_xlsx(path, ['Run', 'File', self.xlabel, 'Measured', 'Fitted', 'Difference',
                                    'Actual area', 'Relative area (%)', '100% reference'], rows, meta)
        except (OSError, ValueError) as ex:
            self.summary.ChangeValue('Export failed: %s' % ex)
            wx.MessageBox(str(ex), 'Export failed', wx.OK | wx.ICON_ERROR, parent=self)
            return None
        return path
