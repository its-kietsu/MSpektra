"""Cards of the polymer analysis in HRMS Analysis: a white card holding controls, an owner drawn table (check
boxes, coloured cells, grey rows, columns that fit the width) and the card of the distribution results."""
import wx
import unilcms as U
import unidec_theme as T
from unidec_theme import C, ui_font

WHITE = "#FFFFFF"
SELECTED = "#FCE6DA"  # (the selection colour of the plots, light)


class Card(wx.Panel):
    """A white rounded card like the plot and table tiles; children sit on white. title: drawn in its first
    row (header height top, DIP); no title: the children fill it."""
    own_refresh = True

    def __init__(self, parent, title="", top=0):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetBackgroundColour(WHITE)
        self.title, self.info, self.top = title, "", top
        self.Bind(wx.EVT_PAINT, self._paint)
        self.Bind(wx.EVT_SIZE, lambda e: (self.Layout(), U.refresh_card(self), e.Skip()))

    def set_info(self, info):
        if info != self.info:
            self.info = info
            self.Refresh()

    def _paint(self, e):
        U.paint_card(self, self.FromDIP(self.top or 14), U.card_frame, self._header,
                     (self.title, self.info, self.top, self.GetDPIScaleFactor()))

    def _header(self, gc, w, font=ui_font):
        if not self.title:
            return
        gc.SetFont(font(10, 700), wx.Colour(C["text"]))
        tw, th = gc.GetTextExtent(self.title)
        gc.DrawText(self.title, self.FromDIP(14), self.FromDIP(9))
        if self.info:
            gc.SetFont(font(9, 400), wx.Colour(C["muted"]))
            sw, sh = gc.GetTextExtent(self.info)
            gc.DrawText(self.info, self.FromDIP(14) + tw + self.FromDIP(10), self.FromDIP(9) + (th - sh) / 2.0)


def label(parent, text, muted=False, bold=False):
    t = wx.StaticText(parent, label=text)
    t.SetBackgroundColour(parent.GetBackgroundColour())
    if bold:
        t.SetFont(ui_font(10, 700))
    if muted:
        t.SetForegroundColour(wx.Colour(C["muted"]))
    return t


class Column(object):
    """key, title, align ('left'/'right'/'center'), least width (DIP), share of the room left, shown, optional
    (may be hidden by the user and is hidden first when the table is narrow); swatch: the colour of a cell is
    a small square before its text (series, as in the plots), not the colour of the text."""

    def __init__(self, key, title, align="right", width=60, grow=1.0, shown=True, optional=False, check=False,
                 swatch=False):
        self.key, self.title, self.align, self.width, self.grow = key, title, align, width, grow
        self.shown, self.optional, self.check, self.swatch = shown, optional, check, swatch


class _Rows(wx.VListBox):
    def __init__(self, table, multi):
        wx.VListBox.__init__(self, table, style=wx.BORDER_NONE | (wx.LB_MULTIPLE if multi else 0))
        self.table = table
        self.SetBackgroundColour(WHITE)
        self.SetFont(ui_font(9, 400))
        self.SetSelectionBackground(wx.Colour(SELECTED))

    def OnMeasureItem(self, n):
        return self.FromDIP(22)

    def OnDrawBackground(self, dc, rect, n):
        selected = self.IsSelected(n)
        dc.SetPen(wx.TRANSPARENT_PEN)
        dc.SetBrush(wx.Brush(wx.Colour(SELECTED if selected else WHITE)))
        dc.DrawRectangle(rect)

    def OnDrawItem(self, dc, rect, n):
        self.table.draw_row(self, dc, rect, n)


class Table(wx.Panel):
    """A table card: title and information line, then column titles and rows drawn here (a check box column,
    coloured cells, grey rows). rows: list of dict(cells=[text per column], colours={column: colour},
    muted=bool, check=True/False/None). Columns share the width; optional ones are left out when the table is
    narrow (or by the user, right click > Columns)."""
    own_refresh = True
    no_tile_menu = True

    def __init__(self, parent, title, columns, multi=False, empty=""):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetBackgroundColour(WHITE)
        self.title, self.info, self.empty = title, "", empty
        self.columns = list(columns)
        self.rows = []
        self.widths = []
        self.on_select = self.on_activate = self.on_check = self.on_cell = self.on_delete = None
        self.extra_menu = None
        self.head = wx.Panel(self, style=wx.BORDER_NONE)
        self.head.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.head.SetMinSize(wx.Size(-1, self.FromDIP(24)))
        self.body = _Rows(self, multi)
        s = wx.BoxSizer(wx.VERTICAL)
        s.AddSpacer(self.FromDIP(30))
        s.Add(self.head, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, self.FromDIP(10))
        s.Add(self.body, 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, self.FromDIP(10))
        self.SetSizer(s)
        self.Bind(wx.EVT_PAINT, self._paint)
        self.Bind(wx.EVT_SIZE, lambda e: (self.Layout(), self.fit_columns(), U.refresh_card(self), e.Skip()))
        self.head.Bind(wx.EVT_PAINT, self._paint_head)
        self.head.Bind(wx.EVT_RIGHT_DOWN, lambda e: self.show_menu())
        self.body.Bind(wx.EVT_LISTBOX, self._selected)
        self.body.Bind(wx.EVT_LISTBOX_DCLICK, lambda e: self.on_activate and self.on_activate(e.GetSelection()))
        self.body.Bind(wx.EVT_LEFT_DOWN, self._left)
        self.body.Bind(wx.EVT_CONTEXT_MENU, lambda e: self.show_menu())
        self.body.Bind(wx.EVT_KEY_DOWN, self._key)
        self.body.Bind(wx.EVT_MOUSEWHEEL, self._wheel)

    # ------------------------------------------------------------ data
    def set_columns(self, columns):
        self.columns = list(columns)
        self.fit_columns()

    def visible_columns(self):
        return [c for c, w in zip(self.columns, self.widths) if c.shown and w > 0]

    def fit_columns(self):
        """Widths: each shown column its least width, the room left by share; optional columns go when the
        least widths do not fit."""
        avail = max(0, self.body.GetClientSize()[0] - self.FromDIP(4))
        use = [c for c in self.columns if c.shown]
        while use and sum(self.FromDIP(c.width) for c in use) > avail and any(c.optional for c in use):
            drop = [c for c in use if c.optional][-1]
            use.remove(drop)
        least = sum(self.FromDIP(c.width) for c in use)
        grow = sum(c.grow for c in use) or 1.0
        extra = max(0, avail - least)
        scale = min(1.0, avail / float(least)) if least else 1.0  # (narrower than all least widths: texts shorten)
        self.widths = [(int(self.FromDIP(c.width) * scale) + int(extra * c.grow / grow)) if c in use else 0
                       for c in self.columns]
        self.head.Refresh()
        self.body.Refresh()

    def set_rows(self, rows, info=""):
        keep = self.selected()
        self.rows = list(rows)
        self.body.SetItemCount(len(self.rows))
        if self.body.HasMultipleSelection():
            self.body.DeselectAll()
        else:
            self.body.SetSelection(wx.NOT_FOUND)
        self.info = info
        self.fit_columns()
        self.Refresh()
        return keep

    def selected(self):
        if not self.body.HasMultipleSelection():
            i = self.body.GetSelection()
            return [i] if 0 <= i < len(self.rows) else []
        return [i for i in range(len(self.rows)) if self.body.IsSelected(i)]

    def select(self, i, only=True, show=True):
        if not 0 <= i < len(self.rows):
            return
        if self.body.HasMultipleSelection():
            if only:
                self.body.DeselectAll()
            self.body.Select(i, True)
        else:
            self.body.SetSelection(i)
        if show:
            self.body.ScrollToRow(i) if not self.body.IsVisible(i) else None
        self.body.Refresh()

    def text(self):
        cols = [k for k, c in enumerate(self.columns) if c.shown and not c.check]
        head = "\t".join(self.columns[k].title for k in cols)
        return "\n".join([head] + ["\t".join(str(r["cells"][k]) for k in cols) for r in self.rows])

    # ------------------------------------------------------------ drawing
    def _paint(self, e):
        U.paint_card(self, self.FromDIP(30), U.table_frame, self._draw_header,
                     (self.title, self.info, self.empty, bool(self.rows), self.GetDPIScaleFactor()))

    def _draw_header(self, gc, w, font=ui_font):
        gc.SetFont(font(10, 700), wx.Colour(C["text"]))
        tw, th = gc.GetTextExtent(self.title)
        gc.DrawText(self.title, self.FromDIP(14), self.FromDIP(9))
        gc.SetFont(font(9, 400), wx.Colour(C["muted"]))
        txt = self.info if self.rows else (self.empty or self.info)
        room = w - tw - self.FromDIP(40)
        if gc.GetTextExtent(txt)[0] > room:
            while txt and gc.GetTextExtent(txt + "…")[0] > room:
                txt = txt[:-1]
            txt = txt + "…" if txt else ""
        sw, sh = gc.GetTextExtent(txt)
        gc.DrawText(txt, self.FromDIP(14) + tw + self.FromDIP(10), self.FromDIP(9) + (th - sh) / 2.0)

    def _cells(self):
        x = self.FromDIP(2)
        for c, w in zip(self.columns, self.widths):
            if w > 0:
                yield c, x, w
                x += w

    def _text_at(self, dc, text, x, w, y, h, align):
        pad = self.FromDIP(6)
        tw, th = dc.GetTextExtent(text)
        if tw > w - 2 * pad:  # (shortened with an ellipsis, never cut)
            while text and dc.GetTextExtent(text + "…")[0] > w - 2 * pad:
                text = text[:-1]
            text += "…"
            tw = dc.GetTextExtent(text)[0]
        tx = x + pad if align == "left" else x + w - pad - tw if align == "right" else x + (w - tw) // 2
        dc.DrawText(text, int(tx), int(y + (h - th) // 2))

    def _paint_head(self, e):
        dc = wx.AutoBufferedPaintDC(self.head)
        dc.SetBackground(wx.Brush(wx.Colour(WHITE)))
        dc.Clear()
        w, h = self.head.GetClientSize()
        dc.SetFont(ui_font(9, 400))  # (as the column titles of the other tables)
        dc.SetTextForeground(wx.Colour(C["text"]))
        for c, x, cw in self._cells():
            if not c.check:
                self._text_at(dc, c.title, x, cw, 0, h, c.align)
        dc.SetPen(wx.Pen(wx.Colour(C["line"])))
        dc.DrawLine(0, h - 1, w, h - 1)

    def draw_row(self, win, dc, rect, n):
        if not 0 <= n < len(self.rows):
            return
        r = self.rows[n]
        dc.SetFont(ui_font(9, 400))
        muted = r.get("muted")
        for k, (c, x, w) in enumerate(self._cells()):
            i = self.columns.index(c)
            if c.check:
                if r.get("check") is None:
                    continue
                s = self.FromDIP(14)
                box = wx.Rect(rect.x + x + (w - s) // 2, rect.y + (rect.height - s) // 2, s, s)
                flags = wx.CONTROL_CHECKED if r.get("check") else 0
                if r.get("check_disabled"):
                    flags |= wx.CONTROL_DISABLED
                wx.RendererNative.Get().DrawCheckBox(win, dc, box, flags)
                continue
            colour = r.get("colours", {}).get(i)
            text = str(r["cells"][i])
            if c.swatch:
                if colour and text and not muted:
                    s = self.FromDIP(8)
                    dc.SetPen(wx.TRANSPARENT_PEN)
                    dc.SetBrush(wx.Brush(wx.Colour(colour)))
                    dc.DrawRectangle(rect.x + x + self.FromDIP(6), rect.y + (rect.height - s) // 2, s, s)
                    x, w = x + s + self.FromDIP(4), w - s - self.FromDIP(4)
                colour = None
            dc.SetTextForeground(wx.Colour("#8A93A0" if muted else colour or C["text"]))
            self._text_at(dc, text, rect.x + x, w, rect.y, rect.height, c.align)

    # ------------------------------------------------------------ input
    def column_at(self, x):
        for c, cx, w in self._cells():
            if cx <= x < cx + w:
                return c
        return None

    def _selected(self, e):
        if self.on_select:
            i = e.GetSelection()
            if 0 <= i < len(self.rows):
                self.on_select(i)

    def _left(self, e):
        n = self.body.VirtualHitTest(e.GetY())
        c = self.column_at(e.GetX())
        if n != wx.NOT_FOUND and c is not None and not (e.ControlDown() or e.ShiftDown()):
            if c.check and self.on_check and self.rows[n].get("check") is not None:
                self.on_check(n)
                return
            if self.on_cell and self.on_cell(n, c.key, e.GetPosition()):
                return
        e.Skip()

    def key_down(self, code, ctrl=False):
        """Delete (or Backspace) removes the selected rows; Ctrl+A selects all; Ctrl+C copies the table."""
        if code in (wx.WXK_DELETE, wx.WXK_BACK) and self.on_delete and self.selected():
            self.on_delete(self.selected())
            return True
        if ctrl and code == ord("A") and self.body.HasMultipleSelection():
            self.body.SelectAll()
            return True
        if ctrl and code == ord("C"):
            self.copy()
            return True
        return False

    def _key(self, e):
        if not self.key_down(e.GetKeyCode(), e.ControlDown()):
            e.Skip()

    def _wheel(self, e):
        # the rows scroll while they can, then the column of tiles does
        first = self.body.GetVisibleRowsBegin()
        last = self.body.GetVisibleRowsEnd()
        down = e.GetWheelRotation() < 0
        if (down and last >= len(self.rows)) or (not down and first <= 0):
            U.forward_wheel(self, e)
        else:
            e.Skip()

    def copy(self):
        U._clip(self.text())
        U._status(self, "Table copied")

    def menu_items(self):
        items = list(self.extra_menu()) if self.extra_menu else []
        if self.rows:
            items += [(None, None), ("Copy table", self.copy)]
        optional = [c for c in self.columns if c.optional]
        if optional:
            def toggle(c):
                c.shown = not c.shown
                self.fit_columns()
            items.append(("Columns", [(c.title, lambda c=c: toggle(c), c.shown) for c in optional]))
        return items

    def show_menu(self):
        items = self.menu_items()
        if not items:
            return
        m = U.build_menu(self, items)
        U.popup_menu(self, m)
        m.Destroy()


class ResultCard(Table):
    """The results of the mass distribution: a Measure / Value table like the other result tables."""

    def __init__(self, parent):
        Table.__init__(self, parent, "Result", [Column("measure", "Measure", "left", 80, 1.0),
                                                Column("value", "Value", "right", 70, 1.0)])
        self.values = []  # [(measure, value)]

    def set_values(self, values, info="", note=""):
        self.values = list(values)
        self.empty = note
        self.set_rows([dict(cells=[name, text]) for name, text in self.values], info)
