"""Cards of the polymer analysis in HRMS Analysis: a white card holding controls, an owner drawn table (check
boxes, coloured cells, grey rows, columns that fit the width) and the card of the distribution results."""
import wx
import ribbon as RB
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

    def __init__(self, parent, title, columns, multi=False, empty="", flat=False, note=False):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetBackgroundColour(WHITE)
        self.title, self.info, self.empty = title, "", empty
        self.flat = flat  # (inside a tile: no title row and no edge; the tile shows title and information)
        self.on_info = None  # callback(text) with the information line (flat)
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
        s.AddSpacer(self.FromDIP(2 if flat else 30))
        side = 0 if flat else self.FromDIP(10)
        s.Add(self.head, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, side)
        s.Add(self.body, 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, side)
        self.note = None  # (note: the empty text where the rows go)
        if note:
            self.note = wx.StaticText(self, label=empty)
            self.note.SetBackgroundColour(WHITE)
            self.note.SetForegroundColour(wx.Colour(C["muted"]))
            s.Add(self.note, 0, wx.LEFT | wx.TOP, self.FromDIP(6))
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
        if self.note is not None:
            self.note.SetLabel(self.empty)
            self.note.Show(not self.rows)
            self.body.Show(bool(self.rows))
            self.Layout()
        self.fit_columns()
        self.Refresh()
        if self.on_info:
            self.on_info(self.info_text())
        return keep

    def info_text(self):
        return self.info if self.rows else (self.empty or self.info)

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
        if self.flat:
            U.paint_card(self, self.FromDIP(14), U.inner_frame, lambda gc, w, font=None: None, ("flat",))
            return
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


class Toggle(RB.RibbonButton):
    """A small button as the views of the ribbon (MS, PDA, Compare): the chosen one filled in blue. icon: a
    ribbon icon name or None (the label alone)."""

    def __init__(self, parent, label, on_click, icon=None, active=False, tip=None):
        state = dict(label=label, active=active)
        item = RB.Item(label, lambda: state["label"], icon, lambda b: on_click(), size="small",
                       tip=lambda: tip or state["label"].split(" · ")[0], active=lambda: state["active"])
        RB.RibbonButton.__init__(self, parent, item, style="view")
        self._state = state
        self.mode = "small"
        self.SetFont(ui_font(9))
        self.SetMinSize(self.best("small"))

    @property
    def active(self):
        return self._state["active"]

    def set_active(self, on):
        if bool(on) != self._state["active"]:
            self._state["active"] = bool(on)
            self.read_state()

    def set_label(self, text):
        """True when its size changed."""
        if text != self._state["label"]:
            self._state["label"] = text
            self.read_state()
            size = self.best("small")
            if size != self.GetMinSize():
                self.SetMinSize(size)
                self.SetSize(size)
                return True
        return False

    def best(self, mode):
        if mode != "small" or self.item.icon:
            return RB.RibbonButton.best(self, mode)
        dc = wx.ClientDC(self)
        dc.SetFont(ui_font(9, 600))
        return wx.Size(dc.GetTextExtent(self._text())[0] + self.FromDIP(16), self.FromDIP(RB.SMALL_H + 2))

    def _paint(self, e):
        if self.item.icon:
            return RB.RibbonButton._paint(self, e)
        dc = wx.AutoBufferedPaintDC(self)
        dc.SetBackground(wx.Brush(self.GetParent().GetBackgroundColour()))
        dc.Clear()
        gc = T.crisp(dc)
        w, h = self.GetClientSize()
        label, _, enabled, active = self.state
        r = 6 * self.FromDIP(10) / 10.0
        if active:
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.SetBrush(wx.Brush(wx.Colour(C["accent_hover"] if self.hover else C["accent"])))
        else:
            gc.SetBrush(wx.Brush(wx.Colour(C["hover"] if (self.hover and enabled) else C["panel"])))
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(C["line2"])).Width(1.0)))
        gc.DrawRoundedRectangle(0.5, 0.5, w - 1, h - 1, r)
        colour = RB.GREY if not enabled else "#FFFFFF" if active else RB.STROKE
        gc.SetFont(ui_font(9, 600 if active else 400), wx.Colour(colour))
        text = self._text()
        tw, th = gc.GetTextExtent(text)
        gc.DrawText(text, (w - tw) / 2.0, (h - th) / 2.0)


class KeyChoice(object):
    """A choice (self.choice, a wx.Choice) whose items are kept by key (the texts saved, such as "+1") and shown by
    label ("+1 · 39 masses"), in an order of their own (indices of keys; keys left out are not offered). Used
    like a choice, by key index."""

    def __init__(self, parent, on_change=None):
        self.choice = wx.Choice(parent)
        self.on_change = on_change
        self.keys, self.labels, self.order = [], [], []
        self.choice.Bind(wx.EVT_CHOICE, self._chosen)

    def _chosen(self, e=None):
        if self.on_change:
            self.on_change()

    def GetParent(self):
        return self.choice.GetParent()

    def SetItems(self, keys, labels=None, order=None):
        old = self.GetStringSelection()
        self.keys = list(keys)
        self.labels = list(labels) if labels is not None else list(self.keys)
        self.order = list(order) if order is not None else list(range(len(self.keys)))
        self.choice.SetItems([self.labels[i] for i in self.order])
        self.choice.InvalidateBestSize()
        if not self.SetStringSelection(old) and self.order:
            self.choice.SetSelection(0)

    def GetSelection(self):
        k = self.choice.GetSelection()
        return self.order[k] if 0 <= k < len(self.order) else -1

    def SetSelection(self, i):
        self.choice.SetSelection(self.order.index(i) if i in self.order else wx.NOT_FOUND)

    def GetStringSelection(self):
        i = self.GetSelection()
        return self.keys[i] if 0 <= i < len(self.keys) else ""

    def SetStringSelection(self, key):
        if key in self.keys and self.keys.index(key) in self.order:
            self.SetSelection(self.keys.index(key))
            return True
        return False

    def GetStrings(self):
        return list(self.keys)

    def GetCount(self):
        return len(self.keys)

    def shown(self):
        """The labels offered, in their order."""
        return [self.labels[i] for i in self.order]


class Line(wx.Panel):
    """A thin line (orient: wx.HORIZONTAL or wx.VERTICAL)."""

    def __init__(self, parent, orient=wx.HORIZONTAL):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundColour(C["line"])
        self.SetMinSize(wx.Size(-1, 1) if orient == wx.HORIZONTAL else wx.Size(1, -1))


class Flow(wx.Panel):
    """Panels from left to right, a thin line between them; when the row is full the next panel starts a new
    line. height_for(width): the height at a width."""

    def __init__(self, parent, gap=14, vgap=8, lines=True):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundColour(parent.GetBackgroundColour())
        self.gap, self.vgap, self.lines = self.FromDIP(gap), self.FromDIP(vgap), lines
        self.items = []
        self.right = None  # a panel at the right end of the first line
        self.on_height = None  # callback() when the height needed changed
        self._lines = []
        self.Bind(wx.EVT_SIZE, lambda e: (self.place(), e.Skip()))
        self.Bind(wx.EVT_PAINT, self._paint)

    def set_items(self, items):
        self.items = list(items)
        self.place()

    def _arrange(self, width):
        """[[item, x, y, w, h, (line top, line height)]], height."""
        rw, rh = self.right.GetBestSize() if self.right is not None else (0, 0)
        full = width
        if self.right is not None:
            width = max(1, width - rw - 2 * self.gap - 1)
        out, x, y, line_h, first = [], 0, 0, 0, True
        for it in self.items:
            if not it.IsShown():
                continue
            w, h = it.GetBestSize()
            if not first and x + 2 * self.gap + 1 + w > width:
                y += line_h + self.vgap
                x, line_h, first = 0, 0, True
            if not first:
                x += 2 * self.gap + 1
            out.append([it, x, y, w, h])
            x += w
            line_h = max(line_h, h)
            first = False
        if self.right is not None:
            out.append([self.right, full - rw, 0, rw, rh])
            line_h = max(line_h, rh) if y == 0 else line_h
        lines = {}
        for o in out:
            lines.setdefault(o[2], []).append(o)
        for y0, line in lines.items():
            lh = max(o[4] for o in line)
            for o in line:  # (each in the middle of its line)
                o[2] = y0 + (lh - o[4]) // 2
                o.append((y0, lh))
        return out, (max(o[5][0] + o[5][1] for o in out) if out else 0)

    def height_for(self, width):
        return self._arrange(width)[1]

    def items_changed(self):
        self.place()

    def place(self):
        w = self.GetClientSize()[0]
        if w <= 0:
            return
        out, h = self._arrange(w)
        for it, x, y, iw, ih, line in out:
            it.SetSize(x, y, iw, ih)
        self._lines = [(x - self.gap - 1, line) for (it, x, y, iw, ih, line) in out
                       if x > 0 and self.lines and it is not self.right]
        self.Refresh()
        if self.GetMinSize()[1] != h:
            self.SetMinSize(wx.Size(-1, h))
            if self.on_height:
                wx.CallAfter(self.on_height)

    def _paint(self, e):
        dc = wx.PaintDC(self)
        dc.SetPen(wx.Pen(wx.Colour(C["line"]), 1))
        for x, (y0, lh) in self._lines:
            dc.DrawLine(x, y0, x, y0 + lh)


class Summary(wx.StaticText):
    """The result of the mass distribution in one line (Mn, Mw, Đ); values: [(measure, value)]."""

    SHOWN = ("Mn", "Mw", "Đ")

    def __init__(self, parent):
        wx.StaticText.__init__(self, parent, label="")
        self.SetBackgroundColour(parent.GetBackgroundColour())
        self.values = []

    def set_values(self, values, info="", note=""):
        self.values = list(values)
        text = " · ".join("%s %s" % (name, value.replace(" Da", "")) for name, value in self.values
                          if name in self.SHOWN)
        self.SetForegroundColour(wx.Colour(C["text"] if text else C["muted"]))
        self.SetLabel(text or note)
        self.SetToolTip(info or "")
        self.InvalidateBestSize()
