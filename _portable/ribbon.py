"""The ribbon of LCMS and HRMS Analysis: a row of tabs and, under it, the commands of the chosen tab in named
groups (large buttons; small or icon only buttons when the window is narrow). The views of the data (Mass
spectrometry, PDA, Compare ...) are large buttons pinned at the start of the ribbon, on every tab.

The ribbon only shows commands: each Item calls code that exists elsewhere (the hidden tool strips of the
views, the hidden top bar buttons, menus), and its state (shown, enabled, active, label) is read again a few
times a second, so every other way to use a command (shortcuts, menus) keeps the ribbon right."""
import time
import wx

import unidec_theme as T
from unidec_theme import C, ui_font

TAB_H, MENU_H, BODY_H, CAP_H = 34, 30, 96, 18
LARGE_H, SMALL_H = 68, 22
ICON_L, ICON_S = 26, 16
GREY = "#A9B1BC"
STROKE = "#3B4450"

ICONS = {
    "ms": '<path d="M4 20v-8M8 20V5M12 20v-9M16 20V8M20 20v-3"/>',
    "pda": '<path d="M3 16c3-9 6-9 9 0s6 9 9 0"/>',
    "compare": '<path d="M3 19c3 0 4-10 6-10s3 10 6 10h6M3 14.5c3 0 4-8 6-8s3 8 6 8h6"/>',
    "save": '<path d="M5 4h11l3 3v13H5z"/><path d="M8 4v5h7V4M8 20v-6h8v6"/>',
    "export": '<path d="M12 3v12M7 10l5 5 5-5M4 21h16"/>',
    "save_as": '<path d="M5 4h11l3 3v5"/><path d="M5 4v16h7"/><path d="M8 4v5h7V4"/>'
               '<path d="M14.5 20.5l1-3.5 4.5-4.5 2.5 2.5-4.5 4.5z"/>',
    "copy": '<rect x="8" y="8" width="12" height="12" rx="1.5"/><path d="M4 16V5.5A1.5 1.5 0 0 1 5.5 4H16"/>',
    "image": '<rect x="3" y="5" width="18" height="14" rx="1.5"/><path d="M3 16l5-5 4 4 3-3 6 6"/>'
             '<circle cx="15.5" cy="9" r="1.3"/>',
    "method": '<path d="M4 6h16M4 12h10M4 18h13"/>',
    "help": '<circle cx="12" cy="12" r="9"/><path d="M9.5 9.5a2.5 2.5 0 1 1 3.5 2.3c-.7.3-1 .8-1 1.5M12 17v.3"/>',
    "table": '<rect x="3" y="5" width="18" height="14" rx="1.5"/><path d="M3 10h18M9 5v14"/>',
    "settings": '<circle cx="12" cy="12" r="3"/><path d="M12 2.5v3M12 18.5v3M2.5 12h3M18.5 12h3M5.3 5.3l2.1 2.1'
                'M16.6 16.6l2.1 2.1M5.3 18.7l2.1-2.1M16.6 7.4l2.1-2.1"/>',
    "folder": '<path d="M3 7h6l2 2h10v10H3z"/>',
    "list": '<path d="M8 6h13M8 12h13M8 18h13M3.5 6h.01M3.5 12h.01M3.5 18h.01"/>',
    "prev": '<path d="M15 5l-7 7 7 7"/>',
    "next": '<path d="M9 5l7 7-7 7"/>',
}


def icon_body(name):
    """The SVG paths of an icon: ribbon icons, the tool icons of the views, the theme icons."""
    if name in ICONS:
        return ICONS[name]
    try:
        import unilcms
        if name in unilcms.TOOL_ICONS:
            return unilcms.TOOL_ICONS[name]
    except Exception:
        pass
    return T.ICONS.get(name, "")


def _bitmap(name, size, colour, width):
    body = icon_body(name)
    if not body:
        return None
    try:
        from wx.svg import SVGimage
        svg = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="%s" '
               'stroke-width="%.2f" stroke-linecap="round" stroke-linejoin="round">%s</svg>'
               % (colour, width, body.replace("{c}", colour)))
        return SVGimage.CreateFromBytes(svg.encode("utf-8")).ConvertToScaledBitmap(wx.Size(size, size))
    except Exception:
        return None


def _call(f, default=None):
    if f is None:
        return default
    if not callable(f):
        return f
    try:
        return f()
    except Exception:
        return default


class Item:
    """One command. run(button) is called on a click; label, tip, shown, enabled and active may be values or
    functions (read again while the window is open). size: "large" or "small" (the size when there is room).
    menu=True: the label gets a small arrow (run opens a menu at the button)."""

    def __init__(self, key, label, icon, run, size="large", tip=None, shown=True, enabled=True, active=False,
                 menu=False):
        self.key, self._label, self.icon, self.run, self.size = key, label, icon, run, size
        self._tip, self._shown, self._enabled, self._active, self.menu = tip, shown, enabled, active, menu

    def label(self):
        return _call(self._label, "") or ""

    def tip(self):
        return _call(self._tip, "") or self.label()

    def shown(self):
        return bool(_call(self._shown, True))

    def enabled(self):
        return bool(_call(self._enabled, True))

    def active(self):
        return bool(_call(self._active, False))


class Group:
    def __init__(self, title, items):
        self.title, self.items = title, items


class Tab:
    def __init__(self, name, groups):
        self.name, self.groups = name, groups


# ---------------------------------------------------------------------------------------------- buttons
class RibbonButton(wx.Control):
    """A large (icon over label), small (icon and label) or icon only button; style "view": the large toggles
    of the views (the active one filled in blue)."""

    def __init__(self, parent, item, style="normal"):
        wx.Control.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.item, self.style = item, style
        self.mode = "large"
        self.hover = self.down = False
        self.state = None
        self._bmps = {}
        for ev, fn in ((wx.EVT_PAINT, self._paint), (wx.EVT_ENTER_WINDOW, lambda e: self._h(True)),
                       (wx.EVT_LEAVE_WINDOW, lambda e: self._h(False)), (wx.EVT_LEFT_DOWN, self._down),
                       (wx.EVT_LEFT_UP, self._up), (wx.EVT_LEFT_DCLICK, self._down)):
            self.Bind(ev, fn)
        self.read_state()

    def read_state(self):
        """(label, tip, enabled, active) now; True when it changed."""
        it = self.item
        st = (it.label(), it.tip(), it.enabled(), it.active())
        if st != self.state:
            self.state = st
            if self.GetToolTipText() != st[1]:
                self.SetToolTip(st[1])
            self.Refresh()
            return True
        return False

    def _h(self, v):
        self.hover = v
        if not v:
            self.down = False
        self.Refresh()

    def _down(self, e):
        if self.state[2]:
            self.down = True
            self.Refresh()

    def _up(self, e):
        was = self.down
        self.down = False
        self.Refresh()
        if was and self.state[2]:
            wx.CallAfter(self._run)

    def _run(self):
        try:
            self.item.run(self)
        except RuntimeError:
            pass
        r = self.GetTopLevelParent().__dict__.get("ribbon")
        if r:
            wx.CallLater(60, r.read_state)

    def _text(self):
        lab = self.state[0]
        return lab + (" ▾" if self.item.menu and lab else "")

    def best(self, mode):
        """The size of the button in a mode ("large", "small", "icon")."""
        dc = wx.ClientDC(self)
        if mode == "large":
            dc.SetFont(ui_font(9 if self.style == "normal" else 9.5, 400 if self.style == "normal" else 600))
            tw = dc.GetTextExtent(self._text())[0]  # (a view button: as wide as its label in bold)
            w = max(self.FromDIP(50 if self.style == "normal" else 88), tw + self.FromDIP(14 if self.style == "normal" else 24))
            return wx.Size(w, self.FromDIP(LARGE_H if self.style == "normal" else 74))
        if mode == "small":
            dc.SetFont(ui_font(9, 400))
            tw = dc.GetTextExtent(self._text())[0]
            return wx.Size(self.FromDIP(6 + ICON_S + 6 + 8) + tw, self.FromDIP(SMALL_H))
        return wx.Size(self.FromDIP(28), self.FromDIP(SMALL_H + 2))

    def _bmp(self, size, colour, width):
        k = (size, colour, width)
        if k not in self._bmps:
            self._bmps[k] = _bitmap(self.item.icon, size, colour, width)
        return self._bmps[k]

    def _paint(self, e):
        dc = wx.AutoBufferedPaintDC(self)
        dc.SetBackground(wx.Brush(self.GetParent().GetBackgroundColour()))
        dc.Clear()
        gc = T.crisp(dc)
        w, h = self.GetClientSize()
        k = self.FromDIP(10) / 10.0
        label, _, enabled, active = self.state
        view = self.style == "view"
        r = 6 * k if view else 4 * k
        if view and active:
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.SetBrush(wx.Brush(wx.Colour(C["accent_hover"] if self.hover else C["accent"])))
            gc.DrawRoundedRectangle(0.5, 0.5, w - 1, h - 1, r)
        elif view:
            gc.SetBrush(wx.Brush(wx.Colour(C["hover"] if (self.hover and enabled) else C["panel"])))
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(C["line2"])).Width(1.0)))
            gc.DrawRoundedRectangle(0.5, 0.5, w - 1, h - 1, r)
        elif active:
            gc.SetBrush(wx.Brush(wx.Colour(C["accent_bg"])))
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(59, 130, 217, 150)).Width(1.0)))
            gc.DrawRoundedRectangle(0.5, 0.5, w - 1, h - 1, r)
        elif enabled and (self.down or self.hover):
            gc.SetBrush(wx.Brush(wx.Colour(C["down"] if self.down else C["hover"])))
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.DrawRoundedRectangle(0.5, 0.5, w - 1, h - 1, r)
        if not enabled:
            col = GREY
        elif view:
            col = "#FFFFFF" if active else STROKE
        else:
            col = C["accent_text"] if (active or self.hover) else STROKE
        tcol = col if (not enabled or view or active) else (C["accent_text"] if self.hover else C["text"])
        text = self._text()
        if self.mode == "large":
            s = self.FromDIP(ICON_L)
            bmp = self._bmp(s, col, 1.55)
            top = self.FromDIP(9 if view else 7)
            if bmp:
                gc.DrawBitmap(bmp, (w - s) / 2.0, top, s, s)
            gc.SetFont(ui_font(9.5 if view else 9, 600 if (view and active) else 400), wx.Colour(tcol))
            tw, th = gc.GetTextExtent(text)
            gc.DrawText(text, (w - tw) / 2.0, top + s + self.FromDIP(6))
        else:
            s = self.FromDIP(ICON_S)
            bmp = self._bmp(s, col, 1.8)
            x = self.FromDIP(6) if self.mode == "small" else (w - s) / 2.0
            if bmp:
                gc.DrawBitmap(bmp, x, (h - s) / 2.0, s, s)
            if self.mode == "small" and text:
                gc.SetFont(ui_font(9, 400), wx.Colour(tcol))
                tw, th = gc.GetTextExtent(text)
                gc.DrawText(text, x + s + self.FromDIP(6), (h - th) / 2.0)


# ---------------------------------------------------------------------------------------------- groups
class GroupPanel(wx.Panel):
    """The buttons of a group and its title under them; level 0: the sizes asked for, 1: all small, 2: icons."""

    def __init__(self, parent, group, style="normal"):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        self.group, self.style = group, style
        self.level = 0
        self.buttons = [RibbonButton(self, it, style) for it in group.items]
        self.shown = None
        self.Bind(wx.EVT_PAINT, self._paint)

    def visible(self):
        return [b for b in self.buttons if b.item.shown()]

    def any_shown(self):
        return any(b.item.shown() for b in self.buttons)

    def _columns(self, level):
        """[(mode, [buttons])] in the order of the items: large buttons alone, small ones three per column."""
        cols = []
        for b in self.visible():
            mode = b.item.size if level == 0 else ("small" if level == 1 else "icon")
            if self.style == "view":
                mode = "large"
            if mode == "large":
                cols.append(("large", [b]))
            elif cols and cols[-1][0] == mode and len(cols[-1][1]) < 3:
                cols[-1][1].append(b)
            else:
                cols.append((mode, [b]))
        return cols

    def width(self, level):
        pad, gap = self.FromDIP(6), self.FromDIP(2 if self.style == "normal" else 6)
        cols = self._columns(level)
        w = sum(max(b.best(m).width for b in bs) for m, bs in cols) + gap * max(0, len(cols) - 1) + 2 * pad
        dc = wx.ClientDC(self)
        dc.SetFont(ui_font(8, 400))
        return max(w, dc.GetTextExtent(self.group.title)[0] + 2 * pad) + self.FromDIP(1)

    def place(self, level):
        self.level = level
        pad, gap = self.FromDIP(6), self.FromDIP(2 if self.style == "normal" else 6)
        area = self.FromDIP(BODY_H - CAP_H - 6)
        x = pad
        vis = set(self.visible())
        for b in self.buttons:
            b.Show(b in vis)
        for mode, bs in self._columns(level):
            cw = max(b.best(mode).width for b in bs)
            if mode == "large":
                sz = bs[0].best("large")
                bs[0].mode = "large"
                bs[0].SetSize(x, self.FromDIP(4), cw, sz.height)
            else:
                hh = self.FromDIP(SMALL_H + (2 if mode == "icon" else 0))
                y0 = self.FromDIP(4) + max(0, (area - 3 * hh) // 2) if len(bs) == 3 else self.FromDIP(4) + (area - len(bs) * hh) // 2
                for i, b in enumerate(bs):
                    b.mode = mode
                    b.SetSize(x, y0 + i * hh, cw if mode == "small" else b.best("icon").width, hh)
            x += cw + gap
        for b in self.buttons:
            b.Refresh()
        self.Refresh()

    def _paint(self, e):
        dc = wx.AutoBufferedPaintDC(self)
        dc.SetBackground(wx.Brush(wx.Colour(C["panel"])))
        dc.Clear()
        gc = T.crisp(dc)
        w, h = self.GetClientSize()
        gc.SetFont(ui_font(8, 400), wx.Colour(C["faint"]))
        tw, th = gc.GetTextExtent(self.group.title)
        gc.DrawText(self.group.title, (w - tw) / 2.0, h - self.FromDIP(CAP_H) + (self.FromDIP(CAP_H) - th) / 2.0 - 1)
        line = C["line2"] if self.style == "view" else C["line"]
        gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(line)).Width(2.0 if self.style == "view" else 1.0)))
        gc.StrokeLine(w - 1, self.FromDIP(6), w - 1, h - self.FromDIP(6))


# ---------------------------------------------------------------------------------------------- tab row
class TabRow(wx.Panel):
    """The tabs (painted) and, at the right, other controls (chips, help) added with add()."""

    def __init__(self, parent, ribbon):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        self.ribbon = ribbon
        self.hover = -1
        self.rects = []
        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.tabs_space = self.sizer.Add(self.FromDIP(10), self.FromDIP(TAB_H))
        self.sizer.AddStretchSpacer(1)
        self.SetSizer(self.sizer)
        self.Bind(wx.EVT_PAINT, self._paint)
        self.Bind(wx.EVT_LEFT_DOWN, self._click)
        self.Bind(wx.EVT_MOTION, self._motion)
        self.Bind(wx.EVT_LEAVE_WINDOW, lambda e: self._set_hover(-1))
        # moved to a display of another scale: the titles measured again (their places were kept in pixels of
        # the old scale: wide gaps between the menus after a move from 200 % to 100 %)
        T.keep_dip(self, then=lambda: (self.measure(), self.Refresh()))

    # the interface of ModernToolbar used by the chips (method_presets.add_chip)
    def add(self, w, gap=4):
        self.sizer.Add(w, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, self.FromDIP(gap))
        self.Layout()

    def separator(self):
        pass

    def stretch(self):
        pass

    def titles(self):
        """The menus of the menu bar, else the tabs."""
        if self.ribbon.menus:
            return [m[0] for m in self.ribbon.menus]
        return [t.name for t in self.ribbon.tabs]

    def measure(self):
        dc = wx.ClientDC(self)
        menus = bool(self.ribbon.menus)
        dc.SetFont(ui_font(9.5, 400 if menus else 600))
        x, self.rects = self.FromDIP(4 if menus else 8), []
        h = self.FromDIP(MENU_H if menus else TAB_H)
        for t in self.titles():
            tw = dc.GetTextExtent(t)[0] + self.FromDIP(20 if menus else 28)
            self.rects.append(wx.Rect(x, 0, tw, h))
            x += tw
        self.tabs_space.AssignSpacer(wx.Size(x + self.FromDIP(8), h))
        self.Layout()

    def _at(self, pos):
        for i, r in enumerate(self.rects):
            if r.Contains(pos):
                return i
        return -1

    def _set_hover(self, i):
        if i != self.hover:
            self.hover = i
            self.Refresh()

    def _motion(self, e):
        self._set_hover(self._at(e.GetPosition()))

    def _click(self, e):
        i = self._at(e.GetPosition())
        if i < 0:
            return
        if self.ribbon.menus:
            # a click on the title of the menu that is open closes it: Windows closes the menu and then
            # still hands the click to this row, which opened the menu again at once
            ci, ct = self.ribbon.closed
            if i == self.ribbon.open_i or (i == ci and time.monotonic() - ct < 0.3):
                return
            wx.CallAfter(self.ribbon.open_menu, i)
        else:
            self.ribbon.select(i)

    def _paint(self, e):
        dc = wx.AutoBufferedPaintDC(self)
        dc.SetBackground(wx.Brush(wx.Colour(C["panel"])))
        dc.Clear()
        gc = T.crisp(dc)
        w, h = self.GetClientSize()
        if self.ribbon.menus:  # a menu bar: the titles, the one open or under the mouse on a light band
            for i, (t, r) in enumerate(zip(self.titles(), self.rects)):
                if i in (self.hover, self.ribbon.open_i):
                    gc.SetPen(wx.TRANSPARENT_PEN)
                    gc.SetBrush(wx.Brush(wx.Colour(C["down"] if i == self.ribbon.open_i else C["hover"])))
                    gc.DrawRoundedRectangle(r.x + 1, r.y + self.FromDIP(3), r.width - 2, r.height - self.FromDIP(6),
                                            self.FromDIP(4))
                gc.SetFont(ui_font(9.5, 400), wx.Colour(C["text"]))
                tw, th = gc.GetTextExtent(t)
                gc.DrawText(t, r.x + (r.width - tw) / 2.0, (h - th) / 2.0)
            return
        for i, (t, r) in enumerate(zip(self.ribbon.tabs, self.rects)):
            on = i == self.ribbon.current
            if self.hover == i and not on:
                gc.SetPen(wx.TRANSPARENT_PEN)
                gc.SetBrush(wx.Brush(wx.Colour(C["hover"])))
                gc.DrawRoundedRectangle(r.x + 2, r.y + self.FromDIP(4), r.width - 4, r.height - self.FromDIP(8),
                                        self.FromDIP(4))
            gc.SetFont(ui_font(9.5, 600 if on else 400), wx.Colour(C["accent_text"] if on else C["muted"]))
            tw, th = gc.GetTextExtent(t.name)
            gc.DrawText(t.name, r.x + (r.width - tw) / 2.0, (h - th) / 2.0)
            if on:
                gc.SetPen(wx.TRANSPARENT_PEN)
                gc.SetBrush(wx.Brush(wx.Colour(C["accent"])))
                bw = tw + self.FromDIP(8)
                gc.DrawRoundedRectangle(r.x + (r.width - bw) / 2.0, h - self.FromDIP(3), bw, self.FromDIP(3), 1.5)


# ---------------------------------------------------------------------------------------------- ribbon
class Ribbon(wx.Panel):
    """menus: [(title, entries())] for a menu bar above one toolbar (the first tab) instead of tabs; entries()
    gives the menu when it is opened (unilcms.build_menu entries)."""

    def __init__(self, parent, tabs, pinned=None, menus=None):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        self.tabs, self.current = tabs, 0
        self.menus, self.open_i = menus, -1
        self.closed = (-1, 0.0)   # the menu closed last and when (time.monotonic)
        self.tabrow = TabRow(self, self)
        self.body = wx.Panel(self, style=wx.BORDER_NONE)
        self.body.SetBackgroundColour(wx.Colour(C["panel"]))
        self.body.SetMinSize(wx.Size(-1, self.FromDIP(BODY_H)))
        self.pinned = GroupPanel(self.body, pinned, "view") if pinned else None
        self.pages = []
        for t in tabs:
            p = wx.Panel(self.body, style=wx.BORDER_NONE)
            p.SetBackgroundColour(wx.Colour(C["panel"]))
            p.groups = [GroupPanel(p, g) for g in t.groups]
            p.Show(False)
            self.pages.append(p)
        line = wx.Panel(self, size=wx.Size(-1, 1))
        line.SetBackgroundColour(wx.Colour(C["line2"]))
        top = wx.Panel(self, size=wx.Size(-1, 1))
        top.SetBackgroundColour(wx.Colour(C["line"]))
        s = wx.BoxSizer(wx.VERTICAL)
        s.Add(self.tabrow, 0, wx.EXPAND)
        s.Add(top, 0, wx.EXPAND)
        s.Add(self.body, 0, wx.EXPAND)
        s.Add(line, 0, wx.EXPAND)
        self.SetSizer(s)
        self.tabrow.measure()
        self.body.Bind(wx.EVT_SIZE, lambda e: (self.fit(), e.Skip()))
        self.select(0)
        self._timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, lambda e: self.read_state(), self._timer)
        self._timer.Start(300)
        self.Bind(wx.EVT_WINDOW_DESTROY, self._gone)
        T.keep_dip(self, then=self.fit)  # another scale: the groups placed again

    def _gone(self, e):
        if e.GetEventObject() is self:
            self._timer.Stop()
        e.Skip()

    def open_menu(self, i):
        import unilcms
        fr = self.GetTopLevelParent()
        try:
            entries = self.menus[i][1]()
        except Exception as ex:
            print("menu %s: %s" % (self.menus[i][0], ex))
            return
        m = unilcms.build_menu(fr, entries)
        r = self.tabrow.rects[i]
        self.open_i = i
        self.tabrow.Refresh()
        self.tabrow.Update()
        try:
            unilcms.popup_menu(self.tabrow, m, wx.Point(r.x, r.y + r.height))
        finally:
            m.Destroy()
            self.open_i = -1
            self.closed = (i, time.monotonic())
            if self:
                self.tabrow.Refresh()
                wx.CallLater(60, self.read_state)

    def select(self, i):
        self.current = i
        for j, p in enumerate(self.pages):
            p.Show(j == i)
        self.tabrow.Refresh()
        self.fit()

    def read_state(self):
        """Labels, tooltips, enabled and active state of the buttons shown, and which buttons are shown (the
        view shown changes the tools): drawn again when changed."""
        if not self or not self.IsShownOnScreen():
            return
        groups = list(self.pages[self.current].groups) + ([self.pinned] if self.pinned else [])
        refit = False
        for g in groups:
            vis = tuple(b.item.shown() for b in g.buttons)
            if vis != g.shown:
                g.shown = vis
                refit = True
            for b in g.buttons:
                if b.IsShown() and b.read_state() and b.item.label() and b.mode != "icon":
                    refit = True
        if refit:
            self.fit()

    def fit(self):
        """The groups of the tab shown in the width of the window: when there is no room, the groups from the
        right get small buttons, then icons only."""
        w, h = self.body.GetClientSize()
        if w <= 0:
            return
        x = 0
        hh = self.FromDIP(BODY_H)
        if self.pinned:
            self.pinned.shown = tuple(b.item.shown() for b in self.pinned.buttons)
            self.pinned.Show(self.pinned.any_shown())  # (one view only: no view buttons)
        if self.pinned and self.pinned.IsShown():
            pw = self.pinned.width(0)
            self.pinned.SetSize(0, 0, pw, hh)
            self.pinned.place(0)
            x = pw + self.FromDIP(4)
        page = self.pages[self.current]
        page.SetSize(x, 0, w - x, hh)
        groups = [g for g in page.groups if g.any_shown()]
        for g in page.groups:
            g.shown = tuple(b.item.shown() for b in g.buttons)
            g.Show(g in groups)
        levels = [0] * len(groups)
        room = w - x - self.FromDIP(4)

        def total():
            return sum(g.width(l) for g, l in zip(groups, levels))
        for target in (1, 2):  # small buttons from the right, then icons only from the right
            for i in range(len(groups) - 1, -1, -1):
                if total() <= room:
                    break
                levels[i] = max(levels[i], target)
        gx = 0
        for g, l in zip(groups, levels):
            gw = g.width(l)
            g.SetSize(gx, 0, gw, hh)
            g.place(l)
            gx += gw
        page.Refresh()
