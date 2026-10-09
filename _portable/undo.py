"""Undo and redo in the analysis windows of MSpektra (LCMS Analysis, HRMS Analysis).

Design
* One undo stack (with redo) per window. Every open file adds its "domains": the
  Mass spectrometry view (integration, mass chromatograms, averaged and background
  ranges, spectrum labels and measurements, side panel settings; HRMS: formula check),
  the PDA view (sliders, ranges, integration, settings, MS detector delay), the
  deconvolution results of the file and their settings, the HRMS calibration; the
  window adds the Compare view (settings, files, labels, colours, shifts, blank, guide
  lines, alignment, time window, reference, graph properties, m/z tool) and the
  "link MS and PDA times" choice. A domain is a pair of functions: get() returns a
  snapshot of its state in plain Python (dicts, lists, numbers, text; heavy results
  such as a deconvolution result are kept by reference, never copied), set(state)
  restores it: the controls are set and the views recomputed as a user change would.
* commit() compares the snapshot of every domain with the one taken at the last
  commit; the domains that changed make one step (before, after). A commit runs after
  every user action: a global event filter sees the end of each mouse click, menu
  choice, button, choice, check box, spin, text entered (Enter or leaving the field)
  and key (outside text fields) in an analysis window or a dialog of it, and schedules a
  commit once the action is done (while a modal dialog is open the commit waits for it
  to close, so a whole dialog is one step). Results that arrive later (deconvolution
  worker) commit when they arrive. Because steps are differences of snapshots, an
  action reached by any path (mouse, menu, side panel, shortcut) is covered, and the
  name of the step comes from the difference (e.g. "delete peak 9.04 min").
* Quick repeats of the same kind (the same spin control, the arrow keys moving a
  slider or the picked scan) within 1.5 s are merged into one step; text fields commit
  when Enter is pressed or the field is left, so typing is one step. At most 100 steps
  per window; the steps of a file are removed when it is closed. Opening and closing
  files, zoom, pan and other view changes, exports and reports are not steps.
* Undo and Redo are refused (with a note in the status bar) while a spectrum is being
  averaged in the window, and for a step of a file whose deconvolution is still
  running; restoring never records a new step.

Hooks for the other parts (merge):
* mass shift finder (deconv_tab): a DeconvPanel (or the view) with shift_state() and
  set_shift_state(d) gets its own domain "shift" (see _shift_domain); nothing to do here.
* method presets (method_presets.apply_preset(window, preset, scope)): call
  preset_applied(window, name) right after applying (or wrap the call in
  `with undo.action(window, "apply preset NAME"):`) so that the step is named.
"""
import time

import wx

MAX_STEPS = 100
COALESCE_S = 1.5
DEBUG_TIMING = False

# icons of the top bar buttons (24 x 24, stroked like the other icons of the theme)
ICON_UNDO = '<path d="M9 14L4 9l5-5"/><path d="M4 9h10.5a5.5 5.5 0 0 1 0 11H11"/>'
ICON_REDO = '<path d="M15 14l5-5-5-5"/><path d="M20 9H9.5a5.5 5.5 0 0 0 0 11H13"/>'


def _log(*a):
    try:
        print("[undo]", *a)
    except Exception:
        pass


# ------------------------------------------------------------------ snapshots
def _copy(x):
    """Copy of the containers of a snapshot (dicts, lists, tuples, sets); other objects (numbers, text,
    arrays, results, files) are shared."""
    if isinstance(x, dict):
        return {k: _copy(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_copy(v) for v in x]
    if isinstance(x, tuple):
        return tuple(_copy(v) for v in x)
    if isinstance(x, (set, frozenset)):
        return set(x)
    return x


def _eq(a, b):
    """Equality of two snapshots (numbers by value, arrays and other objects by identity or value)."""
    if a is b:
        return True
    if isinstance(a, dict):
        if not isinstance(b, dict) or len(a) != len(b):
            return False
        for k, v in a.items():
            if k not in b or not _eq(v, b[k]):
                return False
        return True
    if isinstance(a, (list, tuple)):
        if not isinstance(b, (list, tuple)) or len(a) != len(b):
            return False
        return all(_eq(x, y) for x, y in zip(a, b))
    if isinstance(a, (set, frozenset)):
        return isinstance(b, (set, frozenset)) and set(a) == set(b)
    if isinstance(a, (bool, int, float, str)) or a is None:
        if isinstance(b, (bool, int, float, str)) or b is None:
            try:
                if isinstance(a, float) and isinstance(b, float) and a != a and b != b:
                    return True  # both nan
                return a == b and type(a is None) == type(b is None)
            except Exception:
                return False
        try:  # numpy scalars
            return bool(a == b)
        except Exception:
            return False
    try:
        import numpy as np
        if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
            return isinstance(a, np.ndarray) and isinstance(b, np.ndarray) and a.shape == b.shape and \
                bool(np.array_equal(a, b))
        if isinstance(a, np.generic):
            return bool(a == b)
    except Exception:
        pass
    return False  # other objects: the same object only


class Domain(object):
    """One part of the state of a window that undo restores: get() -> snapshot, set(snapshot)."""

    def __init__(self, key, doc, get, set_, describe=None, page=None, busy=None, order=5, passengers=(), tab=None):
        self.key, self.doc = key, doc
        self.tab = tab  # () -> the view of the domain (the step is named after the domain of the view shown)
        self.get, self.set = get, set_
        self.describe = describe  # (before, after) -> (label, coalesce key or None)
        self.page = page  # () -> callable that shows the view of the domain, or None
        self.busy = busy  # () -> text (why it cannot be restored now) or None
        self.order = order  # restored in this order (calibration before the views)
        self.passengers = tuple(passengers)  # keys of the snapshot dict restored with a step but not compared

    def compare_part(self, st):
        if not self.passengers or not isinstance(st, dict):
            return st
        return {k: v for k, v in st.items() if k not in self.passengers}


class Step(object):
    __slots__ = ("label", "changes", "t", "ckey", "t0", "key")

    def __init__(self, label, changes, ckey, key=None):
        self.label, self.changes, self.ckey = label, changes, ckey
        self.key = key  # the domain that names the step (its view is shown by Undo and Redo)
        self.t = self.t0 = time.time()


def _modal_open():
    try:
        for w in wx.GetTopLevelWindows():
            if isinstance(w, wx.Dialog) and w.IsModal():
                return True
    except Exception:
        pass
    return False


class UndoManager(object):
    def __init__(self, frame):
        self.frame = frame
        self.domains = {}
        self.base = {}
        self.undo_stack, self.redo_stack = [], []
        self.restoring = False
        self._pending = False
        self._retry = None
        self.listeners = []  # callback(manager) when the stacks changed (top bar buttons)
        self.last_commit_ms = 0.0
        self.shift_rev = _shift_rev()  # revision of the mass shift settings this window has seen
        self.link_rev = _link_rev()  # the same for "link MS and PDA times"

    # ------------------------------------------------------------- domains
    def add(self, dom):
        self.domains[dom.key] = dom
        self.base[dom.key] = self._get(dom)
        if dom.key == "shiftset":
            self.shift_rev = _shift_rev()
        if dom.key == "link":
            self.link_rev = _link_rev()

    def remove_doc(self, doc):
        """A file was closed: its domains and its part of every step go."""
        keys = [k for k, d in self.domains.items() if d.doc is doc]
        for k in keys:
            self.domains.pop(k, None)
            self.base.pop(k, None)
        for stack in (self.undo_stack, self.redo_stack):
            for st in list(stack):
                for k in keys:
                    st.changes.pop(k, None)
                if not st.changes:
                    stack.remove(st)
        self._notify()

    def has_doc(self, doc):
        return any(d.doc is doc for d in self.domains.values())

    def _get(self, dom):
        try:
            return dom.get()
        except Exception as ex:
            _log("snapshot of %s failed: %r" % (dom.key, ex))
            return None

    def rebase(self, keys=None):
        """The present state becomes the reference without a step (opening a file, a view made)."""
        for k, d in list(self.domains.items()):
            if keys is None or k in keys:
                self.base[k] = self._get(d)

    def rebase_passengers(self):
        """The reference of every domain that did not change (apart from its passengers) becomes the present state:
        a method preset is then undone to exactly the state before it (deconvolution settings saved by another
        file, calibration choices not applied yet). Busy domains keep their reference (a deconvolution running
        takes its settings with its result)."""
        for k, d in list(self.domains.items()):
            try:
                if d.busy is not None and (d.busy("commit") or d.busy("rebase")):
                    continue
            except Exception:
                pass
            old = self.base.get(k)
            cur = self._get(d)
            if old is not None and cur is not None and _eq(d.compare_part(cur), d.compare_part(old)):
                self.base[k] = cur

    def rebase_doc(self, doc):
        self.rebase([k for k, d in self.domains.items() if d.doc is doc or d.doc is None])

    # ------------------------------------------------------------- commits
    def alive(self):
        try:
            return bool(self.frame)
        except Exception:
            return False

    def schedule(self):
        """A commit once the present event is handled."""
        if self._pending or not self.alive():
            return
        self._pending = True
        wx.CallAfter(self._run)

    def _run(self):
        self._pending = False
        if not self.alive():
            return
        if _modal_open():  # a dialog of the action is still open: one step when it closes
            self._later(250)
            return
        self.commit()

    def _later(self, ms):
        if self._retry is not None:
            return

        def go():
            self._retry = None
            if self.alive():
                self.schedule()
        try:
            self._retry = wx.CallLater(ms, go)
        except Exception:
            self._retry = None

    def commit(self, label=None, ckey=None, full=False):
        """A step from the domains that changed since the last commit; True when one was recorded.
        full: passengers are compared too (a method preset changes the deconvolution settings without a run)."""
        if self.restoring or not self.alive():
            return False
        t_start = time.perf_counter()
        changes, busy = {}, False
        for k, d in list(self.domains.items()):
            try:
                if d.busy is not None and d.busy("commit"):
                    busy = True  # (averaging: committed when it is done)
                    continue
            except Exception:
                pass
            cur = self._get(d)
            old = self.base.get(k)
            if old is None or cur is None:
                self.base[k] = cur  # (a view made or a file shown: nothing to undo yet)
                continue
            if k == "shiftset" and _shift_rev() == self.shift_rev:
                self.base[k] = cur  # changed by another program window (the settings file): not a step here
                continue
            if k == "link" and _link_rev() == self.link_rev:
                self.base[k] = cur  # (the same: the choice is one setting of the program, not of this window)
                continue
            if not (_eq(cur, old) if full else _eq(d.compare_part(cur), d.compare_part(old))):
                changes[k] = (old, cur)
        if busy:
            self._later(300)
        self.last_commit_ms = (time.perf_counter() - t_start) * 1000.0
        if DEBUG_TIMING:
            _log("commit %.1f ms (%d domains)" % (self.last_commit_ms, len(self.domains)))
        if not changes:
            return False
        lab, ck, nkey = self._describe(changes)
        label = label or lab or "change"
        ckey = ckey if ckey is not None else (None if label != lab else ck)
        for k, (old, cur) in changes.items():
            self.base[k] = cur
        top = self.undo_stack[-1] if self.undo_stack else None
        now = time.time()
        if top is not None and ckey is not None and top.ckey == ckey and now - top.t <= COALESCE_S and \
                set(top.changes) == set(changes) and not self.redo_stack:
            for k, (old, cur) in changes.items():
                top.changes[k] = (top.changes[k][0], cur)
            top.t, top.label = now, label
            if all(_eq(self.domains[k].compare_part(a), self.domains[k].compare_part(b))
                   for k, (a, b) in top.changes.items() if k in self.domains):
                self.undo_stack.pop()  # back where it started
        else:
            self.undo_stack.append(Step(label, changes, ckey, nkey))
            del self.undo_stack[:-MAX_STEPS]
        self.redo_stack = []
        if "shiftset" in changes:
            _shift_settings_recorded(self)
        if "link" in changes:
            _link_recorded(self)
        self._notify()
        return True

    def _describe(self, changes):
        doms = [self.domains[k] for k in changes if k in self.domains]
        shown = None
        try:
            import shortcuts
            shown = shortcuts.view_of(self.frame)
        except Exception:
            pass

        def here(d):
            try:
                return d.tab is not None and d.tab() is shown and shown is not None
            except Exception:
                return False
        for d in sorted(doms, key=lambda d: (d.describe is None, not here(d), -_PRIORITY.get(d.key.split(":")[0], 0))):
            if d.describe is None:
                continue
            try:
                r = d.describe(*changes[d.key])
            except Exception as ex:
                _log("describe %s failed: %r" % (d.key, ex))
                r = None
            if r and r[0]:
                return r[0], r[1], d.key
        return None, None, None

    # ------------------------------------------------------------- undo
    def undo_label(self):
        return self.undo_stack[-1].label if self.undo_stack else None

    def redo_label(self):
        return self.redo_stack[-1].label if self.redo_stack else None

    def _blocked(self, step):
        for k, d in self.domains.items():
            if d.busy is None:
                continue
            try:
                why = d.busy("window") if k not in step.changes else (d.busy("window") or d.busy("step"))
            except Exception:
                why = None
            if why:
                return why
        return None

    def _status(self, text):
        try:
            self.frame.SetStatusText(text, 0)
        except Exception:
            pass

    def undo(self):
        return self._move(self.undo_stack, self.redo_stack, 0, "Undone", "Nothing to undo")

    def redo(self):
        return self._move(self.redo_stack, self.undo_stack, 1, "Redone", "Nothing to redo")

    def _move(self, src, dst, which, verb, empty):
        if self.restoring:
            return False
        if _modal_open():
            return False
        # an action not committed yet (its commit still pending) is a step first: undone first, or (redo) it
        # clears the steps that could be redone, as every new edit does (commit makes a new redo list: the
        # lists are taken again)
        self.commit()
        src, dst = (self.undo_stack, self.redo_stack) if which == 0 else (self.redo_stack, self.undo_stack)
        if not src:
            self._status(empty)
            return False
        step = src[-1]
        why = self._blocked(step)
        if why:
            self._status("%s waits: %s" % ("Undo" if which == 0 else "Redo", why))
            return False
        src.pop()
        self._apply(step, which)
        step.t = 0.0  # never merged with a later step
        dst.append(step)
        other = "Ctrl+Y redoes it" if which == 0 else "Ctrl+Z undoes it again"
        self._status("%s: %s (%s)" % (verb, step.label, other))
        self._notify()
        return True

    def _apply(self, step, which):
        t0 = time.perf_counter()
        doms = sorted((self.domains[k] for k in step.changes if k in self.domains), key=lambda d: d.order)
        self._show(doms, step.key)
        self.restoring = True
        try:
            for d in doms:
                st = _copy(step.changes[d.key][which])
                other = step.changes[d.key][1 - which]
                if d.passengers and isinstance(st, dict) and isinstance(other, dict) and \
                        _eq(d.compare_part(st), d.compare_part(other)):
                    st["_ponly"] = True  # a step of the passengers only (a method preset): they are restored
                try:
                    d.set(st)
                except Exception:
                    import traceback
                    _log("restoring %s failed:\n%s" % (d.key, traceback.format_exc()))
        finally:
            self.restoring = False
        for d in doms:
            self.base[d.key] = self._get(d)
        if any(d.key == "shiftset" for d in doms):
            _shift_settings_recorded(self)
        if any(d.key == "link" for d in doms):
            _link_recorded(self)
        _log("%s %r in %.0f ms" % ("undo" if which == 0 else "redo", step.label, (time.perf_counter() - t0) * 1000))

    def _show(self, doms, key=None):
        """The file and the view of the step are shown (the view of the domain that names it)."""
        fr = self.frame
        docs = [d.doc for d in doms if d.doc is not None]
        try:
            if docs and all(x is docs[0] for x in docs) and getattr(fr, "active", None) is not docs[0] and \
                    docs[0] in getattr(fr, "docs", []):
                fr.activate(docs[0])
            for d in sorted(doms, key=lambda d: (d.key != key, -_PRIORITY.get(d.key.split(":")[0], 0))):
                show = d.page() if d.page is not None else None
                if show is not None:
                    show()
                    break
        except Exception as ex:
            _log("show failed: %r" % ex)

    def _notify(self):
        for cb in list(self.listeners):
            try:
                cb(self)
            except RuntimeError:  # its window is gone
                try:
                    self.listeners.remove(cb)
                except ValueError:
                    pass
            except Exception as ex:
                _log("listener failed: %r" % ex)


# which domain names the step when several changed (higher first)
_PRIORITY = {"compare": 6, "dec": 5, "shift": 5, "cal": 4, "ms": 3, "pda": 2, "link": 1}


# ------------------------------------------------------------------ the window
def manager(win):
    """The undo manager of the analysis window of win (a view, a file (_Doc), a dialog or the window)."""
    w = win
    for _ in range(60):
        if w is None:
            return None
        try:
            fr = w.window() if hasattr(w, "window") and callable(getattr(w, "window")) else w
            m = fr.__dict__.get("_undo") if hasattr(fr, "__dict__") else None
        except Exception:
            m = None
        if m is not None:
            return m
        try:
            w = w.GetParent()
        except Exception:
            return None
    return None


def install(frame):
    """Undo for an analysis window (called while it is built, before its first file)."""
    m = UndoManager(frame)
    frame.__dict__["_undo"] = m
    _install_filter()
    _WINDOWS.add(frame)
    frame.Bind(wx.EVT_WINDOW_DESTROY, lambda e: (_closed(frame, e), e.Skip()))
    return m


def _closed(frame, e):
    if e.GetEventObject() is not frame:
        return
    _WINDOWS.discard(frame)
    m = frame.__dict__.get("_undo")
    if m is not None:
        m.listeners = []
    if not _WINDOWS:
        _remove_filter()


def add_buttons(frame, tb):
    """Undo and Redo in the top bar: greyed when there is nothing to undo or redo, the tooltip names the step."""
    import unidec_theme as T
    T.ICONS.setdefault("undo", ICON_UNDO)
    T.ICONS.setdefault("redo", ICON_REDO)
    FB = T._cls("FlatButton")
    bu = FB(tb, wx.ID_ANY, "Undo", icon="undo", kind="ghost", tooltip="Nothing to undo (Ctrl+Z)", height=34, padx=9)
    br = FB(tb, wx.ID_ANY, "Redo", icon="redo", kind="ghost", tooltip="Nothing to redo (Ctrl+Y)", height=34, padx=9)
    bu.Bind(wx.EVT_BUTTON, lambda e: do_undo(frame))
    br.Bind(wx.EVT_BUTTON, lambda e: do_redo(frame))
    tb.add(bu, 2)
    tb.add(br, 2)
    tb.separator()
    frame.__dict__["_undo_buttons"] = (bu, br)
    m = manager(frame)

    def update(mgr):
        try:
            u, r = mgr.undo_label(), mgr.redo_label()
            bu.Enable(u is not None)
            br.Enable(r is not None)
            bu.SetToolTip("Undo: %s (Ctrl+Z)" % u if u else "Nothing to undo (Ctrl+Z)")
            br.SetToolTip("Redo: %s (Ctrl+Y)" % r if r else "Nothing to redo (Ctrl+Y)")
        except RuntimeError:
            raise
    if m is not None:
        m.listeners.append(update)
        update(m)
    return bu, br


def do_undo(win):
    m = manager(win)
    return m.undo() if m is not None else False


def do_redo(win):
    m = manager(win)
    return m.redo() if m is not None else False


def menu_items(win):
    """Undo and Redo for the right click menus of the plots (greyed when empty), with the step named."""
    m = manager(win)
    if m is None:
        return []
    try:
        if wx.GetTopLevelParent(win) is not m.frame:
            return []  # (a plot of a dialog, e.g. the calibration window: the step is made when it closes)
    except Exception:
        pass
    u, r = m.undo_label(), m.redo_label()
    return [("Undo: %s\tCtrl+Z" % _short(u) if u else "Undo\tCtrl+Z", (lambda: m.undo()) if u else None),
            ("Redo: %s\tCtrl+Y" % _short(r) if r else "Redo\tCtrl+Y", (lambda: m.redo()) if r else None)]


def _short(text, n=60):
    text = text or ""
    return text if len(text) <= n else text[:n - 1] + "…"


def commit(win, label=None, full=False):
    m = manager(win)
    return m.commit(label, full=full) if m is not None else False


def schedule(win):
    m = manager(win)
    if m is not None:
        m.schedule()


class action(object):
    """with undo.action(window, "apply preset X"): ...  -> the changes made inside are one step with this name."""

    def __init__(self, win, label, full=False):
        self.m, self.label, self.full = manager(win), label, full

    def __enter__(self):
        if self.m is not None:
            self.m.commit()
            if self.full:
                self.m.rebase_passengers()
        return self

    def __exit__(self, *exc):
        if self.m is not None:
            self.m.commit(self.label, full=self.full)
        return False


def preset_applied(win, name, scope=None):
    """Method presets: method_presets.apply_preset wraps its work in action(window, "apply method NAME", full=True);
    this is the same for a caller that applied settings itself (the settings it changed become one step)."""
    return commit(win, "apply method %s" % name + (" (%s)" % scope if scope else ""), full=True)


# ------------------------------------------------------------------ event filter
_WINDOWS = set()
_FILTER = []


def _trigger_types():
    names = ["wxEVT_LEFT_UP", "wxEVT_TEXT_ENTER", "wxEVT_CHOICE", "wxEVT_CHECKBOX", "wxEVT_CHECKLISTBOX",
             "wxEVT_SPINCTRL", "wxEVT_SPINCTRLDOUBLE", "wxEVT_COLOURPICKER_CHANGED", "wxEVT_BUTTON", "wxEVT_MENU",
             "wxEVT_TOGGLEBUTTON", "wxEVT_RADIOBUTTON", "wxEVT_RADIOBOX", "wxEVT_SLIDER", "wxEVT_COMBOBOX"]
    out = set()
    for n in names:
        v = getattr(wx, n, None)
        if v is not None:
            out.add(v)
    return out, getattr(wx, "wxEVT_KILL_FOCUS", None), getattr(wx, "wxEVT_KEY_UP", None)


def _is_text(w):
    return isinstance(w, (wx.TextCtrl, wx.SpinCtrl, wx.ComboBox)) or \
        isinstance(getattr(w, "GetParent", lambda: None)(), (wx.SpinCtrl, wx.ComboBox))


class _Filter(wx.EventFilter):
    """Sees every event first (cheap: one set lookup); the end of a user action in an analysis window (or a dialog
    of it) schedules a commit of its undo manager."""

    def __init__(self):
        wx.EventFilter.__init__(self)
        self.types, self.kill, self.keyup = _trigger_types()

    def FilterEvent(self, event):
        try:
            t = event.GetEventType()
            if t in self.types or t == self.kill or t == self.keyup:
                obj = event.GetEventObject()
                if isinstance(obj, wx.Menu):  # a popup menu (right click): the window it was shown in
                    obj = _menu_window(obj)
                if t == self.kill and not _is_text(obj):
                    return -1
                if t == self.keyup and _is_text(obj):
                    return -1  # typing: committed on Enter or when the field is left
                m = manager(obj) if isinstance(obj, wx.Window) else None
                if m is not None and not m.restoring:
                    m.schedule()
        except Exception:
            pass
        return -1  # wx.EventFilter.Event_Skip: processed as usual


def _menu_window(menu):
    w = None
    for get in ("GetWindow", "GetInvokingWindow"):
        try:
            w = getattr(menu, get)()
        except Exception:
            w = None
        if w:
            return w
    try:
        w = wx.GetActiveWindow()
    except Exception:
        w = None
    return w or None


def _install_filter():
    if _FILTER:
        return
    try:
        f = _Filter()
        wx.EvtHandler.AddFilter(f)
        _FILTER.append(f)
    except Exception as ex:
        _log("event filter not installed: %r" % ex)


def _remove_filter():
    while _FILTER:
        f = _FILTER.pop()
        try:
            wx.EvtHandler.RemoveFilter(f)
        except Exception:
            pass


# ------------------------------------------------------------------ files
def file_opened(win, doc):
    """A file was read and shown: its domains (nothing to undo yet; the state after the CallAfters of the
    opening, such as the Compare list and an automatic calibration, is the reference)."""
    m = manager(win)
    if m is None:
        return
    if not m.has_doc(doc):
        for d in _domains_of(doc):
            m.add(d)
    m.rebase_doc(doc)
    _ensure_window_domains(m)

    def later():
        if m.alive():
            _ensure_window_domains(m)
            m.rebase_doc(doc)
    wx.CallAfter(later)


def file_closed(win, doc):
    m = manager(win)
    if m is None:
        return
    m.remove_doc(doc)

    def later():
        if m.alive():
            m.rebase([k for k, d in m.domains.items() if d.doc is None])
    wx.CallAfter(later)


def project_restored(win, docs):
    """A saved analysis was restored into the window of win (ms_project.Controller.finish_restore, which set
    restoring while it ran): the restored state is the reference of the domains of these files and of the
    window, and no step holds them (Ctrl+Z right after opening does nothing)."""
    m = manager(win)
    if m is None:
        return
    m.restoring = False
    keys = [k for k, d in m.domains.items() if d.doc is None or any(d.doc is x for x in docs)]
    m.rebase(keys)
    for stack in (m.undo_stack, m.redo_stack):
        stack[:] = [s for s in stack if not set(s.changes) & set(keys)]
    _shift_settings_recorded(m)  # (the mass shift settings restored: no step in the other windows either)
    _link_recorded(m)
    m._notify()


def _ensure_window_domains(m):
    fr = m.frame
    if "compare" not in m.domains and hasattr(fr, "compare_tab"):
        m.add(Domain("compare", None, lambda: cmp_get(fr), lambda st: cmp_set(fr, st),
                     lambda a, b: cmp_describe(a, b, fr),
                     page=lambda: _compare_page(fr), order=8, passengers=("sel",),
                     tab=lambda: fr.__dict__.get("_compare_tab")))
    if "link" not in m.domains and hasattr(fr, "link_menu") and fr.link_menu():
        m.add(Domain("link", None, _link_get, _link_set, _link_describe, order=9))
    if "shiftset" not in m.domains and any(k.startswith("shift:") for k in m.domains):
        m.add(Domain("shiftset", None, _shiftset_get, _shiftset_set, _shiftset_describe, order=3))


def _domains_of(doc):
    out = []
    ms = doc.attrs.get("ms")
    pda = doc.attrs.get("pda")
    n = id(doc)
    if ms is not None and hasattr(ms, "cols"):
        cal = getattr(ms, "cal", None)
        if cal is not None and hasattr(cal, "data_opened"):
            out.append(Domain("cal:%d" % n, doc, lambda: cal_get(cal), lambda st: cal_set(cal, st), cal_describe,
                              page=lambda: _view_page(doc, ms), busy=lambda why: _ms_busy(ms, why), order=1,
                              passengers=("note", "cfg"), tab=lambda: ms))
        out.append(Domain("ms:%d" % n, doc, lambda: ms_get(ms), lambda st: ms_set(ms, st),
                          lambda a, b: ms_describe(ms, a, b), page=lambda: _view_page(doc, ms),
                          busy=lambda why: _ms_busy(ms, why), order=2, passengers=_MS_PASSENGERS, tab=lambda: ms))
        dec = getattr(ms, "dec", None)
        if dec is not None and hasattr(dec, "results"):
            out.append(Domain("dec:%d" % n, doc, lambda: dec_get(dec), lambda st: dec_set(dec, st),
                              dec_describe, page=lambda: _view_page(doc, ms), busy=lambda why: _dec_busy(dec, why),
                              order=3, passengers=("active", "cfg", "sel"), tab=lambda: ms))
            sh = _shift_domain(doc, dec, ms)
            if sh is not None:
                out.append(sh)
            try:  # results arrive from the worker: a commit then
                dec.listeners.append(lambda pn, d=doc: schedule(d))
            except Exception:
                pass
    if pda is not None and hasattr(pda, "map_card"):
        out.append(Domain("pda:%d" % n, doc, lambda: pda_get(pda), lambda st: pda_set(pda, st),
                          lambda a, b: pda_describe(pda, a, b), page=lambda: _view_page(doc, pda),
                          busy=lambda why: None, order=4, passengers=_PDA_PASSENGERS, tab=lambda: pda))
    return out


def _shift_domain(doc, dec, tab):
    """MERGE (agent F, mass shift finder): the state of the mass shift finder of this file as its own undo
    domain when the deconvolution panel (or the view) offers shift_state() and set_shift_state(d)."""
    for obj in (dec, tab):
        get, put = getattr(obj, "shift_state", None), getattr(obj, "set_shift_state", None)
        if callable(get) and callable(put):
            # per file: on or off and the reference of each result; the settings (shared by every file and
            # window) are the window domain "shiftset"
            return Domain("shift:%d" % id(doc), doc, lambda g=get: {"results": g().get("results") or []},
                          lambda st, p=put: p({"results": st.get("results") or []}), _shift_describe,
                          page=lambda: _view_page(doc, tab), busy=lambda why: _dec_busy(dec, why), order=3,
                          tab=lambda: tab)
    return None


def _shift_rev():
    """Revision of the mass shift settings in this program (ms_shifts.SETTINGS_REV: every save adds 1)."""
    try:
        import ms_shifts
        return ms_shifts.SETTINGS_REV[0]
    except Exception:
        return 0


def _shift_settings_recorded(m):
    """A change of the mass shift settings was recorded (or undone) in the window of m: the other windows of
    this program take it as their reference (it is not a step there)."""
    m.shift_rev = _shift_rev()
    for fr in list(_WINDOWS):
        o = fr.__dict__.get("_undo") if hasattr(fr, "__dict__") else None
        if o is None or o is m or not o.alive():
            continue
        if "shiftset" in o.domains:
            o.base["shiftset"] = o._get(o.domains["shiftset"])
        o.shift_rev = m.shift_rev


def _shiftset_get():
    import deconv_shifts
    return deconv_shifts.settings()


def _shiftset_set(st):
    import deconv_shifts
    import ms_shifts
    new = ms_shifts.clean_settings({k: v for k, v in st.items() if k != "_ponly"})
    if new != deconv_shifts.settings():
        deconv_shifts.save_settings(new)
    try:
        deconv_shifts.request_redraw_all()
    except AttributeError:
        pass


def _shiftset_describe(a, b):
    try:
        if {k: v for k, v in a.items() if k != "show_all"} == {k: v for k, v in b.items() if k != "show_all"}:
            return ("every matched pair shown" if b.get("show_all") else "matched pairs: links only"), None
    except Exception:
        pass
    return "mass shift settings", None


def _shift_describe(a, b):
    """Names a step of the mass shift finder (deconv_shifts.ShiftPanel.get_state)."""
    try:
        ra = {r.get("key"): r for r in a.get("results") or []}
        for r in b.get("results") or []:
            o = ra.get(r.get("key"))
            if o is None:
                continue
            if bool(o.get("on")) != bool(r.get("on")):
                return ("mass shifts on" if r.get("on") else "mass shifts off"), None
            if o.get("ref") != r.get("ref"):
                ref = r.get("ref")
                return ("unmodified species %.1f Da" % ref if isinstance(ref, (int, float))
                        else "unmodified species: automatic"), None
        if a.get("settings") is not None and a.get("settings") != b.get("settings"):
            return "mass shift settings", None
    except Exception:
        pass
    return "mass shift finder", None


def _view_page(doc, tab):
    def show():
        try:
            if doc.book.GetPageIndex(tab) >= 0 and doc.book.GetSelection() != doc.book.GetPageIndex(tab):
                doc.show_page(tab)
        except Exception:
            pass
    return show


def _compare_page(fr):
    def show():
        try:
            doc = fr.active
            k = next((i for i, (_, pg) in enumerate(doc.pages) if getattr(pg, "is_compare", False)), -1)
            if k >= 0 and doc.book.GetSelection() != k:
                doc.show_page(doc.pages[k][1])
        except Exception:
            pass
    return show


def _ms_busy(ms, why):
    if getattr(ms, "_averaging", False):
        return "a spectrum is being averaged"
    if why == "step":
        dec = getattr(ms, "dec", None)
        if dec is not None and getattr(dec, "busy", False):
            return "the deconvolution of this file is running"
    return None


def _dec_busy(dec, why):
    if why == "rebase" and getattr(dec, "busy", False):
        return "running"
    if why == "step" and getattr(dec, "busy", False):
        return "the deconvolution of this file is running"
    return None


def _set_linking(doc, on):
    try:
        doc.attrs["_linking"] = bool(on)  # LCMSFrame.link_time: the other view is not moved while restoring
    except Exception:
        pass


# ------------------------------------------------------------------ controls
def _read(c, kind):
    if kind == "sel":
        return c.GetSelection()
    if kind == "str":
        return c.GetStringSelection()
    if kind == "bool":
        return bool(c.GetValue())
    if kind == "int":
        return int(c.GetValue())
    return c.GetValue()


def _write(c, kind, v):
    if v is None:
        return
    if kind == "sel":
        if 0 <= v < c.GetCount() and c.GetSelection() != v:
            c.SetSelection(v)
    elif kind == "str":
        if v and c.GetStringSelection() != v and c.FindString(v) >= 0:
            c.SetStringSelection(v)
    elif kind == "bool":
        if bool(c.GetValue()) != bool(v):
            c.SetValue(bool(v))
    elif kind == "int":
        if int(c.GetValue()) != int(v):
            c.SetValue(int(v))
    elif c.GetValue() != v:
        c.ChangeValue(v)  # (no text event: a list typed is not marked as typed by the user)


def _ctrl_text(c, kind, v):
    if kind == "sel":
        try:
            return c.GetString(v) if 0 <= v < c.GetCount() else str(v)
        except Exception:
            return str(v)
    if kind == "bool":
        return "on" if v else "off"
    return str(v)


def _ctrls_get(obj, spec):
    out = {}
    for name, kind, _ in spec:
        c = getattr(obj, name, None)
        if c is not None:
            try:
                out[name] = _read(c, kind)
            except (RuntimeError, ValueError, TypeError):
                pass
    return out


def _ctrls_set(obj, spec, vals):
    for name, kind, _ in spec:
        c = getattr(obj, name, None)
        if c is not None and name in vals:
            try:
                _write(c, kind, vals[name])
            except (RuntimeError, ValueError, TypeError):
                pass


def _ctrls_describe(obj, spec, a, b, unit=None):
    for name, kind, title in spec:
        if a.get(name) != b.get(name):
            c = getattr(obj, name, None)
            txt = _ctrl_text(c, kind, b.get(name)) if c is not None else str(b.get(name))
            u = (unit or {}).get(name, "")
            return ("%s: %s%s" % (title, txt if txt != "" else "empty", (" " + u) if u and txt else ""),
                    "ctrl:%s" % name)
    return None


def _t(v):
    return "%.2f min" % v


# ------------------------------------------------------------------ peaks (MS and PDA)
def _pkey(p):
    return (p.get("key") or p.get("trace"), round(float(p["t0"]), 5), round(float(p["t1"]), 5))


def _peaks_describe(a, b):
    ka, kb = set(_pkey(p) for p in a), set(_pkey(p) for p in b)
    added = [p for p in b if _pkey(p) not in ka]
    removed = [p for p in a if _pkey(p) not in kb]
    if a and not b:
        return "remove every integrated peak (%d)" % len(a)
    if len(removed) == 1 and len(added) == 2 and all(removed[0]["t0"] - 1e-9 <= q["t0"] and q["t1"] <= removed[0]["t1"]
                                                    + 1e-9 for q in added):
        return "split peak %s" % _t(removed[0]["rt"])
    if len(removed) == 1 and not added:
        return "delete peak %s" % _t(removed[0]["rt"])
    if len(added) == 1 and not removed:
        return "add peak %s" % _t(added[0]["rt"])
    if added or removed:
        tr = (added or removed)[0].get("trace") or ""
        tr = " (%s)" % tr if tr else ""
        if not added:
            return "integrate: no peak found, %d removed%s" % (len(removed), tr)
        return "integrate %d peak%s%s" % (len(added), "" if len(added) == 1 else "s", tr)
    return None  # the same peaks, measured again (the trace changed): named by the change of the trace


# ------------------------------------------------------------------ Mass spectrometry view
MS_CTRLS = [("kind", "sel", "Trace"), ("xic_mode", "sel", "Mass chromatograms shown"), ("smooth", "int", "Smoothing"),
            ("norm", "bool", "Scale each trace to 100 %"), ("link", "bool", "Link the time axes"),
            ("tol", "text", "Default window"), ("tol_unit", "sel", "Default window unit"),
            ("use_bg", "bool", "Subtract background"), ("binw", "text", "Bin width"),
            ("spec_style", "sel", "Spectrum display"), ("use_ev", "sel", "Spectrum to use"),
            ("formula", "text", "Formula"), ("adduct", "sel", "Ion"), ("ppm", "text", "Tolerance"),
            ("int_trace", "str", "Integration trace"), ("int_thr", "text", "Minimum height"),
            ("int_w", "text", "Minimum width")]
MS_UNITS = {"smooth": "points", "binw": "m/z", "ppm": "ppm", "int_thr": "% of top", "int_w": "s"}
_CHROM_CTRLS = ("kind", "xic_mode", "smooth")
_SPEC_CTRLS = ("use_bg", "binw")
_MS_PASSENGERS = ("sel_peak", "fields", "xic_text", "note", "finder", "method")


def ms_get(tab):
    if getattr(tab, "data", None) is None:
        return None
    cols = []
    for c in tab.cols:
        cols.append({"pins": list(c.get("pins") or []), "measure": c.get("measure"), "m1": c.get("m1"),
                     "overlay": c.get("overlay")})
    st = {"ctrl": _ctrls_get(tab, MS_CTRLS),
          "peaks": [dict(p) for p in tab.peaks], "sel_peak": tab.sel_peak,
          "avg": tab.avg, "bg": tab.bg, "pick_t": tab.pick_t,
          "xics": {e: [dict(x) for x in lst] for e, lst in tab.xics.items()},
          "cols": cols,
          "fields": tuple(c.GetValue() for c in (tab.t0, tab.t1, tab.b0, tab.b1)),
          "xic_text": tuple(c.GetValue() for c in tab.xic_fields)}
    note = getattr(tab, "exact_note", None)
    if note is not None:
        st["note"] = note.GetLabel()
    try:  # the name on the Method chip of the file (method presets)
        st["method"] = object.__getattribute__(tab.frame, "attrs").get("method_name")
    except Exception:
        pass
    if callable(getattr(tab, "finder_get", None)):  # HRMS: the formula finder (element limits, ranking)
        try:
            st["finder"] = tab.finder_get()
        except Exception:
            pass
    return st


def ms_set(tab, st):
    if st is None or getattr(tab, "data", None) is None:
        return
    cur = ms_get(tab)
    doc = tab.frame
    _set_linking(doc, True)
    try:
        c0, c1 = cur["ctrl"], st["ctrl"]
        chrom = any(c0.get(k) != c1.get(k) for k in _CHROM_CTRLS) or not _eq(cur["xics"], st["xics"])
        spec = any(c0.get(k) != c1.get(k) for k in _SPEC_CTRLS) or \
            not _eq((cur["avg"], cur["bg"], cur["pick_t"]), (st["avg"], st["bg"], st["pick_t"]))
        peaks = not _eq(cur["peaks"], st["peaks"]) or cur.get("sel_peak") != st.get("sel_peak")
        chroms = chrom or spec or peaks or c0.get("norm") != c1.get("norm")
        style = c0.get("spec_style") != c1.get("spec_style")
        marks = [not _eq(a, b) for a, b in zip(cur["cols"], st["cols"])]
        _ctrls_set(tab, MS_CTRLS, c1)
        tab.xics = {int(e): [dict(x) for x in lst] for e, lst in st["xics"].items()}
        for e in range(len(tab.xic_fields)):
            tab._set_xic_field(e)
        try:
            tab._xic_dirty.clear()
        except AttributeError:
            pass
        if chrom:
            tab.on_update()
            _ctrls_set(tab, MS_CTRLS, c1)  # (the trace of the integration, chosen by its name)
        tab.peaks = [dict(p) for p in st["peaks"]]
        sp = st.get("sel_peak")
        tab.sel_peak = sp if sp is not None and 0 <= sp < len(tab.peaks) else None
        if peaks or chrom:
            tab.table.set_peaks(tab.peaks, tab.units)
        tab.avg, tab.bg, tab.pick_t = st["avg"], st["bg"], st["pick_t"]
        if spec:
            if tab.avg and tab.avg[1] != tab.avg[0]:
                tab.compute_average()
            elif tab.pick_t is not None:
                tab.chrom_pick(tab.pick_t, None, None)
            elif tab.avg:
                tab.chrom_pick(tab.avg[0], None, None)
            else:
                for col in tab.cols:
                    col.update(spec=None, range=None, desc="")
                tab.spectra_changed()
            if callable(getattr(tab, "_remember_spectrum_selection", None)):
                tab._remember_spectrum_selection()  # (a failed averaging later goes back to this, not before it)
        for col, cs in zip(tab.cols, st["cols"]):
            col["pins"] = list(cs.get("pins") or [])
            col["measure"], col["m1"] = cs.get("measure"), cs.get("m1")
            if cs.get("overlay") is not None:
                col["overlay"] = cs["overlay"]
            else:
                col.pop("overlay", None)
        for c, v in zip((tab.t0, tab.t1, tab.b0, tab.b1), st.get("fields") or ()):
            if c.GetValue() != v:
                c.ChangeValue(v)
        for c, v in zip(tab.xic_fields, st.get("xic_text") or ()):
            if c.GetValue() != v:
                c.ChangeValue(v)
        if "method" in st:
            _method_chip(tab, st.get("method"))
        if st.get("_ponly") and st.get("finder") is not None and callable(getattr(tab, "finder_set", None)):
            try:  # a method preset: the limits of the formula finder go back too
                tab.finder_set(st["finder"])
            except Exception:
                pass
        if "note" in st and getattr(tab, "exact_note", None) is not None:
            tab.exact_note.SetLabel(st["note"])
            try:
                tab._rewrap()
            except Exception:
                pass
        # drawn again only where something changed (an undo stays as quick as the edit it undoes)
        if chroms:
            tab.plot_chroms(keep_view=True)
        for k, (col, cs) in enumerate(zip(tab.cols, st["cols"])):
            if style or (k < len(marks) and marks[k]) or (spec and (cs.get("pins") or cs.get("measure") or
                                                                   cs.get("overlay") is not None)):
                tab.plot_spec(col, keep_view=True)
    finally:
        _set_linking(doc, False)


def _method_chip(tab, name):
    """The Method chip of the file back to the name of a step (a method renamed or deleted since: left)."""
    try:
        import method_presets as MP
        doc = tab.frame
        if (object.__getattribute__(doc, "attrs").get("method_name") or None) == (name or None):
            return
        if name and name not in MP.names(MP.kind_of(doc)):
            return
        MP._set_chip(doc, name)
    except Exception as ex:
        _log("method chip: %r" % ex)


def _xic_name(tab, x, e):
    try:
        return ("m/z " + tab.PIN_FMT + " %s") % (x["mz"], tab.ev_short(e))
    except Exception:
        return "m/z %g" % x.get("mz", 0)


_MS_TRACE_CTRLS = [c for c in MS_CTRLS if c[0] in ("kind", "xic_mode", "smooth")]


def ms_describe(tab, a, b):
    # the mass chromatograms and the settings that change the traces first: their peaks are measured again or
    # dropped with them (the step is the change of the trace, not of the peaks)
    if _eq(a["xics"], b["xics"]):
        r = _ctrls_describe(tab, _MS_TRACE_CTRLS, a["ctrl"], b["ctrl"], MS_UNITS)
        if r:
            return r
    r = _xics_describe(tab, a, b)
    if r:
        return r
    if not _eq(a["peaks"], b["peaks"]):
        t = _peaks_describe(a["peaks"], b["peaks"])
        if t:
            return t, None
    return _ms_describe_rest(tab, a, b)


def _xics_describe(tab, a, b):
    if not _eq(a["xics"], b["xics"]):
        for e in sorted(set(a["xics"]) | set(b["xics"])):
            la, lb = a["xics"].get(e, []), b["xics"].get(e, [])
            ka = [(x["mz"], x["w"], x["ppm"]) for x in la]
            kb = [(x["mz"], x["w"], x["ppm"]) for x in lb]
            if ka == kb:
                continue
            add = [x for x in lb if (x["mz"], x["w"], x["ppm"]) not in ka]
            rem = [x for x in la if (x["mz"], x["w"], x["ppm"]) not in kb]
            if add and rem:
                return "edit mass chromatogram %s" % _xic_name(tab, add[0], e), None
            if add:
                more = " and %d more" % (len(add) - 1) if len(add) > 1 else ""
                return "add mass chromatogram %s%s" % (_xic_name(tab, add[0], e), more), None
            if rem:
                return ("remove mass chromatogram %s" % _xic_name(tab, rem[0], e) if len(rem) == 1 else
                        "remove %d mass chromatograms" % len(rem)), None
            return "order of the mass chromatograms", None
    return None


def _ms_describe_rest(tab, a, b):
    if (a["avg"], a["pick_t"]) != (b["avg"], b["pick_t"]):
        if b["avg"] and b["avg"][1] != b["avg"][0]:
            return "average %.2f to %.2f min" % tuple(b["avg"]), None
        if b["pick_t"] is not None:
            return "spectra at %.3f min" % b["pick_t"], "ms-pick"
        return "clear the spectra", None
    if a["bg"] != b["bg"]:
        return ("background range %.2f to %.2f min" % tuple(b["bg"]) if b["bg"] else "clear the background range"), None

    def ev_of(k):
        try:
            t = tab.ev_short(tab.cols[k]["e"])
            return " " + t if t.startswith("(") else " (%s)" % t
        except Exception:
            return ""
    for k, (ca, cb) in enumerate(zip(a["cols"], b["cols"])):
        if ca.get("overlay") is not cb.get("overlay"):  # (before the fields of the formula: the check names it)
            ov = cb.get("overlay")
            if ov:
                name = (ov.get("info") or {}).get("name") or "formula"
                return "check formula %s%s" % (name, ev_of(k)), None
            return "remove the isotope pattern%s" % ev_of(k), None
    r = _ctrls_describe(tab, MS_CTRLS, a["ctrl"], b["ctrl"], MS_UNITS)
    if r:  # (a new spectrum also removes a measurement: named by the setting)
        return r
    for k, (ca, cb) in enumerate(zip(a["cols"], b["cols"])):
        ev = ev_of(k)
        if not _eq(ca["pins"], cb["pins"]):
            add = [m for m in cb["pins"] if m not in ca["pins"]]
            rem = [m for m in ca["pins"] if m not in cb["pins"]]
            if add:
                return ("label m/z " + tab.PIN_FMT + "%s") % (add[0], ev), None
            if rem and not cb["pins"] and len(rem) > 1:
                return "remove every label%s" % ev, None
            if rem:
                return ("remove the label of m/z " + tab.PIN_FMT + "%s") % (rem[0], ev), None
        if not _eq((ca["measure"], ca["m1"]), (cb["measure"], cb["m1"])):
            m = cb["measure"]
            if m and m[1] is not None:
                return ("measure m/z " + tab.PIN_FMT + " to " + tab.PIN_FMT + "%s") % (m[0], m[1], ev), None
            if m:
                return ("measure from m/z " + tab.PIN_FMT + "%s") % (m[0], ev), None
            return "remove the measurement%s" % ev, None
    return ("peak areas measured again", None) if not _eq(a["peaks"], b["peaks"]) else None


# ------------------------------------------------------------------ PDA view
PDA_CTRLS = [("bw", "text", "Bandwidth"), ("cscale", "sel", "Colour scale"), ("wls", "text", "More wavelengths"),
             ("maxplot", "bool", "Max plot"), ("smooth", "int", "Smoothing"),
             ("norm", "bool", "Scale each trace to 100 %"), ("overlay", "bool", "Overlay the MS trace"),
             ("ms_event", "sel", "MS trace"), ("delay", "text", "MS detector delay"),
             ("use_bg", "bool", "Subtract background"), ("int_trace", "str", "Integration trace"),
             ("int_thr", "text", "Minimum height"), ("int_w", "text", "Minimum width")]
PDA_UNITS = {"bw": "nm", "wls": "nm", "smooth": "points", "delay": "min", "int_thr": "% of top", "int_w": "s"}
_PDA_PASSENGERS = ("sel_peak", "fields")
_PDA_TRACE_CTRLS = [c for c in PDA_CTRLS if c[0] in ("bw", "wls", "maxplot", "smooth")]


def pda_get(tab):
    if getattr(tab, "pda", None) is None:
        return None
    return {"ctrl": _ctrls_get(tab, PDA_CTRLS), "t": tab.t_sel, "wl": tab.wl_sel, "avg": tab.avg, "bg": tab.bg,
            "peaks": [dict(p) for p in tab.peaks], "sel_peak": tab.sel_peak, "pins": list(tab.uv_pins),
            "fields": tuple(c.GetValue() for c in (tab.t0, tab.t1, tab.b0, tab.b1))}


def pda_set(tab, st):
    if st is None or getattr(tab, "pda", None) is None:
        return
    cur = pda_get(tab)
    doc = tab.frame
    _set_linking(doc, True)
    c0, c1 = cur["ctrl"], st["ctrl"]
    traces = any(c0.get(k) != c1.get(k) for k in ("bw", "wls", "maxplot", "smooth", "overlay", "ms_event", "delay")) \
        or cur["wl"] != st["wl"]
    uv = cur["t"] != st["t"] or cur["wl"] != st["wl"] or not _eq((cur["avg"], cur["bg"]), (st["avg"], st["bg"])) or \
        cur["pins"] != st["pins"] or c0.get("use_bg") != c1.get("use_bg")
    peaks = not _eq(cur["peaks"], st["peaks"]) or cur.get("sel_peak") != st.get("sel_peak")
    chroms = traces or uv or peaks or c0.get("norm") != c1.get("norm")
    try:
        _ctrls_set(tab, PDA_CTRLS, st["ctrl"])
        tab.t_sel, tab.wl_sel = st["t"], st["wl"]
        tab.avg, tab.bg = st["avg"], st["bg"]
        tab.uv_pins = list(st["pins"])
        tab._last = (tab.t_sel, tab.wl_sel)
        tab._sync_fields()
        if c0.get("cscale") != c1.get("cscale"):
            tab.plot_map()
        elif (cur["t"], cur["wl"]) != (st["t"], st["wl"]):
            tab.map_card.set_cross(tab.t_sel, tab.wl_sel)
        if traces:
            tab.on_update(keep=True)
            _ctrls_set(tab, PDA_CTRLS, st["ctrl"])  # (the trace of the integration, chosen by its name)
        tab.peaks = [dict(p) for p in st["peaks"]]
        sp = st.get("sel_peak")
        tab.sel_peak = sp if sp is not None and 0 <= sp < len(tab.peaks) else None
        if peaks or traces:
            tab.table.set_peaks(tab.peaks, tab.units)
        for c, v in zip((tab.t0, tab.t1, tab.b0, tab.b1), st.get("fields") or ()):
            if c.GetValue() != v:
                c.ChangeValue(v)
        if chroms:
            tab.plot_chroms(keep_view=True)
        if uv:
            tab.compute_uv()
    finally:
        _set_linking(doc, False)


def pda_describe(tab, a, b):
    r = _ctrls_describe(tab, _PDA_TRACE_CTRLS, a["ctrl"], b["ctrl"], PDA_UNITS)
    if r:  # (before the peaks they measure again or drop)
        return r
    if not _eq(a["peaks"], b["peaks"]):
        t = _peaks_describe(a["peaks"], b["peaks"])
        if t:
            return t, None
    if (a["t"], a["wl"]) != (b["t"], b["wl"]):
        if a["wl"] == b["wl"]:
            return "PDA time %.3f min" % b["t"], "pda-cross"
        if a["t"] == b["t"]:
            return "PDA wavelength %.0f nm" % b["wl"], "pda-cross"
        return "PDA %.3f min, %.0f nm" % (b["t"], b["wl"]), "pda-cross"
    if a["avg"] != b["avg"]:
        if b["avg"] and b["avg"][1] != b["avg"][0]:
            return "UV average %.2f to %.2f min" % tuple(b["avg"]), None
        return "UV spectrum at one time", None
    if a["bg"] != b["bg"]:
        return ("UV background %.2f to %.2f min" % tuple(b["bg"]) if b["bg"] else "clear the UV background"), None
    if a["pins"] != b["pins"]:
        add = [w for w in b["pins"] if w not in a["pins"]]
        rem = [w for w in a["pins"] if w not in b["pins"]]
        return ("label %.0f nm" % add[0] if add else "remove the label of %.0f nm" % rem[0] if rem else "labels"), None
    return _ctrls_describe(tab, PDA_CTRLS, a["ctrl"], b["ctrl"], PDA_UNITS) or \
        (("peak areas measured again", None) if not _eq(a["peaks"], b["peaks"]) else None)


# ------------------------------------------------------------------ deconvolution results
_DCFG = {}


def _dec_cfg(pn):
    """The deconvolution settings in effect for pn: the saved ones (pn._saved_cfg(), what the settings window opens
    with; written by any file or window of this kind and by method presets), else pn.cfg. Read again only when the
    settings file changed."""
    try:
        import unilcms as U
        st = U.settings()
    except Exception:
        st = None
    c = _DCFG.get(pn.key)
    if st is None or c is None or c[0] is not st:
        try:
            cfg = pn._saved_cfg()
        except Exception:
            cfg = None
        c = (st, cfg)
        _DCFG[pn.key] = c
    return dict(c[1] or pn.cfg or {})


def dec_get(pn):
    res = []
    for ent in pn.results:
        r = ent.res or {}
        res.append((ent, bool(ent.visible), bool(r.get("label_all")), [dict(q) for q in (r.get("peaks") or [])]))
    return {"results": res, "active": pn.active, "cfg": _dec_cfg(pn),
            "sel": tuple(e.sel for e in pn.results)}


def dec_set(pn, st, save=True):
    """save=False (a saved analysis being opened): the settings go to the panel only, the deconvolution settings
    saved for every file and window are kept."""
    import unidec_theme as T
    cur = list(pn.results)
    # the deconvolution settings go back with a step that adds or removes a result (the run they started) or
    # with a step of the settings alone (a method preset); not with other edits of the results (they are
    # shared by every file and window: a run elsewhere since must not be undone here)
    restore_cfg = bool(st.get("_ponly")) or [r[0] for r in st["results"]] != cur
    ents = [r[0] for r in st["results"]]
    changed_peaks = []
    for ent, vis, lab, peaks in st["results"]:
        ent.visible = bool(vis)
        if lab:
            ent.res["label_all"] = True
        else:
            ent.res.pop("label_all", None)
        if not _eq(ent.res.get("peaks"), peaks):
            changed_peaks.append(ent)
        ent.res["peaks"] = [dict(q) for q in peaks]
        if ent.sel is not None and ent.sel >= len(ent.res["peaks"]):
            ent.sel = None
        if ent.row is None:
            pn._make_row(ent)
    for ent, s in zip(ents, st.get("sel") or ()):
        if s is None or s < len(ent.res.get("peaks") or []):
            ent.sel = s
    gone = [e for e in cur if e not in ents]
    pn.results = list(ents)
    act = st.get("active")
    if act not in ents or not act.visible:
        act = pn._newest_visible()
    pn.active = act
    pn._layout()
    for e in gone:
        pn._drop_row(e)
    for ent in pn.results:
        if ent.visible and ent.row is not None:
            try:
                pn.plot_results(ent)
            except Exception as ex:
                _log("deconvolution result not drawn: %r" % ex)
    if act is not None:
        pn._ensure_table()
        pn._fill_table(act)
    pn._mark_rows()
    events = {e.res["e"] for e in list(cur) + ents if e.res is not None and "e" in e.res}
    pn._replot_spec(events)
    for ent in changed_peaks:  # the mass table written next to the data, as after adding a mass by hand
        if ent.res.get("saved"):
            try:
                pn._rewrite_peaks(ent.res)  # only while the files there are still this result's
            except Exception:
                pass
    cfg = st.get("cfg")
    if restore_cfg and cfg and (not _eq(cfg, pn.cfg) or not _eq(cfg, _dec_cfg(pn))):
        pn.cfg = dict(cfg)
        if save:
            try:
                T._save({pn.key: dict(cfg)})  # the settings window opens with them
            except Exception:
                pass
    pn._notify()


def _res_name(ent):
    r = ent.res or {}
    meth = (r.get("method") or r.get("tag") or "deconvolution").split(" (")[0]
    return "%s, %s" % (meth, r.get("input", "")) if r.get("input") else meth


def dec_describe(a, b):
    ea = [r[0] for r in a["results"]]
    eb = [r[0] for r in b["results"]]
    add = [e for e in eb if e not in ea]
    rem = [e for e in ea if e not in eb]
    if add:
        return "deconvolution (%s)" % _res_name(add[0]), None
    if rem:
        if not eb and len(rem) > 1:
            return "close every deconvolution result", None
        return "close the result (%s)" % _res_name(rem[0]), None
    da = {id(r[0]): r for r in a["results"]}
    for ent, vis, lab, peaks in b["results"]:
        o = da.get(id(ent))
        if o is None:
            continue
        if o[1] != vis:
            return ("show the result (%s)" if vis else "hide the result (%s)") % _res_name(ent), None
        if o[2] != lab:
            return ("label all peaks in view" if lab else "automatic peak labels"), None
        if not _eq(o[3], peaks):
            ma = [(q.get("mass"), q.get("added")) for q in o[3]]
            mb = [(q.get("mass"), q.get("added")) for q in peaks]
            if len(mb) > len(ma):
                q = next((q for q in peaks if (q.get("mass"), q.get("added")) not in ma), peaks[-1])
                return "add the mass %.6g Da" % q["mass"], None
            if len(mb) < len(ma):
                q = next((q for q in o[3] if (q.get("mass"), q.get("added")) not in mb), o[3][-1])
                return "remove the mass %.6g Da" % q["mass"], None
            for qa, qb in zip(o[3], peaks):
                if bool(qa.get("pinned")) != bool(qb.get("pinned")):
                    return ("label %.6g Da" if qb.get("pinned") else "remove the label of %.6g Da") % qb["mass"], None
            return "masses of the result", None
    if ea != eb:
        return "order of the results", None
    return None


# ------------------------------------------------------------------ HRMS calibration
def cal_get(cp):
    d = cp.data
    if d is None:
        return None
    return {"calibs": [tuple(c) for c in (getattr(d, "calibs", None) or [])], "calib": d.calib,
            "info": getattr(d, "calib_info", ""), "events": getattr(d, "calib_events", None),
            "panel": (cp.cal, cp.range, cp.event, cp.spec), "rows": [dict(r) for r in (cp.rows or [])],
            "auto": bool(cp.auto.GetValue()), "note": cp.note.GetLabel(), "cfg": dict(cp.cfg or {})}


def cal_set(cp, st, save=True):
    """save=False (a saved analysis being opened): the calibration choices go to the panel only, the ones saved
    for every file and window are kept."""
    d = cp.data
    if st is None or d is None:
        return

    def put_cfg(**kw):
        if save:
            cp.save_cfg(**kw)
        else:
            cp.cfg = dict(cp.cfg or {}, **kw)
    cur = cal_get(cp)
    if not _eq(cur["calibs"], st["calibs"]) or cur["calib"] is not st["calib"] or \
            not _eq(cur["panel"][:3], st["panel"][:3]):
        d.calibs = [tuple(c) for c in st["calibs"]]
        d.calib, d.calib_info, d.calib_events = st["calib"], st["info"], st["events"]
        cp.cal, cp.range, cp.event, cp.spec = st["panel"]
        cp.rows = [dict(r) for r in st["rows"]]
        cp._set_note(st.get("note", ""))
        cfg = st.get("cfg")
        if cfg and not _eq(cfg, cp.cfg):  # the choices of the calibration window (calibrant, model, references)
            try:
                put_cfg(**{k: v for k, v in cfg.items() if k != "auto"})
            except Exception:
                pass
        cp.frame.calibration_changed()
    else:
        cp._set_note(st.get("note", ""))
        cfg = st.get("cfg")
        if st.get("_ponly") and cfg and not _eq(cfg, cp.cfg):  # a method preset: the choices of the window
            try:
                put_cfg(**{k: v for k, v in cfg.items() if k != "auto"})
            except Exception:
                pass
    if bool(cp.auto.GetValue()) != bool(st["auto"]):
        cp.auto.SetValue(bool(st["auto"]))
        try:
            put_cfg(auto=bool(st["auto"]))
        except Exception:
            pass


def cal_describe(a, b):
    if a["calib"] is not b["calib"] or not _eq(a["calibs"], b["calibs"]) or not _eq(a["panel"][:3], b["panel"][:3]):
        c = b["calib"]
        if c is None and not b["calibs"]:
            return "remove the calibration", None
        try:
            return "calibration (%s)" % c.describe(), None
        except Exception:
            return "calibration", None
    if a["auto"] != b["auto"]:
        return "Calibrate automatically on opening: %s" % ("on" if b["auto"] else "off"), None
    return None


# ------------------------------------------------------------------ link MS and PDA
def _link_rev():
    """Revision of the "link MS and PDA times" choice in this program (unidec_theme._save counts the saves of each
    setting): unchanged while only another program changed it."""
    try:
        import unidec_theme as T
        return T._ST.get("revs", {}).get("link_ms_pda", 0)
    except Exception:
        return 0


def _link_recorded(m):
    """A change of the link choice was recorded (or undone) in the window of m: the other windows of this program
    take it as their reference (as _shift_settings_recorded)."""
    m.link_rev = _link_rev()
    for fr in list(_WINDOWS):
        o = fr.__dict__.get("_undo") if hasattr(fr, "__dict__") else None
        if o is None or o is m or not o.alive():
            continue
        if "link" in o.domains:
            o.base["link"] = o._get(o.domains["link"])
        o.link_rev = m.link_rev


def _link_get():
    import unidec_theme as T
    return {"link": bool(T._load().get("link_ms_pda", True))}


def _link_set(st):
    import unidec_theme as T
    T._save({"link_ms_pda": bool(st["link"])})


def _link_describe(a, b):
    return ("Link MS and PDA times: %s" % ("on" if b["link"] else "off")), None


# ------------------------------------------------------------------ Compare view
_CMP_NAMES = {"signal": "Trace", "wl": "Wavelength", "bw": "Bandwidth", "own_wl": "Each file at its own wavelength",
              "polarity": "MS scans", "mz": "m/z", "mz_win": "m/z window", "mz_ppm": "m/z window in ppm",
              "smooth": "Smoothing", "baseline": "Baseline", "rolling_min": "Rolling window",
              "align_win": "Alignment window", "layout": "Layout", "scale": "Scale", "ref_t": "Reference time",
              "spacing": "Spacing", "skew": "Skew", "off_spacing": "Spacing", "off_skew": "Skew",
              "reverse": "First file at the bottom", "colours": "Colours", "lw": "Line width", "fill": "Fill",
              "labels": "Labels", "label_text": "Label text", "rt_labels": "Retention times", "rt_min": "Threshold",
              "yaxis": "Y axis", "spec_avg": "Spectra averaging", "spec_bg": "Subtract the spectrum before the peak",
              "spec_labels": "Ions labelled"}
_CMP_UNITS = {"wl": "nm", "bw": "nm", "smooth": "points", "rolling_min": "min", "align_win": "min",
              "ref_t": "min", "spacing": "%", "skew": "%", "off_spacing": "%", "off_skew": "%", "lw": "pt",
              "rt_min": "%"}


def _cmp_choices():
    try:
        import lcms_compare_core as K
        import lcms_compare as LC
        return {"signal": K.SIGNALS, "baseline": K.BASELINES, "layout": K.LAYOUTS, "scale": K.SCALES,
                "colours": K.COLOURS, "labels": K.LABELS, "label_text": K.LABEL_TEXT, "rt_labels": K.RT_LABELS,
                "yaxis": K.YAXES, "spec_avg": LC.SPEC_AVG}
    except Exception:
        return {}


_PICK_KEYS = ("kind", "pda", "delay", "ms_apex", "ms_t0", "ms_t1", "trace_apex", "signal", "file_delay", "label")


def cmp_get(fr):
    tab = fr.__dict__.get("_compare_tab")
    if tab is None:
        return None
    pk = tab.mzpick
    pick = None
    if pk:
        pick = {k: pk.get(k) for k in _PICK_KEYS}
        pick["doc"] = pk.get("doc")
    cols = {}
    for pol, c in tab.mzcols.items():
        cols[pol] = {"pins": list(c.get("pins") or []), "measure": c.get("measure"), "m1": c.get("m1"),
                     "pin_doc": c.get("pin_doc")}
    # The tab keys traces by id(doc), which means nothing in another process: the snapshot holds the
    # files themselves (a saved project stores them as references to its raw files).
    files = {id(en["doc"]): en["doc"] for en in tab.entries}
    # the region as applied (not the text typed in its fields until it is used), the reference chosen (also
    # while its trace is left out), the X values and X axis label of the areas shown
    area = {"range": tab.area_range, "reference": files.get(getattr(tab, "_area_ref_key", None)),
            "x": getattr(tab, "_area_x_used", tab.area_x.GetValue()),
            "xlabel": getattr(tab, "_area_xlabel_used", tab.area_xlabel.GetValue()), "y": tab.area_y.GetSelection()}
    props = _copy(tab.props)
    props["traces"] = {files[k]: v for k, v in props.get("traces", {}).items() if k in files}
    return {"s": _copy(tab.s), "area": area,
            "entries": [(en["doc"], bool(en["include"]), en.get("label"), en.get("colour"),
                         float(en.get("shift", 0.0) or 0.0), tuple(en["mzbg"]) if en.get("mzbg") else None)
                        for en in tab.entries],
            "blank": tab.blank_doc, "props": props, "mz_on": bool(tab.mz_on),
            "spec_style": tab._spec_style, "pick": pick, "cols": cols,
            "sel": tab.entries[tab.sel]["doc"] if tab.sel is not None and tab.sel < len(tab.entries) else None}


def _g10(v):
    return "" if v is None else "%.10g" % v


def cmp_set(fr, st):
    tab = fr.__dict__.get("_compare_tab")
    if tab is None or st is None:
        return
    keep_view = (tab.s.get("t0"), tab.s.get("t1")) == (st["s"].get("t0"), st["s"].get("t1")) and \
        (tab.props.get("x0"), tab.props.get("x1")) == (st["props"].get("x0"), st["props"].get("x1"))
    if bool(tab.mz_on) != bool(st["mz_on"]):
        tab.set_mz(st["mz_on"])  # (off: the background ranges and the pick go; restored below)
    tab.s = _copy(st["s"])
    tab._fill_controls()
    tab._busy = True
    try:
        s = tab.s
        for c, v in ((tab.t0, _g10(s.get("t0"))), (tab.t1, _g10(s.get("t1"))), (tab.align, _g10(s.get("align"))),
                     (tab.ref_t, _g10(s.get("ref_t"))),
                     (tab.guides, ", ".join("%.10g" % g for g in (s.get("guides") or [])))):
            if c.GetValue() != v:
                c.ChangeValue(v)
    finally:
        tab._busy = False
    docs = tab._docs()
    ents = []
    for doc, inc, lab, col, sh, mzbg in st["entries"]:
        if doc in docs and not any(e["doc"] is doc for e in ents):
            en = {"doc": doc, "include": inc, "label": lab, "colour": col, "shift": sh}
            if mzbg:
                en["mzbg"] = tuple(mzbg)
            ents.append(en)
    for d in docs:  # files opened after the step: as files_changed adds them
        if not any(e["doc"] is d for e in ents):
            ents.append({"doc": d, "include": True, "label": None, "colour": None, "shift": 0.0})
    tab.entries = ents
    tab.blank_doc = st["blank"] if st["blank"] in docs else None
    tab.props = _copy(st["props"])
    # back to the tab's id(doc) keys; projects saved by 3.4 hold ids of another process: dropped
    tab.props["traces"] = {id(d): _copy(v) for d, v in st["props"].get("traces", {}).items() if d in docs}
    area = st.get("area")
    if area is not None:
        rng = area["range"]
        tab.area_range = tuple(rng) if rng is not None else None
        ref = area.get("reference")
        tab._area_ref_key = id(ref) if ref is not None and ref in docs else None
        for ctrl, value in zip((tab.area_lo, tab.area_hi), ("%.4f" % rng[0], "%.4f" % rng[1]) if rng else ("", "")):
            ctrl.ChangeValue(value)  # (as set_area_range shows the region)
        tab.area_x.ChangeValue(area["x"])
        tab.area_xlabel.ChangeValue(area["xlabel"])
        tab.area_y.SetSelection(area["y"])
    tab._spec_style = st["spec_style"]
    sel = st.get("sel")
    tab.sel = next((i for i, e in enumerate(ents) if e["doc"] is sel), None) if sel is not None else tab.sel
    if tab.sel is not None and tab.sel >= len(ents):
        tab.sel = None
    tab._save()
    tab._enable()
    tab._fill_list()
    pk = st.get("pick")
    if pk and pk.get("doc") in docs and tab.mz_on:
        en = next((e for e in ents if e["doc"] is pk["doc"]), None)
        if en is not None:
            new = dict(pk)
            new.update(entry=en, spectra={}, masses=[], seconds=0.0)
            tab.mzpick = new
            tab._mz_compute()
    else:
        tab.mzpick = None
    for pol, cs in (st.get("cols") or {}).items():
        c = tab.mzcols.get(pol)
        if c is not None:
            c.update(pins=list(cs.get("pins") or []), measure=cs.get("measure"), m1=cs.get("m1"),
                     pin_doc=cs.get("pin_doc"))
    tab.card.ylock = None
    for p in ("+", "-"):
        tab.plot_mz(p, keep_view=True)
    tab.replot(keep_view=keep_view)


def _cmp_label_of(fr, doc):
    tab = fr.__dict__.get("_compare_tab")
    try:
        en = next(e for e in tab.entries if e["doc"] is doc)
        return tab._base_label(en)
    except Exception:
        import os
        try:
            return os.path.splitext(os.path.basename((doc.attrs.get("path") or "").rstrip("\\/")))[0] or "a file"
        except Exception:
            return "a file"


def cmp_describe(a, b, fr=None):
    sa, sb = a["s"], b["s"]
    if a["mz_on"] != b["mz_on"]:
        return "m/z tool %s" % ("on" if b["mz_on"] else "off"), None
    pa, pb = a.get("pick"), b.get("pick")
    if not _eq(pa, pb) and pb:
        who = pb.get("label") or "a run"
        if pb.get("kind") == "range":
            return "m/z: %s averaged %.2f to %.2f min" % (who, pb["ms_t0"], pb["ms_t1"]), None
        if pb.get("kind") == "scan":
            return "m/z: scans of %s at %.2f min" % (who, pb["ms_apex"]), None
        return "m/z: peak of %s at %.2f min" % (who, pb["ms_apex"]), None
    if sa.get("guides") != sb.get("guides"):
        ga, gb = list(sa.get("guides") or []), list(sb.get("guides") or [])
        add = [g for g in gb if g not in ga]
        rem = [g for g in ga if g not in gb]
        if add and not rem:
            return "guide line at %.2f min" % add[0], None
        if rem and not gb and len(rem) > 1:
            return "remove every guide line", None
        if rem and not add:
            return "remove the guide line at %.2f min" % rem[0], None
        return "guide lines", None
    if sa.get("align") != sb.get("align"):
        return ("align on %.2f min" % sb["align"] if sb.get("align") is not None else "undo the alignment"), None
    if (sa.get("t0"), sa.get("t1")) != (sb.get("t0"), sb.get("t1")):
        if sb.get("t0") is None and sb.get("t1") is None:
            return "the whole run (no time window)", None
        return "time window %s to %s min" % (_g10(sb.get("t0")) or "start", _g10(sb.get("t1")) or "end"), None
    if sa.get("ref_t") != sb.get("ref_t") or (sa.get("scale") != sb.get("scale") and sb.get("scale") == "ref"):
        if sb.get("ref_t") is not None:
            return "scale to the peak at %.2f min" % sb["ref_t"], None
    ea, eb = a["entries"], b["entries"]
    if [e[0] for e in ea] != [e[0] for e in eb]:
        if len(ea) == len(eb) and set(id(e[0]) for e in ea) == set(id(e[0]) for e in eb):
            return "order of the files", None
    da = {id(e[0]): e for e in ea}
    for e in eb:
        o = da.get(id(e[0]))
        if o is None:
            continue
        name = _cmp_label_of(fr, e[0]) if fr is not None else "a file"
        if o[1] != e[1]:
            return ("tick %s" if e[1] else "untick %s") % name, None
        if o[2] != e[2]:
            return ("label of %s: %s" % (name, e[2]) if e[2] else "label of %s: automatic" % name), None
        if o[3] != e[3]:
            return "colour of %s" % name, None
        if o[4] != e[4]:
            return "time shift of %s: %g min" % (name, e[4]), None
        if o[5] != e[5]:
            return ("background range of %s: %.2f to %.2f min" % ((name,) + tuple(e[5])) if e[5] else
                    "remove the background range of %s" % name), None
    if a["blank"] is not b["blank"]:
        return ("blank: %s" % _cmp_label_of(fr, b["blank"]) if b["blank"] is not None and fr is not None else
                "no blank"), None
    # the region areas after the files: a file moved or left out also moves the X values in their field
    aa, ab = a.get("area") or {}, b.get("area") or {}
    if aa.get("range") != ab.get("range"):
        rng = ab.get("range")
        return ("integrate %.4f to %.4f min in every trace" % tuple(rng) if rng else "clear the integrated region"), None
    if aa.get("reference") is not ab.get("reference"):
        ref = ab.get("reference")
        return ("100 %% area reference: %s" % _cmp_label_of(fr, ref) if ref is not None and fr is not None else
                "100 % area reference: automatic"), None
    if aa.get("x") != ab.get("x"):
        return "X values: %s" % (ab.get("x") or "none"), None
    if not _eq(aa, ab):
        return "area plot settings", None
    if not _eq(a["props"], b["props"]):
        return "graph properties", None
    ch = _cmp_choices()
    for k in list(_CMP_NAMES) + sorted(k for k in set(sa) | set(sb) if k not in _CMP_NAMES):
        if k in ("guides", "align", "t0", "t1", "ref_t"):
            continue
        if not _eq(sa.get(k), sb.get(k)):
            v = sb.get(k)
            if k in ch:
                v = next((lab for key, lab in ch[k] if key == v), v)
            elif isinstance(v, bool):
                v = "on" if v else "off"
            elif k == "polarity":
                v = "positive" if v == "+" else "negative"
            elif v is None:
                v = "none"
            u = _CMP_UNITS.get(k, "")
            txt = "%g" % v if isinstance(v, float) else str(v)
            return "%s: %s%s" % (_CMP_NAMES.get(k, k), txt, (" " + u) if u else ""), "cmp:%s" % k
    if a["spec_style"] != b["spec_style"]:
        return "spectra as %s" % ("sticks" if b["spec_style"] == 1 else "profile"), None
    for pol in ("+", "-"):
        ca, cb = (a.get("cols") or {}).get(pol), (b.get("cols") or {}).get(pol)
        if ca is None or cb is None:
            continue
        nm = "ESI+" if pol == "+" else "ESI−"
        if ca["pins"] != cb["pins"]:
            add = [m for m in cb["pins"] if m not in ca["pins"]]
            rem = [m for m in ca["pins"] if m not in cb["pins"]]
            if add:
                return "label m/z %.2f (%s)" % (add[0], nm), None
            if rem:
                return ("remove every label (%s)" % nm if not cb["pins"] and len(rem) > 1 else
                        "remove the label of m/z %.2f (%s)" % (rem[0], nm)), None
        if not _eq((ca["measure"], ca["m1"]), (cb["measure"], cb["m1"])):
            m = cb["measure"]
            if m and m[1] is not None:
                return "measure m/z %.2f to %.2f (%s)" % (m[0], m[1], nm), None
            if m:
                return "measure from m/z %.2f (%s)" % (m[0], nm), None
            return "remove the measurement (%s)" % nm, None
    if not _eq(pa, pb):
        return "m/z: no peak picked", None
    return "Compare view", None
