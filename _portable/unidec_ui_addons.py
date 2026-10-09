"""
User-interface add-ons for the portable UniDec package.

1. Right-click menu on every plot:
     Save plot as... (PNG/TIFF at 300 dpi, PDF, SVG, EPS)
     Copy plot to clipboard
     Export plot data as text...
     Reset zoom
     Clear labels
     plus the plot's original right-click action as an explicit item
     ("Integrate peaks" on the mass plot, "Set m/z range to current view" on
     the m/z plot, ...). These used to fire silently on any right-click.
   Ctrl + right-click keeps the original UniDec gestures (e.g. Ctrl + double
   right-click on the m/z plot removes the zoomed region). Plots that use the
   right button for zooming (UniChrom chromatogram) keep their behaviour.

2. Persistent labels: labels added from the peak list (masses, areas, names,
   differences, charge states) stay when the plot is redrawn (integration,
   Ignore/Isolate/Repopulate, colours, zoom). They are removed only by
   "Clear labels" (plot menu or peak-list menu), by choosing another label
   type (replaces the previous set), or when new results are produced
   (opening a file, processing, running UniDec, peak picking).

3. Movable labels: drag a label with the left mouse button to move it (a thin
   leader line then connects it to its peak); double-click it to edit the
   text; right-click it for Edit text / Reset position / Delete. Moved
   positions are kept through redraws, saved figures and re-labelling of the
   same peaks.

UniDec's own code is not modified; everything is patched in at start-up.
"""
import itertools
import os
import sys

import wx

_STATE = {"dir": None, "filter": 0}

_FRIENDLY = {  # UniDec main window plot names used for default file names
    "plot1": "ms_data_fit", "plot2": "mass_distribution", "plot3": "mz_grid",
    "plot4": "individual_peaks", "plot5": "mass_vs_charge", "plot6": "peak_intensities",
    "plot7": "plot7", "plot8": "plot8", "plot9": "mz_cube", "plot10": "mass_cube",
}

_WILDCARD = ("PNG image, 300 dpi (*.png)|*.png|"
             "PDF vector (*.pdf)|*.pdf|"
             "SVG vector (*.svg)|*.svg|"
             "EPS vector (*.eps)|*.eps|"
             "TIFF image, 300 dpi (*.tif)|*.tif")
_EXTS = [".png", ".pdf", ".svg", ".eps", ".tif"]

# Windows whose default right-click posts a selection event that they use
_SELECTION_WINDOWS = ("CDWindow", "ChromWindow")


def _log(*args):
    try:
        print("[add-on]", *args)
    except Exception:
        pass


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _top_and_name(plot):
    """Return (top level window, attribute name of this plot in it)."""
    try:
        top = plot.GetTopLevelParent()
    except Exception:
        return None, None
    name = None
    try:
        for k, v in vars(top).items():
            if v is plot:
                name = k
                break
    except Exception:
        pass
    return top, name


def _data_path(top):
    """Directory and base name of the open data file, if any."""
    try:
        cfg = top.pres.eng.config
        fname = getattr(cfg, "filename", "") or ""
        if fname:
            d = getattr(cfg, "dirname", "") or os.path.dirname(fname)
            return d, os.path.splitext(os.path.basename(fname))[0]
    except Exception:
        pass
    return None, None


def _plot_labels(plot):
    return len(getattr(plot, "text", []) or []) + len(getattr(plot, "labeltexts", []) or [])


def _all_plots(window):
    from unidec.modules.plotting.PlotBase import PlotBase
    plots = []
    try:
        for v in vars(window).values():
            if isinstance(v, PlotBase):
                plots.append(v)
    except Exception:
        pass
    return plots


def _clear_specs(plot, forget_positions=True):
    plot.__dict__["_udp_label_specs"] = []
    plot.__dict__["_udp_pending"] = False
    if forget_positions:
        plot.__dict__["_udp_prev_off"] = {}


# --------------------------------------------------------------------------
# 1. right-click menu
# --------------------------------------------------------------------------
def _save_plot(plot):
    top, name = _top_and_name(plot)
    ddir, base = _data_path(top)
    label = _FRIENDLY.get(name, name or "plot")
    default = ("%s_%s" % (base, label)) if base else label
    start_dir = _STATE["dir"] or ddir or os.path.expanduser("~")
    dlg = wx.FileDialog(plot, "Save plot as", start_dir, default + _EXTS[_STATE["filter"]], _WILDCARD,
                        wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
    dlg.SetFilterIndex(_STATE["filter"])
    try:
        if dlg.ShowModal() != wx.ID_OK:
            return
        path = dlg.GetPath()
        idx = dlg.GetFilterIndex()
    finally:
        dlg.Destroy()
    ext = os.path.splitext(path)[1].lower()
    if ext not in _EXTS + [".tiff", ".jpg", ".jpeg"]:
        path += _EXTS[idx]
        ext = _EXTS[idx]
    _STATE["dir"], _STATE["filter"] = os.path.dirname(path), idx
    kwargs = {"transparent": False}
    if ext in (".png", ".tif", ".tiff", ".jpg", ".jpeg"):
        kwargs["dpi"] = 300
    try:
        plot.figure.savefig(path, facecolor="white", bbox_inches="tight", pad_inches=0.05, **kwargs)
        _log("Saved plot:", path)
    except Exception as e:
        wx.MessageBox("Could not save the plot:\n%s" % e, "Save plot", wx.OK | wx.ICON_ERROR)


def _export_data(plot):
    top, name = _top_and_name(plot)
    ddir, base = _data_path(top)
    label = _FRIENDLY.get(name, name or "plot")
    default = (("%s_%s" % (base, label)) if base else label) + "_data.txt"
    start_dir = _STATE["dir"] or ddir or os.path.expanduser("~")
    dlg = wx.FileDialog(plot, "Export plot data as text", start_dir, default,
                        "Text file (*.txt)|*.txt|All files (*.*)|*.*", wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
    try:
        if dlg.ShowModal() != wx.ID_OK:
            return
        path = dlg.GetPath()
    finally:
        dlg.Destroy()
    _STATE["dir"] = os.path.dirname(path)
    try:
        plot.write_data(path)
        _log("Exported plot data:", path)
    except Exception as e:
        wx.MessageBox("Could not export the data:\n%s" % e, "Export plot data", wx.OK | wx.ICON_ERROR)


def _reset_zoom(plot):
    try:
        zoom = getattr(plot, "zoom", None)
        if zoom is not None and hasattr(zoom, "zoomout"):
            try:
                zoom.zoomout()
            except TypeError:
                zoom.zoomout(None)
            return
    except Exception:
        pass
    try:
        plot.repaint(resetzoom=True)
    except Exception as e:
        _log("reset zoom failed:", e)


def _clear_labels(plot):
    try:
        if getattr(plot, "zoom", None) is not None:
            plot.zoom.lflag = 0
            plot.zoom.kill_labels(repaint=False)
    except Exception:
        pass
    try:
        plot.textremove()  # also clears the stored (persistent) labels
    except Exception as e:
        _log("clear labels failed:", e)
    _clear_specs(plot)


class _FakeEvent(object):
    """Minimal stand-in for the matplotlib event UniDec's handler expects."""

    def __init__(self, x=None, y=None):
        self.dblclick = False
        self.button = 3
        self.xdata = x
        self.ydata = y
        self.inaxes = None


def _original_action(plot):
    """(label, callable) for the plot's original right-click action, or None."""
    orig = _ORIG.get("on_right_click")
    if orig is None:
        return None
    if getattr(plot, "int", 0) == 1:
        return "Integrate peaks", lambda: orig(plot, _FakeEvent())
    smash = getattr(plot, "smash", 0)
    if smash == 1:
        return "Set m/z range to current view", lambda: orig(plot, _FakeEvent())
    if smash == 2:
        return "Use current view as m/z and charge selection", lambda: orig(plot, _FakeEvent())
    top, _ = _top_and_name(plot)
    if top is not None and top.__class__.__name__ in _SELECTION_WINDOWS:
        return "Use current view as selection", lambda: orig(plot, _FakeEvent())
    return None


def _style_module():
    try:
        import unidec_plotstyle
        return unidec_plotstyle
    except Exception:
        return None


def _toggle_style(plot):
    ps = _style_module()
    if ps is None:
        return
    ps.set_modern(not ps.is_modern())
    top, _ = _top_and_name(plot)
    redrawn = False
    try:
        redrawn = ps.restyle_window(top)
    except Exception as e:
        _log("restyle failed:", e)
    if not redrawn:
        wx.MessageBox("The %s plot style applies to plots drawn from now on." %
                      ("modern" if ps.is_modern() else "classic"), "Plot style",
                      wx.OK | wx.ICON_INFORMATION)


def _show_menu(plot, spec=None):
    try:
        if not plot:  # window destroyed meanwhile
            return
    except Exception:
        return
    menu = wx.Menu()
    actions = {}

    def add(label, func, enabled=True, check=None):
        if check is None:
            item = menu.Append(wx.ID_ANY, label)
        else:
            item = menu.AppendCheckItem(wx.ID_ANY, label)
            item.Check(bool(check))
        item.Enable(enabled)
        if func is not None:
            actions[item.GetId()] = func

    if spec is not None:
        first = (str(spec.get("txt", "")).splitlines() or [""])[0]
        if len(first) > 30:
            first = first[:29] + "..."
        add('Label "%s"' % first, None, enabled=False)
        add("Edit label text...", lambda: _edit_label(plot, spec))
        add("Reset label position", lambda: _reset_label(plot, spec), bool(spec.get("off")))
        add("Delete this label", lambda: _delete_label(plot, spec))
        menu.AppendSeparator()

    add("Save plot as...", lambda: _save_plot(plot))
    add("Copy plot to clipboard", lambda: plot.copy_to_clipboard())
    add("Export plot data as text...", lambda: _export_data(plot),
        getattr(plot, "data", None) is not None)
    menu.AppendSeparator()
    add("Reset zoom", lambda: _reset_zoom(plot))
    add("Clear labels", lambda: _clear_labels(plot), _plot_labels(plot) > 0)
    orig = _original_action(plot)
    if orig is not None:
        menu.AppendSeparator()
        add(orig[0], orig[1])
    ps = _style_module()
    if ps is not None:
        menu.AppendSeparator()
        add("Classic plot style", lambda: _toggle_style(plot), check=not ps.is_modern())

    try:
        wx.ToolTip.Enable(False)
        selected = plot.GetPopupMenuSelectionFromUser(menu)
    finally:
        wx.ToolTip.Enable(True)
        menu.Destroy()
    func = actions.get(selected)
    if func is not None:
        try:
            func()
        except Exception as e:
            _log("menu action failed:", e)


def _uses_right_button_zoom(plot):
    try:
        from unidec.modules.plotting.NoZoomSpan import NoZoomSpan
        z = getattr(plot, "zoom", None)
        return isinstance(z, NoZoomSpan) and getattr(z, "zoombutton", None) == 3
    except Exception:
        return False


def _ctrl_down():
    try:
        return wx.GetKeyState(wx.WXK_CONTROL)
    except Exception:
        return False


def _on_press(plot, event):
    button = getattr(event, "button", None)
    if button == 1:
        try:
            _drag_press(plot, event)
        except Exception as e:
            _log("label drag:", e)
        return
    if button != 3:
        return
    if _ctrl_down() or _uses_right_button_zoom(plot):
        return  # original UniDec behaviour
    spec = None
    try:
        art = _hit_label(plot, event)
        if art is not None:
            spec = _spec_of(plot, art)
    except Exception:
        spec = None
    wx.CallAfter(_show_menu, plot, spec)


_ORIG = {}


def _install_plot_menu():
    from unidec.modules.plotting import PlottingWindow as pw
    base = pw.PlottingWindowBase

    orig_init = base.__init__
    orig_right = base.on_right_click
    _ORIG["on_right_click"] = orig_right

    def __init__(self, *args, **kwargs):
        orig_init(self, *args, **kwargs)
        try:
            # connected before any zoom handler, so a press on a label can
            # take the canvas lock before the zoom box sees it
            self.canvas.mpl_connect("button_press_event", lambda ev, p=self: _on_press(p, ev))
            self.canvas.mpl_connect("motion_notify_event", lambda ev, p=self: _on_motion(p, ev))
            self.canvas.mpl_connect("button_release_event", lambda ev, p=self: _on_release(p, ev))
        except Exception as e:
            _log("could not attach plot menu:", e)

    def on_right_click(self, event=None):
        # the hidden actions now live in the menu; keep them for Ctrl+click
        # and for plots that zoom with the right button
        if _ctrl_down() or _uses_right_button_zoom(self):
            return orig_right(self, event)
        return None

    base.__init__ = __init__
    base.on_right_click = on_right_click


# --------------------------------------------------------------------------
# 2. persistent labels
# --------------------------------------------------------------------------
def _install_persistent_labels():
    from unidec.modules.plotting.PlotBase import PlotBase
    from unidec.modules.plotting import PlottingWindow as pw

    orig_addtext = PlotBase.addtext
    orig_textremove = PlotBase.textremove
    orig_clear = PlotBase.clear_plot
    _ORIG["addtext"] = orig_addtext

    def addtext(self, txt, x, y, *args, **kwargs):
        if self.__dict__.get("_udp_reapplying", False):
            return orig_addtext(self, txt, x, y, *args, **kwargs)
        n_text = len(getattr(self, "text", []) or [])
        n_lines = len(getattr(self, "lines", []) or [])
        result = orig_addtext(self, txt, x, y, *args, **kwargs)
        try:
            spec = _new_spec(txt, x, y, args, kwargs)
            prev = (self.__dict__.get("_udp_prev_off") or {}).get(_spec_key(txt, x))
            if prev:
                spec["off"] = tuple(prev)
            self.__dict__.setdefault("_udp_label_specs", []).append(spec)
            _bind_artists(self, spec, n_text, n_lines)
            if prev and not kwargs.get("nopaint", False) and getattr(self, "canvas", None) is not None:
                self.canvas.draw_idle()
        except Exception as e:
            _log("label bookkeeping:", e)
        return result

    def textremove(self, *args, **kwargs):
        # remember moved positions, so re-labelling the same peaks keeps them
        try:
            prev = {}
            for s in self.__dict__.get("_udp_label_specs") or []:
                if s.get("off"):
                    prev[_spec_key(s["txt"], s["x"])] = s["off"]
            self.__dict__["_udp_prev_off"] = prev
        except Exception:
            pass
        _clear_specs(self, forget_positions=False)
        _cancel_drag(self)
        return orig_textremove(self, *args, **kwargs)

    def clear_plot(self, *args, **kwargs):
        if self.__dict__.get("_udp_label_specs"):
            self.__dict__["_udp_pending"] = True
        return orig_clear(self, *args, **kwargs)

    def reapply(self):
        if not self.__dict__.get("_udp_pending") or self.__dict__.get("_udp_reapplying"):
            return
        specs = self.__dict__.get("_udp_label_specs") or []
        ax = getattr(self, "subplot1", None)
        try:
            if not specs or ax is None or ax not in self.figure.axes:
                return  # plot not rebuilt yet; try again at the next repaint
        except Exception:
            return
        self.__dict__["_udp_reapplying"] = True
        try:
            for spec in list(specs):
                kw = dict(spec["kw"])
                kw["nopaint"] = True
                n_text = len(getattr(self, "text", []) or [])
                n_lines = len(getattr(self, "lines", []) or [])
                try:
                    orig_addtext(self, spec["txt"], spec["x"], spec["y"], *spec["args"], **kw)
                    _bind_artists(self, spec, n_text, n_lines)
                except Exception as e:
                    _log("could not restore label", spec.get("txt"), e)
        finally:
            self.__dict__["_udp_reapplying"] = False
            self.__dict__["_udp_pending"] = False

    def wrap_repaint(cls):
        orig = cls.__dict__.get("repaint")
        if orig is None:
            return

        def repaint(self, *args, **kwargs):
            try:
                reapply(self)
            except Exception as e:
                _log("label restore failed:", e)
            return orig(self, *args, **kwargs)

        cls.repaint = repaint

    PlotBase.addtext = addtext
    PlotBase.textremove = textremove
    PlotBase.clear_plot = clear_plot
    wrap_repaint(PlotBase)
    wrap_repaint(pw.PlottingWindowBase)

    # New results invalidate labels (masses change)
    def wrap_invalidate(cls, names):
        for n in names:
            orig = cls.__dict__.get(n)
            if orig is None:
                continue

            def make(orig=orig):
                def f(self, *args, **kwargs):
                    try:
                        for p in _all_plots(getattr(self, "view", None)):
                            _clear_specs(p)
                    except Exception:
                        pass
                    return orig(self, *args, **kwargs)
                f.__name__ = orig.__name__
                return f

            setattr(cls, n, make())

    targets = [("unidec.GUniDec", "UniDecApp"), ("unidec.MetaUniDec", "MetaUniDecBase"),
               ("unidec.UniDecCD", "UniDecCDApp")]
    names = ["on_open_file", "on_dataprep_button", "on_unidec_button", "on_pick_peaks", "on_auto"]
    import importlib
    for mod, cls in targets:
        try:
            wrap_invalidate(getattr(importlib.import_module(mod), cls), names)
        except Exception as e:
            _log("label reset hook not installed for", mod, e)

    # Ignore/Isolate in the peak list: drop the labels of peaks that are hidden
    try:
        gmod = importlib.import_module("unidec.GUniDec")
        orig_delete = gmod.UniDecApp.on_delete

        def on_delete(self, *args, **kwargs):
            try:
                _drop_ignored_labels(self)
            except Exception as e:
                _log("ignored-peak labels:", e)
            return orig_delete(self, *args, **kwargs)

        gmod.UniDecApp.on_delete = on_delete
    except Exception as e:
        _log("ignore hook not installed:", e)

    # "Clear Labels" at the end of the peak-list right-click menu
    try:
        from unidec.modules.gui_elements import peaklistsort
        panel_cls = peaklistsort.PeakListCtrlPanel

        def PopupMenu(self, menu, *args, **kwargs):
            try:
                if not hasattr(self, "_udp_clear_id"):
                    self._udp_clear_id = wx.NewIdRef()
                    self.Bind(wx.EVT_MENU, lambda e, s=self: _clear_window_labels(s), id=self._udp_clear_id)
                menu.AppendSeparator()
                menu.Append(self._udp_clear_id, "Clear Labels")
            except Exception as e:
                _log("peak list menu:", e)
            try:
                wx.ToolTip.Enable(False)
                return wx.Panel.PopupMenu(self, menu, *args, **kwargs)
            finally:
                wx.ToolTip.Enable(True)

        panel_cls.PopupMenu = PopupMenu
    except Exception as e:
        _log("peak list Clear Labels not installed:", e)


def _drop_ignored_labels(pres):
    masses = [float(p.mass) for p in pres.eng.pks.peaks if getattr(p, "ignore", 0)]
    if not masses:
        return
    for plot in _all_plots(pres.view):
        specs = plot.__dict__.get("_udp_label_specs")
        if not specs:
            continue
        keep = []
        for spec in specs:
            try:
                xf = float(spec["x"])
            except Exception:
                keep.append(spec)
                continue
            if any(abs(xf - m) <= 1e-6 * max(1.0, abs(m)) for m in masses):
                continue
            keep.append(spec)
        plot.__dict__["_udp_label_specs"] = keep


def _clear_window_labels(panel):
    try:
        top = panel.GetTopLevelParent()
    except Exception:
        return
    for p in _all_plots(top):
        if _plot_labels(p) or p.__dict__.get("_udp_label_specs"):
            _clear_labels(p)


# --------------------------------------------------------------------------
# 3. movable / editable labels
# --------------------------------------------------------------------------
_IDS = itertools.count(1)
_LEADER = dict(arrowstyle="-", color="#8a8a86", lw=0.8, shrinkA=0, shrinkB=3)
_BOX = dict(boxstyle="round,pad=0.25", fc="white", ec="#c9c9c4", lw=0.6, alpha=0.92)
_LEADER_MIN_PT = 9.0  # no leader line for labels moved less than this


def _new_spec(txt, x, y, args, kwargs):
    kw = dict(kwargs)
    kw.pop("nopaint", None)
    return {"id": next(_IDS), "txt": txt, "x": x, "y": y, "args": tuple(args), "kw": kw, "off": None}


def _spec_key(txt, x):
    try:
        return (str(txt), round(float(x), 4))
    except Exception:
        return (str(txt), str(x))


def _bind_artists(plot, spec, n_text, n_lines):
    """Tag the artists addtext just created; apply a stored offset."""
    texts = (getattr(plot, "text", None) or [])[n_text:]
    for a in texts:
        a._udp_id = spec["id"]
    for ln in (getattr(plot, "lines", None) or [])[n_lines:]:
        try:
            ln._udp_id = spec["id"]
        except Exception:
            pass
    if spec.get("off") and texts:
        _make_moved(plot, texts[-1], spec)


def _specs(plot):
    return plot.__dict__.get("_udp_label_specs") or []


def _spec_of(plot, artist):
    sid = getattr(artist, "_udp_id", None)
    for s in _specs(plot):
        if s["id"] == sid:
            return s
    return None


def _artist_of(plot, spec):
    for a in getattr(plot, "text", None) or []:
        if getattr(a, "_udp_id", None) == spec["id"]:
            return a
    return None


def _renderer(plot):
    try:
        return plot.canvas.get_renderer()
    except Exception:
        try:
            return plot.figure._get_renderer()
        except Exception:
            return None


def _text_bbox(a, renderer):
    from matplotlib.text import Text, Annotation
    if isinstance(a, Annotation):
        if not a._check_xy(renderer):
            return None
        a.update_positions(renderer)
    return Text.get_window_extent(a, renderer)


def _hit_label(plot, event, pad=4):
    if getattr(event, "x", None) is None or not _specs(plot):
        return None
    arts = [a for a in (getattr(plot, "text", None) or []) if getattr(a, "_udp_id", None) is not None]
    if not arts:
        return None
    r = _renderer(plot)
    if r is None:
        return None
    for a in reversed(arts):
        try:
            if a.axes is None or not a.get_visible():
                continue
            bb = _text_bbox(a, r)
            if bb is None:
                continue
            if bb.x0 - pad <= event.x <= bb.x1 + pad and bb.y0 - pad <= event.y <= bb.y1 + pad:
                return a
        except Exception:
            continue
    return None


def _has_vline(spec):
    a = spec.get("args") or ()
    if len(a) > 0:
        return bool(a[0])
    return bool(spec.get("kw", {}).get("vlines", True))


def _anchor(ax, spec, x, y):
    """Point the leader line goes to: the peak top under the label (mass
    labels) or the top of the label's dashed line (charge states)."""
    import numpy as np
    if _has_vline(spec):
        return x, (y * 0.95 if y > 0 else y)
    for ln in ax.lines:
        try:
            if ln.get_linestyle() in ("None", "none", "", " "):
                continue
            xd = np.asarray(ln.get_xdata(), dtype=float)
            yd = np.asarray(ln.get_ydata(), dtype=float)
            if len(xd) < 3:
                continue
            if xd[0] > xd[-1]:
                xd, yd = xd[::-1], yd[::-1]
            if xd[0] <= x <= xd[-1]:
                return x, min(float(np.interp(x, xd, yd)), y)
            break
        except Exception:
            break
    return x, y


def _offset_transform(ax, dx, dy):
    from matplotlib.transforms import offset_copy
    return offset_copy(ax.transData, fig=ax.figure, x=dx, y=dy, units="points")


def _set_offset(ann, dx, dy):
    ann.anncoords = _offset_transform(ann.axes, dx, dy)
    try:
        ann.arrow_patch.set_visible((dx * dx + dy * dy) ** 0.5 > _LEADER_MIN_PT)
    except Exception:
        pass


def _make_moved(plot, t, spec):
    """Replace a plain label by an annotation drawn at an offset (points)
    from its original place, with a white box and a leader line."""
    from matplotlib.text import Annotation
    ax = t.axes
    if ax is None:
        return t
    dx, dy = spec.get("off") or (0.0, 0.0)
    if isinstance(t, Annotation):
        _set_offset(t, dx, dy)
        return t
    x, y = t.get_position()
    xa, ya = _anchor(ax, spec, x, y)
    ann = ax.annotate(t.get_text(), xy=(xa, ya), xycoords="data", xytext=(x, y),
                      textcoords=_offset_transform(ax, dx, dy),
                      ha=t.get_horizontalalignment(), va=t.get_verticalalignment(),
                      color=t.get_color(), fontproperties=t.get_fontproperties(),
                      zorder=max(4, t.get_zorder() + 1), bbox=dict(_BOX), arrowprops=dict(_LEADER),
                      annotation_clip=True)
    ann._udp_id = spec["id"]
    ann._udp_plain = {"x": x, "y": y, "ha": t.get_horizontalalignment(), "va": t.get_verticalalignment(),
                      "color": t.get_color(), "fp": t.get_fontproperties(), "z": t.get_zorder()}
    ann.set_visible(t.get_visible())
    _set_offset(ann, dx, dy)
    try:
        plot.text[plot.text.index(t)] = ann
    except ValueError:
        plot.text.append(ann)
    try:
        t.remove()
    except Exception:
        pass
    return ann


def _make_plain(plot, ann, spec):
    p = getattr(ann, "_udp_plain", None)
    ax = ann.axes
    if p is None or ax is None:
        return ann
    t = ax.text(p["x"], p["y"], ann.get_text(), horizontalalignment=p["ha"], verticalalignment=p["va"],
                color=p["color"], fontproperties=p["fp"], zorder=p["z"])
    t._udp_id = spec["id"]
    try:
        plot.text[plot.text.index(ann)] = t
    except ValueError:
        plot.text.append(t)
    try:
        ann.remove()
    except Exception:
        pass
    return t


def _drag_allowed(plot):
    z = getattr(plot, "zoom", None)
    return z is None or z.__class__.__name__ == "ZoomBox"


def _modifier_down():
    try:
        return (wx.GetKeyState(wx.WXK_CONTROL) or wx.GetKeyState(wx.WXK_SHIFT) or
                wx.GetKeyState(wx.WXK_ALT))
    except Exception:
        return False


def _cancel_drag(plot):
    st = plot.__dict__.get("_udp_drag")
    plot.__dict__["_udp_drag"] = None
    if st is not None:
        try:
            plot.canvas.widgetlock.release(plot)
        except Exception:
            pass


def _drag_press(plot, event):
    _cancel_drag(plot)
    if event.inaxes is None or _modifier_down() or not _drag_allowed(plot):
        return
    art = _hit_label(plot, event)
    if art is None:
        return
    spec = _spec_of(plot, art)
    if spec is None:
        return
    lock = plot.canvas.widgetlock
    if not lock.available(plot):
        return
    lock(plot)  # the zoom box ignores events while the canvas is locked
    plot.__dict__["_udp_drag"] = {"a": art, "spec": spec, "x0": event.x, "y0": event.y,
                                  "off0": tuple(spec.get("off") or (0.0, 0.0)), "moved": False,
                                  "edit": bool(getattr(event, "dblclick", False))}


def _set_hover(plot, over):
    if over == plot.__dict__.get("_udp_hover", False):
        return
    plot.__dict__["_udp_hover"] = over
    try:
        from matplotlib.backend_tools import Cursors
        plot.canvas.set_cursor(Cursors.MOVE if over else Cursors.POINTER)
    except Exception:
        pass


def _on_motion(plot, event):
    try:
        st = plot.__dict__.get("_udp_drag")
        if st is None:
            if _specs(plot) or plot.__dict__.get("_udp_hover"):
                _set_hover(plot, event.inaxes is not None and _hit_label(plot, event) is not None)
            return
        if getattr(event, "x", None) is None:
            return
        ddx, ddy = event.x - st["x0"], event.y - st["y0"]
        if not st["moved"] and ddx * ddx + ddy * ddy < 9:
            return
        st["moved"] = True
        k = 72.0 / plot.figure.dpi
        off = (st["off0"][0] + ddx * k, st["off0"][1] + ddy * k)
        st["spec"]["off"] = off
        from matplotlib.text import Annotation
        if isinstance(st["a"], Annotation):
            _set_offset(st["a"], *off)
        else:
            st["a"] = _make_moved(plot, st["a"], st["spec"])
        plot.canvas.draw_idle()
    except Exception as e:
        _log("label drag:", e)
        _cancel_drag(plot)


def _on_release(plot, event):
    st = plot.__dict__.get("_udp_drag")
    if st is None:
        return
    _cancel_drag(plot)
    if st["moved"]:
        try:
            plot.canvas.draw_idle()
        except Exception:
            pass
    elif st["edit"]:
        wx.CallAfter(_edit_label, plot, st["spec"])


def _ask_label_text(parent, value):
    dlg = wx.Dialog(parent, title="Edit label", style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
    sizer = wx.BoxSizer(wx.VERTICAL)
    sizer.Add(wx.StaticText(dlg, label="Label text (empty deletes the label):"),
              0, wx.ALL, 8)
    txt = wx.TextCtrl(dlg, value=value, style=wx.TE_MULTILINE, size=dlg.FromDIP(wx.Size(340, 72)))
    sizer.Add(txt, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 8)
    buttons = dlg.CreateStdDialogButtonSizer(wx.OK | wx.CANCEL)
    sizer.Add(buttons, 0, wx.EXPAND | wx.ALL, 8)
    dlg.SetSizerAndFit(sizer)
    dlg.CentreOnParent()
    txt.SetFocus()
    txt.SelectAll()
    try:
        if dlg.ShowModal() != wx.ID_OK:
            return None
        return txt.GetValue()
    finally:
        dlg.Destroy()


def _edit_label(plot, spec):
    a = _artist_of(plot, spec)
    if a is None:
        return
    new = _ask_label_text(plot, a.get_text())
    if new is None:
        return
    new = new.replace("\r\n", "\n").rstrip("\n")
    if not new.strip():
        _delete_label(plot, spec)
        return
    spec["txt"] = new
    a.set_text(new)
    plot.canvas.draw_idle()


def _reset_label(plot, spec):
    spec["off"] = None
    a = _artist_of(plot, spec)
    from matplotlib.text import Annotation
    if isinstance(a, Annotation):
        _make_plain(plot, a, spec)
    plot.canvas.draw_idle()


def _delete_label(plot, spec):
    sid = spec["id"]
    for name in ("text", "lines"):
        lst = getattr(plot, name, None)
        if not lst:
            continue
        for a in [a for a in lst if getattr(a, "_udp_id", None) == sid]:
            try:
                a.remove()
            except Exception:
                pass
            lst.remove(a)
    try:
        _specs(plot).remove(spec)
    except ValueError:
        pass
    plot.canvas.draw_idle()


# --------------------------------------------------------------------------
# 4. taskbar identity (pinned UniDec icon groups with UniDec windows)
# --------------------------------------------------------------------------
APP_ID = "UniDec.Portable.Launcher"


def set_process_app_id():
    if os.name != "nt":
        return
    try:
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(APP_ID)
    except Exception as e:
        _log("AppUserModelID not set:", e)


def _set_window_relaunch(hwnd, exe_path):
    from win32com.propsys import propsys, pscon
    import pythoncom
    store = propsys.SHGetPropertyStoreForWindow(hwnd, propsys.IID_IPropertyStore)
    store.SetValue(pscon.PKEY_AppUserModel_ID, propsys.PROPVARIANTType(APP_ID, pythoncom.VT_LPWSTR))
    store.SetValue(pscon.PKEY_AppUserModel_RelaunchCommand,
                   propsys.PROPVARIANTType('"%s"' % exe_path, pythoncom.VT_LPWSTR))
    store.SetValue(pscon.PKEY_AppUserModel_RelaunchDisplayNameResource,
                   propsys.PROPVARIANTType("MS Analysis", pythoncom.VT_LPWSTR))
    store.SetValue(pscon.PKEY_AppUserModel_RelaunchIconResource,
                   propsys.PROPVARIANTType("%s,0" % exe_path, pythoncom.VT_LPWSTR))
    store.Commit()


def _install_window_relaunch(exe_path):
    if os.name != "nt" or not os.path.isfile(exe_path):
        return
    try:
        import win32com.propsys.propsys  # noqa: F401  (pywin32 available?)
    except Exception as e:
        _log("taskbar relaunch properties not available:", e)
        return
    orig_show = wx.Frame.Show
    done = set()

    def Show(self, show=True):
        result = orig_show(self, show)
        try:
            h = int(self.GetHandle())
            if show and h and h not in done:
                done.add(h)
                _set_window_relaunch(h, exe_path)
        except Exception as e:
            if not done:
                _log("taskbar relaunch properties not set:", e)
            done.add(-1)
        return result

    wx.Frame.Show = Show


# --------------------------------------------------------------------------
_RELAUNCH = {"done": False}


def install_relaunch(root=None, frames=()):
    """Taskbar pinning of every window (light: no UniDec modules). frames:
    windows already shown, which get the properties at once."""
    if _RELAUNCH["done"] or not root:
        return
    _RELAUNCH["done"] = True
    try:
        # MS Analysis.exe only (no relaunch properties without it): the old
        # UniDec.exe is not in the root any more and cannot start from elsewhere
        exe = os.path.join(root, "MS Analysis.exe")
        if not os.path.isfile(exe):
            return
        _install_window_relaunch(exe)
        for f in frames:
            try:
                _set_window_relaunch(int(f.GetHandle()), exe)
            except Exception:
                pass
    except Exception as e:
        _log("taskbar relaunch hook failed:", e)


def install(root=None):
    for name, func in (("plot menu", _install_plot_menu), ("persistent labels", _install_persistent_labels)):
        try:
            func()
            _log(name, "enabled")
        except Exception as e:
            _log(name, "NOT enabled:", e)
    install_relaunch(root)
