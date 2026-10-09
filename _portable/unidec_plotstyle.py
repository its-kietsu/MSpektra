"""
Modern plot style add-on for the portable UniDec package.

* Plots are drawn in a closed box (all four sides) with ticks and minor
  ticks outside the frame on the bottom and left axes.
* Intensity axes read "Relative intensity (%)" with 0-100 ticks (the data
  themselves are not changed; only the tick labels are scaled).
* Some head room above the tallest peak, so peaks and labels do not touch the
  frame.
* Calmer colours: dark grey data line, orange fit, a colour-blind friendly
  peak palette (used only while the peak colour map is UniDec's default
  "rainbow"; any other colour map chosen in UniDec is respected), white-edged
  peak markers, lighter integration shading.
* Axis numbers never switch to offset notation ("+1.83e4") when zoomed in.

"Classic UniDec style" in the plot right-click menu switches all of this off
(the choice is remembered). UniDec's own code is not modified.
"""
import os

_SETTINGS = {"path": None, "modern": True}
_ORIG = {}

INK = "#000000"  # axes, ticks and labels
TRACE = "#1D4ED8"  # data drawn in black by UniDec
FIT = "#eb6834"
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#8e5bd0", "#d4a020", "#e0457b", "#5b6b7c", "#9c6b43"]
HEADROOM = 0.08  # fraction of the y range added above the highest point


def _log(*args):
    try:
        print("[plot style]", *args)
    except Exception:
        pass


# --------------------------------------------------------------------------
# settings
# --------------------------------------------------------------------------
def _load():
    p = _SETTINGS["path"]
    if not p or not os.path.isfile(p):
        return
    try:
        import unidec_theme
        data = unidec_theme.read_settings(p)
        _SETTINGS["modern"] = data.get("plot_style", "modern") != "classic"
    except Exception as e:
        _log("settings not read:", e)


def _save():
    p = _SETTINGS["path"]
    if not p:
        return
    try:
        # as unidec_theme._save: the file holds every setting of the program, shared by its processes; it is
        # changed under their common lock and swapped in whole, and kept as it is when it cannot be read
        import unidec_theme
        unidec_theme.update_settings(p, {"plot_style": "modern" if _SETTINGS["modern"] else "classic"})
    except Exception as e:
        _log("settings not saved:", e)


def is_modern():
    return bool(_SETTINGS["modern"])


def set_modern(flag):
    _SETTINGS["modern"] = bool(flag)
    _save()


# --------------------------------------------------------------------------
# tick helpers
# --------------------------------------------------------------------------
def _ticker_classes():
    from matplotlib.ticker import Locator, MaxNLocator, Formatter

    class PercentLocator(Locator):
        """Nice 0/25/50/75/100 style ticks on a data axis whose 100 % is ref."""

        def __init__(self, ref):
            self.ref = float(ref)
            self._m = MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10])

        def __call__(self):
            vmin, vmax = self.axis.get_view_interval()
            return self.tick_values(vmin, vmax)

        def tick_values(self, vmin, vmax):
            r = self.ref
            if vmin > vmax:
                vmin, vmax = vmax, vmin
            ticks = self._m.tick_values(vmin / r * 100.0, vmax / r * 100.0)
            return [t * r / 100.0 for t in ticks]

    class PercentFormatter(Formatter):
        def __init__(self, ref, hide_negative=True):
            self.ref = float(ref)
            self.hide_negative = hide_negative

        def __call__(self, v, pos=None):
            p = v / self.ref * 100.0
            if abs(p) < 1e-7:
                p = 0.0
            if self.hide_negative and p < 0:
                return ""
            return ("%.3f" % p).rstrip("0").rstrip(".")

    return PercentLocator, PercentFormatter


# --------------------------------------------------------------------------
# style pass (runs before every repaint; idempotent)
# --------------------------------------------------------------------------
def _width_locator():
    """Round-number x ticks, at most one per 0.55 inch of axis width, so four
    and five digit labels do not run into each other on small plots."""
    from matplotlib.ticker import MaxNLocator

    class WidthLocator(MaxNLocator):
        def __call__(self):
            try:
                ax = self.axis.axes
                w_in = ax.get_position().width * ax.figure.get_size_inches()[0]
                self.set_params(nbins=int(max(3, min(10, w_in / 0.55))))
            except Exception:
                pass
            return MaxNLocator.__call__(self)

    return WidthLocator(nbins=6, steps=[1, 2, 5, 10], min_n_ticks=2)


def _is_color(c, rgba):
    try:
        from matplotlib.colors import to_rgba
        return all(abs(a - b) < 1e-3 for a, b in zip(to_rgba(c), rgba))
    except Exception:
        return False


def _style_box(ax, fig):
    from matplotlib.ticker import AutoMinorLocator
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(0.9)
        s.set_color(INK)
    # closed frame, but ticks only where the numbers are (bottom and left)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.tick_params(which="both", direction="out", top=False, right=False, bottom=True, left=True,
                   color=INK, width=0.9)
    ax.tick_params(which="major", length=5, labelsize=11, labelcolor=INK, pad=4)
    ax.tick_params(which="minor", length=2.5)
    try:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    except Exception:
        pass
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    ax.xaxis.label.set_fontsize(12)
    ax.yaxis.label.set_fontsize(12)
    # title: small, left aligned, grey
    t = ax.get_title()
    if t:
        ax.set_title("")
        ax.set_title(t, loc="left", fontsize=11, color="#555555", pad=6)
    # margins in inches (fit the labels at any plot size), updated on resize;
    # _place_axes itself runs at the end of the style functions
    try:
        if getattr(fig, "_udp_resize_cid", None) is None and fig.canvas is not None:
            fig._udp_resize_cid = fig.canvas.mpl_connect("resize_event", lambda e, f=fig: _on_resize(f))
    except Exception:
        pass


def _place_axes(ax):
    """Axes position from fixed margins in inches (labels never clipped,
    no wasted space on large plots)."""
    try:
        fig = ax.figure
        if len([a for a in fig.axes if a.get_visible()]) > 2:
            return  # multi-panel figures keep their own layout
        w, h = fig.get_size_inches()
        if w <= 0 or h <= 0:
            return
        meta = getattr(ax, "_udp_meta", {}) or {}
        left = 0.68 if (meta.get("relative") or meta.get("kind") == "bar") else 0.86
        if not ax.get_ylabel():
            left -= 0.2
        bottom = 0.54 if ax.get_xlabel() else 0.36
        if meta.get("kind") == "bar":
            bottom += 0.12  # rotated peak letters
        top = 0.30 if ax.get_title(loc="left") else 0.12
        right = 0.14
        x0, y0 = left / w, bottom / h
        ax.set_position([x0, y0, max(0.2, 1 - x0 - right / w), max(0.2, 1 - y0 - top / h)])
    except Exception:
        pass


def _on_resize(fig):
    if not is_modern():
        return
    for ax in fig.axes:
        if getattr(ax, "_udp_meta", None):
            _place_axes(ax)


def _style_line_axes(plot, ax, meta):
    from matplotlib.ticker import ScalarFormatter
    fig = ax.figure
    _style_box(ax, fig)

    # x axis: plain numbers (no "+1.83e4" offsets), round tick steps, italic m/z
    try:
        from matplotlib.ticker import AutoLocator
        fmt = ScalarFormatter(useOffset=False)
        fmt.set_scientific(False)
        ax.xaxis.set_major_formatter(fmt)
        if type(ax.xaxis.get_major_locator()) is AutoLocator:
            ax.xaxis.set_major_locator(_width_locator())
    except Exception:
        pass
    xl = ax.get_xlabel()
    if xl.strip().lower().startswith("m/z"):
        ax.set_xlabel(r"$m/z$")

    # y axis: relative intensity in %
    PercentLocator, PercentFormatter = _ticker_classes()
    ref = meta.get("ref") or 0
    if meta.get("relative") and ref > 0:
        ax.yaxis.set_major_locator(PercentLocator(ref))
        ax.yaxis.set_major_formatter(PercentFormatter(ref, hide_negative=True))
        ax.set_ylabel("Relative intensity (%)")
    elif meta.get("pub") and not ax.get_ylabel():
        ax.set_ylabel("Intensity")

    # lines
    for ln in ax.lines:
        if getattr(ln, "_udp_st", False):
            continue
        ln._udp_st = True
        try:
            ls = ln.get_linestyle()
            if ls in ("None", "none", "", " ") and ln.get_marker() not in (None, "None", "none", ""):
                if _is_color(ln.get_markeredgecolor(), (0, 0, 0, 1)):
                    ln.set_markeredgecolor("white")
                    ln.set_markeredgewidth(0.9)
                continue
            if _is_color(ln.get_color(), (0, 0, 0, 1)):
                ln.set_color(TRACE)
                ln.set_linewidth(min(float(ln.get_linewidth() or 1.0), 0.8))  # thinner data lines
            elif _is_color(ln.get_color(), (1, 0, 0, 1)) and ln.get_label() == "Fit Data":
                ln.set_color(FIT)
                ln.set_alpha(0.9)
        except Exception:
            pass

    # integration shading from UniDec (alpha 0.75, black edge) -> lighter
    for c in list(ax.collections):
        if getattr(c, "_udp_st", False):
            continue
        c._udp_st = True
        try:
            if c.__class__.__name__ in ("PolyCollection", "FillBetweenPolyCollection") and \
                    abs((c.get_alpha() or 0) - 0.75) < 1e-6:
                c.set_alpha(0.35)
                c.set_edgecolor("none")
        except Exception:
            pass

    # a faint fill under the zero-charge mass spectrum
    if meta.get("mass") and not getattr(ax, "_udp_fill", False) and ax.lines:
        try:
            import numpy as np
            main = ax.lines[0]
            x = np.asarray(main.get_xdata(), dtype=float)
            y = np.asarray(main.get_ydata(), dtype=float)
            if len(x) > 2 and np.all(np.isfinite(y)) and np.nanmin(y) >= 0:
                f = ax.fill_between(x, y, 0, color=PALETTE[0], alpha=0.07, linewidth=0, zorder=0.8)
                f._udp_st = True
            ax._udp_fill = True
        except Exception:
            ax._udp_fill = True

    _place_axes(ax)

    # legend without frame
    leg = ax.get_legend()
    if leg is not None and not getattr(leg, "_udp_st", False):
        leg._udp_st = True
        try:
            leg.get_frame().set_visible(False)
            handles = getattr(leg, "legend_handles", None) or getattr(leg, "legendHandles", [])
            by_label = {ln.get_label(): ln for ln in ax.lines}
            for h, t in zip(handles, leg.get_texts()):
                t.set_fontsize(10)
                t.set_color(INK)
                src = by_label.get(t.get_text())
                if src is not None and hasattr(h, "set_color"):
                    h.set_color(src.get_color())
        except Exception:
            pass


def _style_bar_axes(plot, ax, meta):
    _style_box(ax, ax.figure)
    for ln in ax.lines:
        if getattr(ln, "_udp_st", False):
            continue
        ln._udp_st = True
        try:
            if _is_color(ln.get_markeredgecolor(), (0, 0, 0, 1)):
                ln.set_markeredgecolor("white")
                ln.set_markeredgewidth(0.9)
        except Exception:
            pass
    for p in ax.patches:
        if getattr(p, "_udp_st", False):
            continue
        p._udp_st = True
        try:
            if p.__class__.__name__ == "Rectangle" and p.get_visible() and p.get_width() == 1:
                p.set_edgecolor("white")
                p.set_linewidth(0.6)
        except Exception:
            pass
    _place_axes(ax)


def style_plot(plot):
    """Apply the modern style to the plot's main axes (no drawing)."""
    if not is_modern():
        return
    ax = getattr(plot, "subplot1", None)
    if ax is None:
        return
    meta = getattr(ax, "_udp_meta", None)
    if not meta:
        return
    try:
        if ax not in plot.figure.axes:
            return
    except Exception:
        return
    # idempotent: artists already styled are skipped
    try:
        if meta.get("kind") == "line":
            _style_line_axes(plot, ax, meta)
        elif meta.get("kind") == "bar":
            _style_bar_axes(plot, ax, meta)
    except Exception as e:
        _log("style failed:", e)


# --------------------------------------------------------------------------
# patches
# --------------------------------------------------------------------------
def _install_plot_hooks():
    import numpy as np
    from unidec.modules.plotting.plot1d import Plot1dBase
    from unidec.modules.plotting.PlotBase import PlotBase

    orig_top = Plot1dBase.plotrefreshtop
    orig_bar = Plot1dBase.barplottop

    def plotrefreshtop(self, xvals, yvals, *args, **kwargs):
        # draw without painting, tag the axes, then paint the way the
        # original would (so the zoom set-up already sees the tag)
        args = list(args)
        if len(args) > 8:
            nopaint = bool(args[8])
            args[8] = True
        else:
            nopaint = bool(kwargs.get("nopaint", False))
            kwargs["nopaint"] = True
        zoomout = bool(args[11]) if len(args) > 11 else bool(kwargs.get("zoomout", False))
        result = orig_top(self, xvals, yvals, *args, **kwargs)
        try:
            ax = self.subplot1
            y = np.asarray(yvals, dtype=float)
            y = y[np.isfinite(y)]
            config = kwargs.get("config")
            if config is None and len(args) >= 5:
                config = args[4]
            datanorm = getattr(config, "datanorm", 1) if config is not None else 1
            pub = bool(getattr(config, "publicationmode", 0)) if config is not None else False
            ymin = float(np.amin(y)) if len(y) else 0.0
            ref = float(np.amax(y)) if len(y) else 0.0
            xl = ax.get_xlabel().strip().lower()
            ax._udp_meta = {
                "kind": "line",
                "ref": ref,
                "relative": (datanorm != 0 or pub) and ref > 0 and config is not None,
                "negative": ymin < 0,
                "pub": pub,
                "mass": xl.startswith("mass") and ymin >= 0,
            }
        except Exception as e:
            _log("meta failed:", e)
        if not nopaint:
            if zoomout:
                self.repaint(setupzoom=True, resetzoom=True)
            else:
                self.repaint(setupzoom=True)
        return result

    def barplottop(self, *args, **kwargs):
        repaint = kwargs.get("repaint", True)
        kwargs["repaint"] = False
        result = orig_bar(self, *args, **kwargs)
        try:
            self.subplot1._udp_meta = {"kind": "bar"}
        except Exception:
            pass
        if repaint:
            self.repaint()
        return result

    Plot1dBase.plotrefreshtop = plotrefreshtop
    Plot1dBase.barplottop = barplottop

    def wrap_repaint(cls):
        orig = cls.__dict__.get("repaint")
        if orig is None:
            return

        def repaint(self, *args, **kwargs):
            try:
                style_plot(self)
            except Exception as e:
                _log("style pass failed:", e)
            return orig(self, *args, **kwargs)

        cls.repaint = repaint

    wrap_repaint(PlotBase)
    try:
        from unidec.modules.plotting import PlottingWindow as pw
        wrap_repaint(pw.PlottingWindowBase)
    except Exception as e:
        _log("window repaint hook:", e)


def _headroom_ok(axes_list):
    if not is_modern():
        return False
    try:
        return all(getattr(a, "_udp_meta", None) for a in axes_list)
    except Exception:
        return False


def _install_headroom():
    from unidec.modules.plotting import ZoomBox as zb
    from unidec.modules.plotting import ZoomCommon as zc

    orig_init = zb.ZoomBox.__init__

    def __init__(self, *args, **kwargs):
        orig_init(self, *args, **kwargs)
        try:
            if _headroom_ok(self.axes):
                xmin, ymin, xmax, ymax = self.data_lims
                if ymax > ymin:
                    if 0 < ymin < 0.05 * (ymax - ymin):
                        ymin = 0.0  # start the intensity axis at 0 (shows the "0" tick)
                    self.data_lims = [xmin, ymin, xmax, ymax + HEADROOM * (ymax - ymin)]
                    self.initialize()
        except Exception as e:
            _log("head room:", e)

    orig_auto = zc.ZoomCommon.set_auto_ylim

    def set_auto_ylim(self, xmin, xmax, draw=True):
        orig_auto(self, xmin, xmax, draw=False)
        try:
            if _headroom_ok(self.axes):
                for axes in self.axes:
                    y0, y1 = axes.get_ylim()
                    if y1 > y0:
                        axes.set_ylim((y0, y1 + HEADROOM * (y1 - y0)))
                        self.yzoom = [y0, y1 + HEADROOM * (y1 - y0)]
        except Exception as e:
            _log("head room (zoom):", e)
        if draw:
            self.canvas.draw()

    zb.ZoomBox.__init__ = __init__
    zc.ZoomCommon.set_auto_ylim = set_auto_ylim

    # plots made without a window (HTML report, scripts)
    from unidec.modules.plotting.PlotBase import PlotBase
    orig_setup = PlotBase.setup_zoom

    def setup_zoom(self, *args, **kwargs):
        result = orig_setup(self, *args, **kwargs)
        try:
            ax = self.subplot1
            if _headroom_ok([ax]):
                y0, y1 = ax.get_ylim()
                if y1 > y0:
                    ax.set_ylim((y0, y1 + HEADROOM * (y1 - y0)))
        except Exception:
            pass
        return result

    PlotBase.setup_zoom = setup_zoom


def _install_palette():
    import numpy as np
    from matplotlib.colors import to_rgba
    from unidec.modules import peakstructure

    orig = peakstructure.Peaks.default_params
    rgba = [np.array(to_rgba(c)) for c in PALETTE]

    def default_params(self, cmap="rainbow", *args, **kwargs):
        result = orig(self, cmap, *args, **kwargs)
        try:
            name = cmap
            if isinstance(name, bytes):
                name = name.decode("utf-8", "ignore")
            name = str(name)
            if name.startswith("b'"):
                name = name[2:-1]
            if is_modern() and name == "rainbow":
                self.peakcolors = np.array([rgba[i % len(rgba)] for i in range(len(self.peaks))])
                for i, p in enumerate(self.peaks):
                    p.color = self.peakcolors[i]
        except Exception as e:
            _log("palette:", e)
        return result

    peakstructure.Peaks.default_params = default_params


def restyle_window(top):
    """Redraw a UniDec window after switching style (recolours peaks)."""
    pres = getattr(top, "pres", None)
    eng = getattr(pres, "eng", None) if pres is not None else None
    try:
        pks = getattr(eng, "pks", None)
        if pks is not None and getattr(pks, "plen", 0) > 0:
            pks.default_params(cmap=eng.config.peakcmap)
            panel = getattr(getattr(pres, "view", None), "peakpanel", None)
            if panel is not None:
                panel.add_data(pks)
    except Exception as e:
        _log("peak recolour:", e)
    try:
        data2 = getattr(getattr(eng, "data", None), "data2", None)
        has_data = data2 is not None and len(data2) > 1
    except Exception:
        has_data = False
    if pres is not None and hasattr(pres, "on_replot") and has_data:
        try:
            pres.on_replot()
            return True
        except Exception as e:
            _log("replot failed:", e)
    return False


def install(settings_path=None):
    _SETTINGS["path"] = settings_path
    _load()
    for name, func in (("plot hooks", _install_plot_hooks), ("head room", _install_headroom),
                       ("peak palette", _install_palette)):
        try:
            func()
        except Exception as e:
            _log(name, "NOT enabled:", e)
    _log("style:", "modern" if is_modern() else "classic")
