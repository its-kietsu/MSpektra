# -*- coding: utf-8 -*-
"""
Mass shifts on the zero charge mass spectrum of a deconvolution result (MS Analysis 3.35): the window part
of the mass shift finder (the calculation is ms_shifts.find_shifts, in the C++ library).

Right click a result > Mass shifts (or Mass shifts... in the Deconvolution section of the side panel):
brackets above the peaks of the zero charge mass spectrum join the masses that differ by a known shift
(Na, oxidation, water, carbamylation, ...) or by k tags (a dye or another label) plus at most two other
shifts, each labelled with its name and the difference; a mass whose difference from the reference
matches nothing is marked as unknown above its peak (e.g. "+429.1 unknown"). By default each mass is
joined to the reference or to its neighbouring species (the nearest mass it follows from); "Show every
matched pair" draws them all. Next to the plot, in one tile, a table lists the pairs (from, to,
difference, match, error, what the mass is) and, with a tag, a second table below it the degree of
conjugation (the share of the species with 0, 1, 2, ... tags by height and by area, and the tags per
molecule). The reference is chosen automatically (the tallest mass; with a tag, the lightest of a tag
series below it when the tallest is a conjugate); right click a peak > Use as the unmodified species
changes it. The shift list (name, average and monoisotopic difference, on or off), the tags, the tolerance
and the display are set in the Mass shifts window and kept in the settings ("mass_shifts").

Proteins (isotope envelopes) are compared by average masses, isotope resolved results by monoisotopic
differences (IsoDec masses; resolved species by the monoisotopic masses estimated from their most abundant
isotope with the averagine model, one isotope more or less allowed, shown by the isotope of the mass table;
the isotope peaks of a resolved spectrum are grouped first).

The state of the finder (settings and, per result, on or off and the reference) is a plain dict:
DeconvPanel.shift_state() and DeconvPanel.set_shift_state(d) (ShiftPanel.get_state / set_state).
"""
import math
import uuid

import numpy as np
import wx

import unidec_theme as T
from unidec_theme import C, ui_font
import unilcms as U
import ms_shifts as MS

KEY = "mass_shifts"
ROW_PT = 10.5  # height of a bracket with its label (points)
MAX_UNKNOWN = 8  # unknown masses marked in the view (the tallest ones; all of them in the table)
MAX_BRACKETS = 40  # brackets drawn at most (Show every matched pair; the others in the table)
COL_SHIFT, COL_TAG, COL_UNKNOWN, COL_GUIDE = "#0E7C5A", "#7B4FC4", "#C8324F", "#9AA3AE"


def settings():
    """The finder's settings (the program's settings file, made valid)."""
    try:
        return MS.clean_settings(U.settings().get(KEY))
    except Exception:
        return MS.default_settings()


def save_settings(s):
    T._save({KEY: MS.clean_settings(s)})
    MS.SETTINGS_REV[0] += 1  # (undo: the settings changed)


def _fmt_num(v, d):
    return ("%%.%df" % d) % v if v is not None and math.isfinite(v) else ""


def _key(res):
    """The name of a result for set_state: an id of its own (two runs of one method on one spectrum are saved
    under the same name; kept with the result)."""
    k = res.get("shift_id")
    if not isinstance(k, str) or not k:
        k = res["shift_id"] = uuid.uuid4().hex
    return k


_REDRAW = {"pending": False}


def _tallest(r):
    """Index of the tallest species of a result of ms_shifts.find_shifts (the first of equal heights), or -1."""
    best = -1
    for j, q in enumerate(r["species"]):
        if q["status"] != 4 and (best < 0 or q["height"] > r["species"][best]["height"]):
            best = j
    return best


class ShiftPanel(object):
    """The mass shift finder of a DeconvPanel: per result its state (ent.shift: on, reference), its two
    tables (ent.card_shift and ent.card_conj in the tile ent.shift_box, in the row of the result next to
    the mass spectrum), the brackets on the zero charge mass spectrum, the menu entries and the settings
    window."""

    def __init__(self, dec, sp=None):
        self.dec = dec
        self.btn = None
        if sp is not None:
            self.btn = U.flat(sp, "Mass shifts…", icon=None, handler=lambda e: self.open_dialog(),
                              tooltip="Mass differences, modifications and tags")
            sp.buttons(self.btn)

    # ----------------------------------------------------------- state
    @staticmethod
    def state(ent):
        st = getattr(ent, "shift", None)
        if not isinstance(st, dict):
            st = ent.shift = {"on": False, "ref": None}
        return st

    def is_on(self, ent):
        return ent is not None and bool(self.state(ent)["on"])

    def get_state(self):
        """Plain dict of everything the finder shows (for undo): the settings and per result whether it is
        on and its reference (a mass, or None: automatic)."""
        return {"settings": settings(),
                "results": [{"key": _key(e.res), "on": bool(self.state(e)["on"]), "ref": self.state(e)["ref"]}
                            for e in self.dec.results]}

    def set_state(self, d):
        """Back to a state of get_state (the settings are saved; every result drawn again)."""
        if not isinstance(d, dict):
            return
        if isinstance(d.get("settings"), dict):
            new = MS.clean_settings(d["settings"])
            if new != settings():
                save_settings(new)
        items = d.get("results") if isinstance(d.get("results"), list) else []
        by_key = {it.get("key"): it for it in items if isinstance(it, dict)}
        for n, ent in enumerate(self.dec.results):
            it = by_key.get(_key(ent.res))
            if it is None and n < len(items) and isinstance(items[n], dict) and "key" not in items[n]:
                it = items[n]
            st = self.state(ent)
            on = bool(it.get("on")) if it else False
            ref = it.get("ref") if it else None
            ref = float(ref) if isinstance(ref, (int, float)) and math.isfinite(ref) else None
            changed = on != st["on"]
            st["on"], st["ref"] = on, ref
            if changed:
                self.dec._apply_z(ent)
        self.request_redraw_all()  # (undo restores every file at once: one redraw)

    # ----------------------------------------------------------- tables in the row
    def make_cards(self, ent):
        """The tile of the finder in the row of a result: the table of the pairs above the degree of
        conjugation (shown with a tag), hidden until the finder is turned on for the result."""
        import deconv_tab as D
        ent.shift_box = box = U.SplitBox(ent.row, wx.VERTICAL)
        box.title = "Mass shifts"
        ent.card_shift = D.TableCard(box, "Mass shifts", self._cols("average"),
                                     on_select=lambda i, r=ent: self._row_selected(r, i),
                                     empty="Turn on with right click > Mass shifts")
        ent.card_conj = D.TableCard(box, "Degree of conjugation", self._conj_cols(),
                                    empty="Set a tag in the Mass shifts window")
        for card in (ent.card_shift, ent.card_conj):
            card.extra_menu = lambda r=ent, c=card: self._table_menu(r, c)
            card.Bind(wx.EVT_MOUSEWHEEL, lambda e, c=card: U.forward_wheel(c, e))
        box.set_panes([ent.card_shift], [1.0])
        ent.card_conj.Hide()  # (a pane left out of a SplitBox is not hidden by it)
        ent._shift_filled = ent._shift_cache = None  # new tables (a row made again: undo of Close this result)
        box.on_close_tile = lambda pane, r=ent: self.close_tile(r, pane)
        box.Hide()
        ent.card_shift._shift_rows = []

    def panes(self, ent):
        """(panes, weights) of the finder in the row of ent: none while it is off."""
        box = getattr(ent, "shift_box", None)
        if box is None:
            return [], []
        if not self.is_on(ent):
            try:
                box.Hide()
            except RuntimeError:
                pass
            return [], []
        try:
            if settings()["tags"]:
                w = [box.weight_of(ent.card_shift, 2.0), box.weight_of(ent.card_conj, 1.25)]
                box.set_panes([ent.card_shift, ent.card_conj], w if len(box.panes) == 2 else [2.0, 1.25])
            else:
                box.set_panes([ent.card_shift], [1.0])
                ent.card_conj.Hide()
        except RuntimeError:
            return [], []
        return [box], [3.0]

    def close_tile(self, ent, pane):
        """Close this tile on a table: the finder is turned off for the result."""
        if pane in (getattr(ent, "card_shift", None), getattr(ent, "card_conj", None), getattr(ent, "shift_box", None)):
            self.set_on(ent, False)
            return True
        return False

    # ----------------------------------------------------------- turning on and off
    def set_on(self, ent, flag):
        if ent is None or ent not in self.dec.results:
            return
        st = self.state(ent)
        flag = bool(flag)
        if st["on"] == flag:
            return
        st["on"] = flag
        self.dec._apply_z(ent)
        self.refresh(ent)
        if flag:
            v = self.view(ent)
            n_unk = sum(1 for q in v["r"]["species"] if q["status"] == 3) if v else 0
            self.dec.status("Mass shifts: %s; %d unknown difference%s from the reference" % (
                self._info(v) if v else "no masses", n_unk, "" if n_unk == 1 else "s"))
        else:
            self.dec.status("Mass shifts hidden")

    def set_ref(self, ent, i):
        """Peak i of the result as the unmodified species (None: automatic)."""
        if ent is None:
            return
        st = self.state(ent)
        peaks = ent.res.get("peaks") or []
        if i is None or not (0 <= i < len(peaks)):
            st["ref"] = None
        else:
            st["ref"] = float(peaks[i]["mass"])
        if not st["on"]:
            st["on"] = True
            self.dec._apply_z(ent)
        self.refresh(ent)
        v = self.view(ent)
        if v:
            self.dec.status("Mass shifts relative to %s Da%s" % (self._mass_text(v, v["r"]["ref"]),
                                                                 self._ref_how(ent, v)))

    def _ref_index(self, ent):
        ref = self.state(ent)["ref"]
        if ref is None:
            return -1
        peaks = ent.res.get("peaks") or []
        for i, q in enumerate(peaks):
            if abs(float(q["mass"]) - ref) <= 1e-6 * max(1.0, abs(ref)):
                return i
        return -1  # the mass is gone (removed by hand): automatic

    # ----------------------------------------------------------- calculation
    def view(self, ent):
        """The finder's result for ent (cached while the masses and settings stay the same): {"r": the
        dict of ms_shifts.find_shifts, "basis", "names", "tags", "tol" (value, unit), "mass", "dm", "err"
        (decimals), ...}, or None (no masses, or the calculation failed)."""
        res = ent.res
        peaks = res.get("peaks") or []
        if not peaks:
            return None
        s = settings()
        b = MS.result_basis(res)
        col = "avg" if b["basis"] == "average" else "mono"
        on = [e for e in s["shifts"] if e["on"]]
        vals = [float(e[col]) for e in on]
        tags = [float(t[col]) for t in s["tags"]]
        tol = s["tol_envelope"] if b["kind"] == "envelope" else s["tol_resolved"]
        auto = tol is None
        tol = list(tol) if tol else list(b["tol"])
        td, tp = (float(tol[0]), 0.0) if tol[1] == "Da" else (0.0, float(tol[0]))
        ref = self._ref_index(ent)
        key = (tuple(b["masses"]), tuple(b["shown"]), tuple(b["heights"]), tuple(b["areas"]), tuple(vals), tuple(tags), td, tp, ref,
               s["kmax"], s["max_extra"], b["flags"], tuple(e["name"] for e in on), tuple(t["name"] for t in s["tags"]),
               s["show_all"], s["min_rel"])
        cache = getattr(ent, "_shift_cache", None)
        if cache is not None and cache[0] == key:
            return cache[1]
        try:
            r = MS.find_shifts(b["masses"], b["heights"], b["areas"], vals, tags, td, tp, ref, s["kmax"],
                               s["max_extra"], b["flags"], s["min_rel"])
        except Exception as ex:
            print("[mass shifts] %s" % ex)
            return None
        env = b["basis"] == "average"
        v = {"r": r, "basis": b["basis"], "kind": b["kind"], "names": [e["name"] for e in on],
             "tags": [t["name"] for t in s["tags"]], "tol": tol, "auto": auto, "masses": b["shown"],
             "dm": 1 if env else 4, "dd": 1 if env else 3, "de": 2 if env else 4, "show_all": s["show_all"],
             "kmax": s["kmax"], "min_rel": s["min_rel"], "auto_ref": ref < 0, "real_areas": b["real_areas"]}
        ent._shift_cache = (key, v)
        return v

    def _mass_text(self, v, i):
        return _fmt_num(v["masses"][i], v["dm"]) if 0 <= i < len(v["masses"]) else ""

    def _info(self, v):
        r = v["r"]
        what = "average masses" if v["basis"] == "average" else (
            "monoisotopic differences (estimated)" if v["kind"] == "species" else
            "monoisotopic differences")
        tol = "%g %s%s" % (v["tol"][0], v["tol"][1], " (automatic)" if v["auto"] else "")
        ref = ("reference %s Da%s" % (self._mass_text(v, r["ref"]), self._ref_how(None, v)) if r["ref"] >= 0
               else "no reference")
        txt = "%s, %s, tolerance %s" % (ref, what, tol)
        low = sum(1 for q in r["species"] if q["status"] == 5)
        if low:
            txt += "; %d mass%s below %g %% of the tallest left out" % (low, "" if low == 1 else "es",
                                                                      100.0 * v["min_rel"])
        return txt

    def _ref_how(self, ent, v):
        """How the reference was chosen: " (tallest)", " (the lightest of the tag series)" (automatic) or
        " (chosen)"."""
        if not v.get("auto_ref", True):
            return " (chosen)"
        return " (tallest)" if v["r"]["ref"] == _tallest(v["r"]) else " (the lightest of the tag series)"

    # ----------------------------------------------------------- text of a pair
    def comp_text(self, v, comp, empty="reference"):
        return MS.comp_text(comp, v["names"], v["tags"], empty)

    def match_text(self, v, pr, alt=True):
        if pr["kind"] == 2:
            return "unknown"
        nt, ns = len(v["tags"]), len(v["names"])
        txt = self.comp_text(v, MS.pair_comp(pr, nt, ns))
        if alt and pr.get("alt"):
            txt += " (or %s)" % self.comp_text(v, MS.pair_comp(pr["alt"], nt, ns))
        return txt

    def label_text(self, v, pr):
        d = ("%%+.%df" % v["dd"]) % pr["delta"]
        if pr["kind"] == 2:
            return "%s unknown" % d
        return "%s %s" % (self.match_text(v, pr, alt=False), d)

    def error_text(self, v, pr):
        if pr["kind"] == 2 or not math.isfinite(pr["error"]):
            return ""
        if v["basis"] == "average":
            return "%+.2f Da" % pr["error"]
        heavier = max(v["masses"][pr["from"]], v["masses"][pr["to"]])
        return "%+.1f ppm" % (1e6 * pr["error"] / heavier)

    def shown_pairs(self, v):
        """The pairs drawn: the links (each mass to the reference or to its neighbouring species, unknown
        masses from the reference), or every explained pair and the unknown links."""
        r = v["r"]
        if v["show_all"]:
            return [p for p in r["pairs"] if p["kind"] in (0, 1) or (p["kind"] == 2 and p["link"])]
        return [p for p in r["pairs"] if p["link"] and p["kind"] in (0, 1, 2)]

    def table_pairs(self, v):
        """The rows of the table: every pair with the reference, the links, and with Show every matched
        pair every explained pair; by the mass of To, then From."""
        r = v["r"]
        rows = [p for p in r["pairs"] if p["kind"] in (0, 1, 2) and (p["ref"] or p["link"] or v["show_all"])]
        return sorted(rows, key=lambda p: (v["masses"][p["to"]], v["masses"][p["from"]]))

    # ----------------------------------------------------------- tables
    @staticmethod
    def _cols(basis):
        return [("From (Da)", 68 if basis == "average" else 84), ("To (Da)", 68 if basis == "average" else 84),
                ("Δ (Da)", 60 if basis == "average" else 72), ("Match", 140),
                ("Error", 62 if basis == "average" else 66), ("To is", 130)]

    @staticmethod
    def _conj_cols():
        return [("Tags", 124), ("Height %", 66), ("Area %", 62), ("Masses (Da)", 160)]

    def table_rows(self, v):
        rows = []
        r = v["r"]
        for p in self.table_pairs(v):
            q = r["species"][p["to"]]
            what = "unknown" if q["status"] == 3 else "reference"
            if q["status"] in (1, 2):
                txt = self.comp_text(v, q["comp"])
                what = "reference \u2212 " + txt[1:] if txt.startswith("\u2212") else "reference + " + txt
            rows.append((self._mass_text(v, p["from"]), self._mass_text(v, p["to"]),
                         ("%%+.%df" % v["dd"]) % p["delta"], self.match_text(v, p), self.error_text(v, p), what))
        return rows

    def conj_rows(self, v):
        """Rows of the degree of conjugation: per tag one row per k (up to the largest k found), then
        the average number of tags per molecule."""
        r = v["r"]
        out = []
        ar = v.get("real_areas", True)  # no areas from the method (UniDec envelopes, IsoDec): left empty
        for ti, c in enumerate(r["conj"]):
            name = v["tags"][ti]
            classes = {}
            for j, q in enumerate(r["species"]):
                if q["status"] not in (0, 1, 2):
                    continue
                tc = q["comp"][0]
                kk = tc[ti]
                if kk < 0 or kk > v["kmax"] or any(tc[x] for x in range(len(tc)) if x != ti):
                    continue
                classes.setdefault(kk, []).append(v["masses"][j])
            top = max([k for k in classes] + [1])
            for k in range(0, top + 1):
                ms_ = ", ".join(_fmt_num(m, v["dm"]) for m in sorted(classes.get(k, [])))
                out.append(("%d %s" % (k, name), _fmt_num(c["height"][k], 1),
                            _fmt_num(c["area"][k], 1) if ar else "", ms_))
            out.append(("Tags per molecule", _fmt_num(c["avg_height"], 2), _fmt_num(c["avg_area"], 2) if ar else "",
                        "%d mass%s counted" % (c["n"], "" if c["n"] == 1 else "es")))
        return out

    def conj_info(self, v):
        r = v["r"]
        unk = sum(1 for q in r["species"] if q["status"] == 3)
        txt = "by height and area" if v.get("real_areas", True) else "by height (no peak areas)"
        if unk:
            txt += "; %d unknown mass%s left out" % (unk, "" if unk == 1 else "es")
        neg = sum(1 for q in r["species"] if q["status"] in (1, 2) and any(c < 0 for c in q["comp"][0]))
        if neg:
            txt += "; %d mass%s lighter than the reference by a tag" % (neg, "" if neg == 1 else "es")
        return txt

    def fill(self, ent):
        """The tables of ent from the finder's result."""
        cs, cc = getattr(ent, "card_shift", None), getattr(ent, "card_conj", None)
        if cs is None or not self.is_on(ent):
            return
        v = self.view(ent)
        try:
            if v is None:
                cs._shift_rows = []
                cs.set_rows([], "No masses in this result")
                cc.set_rows([], "")
                return
            cs.set_columns(self._cols(v["basis"]))
            rows = self.table_rows(v)
            cs._shift_rows = self.table_pairs(v)
            n_unk = sum(1 for p in cs._shift_rows if p["kind"] == 2)
            cs.set_rows(rows, self._info(v) + ("; %d unknown" % n_unk if n_unk else ""))
            if v["tags"]:
                cc.set_rows(self.conj_rows(v), self.conj_info(v))
        except RuntimeError:
            pass

    def _row_selected(self, ent, i):
        rows = getattr(ent.card_shift, "_shift_rows", [])
        if 0 <= i < len(rows):
            j = rows[i]["to"]
            if ent is not self.dec.active:
                self.dec.activate(ent)
            if self.dec.sel != j:
                self.dec.select_peak(j)

    def _table_menu(self, ent, card):
        items = []
        rows = getattr(ent.card_shift, "_shift_rows", [])
        sel = ent.card_shift.selected() if card is ent.card_shift else []
        v = self.view(ent)
        if sel and v is not None and sel[0] < len(rows):
            j = rows[sel[0]]["to"]
            items.append(("Use %s Da as the unmodified species" % self._mass_text(v, j),
                          lambda: self.set_ref(ent, j)))
        if self.state(ent)["ref"] is not None:
            items.append(("Automatic reference", lambda: self.set_ref(ent, None)))
        items += [("Mass shift settings…", lambda: self.open_dialog(ent)),
                  ("Hide the mass shifts", lambda: self.set_on(ent, False)), (None, None)]
        return items

    # ----------------------------------------------------------- menus of the plot
    def menu_items(self, ent, x):
        """Entries of the right click menu of the zero charge mass spectrum."""
        on = self.is_on(ent)
        items = [("Mass shifts", lambda: self.set_on(ent, not on), on)]
        if on:
            i = self.dec._peak_near(x) if x is not None else None
            v = self.view(ent)
            if i is not None and v is not None and i != v["r"]["ref"]:
                items.append(("Use %s Da as the unmodified species" % self._mass_text(v, i),
                              lambda: self.set_ref(ent, i)))
            if self.state(ent)["ref"] is not None:
                items.append(("Automatic reference", lambda: self.set_ref(ent, None)))
            s = settings()
            items.append(("Show every matched pair", lambda: self.set_show_all(not s["show_all"]), s["show_all"]))
        items += [("Mass shift settings…", lambda: self.open_dialog(ent)), (None, None)]
        return items

    def set_show_all(self, flag):
        s = settings()
        s["show_all"] = bool(flag)
        save_settings(s)
        self.redraw_all()

    # ----------------------------------------------------------- drawing
    def refresh(self, ent):
        """ent drawn again; a result without masses is not drawn: its table says so."""
        self.dec.plot_results(ent)
        if self.is_on(ent) and self.view(ent) is None:
            self.fill(ent)

    def request_redraw_all(self):
        """redraw_all once, after the present event (several calls in one event: one redraw)."""
        if _REDRAW["pending"]:
            return
        _REDRAW["pending"] = True

        def go():
            _REDRAW["pending"] = False
            self.redraw_all()
        try:
            wx.CallAfter(go)
        except Exception:  # no wx application (tests)
            go()

    def redraw_all(self):
        """Every result of every open file with the finder on drawn again (the settings changed)."""
        import deconv_tab as D
        for pn in list(D.DeconvPanel.instances):
            try:
                for ent in pn.results:
                    if ent.row is not None:
                        pn._apply_z(ent)
                        if pn.shifts.is_on(ent) or getattr(ent, "_shift_drawn", False):
                            pn.shifts.refresh(ent)
            except RuntimeError:
                pass

    def add_room(self, ent):
        """Room above the peaks for the brackets (points), added to the label room of the plot before its y
        range is made (plot_results): the room the last drawing needed (draw() adjusts it)."""
        cm = ent.card_mass if ent is not None else None
        if cm is None:
            return 0.0
        if not self.is_on(ent):
            cm._shift_room = 0.0
        return float(getattr(cm, "_shift_room", 0.0))

    def _brackets(self, ent, v):
        """The brackets of the view: [{"pair", "xa", "xb" (data), "f0", "f1" (extent with the label, axes
        fraction), "fc" (label centre), "text"}]."""
        import deconv_tab as D
        cm = ent.card_mass
        x0v, x1v = cm.ax.get_xlim()
        wv = (x1v - x0v) or 1.0
        sw, sh = D._label_room(cm.ax, screen=True)
        res = ent.res
        peaks = res.get("peaks") or []
        mass = res.get("mass")
        iso = res.get("tag") == "isodec"
        resolved = D._grouped(res)
        out = []
        unknown = []
        for pr in self.shown_pairs(v):
            a, b = pr["from"], pr["to"]
            if not (0 <= a < len(peaks) and 0 <= b < len(peaks)):
                continue
            xb = D.DeconvPanel._marker_xy(peaks[b], mass, iso, resolved)[0]
            # an unknown mass is marked above its own peak (a bracket from the reference crossed the plot)
            xa = xb if pr["kind"] == 2 else D.DeconvPanel._marker_xy(peaks[a], mass, iso, resolved)[0]
            lo, hi = min(xa, xb), max(xa, xb)
            if hi < x0v or lo > x1v:
                continue
            text = self.label_text(v, pr)
            fa, fb = (max(lo, x0v) - x0v) / wv, (min(hi, x1v) - x0v) / wv
            fc = 0.5 * (fa + fb)
            hw = 0.5 * len(text) * 0.55 * 7.0 / sw
            bk = {"pair": pr, "xa": xa, "xb": xb, "f0": min(fa, fc - hw) - 2.0 / sw,
                  "f1": max(fb, fc + hw) + 2.0 / sw, "fc": fc, "text": text,
                  "h": float(peaks[b].get("height", 0.0) or 0.0)}
            (unknown if pr["kind"] == 2 else out).append(bk)
        unknown = sorted(unknown, key=lambda b: -b["h"])[:MAX_UNKNOWN]
        # the most important first when there are many: links of the reference, then the tallest masses
        out = sorted(out, key=lambda b: (0 if b["pair"]["ref"] else 1, -b["h"]))[:MAX_BRACKETS]
        return out + unknown, (sw, sh)

    def _obstacles(self, ent, sw, h_pt):
        """What the brackets must stay above, in points from the bottom of the plot: the drawn mass spectrum
        (a function of the axes fraction) and the boxes of the mass labels [(f0, f1, top)]."""
        from matplotlib.text import Annotation
        cm = ent.card_mass
        ax = cm.ax
        x0v, x1v = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        wv, hy = (x1v - x0v) or 1.0, (y1 - y0) or 1.0
        res = ent.res
        if res.get("tag") == "isodec":
            pk = res.get("peaks") or []
            xs = np.array(sorted(float(q["mass"]) for q in pk)) if pk else np.zeros(0)
            ys = np.array([float(q["height"]) for q in sorted(pk, key=lambda q: float(q["mass"]))]) if pk else np.zeros(0)
        else:
            mass = res.get("mass")
            xs = np.asarray(mass[:, 0], float) if mass is not None and len(mass) else np.zeros(0)
            ys = np.asarray(mass[:, 1], float) if mass is not None and len(mass) else np.zeros(0)

        def curve(f0, f1):
            if not len(xs):
                return 0.0
            a, b = np.searchsorted(xs, [x0v + f0 * wv, x0v + f1 * wv])
            a, b = max(0, a - 1), min(len(xs), b + 1)
            if b <= a:
                return 0.0
            return (float(np.nanmax(ys[a:b])) - y0) / hy * h_pt
        boxes = []
        for t in ax.texts:
            if not isinstance(t, Annotation) or getattr(t, "_shift_art", False) or not getattr(t, "_dec_label", False):
                continue
            try:
                fx = (t.xy[0] - x0v) / wv
                dy = float(t.xyann[1]) if t.get_anncoords() == "offset points" else 5.0
                fs = float(t.get_fontsize())
                yb = (t.xy[1] - y0) / hy * h_pt + dy
                hw = 0.5 * len(t.get_text()) * 0.6 * fs / sw
                ha = t.get_ha()
                f0, f1 = (fx - 2 * hw, fx) if ha == "right" else (fx, fx + 2 * hw) if ha == "left" else (fx - hw, fx + hw)
                boxes.append((f0, f1, yb + 1.25 * fs))
            except Exception:
                continue
        return curve, boxes

    def _place(self, ent, brackets, sw, h_pt):
        """Each bracket just above what lies below it (spectrum, mass labels, brackets placed before it),
        the shortest first: sets b["y"] (points, the bracket line) and returns the highest top needed. When
        a bracket of tags finds no room, the tags are placed first and a shorter one is left out instead."""
        need = self._place_order(ent, brackets, sw, h_pt, False)
        if any(b["pair"]["kind"] == 1 and b["y"] + ROW_PT - 3.0 > h_pt + 1.0 for b in brackets):
            need = self._place_order(ent, brackets, sw, h_pt, True)
        return need

    def _place_order(self, ent, brackets, sw, h_pt, tags_first):
        curve, boxes = self._obstacles(ent, sw, h_pt)
        placed = []  # (f0, f1, bottom, top)
        need = 0.0
        for b in sorted(brackets, key=lambda b: ((0 if b["pair"]["kind"] == 1 else 1) if tags_first else 0,
                                                 b["f1"] - b["f0"], b["fc"])):
            base = curve(b["f0"], b["f1"]) + 3.0
            for f0, f1, top in boxes:
                if f1 > b["f0"] and f0 < b["f1"]:
                    base = max(base, top + 1.0)
            y = base + 3.0  # the line (its ticks reach 3 pt down)
            moved = True
            while moved:
                moved = False
                for f0, f1, bot, top in placed:
                    if f1 > b["f0"] and f0 < b["f1"] and y - 3.0 < top and bot < y + ROW_PT - 3.0:
                        y = top + 3.5
                        moved = True
            b["y"] = y
            placed.append((b["f0"], b["f1"], y - 3.0, y + ROW_PT - 3.0))
            need = max(need, y + ROW_PT - 3.0)
        return need

    def draw(self, ent):
        """Brackets of the finder on the zero charge mass spectrum of ent (called after its mass labels
        were placed: on a new plot, after a zoom and for copied, saved and report images)."""
        cm = ent.card_mass if ent is not None else None
        if cm is None:
            return
        ax = cm.ax
        for a in list(ax.lines) + list(ax.texts):
            if getattr(a, "_shift_art", False):
                try:
                    a.remove()
                except Exception:
                    pass
        ent._shift_drawn = False
        if not self.is_on(ent):
            return
        v = self.view(ent)
        if v is None:
            return
        brackets, (sw, sh) = self._brackets(ent, v)

        def h_now():
            try:
                return max(20.0, ax.get_window_extent().height / ax.figure.dpi * 72.0)
            except Exception:
                return sh
        h_pt = h_now()
        need = self._place(ent, brackets, sw, h_pt) if brackets else 0.0
        # more room above the peaks when the brackets need it (the y range is made again), less when they
        # need less; at most what the plot gives (autoscale keeps 40 % of the height at most)
        for _ in range(3):
            if cm.ylock is not None:
                break
            over = need - (h_pt - 2.0)
            have = float(getattr(cm, "_shift_room", 0.0))
            base = cm.label_room_pt - have
            if over > 0.5 or (over < -ROW_PT and have > 0):
                # no more than the plot gives (autoscale_y keeps 40 % of the height for labels at most)
                new_room = min(max(0.0, have + over + (1.0 if over > 0 else 0.0)), max(0.0, 0.4 * h_pt - base))
                if abs(new_room - have) < 0.25:
                    break
                cm.label_room_pt = base + new_room
                cm._shift_room = new_room
                pad = cm._auto[0] if getattr(cm, "_auto", None) else 0.16
                cm.autoscale_y(pad=pad)
                h_pt = h_now()
                need = self._place(ent, brackets, sw, h_pt) if brackets else 0.0
            else:
                break
        y0, y1 = ax.get_ylim()
        hy = (y1 - y0) or 1.0

        def yd(pt):
            return y0 + pt / h_pt * hy
        import deconv_tab as D
        res = ent.res
        peaks = res.get("peaks") or []
        mass = res.get("mass")
        iso = res.get("tag") == "isodec"
        resolved = D._grouped(res)
        left = 0
        for b in brackets:
            if b["y"] + ROW_PT - 3.0 > h_pt + 1.0:
                left += 1  # no room in the plot: in the table only
                continue
            pr = b["pair"]
            kind = pr["kind"]
            col = COL_UNKNOWN if kind == 2 else COL_TAG if kind == 1 else COL_SHIFT
            y_line, y_tick = yd(b["y"]), yd(b["y"] - 3.0)
            xa, xb = b["xa"], b["xb"]
            if kind == 2:  # unknown: a short stroke under its label, above its peak
                ln = ax.plot([xb, xb], [y_tick, y_line], color=col, lw=0.9, zorder=7, solid_capstyle="butt")[0]
            else:
                ln = ax.plot([xa, xa, xb, xb], [y_tick, y_line, y_line, y_tick], color=col, lw=0.9, zorder=7,
                             solid_capstyle="butt")[0]
            ln._shift_art = ln._no_scale = True
            for i, x in (((pr["to"], xb),) if kind == 2 else ((pr["from"], xa), (pr["to"], xb))):
                if not (0 <= i < len(peaks)):
                    continue
                ty = D.DeconvPanel._marker_xy(peaks[i], mass, iso, resolved)[1]
                if (b["y"] - 3.0) - (ty - y0) / hy * h_pt > 6.0:
                    g = ax.plot([x, x], [ty + 2.5 / h_pt * hy, y_tick], color=COL_GUIDE, lw=0.6, ls=(0, (1, 2)),
                                zorder=2)[0]
                    g._shift_art = g._no_scale = True
            xc = ax.get_xlim()[0] + b["fc"] * (ax.get_xlim()[1] - ax.get_xlim()[0])
            t = ax.annotate(b["text"], (xc, y_line), xytext=(0, 1.0), textcoords="offset points", ha="center",
                            va="bottom", fontsize=7, color=col, zorder=8, annotation_clip=False,
                            bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none", alpha=0.8))
            t._shift_art = t._dec_label = True
        if left:
            t = ax.text(0.995, 0.995, "%d more in the table" % left, transform=ax.transAxes, ha="right", va="top",
                        fontsize=6.5, color=U.MUTED, zorder=8)
            t._shift_art = True
        ent._shift_drawn = True
        if getattr(ent, "_shift_filled", None) is not v:
            ent._shift_filled = v
            self.fill(ent)

    # ----------------------------------------------------------- report
    def report_tables(self, ent):
        """For a report: {"title", "info", "pairs": (header, rows, widths), "conj": (header, rows, widths) or
        None, "conj_info"}, or None while the finder is off for ent."""
        if ent is None or not self.is_on(ent):
            return None
        v = self.view(ent)
        if v is None:
            return None
        hdr = [c[0] for c in self._cols(v["basis"])]
        rows = [list(r) for r in self.table_rows(v)]
        env = v["basis"] == "average"
        wid = [2.0, 2.0, 1.6, 4.4, 1.7, 4.0] if env else [2.3, 2.3, 1.8, 4.2, 1.7, 3.6]
        out = {"title": "Mass shifts", "info": self._info(v), "pairs": (hdr, rows, wid), "conj": None,
               "conj_info": ""}
        if v["tags"]:
            out["conj"] = ([c[0] for c in self._conj_cols()], [list(r) for r in self.conj_rows(v)],
                           [3.4, 2.0, 2.0, 7.0])
            out["conj_info"] = self.conj_info(v)
        return out

    # ----------------------------------------------------------- settings window
    def open_dialog(self, ent=None):
        ent = ent if ent is not None else self.dec.active
        dlg = ShiftDialog(self.dec.tab, self, ent)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return
            s, ref = dlg.result()
        finally:
            dlg.Destroy()
        save_settings(s)
        if ent is not None and ent in self.dec.results:
            st = self.state(ent)
            st["ref"] = ref
            if not st["on"]:
                st["on"] = True
        self.redraw_all()
        if ent is not None and ent in self.dec.results:
            v = self.view(ent)
            if v is not None:
                self.dec.status("Mass shifts: " + self._info(v))


class ShiftDialog(wx.Dialog):
    """The Mass shifts window: tolerance, reference, tags, the shift list and the display."""

    def __init__(self, parent, panel, ent):
        wx.Dialog.__init__(self, wx.GetTopLevelParent(parent), title="Mass shifts",
                           style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        self.panel, self.ent = panel, ent
        s = settings()
        self.shifts = [dict(e) for e in s["shifts"]]
        fm = self.form = U.Form(self, label_w=130, note_w=470)
        v = panel.view(ent) if ent is not None else None
        fm.section("Comparison", T.GROUP["yellow"])
        if v is not None:
            what = ("average masses" if v["basis"] == "average" else
                    "monoisotopic differences (estimated)" if v["kind"] == "species" else "monoisotopic differences")
            fm.note("This result: %s." % what)
        self.tol = {}
        for key, lab, default in (("tol_envelope", "Tolerance, envelopes", "1.5 Da, at least 1.5 mass steps"),
                                  ("tol_resolved", "Tolerance, resolved", "10 ppm")):
            t = s[key]
            txt = fm.text("" if t is None else "%g" % t[0], 64)
            unit = fm.choice(["Da", "ppm"], 64)
            unit.SetSelection(1 if (t and t[1] == "ppm") or (t is None and key == "tol_resolved") else 0)
            fm.row(lab, txt, unit, "empty: %s" % default)
            self.tol[key] = (txt, unit)
        self.min_rel = fm.text("%g" % (100.0 * s["min_rel"]), 64)
        fm.row("Masses used", self.min_rel, "% of the tallest at least")
        refs = ["Automatic"]
        self.ref_masses = [None]
        if ent is not None:
            for q in ent.res.get("peaks") or []:
                refs.append("%s Da" % U._g(float(q["mass"])) if v is None else
                            "%s Da" % _fmt_num(float(q["mass"]), max(v["dm"], 1)))
                self.ref_masses.append(float(q["mass"]))
        self.ref = fm.choice(refs, 200)
        cur = panel.state(ent)["ref"] if ent is not None else None
        if cur is not None:
            for n, m in enumerate(self.ref_masses):
                if m is not None and abs(m - cur) <= 1e-6 * max(1.0, abs(cur)):
                    self.ref.SetSelection(n)
        fm.row("Unmodified species", self.ref)
        fm.section("Tags (degree of conjugation)", T.GROUP["green"])
        self.tags = []
        tg = s["tags"] + [{"name": "", "avg": None, "mono": None}] * (MS.MAX_TAGS - len(s["tags"]))
        for n, t in enumerate(tg[:MS.MAX_TAGS]):
            name = fm.text(t["name"], 110)
            avg = fm.text("" if t["avg"] is None else "%.4f" % t["avg"], 84)
            mono = fm.text("" if t["mono"] is None or t["mono"] == t["avg"] else "%.4f" % t["mono"], 84,
                           tooltip="Empty: the average mass")
            fm.row("Tag %d" % (n + 1), name, "average", avg, "monoisotopic", mono, "Da")
            self.tags.append((name, avg, mono))
        self.kmax = fm.text(str(s["kmax"]), 48, tooltip="1 to 20")
        self.extra = fm.choice(["0", "1", "2"], 48, tooltip="Other shifts in one difference with the tags")
        self.extra.SetSelection(int(s["max_extra"]))
        fm.row("At most", self.kmax, "tags per molecule, with up to", self.extra, "other shifts")
        fm.section("Shifts", T.GROUP["red"])
        self.list = wx.ListCtrl(fm, style=wx.LC_REPORT | wx.BORDER_THEME, size=wx.Size(-1, self.FromDIP(150)))
        self.list.SetFont(ui_font(9, 400))
        try:
            self.list.EnableCheckBoxes(True)
        except Exception:
            pass
        for n, (lab, w) in enumerate((("Name (ticked: in use)", 210), ("Δ average (Da)", 120),
                                      ("Δ monoisotopic (Da)", 140))):
            self.list.InsertColumn(n, lab, wx.LIST_FORMAT_LEFT if n == 0 else wx.LIST_FORMAT_RIGHT,
                                   width=self.FromDIP(w))
        fm.full(self.list)
        self.list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, lambda e: self.edit(e.GetIndex()))
        fm.buttons(U.flat(fm, "Add…", handler=lambda e: self.edit(None)),
                   U.flat(fm, "Edit…", handler=lambda e: self.edit(self.list.GetFirstSelected())),
                   U.flat(fm, "Remove", handler=lambda e: self.remove()),
                   U.flat(fm, "Defaults", handler=lambda e: self.defaults()))
        self.show_all = fm.check("Show every matched pair", s["show_all"])
        self._fill_list()
        vs = wx.BoxSizer(wx.VERTICAL)
        vs.Add(fm, 1, wx.EXPAND)
        vs.Add(U.dialog_buttons(self, U.flat(self, "Cancel", handler=lambda e: self.EndModal(wx.ID_CANCEL)),
                                U.flat(self, "Show", "primary", handler=lambda e: self._ok())), 0, wx.EXPAND)
        self.SetSizerAndFit(vs)
        try:  # never taller than the screen (the buttons were cut off): the list gets shorter
            room = T.display_area(self).GetHeight() - self.FromDIP(40)
            if self.GetSize()[1] > room:
                lh = max(self.FromDIP(80), self.list.GetSize()[1] - (self.GetSize()[1] - room))
                self.list.SetMinSize(wx.Size(-1, lh))
                self.list.SetSize(wx.Size(-1, lh))
                self.Fit()
        except Exception:
            pass
        self.CentreOnParent()

    def _fill_list(self):
        self.list.DeleteAllItems()
        for n, e in enumerate(self.shifts):
            i = self.list.InsertItem(n, e["name"])
            self.list.SetItem(i, 1, "%.5f" % e["avg"])
            self.list.SetItem(i, 2, "%.6f" % e["mono"])
            try:
                self.list.CheckItem(i, bool(e["on"]))
            except Exception:
                pass

    def _sync_checks(self):
        for n, e in enumerate(self.shifts):
            try:
                e["on"] = bool(self.list.IsItemChecked(n))
            except Exception:
                pass

    def edit(self, i):
        self._sync_checks()
        new = i is None or i < 0 or i >= len(self.shifts)
        e = {"name": "", "avg": "", "mono": "", "on": True} if new else self.shifts[i]
        vals = U.ask_form(self, "Add a shift" if new else "Edit the shift", [
            {"key": "name", "label": "Name", "value": e["name"], "width": 180},
            {"key": "avg", "label": "Δ average mass", "value": "" if e["avg"] == "" else "%.5f" % e["avg"],
             "unit": "Da", "width": 100, "tip": "A loss is negative"},
            {"key": "mono", "label": "Δ monoisotopic mass", "value": "" if e["mono"] == "" else "%.6f" % e["mono"],
             "unit": "Da", "width": 100, "tip": "Empty: the average mass"}], note=None, ok="OK")
        if vals is None:
            return
        avg, mono = MS._num(vals["avg"]), MS._num(vals["mono"])
        avg = mono if avg is None else avg
        mono = avg if mono is None else mono
        if avg is None or avg == 0 or mono == 0:
            wx.MessageBox("The shift needs a mass difference other than 0.", "Mass shifts", wx.ICON_INFORMATION)
            return
        entry = {"name": vals["name"].strip() or "shift", "avg": avg, "mono": mono, "on": True if new else e["on"]}
        if new:
            if len(self.shifts) >= MS.MAX_SHIFTS:
                wx.MessageBox("At most %d shifts." % MS.MAX_SHIFTS, "Mass shifts", wx.ICON_INFORMATION)
                return
            self.shifts.append(entry)
        else:
            self.shifts[i] = entry
        self._fill_list()

    def remove(self):
        self._sync_checks()
        sel = []
        i = self.list.GetFirstSelected()
        while i >= 0:
            sel.append(i)
            i = self.list.GetNextSelected(i)
        for i in reversed(sel):
            del self.shifts[i]
        self._fill_list()

    def defaults(self):
        self.shifts = MS.default_shifts()
        self._fill_list()

    def result(self):
        """(settings, reference mass or None)."""
        self._sync_checks()
        s = settings()
        s["shifts"] = [dict(e) for e in self.shifts]
        tags = []
        for name, avg, mono in self.tags:
            a, m = MS._num(avg.GetValue()), MS._num(mono.GetValue())
            if a is None and m is None:
                continue
            a = m if a is None else a
            m = a if m is None else m
            tags.append({"name": name.GetValue().strip() or "tag", "avg": a, "mono": m})
        s["tags"] = tags
        k = MS._num(self.kmax.GetValue())
        s["kmax"] = int(min(MS.MAX_K, max(1, int(k)))) if k is not None else s["kmax"]
        s["max_extra"] = self.extra.GetSelection()
        s["show_all"] = self.show_all.GetValue()
        mr = MS._num(self.min_rel.GetValue())
        if mr is not None:
            s["min_rel"] = min(1.0, max(0.0, mr / 100.0))
        for key, (txt, unit) in self.tol.items():
            val = MS._num(txt.GetValue())
            s[key] = [val, unit.GetStringSelection()] if val is not None and val > 0 else None
        ref = self.ref_masses[self.ref.GetSelection()] if self.ref.GetSelection() >= 0 else None
        return MS.clean_settings(s), ref

    def _ok(self):
        for name, avg, mono in self.tags:
            for c in (avg, mono):
                t = c.GetValue().strip()
                if t and (MS._num(t) is None or MS._num(t) <= 0):
                    wx.MessageBox("A tag mass must be a number above 0.", "Mass shifts", wx.ICON_INFORMATION)
                    c.SetFocus()
                    return
        self.EndModal(wx.ID_OK)
