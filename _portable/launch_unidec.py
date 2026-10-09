"""
Launcher of MSpektra, a portable build of UniDec 8.2.1 (unmodified code
from PyPI) with add-ons for LC-MS and HRMS analysis.

UniDec: Marty et al., Anal. Chem. 2015, DOI 10.1021/acs.analchem.5b00140.
Please cite this paper in any publication that uses UniDec, including the
deconvolution in LCMS Analysis and HRMS Analysis.

What this launcher adds on top of a normal "pip install unidec":
  * runs from a self contained, isolated Python (no Anaconda / PATH / user
    site-packages interference),
  * removes the Windows "downloaded from the internet" block (Mark of the Web)
    from all files on first start (this block is what breaks Thermo .raw
    loading and some DLLs in the official zip),
  * keeps matplotlib and numba settings/caches inside this folder,
  * writes a log file (logs\\) when started without a console and records
    native crashes there, so problems can be diagnosed,
  * add-ons (all in this folder, UniDec's own code is unchanged):
      unidec_theme.py      start screen (LCMS Analysis, HRMS Analysis,
                           Deconvolute), modern interface, sharp on scaled
                           displays, window memory,
      unilcms.py           LCMS Analysis: Shimadzu LabSolutions .lcd files, MS
                           and PDA (lcms_data.py, lcms_pda.py,
                           lcms_integrate.py),
      hrms.py              HRMS Analysis: Bruker .d folders (Baf2Sql library)
                           and mzML; internal calibration (hrms_calib.py),
                           exact mass tools (ms_formula.py), reader
                           (hrms_data.py),
      deconv_tab.py        deconvolution inside both analysis windows
                           (ms_deconv.py: UniDec, maximum entropy, IsoDec),
      unidec_jcamp.py      opens JCAMP-DX (.jdx/.dx/.jcamp) mass spectra,
      unidec_ui_addons.py  right-click menu on plots, persistent and movable
                           labels, taskbar icon grouping,
      unidec_plotstyle.py  modern boxed plot style (Classic style switch in
                           the plot menu),
      unidec_fast.py       loads the Thermo/.NET reader only when needed,
  * a data file given as argument (drag and drop onto MSpektra.exe, or
    "Open with") opens directly: .lcd in LCMS Analysis, a Bruker .d folder or
    .mzML in HRMS Analysis, anything else in the Deconvolute window,
  * works from write-protected locations (e.g. C:\\Program Files, or a folder
    copied to C:\\ with administrator permission): UniDec normally writes its
    recent-file list and default settings into its own folder, which fails
    there and leaves the Deconvolute window unusable. In that case these
    files, the logs and the caches are kept in %LOCALAPPDATA%\\UniDecPortable.
"""
import os
import sys
import time
import shutil
import tempfile
import traceback
import faulthandler

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _app_version():
    """Version of MSpektra (APP_VERSION in unidec_theme.py), for the log."""
    try:
        import re
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "unidec_theme.py"),
                  encoding="utf-8") as fh:
            m = re.search(r'^APP_VERSION = "([^"]+)"', fh.read(), re.M)
        return m.group(1) if m else "?"
    except Exception:
        return "?"
PYDIR = os.path.join(ROOT, "python")
UNIDEC_BIN = os.path.join(PYDIR, "Lib", "site-packages", "unidec", "bin")


def _is_writable(path):
    test = os.path.join(path, ".write_test_%d" % os.getpid())
    try:
        with open(test, "w") as fh:
            fh.write("x")
        os.remove(test)
        return True
    except OSError:
        return False


READONLY = not _is_writable(ROOT)
if READONLY:
    _base = os.environ.get("LOCALAPPDATA") or tempfile.gettempdir()
    USERDIR = os.path.join(_base, "UniDecPortable")
else:
    USERDIR = ROOT
LOGDIR = os.path.join(USERDIR, "logs")
CONFDIR = os.path.join(USERDIR, "config")

_LOGFILE = None
_LOGPATH = None


def _setup_logging():
    """pythonw.exe has no console: redirect output to a log file."""
    global _LOGFILE, _LOGPATH
    if sys.stdout is not None and sys.stderr is not None:
        return
    for logdir in (LOGDIR, os.path.join(tempfile.gettempdir(), "UniDecPortable", "logs")):
        try:
            os.makedirs(logdir, exist_ok=True)
            old = sorted(f for f in os.listdir(logdir) if f.endswith(".log"))
            for name in old[:-20]:  # keep the 20 most recent logs
                try:
                    os.remove(os.path.join(logdir, name))
                except OSError:
                    pass
            # the process number in the name: the start screen and the windows it
            # starts as their own process (same second) must not share a file
            _LOGPATH = os.path.join(logdir, time.strftime("unidec_%Y%m%d_%H%M%S") + "_%d.log" % os.getpid())
            _LOGFILE = open(_LOGPATH, "w", encoding="utf-8", errors="replace", buffering=1)
            sys.stdout = _LOGFILE
            sys.stderr = _LOGFILE
            return
        except OSError:
            continue
    # Nothing writable at all: fall back to a sink so print() never fails
    _LOGFILE = open(os.devnull, "w")
    sys.stdout = sys.stderr = _LOGFILE


def _unblock_files():
    """Strip the Zone.Identifier stream (Mark of the Web) once."""
    if os.name != "nt":
        return
    marker = os.path.join(CONFDIR if READONLY else PYDIR, ".unblocked")
    if os.path.exists(marker):
        return
    count = 0
    # only files Windows checks for the block (DLLs and .NET assemblies,
    # programs, help): trying all 11 000 files took seconds on a new computer
    kinds = (".dll", ".exe", ".pyd", ".chm", ".bat", ".cmd", ".vbs", ".ps1", ".ocx", ".config", ".manifest")
    for dirpath, dirs, files in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d != "__pycache__" and not d.endswith(".dist-info")]
        for name in files:
            if not name.lower().endswith(kinds):
                continue
            try:
                os.remove(os.path.join(dirpath, name) + ":Zone.Identifier")
                count += 1
            except OSError:
                pass
    try:
        with open(marker, "w") as fh:
            fh.write(str(count))
    except OSError:
        pass
    print("First start: removed download block from %d files" % count)


def _apply_pending_updates():
    """Files of an update that could not be written while MSpektra was running (Windows locks a
    library in use) are delivered as <name>.new next to it: they replace the old file at the next
    start, before anything loads it. If another MSpektra window still uses it, the next start
    tries again (the program works with the older library meanwhile)."""
    folder = os.path.join(ROOT, "_portable", "msengine")
    try:
        names = [n for n in os.listdir(folder) if n.lower().endswith(".new")]
    except OSError:
        return
    for n in names:
        new, old = os.path.join(folder, n), os.path.join(folder, n[:-4])
        try:
            os.replace(new, old)
            print("Update: %s replaced by the new version" % n[:-4])
            try:
                os.remove(old + ":Zone.Identifier")  # a download block travels with the file
            except OSError:
                pass
        except OSError as ex:
            print("Update: %s still in use (%s); replaced at a later start" % (n[:-4], ex))


def _retire_old_launcher():
    """MSpektra.exe is the only program file of the folder: the old launcher UniDec.exe of the
    UniDec named versions (the same launcher under its old name, it starts nothing else) is moved
    to backups\\old_launcher once, so that it cannot be confused with the program."""
    if os.name != "nt" or READONLY:
        return
    old = os.path.join(ROOT, "UniDec.exe")
    if not (os.path.isfile(old) and os.path.isfile(os.path.join(ROOT, "MSpektra.exe"))):
        return
    dst_dir = os.path.join(ROOT, "backups", "old_launcher")
    try:
        os.makedirs(dst_dir, exist_ok=True)
        os.replace(old, os.path.join(dst_dir, "UniDec.exe"))
        print("Moved the old launcher UniDec.exe to backups\\old_launcher (MSpektra.exe starts the program)")
    except OSError as ex:  # e.g. in use: tried again at the next start
        print("Old launcher UniDec.exe not moved yet: %s" % ex)


def _isolate_sys_path():
    """Drop anything outside this folder (e.g. a per-user site-packages)."""
    root = os.path.normcase(os.path.abspath(ROOT))
    keep = []
    for p in sys.path:
        full = os.path.normcase(os.path.abspath(p)) if p else p
        if not p or full == root or full.startswith(root + os.sep):
            keep.append(p)
        else:
            print("Ignoring external path:", p)
    sys.path[:] = keep


def _setup_environment():
    for d in (os.path.join(CONFDIR, "matplotlib"), os.path.join(CONFDIR, "numba_cache")):
        try:
            os.makedirs(d, exist_ok=True)
        except OSError:
            pass
    # Keep settings isolated from any other Python/matplotlib on this PC
    os.environ["MPLCONFIGDIR"] = os.path.join(CONFDIR, "matplotlib")
    os.environ.pop("MPLBACKEND", None)
    os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(CONFDIR, "numba_cache"))
    # .NET Framework runtime for the Thermo RawFileReader (Windows default)
    os.environ.setdefault("PYTHONNET_RUNTIME", "netfx")
    # The package ships the compiled Python files (__pycache__, hash based,
    # so they stay valid after copying). No sys.pycache_prefix for a
    # write-protected folder: with a prefix Python looks only there and
    # ignores the shipped files (every module would be compiled again on the
    # first start). Where __pycache__ cannot be written, Python just skips
    # the write.


def _redirect_unidec_user_files():
    """UniDec writes recent.txt, recentCD.txt and default_conf.dat into its
    own bin folder; move them to a writable per-user folder."""
    udir = os.path.join(USERDIR, "unidec_user")
    os.makedirs(udir, exist_ok=True)
    recent = os.path.join(udir, "recent.txt")
    recent_cd = os.path.join(udir, "recentCD.txt")
    dconf = os.path.join(udir, "default_conf.dat")
    for p in (recent, recent_cd):
        if not os.path.exists(p):
            open(p, "a").close()
    if not os.path.exists(dconf):
        src = os.path.join(UNIDEC_BIN, "default_conf.dat")
        if os.path.exists(src):
            shutil.copyfile(src, dconf)

    from unidec.modules import unidecstructure

    original = unidecstructure.UniDecConfig.initialize_system_paths

    def initialize_system_paths(self, *args, **kwargs):
        result = original(self, *args, **kwargs)
        self.recentfile = recent
        self.recentfileCD = recent_cd
        self.defaultconfig = dconf
        return result

    unidecstructure.UniDecConfig.initialize_system_paths = initialize_system_paths
    print("UniDec settings and recent files:", udir)


def _fatal(msg):
    print(msg)
    try:
        import ctypes
        ctypes.windll.user32.MessageBoxW(None, msg, "MSpektra", 0x10)
    except Exception:
        pass


def _is_project(arg):
    """A saved analysis: a .msanalysis file or a folder holding session.msanalysis (as
    ms_project_store.is_project, without importing it)."""
    a = arg.rstrip("\\/")
    return a.lower().endswith(".msanalysis") or \
        (os.path.isdir(a) and os.path.isfile(os.path.join(a, "session.msanalysis")))


def _addon(label, func):
    try:
        func()
    except Exception:
        print("%s could not be enabled:\n%s" % (label, traceback.format_exc()))


def main():
    t0 = time.perf_counter()
    if READONLY:
        try:
            os.makedirs(USERDIR, exist_ok=True)
        except OSError:
            pass
    _setup_logging()
    try:
        faulthandler.enable(file=sys.stderr, all_threads=True)
    except Exception:
        pass
    print("MSpektra %s | Python %s | %s" % (_app_version(), sys.version.split()[0], ROOT))
    print("Started", time.strftime("%Y-%m-%d %H:%M:%S"))
    if READONLY:
        print("Program folder is write-protected; using", USERDIR,
              "for settings, logs and caches")
    _isolate_sys_path()
    _unblock_files()
    _apply_pending_updates()
    _retire_old_launcher()
    _setup_environment()
    sys.path.insert(0, os.path.join(ROOT, "_portable"))

    # Add-ons (each one is optional; a failure only disables that add-on).
    # Only wx is loaded before the first window. The libraries of the analysis
    # windows (NumPy, Matplotlib) load in the background while the start
    # screen is shown; the UniDec interface and its libraries (about 10 s on
    # a cold start) load only in the process of the Deconvolute window.
    def dpi():
        import unidec_theme
        scale = unidec_theme.enable_dpi_awareness()
        print("Display scaling: %d %%" % round(scale * 100))
    _addon("Display scaling", dpi)

    def fast():
        import unidec_fast
        if unidec_fast.install_lazy_thermo():
            print("Thermo reader will load on first use")
    _addon("Start-up speed add-on", fast)

    def app_id():
        import unidec_ui_addons
        unidec_ui_addons.set_process_app_id()
    _addon("Taskbar identity", app_id)

    settings_path = os.path.join(CONFDIR, "portable_settings.json")

    def light():
        import unidec_theme
        unidec_theme.configure(settings_path)
    _addon("Modern interface add-on", light)
    print("Ready for the first window after %.1f s" % (time.perf_counter() - t0))

    # A single data file as argument (drag and drop onto MSpektra.exe, or
    # "Open with"): open it straight in the right window.
    args = [a for a in sys.argv[1:] if a.strip()]
    # Saved analyses (.msanalysis files, analysis folders) open in the window of their kind, together with
    # the raw data files given with them. The project reader (and NumPy) loads only when there is one.
    saved = [a for a in args if _is_project(a)]
    projects = {"lcms": [], "hrms": []}
    if saved:
        args = [a for a in args if a not in saved]
        try:
            import ms_project_store
            reader = ms_project_store.peek
        except Exception as ex:
            def reader(path, ex=ex):
                raise ex
        for a in saved:
            try:
                kind = reader(a)["kind"]
            except Exception as ex:
                print("Saved analysis not opened: %s\n%s" % (a, traceback.format_exc()))
                _fatal("Could not open %s\n\n%s" % (a, ex))
                continue
            projects[kind].append(os.path.abspath(a))
        if not projects["lcms"] and not projects["hrms"] and not args:
            return 1
    # LCMS Analysis: "--lcms" or a Shimadzu .lcd file
    lcd = [a for a in args if a.lower().endswith(".lcd") and os.path.isfile(a)]
    # HRMS Analysis: "--hrms", a Bruker .d folder (or a file inside it) or .mzML
    hr = []
    try:
        if args:  # hrms_data imports NumPy: only when there is a file argument
            import hrms_data
            hr = [a for a in args if hrms_data.find_d_folder(a) or a.lower().endswith(".mzml")]
    except Exception:
        print("HRMS reader not available:\n" + traceback.format_exc())
    # every file dropped onto MSpektra.exe opens (one window per kind, a file each); projects first,
    # so a raw file of a project being opened is not asked about again
    lcms_paths = projects["lcms"] + [os.path.abspath(a) for a in lcd]
    hrms_paths = list(projects["hrms"])
    for a in hr:
        t = hrms_data.find_d_folder(a) or os.path.abspath(a)
        if t not in hrms_paths:
            hrms_paths.append(t)
    want_lcms = "--lcms" in args or bool(lcms_paths)
    want_hrms = "--hrms" in args or bool(hrms_paths)
    if want_lcms or want_hrms:
        sys.argv = [sys.argv[0]]
        _relaunch_hook()
        if want_lcms and want_hrms:
            import wx
            import unilcms
            import hrms
            app = wx.App(False)
            app.SetAppName(unilcms.APP)
            unilcms.open_window(lcms_paths or None)
            hrms.open_window(hrms_paths or None)
            app.MainLoop()
        elif want_lcms:
            import unilcms
            unilcms.run(lcms_paths or None)
        else:
            import hrms
            hrms.run(hrms_paths or None)
        return 0
    if len(args) == 1 and os.path.isfile(args[0]):
        path = os.path.abspath(args[0])
        sys.argv = [sys.argv[0]]
        if not ensure_unidec():
            return 1
        import wx
        from unidec.GUniDec import UniDecApp
        app = UniDecApp()
        wx.CallAfter(app.on_open_file, os.path.basename(path), os.path.dirname(path))
        app.start()
        return 0
    # The Deconvolute window ("--deconv") or a classic tool ("--classic N",
    # N as in unidec_theme.CLASSIC_TOOLS): started this way by the start
    # screen, so they run in their own process (see spawn below)
    if "--deconv" in args:
        sys.argv = [sys.argv[0]]
        if not ensure_unidec():
            return 3  # reported to the user already; the start screen does not retry
        from unidec.GUniDec import UniDecApp
        print("Launching UniDec")
        app = UniDecApp()
        app.start()
        return 0
    if "--classic" in args:
        try:
            n = int(args[args.index("--classic") + 1])
        except (IndexError, ValueError):
            n = 0
        sys.argv = [sys.argv[0]]
        if not ensure_unidec():
            return 3
        import importlib
        L = importlib.import_module("unidec.Launcher")
        print("Launching classic tool %d" % n)
        getattr(L.Lview, "button%d" % n)(None)
        return 0
    if not any(a.startswith("-") for a in args):
        try:
            import unidec_theme
            return unidec_theme.start_screen(ensure_unidec, ROOT, spawn=spawn)
        except Exception:
            print("Start screen failed, opening UniDec's launcher:\n" + traceback.format_exc())
    # UniDec's own launcher (its command line options, e.g. --meta, --chrom)
    if not ensure_unidec():
        return 1
    try:
        import unidec_theme
        unidec_theme._install_launcher()
    except Exception:
        pass
    from unidec.Launcher import run_launcher
    run_launcher()
    return 0


def _relaunch_hook():
    try:
        import unidec_ui_addons
        unidec_ui_addons.install_relaunch(ROOT)
    except Exception:
        pass


def spawn(args):
    """Starts MSpektra again in a new process with the given arguments
    (pythonw.exe, no console). The start screen opens the Deconvolute window
    and the classic tools this way, so UniDec's interface and its libraries
    (about 10 s on a cold start) load on another processor core and never
    into the process of the start screen and the analysis windows; LCMS and
    HRMS Analysis do the same for the deconvolution (deconv_worker.py)."""
    import subprocess
    exe = sys.executable
    if os.name == "nt" and exe.lower().endswith("python.exe"):
        w = exe[:-10] + "pythonw.exe"
        if os.path.isfile(w):
            exe = w
    cmd = [exe, "-s", os.path.abspath(__file__)] + [str(a) for a in args]
    return subprocess.Popen(cmd, cwd=ROOT, close_fds=True)


_UNIDEC = {"ready": None}


def ensure_unidec():
    """UniDec's interface and the add-ons for it; done once, when the first
    UniDec window (Deconvolute, classic tools) is opened."""
    if _UNIDEC["ready"] is not None:
        return _UNIDEC["ready"]
    t0 = time.perf_counter()
    settings_path = os.path.join(CONFDIR, "portable_settings.json")

    def jcamp():
        import unidec_jcamp
        unidec_jcamp.install()
    _addon("JCAMP-DX support", jcamp)
    try:
        if READONLY:
            _redirect_unidec_user_files()
        import unidec.GUniDec  # noqa: F401  (the Deconvolute window)
    except Exception:
        tb = traceback.format_exc()
        print(tb)
        where = ("\n\nLog: " + _LOGPATH) if _LOGPATH else ""
        _fatal("The Deconvolute window could not be loaded: " + tb.strip().splitlines()[-1] + where)
        _UNIDEC["ready"] = False
        return False

    def style():
        import unidec_plotstyle
        unidec_plotstyle.install(settings_path)
    _addon("Plot style add-on", style)

    def ui():
        import unidec_ui_addons
        unidec_ui_addons.install(ROOT)
    _addon("User-interface add-ons", ui)

    def theme():
        import unidec_theme
        unidec_theme.install(settings_path)
    _addon("Modern interface add-on", theme)
    print("UniDec interface loaded in %.1f s" % (time.perf_counter() - t0))
    _UNIDEC["ready"] = True
    return True


if __name__ == "__main__":
    sys.exit(main())
