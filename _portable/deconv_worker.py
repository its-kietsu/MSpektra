"""
Worker process for the deconvolutions of MSpektra (LCMS Analysis and HRMS
Analysis).

The deconvolution runs here, in its own Python process, and not in the
window's process: the window stays responsive whatever the settings, a run
can be cancelled at once (the process is ended), and the numerical work does
not compete with the interface for Python's global lock. The process is
started at the first deconvolution and kept for the next ones, so the
libraries are loaded only once.

Protocol (binary, over the standard input and output of this process):
8 byte little endian length, then a pickled Python object.
  parent -> worker   ("run", job) | ("quit", None)
  worker -> parent   ("ready", pid) | ("progress", i, n, text)
                     | ("done", result) | ("error", message, traceback)
Everything else the libraries print goes to the standard error stream, which
the window sends to logs\\deconvolution_worker.log. A line there marks the
start of the process and each job; a crash of native code (the engines)
leaves the Python stack of every thread there (faulthandler).
"""
import os
import sys
import time
import pickle
import struct
import threading
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))


def _read_exact(stream, n):
    buf = b""
    while len(buf) < n:
        chunk = stream.read(n - len(buf))
        if not chunk:
            return None
        buf += chunk
    return buf


def _log(text):
    """A line of the worker log (standard error), with the time and process."""
    try:
        sys.stderr.write("[deconvolution worker %d %s] %s\n" % (os.getpid(), time.strftime("%Y-%m-%d %H:%M:%S"),
                                                                text))
        sys.stderr.flush()
    except Exception:
        pass


def _version():
    try:
        import re
        with open(os.path.join(HERE, "unidec_theme.py"), encoding="utf-8") as fh:
            m = re.search(r'^APP_VERSION = "([^"]+)"', fh.read(), re.M)
        return m.group(1) if m else "?"
    except Exception:
        return "?"


def _job_text(job):
    """What a job is, for the log: method, number of points, m/z range."""
    try:
        spec = job.get("spec")
        n = len(spec) if spec is not None else 0
        txt = "%s on %d points" % (job.get("method"), n)
        if n:
            txt += ", m/z %.4g to %.4g" % (float(spec[0][0]), float(spec[-1][0]))
        if job.get("name"):
            txt += " (%s)" % job["name"]
        return txt
    except Exception:
        return str(job.get("method") if isinstance(job, dict) else job)


def main():
    # a crash in native code (UniDec engine, the C++ core) writes the Python
    # stack of every thread to the log before the process ends
    try:
        import faulthandler
        faulthandler.enable(file=sys.stderr, all_threads=True)
    except Exception:
        pass
    # protocol channel: a private copy of the standard output; the standard
    # output itself then points at the error stream, so prints of UniDec or
    # of its engine cannot corrupt the protocol
    out = os.fdopen(os.dup(1), "wb")
    os.dup2(2, 1)
    sys.stdout = sys.stderr
    inp = sys.stdin.buffer
    lock = threading.Lock()  # one message at a time, whatever thread reports progress

    def send(obj):
        data = pickle.dumps(obj, protocol=4)
        with lock:
            out.write(struct.pack("<Q", len(data)))
            out.write(data)
            out.flush()

    t_start = time.time()
    _log("started (MSpektra %s, Python %s)" % (_version(), sys.version.split()[0]))
    sys.path.insert(0, HERE)
    try:
        import launch_unidec as L
        L._isolate_sys_path()
        L._setup_environment()
        if getattr(L, "READONLY", False):
            L._redirect_unidec_user_files()
    except Exception:
        traceback.print_exc()
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    import numpy as np  # noqa: F401  (loaded once, for every job)
    import ms_deconv as D
    try:
        import msengine_py
        core = msengine_py.lib() is not None  # UniDec, IsoDec and maximum entropy in C++
    except Exception:
        core = False
    if not core:
        try:
            from unidec import engine  # noqa: F401  (slow import, done before the first job)
        except Exception:
            traceback.print_exc()
    _log("ready after %.1f s (%s)" % (time.time() - t_start, "C++ core" if core else "Python code"))
    send(("ready", os.getpid()))

    while True:
        hdr = _read_exact(inp, 8)
        if hdr is None:
            break
        n = struct.unpack("<Q", hdr)[0]
        body = _read_exact(inp, n)
        if body is None:
            break
        try:
            kind, job = pickle.loads(body)
        except Exception:
            # the stream is out of step: nothing that follows can be trusted
            _log("unreadable message from the window, stopping:\n" + traceback.format_exc())
            break
        if kind == "quit":
            break
        if kind != "run":
            continue
        t0 = time.time()
        _log("job: " + _job_text(job))

        def progress(i, n_, text, _t0=t0):
            send(("progress", int(i), int(n_), str(text)))

        try:
            res = D.run_method(job["method"], job["spec"], job.get("folder"), job.get("name"), job["params"],
                               progress=progress)
            res["seconds_worker"] = time.time() - t0
            send(("done", res))
            _log("done in %.1f s" % (time.time() - t0))
        except (Exception, SystemExit) as ex:  # a library calling sys.exit() must not end the worker silently
            tb = traceback.format_exc()
            _log("failed after %.1f s: %s\n%s" % (time.time() - t0, str(ex) or ex.__class__.__name__, tb))
            send(("error", str(ex) or ex.__class__.__name__, tb))
    _log("stopped")


if __name__ == "__main__":
    main()
