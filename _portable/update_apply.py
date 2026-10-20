"""Installs a downloaded update of MSpektra (updater.py) after MSpektra has closed: the changed parts in
_update\\staging take the place of the installed ones, which move to _update\\old; every move is written to
_update\\old\\journal.json, so that it can be undone. Then MSpektra starts again; when its first window does
not come up (it removes _update\\pending.json), the previous version is put back and started.

    update_apply.py --root <MSpektra folder> --pid <MSpektra process>             install the update
    update_apply.py --root <MSpektra folder> --pid <MSpektra process> --rollback  back to _update\\old

Runs from a copy in the temp folder with the Python of the MSpektra folder (which an update never replaces),
standard library only. Log: logs\\update.log in the MSpektra folder."""
import ctypes
import json
import os
import shutil
import subprocess
import sys
import time
from ctypes import wintypes as W

ARGS = sys.argv[1:]
ROOT = os.path.abspath(ARGS[ARGS.index("--root") + 1])
PID = int(ARGS[ARGS.index("--pid") + 1]) if "--pid" in ARGS else 0
UPD = os.path.join(ROOT, "_update")
STAGING = os.path.join(UPD, "staging")
OLD = os.path.join(UPD, "old")
PENDING = os.path.join(UPD, "pending.json")
K32 = ctypes.windll.kernel32
U32 = ctypes.windll.user32
K32.OpenProcess.restype = W.HANDLE
K32.OpenProcess.argtypes = [W.DWORD, W.BOOL, W.DWORD]
K32.CloseHandle.argtypes = [W.HANDLE]
K32.WaitForSingleObject.argtypes = [W.HANDLE, W.DWORD]
K32.TerminateProcess.argtypes = [W.HANDLE, W.UINT]
K32.QueryFullProcessImageNameW.argtypes = [W.HANDLE, W.DWORD, W.LPWSTR, ctypes.POINTER(W.DWORD)]
U32.GetWindowThreadProcessId.argtypes = [W.HWND, ctypes.POINTER(W.DWORD)]
U32.IsWindowVisible.argtypes = [W.HWND]
U32.GetWindowTextLengthW.argtypes = [W.HWND]
U32.GetWindowTextW.argtypes = [W.HWND, W.LPWSTR, ctypes.c_int]
LOCKED = (5, 32, 33)  # access denied, sharing violation, lock violation: a file held by another program


def log(*a):
    try:
        os.makedirs(os.path.join(ROOT, "logs"), exist_ok=True)
        with open(os.path.join(ROOT, "logs", "update.log"), "a", encoding="utf-8") as f:
            f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + " ".join(str(x) for x in a) + "\n")
    except OSError:
        pass


def message(text, retry=False):
    """A Windows message box (the window of MSpektra is closed); True for Retry."""
    flags = (0x5 if retry else 0x0) | 0x40 | 0x10000  # RETRYCANCEL / OK, information, foreground
    return U32.MessageBoxW(None, text, "MSpektra update", flags) == 4


def read_json(path):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def replace(src, dst, tries=240):
    """os.replace, tried again for up to a minute while another program holds a file for a moment: OneDrive
    syncing a folder on the OneDrive Desktop, or a virus scanner (15.10.2026: "Access is denied" on journal.json
    at the 102nd move of the update to 4.36, which then was undone; a sharing violation is an OSError too)."""
    for i in range(tries):
        try:
            os.replace(src, dst)
            return
        except OSError as ex:
            if i == tries - 1 or not (isinstance(ex, PermissionError) or getattr(ex, "winerror", None) in LOCKED):
                raise
            if i in (8, 80):
                log("waiting for a file held by another program:", dst, ex)
            time.sleep(0.25)


def write_json(path, data):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=1)
    replace(tmp, path)


def processes_in_folder():
    """Processes started from the MSpektra folder (its Python, MSpektra.exe), this one left out."""
    psapi = ctypes.windll.psapi
    arr = (W.DWORD * 4096)()
    needed = W.DWORD()
    if not psapi.EnumProcesses(ctypes.byref(arr), ctypes.sizeof(arr), ctypes.byref(needed)):
        return []
    me = os.getpid()
    root = os.path.normcase(ROOT) + os.sep
    out = []
    for pid in arr[:needed.value // ctypes.sizeof(W.DWORD)]:
        if not pid or pid == me:
            continue
        h = K32.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
        if not h:
            continue
        try:
            buf = ctypes.create_unicode_buffer(1024)
            size = W.DWORD(1024)
            if K32.QueryFullProcessImageNameW(h, 0, buf, ctypes.byref(size)):
                if os.path.normcase(buf.value).startswith(root):
                    out.append(pid)
        finally:
            K32.CloseHandle(h)
    return out


def image_name(pid):
    h = K32.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
    if not h:
        return "?"
    try:
        buf = ctypes.create_unicode_buffer(1024)
        size = W.DWORD(1024)
        return buf.value if K32.QueryFullProcessImageNameW(h, 0, buf, ctypes.byref(size)) else "?"
    finally:
        K32.CloseHandle(h)


def windows_of(pids):
    """{pid: [titles of its visible windows]} for the processes pids."""
    out = {p: [] for p in pids}
    proc = ctypes.WINFUNCTYPE(W.BOOL, W.HWND, W.LPARAM)

    def each(hwnd, lp):
        try:
            if U32.IsWindowVisible(hwnd):
                pid = W.DWORD()
                U32.GetWindowThreadProcessId(hwnd, ctypes.byref(pid))
                if pid.value in out:
                    n = U32.GetWindowTextLengthW(hwnd)
                    buf = ctypes.create_unicode_buffer(n + 1)
                    U32.GetWindowTextW(hwnd, buf, n + 1)
                    out[pid.value].append(buf.value or "MSpektra")
        except Exception:
            pass
        return True
    U32.EnumWindows(proc(each), 0)
    return out


def end_process(pid):
    """Ends a process of the MSpektra folder that has no window left (the MSpektra that asked for the update
    and did not end, or a helper of it): nothing of the user's is open in it."""
    h = K32.OpenProcess(0x0001 | 0x00100000, False, pid)  # PROCESS_TERMINATE | SYNCHRONIZE
    if not h:
        return False
    try:
        ok = bool(K32.TerminateProcess(h, 1))
        K32.WaitForSingleObject(h, 5000)
        return ok
    finally:
        K32.CloseHandle(h)


def wait_closed(grace=15.0, poll=0.5):
    """Waits until MSpektra and its other processes (deconvolution, other windows) have closed: their files
    are locked while they run. Processes of the MSpektra folder without any window (the MSpektra that asked for
    the update but did not end, a helper) are ended after grace seconds; MSpektra windows still open are named
    in a message (Retry: wait again). False when the user gives up (Cancel).
    (20.10.2026: on another PC the window closed, nothing came back and "Close every MSpektra window" showed,
    with no window to close: a process left without a window kept the installer waiting.)"""
    if PID:
        h = K32.OpenProcess(0x00100000, False, PID)  # SYNCHRONIZE
        if h:
            K32.WaitForSingleObject(h, 10000)
            K32.CloseHandle(h)
    t0 = time.time()
    while True:
        procs = processes_in_folder()
        if not procs:
            time.sleep(0.5)  # (files are released a moment after the process ends)
            return True
        wins = windows_of(procs)
        hidden = [p for p in procs if not wins.get(p)]
        shown = sorted({t for p in procs for t in wins.get(p) or []})
        waited = time.time() - t0
        if hidden and waited > grace:
            for p in hidden:
                log("ended a process without a window:", p, image_name(p), end_process(p))
            continue
        if shown and waited > grace:
            log("windows still open:", shown)
            if not message("These MSpektra windows are still open:\n\n    %s\n\nClose them, then press Retry to finish "
                           "the update. Cancel starts MSpektra again without the update (it is installed the next "
                           "time)." % "\n    ".join(shown[:8]), retry=True):
                return False
            t0 = time.time()
            continue
        time.sleep(poll)


def move(src, dst, journal):
    """Moves src to dst (paths in the MSpektra folder) and writes it to the journal."""
    s, d = os.path.join(ROOT, src), os.path.join(ROOT, dst)
    os.makedirs(os.path.dirname(d), exist_ok=True)
    replace(s, d)
    journal.append([src, dst])
    write_json(os.path.join(OLD, "journal.json"), journal)


def undo(journal):
    """Puts every moved file or folder back, last move first."""
    for src, dst in reversed(journal):
        s, d = os.path.join(ROOT, src), os.path.join(ROOT, dst)
        try:
            if os.path.exists(d):
                os.makedirs(os.path.dirname(s), exist_ok=True)
                if os.path.exists(s):  # (an unfinished move)
                    shutil.rmtree(s, ignore_errors=True) if os.path.isdir(s) else os.remove(s)
                replace(d, s)
        except OSError as ex:
            log("undo failed:", src, ex)


def paths(manifest, key):
    it = (manifest or {}).get("items", {}).get(key) or {}
    return list(it.get("paths") or [])


def install():
    ready = read_json(os.path.join(STAGING, "ready.json"))
    new = read_json(os.path.join(STAGING, "install.json"))
    cur = read_json(os.path.join(ROOT, "install.json"))
    if not ready or not new:
        log("nothing to install")
        return None
    if os.path.isdir(OLD):
        shutil.rmtree(OLD, ignore_errors=True)  # the version before the previous one
    os.makedirs(OLD, exist_ok=True)
    journal = []
    try:
        for key in list(ready["changed"]) + list(ready["removed"]):
            for p in paths(cur, key):
                if os.path.exists(os.path.join(ROOT, p)):
                    move(p, "_update/old/" + p, journal)
            if key in ready["changed"]:
                for p in paths(new, key):
                    move("_update/staging/" + p, p, journal)
        if os.path.isdir(os.path.join(ROOT, "_portable", "pylibs")):  # unpacked again from the wheels
            move("_portable/pylibs", "_update/old/_portable/pylibs", journal)
        move("install.json", "_update/old/install.json", journal)
        move("_update/staging/install.json", "install.json", journal)
    except Exception as ex:
        log("install failed:", repr(ex))
        undo(journal)
        message("The update to MSpektra %s could not be installed (%s). MSpektra %s is unchanged."
                % (ready["version"], ex, ready.get("from", "")))
        return False
    write_json(os.path.join(OLD, "info.json"), {"from": ready.get("from"), "to": ready["version"]})
    shutil.rmtree(STAGING, ignore_errors=True)
    log("installed %s (was %s): %s" % (ready["version"], ready.get("from"), ", ".join(ready["changed"])))
    return ready


def start_and_watch(version):
    """Starts MSpektra; True when its first window came up (pending.json removed), or it is still
    starting on a slow computer; False when it ended without a window."""
    write_json(PENDING, {"version": version, "time": time.time()})
    exe = os.path.join(ROOT, "MSpektra.exe")
    subprocess.Popen([exe], cwd=ROOT, close_fds=True)
    t0 = time.time()
    while time.time() - t0 < 180:
        time.sleep(1)
        if not os.path.exists(PENDING):
            return True
        if time.time() - t0 > 20 and not processes_in_folder():
            log("%s ended without a window" % version)
            return False
    alive = bool(processes_in_folder())
    log("%s not confirmed after 180 s (%s)" % (version, "still running" if alive else "not running"))
    return alive


def rollback(reason):
    """Back to the version in _update\\old; that version is not offered again."""
    info = read_json(os.path.join(OLD, "info.json")) or {}
    journal = read_json(os.path.join(OLD, "journal.json"))
    if not journal:
        message("There is no previous version to go back to.")
        return False
    undo(journal)
    shutil.rmtree(STAGING, ignore_errors=True)
    shutil.rmtree(OLD, ignore_errors=True)
    write_json(os.path.join(UPD, "skip.json"), {"version": info.get("to"), "reason": reason})
    try:
        os.remove(PENDING)
    except OSError:
        pass
    log("back to %s (%s)" % (info.get("from"), reason))
    return info.get("from")


def main():
    if not wait_closed():
        log("cancelled: MSpektra still open; started again without the update (it stays ready)")
        try:
            subprocess.Popen([os.path.join(ROOT, "MSpektra.exe")], cwd=ROOT, close_fds=True)
        except OSError as ex:
            log("could not start MSpektra:", ex)
        return
    if "--rollback" in ARGS:
        if rollback("chosen in About"):
            subprocess.Popen([os.path.join(ROOT, "MSpektra.exe")], cwd=ROOT, close_fds=True)
        return
    ready = install()
    if not ready:
        if ready is False:
            subprocess.Popen([os.path.join(ROOT, "MSpektra.exe")], cwd=ROOT, close_fds=True)
        return
    if start_and_watch(ready["version"]):
        return
    back = rollback("did not start")
    subprocess.Popen([os.path.join(ROOT, "MSpektra.exe")], cwd=ROOT, close_fds=True)
    message("MSpektra %s did not start, so MSpektra %s was put back." % (ready["version"], back or "the previous"))


if __name__ == "__main__":
    try:
        main()
    except Exception as ex:
        log("error:", repr(ex))
    finally:
        try:
            os.remove(os.path.abspath(__file__))  # the copy in the temp folder
        except OSError:
            pass
