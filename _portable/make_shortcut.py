"""Create 'MSpektra' shortcuts (desktop and Start menu) for this portable
folder.

The shortcuts point to MSpektra.exe and carry the same taskbar identity
(AppUserModelID) as the running windows, so a pinned shortcut and the open
windows share one taskbar icon. Older 'UniDec' shortcuts made by earlier
versions are replaced.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Only MSpektra.exe in the root is used. The old launcher UniDec.exe no longer
# lives in the root (it is kept in _portable\tools); it finds python\ next to
# itself, so from there it cannot start the program and is never a target.
EXE = os.path.join(ROOT, "MSpektra.exe")
PYTHONW = os.path.join(ROOT, "python", "pythonw.exe")
LAUNCHER = os.path.join(ROOT, "_portable", "launch_unidec.py")
ICON = os.path.join(ROOT, "_portable", "msanalysis.ico")
APP_ID = "UniDec.Portable.Launcher"  # must match unidec_ui_addons.APP_ID


def make_link(path):
    import pythoncom
    from win32com.shell import shell
    from win32com.propsys import propsys, pscon

    link = pythoncom.CoCreateInstance(shell.CLSID_ShellLink, None, pythoncom.CLSCTX_INPROC_SERVER,
                                      shell.IID_IShellLink)
    if os.path.isfile(EXE):
        link.SetPath(EXE)
        link.SetIconLocation(EXE, 0)
    else:  # MSpektra.exe missing (blocked or deleted): start the launcher script
        link.SetPath(PYTHONW)
        link.SetArguments('-s "%s"' % LAUNCHER)
        link.SetIconLocation(ICON, 0)
    link.SetWorkingDirectory(ROOT)
    link.SetDescription("MSpektra")
    store = link.QueryInterface(propsys.IID_IPropertyStore)
    store.SetValue(pscon.PKEY_AppUserModel_ID, propsys.PROPVARIANTType(APP_ID, pythoncom.VT_LPWSTR))
    store.Commit()
    link.QueryInterface(pythoncom.IID_IPersistFile).Save(path, 0)


def main():
    from win32com.shell import shell, shellcon
    made = []
    for csidl, label in ((shellcon.CSIDL_DESKTOPDIRECTORY, "desktop"), (shellcon.CSIDL_PROGRAMS, "Start menu")):
        try:
            folder = shell.SHGetFolderPath(0, csidl, None, 0)
            path = os.path.join(folder, "MSpektra.lnk")
            make_link(path)
            old = os.path.join(folder, "UniDec.lnk")
            if os.path.isfile(old):
                try:
                    os.remove(old)
                except OSError:
                    pass
            made.append(path)
            print("Shortcut created (%s): %s" % (label, path))
        except Exception as e:
            print("Could not create the %s shortcut: %s" % (label, e))
    if made:
        print("\nTo pin MSpektra to the taskbar: start it from one of these shortcuts,")
        print("right-click its taskbar icon and choose 'Pin to taskbar'.")
        print("If you move the folder, run this again.")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
