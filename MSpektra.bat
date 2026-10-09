@echo off
setlocal
set "UD_ROOT=%~dp0"
start "" "%UD_ROOT%python\pythonw.exe" -s "%UD_ROOT%_portable\launch_unidec.py" %*
