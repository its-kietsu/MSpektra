@echo off
setlocal
set "UD_ROOT=%~dp0"
title MSpektra (console)
"%UD_ROOT%python\python.exe" -s "%UD_ROOT%_portable\launch_unidec.py" %*
echo.
echo MSpektra exited with code %errorlevel%. Press any key to close this window.
pause >nul
