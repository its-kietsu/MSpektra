@echo off
setlocal
title JCAMP-DX to UniDec text converter
rem Drag .jdx files (or a folder containing them) onto this file.
rem Each file is saved as <name>_unidec.txt (tab separated m/z, intensity)
rem next to the original. UniDec itself can also open .jdx files directly.
"%~dp0python\python.exe" -s "%~dp0_portable\unidec_jcamp.py" %*
echo.
pause
