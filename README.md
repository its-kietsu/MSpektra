# MSpektra

Windows program for LC-MS and HRMS analysis with deconvolution.

- **LCMS Analysis:** Shimadzu LabSolutions .lcd files (MS and PDA), peak integration, Compare view
- **HRMS Analysis:** Bruker .d folders and mzML, calibration, exact mass, formula finder
- **Deconvolution:** Bayesian (UniDec), maximum entropy and IsoDec

## Download

Get `MSpektra-<version>.zip` from [Releases](../../releases), extract it to a short folder path (for example the Desktop) and start `MSpektra.exe`. Nothing is installed; the zip includes its own Python. See `README.txt` for the manual.

## Deconvolution engine: UniDec

The deconvolution uses UniDec by Michael T. Marty (https://github.com/michaelmarty/UniDec), partly compiled into `msengine.dll` from UniDec's C source. If you publish results obtained with it, please cite:

> Marty et al., Anal. Chem. 2015, 87 (8), 4370-4376. DOI: 10.1021/acs.analchem.5b00140

## Building from source

Toolchain: MinGW-w64 GCC (WinLibs, posix threads, msvcrt), CMake and Ninja.

```
cmake -S cpp\msengine -B cpp\msengine\build-win -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ "-DCMAKE_SHARED_LINKER_FLAGS=-static -static-libgcc -static-libstdc++ -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive"
cmake --build cpp\msengine\build-win
cpp\maxent_native\build.ps1
cpp\kinetics\build.ps1
cpp\polymer\build.ps1
```

Copy the DLLs into `_portable\msengine\`. To run from source, copy `python\`, `_portable\baf2sql\` and `_portable\timsdata\` from a release zip next to the sources and start `MSpektra.bat`.

## License

MSpektra is released under the [MIT License](LICENSE). Third-party components keep their own licenses (see `LICENSES\`): UniDec (BSD-style, including the code compiled from it in `cpp\msengine\src\udengine` and `unidec.cpp`), Python, OpenSZRaw, olefile, the Bruker, Thermo and Waters libraries.
