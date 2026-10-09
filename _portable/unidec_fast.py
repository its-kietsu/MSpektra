"""
Start-up speed add-on for the portable UniDec package.

UniDec imports its Thermo .raw reader at start, which boots the .NET runtime
(pythonnet) even when no Thermo file is ever opened. This replaces that module
with a small stand-in that loads the real reader the first time a .raw file is
opened. Nothing else changes.
"""
import importlib
import sys
import time
import types

_THERMO = "unidec.UniDecImporter.Thermo.Thermo"


def install_lazy_thermo():
    if _THERMO in sys.modules:
        return False
    stub = types.ModuleType(_THERMO)
    stub.__file__ = "<Thermo reader, loaded on first use (portable add-on)>"
    cache = {}

    def ThermoImporter(*args, **kwargs):
        cls = cache.get("cls")
        if cls is None:
            if sys.modules.get(_THERMO) is stub:
                del sys.modules[_THERMO]
            t = time.perf_counter()
            module = importlib.import_module(_THERMO)
            cls = cache["cls"] = module.ThermoImporter
            print("Thermo reader (.NET) loaded on first use in %.1f s" % (time.perf_counter() - t))
        return cls(*args, **kwargs)

    stub.ThermoImporter = ThermoImporter
    sys.modules[_THERMO] = stub
    return True
