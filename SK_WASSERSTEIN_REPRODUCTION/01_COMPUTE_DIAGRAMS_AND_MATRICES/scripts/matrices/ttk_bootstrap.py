#!/usr/bin/env python3


from __future__ import annotations

import os
from pathlib import Path

_LOADED = False


def ensure_ttk_plugin_loaded() -> None:
    global _LOADED
    if _LOADED:
        return

    plugin = os.environ.get("TTK_PLUGIN", "").strip()
    if plugin:
        plugin_path = Path(plugin).expanduser().resolve()
        if not plugin_path.is_file():
            raise FileNotFoundError(f"TTK_PLUGIN points to a missing file: {plugin_path}")
        import paraview.simple as pvs

        try:
            pvs.LoadPlugin(str(plugin_path), remote=False, ns=pvs.__dict__)
        except TypeError:
            pvs.LoadPlugin(str(plugin_path), remote=False)

    _LOADED = True
