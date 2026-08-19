"""
launcher.py — makes `pharmol-gui` (or whatever entry-point name is configured
in pyproject.toml) launch PyMOL with PhArMol already open, in one command —
matching the UX of tools like OpenPharmaco's `openph` command.

This does exactly two things, in order:
  1. Launches PyMOL's real, normal windowed GUI (pymol.finish_launching)
  2. Registers this plugin's menu item and immediately opens its dialog

After that, PyMOL behaves completely normally — the plugin is also still
available from Plugin > PhArMol for the rest of the session, exactly as if
it had been loaded through the Plugin Manager.
"""
import sys


def main():
    import pymol
    # Launch PyMOL's actual windowed GUI (not -cq/quiet mode) — this is a
    # real, normal PyMOL session, not a special headless mode.
    pymol.finish_launching(['pymol'])

    # Import here, after PyMOL/Qt are fully initialised.
    from . import __init_plugin__, run_plugin_gui

    __init_plugin__()   # registers the menu item, so it's there for later too
    run_plugin_gui()    # and open the dialog immediately


if __name__ == "__main__":
    sys.exit(main())
