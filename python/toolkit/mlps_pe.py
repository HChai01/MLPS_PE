#!/usr/bin/env python3
"""
Script Name:  mlps_pe.py
Description:  Dispatch registered MLPS Python toolkit modules from one CLI.
Author:       Dr. Haiting Chai
Last update:  2026-AUG-04
Version:      ver.2026-08-04

DEPENDENCIES:
    - Python 3.10
    - __future__ (standard library)
    - csv (standard library)
    - os (standard library)
    - pathlib (standard library)
    - subprocess (standard library)
    - sys (standard library)
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path


# ========================================================================================
# Configuration
# ========================================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]
TOOLKIT_DIR = Path(__file__).resolve().parent
TOOLKIT_REGISTRY = TOOLKIT_DIR / "toolkit.list.csv"
TOOLKIT_SCRIPT_DIR = PROJECT_ROOT / "python"


# ========================================================================================
# Helper Functions
# ========================================================================================


def load_modules() -> dict[str, dict[str, str]]:
    """Load and validate module definitions from the toolkit registry CSV."""
    if not TOOLKIT_REGISTRY.is_file():
        raise SystemExit(f"Toolkit registry not found: {TOOLKIT_REGISTRY}")

    with TOOLKIT_REGISTRY.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))

    modules: dict[str, dict[str, str]] = {}
    for row in rows:
        name = row.get("module", "").strip()
        script = row.get("script", "").strip()
        if name and script:
            modules[name] = row
    return modules


def print_help(modules: dict[str, dict[str, str]]) -> None:
    """Print the top-level command help and all registered module summaries."""
    print("MLPS_PE toolkit")
    print("Run the registered MLPS Python utilities from one terminal command.\n")
    print("Usage:")
    print("  mlps_pe <module_name> [module_parameters]")
    print("  mlps_pe --help")
    print("  mlps_pe <module_name> --help\n")
    print("Available modules:")
    for name, row in modules.items():
        print(f"\n  {name}")
        print(f"      {row.get('description', '').strip()}")
        dependencies = row.get("external_dependency", "na").strip() or "na"
        print(f"      Dependencies: {dependencies}")
        parameters = row.get("parameters", "na").strip() or "na"
        print(f"      Parameters: {parameters.replace(';', '; ')}")
    print("\nFor complete option details, run: mlps_pe <module_name> --help")


def main(argv: list[str] | None = None) -> int:
    """Dispatch one registered module or print the top-level toolkit help."""
    modules = load_modules()
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args in (["-h"], ["--help"]):
        print_help(modules)
        return 0
    module_name = args.pop(0)
    if module_name not in modules:
        print(f"Unknown module: {module_name}\n", file=sys.stderr)
        print_help(modules)
        return 2

    # Resolve the script from the registry
    script = TOOLKIT_SCRIPT_DIR / modules[module_name]["script"]
    if not script.is_file():
        print(f"Registered script not found: {script}", file=sys.stderr)
        return 2

    # Pass the project location to child processes for future runtime-aware modules
    environment = os.environ.copy()
    environment["MLPS_PE_ROOT"] = str(PROJECT_ROOT)
    completed = subprocess.run([sys.executable, str(script), *args], env=environment)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
