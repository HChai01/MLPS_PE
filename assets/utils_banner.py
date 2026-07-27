#!/usr/bin/env python3
"""
Script Name:  utils_banner.py
Description:  Utils banners
Author:       Dr. Haiting Chai
Last update:  2025-OCT-28
Version:      ver.2025-10-28

DEPENDENCIES:
    - Python 3.10
    - sys (standard library)
    - time (standard library)
    - logging (standard library)
    - socket (standard library)
    - getpass (standard library)
    - pathlib (standard library)
"""

from __future__ import annotations

# ========================================================================================
# Imports
# ========================================================================================

import getpass
import logging
import socket
import sys
import time
from pathlib import Path

# ========================================================================================
# Constants
# ========================================================================================

class UBColors:
    """ANSI escape sequences used by the shared CLI logging helpers."""

    NC = "\033[0m"
    BLACK = "\033[0;30m"
    RED = "\033[0;31m"
    GREEN = "\033[0;32m"
    YELLOW = "\033[0;33m"
    BLUE = "\033[0;34m"
    MAGENTA = "\033[0;35m"
    CYAN = "\033[0;36m"
    WHITE = "\033[0;37m"
    GREY = "\033[1;30m"

# ========================================================================================
# Helper Functions
# ========================================================================================

def setup_logging(level: int = logging.INFO) -> logging.Logger:
    """Configure a consistent logger for the project scripts."""
    root = logging.getLogger()
    if root.handlers:
        for handler in list(root.handlers):
            root.removeHandler(handler)
    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger(__name__)


def get_project_root(current_file: str | Path) -> Path:
    """Return the repository root for a script living under `python/` or `assets/`."""
    return Path(current_file).resolve().parents[1]


def print_header(script_name: str = "MLPS utility script") -> None:
    """Log a standard execution banner at the start of a helper script."""
    logging.info("=" * 60)
    logging.info("        %s", script_name)
    logging.info("=" * 60)
    logging.info("User:        %s", getpass.getuser())
    logging.info("Machine:     %s", socket.gethostname())
    logging.info("Python ver.: %s", sys.version.split()[0])
    logging.info("Working dir: %s", Path.cwd())
    logging.info("=" * 60)


def print_summary(start_time: float, script_name: str = "MLPS utility script") -> None:
    """Log the elapsed runtime summary for a helper script."""
    runtime = time.time() - start_time
    logging.info("")
    logging.info("=" * 60)
    logging.info("%s completed in %.2f seconds", script_name, runtime)
    logging.info("=" * 60)


def parse_yes_no(value: str | bool, *, param_name: str = "value") -> bool:
    """Parse yes/no style CLI flags used across the legacy helper scripts."""
    if isinstance(value, bool):
        return value
    normalised = str(value).strip().lower()
    if normalised in {"yes", "true", "1"}:
        return True
    if normalised in {"no", "false", "0"}:
        return False
    raise ValueError(f"{param_name} must be yes/no or true/false, received: {value}")


# ========================================================================================
# Public Exports
# ========================================================================================

__all__ = [
    "UBColors",
    "get_project_root",
    "parse_yes_no",
    "print_header",
    "print_summary",
    "setup_logging",
]
