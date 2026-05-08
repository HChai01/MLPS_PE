#!/usr/bin/env python3
"""
Script Name:  export_command_md.py
Description:  Export Nextflow task command scripts into a readable Markdown report.
Author:       Dr. Haiting Chai
Last update:  2026-APR-12
Version:      ver.2026-04-12
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from assets.utils_banner import print_header, print_summary, setup_logging


logger = setup_logging()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export Nextflow .command.sh files into one Markdown document.",
    )
    parser.add_argument("--trace", type=Path, required=True, help="Nextflow trace.txt path")
    parser.add_argument("--workdir", type=Path, required=True, help="Nextflow work directory")
    parser.add_argument("--outputmd", type=Path, required=True, help="Output Markdown path")
    parser.add_argument("--project", required=True, help="Project name for the report title")
    return parser


def read_trace_rows(trace_path: Path) -> list[dict[str, str]]:
    if not trace_path.exists():
        logger.warning("Trace file not found: %s", trace_path)
        return []

    with trace_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [row for row in reader if row]

    if not rows:
        logger.warning("Trace file had no task rows: %s", trace_path)
    return rows


def build_trace_entries(rows: list[dict[str, str]], workdir: Path) -> list[dict[str, str | Path]]:
    entries: list[dict[str, str | Path]] = []
    for row in rows:
        task_hash = (row.get("hash") or "").strip()
        if "/" not in task_hash:
            continue

        shard, task_id = task_hash.split("/", 1)
        task_dir = workdir / shard / task_id
        command_path = task_dir / ".command.sh"
        if not command_path.exists():
            continue

        entries.append(
            {
                "task_id": (row.get("task_id") or "").strip(),
                "hash": task_hash,
                "name": (row.get("name") or "Unknown task").strip(),
                "status": (row.get("status") or "UNKNOWN").strip(),
                "exit": (row.get("exit") or "-").strip(),
                "submit": (row.get("submit") or "-").strip(),
                "duration": (row.get("duration") or "-").strip(),
                "task_dir": task_dir,
                "command_path": command_path,
            }
        )
    return entries


def build_fallback_entries(workdir: Path) -> list[dict[str, str | Path]]:
    command_paths = sorted(
        workdir.glob("*/*/.command.sh"),
        key=lambda path: path.stat().st_mtime,
    )
    entries: list[dict[str, str | Path]] = []
    for index, command_path in enumerate(command_paths, start=1):
        task_dir = command_path.parent
        entries.append(
            {
                "task_id": str(index),
                "hash": f"{task_dir.parent.name}/{task_dir.name}",
                "name": "Unknown task",
                "status": "UNKNOWN",
                "exit": ((task_dir / ".exitcode").read_text(encoding="utf-8").strip() if (task_dir / ".exitcode").exists() else "-"),
                "submit": "-",
                "duration": "-",
                "task_dir": task_dir,
                "command_path": command_path,
            }
        )
    return entries


def load_command_text(command_path: Path) -> str:
    return command_path.read_text(encoding="utf-8").rstrip() or "# Empty .command.sh"


def write_markdown(output_path: Path, project_name: str, entries: list[dict[str, str | Path]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        f"# {project_name} command audit",
        "",
        "This file collects the Nextflow task `.command.sh` scripts that were available when the workflow ended.",
        "",
        f"Tasks captured: {len(entries)}",
        "",
    ]

    for entry in entries:
        command_text = load_command_text(entry["command_path"])  # type: ignore[index]
        lines.extend(
            [
                f"## Task {entry['task_id']}: {entry['name']}",
                "",
                f"- Status: `{entry['status']}`",
                f"- Exit code: `{entry['exit']}`",
                f"- Submitted: `{entry['submit']}`",
                f"- Duration: `{entry['duration']}`",
                f"- Hash: `{entry['hash']}`",
                f"- Work dir: `{entry['task_dir']}`",
                "",
                "```bash",
                command_text,
                "```",
                "",
            ]
        )

    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = build_parser().parse_args()

    trace_rows = read_trace_rows(args.trace)
    entries = build_trace_entries(trace_rows, args.workdir)
    if not entries:
        logger.warning("Falling back to scanning work directories for .command.sh files")
        entries = build_fallback_entries(args.workdir)

    write_markdown(args.outputmd, args.project, entries)
    logger.info("Command audit written to %s", args.outputmd)


if __name__ == "__main__":
    start_time = time.time()
    print_header("export_command_md.py")
    main()
    print_summary(start_time, "export_command_md.py")
