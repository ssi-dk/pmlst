#!/usr/bin/env python3

import argparse
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any


def version_from_database(database: Path) -> str:
    version_file = database / "version.txt"
    if not version_file.exists():
        return "unknown"

    for line in version_file.read_text(encoding="utf-8").splitlines():
        if line.startswith("version:"):
            return line.split(":", 1)[1].strip()
    return "unknown"


def install_database(destination: Path, fixture_db: Path | None) -> None:
    if fixture_db is not None:
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(fixture_db, destination)
        return

    subprocess.run(
        ["pmlst-download-db", str(destination)],
        check=True,
    )


def write_data_manager_json(
    out_file: Path,
    db_key: str,
    db_name: str,
    database: Path,
) -> None:
    data_manager_json: dict[str, Any] = {
        "data_tables": {
            "pmlst_db": [
                {
                    "value": db_key,
                    "name": db_name,
                    "path": str(database),
                    "version": version_from_database(database),
                }
            ]
        }
    }
    out_file.write_text(json.dumps(data_manager_json, sort_keys=True), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Install a pMLST database for Galaxy.")
    parser.add_argument("--out-file", required=True, type=Path)
    parser.add_argument("--db-key", required=True)
    parser.add_argument("--db-name", required=True)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument(
        "--fixture-db",
        type=Path,
        help="Copy this local fixture database instead of downloading upstream data.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    destination = args.destination.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    fixture_db = args.fixture_db.resolve() if args.fixture_db else None
    install_database(destination, fixture_db)
    write_data_manager_json(args.out_file, args.db_key, args.db_name, destination)


if __name__ == "__main__":
    main()
