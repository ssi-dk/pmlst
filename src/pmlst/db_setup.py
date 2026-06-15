#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path

from .database import default_database_destination

DEFAULT_REPO_URL = "https://bitbucket.org/genomicepidemiology/pmlst_db.git"
DEFAULT_BRANCH = "master"
VERSION_FILE = "version.txt"


def run_command(command: list[str], dry_run: bool = False) -> str:
    if dry_run:
        print("+ " + " ".join(command))
        return ""

    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def ensure_database(
    destination: Path, repo_url: str, branch: str, dry_run: bool
) -> None:
    if destination.exists():
        if (destination / ".git").is_dir():
            run_command(
                ["git", "-C", str(destination), "fetch", "origin", branch], dry_run
            )
            run_command(["git", "-C", str(destination), "checkout", branch], dry_run)
            run_command(
                ["git", "-C", str(destination), "pull", "--ff-only", "origin", branch],
                dry_run,
            )
            return

        if any(destination.iterdir()):
            raise SystemExit(
                f"Destination exists and is not an empty directory or git checkout: "
                f"{destination}"
            )

    if dry_run:
        if not destination.parent.exists():
            print(f"Would create directory: {destination.parent}")
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)

    run_command(
        [
            "git",
            "clone",
            "--branch",
            branch,
            "--single-branch",
            repo_url,
            str(destination),
        ],
        dry_run,
    )


def install_database(destination: Path, dry_run: bool) -> None:
    install_script = destination / "INSTALL.py"
    if not dry_run and not install_script.exists():
        raise SystemExit(f"Database INSTALL.py was not found: {install_script}")

    run_command([sys.executable, str(install_script), "kma_index"], dry_run)


def commit_metadata(destination: Path, dry_run: bool) -> tuple[str, str]:
    commit_date = run_command(
        ["git", "-C", str(destination), "log", "-1", "--format=%cs"],
        dry_run,
    )
    commit_hash = run_command(
        ["git", "-C", str(destination), "rev-parse", "HEAD"],
        dry_run,
    )
    if dry_run:
        commit_date = "<latest-commit-date>"
        commit_hash = "<latest-commit-hash>"

    return commit_date, commit_hash


def write_version_file(
    destination: Path,
    repo_url: str,
    branch: str,
    commit_date: str,
    commit_hash: str,
    dry_run: bool,
) -> Path:
    version_file = destination / VERSION_FILE
    content = "\n".join(
        [
            f"version: {commit_date}",
            f"commit: {commit_hash}",
            f"source: {repo_url}",
            f"branch: {branch}",
            "",
        ]
    )

    if dry_run:
        print(f"Would write {version_file}:")
        print(content, end="")
        return version_file

    version_file.write_text(content, encoding="utf-8")
    return version_file


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download or update the pMLST database after tool installation."
    )
    parser.add_argument(
        "destination",
        nargs="?",
        type=Path,
        help=(
            "Database destination directory. Defaults to PMLST_DB, then "
            "$CONDA_PREFIX/share/pmlst/db, then the platformdirs user data path."
        ),
    )
    parser.add_argument(
        "--repo-url",
        default=DEFAULT_REPO_URL,
        help="Upstream pMLST database git repository URL.",
    )
    parser.add_argument(
        "--branch",
        default=DEFAULT_BRANCH,
        help="Upstream pMLST database branch to clone or update.",
    )
    parser.add_argument(
        "--skip-install",
        action="store_true",
        help=(
            "Fetch/update the database and write version.txt without running "
            "INSTALL.py."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned actions without changing files.",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    destination = (args.destination or default_database_destination()).expanduser()
    destination = destination.resolve()

    ensure_database(destination, args.repo_url, args.branch, args.dry_run)
    if not args.skip_install:
        install_database(destination, args.dry_run)

    commit_date, commit_hash = commit_metadata(destination, args.dry_run)
    version_file = write_version_file(
        destination,
        args.repo_url,
        args.branch,
        commit_date,
        commit_hash,
        args.dry_run,
    )

    print(f"pMLST database is ready at {destination}")
    print(f"Database version: {commit_date}")
    print(f"Version file: {version_file}")


if __name__ == "__main__":
    main()
