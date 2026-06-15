import os
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from .utils import PmlstError

ROOT_DATABASE_FILES = ("config", "INSTALL.py", "version.txt")
SCHEME_DATABASE_FILES = ("fsa", "txt.clean", "clpx", "name")
KMA_DATABASE_FILES = ("comp.b", "length.b", "seq.b")


@dataclass(frozen=True)
class DatabaseCandidate:
    source: str
    path: Path


def conda_database_path() -> Path | None:
    conda_prefix = os.getenv("CONDA_PREFIX")
    if not conda_prefix:
        return None
    return Path(conda_prefix) / "share" / "pmlst" / "db"


def user_database_path() -> Path:
    # Keep this lazy so --help/--version work in lightweight no-deps smoke checks.
    from platformdirs import user_data_path  # noqa: PLC0415

    return user_data_path("pmlst") / "db"


def default_database_destination() -> Path:
    env_database = os.getenv("PMLST_DB")
    if env_database:
        return Path(env_database)

    conda_database = conda_database_path()
    if conda_database is not None:
        return conda_database

    return user_database_path()


def resolve_database_path(database_option: str | None) -> Path:
    if database_option:
        return _resolve_authoritative_database(
            DatabaseCandidate("-p/--database", Path(database_option))
        )

    env_database = os.getenv("PMLST_DB")
    if env_database:
        return _resolve_authoritative_database(
            DatabaseCandidate("PMLST_DB", Path(env_database))
        )

    checked = default_database_candidates()
    for candidate in checked:
        if candidate.path.exists():
            return candidate.path.expanduser().resolve()

    raise PmlstError(
        "No usable pMLST database found. "
        f"Checked: {_format_checked_paths(checked)}. "
        "Run pmlst-download-db DESTINATION, set PMLST_DB, or pass -p/--database."
    )


def default_database_candidates() -> list[DatabaseCandidate]:
    candidates: list[DatabaseCandidate] = []
    conda_database = conda_database_path()
    if conda_database is not None:
        candidates.append(
            DatabaseCandidate("$CONDA_PREFIX/share/pmlst/db", conda_database)
        )
    candidates.append(DatabaseCandidate("platformdirs user data", user_database_path()))
    return candidates


def validate_database_root(database: Path) -> None:
    missing_files = _missing_files(database, ROOT_DATABASE_FILES)
    if missing_files:
        raise PmlstError(
            f"Database at {database} is incomplete. "
            f"Missing root files: {', '.join(missing_files)}. "
            f"Run pmlst-download-db {database} to install or update the database."
        )


def validate_database_schemes(
    database: Path, schemes: Sequence[str], require_kma: bool
) -> None:
    for scheme in schemes:
        missing_scheme_files = _missing_scheme_files(
            database, scheme, SCHEME_DATABASE_FILES
        )
        if missing_scheme_files:
            raise PmlstError(
                f"Database at {database} is incomplete for scheme '{scheme}'. "
                f"Missing files: {', '.join(missing_scheme_files)}. "
                f"Run pmlst-download-db {database} to install or update the database."
            )

        if require_kma:
            missing_kma_files = _missing_scheme_files(
                database, scheme, KMA_DATABASE_FILES
            )
            if missing_kma_files:
                raise PmlstError(
                    f"Database at {database} is incomplete for FASTQ/KMA scheme "
                    f"'{scheme}'. Missing files: {', '.join(missing_kma_files)}. "
                    f"Run pmlst-download-db {database} to install or update the "
                    "database."
                )


def _resolve_authoritative_database(candidate: DatabaseCandidate) -> Path:
    database = candidate.path.expanduser()
    if not database.exists():
        raise PmlstError(
            f"Database path from {candidate.source} does not exist: {database}. "
            f"Run pmlst-download-db {database} to install the database, "
            "or provide a different path."
        )
    return database.resolve()


def _missing_files(directory: Path, filenames: Sequence[str]) -> list[str]:
    return [filename for filename in filenames if not (directory / filename).exists()]


def _missing_scheme_files(
    database: Path, scheme: str, suffixes: Sequence[str]
) -> list[str]:
    filenames = [f"{scheme}.{suffix}" for suffix in suffixes]
    return _missing_files(database, filenames)


def _format_checked_paths(candidates: Sequence[DatabaseCandidate]) -> str:
    return "; ".join(
        f"{candidate.source}: {candidate.path}" for candidate in candidates
    )
