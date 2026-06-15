import sys
from pathlib import Path


def _remove_script_dir_from_path() -> None:
    script_dir = Path(sys.argv[0]).resolve().parent
    if sys.path and Path(sys.path[0]).resolve() == script_dir:
        sys.path.pop(0)


def main() -> None:
    _remove_script_dir_from_path()

    from pmlst.cli import main as package_main

    package_main()


def db_setup_main() -> None:
    _remove_script_dir_from_path()

    from pmlst.db_setup import main as package_main

    package_main()


if __name__ == "__main__":
    main()
